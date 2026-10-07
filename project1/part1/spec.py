# -*- coding: utf-8 -*-
"""Part 1 스펙 시트 8항목 측정 하네스.

과제 29쪽 표를 그대로 구현한다. 추출 규칙을 여기 한 곳에 모아 두어야
탐색·ablation·최종 보고가 같은 숫자를 말한다.

  Vth        constant-current Id = 1e-7 A/um, Vd = 0.05 V      0.45 +- 0.05 V
  Ion        Vg = Vd = 2 V                                     >= 450 uA/um
  Ioff       Vg = 0, Vd = 2 V, 300 K                           <= 1 pA/um
  SS         Vd = 0.05 V, Id 1e-10~1e-8 구간 평균               <= 75 mV/dec
  고온 Ioff   Vg = 0, Vd = 2 V, 398 K                           <= 100 pA/um
  DIBL       (Vth@Vd0.05 - Vth@Vd2) / 1.95 V                   <= 30 mV/V
  Body       Vth 상승, Vb = -0.5 V                              <= 0.08 V
  E_ox       VDD / tox                                          <= 5 MV/cm
"""

from __future__ import annotations

import contextlib
import dataclasses
import io
import math
import time

import numpy as np

from .device import Part1Device
from .simulator import Part1Simulator

VDD = 2.0
VD_LIN = 0.05
VB_BODY = -0.5
T_COLD = 300.0
T_HOT = 398.0

VTH_CURRENT = 1.0e-7      # A/um
SS_LO, SS_HI = 1.0e-10, 1.0e-8
STEP_COARSE = 0.1
STEP_FINE = 0.01
BRACKET_MAX = 0.0105      # 이보다 넓은 구간에서 보간한 Vth 는 신뢰하지 않는다
SS_MIN_POINTS = 5
# 배정밀도 바닥. 이보다 작은(또는 음수인) Ioff 는 측정값이 아니라 반올림 잡음이다.
# 128비트 확장 정밀도를 켜도 바닥 자체는 남으므로 값이 아니라 게이트로 막는다.
IOFF_FLOOR_pA = 1.0e-4   # = 1e-16 A/um

LIMITS = {
    "vth": (0.40, 0.50),
    "ion_uA_per_um": (450.0, None),
    "ioff_pA_per_um": (None, 1.0),
    "ss_mV_per_dec": (None, 75.0),
    "ioff_hot_pA_per_um": (None, 100.0),
    "dibl_mV_per_V": (None, 30.0),
    "body_effect_mV": (None, 80.0),
    "eox_MV_per_cm": (None, 5.0),
}


# ---------------------------------------------------------------------------
# 추출기
# ---------------------------------------------------------------------------
def vg_at_current(vg, idd, target):
    """log10(Id) 선형 보간으로 target 전류가 흐르는 Vg 를 찾는다.

    Id 는 지수적으로 변하므로 선형 보간은 서브스레숄드에서 수십 mV 를 틀린다.
    반환값은 (Vg, 보간에 쓴 구간 폭).
    """
    vg = np.asarray(vg, dtype=float)
    idd = np.asarray(idd, dtype=float)
    order = np.argsort(vg)
    vg, idd = vg[order], idd[order]
    good = idd > 0
    vg, idd = vg[good], idd[good]
    if vg.size < 2 or idd[0] >= target or idd[-1] < target:
        return float("nan"), float("nan")
    k = int(np.argmax(idd >= target))
    l0, l1 = math.log10(idd[k - 1]), math.log10(idd[k])
    lt = math.log10(target)
    if l1 == l0:
        return float(vg[k]), float(vg[k] - vg[k - 1])
    return (float(vg[k - 1] + (lt - l0) * (vg[k] - vg[k - 1]) / (l1 - l0)),
            float(vg[k] - vg[k - 1]))


def subthreshold_swing(vg, idd):
    """Id 1e-10~1e-8 구간의 SS [mV/dec].

    "평균"의 정의를 두 decade 양 끝점 사이 기울기로 고정한다. 점별 차분의 평균은
    격자 잡음에 민감하고, 회귀는 구간 밖 점이 섞이면 흔들린다. 양 끝점 방식은
    정의가 유일하고 재현된다. 회귀값은 교차 확인용으로만 함께 돌려준다.
    """
    lo, _ = vg_at_current(vg, idd, SS_LO)
    hi, _ = vg_at_current(vg, idd, SS_HI)
    vg = np.asarray(vg, dtype=float)
    idd = np.asarray(idd, dtype=float)
    mask = (idd >= SS_LO) & (idd <= SS_HI)
    n = int(mask.sum())
    ss = (hi - lo) / 2.0 * 1e3 if math.isfinite(hi) and math.isfinite(lo) else float("nan")
    ss_fit = float("nan")
    if n >= 2:
        ss_fit = float(np.polyfit(np.log10(idd[mask]), vg[mask], 1)[0] * 1e3)
    return ss, ss_fit, n


# ---------------------------------------------------------------------------
# 스윕
# ---------------------------------------------------------------------------
def _scan(sim, vgs):
    """주어진 Vg 목록을 순서대로 인가하며 드레인 전류를 읽는다."""
    out = []
    for vg in vgs:
        sim.set_bias("gate", float(vg))
        out.append(sim.drain_current())
    return out


def _grid(lo, hi, step):
    n = max(1, int(round(abs(hi - lo) / step)))
    return [round(lo + (hi - lo) * i / n, 6) for i in range(n + 1)]


def sweep_for_vth(sim, vd, vg_hi=VDD, descending=False, extra_points=()):
    """Vth·SS 추출용 2단 스윕. 굵게 훑어 구간을 잡고, 그 근처만 촘촘히 다시 잰다.

    descending=True 면 게이트를 먼저 vg_hi 까지 올리고 드레인을 걸고 내려온다.
    Vd = 2 V 에서는 이 순서가 아니면 수렴하지 않는다 (off 상태에서 드레인을
    올리면 채널이 없어 공핍층이 급격히 퍼진다).
    """
    if descending:
        sim.set_bias("gate", vg_hi)
        sim.set_bias("drain", vd)
        coarse_vgs = _grid(vg_hi, 0.0, STEP_COARSE)
    else:
        sim.set_bias("drain", vd)
        coarse_vgs = _grid(0.0, vg_hi, STEP_COARSE)

    vgs = list(coarse_vgs)
    ids = _scan(sim, vgs)

    # 관심 구간(1e-10 ~ 1e-7)을 촘촘히 다시 잰다
    lo_c, _ = vg_at_current(vgs, ids, SS_LO)
    hi_c, _ = vg_at_current(vgs, ids, VTH_CURRENT)
    if math.isfinite(lo_c) and math.isfinite(hi_c):
        f_lo = max(min(lo_c, hi_c) - 2 * STEP_COARSE, min(vgs))
        f_hi = min(max(lo_c, hi_c) + 2 * STEP_COARSE, max(vgs))
        fine = _grid(f_lo, f_hi, STEP_FINE)
        # 현재 바이어스에서 가까운 쪽부터 찍어 램프 이동을 줄인다
        if abs(sim.bias["gate"] - fine[-1]) < abs(sim.bias["gate"] - fine[0]):
            fine = fine[::-1]
        vgs += fine
        ids += _scan(sim, fine)

    for vg in extra_points:
        vgs.append(float(vg))
        ids.append(_scan(sim, [vg])[0])
    return vgs, ids


# ---------------------------------------------------------------------------
# 전체 측정
# ---------------------------------------------------------------------------
def _new(device, name, **kw):
    sim = Part1Simulator(device, name=name, **kw)
    sim.build()
    sim.solve_equilibrium()
    sim.enable_transport()
    return sim


def measure(device: Part1Device, quiet: bool = True, curves: dict | None = None,
            **sim_kw) -> dict:
    """8항목 전부를 재서 dict 로 돌려준다. 실패한 항목은 nan 과 error 문자열.

    sim_kw 는 Part1Simulator 로 그대로 넘어간다 (use_ni_temperature 등).
    """
    t0 = time.time()
    res: dict = {"device": device.summary(), "errors": {}, "retries": 0}
    sink = io.StringIO()
    ctx = contextlib.redirect_stdout(sink) if quiet else contextlib.nullcontext()

    with ctx:
        # --- (1) Vd = 0.05 V, 300 K : Vth_lin, SS -------------------------
        try:
            sim = _new(dataclasses.replace(device, temperature_k=T_COLD), "lin", **sim_kw)
            vgs, ids = sweep_for_vth(sim, VD_LIN)
            vth, bw = vg_at_current(vgs, ids, VTH_CURRENT)
            ss, ss_fit, npts = subthreshold_swing(vgs, ids)
            res.update(vth=vth, vth_bracket_V=bw, ss_mV_per_dec=ss,
                       ss_fit_mV_per_dec=ss_fit, ss_points=npts)
            res["retries"] += sim.ramp_retries
            if curves is not None:
                curves["lin_300K"] = (list(vgs), list(ids))
        except Exception as exc:
            res["errors"]["lin"] = f"{type(exc).__name__}: {exc}"
            res.update(vth=float("nan"), ss_mV_per_dec=float("nan"))

        # --- (2) Vd = 2 V, 300 K : Ion, Ioff, Vth_sat -> DIBL -------------
        try:
            sim = _new(dataclasses.replace(device, temperature_k=T_COLD), "sat", **sim_kw)
            vgs, ids = sweep_for_vth(sim, VDD, descending=True, extra_points=(-0.5,))
            byvg = dict(zip(vgs, ids))
            vth_sat, _ = vg_at_current(vgs, ids, VTH_CURRENT)
            res.update(
                ion_uA_per_um=byvg[2.0] * 1e6,
                ioff_pA_per_um=byvg[0.0] * 1e12,
                ioff_neg_pA_per_um=byvg[-0.5] * 1e12,
                vth_sat=vth_sat,
                dibl_mV_per_V=(res["vth"] - vth_sat) / (VDD - VD_LIN) * 1e3,
            )
            res["retries"] += sim.ramp_retries
            if curves is not None:
                curves["sat_300K"] = (list(vgs), list(ids))
        except Exception as exc:
            res["errors"]["sat"] = f"{type(exc).__name__}: {exc}"
            res.update(ion_uA_per_um=float("nan"), ioff_pA_per_um=float("nan"),
                       dibl_mV_per_V=float("nan"))

        # --- (3) body effect : Vd = 0.05 V, Vb = -0.5 V -------------------
        try:
            sim = _new(dataclasses.replace(device, temperature_k=T_COLD), "body", **sim_kw)
            sim.set_bias("body", VB_BODY)       # 드레인보다 먼저 건다
            vgs, ids = sweep_for_vth(sim, VD_LIN)
            vth_b, _ = vg_at_current(vgs, ids, VTH_CURRENT)
            res.update(vth_body=vth_b, body_effect_mV=(vth_b - res["vth"]) * 1e3)
            res["retries"] += sim.ramp_retries
        except Exception as exc:
            res["errors"]["body"] = f"{type(exc).__name__}: {exc}"
            res["body_effect_mV"] = float("nan")

        # --- (4) 고온 Ioff : Vd = 2 V, 398 K ------------------------------
        try:
            sim = _new(dataclasses.replace(device, temperature_k=T_HOT), "hot", **sim_kw)
            sim.set_bias("gate", VDD)
            sim.set_bias("drain", VDD)
            hot = _scan(sim, [0.0, -0.5])
            res.update(ioff_hot_pA_per_um=hot[0] * 1e12,
                       ioff_hot_neg_pA_per_um=hot[1] * 1e12)
            # 두 값이 2배 안이면 접합 누설 바닥에 앉은 것 -> Vth 를 내려도 거의 공짜
            res["hot_on_floor"] = bool(hot[0] > 0 and hot[1] > 0
                                       and hot[0] / hot[1] < 2.0)
            res["retries"] += sim.ramp_retries
        except Exception as exc:
            res["errors"]["hot"] = f"{type(exc).__name__}: {exc}"
            res["ioff_hot_pA_per_um"] = float("nan")

    res["eox_MV_per_cm"] = VDD / (device.oxide_thickness_nm * 1e-7) / 1e6
    res["seconds"] = time.time() - t0
    res["pass"] = verdict(res)
    res["all_pass"] = all(res["pass"].values())
    res["quality"] = quality(res)
    res["trustworthy"] = all(res["quality"].values())
    return res


def quality(res: dict) -> dict:
    """측정 자체를 믿을 수 있는지. 스펙 합격과 따로 본다.

    음수이거나 1e-16 A/um 보다 작은 Ioff 는 "아주 좋은 소자"가 아니라 반올림 잡음이다.
    그걸 기준(1 pA) 이하라고 합격시키면 탐색이 쓰레기 값을 향해 굴러간다.
    """
    q = {}
    bw = res.get("vth_bracket_V", float("nan"))
    q["vth_bracket"] = bool(isinstance(bw, float) and math.isfinite(bw) and bw <= BRACKET_MAX)
    q["ss_points"] = bool(res.get("ss_points", 0) >= SS_MIN_POINTS)
    for key in ("ioff_pA_per_um", "ioff_hot_pA_per_um"):
        v = res.get(key, float("nan"))
        q[key + "_above_floor"] = bool(
            isinstance(v, (int, float)) and math.isfinite(v) and v > IOFF_FLOOR_pA)
    q["no_errors"] = not res.get("errors")
    return q


def verdict(res: dict) -> dict:
    """항목별 합격 여부."""
    out = {}
    for key, (lo, hi) in LIMITS.items():
        v = res.get(key, float("nan"))
        if not isinstance(v, (int, float)) or not math.isfinite(v):
            out[key] = False
            continue
        ok = True
        if lo is not None:
            ok = ok and v >= lo
        if hi is not None:
            ok = ok and v <= hi
        out[key] = bool(ok)
    return out


LABEL = {
    "vth": ("Vth [V]", "0.45 +- 0.05"),
    "ion_uA_per_um": ("Ion [uA/um]", ">= 450"),
    "ioff_pA_per_um": ("Ioff 300K [pA/um]", "<= 1"),
    "ss_mV_per_dec": ("SS [mV/dec]", "<= 75"),
    "ioff_hot_pA_per_um": ("Ioff 398K [pA/um]", "<= 100"),
    "dibl_mV_per_V": ("DIBL [mV/V]", "<= 30"),
    "body_effect_mV": ("Body effect [mV]", "<= 80"),
    "eox_MV_per_cm": ("E_ox [MV/cm]", "<= 5"),
}


def report(res: dict) -> str:
    lines = [f"{'item':22s} {'value':>12s}  {'limit':>14s}  result"]
    for key, (name, lim) in LABEL.items():
        v = res.get(key, float("nan"))
        mark = "PASS" if res["pass"].get(key) else "FAIL"
        lines.append(f"{name:22s} {v:12.4g}  {lim:>14s}  {mark}")
    extra = (f"\nVth bracket {res.get('vth_bracket_V', float('nan')):.4f} V | "
             f"SS points {res.get('ss_points', 0)} | SS fit {res.get('ss_fit_mV_per_dec', float('nan')):.2f} | "
             f"ramp retries {res.get('retries', 0)} | {res.get('seconds', 0):.0f} s")
    bad = [k for k, ok in res.get("quality", {}).items() if not ok]
    if bad:
        extra += "\n측정 신뢰도 경고: " + ", ".join(bad)
    if res.get("hot_on_floor"):
        extra += "\n고온 Ioff 는 접합 누설 바닥 (Vg=0 과 -0.5 V 가 2배 안) -> Vth 를 내려도 거의 공짜"
    if res.get("errors"):
        extra += "\nerrors: " + str(res["errors"])
    return "\n".join(lines) + extra
