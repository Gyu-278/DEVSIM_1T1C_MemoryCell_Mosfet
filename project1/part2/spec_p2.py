# -*- coding: utf-8 -*-
"""[Part 2 / 2026-09-27] 셀 하나를 개정 스펙 전 항목으로 평가한다.

    python -m part2.spec_p2 --cell A

판정 항목
  읽기 '0' / '1'   : WL 상승 0.5 ns 뒤 |V_BL - 1 V| >= 60 mV
  64 ms 읽기       : 64 ms 유지 후 그 셀 전압으로 다시 읽어 >= 60 mV
  리텐션 '0' / '1' : 마진이 60 mV 아래로 내려가는 시각 >= 64 ms
  C_STORE          : <= 20 fF
  E_ox             : 2.5 V / tox <= 5 MV/cm
  E_diel           : 1 V / t_diel <= 4 MV/cm
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import math
import time
from pathlib import Path

from part2 import ranges, transient as T
from part2.cell import CELLS, FINAL, Part2Simulator

OUT = Path(__file__).resolve().parent.parent / "results"
T_RET_REQ = 64.0e-3


def measure(device, quiet: bool = True, dv_step: float = 5.0e-3) -> dict:
    t0 = time.time()
    res: dict = {"errors": {}}
    d = device
    res["device"] = dict(
        Lg_um=d.gate_length_um, tox_nm=d.oxide_thickness_nm, xj_um=d.junction_depth_um,
        tSi_um=d.silicon_thickness_um, Ls_um=d.source_length_um, Ld_um=d.drain_length_um,
        NA_cm3=d.body_doping_cm3, ND_cm3=d.sd_doping_cm3, gate_metal=d.gate_material,
        cap_h_um=d.cap_height_um, cap_t_nm=d.cap_dielectric_nm,
        pillar_nm=d.cap_pillar_nm, cap_material=d.cap_material,
        cap_center_um=d.cap_center_um)
    res["range_violations"] = ranges.check(d, v_wl=T.V_WL)
    sink = io.StringIO()
    ctx = contextlib.redirect_stdout(sink) if quiet else contextlib.nullcontext()

    with ctx:
        # --- C_STORE : 평형에서 plate dQ/dV ------------------------------
        sim0 = Part2Simulator(d, name="cs"); sim0.build(); sim0.solve_equilibrium()
        _, c_fF = sim0.c_store(v0=T.V_PLATE)
        res["c_store_fF"] = c_fF
        res["c_store_closed_fF"] = d.c_store_closed_form_fF
        c_F = c_fF * 1e-15

        # --- 읽기 (fresh) -------------------------------------------------
        for state, v0 in (("0", 0.0), ("1", 2.0)):
            try:
                sim = T.build(d, name=f"rd{state}")
                r = T.read_transient(sim, v0, c_F)
                res[f"read_{state}_mV"] = abs(r["dv_bl"]) * 1e3
                res[f"read_{state}_limit_mV"] = abs(
                    T.charge_share_limit(v0, c_F)) * 1e3
            except Exception as exc:
                res["errors"][f"read{state}"] = f"{type(exc).__name__}: {exc}"
                res[f"read_{state}_mV"] = float("nan")

        # --- 리텐션 : V_crit 을 찾고 거기까지 적분 -------------------------
        for state in ("1", "0"):
            try:
                sim = T.build(d, name=f"rt{state}")
                vc = T.v_crit(sim, state, c_F)
                res[f"vcrit_{state}_V"] = vc
                sim2 = T.build(d, name=f"rh{state}")
                t_ret, (ts, vs, cur), stall = T.retention(
                    sim2, state, c_F, vc, dv_step=dv_step)
                res[f"tret_{state}_ms"] = t_ret * 1e3
                res[f"stall_{state}_V"] = stall
                res[f"ileak_{state}_pA_per_um"] = cur[0] * 1e12 if cur else float("nan")
                # 64 ms 시점의 셀 전압으로 다시 읽기 (과제 정의). 멈춘 셀도 궤적이 64 ms 를
                # 지나 있으므로 그 시점 값을 쓴다. (10-05 이전 판은 정지 전압에서 읽었다.)
                v64 = _v_at(ts, vs, T_RET_REQ)
                res[f"v64_{state}_V"] = v64
                sim3 = T.build(d, name=f"r6{state}")
                res[f"read64_{state}_mV"] = abs(
                    T.read_transient(sim3, v64, c_F)["dv_bl"]) * 1e3
                if stall is not None:            # 셀이 끝내 다가가는 전압에서 읽은 최악값
                    sim4 = T.build(d, name=f"rs{state}")
                    res[f"read_stall_{state}_mV"] = abs(
                        T.read_transient(sim4, stall, c_F)["dv_bl"]) * 1e3
            except Exception as exc:
                res["errors"][f"ret{state}"] = f"{type(exc).__name__}: {exc}"
                res[f"tret_{state}_ms"] = float("nan")
                res[f"read64_{state}_mV"] = float("nan")

    res["eox_MV_per_cm"] = T.V_WL / (d.oxide_thickness_nm * 1e-7) / 1e6
    res["ediel_MV_per_cm"] = d.cap_field_MV_per_cm
    res["seconds"] = time.time() - t0
    res["pass"] = verdict(res)
    res["all_pass"] = all(res["pass"].values()) and not res["range_violations"]
    return res


def _v_at(ts, vs, t_target):
    """궤적에서 t_target 시점의 셀 전압. 끝까지 못 가면 마지막 값."""
    if ts[-1] <= t_target:
        return vs[-1]
    for i in range(1, len(ts)):
        if ts[i] >= t_target:
            f = (t_target - ts[i - 1]) / (ts[i] - ts[i - 1])
            return vs[i - 1] + f * (vs[i] - vs[i - 1])
    return vs[-1]


LIMITS = {
    "read_0_mV": (60.0, None), "read_1_mV": (60.0, None),
    "read64_0_mV": (60.0, None), "read64_1_mV": (60.0, None),
    "tret_0_ms": (64.0, None), "tret_1_ms": (64.0, None),
    "c_store_fF": (None, 20.0),
    "eox_MV_per_cm": (None, 5.0), "ediel_MV_per_cm": (None, 4.0),
}
LABEL = {
    "read_0_mV": ("read '0' [mV]", ">= 60"), "read_1_mV": ("read '1' [mV]", ">= 60"),
    "read64_0_mV": ("64ms read '0' [mV]", ">= 60"), "read64_1_mV": ("64ms read '1' [mV]", ">= 60"),
    "tret_0_ms": ("t_ret '0' [ms]", ">= 64"), "tret_1_ms": ("t_ret '1' [ms]", ">= 64"),
    "c_store_fF": ("C_STORE [fF]", "<= 20"),
    "eox_MV_per_cm": ("E_ox [MV/cm]", "<= 5"), "ediel_MV_per_cm": ("E_diel [MV/cm]", "<= 4"),
}


def verdict(res: dict) -> dict:
    out = {}
    for k, (lo, hi) in LIMITS.items():
        v = res.get(k, float("nan"))
        if not isinstance(v, (int, float)) or (isinstance(v, float) and math.isnan(v)):
            out[k] = False; continue
        ok = True
        if lo is not None: ok = ok and v >= lo
        if hi is not None: ok = ok and v <= hi
        out[k] = bool(ok)
    return out


def report(res: dict) -> str:
    lines = [f"{'item':22s} {'value':>12s}  {'limit':>7s}  result"]
    for k, (n, lim) in LABEL.items():
        v = res.get(k, float("nan"))
        lines.append(f"{n:22s} {v:12.4g}  {lim:>7s}  {'PASS' if res['pass'].get(k) else 'FAIL'}")
    lines.append("")
    lines.append(f"충전분배 한계 '0' {res.get('read_0_limit_mV', float('nan')):.1f} / "
                 f"'1' {res.get('read_1_limit_mV', float('nan')):.1f} mV  |  "
                 f"0.5ns 완료율 '1' {100*res.get('read_1_mV',0)/max(res.get('read_1_limit_mV',1),1e-9):.0f}%")
    lines.append(f"V_crit '1' {res.get('vcrit_1_V', float('nan')):.4f} V  '0' {res.get('vcrit_0_V', float('nan')):.4f} V  |  "
                 f"64ms 셀전압 '1' {res.get('v64_1_V', float('nan')):.4f} V")
    lines.append(f"유지 누설 '1' {res.get('ileak_1_pA_per_um', float('nan')):.4g}  "
                 f"'0' {res.get('ileak_0_pA_per_um', float('nan')):.4g} pA/um  |  {res.get('seconds',0):.0f} s")
    for st in ("0", "1"):
        s = res.get(f"stall_{st}_V")
        if s is not None:
            lines.append(f"  '{st}' 셀은 {s:.4f} V 에서 멈춘다 -> 마진 소진 없음")
    if res.get("range_violations"):
        lines.append("범위 위반: " + "; ".join(res["range_violations"]))
    if res.get("errors"):
        lines.append("errors: " + str(res["errors"]))
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", default=FINAL)
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    d = CELLS[args.cell]
    r = measure(d)
    tag = args.tag or args.cell
    print(f"===== {tag} =====")
    print(report(r))
    OUT.mkdir(exist_ok=True)
    (OUT / f"p2_{tag}.json").write_text(json.dumps(r, indent=2, default=str), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
