# -*- coding: utf-8 -*-
"""[Part 2 / 2026-09-27 개정 스펙] Δt 스텝 과도 READ 엔진.

Slide procedure (p.4 "Read Check"):
  Init: V_BL = 1 V ; V_cell = 0 V ("0") or 2 V ("1") ; WL 2.5 V ; V_B -0.5 V ; 125 C.
  Each dt = 10 ps: DC solve -> I -> dQ = I*dt -> V_cell += dQ/C_STORE, V_BL -= dQ/C_BL.
  Stop at 0.5 ns: |V_BL - 1 V| >= 60 mV -> pass.

Sign convention used here (and proved in the harness):
  devsim get_contact_current is POSITIVE for conventional current flowing INTO the
  device through that contact.  A node that injects current into the device loses
  that charge, hence for BOTH nodes the update is a MINUS:

      V_cell -= I_source * W * dt / C_STORE
      V_BL   -= I_drain  * W * dt / C_BL

  which reproduces the slide's single-I form when I_source = -I_drain.
"""
from __future__ import annotations

import math
import dataclasses
import numpy as np
import devsim
from devsim.python_packages.simple_physics import GetContactBiasName
from mosfet_tool.simulator import UM

from part2.cell import Part2Device, Part2Simulator

# ---- operating point from the 2026-09-27 slide -----------------------------
VDD        = 2.0
V_BL_PRE   = 1.0          # VDD/2 pre-charge
V_PLATE    = 1.0          # plate held at VDD/2
V_WL       = 2.5          # word line high
V_BODY     = -0.5
T_HOT      = 398.0        # 125 C
C_BL_F     = 100.0e-15
V_OFFSET   = 0.060        # 60 mV
T_READ     = 0.5e-9       # decision time
W_UM       = 0.1          # W = 100 nm, currents come out in A/um

NODE_CONTACTS = ("source", "storage_l", "storage_r")
PLATE_CONTACTS = ("plate_l", "plate_r")


# ---------------------------------------------------------------------------
def build(device: Part2Device, t_k: float = T_HOT, name: str = "cell"):
    sim = Part2Simulator(dataclasses.replace(device, temperature_k=t_k), name=name)
    sim.build()
    sim.solve_equilibrium()
    sim.enable_transport()
    return sim


def terminal_currents(sim) -> dict:
    """All four terminal conduction currents [A/um], positive = INTO the device."""
    out = {}
    for c in ("source", "drain", "body"):
        e = devsim.get_contact_current(device=sim.name, contact=c,
                                       equation="ElectronContinuityEquation")
        h = devsim.get_contact_current(device=sim.name, contact=c,
                                       equation="HoleContinuityEquation")
        out[c] = (e + h) * UM
    return out


def _set(sim, contact, v):
    devsim.set_parameter(device=sim.name, name=GetContactBiasName(contact), value=float(v))
    sim.bias[contact] = float(v)


def _solve(sim):
    devsim.solve(type="dc", absolute_error=1.0e30, relative_error=1.0e-6,
                 maximum_iterations=80)


def node_ramp(sim, volts, step=0.025):
    """Ramp source + storage_l + storage_r TOGETHER.

    Driving `source` alone while storage_* stay put puts volts across the 5 nm
    dielectric and the Newton solve freezes (Iteration 0 / RelError 0).  This was
    a real bug in the previous round; the harness re-checks it.
    """
    start = sim.bias[NODE_CONTACTS[0]]
    if abs(volts - start) < 1e-12:
        return
    n = max(1, int(math.ceil(abs(volts - start) / step - 1e-9)))
    for i in range(1, n + 1):
        v = start + (volts - start) * i / n
        for c in NODE_CONTACTS:
            _set(sim, c, v)
        _solve(sim)


def plate_ramp(sim, volts, step=0.1):
    start = sim.bias[PLATE_CONTACTS[0]]
    if abs(volts - start) < 1e-12:
        return
    n = max(1, int(math.ceil(abs(volts - start) / step - 1e-9)))
    for i in range(1, n + 1):
        v = start + (volts - start) * i / n
        for c in PLATE_CONTACTS:
            _set(sim, c, v)
        _solve(sim)


def surface_potential_probe(sim):
    """A scalar that must move whenever the Newton solve actually re-converges."""
    p = devsim.get_node_model_values(device=sim.name, region="bulk", name="Potential")
    return float(np.sum(np.abs(p)))


# ---------------------------------------------------------------------------
def init_read(sim, v_cell0: float, v_bl0: float = V_BL_PRE, wl_on: bool = True):
    """Bring the device to the t = 0- state, then raise WL."""
    sim.set_bias("body", V_BODY)
    plate_ramp(sim, V_PLATE)
    sim.set_bias("drain", v_bl0)      # BL pre-charge
    node_ramp(sim, v_cell0)           # stored level, WL still 0 V
    if wl_on:
        sim.set_bias("gate", V_WL)    # WL rises -> t = 0


def read_transient(sim, v_cell0, c_store_F, dt=10e-12, t_stop=T_READ,
                   v_bl0=V_BL_PRE, use_split_currents=True, probe=False):
    """The Delta-t stepped read.  Returns a dict of traces."""
    init_read(sim, v_cell0, v_bl0)

    t = 0.0
    v_cell, v_bl = v_cell0, v_bl0
    ts, i_s, i_d, i_b, vb, vc, probes = [], [], [], [], [], [], []
    nsteps = int(round(t_stop / dt))
    for k in range(nsteps + 1):
        cur = terminal_currents(sim)
        ts.append(t); vb.append(v_bl); vc.append(v_cell)
        i_s.append(cur["source"]); i_d.append(cur["drain"]); i_b.append(cur["body"])
        if probe:
            probes.append(surface_potential_probe(sim))
        if k == nsteps:
            break
        # charge update: each node loses the charge it injects into the device
        if use_split_currents:
            dq_cell = cur["source"] * W_UM * dt       # C leaving the cell node
            dq_bl   = cur["drain"]  * W_UM * dt       # C leaving the BL node
        else:                                          # slide's single-I form
            i_one = cur["drain"]
            dq_cell = -i_one * W_UM * dt
            dq_bl   =  i_one * W_UM * dt
        v_cell -= dq_cell / c_store_F
        v_bl   -= dq_bl / C_BL_F
        # one DC solve per dt, node contacts driven together
        for c in NODE_CONTACTS:
            _set(sim, c, v_cell)
        _set(sim, "drain", v_bl)
        _solve(sim)
        t += dt

    ts = np.asarray(ts); vb = np.asarray(vb); vc = np.asarray(vc)
    return dict(t=ts, v_bl=vb, v_cell=vc,
                i_source=np.asarray(i_s), i_drain=np.asarray(i_d),
                i_body=np.asarray(i_b), probe=np.asarray(probes),
                dv_bl=vb[-1] - v_bl0, dv_cell=vc[-1] - v_cell0)


def charge_share_limit(v_cell0, c_store_F, v_bl0=V_BL_PRE):
    """dV_BL(inf) for an ideal switch: (V_cell0 - V_BL0) * C_S/(C_S + C_BL)."""
    return (v_cell0 - v_bl0) * c_store_F / (c_store_F + C_BL_F)


# ===========================================================================
# 리텐션 (2026-09-27 개정 스펙)
# ===========================================================================
# 슬라이드 절차:
#   data "1" -> WL 0 V, BL 0 V, cell 2 V ;  data "0" -> WL 0 V, BL 2 V, cell 0 V
#   (V_B -0.5 V, 125 C).  매 Δt 마다 DC 풀어 I_leak 을 얻고
#   V_cell -= I_leak*W*Δt/C_STORE.  셀이 한 스텝에 몇 mV 만 움직이도록 Δt 를 키운다.
#   64 ms 지점에서 그 V_cell 로 0.5 ns 읽기를 다시 돌려 마진이 60 mV 이상인지 본다.
#   t_ret = 그 마진이 60 mV 아래로 내려가는 시각.
#
# 종료 판정이 셀 전압 임계값이 아니라 **읽기 마진**이라는 점이 개정 스펙의 핵심이다.
# 그래서 먼저 "마진이 정확히 60 mV 가 되는 셀 전압" V_crit 을 읽기 엔진으로 찾고,
# 거기까지 가는 시간을 적분한다.

HOLD_BIAS = {           # state -> (V_BL during hold, 초기 셀 전압)
    "1": (0.0, 2.0),
    "0": (2.0, 0.0),
}


def hold_current(sim, state: str, v_cell: float) -> float:
    """유지 상태에서 저장 노드를 떠나는 전류 [A/um].

    저장 노드를 충전·방전하는 것은 **소스 단자** 전류다. 비트라인 단자 전류에는
    body 로 빠지는 접합 누설이 섞여 있고 그건 셀에 오지 않는다. 읽기 때는 두
    전류가 9자리까지 같아 아무거나 써도 되지만, 유지 상태에서는 body 몫이
    0.7% 이고 BL 을 1 V 로 두면 body 쪽이 오히려 더 커진다.
    """
    v_bl, _ = HOLD_BIAS[state]
    sim.set_bias("body", V_BODY)
    sim.set_bias("gate", 0.0)
    plate_ramp(sim, V_PLATE)
    sim.set_bias("drain", v_bl)
    node_ramp(sim, v_cell)
    return terminal_currents(sim)["source"]


def v_crit(sim, state: str, c_store_F: float, tol: float = 5.0e-4,
           dt: float = 10e-12) -> float:
    """읽기 마진이 정확히 V_OFFSET 이 되는 셀 전압을 이분법으로 찾는다.

    '1' 은 2 V 에서 내려오고 '0' 은 0 V 에서 올라오므로 탐색 구간이 다르다.
    """
    v0 = HOLD_BIAS[state][1]
    lo, hi = (V_BL_PRE, v0) if state == "1" else (v0, V_BL_PRE)
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        m = abs(read_transient(sim, mid, c_store_F, dt=dt)["dv_bl"])
        if state == "1":
            if m >= V_OFFSET: hi = mid
            else:             lo = mid
        else:
            if m >= V_OFFSET: lo = mid
            else:             hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def retention(sim, state: str, c_store_F: float, v_target: float,
              dv_step: float = 5.0e-3, t_max: float = 10.0):
    """유지 상태에서 셀 전압이 v_target 에 닿는 시각 [s] 과 궤적.

    적응 스텝: 셀이 한 번에 dv_step 만큼만 움직이도록 Δt 를 매번 다시 잡는다.
    마지막 스텝은 v_target 에서 정확히 멈추도록 Δt 를 비례해서 줄인다.
    (줄이지 않고 전체 Δt 를 더하면 일어나지도 않은 이동을 시간에 넣게 되어
     t_ret 이 과대평가된다. 작은 누설 설계에서 최대 19% 까지 틀린다.)

    전류 부호가 도중에 뒤집히면 셀이 그 지점에서 멈춘다 -> (inf, 궤적, 정지전압).
    t_max 안에 v_target 에 못 닿아도 (inf, 궤적, 마지막 전압) 이다. 정지점 바로 아래에서는
    전류가 0 에 가까워 Δt 가 끝없이 커지기 때문이다. (이전 판은 이때 t_max 를 넘긴
    시각을 t_ret 처럼 돌려줬다. F5·G1 의 '0' 이 약 10 s 로 나온 것이 그 경우다.)
    """
    v_bl, v0 = HOLD_BIAS[state]
    sgn = 1.0 if v_target > v0 else -1.0
    v, t = v0, 0.0
    vs, ts, cur = [v0], [0.0], []
    while (v_target - v) * sgn > 1e-12 and t < t_max:
        i_leave = hold_current(sim, state, v)      # A/um, 노드를 떠나면 양수
        dv_dt = -i_leave * W_UM / c_store_F        # V/s
        cur.append(i_leave)
        if dv_dt * sgn <= 0:                       # 더 이상 목표 쪽으로 안 간다
            return float("inf"), (ts, vs, cur), v
        dv = sgn * dv_step
        if (v + dv - v_target) * sgn >= 0.0:       # 마지막 스텝
            dv = v_target - v
        t += dv / dv_dt
        v += dv
        ts.append(t); vs.append(v)
    if (v_target - v) * sgn > 1e-12:               # t_max 안에 v_target 에 못 닿음
        return float("inf"), (ts, vs, cur), v
    return t, (ts, vs, cur), None
