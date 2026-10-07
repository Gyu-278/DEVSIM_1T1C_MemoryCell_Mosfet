# -*- coding: utf-8 -*-
"""Part 2 제출 셀(G2)의 탭 구조와 최종 측정.

body 접점을 소스 왼쪽 윗면의 p+ 탭 위에 둔다(이유는 part1.body_tap 참고). 탭 기하·도핑은
part1.body_tap.TapMixin 과 같고 간격만 넓다: '1' 유지 상태(저장 노드 2 V, body −0.5 V,
125 °C)에서 소스 접합 공핍층이 p 기판(9.5e16) 쪽으로 약 0.21 µm 뻗으므로 0.30 µm 로 둔다.
C~G 라운드 탐색은 배포 구조(바닥 접점)로 했다.

    python -m part2.body_tap        # G2 를 탭 구조로 측정 + 바닥 접점 결과와 비교 + 검사
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import time
from pathlib import Path

import devsim

from mosfet_tool.simulator import UM
from part1.body_tap import TapMixin, check_file

from .cell import CELLS, FINAL, Part2Simulator


class TapPart2Simulator(TapMixin, Part2Simulator):
    """Part2Simulator 와 같고, 영역을 왼쪽으로 늘려 탭 위에 body 접점을 둔다."""

    TAP_GAP = 0.30 * UM

    def _build_mesh(self) -> None:
        pad = 1.0e-8
        d = self.dev
        x_min, x_max = self.x_tap_l - pad, self.x_right + pad
        y_top = min(self.y_cap_top, self.y_metal_top)
        y_min, y_max = y_top - pad, self.y_bottom + pad

        devsim.create_2d_mesh(mesh=self.mesh)
        for pos, ns, ps in (
            (x_min, self.DX_CHANNEL, self.DX_CHANNEL),
            (self.x_tap_l, self.DX_CHANNEL, self.DX_CHANNEL),
            (self.x_tap_r, self.DX_CHANNEL, self.DX_CHANNEL),
            (0.0, self.DX_CHANNEL, self.DX_CHANNEL),
            (self.x_src_contact, self.DX_CHANNEL, self.DX_CHANNEL),
            (self.x_hk_l_lo, self.DX_CHANNEL, self.DX_CAP),
            (self.x_hk_l_hi, self.DX_CAP, self.DX_CHANNEL),
            (self.x_hk_r_lo, self.DX_CHANNEL, self.DX_CAP),
            (self.x_hk_r_hi, self.DX_CAP, self.DX_CHANNEL),
            (self.x_gate_left, self.DX_CHANNEL, self.DX_CHANNEL),
            (self.x_gate_right, self.DX_CHANNEL, self.DX_CHANNEL),
            (self.x_right, self.DX_CHANNEL, self.DX_CHANNEL),
            (x_max, self.DX_CHANNEL, self.DX_CHANNEL),
        ):
            devsim.add_2d_mesh_line(mesh=self.mesh, dir="x", pos=pos, ns=ns, ps=ps)
        for pos, spacing in (
            (y_min, self.DY_CAP),
            (self.y_cap_top, self.DY_CAP),
            (self.y_metal_top, self.DY_METAL),
            (self.y_oxide_top, self.DY_OXIDE),
            (0.0, min(self.DY_OXIDE, self.DY_JUNCTION)),
            (self.y_junction, self.DY_JUNCTION),
            (self.y_bottom, self.DY_BULK),
            (y_max, self.DY_BULK),
        ):
            devsim.add_2d_mesh_line(mesh=self.mesh, dir="y", pos=pos, ps=spacing)

        devsim.add_2d_region(mesh=self.mesh, material="Air", region="air")
        devsim.add_2d_region(mesh=self.mesh, material="Silicon", region="bulk",
                             xl=self.x_tap_l, xh=self.x_right, yl=self.y_bottom, yh=0.0)
        devsim.add_2d_region(mesh=self.mesh, material="Oxide", region="oxide",
                             xl=self.x_gate_left, xh=self.x_gate_right,
                             yl=0.0, yh=self.y_oxide_top)
        devsim.add_2d_region(mesh=self.mesh, material=d.gate_material, region="gate_metal",
                             xl=self.x_gate_left, xh=self.x_gate_right,
                             yl=self.y_oxide_top, yh=self.y_metal_top)
        devsim.add_2d_region(mesh=self.mesh, material=d.cap_material, region="hk_l",
                             xl=self.x_hk_l_lo, xh=self.x_hk_l_hi,
                             yl=self.y_cap_top, yh=0.0)
        devsim.add_2d_region(mesh=self.mesh, material=d.cap_material, region="hk_r",
                             xl=self.x_hk_r_lo, xh=self.x_hk_r_hi,
                             yl=self.y_cap_top, yh=0.0)
        devsim.add_2d_region(mesh=self.mesh, material="TiN", region="metal_pillar",
                             xl=self.x_hk_l_hi, xh=self.x_hk_r_lo,
                             yl=self.y_cap_top, yh=0.0)

        devsim.add_2d_contact(mesh=self.mesh, name="gate", region="oxide", material="metal",
                              xl=self.x_gate_left, xh=self.x_gate_right,
                              yl=self.y_oxide_top, yh=self.y_oxide_top)
        devsim.add_2d_contact(mesh=self.mesh, name="source", region="bulk", material="metal",
                              xl=0.0, xh=self.x_src_contact, yl=0.0, yh=0.0)
        devsim.add_2d_contact(mesh=self.mesh, name="drain", region="bulk", material="metal",
                              xl=self.x_gate_right, xh=x_max, yl=0.0, yh=0.0)
        # body: 바닥이 아니라 p+ 탭 윗면
        devsim.add_2d_contact(mesh=self.mesh, name="body", region="bulk", material="metal",
                              xl=x_min, xh=self.x_tap_r, yl=0.0, yh=0.0)
        devsim.add_2d_contact(mesh=self.mesh, name="plate_l", region="hk_l", material="metal",
                              xl=self.x_hk_l_lo, xh=self.x_hk_l_lo,
                              yl=self.y_cap_top, yh=0.0)
        devsim.add_2d_contact(mesh=self.mesh, name="storage_l", region="hk_l", material="metal",
                              xl=self.x_hk_l_hi, xh=self.x_hk_l_hi,
                              yl=self.y_cap_top, yh=0.0)
        devsim.add_2d_contact(mesh=self.mesh, name="storage_r", region="hk_r", material="metal",
                              xl=self.x_hk_r_lo, xh=self.x_hk_r_lo,
                              yl=self.y_cap_top, yh=0.0)
        devsim.add_2d_contact(mesh=self.mesh, name="plate_r", region="hk_r", material="metal",
                              xl=self.x_hk_r_hi, xh=self.x_hk_r_hi,
                              yl=self.y_cap_top, yh=0.0)
        devsim.add_2d_interface(mesh=self.mesh, name="bulk_oxide",
                                region0="bulk", region1="oxide")
        devsim.add_2d_interface(mesh=self.mesh, name="bulk_hk_l",
                                region0="bulk", region1="hk_l")
        devsim.add_2d_interface(mesh=self.mesh, name="bulk_hk_r",
                                region0="bulk", region1="hk_r")
        devsim.finalize_mesh(mesh=self.mesh)
        devsim.create_device(mesh=self.mesh, device=self.name)


def use_tap() -> None:
    """spec_p2 · transient · figures2 가 만드는 셀을 탭 구조로 바꾼다 (최종 측정·그림용)."""
    from . import figures2, spec_p2
    from . import transient as T
    spec_p2.Part2Simulator = TapPart2Simulator
    T.Part2Simulator = TapPart2Simulator
    figures2.Part2Simulator = TapPart2Simulator


def hold_split(sim_cls, device, v_cells, name: str) -> dict:
    """'1' 유지 바이어스에서 단자별 전류 [pA/um]. 저장 노드(source) 누설이 BL 과 body
    로 어떻게 나뉘는지 본다."""
    from . import transient as T
    T.Part2Simulator = sim_cls
    out = {}
    with contextlib.redirect_stdout(io.StringIO()):
        sim = T.build(device, name=name)
        for v in v_cells:
            T.hold_current(sim, "1", v)
            out[f"{v:.4f}"] = {k: c * 1e12 for k, c in T.terminal_currents(sim).items()}
    return out


def main() -> int:
    from . import spec_p2
    from . import transient as T
    t0 = time.time()
    root = Path(__file__).resolve().parent.parent
    out_dir = root / "results" / "bodytap"
    out_dir.mkdir(parents=True, exist_ok=True)
    g2 = CELLS[FINAL]

    devfile = out_dir / "part2_bodytap.devsim"
    with contextlib.redirect_stdout(io.StringIO()):
        TapPart2Simulator(g2, name="tapfile2", structure_file=str(devfile)).build()
    checker = check_file(devfile)
    print(checker[-500:], flush=True)

    # spec_p2.measure 가 만드는 모든 소자를 탭 구조로, 리텐션 궤적은 가로채 저장
    use_tap()
    traj: dict = {}
    orig_retention = T.retention

    def capture(sim, state, c_store_F, v_target, **kw):
        r = orig_retention(sim, state, c_store_F, v_target, **kw)
        t_ret, (ts, vs, cur), stall = r
        traj[state] = dict(t=list(ts), v=list(vs), i=list(cur), stall=stall,
                           t_ret=t_ret, vcrit=v_target)
        return r

    T.retention = capture
    res = spec_p2.measure(g2)
    T.retention = orig_retention
    print("measure done", round(time.time() - t0), "s", flush=True)

    base = json.loads((root / "results" / "p2_G2.json").read_text(encoding="utf-8"))
    out = {"tap_um": dict(width=TapPart2Simulator.TAP_W / UM, gap=TapPart2Simulator.TAP_GAP / UM),
           "tap_doping_cm3": TapPart2Simulator.TAP_DOPING,
           "bodytap": res, "bottom": base, "retention_bodytap": traj, "checker": checker}
    try:
        out["split_bodytap"] = hold_split(TapPart2Simulator, g2, (2.0, res["vcrit_1_V"]), "spt")
        out["split_bottom"] = hold_split(Part2Simulator, g2, (2.0, base["vcrit_1_V"]), "spb")
    except Exception as exc:                          # 보조 측정이라 실패해도 본 결과는 남긴다
        out["split_error"] = f"{type(exc).__name__}: {exc}"
    out["seconds"] = time.time() - t0
    (out_dir / "p2_G2_bodytap.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=str), encoding="utf-8")
    keys = ("c_store_fF", "read_0_mV", "read_1_mV", "vcrit_1_V", "vcrit_0_V", "tret_1_ms",
            "tret_0_ms", "stall_0_V", "v64_1_V", "read64_1_mV", "read64_0_mV",
            "ileak_1_pA_per_um", "all_pass", "errors")
    for k in keys:
        print(f"{k:>18}: bottom {base.get(k)!s:>22}   tap {res.get(k)!s:>22}")
    for k in ("split_bottom", "split_bodytap"):
        print(k, out.get(k))
    print("seconds", round(out["seconds"], 1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
