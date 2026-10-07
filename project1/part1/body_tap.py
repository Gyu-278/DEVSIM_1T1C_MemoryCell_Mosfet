# -*- coding: utf-8 -*-
"""제출 구조의 body 접점: 소스 왼쪽 윗면의 p+ 탭 위.

HW#1 배포 코드(step2)는 기판 바닥 전체를 body 접점으로 쓴다. 실제 CMOS 는 윗면에서만
만들 수 있어 바닥 접점이 불가능하고, 금속을 옅은 p 실리콘에 바로 대면 쇼트키 장벽이
생길 수 있어 p+ 를 깔아 오믹 접점을 만든다(과제 Background-2). 과제 예시 그림
(Part 1 2쪽, Part 2 4쪽, 변수명 43쪽)도 body 접점을 p+ 탭 위에 그렸다.

설계 탐색(search·refine·ablation, Part 2 의 C~G 라운드)은 배포 구조(바닥 접점)로
했고, 제출 파일과 최종 결과는 이 탭 구조로 만든다. 두 구조의 차이는 접합 누설뿐이다.

구조 (x 는 오른쪽 +, y 는 아래 +):
    x_tap_l ─ p+ 탭(TAP_W) ─ 간격(TAP_GAP, p 기판) ─ x = 0 ─ n+ 소스 ─ 게이트 ─ n+ 드레인
    body 접점 = 탭 윗면, 바닥은 접점 없는 경계(전류가 나가지 않는 면)

간격은 소스 접합 공핍층이 탭 쪽으로 뻗는 거리보다 넓게 잡는다. Part 1 에서 가장
멀리 뻗는 경우는 body −0.5 V 측정이며 약 0.16 µm 다 (NA 7.5e16, 역바이어스 0.5 V).
더 좁으면 공핍층이 탭에 닿아 탭 위치 자체가 결과를 좌우한다.

    python -m part1.body_tap        # 최종 소자 8항목(탭 구조) + 바닥 접점 결과와 비교 + 검사
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import devsim

from mosfet_tool.simulator import UM

from .simulator import Part1Simulator


class TapMixin:
    """탭 기하와 도핑. _build_mesh 는 Part 1 / Part 2 클래스가 각자 쓴다."""

    TAP_W = 0.15 * UM           # p+ 탭 폭
    TAP_GAP = 0.25 * UM         # 탭 오른쪽 끝 ~ n+ 소스 왼쪽 끝(x = 0)
    TAP_DOPING = 1.0e20         # p+ 탭 농도 [cm^-3], 깊이는 S/D 와 같은 xj

    @property
    def x_tap_l(self) -> float:
        return -(self.TAP_GAP + self.TAP_W)

    @property
    def x_tap_r(self) -> float:
        return -self.TAP_GAP

    def _build_doping(self) -> None:
        """배포 도핑 식 + 두 가지.

        1) n+ 소스 식에 x = 0 에서 끊는 erfc 를 곱한다. 배포 식은 왼쪽으로 끝없이
           n+ 인 반평면이라, 영역을 왼쪽으로 늘리면 탭까지 n+ 가 이어진다.
        2) p+ 탭을 같은 erfc 꼴로 빼 준다 (NetDoping 은 n 형 +, p 형 −).
        """
        nd, na = self.dev.sd_doping_cm3, self.dev.body_doping_cm3
        dx = 0.5 * self.DX_CHANNEL
        dy = 0.5 * self.DY_JUNCTION
        yj = self.y_junction
        devsim.node_model(
            device=self.name, region="bulk", name="SourceDoping",
            equation=(f"0.25*{nd:.6e}*erfc((x-{self.x_gate_left:.6e})/{dx:.6e})"
                      f"*erfc((y-{yj:.6e})/{dy:.6e})*0.5*erfc(-x/{dx:.6e})"),
        )
        devsim.node_model(
            device=self.name, region="bulk", name="DrainDoping",
            equation=(f"0.25*{nd:.6e}*erfc(-(x-{self.x_gate_right:.6e})/{dx:.6e})"
                      f"*erfc((y-{yj:.6e})/{dy:.6e})"),
        )
        devsim.node_model(
            device=self.name, region="bulk", name="TapDoping",
            equation=(f"0.25*{self.TAP_DOPING:.6e}*erfc((x-{self.x_tap_r:.6e})/{dx:.6e})"
                      f"*erfc((y-{yj:.6e})/{dy:.6e})"),
        )
        devsim.node_model(
            device=self.name, region="bulk", name="NetDoping",
            equation=f"SourceDoping + DrainDoping - TapDoping - {na:.6e}",
        )


class TapPart1Simulator(TapMixin, Part1Simulator):
    """Part1Simulator 와 같고, 영역을 왼쪽으로 늘려 탭 위에 body 접점을 둔다."""

    def _build_mesh(self) -> None:
        pad = 1.0e-8
        x_min, x_max = self.x_tap_l - pad, self.x_right + pad
        y_min, y_max = self.y_metal_top - pad, self.y_bottom + pad

        devsim.create_2d_mesh(mesh=self.mesh)
        for pos in (x_min, self.x_tap_l, self.x_tap_r, 0.0, self.x_gate_left,
                    self.x_gate_right, self.x_right, x_max):
            devsim.add_2d_mesh_line(mesh=self.mesh, dir="x", pos=pos, ps=self.DX_CHANNEL)
        for pos, spacing in (
            (y_min, self.DY_METAL),
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
        devsim.add_2d_region(mesh=self.mesh, material=self.dev.gate_material,
                             region="gate_metal",
                             xl=self.x_gate_left, xh=self.x_gate_right,
                             yl=self.y_oxide_top, yh=self.y_metal_top)

        devsim.add_2d_contact(mesh=self.mesh, name="gate", region="oxide", material="metal",
                              xl=self.x_gate_left, xh=self.x_gate_right,
                              yl=self.y_oxide_top, yh=self.y_oxide_top)
        devsim.add_2d_contact(mesh=self.mesh, name="source", region="bulk", material="metal",
                              xl=0.0, xh=self.x_gate_left, yl=0.0, yh=0.0)
        devsim.add_2d_contact(mesh=self.mesh, name="drain", region="bulk", material="metal",
                              xl=self.x_gate_right, xh=x_max, yl=0.0, yh=0.0)
        # body: 바닥이 아니라 p+ 탭 윗면
        devsim.add_2d_contact(mesh=self.mesh, name="body", region="bulk", material="metal",
                              xl=x_min, xh=self.x_tap_r, yl=0.0, yh=0.0)
        devsim.add_2d_interface(mesh=self.mesh, name="bulk_oxide",
                                region0="bulk", region1="oxide")
        devsim.finalize_mesh(mesh=self.mesh)
        devsim.create_device(mesh=self.mesh, device=self.name)


def use_tap() -> None:
    """spec.measure 가 만드는 소자를 탭 구조로 바꾼다 (최종 측정·그림용)."""
    from . import spec
    spec.Part1Simulator = TapPart1Simulator


def check_file(path: Path) -> str:
    """제출 검사기(check_structure_file.py)를 돌려 마지막 몇 줄을 돌려준다."""
    root = Path(__file__).resolve().parent.parent
    p = subprocess.run([sys.executable, str(root / "check_structure_file" / "check_structure_file.py"), str(path)],
                       capture_output=True, text=True, cwd=str(root))
    return (p.stdout + p.stderr).strip()


def main() -> int:
    from . import spec
    from .designs import FINAL_DEVICE
    t0 = time.time()
    root = Path(__file__).resolve().parent.parent
    out_dir = root / "results" / "bodytap"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 구조 파일을 따로 써서 검사기에 넣어 본다 (제출 파일 이름은 쓰지 않는다)
    devfile = out_dir / "part1_bodytap.devsim"
    TapPart1Simulator(FINAL_DEVICE, name="tapfile", structure_file=str(devfile)).build()
    checker = check_file(devfile)

    use_tap()
    res = spec.measure(FINAL_DEVICE, quiet=True)
    base = json.loads((root / "results" / "search.json").read_text(encoding="utf-8"))
    out = {"tap_um": dict(width=TapMixin.TAP_W / UM, gap=TapMixin.TAP_GAP / UM),
           "tap_doping_cm3": TapMixin.TAP_DOPING,
           "bodytap": res, "bottom": base.get("final", base),
           "checker": checker, "seconds": time.time() - t0}
    (out_dir / "p1_final_bodytap.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=str), encoding="utf-8")
    keys = ("vth", "ion_uA_per_um", "ioff_pA_per_um", "ss_mV_per_dec", "ioff_hot_pA_per_um",
            "ioff_hot_neg_pA_per_um", "dibl_mV_per_V", "body_effect_mV", "all_pass", "errors")
    for k in keys:
        print(f"{k:>24}: bottom {out['bottom'].get(k)!s:>24}   tap {res.get(k)!s:>24}")
    print(checker[-600:])
    print("seconds", round(out["seconds"], 1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
