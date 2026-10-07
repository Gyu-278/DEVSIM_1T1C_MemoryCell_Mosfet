# -*- coding: utf-8 -*-
"""제출용 part2_<학번>.devsim 저장 + 공식 검사기 실행.

    python -m part2.export_devsim2

제출 구조는 body 접점을 소스 왼쪽 p+ 탭 위에 둔다(part2/body_tap.py).
저장 시점은 Part1Simulator.build() 가 처리한다. 메시·영역·접점·계면·도핑까지
만든 직후, 방정식과 Potential 을 만들기 전이다. Part 2 는 이 규칙이 더 엄격하다.
hk 영역에 Potential 을 만들기만 해도 검사기가 REJECT 한다.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import subprocess
import sys
from pathlib import Path

import devsim

from part2.body_tap import TapPart2Simulator
from part2.cell import CELLS, FINAL
from part2 import ranges

STUDENT_ID = "2021144048"
ROOT = Path(__file__).resolve().parent.parent
CHECKER = ROOT / "check_structure_file" / "check_structure_file.py"


def write(design_key: str = FINAL, out_path: Path | None = None,
          device_name: str = "cell", quiet: bool = True) -> Path:
    device = CELLS[design_key]
    ranges.assert_ok(device)            # 조교가 역산해서 보는 범위를 먼저 막는다
    out = Path(out_path) if out_path else ROOT / f"part2_{STUDENT_ID}.devsim"
    sink = io.StringIO()
    ctx = contextlib.redirect_stdout(sink) if quiet else contextlib.nullcontext()
    with ctx:
        sim = TapPart2Simulator(device, name=device_name, structure_file=str(out))
        sim.build()
        # 검사기는 get_device_list()[0] 만 본다. 소자가 둘 이상 살아 있으면
        # 엉뚱한 쪽이 채점될 수 있으므로 하나뿐임을 확인한다.
        assert len(devsim.get_device_list()) == 1, devsim.get_device_list()
        sim.solve_equilibrium()          # 저장 후에도 계속 풀리는지 확인
    return out


def check(path: Path) -> tuple[int, str]:
    proc = subprocess.run([sys.executable, str(CHECKER), str(path)],
                          capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default=FINAL, choices=sorted(CELLS))
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    out = write(args.design, args.out)
    d = CELLS[args.design]
    print(f"design   : {args.design}")
    print(f"  TRPASS : Lg {d.gate_length_um} um, tox {d.oxide_thickness_nm} nm, "
          f"xj {d.junction_depth_um} um, NA {d.body_doping_cm3:.1e}, "
          f"ND {d.sd_doping_cm3:.1e}, gate {d.gate_material}")
    print(f"  CAP    : h {d.cap_height_um} um, {d.cap_material} {d.cap_dielectric_nm} nm, "
          f"pillar {d.cap_pillar_nm} nm  -> 해석식 {d.c_store_closed_form_fF:.3f} fF")
    print(f"written  : {out}  ({out.stat().st_size} bytes)")
    code, text = check(out)
    for line in text.splitlines():
        if line.startswith(("OK:", "REJECT:", "  ")) or "rror" in line:
            print(line)
    print("checker exit:", code, "(0 = OK)")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
