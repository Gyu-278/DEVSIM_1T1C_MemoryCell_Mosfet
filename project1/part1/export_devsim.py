# -*- coding: utf-8 -*-
"""제출용 구조 파일 part1_<학번>.devsim 을 쓰고, 공식 검사기를 바로 돌린다.

    python -m part1.export_devsim            # 최종 설계로 저장 + 검사
    python -m part1.export_devsim --name 임시  # 다른 이름으로

제출 구조는 body 접점을 소스 왼쪽 p+ 탭 위에 둔다(part1/body_tap.py).
저장 시점은 Part1Simulator.build() 안에서 처리된다. 메시·영역·접점·계면·도핑까지
만든 직후, 방정식을 등록하기 전이다. 이 순서가 아니면 검사기가 REJECT 한다.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import subprocess
import sys
from pathlib import Path

from .designs import DESIGNS, FINAL
from .body_tap import TapPart1Simulator

STUDENT_ID = "2021144048"
ROOT = Path(__file__).resolve().parent.parent          # project1/
CHECKER = ROOT / "check_structure_file" / "check_structure_file.py"


def write(design_key: str = FINAL, out_path: Path | None = None,
          device_name: str = "cell", quiet: bool = True) -> Path:
    device = DESIGNS[design_key]
    out = Path(out_path) if out_path else ROOT / f"part1_{STUDENT_ID}.devsim"
    sink = io.StringIO()
    ctx = contextlib.redirect_stdout(sink) if quiet else contextlib.nullcontext()
    with ctx:
        sim = TapPart1Simulator(device, name=device_name, structure_file=str(out))
        sim.build()          # 도핑 직후 저장 -> 그 다음 물리 등록
        sim.solve_equilibrium()
        sim.enable_transport()   # 저장 후에도 계속 돌아가는지 함께 확인
    return out


def check(path: Path) -> tuple[int, str]:
    proc = subprocess.run([sys.executable, str(CHECKER), str(path)],
                          capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default=FINAL, choices=sorted(DESIGNS))
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    out = write(args.design, args.out)
    print(f"design   : {args.design}")
    print(f"written  : {out}  ({out.stat().st_size} bytes)")
    code, text = check(out)
    tail = [l for l in text.splitlines()
            if l.startswith(("OK:", "REJECT:", "  ")) or "Error" in l or "error" in l]
    print("\n".join(tail))
    print("checker exit:", code, "(0 = OK)")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
