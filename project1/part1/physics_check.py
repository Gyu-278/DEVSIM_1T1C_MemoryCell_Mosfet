# -*- coding: utf-8 -*-
"""보고서 ① 용: 추가한 물리가 실제로 무엇을 바꾸는지 전후 비교표.

    python -m part1.physics_check --design baseline

배포 코드 상태(두 물리 모두 off)에서 시작해 하나씩 켜 가며 8항목을 다시 잰다.
"추가한 물리가 없으면 이 스펙은 측정 자체가 안 된다"는 주장을 숫자로 보인다.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import csvout, spec
from .designs import DESIGNS

OUT = Path(__file__).resolve().parent.parent / "results"

# 추가한 물리 세 가지를 하나씩 켠다. 기준 행은 셋 다 꺼서 배포 코드와 같게 만든다.
_OFF = dict(use_work_function=False, use_ni_temperature=False, use_mobility_t=False)
CASES = [
    ("배포 코드 그대로 (midgap · n_i 동결 · mu 상수)", dict(_OFF)),
    ("+ 게이트 일함수만",       dict(_OFF, use_work_function=True)),
    ("+ n_i(T) 만",             dict(_OFF, use_ni_temperature=True)),
    ("+ mu(T) 만",              dict(_OFF, use_mobility_t=True)),
    ("셋 다 (최종)",            dict(use_work_function=True, use_ni_temperature=True,
                                     use_mobility_t=True)),
]

KEYS = [("vth", "Vth [V]"), ("ss_mV_per_dec", "SS [mV/dec]"),
        ("ioff_pA_per_um", "Ioff300 [pA]"), ("ioff_hot_pA_per_um", "Ioff398 [pA]"),
        ("ion_uA_per_um", "Ion [uA]"), ("dibl_mV_per_V", "DIBL [mV/V]"),
        ("body_effect_mV", "Body [mV]")]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default="baseline", choices=sorted(DESIGNS))
    args = ap.parse_args()
    dev = DESIGNS[args.design]

    rows = []
    for label, kw in CASES:
        r = spec.measure(dev, **kw)
        rows.append((label, r))
        print(f"[done] {label}", flush=True)

    head = f"{'case':42s} " + " ".join(f"{n:>13s}" for _, n in KEYS)
    lines = [f"소자: {args.design}", head, "-" * len(head)]
    for label, r in rows:
        lines.append(f"{label:42s} " +
                     " ".join(f"{r.get(k, float('nan')):13.4g}" for k, _ in KEYS))
    base = rows[0][1]
    full = rows[-1][1]
    lines += ["", "핵심 변화:",
              f"  고온 Ioff 비율 (398K/300K):  물리 없음 {base.get('ioff_hot_pA_per_um',float('nan'))/max(base.get('ioff_pA_per_um',1),1e-30):9.1f}x"
              f"   ->  추가 후 {full.get('ioff_hot_pA_per_um',float('nan'))/max(full.get('ioff_pA_per_um',1),1e-30):9.1f}x",
              f"  Vth 온도계수는 physics_check 의 별도 항목 참조 (n_i 동결 시 부호가 거꾸로)"]
    text = "\n".join(lines)
    print("\n" + text)
    OUT.mkdir(exist_ok=True)
    (OUT / f"physics_check_{args.design}.txt").write_text(text, encoding="utf-8")
    (OUT / f"physics_check_{args.design}.json").write_text(
        json.dumps({l: r for l, r in rows}, indent=2, default=str), encoding="utf-8")
    csv_path = csvout.write_specs(rows, OUT / f"physics_check_{args.design}.csv")
    print(f"\nsaved: {OUT/f'physics_check_{args.design}.txt'} , {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
