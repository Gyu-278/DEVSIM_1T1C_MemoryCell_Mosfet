# -*- coding: utf-8 -*-
"""설계 탐색. 후보를 차례로 재서 탐색 기록표(보고서 ②)를 만든다.

    python -m part1.search                 # designs.py 의 후보 전부
    python -m part1.search --only A B      # 일부만
    python -m part1.search --metals        # 최종 구조에서 게이트 금속만 훑기

게이트 금속만 바꿀 때는 8항목을 다시 재지 않는다. 일함수 오프셋은 게이트 전압축을
평행이동시키는 항등식이라 SS · DIBL · body effect 는 전혀 변하지 않고,
고정된 절대 Vg 에서 읽는 Ion · Ioff · 고온 Ioff 만 변한다.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path

from . import csvout, spec
from .designs import DESIGNS
from .gate_physics import GATE_METALS

OUT = Path(__file__).resolve().parent.parent / "results"

COLS = ["vth", "ion_uA_per_um", "ioff_pA_per_um", "ss_mV_per_dec",
        "ioff_hot_pA_per_um", "dibl_mV_per_V", "body_effect_mV", "eox_MV_per_cm"]
HEAD = ["Vth", "Ion", "Ioff300", "SS", "Ioff398", "DIBL", "Body", "Eox"]


def table(rows: list[tuple[str, dict]]) -> str:
    out = [f"{'design':10s} " + " ".join(f"{h:>9s}" for h in HEAD) + "   verdict"]
    for name, r in rows:
        cells = []
        for c in COLS:
            v = r.get(c, float("nan"))
            cells.append(f"{v:9.3g}" if r["pass"].get(c) else f"{v:8.3g}*")
        out.append(f"{name:10s} " + " ".join(cells) +
                   ("   ALL PASS" if r["all_pass"] else "   fail: " +
                    ",".join(k for k, ok in r["pass"].items() if not ok)))
    out.append("(* = 기준 미달)")
    return "\n".join(out)


def run(keys, metals=False, quiet=True):
    rows = []
    for key in keys:
        dev = DESIGNS[key]
        res = spec.measure(dev, quiet=quiet)
        rows.append((key, res))
        print(f"\n===== {key} =====")
        print(spec.report(res), flush=True)

        if metals:
            base_ewf = dev.gate_work_function_ev
            print(f"\n  -- gate metal sweep (구조 고정, Vth 는 1 V/eV 로 이동) --")
            print(f"  {'metal':8s} {'EWF':>5s} {'Vth pred':>9s} {'Ion':>9s} {'Ioff300':>9s} {'Ioff398':>9s}")
            for m, ewf in sorted(GATE_METALS.items(), key=lambda kv: kv[1]):
                pred = res["vth"] + (ewf - base_ewf)
                if not (0.35 <= pred <= 0.55):
                    continue          # 스펙 창 근처만 실제로 재 본다
                r2 = spec.measure(dataclasses.replace(dev, gate_material=m), quiet=quiet)
                print(f"  {m:8s} {ewf:5.2f} {pred:9.4f} {r2['ion_uA_per_um']:9.1f} "
                      f"{r2['ioff_pA_per_um']:9.3g} {r2['ioff_hot_pA_per_um']:9.3g}"
                      f"  measured Vth {r2['vth']:.4f}"
                      + ("  ALL PASS" if r2["all_pass"] else ""))
                rows.append((f"{key}/{m}", r2))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--metals", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    keys = args.only or list(DESIGNS)
    rows = run(keys, metals=args.metals, quiet=not args.verbose)

    print("\n" + "=" * 100)
    print(table(rows))
    OUT.mkdir(exist_ok=True)
    (OUT / "search.json").write_text(
        json.dumps({n: r for n, r in rows}, indent=2, default=str), encoding="utf-8")
    (OUT / "search_table.txt").write_text(table(rows), encoding="utf-8")
    csv_path = csvout.write_specs(rows, OUT / "search.csv")
    print(f"\nsaved: {OUT/'search.json'} , {OUT/'search_table.txt'} , {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
