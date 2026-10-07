# -*- coding: utf-8 -*-
"""탐색 결과를 한 표로 모은다 (보고서 ② 탐색 기록표).

    python -m part1.collect
"""

from __future__ import annotations

import json
from pathlib import Path

from . import csvout

OUT = Path(__file__).resolve().parent.parent / "results"

COLS = [("vth", "Vth", "0.45±0.05", 0.45, 0.05),
        ("ion_uA_per_um", "Ion", ">=450", None, None),
        ("ioff_pA_per_um", "Ioff300", "<=1", None, None),
        ("ss_mV_per_dec", "SS", "<=75", None, None),
        ("ioff_hot_pA_per_um", "Ioff398", "<=100", None, None),
        ("dibl_mV_per_V", "DIBL", "<=30", None, None),
        ("body_effect_mV", "Body", "<=80", None, None),
        ("eox_MV_per_cm", "Eox", "<=5", None, None)]

# 허용 예산을 얼마나 썼는지 (0 = 여유 만점, 1 = 기준선, >1 = 미달).
# 항목마다 단위가 달라서 비율을 그대로 비교하면 안 된다. 전부 "예산 소진율"로 환산한다.
MARGIN = {"ion_uA_per_um": ("lo", 450.0), "ioff_pA_per_um": ("hi", 1.0),
          "ss_mV_per_dec": ("hi", 75.0), "ioff_hot_pA_per_um": ("hi", 100.0),
          "dibl_mV_per_V": ("hi", 30.0), "body_effect_mV": ("hi", 80.0),
          "eox_MV_per_cm": ("hi", 5.0)}


def margins(r: dict) -> dict:
    """항목별 예산 소진율. 작을수록 안전하다."""
    out = {}
    for k, (kind, lim) in MARGIN.items():
        v = r.get(k)
        if not isinstance(v, (int, float)) or v != v or v <= 0:
            out[k] = float("nan"); continue
        out[k] = (lim / v) if kind == "lo" else (v / lim)
    vth = r.get("vth", float("nan"))
    out["vth"] = abs(vth - 0.45) / 0.05 if vth == vth else float("nan")
    return out


def main() -> int:
    rows = []
    for p in sorted((OUT / "refine").glob("*.json")):
        rows.append(json.loads(p.read_text(encoding="utf-8")))
    sj = OUT / "search.json"
    if sj.exists():
        for name, r in json.loads(sj.read_text(encoding="utf-8")).items():
            r["tag"] = name
            rows.append(r)

    head = f"{'variant':16s} " + " ".join(f"{c[1]:>8s}" for c in COLS) + "  metal  verdict"
    lines = [head, "-" * len(head)]
    scored = []
    for r in rows:
        cells = []
        for key, _, _, _, _ in COLS:
            v = r.get(key, float("nan"))
            ok = r.get("pass", {}).get(key, False)
            cells.append(f"{v:8.3g}" if ok else f"{v:7.3g}*")
        m = margins(r)
        worst = max((v for v in m.values() if v == v), default=float("nan"))
        allp = r.get("all_pass", False)
        trust = r.get("trustworthy", True)
        verdict = ("ALL PASS" if allp else "fail:" + ",".join(
            k for k, ok in r.get("pass", {}).items() if not ok))
        if not trust:
            verdict += " (측정 경고)"
        lines.append(f"{r.get('tag','?'):16s} " + " ".join(cells) +
                     f"  {r.get('device',{}).get('gate_material','?'):5s}  {verdict}")
        if allp and trust:
            scored.append((worst, r))

    lines.append("(* = 기준 미달)")
    lines.append("")
    lines.append("전 항목 통과 후보를 '가장 빡빡한 항목의 예산 소진율'이 낮은 순서로")
    lines.append("(0 = 여유 만점, 1.00 = 기준선에 딱 걸침):")
    lines.append(f"  {'variant':16s} {'최악':>6s}  {'그 항목':<20s} 두 번째로 빡빡한 항목")
    for worst, r in sorted(scored, key=lambda t: t[0]):
        m = {k: v for k, v in margins(r).items() if v == v}
        rank = sorted(m.items(), key=lambda kv: -kv[1])
        second = f"{rank[1][0]} {rank[1][1]:.2f}" if len(rank) > 1 else ""
        lines.append(f"  {r.get('tag','?'):16s} {worst:6.2f}  {rank[0][0]:<20s} {second}")

    text = "\n".join(lines)
    print(text)
    (OUT / "refine_table.txt").write_text(text, encoding="utf-8")
    csv_path = csvout.write_specs([(r.get("tag", "?"), r) for r in rows],
                                  OUT / "search_record.csv")
    print(f"\nsaved: {OUT/'refine_table.txt'} , {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
