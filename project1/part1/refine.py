# -*- coding: utf-8 -*-
"""정밀 탐색. 변형 하나를 별 프로세스에서 재서 결과 JSON 을 남긴다.

DEVSIM 은 한 프로세스에서 순차로만 풀리므로, 변형마다 프로세스를 띄워 병렬화한다.

    python -m part1.refine --spec "oxide_thickness_nm=5.0,body_doping_cm3=6e16" --tag A_tox5_na6
"""

from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path

from . import spec
from .designs import DESIGNS
from .gate_physics import GATE_METALS

OUT = Path(__file__).resolve().parent.parent / "results" / "refine"

FIELDS = {f.name: f.type for f in dataclasses.fields(DESIGNS["A"])}


def parse(text: str) -> dict:
    out = {}
    for part in filter(None, (p.strip() for p in text.split(","))):
        k, v = part.split("=", 1)
        k = k.strip()
        if k not in FIELDS:
            raise SystemExit(f"unknown field {k}; known: {sorted(FIELDS)}")
        out[k] = v.strip() if k == "gate_material" else float(v)
    return out


def best_metal(vth_measured: float, current_metal: str) -> tuple[str, float]:
    """구조를 그대로 두고 Vth 를 0.45 V 에 가장 가깝게 만드는 금속.

    일함수 오프셋은 게이트 전압축 평행이동이라 Vth 가 정확히 1 V/eV 로 움직인다.
    그래서 다시 시뮬레이션하지 않고 표에서 고를 수 있다.
    """
    base = GATE_METALS[current_metal]
    pick = min(GATE_METALS.items(),
               key=lambda kv: abs(vth_measured + (kv[1] - base) - 0.45))
    return pick[0], vth_measured + (pick[1] - base)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="A", choices=sorted(DESIGNS))
    ap.add_argument("--spec", default="")
    ap.add_argument("--tag", required=True)
    args = ap.parse_args()

    dev = dataclasses.replace(DESIGNS[args.base], **parse(args.spec))
    res = spec.measure(dev)
    res["tag"] = args.tag
    metal, vth_pred = best_metal(res.get("vth", float("nan")), dev.gate_material)
    res["best_metal"] = metal
    res["vth_with_best_metal"] = vth_pred

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{args.tag}.json").write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    print(f"===== {args.tag} =====")
    print(f"spec: {args.spec or '(base)'}   gate {dev.gate_material}")
    print(spec.report(res))
    print(f"best metal for Vth->0.45 : {metal}  (Vth would be {vth_pred:.4f} V)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
