# -*- coding: utf-8 -*-
"""Gate oxide thickness(tox)에 따른 Cgg-Vg 특성을 비교한다.

config.yaml의 device: 항목을 기본 소자로 사용하고,
compare: 아래에 지정된 변경값만 소자별로 덮어쓴다.
각 소자에 대해 Step 2에서 구현한 run_cv()를 실행한 뒤
Cgg-Vg 곡선을 하나의 그래프에 겹쳐 표시한다.

실행: STEP2_COMPARE_CV.bat 더블클릭
"""

from pathlib import Path

import pandas as pd
import plotly.express as px
import yaml

from mosfet_tool.config import Device
from mosfet_tool.workflows import run_cv


CONFIG_PATH = Path(__file__).resolve().parent / "config.yaml"


def main():
    # ------------------------------------------------------------------
    # 1) config.yaml 읽기
    # ------------------------------------------------------------------
    with CONFIG_PATH.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    base_device_cfg = dict(cfg["device"])
    compare_cfg = cfg["compare"]
    cv_cfg = cfg["sweeps"]["cv"]

    # ------------------------------------------------------------------
    # 2) compare:에 적힌 각 소자를 차례대로 시뮬레이션
    # ------------------------------------------------------------------
    curves = []

    for name, overrides in compare_cfg.items():
        # 기본 소자 설정을 복사한 뒤 compare:의 값만 변경한다.
        device_cfg = dict(base_device_cfg)
        device_cfg.update(overrides)

        device = Device(**device_cfg)

        print(f"[C-V] running: {name} -> {overrides}")

        curve = run_cv(
            device,
            start_v=float(cv_cfg["start_v"]),
            stop_v=float(cv_cfg["stop_v"]),
            step_v=float(cv_cfg["step_v"]),
        )

        curve["device"] = name
        curves.append(curve)

    result = pd.concat(curves, ignore_index=True)

    # 원래 계산 단위(F/um)는 CSV에 그대로 저장한다.
    result.to_csv("compare_cv.csv", index=False)
    print("saved: compare_cv.csv")

    # ------------------------------------------------------------------
    # 3) 그래프 표시
    # ------------------------------------------------------------------
    # 그래프에서는 숫자를 읽기 쉽게 fF/um 단위로 변환한다.
    plot_data = result.copy()
    plot_data["Cgg_fF_per_um"] = plot_data["Cgg_F_per_um"] * 1.0e15

    fig = px.line(
        plot_data,
        x="Vg_V",
        y="Cgg_fF_per_um",
        color="device",
        markers=True,
        title="Gate Oxide Thickness에 따른 Cgg-Vg 비교",
        labels={
            "Vg_V": "Vg (V)",
            "Cgg_fF_per_um": "Cgg (fF/µm)",
            "device": "Device",
        },
    )

    fig.update_layout(
        height=520,
        margin=dict(l=50, r=30, t=70, b=50),
    )

    fig.show()


if __name__ == "__main__":
    main()
