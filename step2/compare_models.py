# -*- coding: utf-8 -*-
"""Step 1 square-law model과 Step 2 DEVSIM의 Id-Vg를 동일 조건에서 비교한다.

조건:
- Vd = 0.05 V
- Vg = 0.0 ~ 2.0 V
- Step 1: SimpleMosfet(vth=0.6, k=1.4e-4)
- Step 2: config.yaml의 기본 device 설정

실행: STEP2_MODEL_COMPARE.bat 더블클릭
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import plotly.express as px
import yaml

STEP1_DIR = Path(__file__).resolve().parent.parent / "step1"
sys.path.insert(0, str(STEP1_DIR))
from simulator_practice import SimpleMosfet

from mosfet_tool.config import Device
from mosfet_tool.workflows import run_idvg

CONFIG_PATH = Path(__file__).resolve().parent / "config.yaml"


def main():
    with CONFIG_PATH.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    idvg_cfg = cfg["sweeps"]["idvg"]
    vd = float(idvg_cfg["drain_v"])
    start_v = float(idvg_cfg["start_v"])
    stop_v = float(idvg_cfg["stop_v"])
    step_v = float(idvg_cfg["step_v"])

    voltages = np.round(
        np.arange(start_v, stop_v + 0.5 * step_v, step_v),
        9,
    )

    simple = SimpleMosfet(vth=0.6, k=1.4e-4)
    simple_curve = simple.sweep_idvg(vd=vd, voltages=voltages)
    simple_curve["model"] = "Step 1: Square-law"

    device = Device(**cfg["device"])
    devsim_curve = run_idvg(
        device,
        drain_v=vd,
        start_v=start_v,
        stop_v=stop_v,
        step_v=step_v,
    )
    devsim_curve["model"] = "Step 2: DEVSIM"

    result = pd.concat([simple_curve, devsim_curve], ignore_index=True)
    result.to_csv("compare_simple_devsim_idvg.csv", index=False)
    print("saved: compare_simple_devsim_idvg.csv")

    fig = px.line(
        result,
        x="Vg_V",
        y="Id_A_per_um",
        color="model",
        markers=True,
        title=f"Square-law vs DEVSIM Id-Vg Comparison (Vd = {vd:.2f} V)",
        labels={
            "Vg_V": "Vg (V)",
            "Id_A_per_um": "Id (A/µm)",
            "model": "Model",
        },
    )
    fig.update_layout(height=520, margin=dict(l=50, r=30, t=70, b=50))
    fig.show()

    log_data = result.copy()
    log_data.loc[log_data["Id_A_per_um"] <= 0.0, "Id_A_per_um"] = np.nan

    fig_log = px.line(
        log_data,
        x="Vg_V",
        y="Id_A_per_um",
        color="model",
        markers=True,
        log_y=True,
        title=f"Square-law vs DEVSIM Id-Vg Comparison - Log Scale (Vd = {vd:.2f} V)",
        labels={
            "Vg_V": "Vg (V)",
            "Id_A_per_um": "Id (A/µm)",
            "model": "Model",
        },
    )
    fig_log.update_layout(height=520, margin=dict(l=50, r=30, t=70, b=50))
    fig_log.show()


if __name__ == "__main__":
    main()
