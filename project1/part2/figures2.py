# -*- coding: utf-8 -*-
"""Part 2 보고서용 그림: 1T1C 셀 단면과 쓰기 궤적.

    python -m part2.figures2 --design A
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import io
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

import devsim

from part2 import spec2
from part2.cell import CELLS, FINAL, Part2Simulator

OUT = Path(__file__).resolve().parent.parent / "results"
UM = 1.0e-4

COLORS = {"bulk": "#c8d6e5", "oxide": "#ee5253", "gate_metal": "#341f97",
          "hk_l": "#ff9f43", "hk_r": "#ff9f43", "metal_pillar": "#222f3e",
          "air": "#f5f6fa"}


def structure_figure(device, name="fig2"):
    with contextlib.redirect_stdout(io.StringIO()):
        sim = Part2Simulator(device, name=name)
        sim._clear_session(); sim._build_mesh(); sim._build_doping()
    fig = go.Figure()
    for region in devsim.get_region_list(device=name):
        x = np.asarray(devsim.get_node_model_values(device=name, region=region, name="x")) / UM
        y = np.asarray(devsim.get_node_model_values(device=name, region=region, name="y")) / UM
        mat = devsim.get_material(device=name, region=region)
        fig.add_trace(go.Scatter(x=x, y=y, mode="markers", name=f"{region} ({mat})",
                                 marker=dict(size=2.5, color=COLORS.get(region, "#8395a7"))))
    xb = np.asarray(devsim.get_node_model_values(device=name, region="bulk", name="x")) / UM
    yb = np.asarray(devsim.get_node_model_values(device=name, region="bulk", name="y")) / UM
    nd = np.asarray(devsim.get_node_model_values(device=name, region="bulk", name="NetDoping"))
    m = nd > 0
    fig.add_trace(go.Scatter(x=xb[m], y=yb[m], mode="markers", name="n+ (NetDoping > 0)",
                             marker=dict(size=2.5, color="#10ac84")))
    d = device
    fig.update_layout(
        title=(f"Part 2 1T1C 셀 — TRPASS Lg {d.gate_length_um} µm / tox {d.oxide_thickness_nm} nm / "
               f"NA {d.body_doping_cm3:.0e} / gate {d.gate_material} · "
               f"CAP h {d.cap_height_um} µm / {d.cap_material} {d.cap_dielectric_nm} nm / "
               f"기둥 {d.cap_pillar_nm} nm"),
        xaxis_title="x (µm)", yaxis_title="y (µm) · 아래로 증가, 실리콘 표면 y = 0",
        height=700, legend=dict(orientation="h", y=-0.14))
    fig.update_yaxes(autorange="reversed", scaleanchor="x", scaleratio=1.0)
    return fig


def write_figure(device, c_store_fF):
    """쓰기 궤적: 셀 전압 대 시간, 5 ns 지점 표시."""
    with contextlib.redirect_stdout(io.StringIO()):
        sim = Part2Simulator(dataclasses.replace(device, temperature_k=spec2.T_HOT), name="wf")
        sim.build(); sim.solve_equilibrium(); sim.enable_transport()
        v_high, (vs, is_) = spec2.write_level(sim, c_store_fF)
    c_f = c_store_fF * 1e-15
    i_a = np.asarray(is_) * spec2.W_UM
    inv = 1.0 / i_a
    t_ns = c_f * np.concatenate(([0.0], np.cumsum(np.diff(vs) * 0.5 * (inv[1:] + inv[:-1])))) * 1e9
    df = pd.DataFrame({"t_ns": t_ns, "V_cell_V": vs, "I_node_A_per_um": is_})
    fig = px.line(df, x="t_ns", y="V_cell_V", markers=True,
                  title=f"쓰기 궤적 (WL 2.5 V, 비트라인 2 V, 125 °C) — 5 ns 에서 V_high = {v_high:.4f} V",
                  labels={"t_ns": "시간 (ns)", "V_cell_V": "셀 전압 (V)"})
    fig.add_vline(x=5.0, line_dash="dash", annotation_text="WL 펄스 끝 5 ns")
    fig.add_hline(y=v_high, line_dash="dot", annotation_text=f"V_high {v_high:.3f} V")
    fig.update_layout(height=520, margin=dict(l=60, r=30, t=70, b=50))
    return fig, df, v_high


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default=FINAL, choices=sorted(CELLS))
    args = ap.parse_args()
    device = CELLS[args.design]
    OUT.mkdir(exist_ok=True)
    tag = args.design

    f = structure_figure(device)
    p = OUT / f"fig2_structure_{tag}.html"; f.write_html(str(p)); print("saved:", p)

    with contextlib.redirect_stdout(io.StringIO()):
        sim = Part2Simulator(device, name="cs2"); sim.build(); sim.solve_equilibrium()
        _, c_fF = sim.c_store()
    print(f"C_STORE = {c_fF:.4f} fF")

    fig, df, v_high = write_figure(device, c_fF)
    p = OUT / f"fig2_write_{tag}.html"; fig.write_html(str(p)); print("saved:", p)
    p = OUT / f"write_trajectory_{tag}.csv"; df.to_csv(p, index=False, encoding="utf-8-sig")
    print("saved:", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
