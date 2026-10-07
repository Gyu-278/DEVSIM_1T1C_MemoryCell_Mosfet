# -*- coding: utf-8 -*-
"""[Part 2 / 2026-09-27] 보고서용 그림.

과제 2-(b)가 요구하는 것: "ΔVBL(t) curves for '0' and '1'; value at 0.5 ns vs 60 mV".
2-(a)는 1T1C 구조 그림, 2-(d)는 최종 셀 그림이다.

    python -m part2.figures_p2 --cell C4
"""

from __future__ import annotations

import argparse
import contextlib
import io
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from part2 import transient as T
from part2.body_tap import TapPart2Simulator as Part2Simulator, use_tap   # 제출 구조(p+ 탭)
from part2.cell import CELLS, FINAL
from part2.figures2 import structure_figure

OUT = Path(__file__).resolve().parent.parent / "results"


def read_curves(device, c_store_F):
    """'0' 과 '1' 의 과도 읽기 궤적을 모은다."""
    frames = []
    for state, v0 in (("0", 0.0), ("1", 2.0)):
        with contextlib.redirect_stdout(io.StringIO()):
            sim = T.build(device, name=f"fig{state}")
            r = T.read_transient(sim, v0, c_store_F)
        frames.append(pd.DataFrame({
            "t_ps": r["t"] * 1e12,
            "V_BL_V": r["v_bl"],
            "dV_BL_mV": (r["v_bl"] - T.V_BL_PRE) * 1e3,
            "V_cell_V": r["v_cell"],
            "I_source_uA_per_um": r["i_source"] * 1e6,
            "state": f'data "{state}"',
        }))
    return pd.concat(frames, ignore_index=True)


def read_figure(df, device):
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.09,
                        subplot_titles=("비트라인 전압 변화 ΔV_BL(t)",
                                        "저장 노드 전류 I(t)"))
    colors = {'data "0"': "#2e86de", 'data "1"': "#ee5253"}
    for st, g in df.groupby("state"):
        fig.add_trace(go.Scatter(x=g["t_ps"], y=g["dV_BL_mV"], name=st, mode="lines+markers",
                                 line=dict(color=colors[st]), marker=dict(size=4)), row=1, col=1)
        fig.add_trace(go.Scatter(x=g["t_ps"], y=g["I_source_uA_per_um"], name=st + " I",
                                 mode="lines", line=dict(color=colors[st], dash="dot"),
                                 showlegend=False), row=2, col=1)
    # 판정선: 0.5 ns 와 ±60 mV
    for y in (60.0, -60.0):
        fig.add_hline(y=y, line_dash="dash", line_color="#8395a7", row=1, col=1)
    fig.add_vline(x=500.0, line_dash="dash", line_color="#222f3e")
    fig.add_annotation(x=500, y=60, text="판정 시점 0.5 ns · 기준 ±60 mV",
                       showarrow=False, yshift=14, row=1, col=1)
    d = device
    fig.update_layout(
        title=(f"Part 2 과도 읽기 — TRPASS Lg {d.gate_length_um} µm / tox {d.oxide_thickness_nm} nm / "
               f"NA {d.body_doping_cm3:.1e} / ND {d.sd_doping_cm3:.0e} / gate {d.gate_material} · "
               f"CAP {d.cap_material} h {d.cap_height_um} µm / t {d.cap_dielectric_nm} nm · "
               f"125 °C, V_B = −0.5 V"),
        height=760, legend=dict(orientation="h", y=-0.1))
    fig.update_xaxes(title_text="시간 (ps)", row=2, col=1)
    fig.update_yaxes(title_text="ΔV_BL (mV)", row=1, col=1)
    fig.update_yaxes(title_text="I (µA/µm)", row=2, col=1)
    return fig


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", default=FINAL)
    args = ap.parse_args()
    device = CELLS[args.cell]
    OUT.mkdir(exist_ok=True)
    tag = args.cell
    use_tap()                   # T.build · structure_figure 도 탭 구조로

    with contextlib.redirect_stdout(io.StringIO()):
        sim = Part2Simulator(device, name="cs3"); sim.build(); sim.solve_equilibrium()
        _, c_fF = sim.c_store(v0=T.V_PLATE)
    print(f"C_STORE = {c_fF:.4f} fF  (해석식 {device.c_store_closed_form_fF:.4f})")

    f = structure_figure(device, name="p2fig")
    p = OUT / f"fig_p2_structure_{tag}.html"; f.write_html(str(p)); print("saved:", p)

    df = read_curves(device, c_fF * 1e-15)
    fig = read_figure(df, device)
    p = OUT / f"fig_p2_read_{tag}.html"; fig.write_html(str(p)); print("saved:", p)
    p = OUT / f"p2_read_curves_{tag}.csv"
    df.to_csv(p, index=False, encoding="utf-8-sig"); print("saved:", p)

    for st, g in df.groupby("state"):
        last = g.iloc[-1]
        print(f"  {st}: 0.5 ns 에서 ΔV_BL = {last['dV_BL_mV']:+.2f} mV "
              f"({'PASS' if abs(last['dV_BL_mV']) >= 60 else 'FAIL'}), "
              f"셀 {last['V_cell_V']:.4f} V")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
