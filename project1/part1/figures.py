# -*- coding: utf-8 -*-
"""보고서용 그림. 최종 소자 구조 단면과 특성 곡선을 HTML 로 뽑는다.

    python -m part1.figures --design A

fig.show() 대신 write_html 을 쓴다. 리눅스 헤드리스에서는 브라우저가 안 열리고,
보고서에 붙일 때도 파일이 남아 있는 편이 낫다.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import io
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

import devsim

from . import spec
from .designs import DESIGNS, FINAL
from .body_tap import TapPart1Simulator as Part1Simulator   # 제출 구조(p+ 탭)

OUT = Path(__file__).resolve().parent.parent / "results"
UM = 1.0e-4   # cm per um


def structure_figure(device, name="fig"):
    """영역별 메시 노드를 뿌려 소자 단면을 그린다. y 는 아래로 증가(코드 규약)."""
    with contextlib.redirect_stdout(io.StringIO()):
        sim = Part1Simulator(device, name=name)
        sim._clear_session()
        sim._build_mesh()
        sim._build_doping()

    fig = go.Figure()
    colors = {"bulk": "#c8d6e5", "oxide": "#ee5253", "gate_metal": "#341f97", "air": "#f5f6fa"}
    for region in devsim.get_region_list(device=name):
        x = np.asarray(devsim.get_node_model_values(device=name, region=region, name="x")) / UM
        y = np.asarray(devsim.get_node_model_values(device=name, region=region, name="y")) / UM
        mat = devsim.get_material(device=name, region=region)
        fig.add_trace(go.Scatter(
            x=x, y=y, mode="markers", name=f"{region} ({mat})",
            marker=dict(size=2.5, color=colors.get(region, "#8395a7"))))

    # 도핑 부호로 n+ 영역을 덧칠한다
    xb = np.asarray(devsim.get_node_model_values(device=name, region="bulk", name="x")) / UM
    yb = np.asarray(devsim.get_node_model_values(device=name, region="bulk", name="y")) / UM
    nd = np.asarray(devsim.get_node_model_values(device=name, region="bulk", name="NetDoping"))
    m = nd > 0
    fig.add_trace(go.Scatter(x=xb[m], y=yb[m], mode="markers", name="n+ (NetDoping > 0)",
                             marker=dict(size=2.5, color="#10ac84")))

    d = device
    fig.update_layout(
        title=(f"Part 1 최종 소자 단면 — Lg {d.gate_length_um} µm, tox {d.oxide_thickness_nm} nm, "
               f"xj {d.junction_depth_um} µm, t_Si {d.silicon_thickness_um} µm, "
               f"NA {d.body_doping_cm3:.1e}, ND {d.sd_doping_cm3:.1e}, gate {d.gate_material}"),
        xaxis_title="x (µm)", yaxis_title="y (µm)  ·  아래로 증가, 실리콘 표면 y = 0",
        height=620, legend=dict(orientation="h", y=-0.18))
    fig.update_yaxes(autorange="reversed", scaleanchor="x", scaleratio=1.0)
    return fig


def curve_figures(device):
    """Id-Vg 를 선형/로그 두 장으로. Vd = 0.05 V 와 2 V 를 함께 그린다."""
    frames = []
    with contextlib.redirect_stdout(io.StringIO()):
        for label, vd, desc in (("Vd = 0.05 V", spec.VD_LIN, False), ("Vd = 2 V", spec.VDD, True)):
            sim = Part1Simulator(dataclasses.replace(device, temperature_k=spec.T_COLD), name="cv")
            sim.build(); sim.solve_equilibrium(); sim.enable_transport()
            vgs, ids = spec.sweep_for_vth(sim, vd, descending=desc)
            frames.append(pd.DataFrame({"Vg_V": vgs, "Id_A_per_um": ids, "bias": label}))
        # 고온
        sim = Part1Simulator(dataclasses.replace(device, temperature_k=spec.T_HOT), name="ch")
        sim.build(); sim.solve_equilibrium(); sim.enable_transport()
        vgs, ids = spec.sweep_for_vth(sim, spec.VDD, descending=True)
        frames.append(pd.DataFrame({"Vg_V": vgs, "Id_A_per_um": ids, "bias": "Vd = 2 V, 125 °C"}))

    df = pd.concat(frames, ignore_index=True).sort_values(["bias", "Vg_V"])
    lin = px.line(df, x="Vg_V", y="Id_A_per_um", color="bias", markers=True,
                  title="Id–Vg (선형)", labels={"Vg_V": "Vg (V)", "Id_A_per_um": "Id (A/µm)"})
    log_df = df[df["Id_A_per_um"] > 0]
    log = px.line(log_df, x="Vg_V", y="Id_A_per_um", color="bias", markers=True, log_y=True,
                  title="Id–Vg (로그) — SS · Ioff · 고온 Ioff 판정용",
                  labels={"Vg_V": "Vg (V)", "Id_A_per_um": "Id (A/µm)"})
    for f in (lin, log):
        f.update_layout(height=520, margin=dict(l=60, r=30, t=60, b=50))
    return lin, log, df


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", default=FINAL, choices=sorted(DESIGNS))
    args = ap.parse_args()
    device = DESIGNS[args.design]
    OUT.mkdir(exist_ok=True)

    tag = args.design
    f = structure_figure(device)
    p = OUT / f"fig_structure_{tag}.html"
    f.write_html(str(p)); print("saved:", p)

    lin, log, df = curve_figures(device)
    for nm, fig in (("idvg_lin", lin), ("idvg_log", log)):
        p = OUT / f"fig_{nm}_{tag}.html"
        fig.write_html(str(p)); print("saved:", p)
    p = OUT / f"curves_{tag}.csv"
    df.to_csv(p, index=False); print("saved:", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
