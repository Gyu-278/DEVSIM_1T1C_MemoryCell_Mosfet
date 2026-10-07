# -*- coding: utf-8 -*-
"""[3단계] DEVSIM TCAD 시뮬레이터를 Streamlit GUI 툴로 포장한다.

Device 1 / Device 2의 소자 파라미터를 나란히 설정하고,
Id-Vg / Id-Vd / C-V 해석을 선택해 두 소자의 특성을 비교한다.

2단계에서 구현한 mosfet_tool의 Device와 workflow를 그대로 재사용한다.
실행: MY_APP.bat 더블클릭
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

# ---------------------------------------------------------------------------
# 2단계에서 만든 TCAD 엔진을 재사용한다.
# step3/my_app.py 기준으로 상위 폴더의 step2를 Python import 경로에 추가한다.
# ---------------------------------------------------------------------------
STEP2_DIR = Path(__file__).resolve().parent.parent / "step2"
sys.path.insert(0, str(STEP2_DIR))

from mosfet_tool.config import Device
from mosfet_tool.workflows import run_cv, run_idvd, run_idvg


# ---------------------------------------------------------------------------
# 페이지 설정
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="My TCAD MOSFET Tool",
    page_icon="🔬",
    layout="wide",
)

st.title("My TCAD MOSFET Tool")
st.caption(
    "DEVSIM 2-D nMOS TCAD · Device 1 / Device 2의 소자 파라미터와 "
    "Id–Vg / Id–Vd / C–V 특성을 비교합니다."
)


# ---------------------------------------------------------------------------
# 소자 설정 패널
# ---------------------------------------------------------------------------
def device_panel(column, label: str, default_tox: float) -> Device:
    """한쪽 컬럼에 소자 파라미터 입력 GUI를 만들고 Device 객체를 반환한다."""

    with column:
        st.subheader(label)

        gate_length_um = st.slider(
            "Gate length (µm)",
            min_value=0.5,
            max_value=2.0,
            value=1.0,
            step=0.1,
            key=f"{label}-lg",
        )

        oxide_thickness_nm = st.slider(
            "Oxide thickness (nm)",
            min_value=5.0,
            max_value=30.0,
            value=default_tox,
            step=1.0,
            key=f"{label}-tox",
        )

        # 도핑 농도는 값의 범위가 매우 크므로 log10(N) 형태로 입력한다.
        body_doping_log = st.slider(
            "Body doping log10(N / cm³)",
            min_value=15.0,
            max_value=18.0,
            value=16.0,
            step=0.1,
            key=f"{label}-body-doping",
        )

        sd_doping_log = st.slider(
            "Source/Drain doping log10(N / cm³)",
            min_value=18.0,
            max_value=20.0,
            value=19.0,
            step=0.1,
            key=f"{label}-sd-doping",
        )

        with st.expander("More device settings"):
            source_length_um = st.slider(
                "Source length (µm)",
                min_value=0.2,
                max_value=1.0,
                value=0.5,
                step=0.1,
                key=f"{label}-ls",
            )

            drain_length_um = st.slider(
                "Drain length (µm)",
                min_value=0.2,
                max_value=1.0,
                value=0.5,
                step=0.1,
                key=f"{label}-ld",
            )

            silicon_thickness_um = st.slider(
                "Silicon thickness (µm)",
                min_value=0.2,
                max_value=1.0,
                value=0.5,
                step=0.1,
                key=f"{label}-tsi",
            )

            junction_depth_um = st.slider(
                "Junction depth (µm)",
                min_value=0.05,
                max_value=0.30,
                value=0.10,
                step=0.01,
                key=f"{label}-xj",
            )

            temperature_k = st.slider(
                "Temperature (K)",
                min_value=250.0,
                max_value=400.0,
                value=300.0,
                step=10.0,
                key=f"{label}-temp",
            )

        return Device(
            gate_length_um=gate_length_um,
            source_length_um=source_length_um,
            drain_length_um=drain_length_um,
            oxide_thickness_nm=oxide_thickness_nm,
            silicon_thickness_um=silicon_thickness_um,
            junction_depth_um=junction_depth_um,
            body_doping_cm3=10.0 ** body_doping_log,
            sd_doping_cm3=10.0 ** sd_doping_log,
            temperature_k=temperature_k,
        )


col1, col2 = st.columns(2)

# 기본값을 다르게 두어 처음 실행했을 때에도 두 소자의 차이가 보이도록 한다.
device1 = device_panel(col1, "Device 1", default_tox=10.0)
device2 = device_panel(col2, "Device 2", default_tox=20.0)


# ---------------------------------------------------------------------------
# 해석 종류 및 sweep 조건
# ---------------------------------------------------------------------------
st.divider()
st.subheader("Analysis")

analysis = st.selectbox(
    "Analysis type",
    ("Id–Vg", "Id–Vd", "C–V"),
)

set_col1, set_col2, set_col3, set_col4 = st.columns(4)

if analysis == "Id–Vg":
    with set_col1:
        fixed_voltage = st.number_input(
            "Drain voltage Vd (V)",
            min_value=0.0,
            max_value=2.0,
            value=0.05,
            step=0.05,
        )
    default_start, default_stop, default_step = 0.0, 2.0, 0.1

elif analysis == "Id–Vd":
    with set_col1:
        fixed_voltage = st.number_input(
            "Gate voltage Vg (V)",
            min_value=0.0,
            max_value=3.0,
            value=2.0,
            step=0.1,
        )
    default_start, default_stop, default_step = 0.0, 2.0, 0.1

else:
    fixed_voltage = None
    with set_col1:
        st.info("C–V: gate voltage sweep")
    default_start, default_stop, default_step = -1.0, 2.0, 0.1

with set_col2:
    start_v = st.number_input(
        "Sweep start (V)",
        value=float(default_start),
        step=0.1,
        key=f"{analysis}-start",
    )

with set_col3:
    stop_v = st.number_input(
        "Sweep stop (V)",
        value=float(default_stop),
        step=0.1,
        key=f"{analysis}-stop",
    )

with set_col4:
    step_v = st.number_input(
        "Sweep step (V)",
        min_value=0.01,
        max_value=1.0,
        value=float(default_step),
        step=0.01,
        key=f"{analysis}-step",
    )


# ---------------------------------------------------------------------------
# 시뮬레이션
#
# Streamlit은 slider를 움직일 때마다 파일 전체를 다시 실행한다.
# TCAD 해석은 시간이 걸리므로 버튼을 눌렀을 때만 계산하고,
# 결과는 session_state에 저장해 그래프가 사라지지 않도록 한다.
# ---------------------------------------------------------------------------
run_clicked = st.button("Run simulation", type="primary", use_container_width=True)

if run_clicked:
    if start_v == stop_v:
        st.error("Sweep start와 stop은 서로 다른 값이어야 합니다.")
    else:
        curves = []

        try:
            with st.spinner(f"{analysis} TCAD simulation running..."):
                for name, device in {
                    "Device 1": device1,
                    "Device 2": device2,
                }.items():

                    if analysis == "Id–Vg":
                        curve = run_idvg(
                            device,
                            drain_v=float(fixed_voltage),
                            start_v=float(start_v),
                            stop_v=float(stop_v),
                            step_v=float(step_v),
                        )

                    elif analysis == "Id–Vd":
                        curve = run_idvd(
                            device,
                            gate_v=float(fixed_voltage),
                            start_v=float(start_v),
                            stop_v=float(stop_v),
                            step_v=float(step_v),
                        )

                    else:
                        # quasi-static C-V에서는 전류가 필요하지 않으므로
                        # run_cv() 내부에서 enable_transport()를 호출하지 않는다.
                        curve = run_cv(
                            device,
                            start_v=float(start_v),
                            stop_v=float(stop_v),
                            step_v=float(step_v),
                        )

                    curve["device"] = name
                    curves.append(curve)

            result = pd.concat(curves, ignore_index=True)

            # 계산 당시의 해석 종류/조건과 함께 저장한다.
            st.session_state["tcad_result"] = result
            st.session_state["tcad_analysis"] = analysis
            st.session_state["tcad_fixed_voltage"] = fixed_voltage
            st.session_state["tcad_device1"] = device1
            st.session_state["tcad_device2"] = device2

        except Exception as exc:
            st.error("Simulation failed.")
            st.exception(exc)


# ---------------------------------------------------------------------------
# 결과 표시
# ---------------------------------------------------------------------------
result = st.session_state.get("tcad_result")
result_analysis = st.session_state.get("tcad_analysis")

if result is not None and result_analysis is not None:
    st.divider()
    st.subheader(f"Result — {result_analysis}")

    plot_data = result.copy()

    if result_analysis == "Id–Vg":
        x_col = "Vg_V"
        y_col = "Id_A_per_um"
        x_label = "Vg (V)"
        y_label = "Id (A/µm)"

        fixed = st.session_state.get("tcad_fixed_voltage")
        st.caption(f"Simulation condition: Vd = {fixed:g} V")

    elif result_analysis == "Id–Vd":
        x_col = "Vd_V"
        y_col = "Id_A_per_um"
        x_label = "Vd (V)"
        y_label = "Id (A/µm)"

        fixed = st.session_state.get("tcad_fixed_voltage")
        st.caption(f"Simulation condition: Vg = {fixed:g} V")

    else:
        x_col = "Vg_V"

        # F/µm보다 fF/µm가 그래프에서 읽기 쉬우므로 표시할 때만 단위를 바꾼다.
        plot_data["Cgg_fF_per_um"] = plot_data["Cgg_F_per_um"] * 1.0e15
        y_col = "Cgg_fF_per_um"
        x_label = "Vg (V)"
        y_label = "Cgg (fF/µm)"

    fig = px.line(
        plot_data,
        x=x_col,
        y=y_col,
        color="device",
        markers=True,
        labels={
            x_col: x_label,
            y_col: y_label,
            "device": "Device",
        },
    )
    fig.update_layout(
        height=470,
        margin=dict(l=30, r=20, t=20, b=30),
    )
    st.plotly_chart(fig, use_container_width=True)

    # 수치 확인이 가능하도록 원래 단위의 DataFrame도 함께 표시한다.
    with st.expander("Show simulation data"):
        st.dataframe(result, use_container_width=True)

        csv_data = result.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download result CSV",
            data=csv_data,
            file_name=f"{result_analysis.replace('–', '').replace('-', '').lower()}_compare.csv",
            mime="text/csv",
        )
