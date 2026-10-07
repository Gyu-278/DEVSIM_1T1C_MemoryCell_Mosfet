# -*- coding: utf-8 -*-
"""보고서 그림 생성 v2 (라벨 겹침 정리, 구조 그림 재작성).

    cd project1 && /hai/anaconda3/bin/python make_figures2.py
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import matplotlib.tri as mtri
import numpy as np
from matplotlib import font_manager as fm

fm.fontManager.addfont("/usr/share/fonts/google-noto-cjk/NotoSansCJK-Regular.ttc")
fm.fontManager.addfont("/usr/share/fonts/google-noto-cjk/NotoSansCJK-Bold.ttc")
plt.rcParams.update({
    "font.family": "Noto Sans CJK JP", "axes.unicode_minus": False,
    "font.size": 9.5, "axes.titlesize": 10.5, "axes.labelsize": 9.5, "legend.fontsize": 8.5,
    "savefig.dpi": 200, "axes.grid": True, "grid.alpha": 0.25,
    "mathtext.fontset": "dejavusans",
})

R, D, F = Path("results"), Path("report/data"), Path("report/fig")
F.mkdir(parents=True, exist_ok=True)
C0, C1 = "#1f6fb4", "#d1352b"
PASS, FAIL = "#2e9e5b", "#d1352b"
GATE_COL = {"TaN": "#8c6bb1", "TiN": "#2b8cbe", "W": "#e6550d", "Mo": "#31a354"}
NPLUS, PSUB, OXI, GM, HK, PIL = "#9fd8a8", "#dbe9f6", "#f4a6a0", "#b39ddb", "#f8c471", "#7f8c8d"
PTAP, TAPC = "#c9a7e0", "#6a1b9a"   # p+ 탭, body 접점 색


def save(fig, name):
    p = F / name
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print("saved", p)


def read_csvout(path):
    with open(path, encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    hdr, out = rows[0], {}
    for r in rows[1:]:
        if r and r[0] != "기준":
            out[r[0]] = dict(zip(hdr, r))
    return out


def fnum(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def ann(ax, text, xy, xytext, fs=8.5, ha="left", va="center", color="k", **kw):
    ax.annotate(text, xy=xy, xytext=xytext, fontsize=fs, ha=ha, va=va, color=color,
                arrowprops=dict(arrowstyle="->", lw=0.8, color=color, shrinkA=2, shrinkB=2), **kw)


def dim_h(ax, x0, x1, y, text, fs=8.5, above=True):
    ax.annotate("", xy=(x0, y), xytext=(x1, y), arrowprops=dict(arrowstyle="<->", lw=0.9))
    ax.text((x0 + x1) / 2, y + (-0.012 if above else 0.012), text, ha="center",
            va="bottom" if above else "top", fontsize=fs)


def dim_v(ax, x, y0, y1, text, fs=8.5, side="right"):
    ax.annotate("", xy=(x, y0), xytext=(x, y1), arrowprops=dict(arrowstyle="<->", lw=0.9))
    ax.text(x + (0.02 if side == "right" else -0.02), (y0 + y1) / 2, text,
            ha="left" if side == "right" else "right", va="center", fontsize=fs)


def doping_map(ax, struct_csv, title, xlim, ylim, tap=False):
    xs, ys, ns = [], [], []
    with open(struct_csv) as fh:
        for r in csv.DictReader(fh):
            if r["region"] == "bulk":
                xs.append(float(r["x_um"])); ys.append(float(r["y_um"])); ns.append(float(r["NetDoping"]))
    xs, ys, ns = map(np.asarray, (xs, ys, ns))
    v = np.sign(ns) * np.log10(np.maximum(np.abs(ns), 1e14))
    v = np.where(np.abs(v) < 14.0, 0.0, v)
    tri = mtri.Triangulation(xs, ys)
    cs = ax.tricontourf(tri, v, levels=np.linspace(-21, 21, 43), cmap="RdBu_r", extend="both")
    ax.tricontour(tri, ns, levels=[0.0], colors="k", linewidths=0.9)
    if tap:                                   # p+ 탭 경계 (순 도핑 −1e19)
        ax.tricontour(tri, ns, levels=[-1e19], colors=TAPC, linewidths=1.4)
    ax.set_aspect("equal"); ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)"); ax.set_title(title); ax.grid(False)
    return cs


def doping_colorbar(fig, cs, ax):
    cb = fig.colorbar(cs, ax=ax, shrink=0.85, ticks=[-20, -17, 0, 18, 20])
    cb.ax.set_yticklabels([r"p $10^{20}$", r"p $10^{17}$", "0", r"n $10^{18}$", r"n $10^{20}$"], fontsize=8)
    cb.set_label(r"순 도핑 (cm$^{-3}$), 파랑 = p형 · 빨강 = n형", fontsize=8.5)


# ---------------------------------------------------------------------------
def fig01_p1_structure():
    Ls, Lg, Ld, tsi = 0.50, 0.50, 0.50, 0.50
    xjm = 0.061                       # 금속학적 접합 깊이 (erfc 중심 0.05 + 꼬리)
    tox_d, gm_d = 0.03, 0.07          # 그림용 과장 두께
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.0, 4.9), gridspec_kw=dict(width_ratios=[1.12, 1]))
    X = Ls + Lg + Ld
    tl, tr_ = -0.40, -0.25            # p+ 탭 (폭 0.15 µm, 소스와 간격 0.25 µm)
    a.add_patch(mp.Rectangle((tl, 0), X - tl, tsi, fc=PSUB, ec="k", lw=0.8))
    a.add_patch(mp.Rectangle((tl, 0), tr_ - tl, xjm, fc=PTAP, ec=TAPC, lw=0.8))
    a.add_patch(mp.Rectangle((0, 0), Ls, xjm, fc=NPLUS, ec="#2e7d32", lw=0.8))
    a.add_patch(mp.Rectangle((Ls + Lg, 0), Ld, xjm, fc=NPLUS, ec="#2e7d32", lw=0.8))
    a.add_patch(mp.Rectangle((Ls, -tox_d), Lg, tox_d, fc=OXI, ec="#b71c1c", lw=0.8))
    a.add_patch(mp.Rectangle((Ls, -tox_d - gm_d), Lg, gm_d, fc=GM, ec="#4a148c", lw=0.8))
    a.text(Ls + Lg / 2, -tox_d - gm_d / 2, "게이트 (TaN)", ha="center", va="center", fontsize=8.5,
           color="white", weight="bold")
    for x0, x1 in ((0, Ls), (Ls + Lg, X)):
        a.plot([x0, x1], [-0.003, -0.003], color="#333", lw=3.2)
    a.plot([tl, tr_], [-0.003, -0.003], color=TAPC, lw=3.2)
    a.text(0.22, -0.025, "source 접점", ha="center", va="bottom", fontsize=8)
    a.text(1.28, -0.025, "drain 접점", ha="center", va="bottom", fontsize=8)
    a.text((tl + tr_) / 2, -0.025, "body 접점", ha="center", va="bottom", fontsize=8, color=TAPC)
    a.text(X / 2 - 0.1, tsi + 0.03, "바닥: 접점 없음", ha="center", va="top", fontsize=8)
    ann(a, "p+ 탭 1e20\n(폭 0.15 µm)", ((tl + tr_) / 2, 0.035), (-0.64, 0.24), fs=7.8, color=TAPC)
    a.annotate("", xy=(tr_, 0.11), xytext=(0.0, 0.11), arrowprops=dict(arrowstyle="<->", lw=0.8))
    a.text(tr_ / 2, 0.125, "간격 0.25 µm", ha="center", va="top", fontsize=7.2)
    ann(a, r"n+ 소스 (ND $10^{20}$, xj 0.05 µm)", (0.30, 0.03), (0.06, 0.155), fs=8)
    ann(a, r"n+ 드레인 (ND $10^{20}$)", (1.20, 0.03), (1.02, 0.12), fs=8)
    a.text(X / 2, 0.40, "p형 기판 (영역 bulk, Silicon)\n" + r"NA = 7.5$\times10^{16}$ cm$^{-3}$",
           ha="center", va="center", fontsize=9)
    dim_h(a, Ls, Ls + Lg, -0.17, "Lg = 0.50 µm")
    dim_h(a, 0, Ls, 0.29, "Ls = 0.50 µm")
    dim_h(a, Ls + Lg, X, 0.29, "Ld = 0.50 µm")
    dim_v(a, X + 0.07, 0, tsi, "t_Si =\n0.50 µm")
    ann(a, "게이트 산화막 SiO$_2$, tox = 5.0 nm\n(얇아서 그림에서는 과장)", (Ls + Lg * 0.95, -tox_d / 2),
        (1.12, -0.26), fs=8, va="top")
    ann(a, "gate_metal 영역: 재료 이름 TaN\n(유효 일함수 4.45 eV)", (Ls + 0.05, -tox_d - gm_d * 0.6),
        (-0.62, -0.26), fs=8, va="top")
    a.set_xlim(-0.66, 2.10); a.set_ylim(0.62, -0.32)
    a.set_aspect("equal"); a.grid(False)
    a.set_xticks([0, 0.5, 1.0, 1.5]); a.set_xlabel("x (µm)"); a.set_ylabel("y (µm)")
    a.set_title("(a) 구조와 설계 치수 (모식도)")
    cs = doping_map(b, D / "struct_p1_tap.csv", "(b) 실제 메시의 순 도핑 (보라 선 = p+ 탭)", (-0.4, 1.5), (0.5, -0.01), tap=True)
    doping_colorbar(fig, cs, b)
    b.text(0.75, 0.43, "검은 선 = p-n 접합 (순 도핑 0)", ha="center", fontsize=8,
           bbox=dict(fc="white", ec="none", alpha=0.8))
    save(fig, "fig01_p1_structure.png")


# ---------------------------------------------------------------------------
def fig02_p1_idvg():
    p = R / "curves_final.csv"
    if not p.exists():
        return
    cur = {}
    with open(p, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            cur.setdefault(r["bias"], []).append((float(r["Vg_V"]), float(r["Id_A_per_um"])))
    sj = json.load(open(R / "bodytap" / "p1_final_bodytap.json", encoding="utf-8"))["bodytap"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.8, 4.5))
    style = {"Vd = 0.05 V": ("#1f6fb4", "Vd = 0.05 V, 27 °C"), "Vd = 2 V": ("#d1352b", "Vd = 2 V, 27 °C"),
             "Vd = 2 V, 125 °C": ("#ff8c00", "Vd = 2 V, 125 °C")}
    for k, pts in cur.items():
        pts = sorted(pts)
        vg = np.array([x for x, _ in pts]); idd = np.array([y for _, y in pts])
        col, lab = style.get(k, ("gray", k))
        m = idd > 0
        a.semilogy(vg[m], idd[m], "-", color=col, lw=1.7, label=lab)
        b.plot(vg, idd * 1e6, "-", color=col, lw=1.7, label=lab)
    a.axhline(1e-7, color="k", ls="--", lw=0.8)
    a.text(1.98, 1.5e-7, r"Vth 정의: Id = $10^{-7}$ A/µm", fontsize=8, ha="right")
    a.axhspan(1e-10, 1e-8, color="#999", alpha=0.15)
    a.text(1.98, 1.4e-10, "SS 평균 구간\n" + r"$10^{-10}$ ~ $10^{-8}$ A/µm", fontsize=8, ha="right")
    a.axvline(sj["vth"], color="#1f6fb4", ls=":", lw=1)
    a.text(sj["vth"] + 0.03, 1e-15, f"Vth = {sj['vth']:.3f} V", color="#1f6fb4", fontsize=8.5)
    i27, i125 = sj["ioff_pA_per_um"] * 1e-12, sj["ioff_hot_pA_per_um"] * 1e-12
    a.plot([0], [i27], "o", color="#d1352b", zorder=5); a.plot([0], [i125], "o", color="#ff8c00", zorder=5)
    ann(a, f"Ioff(27 °C) = {sj['ioff_pA_per_um']:.3f} pA/µm", (0.0, i27), (0.55, 3e-14), fs=8)
    ann(a, f"Ioff(125 °C) = {sj['ioff_hot_pA_per_um']:.1f} pA/µm", (0.0, i125), (0.55, 4e-12), fs=8)
    a.set_xlabel("Vg (V)"); a.set_ylabel("Id (A/µm)"); a.set_ylim(1e-16, 3e-3)
    a.set_title("(a) 로그 축: Vth · SS · Ioff 판정"); a.legend(loc="lower right")
    b.plot([2.0], [sj["ion_uA_per_um"]], "o", color="#d1352b", zorder=5)
    ann(b, f"Ion = {sj['ion_uA_per_um']:.0f} µA/µm\n(Vg = Vd = 2 V)", (2.0, sj["ion_uA_per_um"]),
        (1.05, 560), fs=8.5)
    b.axhline(450, color="k", ls="--", lw=0.8); b.text(0.03, 462, "Ion 기준 450 µA/µm", fontsize=8)
    b.set_xlabel("Vg (V)"); b.set_ylabel("Id (µA/µm)")
    b.set_title("(b) 선형 축: Ion 판정"); b.legend(loc="upper left")
    save(fig, "fig02_p1_idvg.png")


# ---------------------------------------------------------------------------
def fig03_p1_physics():
    p = R / "physics_check_baseline.csv"
    if not p.exists():
        return
    rows = read_csvout(p)
    names = list(rows.keys())
    short = ["배포 코드\n그대로", "+ 게이트\n일함수", "+ n_i(T)", "+ μ(T)", "셋 다\n(최종)"]
    have_m, have_t = (D / "verify_metal.json").exists(), (D / "verify_temp.json").exists()
    ncol = 1 + int(have_m) + int(have_t)
    fig, axs = plt.subplots(1, ncol, figsize=(4.7 * ncol + 0.6, 4.3))
    axs = np.atleast_1d(axs)
    a = axs[0]
    x = np.arange(len(names))
    i3 = np.array([fnum(rows[n]["Ioff300K_pA_per_um"]) for n in names])
    i4 = np.array([fnum(rows[n]["Ioff398K_pA_per_um"]) for n in names])
    a.bar(x - 0.2, i3, 0.4, color="#1f6fb4", label="27 °C")
    a.bar(x + 0.2, i4, 0.4, color="#ff8c00", label="125 °C")
    a.set_yscale("log")
    for xi, u, w in zip(x, i3, i4):
        a.text(xi + 0.2, w * 1.7, f"×{w/u:,.0f}" if w / u >= 10 else f"×{w/u:.1f}", ha="center",
               fontsize=8.5, weight="bold", color="#b35400")
    a.set_xticks(x); a.set_xticklabels(short[:len(names)], fontsize=8.5)
    a.set_ylabel("Ioff (pA/µm) · Vg = 0, Vd = 2 V"); a.set_ylim(5e-4, 8e2)
    a.set_title("(a) 고온/상온 누설 비 (배포 기본 소자)"); a.legend(loc="upper left")
    k = 1
    if have_m:
        m = json.load(open(D / "verify_metal.json", encoding="utf-8"))
        b = axs[k]; k += 1
        pts = sorted((v["ewf"], v["vth"], n) for n, v in m.items() if v.get("ewf"))
        e = np.array([p[0] for p in pts]); vt = np.array([p[1] for p in pts])
        slope, icpt = np.polyfit(e, vt, 1)
        ee = np.linspace(e.min() - 0.05, e.max() + 0.05, 10)
        b.plot(ee, slope * ee + icpt, "--", color="#5b2a86", lw=1)
        b.plot(e, vt, "o", color="#5b2a86", ms=7)
        for ei, vi, n in pts:
            b.annotate(n, (ei, vi), xytext=(7, -9), textcoords="offset points", fontsize=8.5)
        mg = m.get("midgap(배포)")
        if mg:
            b.axhline(mg["vth"], color="gray", ls=":", lw=1)
            b.text(e.max(), mg["vth"] - 0.015, f"배포 코드 (일함수 항 없음) {mg['vth']:.3f} V",
                   fontsize=7.8, color="gray", ha="right", va="top")
        b.set_xlabel("게이트 금속 유효 일함수 φm (eV)"); b.set_ylabel("Vth (V) · Vd = 0.05 V")
        b.set_title(f"(b) 일함수 1 eV 당 Vth 변화 = {slope:.4f} V")
    if have_t:
        t = json.load(open(D / "verify_temp.json", encoding="utf-8"))
        c = axs[k]
        T = np.array(sorted(float(x) for x in t if x.replace('.', '', 1).isdigit()))
        key = lambda x: next(kk for kk in t if kk.replace('.', '', 1).isdigit() and abs(float(kk) - x) < 1e-6)
        v0 = np.array([t[key(x)]["vth_deployed"] for x in T]); v1 = np.array([t[key(x)]["vth_final"] for x in T])
        s0, s1 = np.polyfit(T, v0, 1)[0] * 1e3, np.polyfit(T, v1, 1)[0] * 1e3
        c.plot(T - 273, v0, "s--", color="gray", label=f"배포 코드: {s0:+.2f} mV/K\n(온도가 오르면 Vth 상승: 틀린 부호)")
        c.plot(T - 273, v1, "o-", color="#d1352b", label=f"물리 셋 다 추가: {s1:+.2f} mV/K\n(온도가 오르면 Vth 감소: 정상)")
        mo = t.get("mu_off")
        if mo:
            c.plot([27, 125], [mo["300"], mo["398"]], "^:", color="#7b3294",
                   label=f"μ(T) 만 끔: {mo['slope_mV_per_K']:+.2f} mV/K")
        c.set_xlabel("온도 (°C)"); c.set_ylabel("Vth (V) · Vd = 0.05 V")
        c.set_title("(c) Vth 의 온도 계수"); c.legend(loc="center left", fontsize=8)
    save(fig, "fig03_p1_physics.png")


# ---------------------------------------------------------------------------
def fig04_p1_search():
    p = R / "search_record.csv"
    if not p.exists():
        return
    rows = read_csvout(p)
    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    ax.axvspan(450, 720, color=PASS, alpha=0.06)
    ax.axvline(450, color="k", ls="--", lw=0.8); ax.axhline(100, color="k", ls="--", lw=0.8)
    ax.text(456, 116, "Ion ≥ 450 µA/µm", fontsize=8, va="center")
    ax.text(105, 85, "고온 Ioff ≤ 100 pA/µm", fontsize=8, va="top")
    # (라벨, 점에서 띄운 거리 [pt], 지시선 여부). 겹치지 않게 손으로 배치했다.
    lab = {"baseline": ("배포 기본 소자", (6, 4), False), "A": ("후보 A", (8, 6), False),
           "B": ("후보 B", (7, -3), False), "final": ("최종", (10, 6), False),
           "A1_na5e16": ("A, NA 5e16", (-62, -2), False), "A2_na6e16": ("A, NA 6e16", (7, 3), False),
           "A4_tox5_na6": ("A, tox 5 + NA 6e16", (14, 30), True, "right"), "A5_lg06": ("A, Lg 0.6", (-62, -3), False),
           "A6_tox5_lg06": ("A, tox 5 + Lg 0.6", (-8, -32), True), "A7_xj03": ("A, xj 0.03", (7, -10), False)}
    for n, r in rows.items():
        ion, ih = fnum(r["Ion_uA_per_um"]), fnum(r["Ioff398K_pA_per_um"])
        ok = r.get("all_pass", "") in ("True", "true", "1")
        if n == "final":
            ax.plot(ion, ih, "*", ms=18, color="#f2b100", mec="k", zorder=5)
        else:
            ax.plot(ion, ih, "o", ms=7, color=PASS if ok else FAIL, mec="k", mew=0.4)
        if n in lab:
            t, off, leader = lab[n][:3]
            ha = lab[n][3] if len(lab[n]) > 3 else "left"
            ax.annotate(t, (ion, ih), xytext=off, textcoords="offset points", fontsize=7.8, ha=ha,
                        arrowprops=dict(arrowstyle="-", lw=0.6, color="0.35", shrinkB=4) if leader else None)
    ax.set_yscale("log"); ax.set_xlim(100, 720); ax.set_ylim(5, 130)
    ax.set_xlabel("Ion (µA/µm) · Vg = Vd = 2 V"); ax.set_ylabel("고온 Ioff (pA/µm) · 125 °C")
    ax.set_title("Part 1 설계 탐색 (탐색 구조): Ion 과 고온 누설의 맞교환\n초록 = 8항목 전부 통과 · 빨강 = 미달 있음 · 별 = 최종 설계")
    save(fig, "fig04_p1_search.png")


# ---------------------------------------------------------------------------
def fig05_p1_ablation():
    p = R / "ablation_final.csv"
    if not p.exists():
        return
    rows = read_csvout(p)
    cols = [("Vth_V", "Vth"), ("Ion_uA_per_um", "Ion"), ("Ioff300K_pA_per_um", "Ioff 27°C"),
            ("SS_mV_per_dec", "SS"), ("Ioff398K_pA_per_um", "Ioff 125°C"), ("DIBL_mV_per_V", "DIBL"),
            ("BodyEffect_mV", "Body"), ("Eox_MV_per_cm", "E_ox")]
    ref = rows["reference"]
    knobs = [n for n in rows if n != "reference"]
    def pretty(n):
        n = n.replace(" [nm]=", " ").replace(" [um]=", " ").replace(" [cm^-3]=", " ")
        unit = {"tox": " nm", "Lg": " µm", "t_Si": " µm", "xj": " µm"}
        for k_, u in unit.items():
            if n.startswith(k_):
                return n + u
        return n.replace("e+", "e")
    M = np.zeros((len(knobs), len(cols))); T = [["" for _ in cols] for _ in knobs]
    ok = np.ones((len(knobs), len(cols)), bool)
    for i, n in enumerate(knobs):
        for j, (c, _) in enumerate(cols):
            v, v0 = fnum(rows[n][c]), fnum(ref[c])
            ch = (v / v0 - 1) * 100 if v0 else 0
            M[i, j] = np.sign(ch) * np.log10(1 + abs(ch))
            T[i][j] = f"{ch:+.0f}%" if abs(ch) < 1000 else f"×{v/v0:.0f}"
            ok[i, j] = rows[n].get(c + "_pass", "PASS") == "PASS"
    fig, ax = plt.subplots(figsize=(10.5, 0.43 * len(knobs) + 2.0))
    ax.imshow(M, cmap="PuOr_r", vmin=-3.5, vmax=3.5, aspect="auto")
    for i in range(len(knobs)):
        for j in range(len(cols)):
            ax.text(j, i, T[i][j], ha="center", va="center", fontsize=8,
                    color="white" if abs(M[i, j]) > 2.2 else "black", weight=None if ok[i, j] else "bold")
            if not ok[i, j]:
                ax.add_patch(mp.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=FAIL, lw=2.3))
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[1] for c in cols])
    ax.set_yticks(range(len(knobs))); ax.set_yticklabels([pretty(n) for n in knobs])
    ax.grid(False)
    refline = "  ".join(f"{lab} {fnum(ref[c]):.3g}" for c, lab in cols)
    ax.set_title("Ablation (탐색 구조): 최종 소자에서 설계 변수 하나만 바꿨을 때의 변화율\n"
                 f"기준값  {refline}\n주황 = 증가 · 보라 = 감소 · 빨간 테두리 = 그 항목이 기준 미달로 바뀜", fontsize=9.5)
    save(fig, "fig05_p1_ablation.png")


# ---------------------------------------------------------------------------
def fig06_p2_structure():
    Ls = Ld = 0.25; Lg = 0.30; tsi = 1.0; xjm = 0.113
    c, w, t, h = 0.175, 0.040, 0.0045, 1.40
    X = Ls + Lg + Ld
    tl, tr_ = -0.45, -0.30            # p+ 탭 (폭 0.15 µm, 소스와 간격 0.30 µm)
    fig = plt.figure(figsize=(13.4, 6.6))
    a = fig.add_axes([0.02, 0.07, 0.46, 0.86])
    a.add_patch(mp.Rectangle((tl, 0), X - tl, tsi, fc=PSUB, ec="k", lw=0.8))
    a.add_patch(mp.Rectangle((tl, 0), tr_ - tl, xjm, fc=PTAP, ec=TAPC, lw=0.8))
    a.plot([tl, tr_], [-0.003, -0.003], color=TAPC, lw=3.2)
    a.add_patch(mp.Rectangle((0, 0), Ls, xjm, fc=NPLUS, ec="#2e7d32", lw=0.8))
    a.add_patch(mp.Rectangle((Ls + Lg, 0), Ld, xjm, fc=NPLUS, ec="#2e7d32", lw=0.8))
    a.add_patch(mp.Rectangle((Ls, -0.03), Lg, 0.03, fc=OXI, ec="#b71c1c", lw=0.8))
    a.add_patch(mp.Rectangle((Ls, -0.10), Lg, 0.07, fc=GM, ec="#4a148c", lw=0.8))
    td = 0.012
    a.add_patch(mp.Rectangle((c - w / 2 - td, -h), td, h, fc=HK, ec="#b9770e", lw=0.8))
    a.add_patch(mp.Rectangle((c + w / 2, -h), td, h, fc=HK, ec="#b9770e", lw=0.8))
    a.add_patch(mp.Rectangle((c - w / 2, -h), w, h, fc=PIL, ec="#2c3e50", lw=0.8))
    a.plot([0, 0.10], [-0.003, -0.003], color="#333", lw=3.2)
    a.plot([Ls + Lg, X], [-0.003, -0.003], color="#333", lw=3.2)
    ann(a, "저장 커패시터 (C_STORE)\nZrO$_2$ 유전막 두 장 + 금속 기둥\n높이 h = 1.40 µm", (c + w / 2 + td, -1.05),
        (0.42, -1.25), fs=8.5)
    ann(a, "plate 접점 = 1 V (VDD/2)\n(유전막 바깥 면)", (c - w / 2 - td, -1.10), (-0.95, -1.32), fs=8)
    ann(a, "워드라인 WL = 게이트 (Mo)\nLg = 0.30 µm, tox = 5.1 nm", (Ls + Lg / 2, -0.10), (0.48, -0.55), fs=8.5)
    ann(a, "비트라인 BL = drain 접점", (Ls + Lg + Ld * 0.6, -0.003), (0.66, -0.22), fs=8.5)
    ann(a, "source 접점 (0~0.1 µm)", (0.05, -0.003), (-0.95, -0.52), fs=8)
    ann(a, "body 접점 (V_B = −0.5 V)\np+ 탭 위, 간격 0.30 µm", ((tl + tr_) / 2, -0.003), (-0.95, -0.22), fs=8, color=TAPC)
    ann(a, r"n+ 저장 노드 (ND 9$\times10^{20}$)", (0.08, 0.06), (-0.40, 0.36), fs=8)
    ann(a, "p+ 탭 1e20", ((tl + tr_) / 2, 0.06), (-0.95, 0.16), fs=8, color=TAPC)
    ann(a, r"n+ 드레인 (ND 9$\times10^{20}$)", (0.70, 0.06), (0.62, 0.30), fs=8)
    a.text(X / 2, 0.66, "p형 기판\n" + r"NA = 9.5$\times10^{16}$ cm$^{-3}$", ha="center", va="center", fontsize=8.5)
    a.text((X + tl) / 2, tsi + 0.04, "바닥: 접점 없음", ha="center", va="top", fontsize=8)
    dim_v(a, 0.07, -h, 0, "h\n1.40 µm", fs=8, side="left")
    dim_v(a, X + 0.05, 0, tsi, "t_Si\n1.0 µm", fs=8)
    dim_h(a, Ls, Ls + Lg, 0.42, "Lg 0.30 µm", fs=8)
    a.set_xlim(-0.98, 1.18); a.set_ylim(1.17, -1.52)
    a.set_aspect("equal"); a.grid(False); a.set_xlabel("x (µm)"); a.set_ylabel("y (µm)")
    a.set_title("(a) 1T1C 셀 전체 (실제 비율, 유전막 두께만 과장)")
    # 확대
    z = fig.add_axes([0.50, 0.52, 0.21, 0.41])
    z.add_patch(mp.Rectangle((c - w / 2 - t, -0.12), t, 0.12, fc=HK, ec="#b9770e"))
    z.add_patch(mp.Rectangle((c + w / 2, -0.12), t, 0.12, fc=HK, ec="#b9770e"))
    z.add_patch(mp.Rectangle((c - w / 2, -0.12), w, 0.12, fc=PIL, ec="#2c3e50"))
    z.add_patch(mp.Rectangle((c - 0.065, 0), 0.13, 0.04, fc=NPLUS, ec="#2e7d32"))
    for xx, col in ((c - w / 2 - t, "#1f6fb4"), (c - w / 2, C1), (c + w / 2, C1), (c + w / 2 + t, "#1f6fb4")):
        z.plot([xx, xx], [-0.12, 0], color=col, lw=2.4)
    z.text(c - w / 2 - t - 0.003, -0.105, "plate_l", ha="right", fontsize=7.5, color="#1f6fb4")
    z.text(c + w / 2 + t + 0.003, -0.105, "plate_r", ha="left", fontsize=7.5, color="#1f6fb4")
    z.text(c - w / 2 + 0.002, -0.02, "storage_l", ha="left", fontsize=6.5, color=C1, rotation=90)
    z.text(c + w / 2 - 0.002, -0.02, "storage_r", ha="right", fontsize=6.5, color=C1, rotation=90)
    z.text(c, -0.065, "metal_pillar\n(폭 40 nm)", ha="center", va="center", fontsize=7, color="white")
    z.text(c, 0.022, "n+ 소스 = 저장 노드", ha="center", va="center", fontsize=7.5)
    ann(z, "hk_l / hk_r: ZrO$_2$\n두께 4.5 nm", (c + w / 2 + t / 2, -0.09), (c + 0.035, -0.142), fs=7, va="center")
    ann(z, "계면 bulk_hk_l (유전막 바닥)", (c - w / 2 - t / 2, 0.0), (c - 0.062, 0.062), fs=7)
    z.set_xlim(c - 0.065, c + 0.065); z.set_ylim(0.075, -0.16)
    z.set_aspect("equal"); z.grid(False); z.set_xticks([]); z.set_yticks([])
    z.set_title("(b) 커패시터 아래쪽 확대 (실제 비율)", fontsize=9)
    # 설계값 표
    tb = fig.add_axes([0.74, 0.52, 0.25, 0.41]); tb.axis("off")
    txt = ("최종 셀 설계값\n\n"
           "패스 트랜지스터\n"
           "  Lg 0.30 µm (고정) · W 100 nm\n"
           "  tox 5.1 nm · 게이트 Mo (4.70 eV)\n"
           r"  NA 9.5$\times10^{16}$ · ND 9$\times10^{20}$ cm$^{-3}$" "\n"
           "  xj 0.10 µm · t_Si 1.0 µm\n"
           "  S/D 길이 0.25 µm\n"
           "  body 접점: p+ 탭 위 (탭 0.15 µm)\n\n"
           "저장 커패시터\n"
           "  ZrO$_2$ (εr 35) 두 장\n"
           "  높이 1.40 µm · 두께 4.5 nm\n"
           "  기둥 폭 40 nm\n"
           "  C_STORE = 19.24 fF (≤ 20 fF)")
    tb.text(0.02, 0.98, txt, va="top", ha="left", fontsize=8.6, linespacing=1.45,
            bbox=dict(fc="#fbfbfb", ec="#999", boxstyle="round,pad=0.6"))
    m = fig.add_axes([0.53, 0.07, 0.44, 0.38])
    cs = doping_map(m, D / "struct_p2_tap.csv", "(c) 실제 메시의 순 도핑 (p+ 탭 ~ 드레인, 보라 선 = p+ 경계)", (-0.45, 0.8), (0.35, -0.01), tap=True)
    save(fig, "fig06_p2_structure.png")


# ---------------------------------------------------------------------------
def fig07_p2_read():
    p = R / "p2_read_curves_G2.csv"
    if not p.exists():
        return
    by = {}
    with open(p, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            by.setdefault(r["state"], []).append(r)
    cs_fF = json.load(open(R / "bodytap" / "p2_G2_bodytap.json"))["bodytap"]["c_store_fF"]
    lim = 1000.0 * cs_fF / (cs_fF + 100.0)
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(14.0, 4.4), gridspec_kw=dict(width_ratios=[1.35, 1, 1]))
    for st, col in (('data "0"', C0), ('data "1"', C1)):
        rs = by.get(st, [])
        t = np.array([float(r["t_ps"]) for r in rs]); dv = np.array([float(r["dV_BL_mV"]) for r in rs])
        vc = np.array([float(r["V_cell_V"]) for r in rs])
        ii = np.abs(np.array([float(r["I_source_uA_per_um"]) for r in rs])) * 0.1
        is1 = "1" in st
        lab = "저장 '1' (시작 V_cell = 2 V)" if is1 else "저장 '0' (시작 V_cell = 0 V)"
        a.plot(t, dv, "-o", ms=2.0, color=col, label=lab)
        a.plot([500], [dv[-1]], "s", ms=7, color=col, mec="k")
        frac = abs(dv[-1]) / lim * 100
        a.annotate(f"0.5 ns: {dv[-1]:+.2f} mV\n(한계의 {frac:.0f}%)", (500, dv[-1]),
                   xytext=(330, dv[-1] + (28 if is1 else -14)), fontsize=8.5, color=col, weight="bold",
                   va="center", arrowprops=dict(arrowstyle="->", lw=0.7, color=col))
        b.plot(t, ii, "-", color=col, label=lab)
        c.plot(t, vc, "-", color=col, label=lab)
    for y in (60, -60):
        a.axhline(y, color="k", ls="--", lw=0.9)
    for y in (lim, -lim):
        a.axhline(y, color="gray", ls=":", lw=1)
    a.text(8, 63, "판정선 +60 mV", fontsize=7.8); a.text(8, -57, "판정선 −60 mV", fontsize=7.8, va="bottom")
    a.text(8, lim - 3, f"전하 분배 한계 ±{lim:.0f} mV (시간이 무한히 흐르면)", fontsize=7.5, color="gray", va="top")
    a.axvline(500, color="gray", ls=":", lw=1)
    a.set_ylim(-lim - 20, lim + 20)
    a.set_xlabel("WL 이 켜진 뒤 시간 (ps)"); a.set_ylabel("ΔV_BL = V_BL − 1 V (mV)")
    a.set_title("(a) 비트라인 전압 변화 · 판정: 0.5 ns 에서 |ΔV_BL| ≥ 60 mV")
    fig.suptitle("최종 셀 G2 · 제출 구조 (p+ 탭 body 접점) · 125 °C, V_B = −0.5 V", fontsize=10.5, y=1.02)
    a.legend(loc="center right", fontsize=8)
    b.set_xlabel("시간 (ps)"); b.set_ylabel("|I| (µA), W = 100 nm")
    b.set_title("(b) 패스 트랜지스터 전류"); b.legend(fontsize=8)
    c.axhline(1.0, color="gray", ls=":", lw=1); c.text(5, 1.03, "BL 예비충전 1 V", fontsize=7.5, color="gray")
    c.set_xlabel("시간 (ps)"); c.set_ylabel("V_cell (V)")
    c.set_title("(c) 저장 노드 전압"); c.legend(fontsize=8, loc="center right")
    save(fig, "fig07_p2_read.png")


# ---------------------------------------------------------------------------
def fig08_p2_retention():
    p = D / "retention_G2_tap.json"
    if not p.exists():
        return
    d = json.load(open(p)); g2 = json.load(open(R / "bodytap" / "p2_G2_bodytap.json"))["bodytap"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.8, 4.4))
    r1 = d["1"]
    t1 = np.array(r1["t"][1:], float) * 1e3; v1 = np.array(r1["v"][1:], float)
    vcrit1, tr = float(r1["vcrit"]), float(r1["t_ret"]) * 1e3
    a.semilogx(t1, v1, "-", color=C1, lw=1.9, label="저장 '1' (유지 중 BL = 0 V)")
    a.axhline(vcrit1, color="k", ls="--", lw=0.9)
    a.text(t1.min() * 1.3, vcrit1 + 0.006, f"V_crit = {vcrit1:.3f} V : 이보다 낮으면 읽기 마진 < 60 mV",
           fontsize=8, va="bottom")
    a.axvline(64, color="gray", ls=":", lw=1.2); a.text(60, 1.99, "64 ms", fontsize=8.5, ha="right", va="top")
    a.axvline(tr, color=C1, ls="-.", lw=1)
    a.text(tr * 0.93, 1.94, f"t_ret = {tr:.0f} ms", fontsize=9, color=C1, weight="bold", ha="right")
    a.plot([64], [g2["v64_1_V"]], "o", color="k", zorder=5)
    ann(a, f"64 ms 지점: V_cell = {g2['v64_1_V']:.3f} V\n여기서 다시 읽으면 {g2['read64_1_mV']:.1f} mV (통과)",
        (64, g2["v64_1_V"]), (t1.min() * 1.3, 1.80), fs=8.5)
    a.set_ylim(vcrit1 - 0.03, 2.01)
    a.set_xlabel("유지 시간 (ms, 로그 축)"); a.set_ylabel("V_cell (V)")
    a.set_title("(a) 저장 '1' : 누설로 전압이 내려간다 (WL 0 V, 125 °C)")
    a.legend(loc="lower left", fontsize=8)
    r0 = d["0"]
    t0 = np.array(r0["t"][1:], float) * 1e3; v0 = np.array(r0["v"][1:], float)
    vcrit0 = float(r0["vcrit"])
    b.semilogx(t0, v0 * 1e3, "-", color=C0, lw=1.9, label="저장 '0' (유지 중 BL = 2 V)")
    v10 = float(np.interp(1.0e4, t0, v0))
    b.plot([1.0e4], [v10 * 1e3], "s", color=C0, ms=5, zorder=5)
    b.text(t0[-1] * 0.8, v0[-1] * 1e3 + 10,
           f"10 s 에 {v10*1e3:.0f} mV (적분이 멈춘 {t0[-1]/1e3:.0f} s 에 {v0[-1]*1e3:.0f} mV)\n계속 아주 천천히 오름",
           fontsize=8.0, color=C0, va="bottom", ha="right")
    b.axvline(64, color="gray", ls=":", lw=1.2)
    b.plot([64], [g2["v64_0_V"] * 1e3], "o", color="k", zorder=5)
    ann(b, f"64 ms: V_cell = {g2['v64_0_V']*1e3:.0f} mV\n다시 읽으면 {g2['read64_0_mV']:.1f} mV",
        (64, g2["v64_0_V"] * 1e3), (t0.min() * 1.3, 300), fs=8.3)
    b.axhline(vcrit0 * 1e3, color="k", ls="--", lw=0.9)
    b.text(t0.min() * 1.3, vcrit0 * 1e3 + 8, f"V_crit = {vcrit0*1e3:.0f} mV : 여기까지 올라가야 실패",
           fontsize=8, va="bottom")
    b.set_ylim(-20, vcrit0 * 1e3 * 1.18)
    b.set_xlabel("유지 시간 (ms, 로그 축)"); b.set_ylabel("V_cell (mV)")
    b.set_title("(b) 저장 '0' : 아주 천천히 오르고 V_crit 에 닿지 않는다")
    b.legend(loc="upper left", bbox_to_anchor=(0.0, 0.82), fontsize=8)
    save(fig, "fig08_p2_retention.png")


# ---------------------------------------------------------------------------
def fig09_p2_tradeoff():
    pts = []
    for p in sorted(R.glob("p2_*.json")):
        k = p.stem[3:]
        if k in ("A", "B", "start", "G1") or k.startswith("G1"):
            continue
        r = json.load(open(p, encoding="utf-8"))
        if "read64_1_mV" in r:
            pts.append((k, r))
    if not pts:
        return
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.2, 5.0))
    name = {"G2": "최종 G2", "P2start": "Part 1 소자 그대로", "P2oldA": "이전(9월) 설계",
            "C2": "C2 (TiN, tox 5.5)", "E3": "E3 (W)", "C1": "C1 (TiN, tox 5.0)"}
    offa = {"G2": (10, 10), "P2start": (6, -12), "P2oldA": (6, 6), "C2": (6, -12), "E3": (-10, 10), "C1": (-60, 10)}
    offb = {"G2": (10, 8), "P2start": (6, -12), "P2oldA": (6, 6), "C2": (6, -12), "E3": (6, -14), "C1": (6, -12)}
    for k, r in pts:
        g = r["device"]["gate_metal"]
        tr = r.get("tret_1_ms", 0.1)
        tr = 0.1 if (not isinstance(tr, (int, float)) or tr != tr or tr < 0.1) else tr
        col = GATE_COL.get(g, "gray")
        mk, ms = ("*", 18) if k == "G2" else ("o", 7)
        early = tr < 64 and r["read_1_mV"] >= 60       # 64 ms 전에 V_crit 도달: 64 ms 값은 60 mV 미만
        y64 = 58.6 if early else r["read64_1_mV"]
        if early:
            a.plot(r["read_1_mV"], y64, "v", color=col, ms=9, mfc="white", mec=col, mew=1.6, zorder=3)
            b.plot(tr, y64, "v", color=col, ms=9, mfc="white", mec=col, mew=1.6, zorder=3)
        else:
            a.plot(r["read_1_mV"], y64, mk, color=col, ms=ms, mec="k", mew=0.5, zorder=4 if k == "G2" else 3)
            b.plot(tr, y64, mk, color=col, ms=ms, mec="k", mew=0.5, zorder=4 if k == "G2" else 3)
        r = dict(r); r["read64_1_mV"] = y64
        if k in name:
            a.annotate(name[k], (r["read_1_mV"], r["read64_1_mV"]), xytext=offa[k], textcoords="offset points",
                       fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.5, color="gray"))
            b.annotate(name[k], (tr, r["read64_1_mV"]), xytext=offb[k], textcoords="offset points",
                       fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.5, color="gray"))
    for g, col in GATE_COL.items():
        a.plot([], [], "o", color=col, mec="k", mew=0.5, label=f"게이트 {g}")
    a.axhline(60, color="k", ls="--", lw=0.8); a.axvline(60, color="k", ls="--", lw=0.8)
    lo = min(min(r["read_1_mV"] for _, r in pts), 50)
    a.plot([lo, 92], [lo, 92], ":", color="gray", lw=0.9)
    a.text(84, 86.5, "대각선 = 64 ms 동안 손실 없음", fontsize=7.5, color="gray", rotation=33, ha="center")
    a.plot([], [], "v", color="gray", mfc="white", mec="gray", mew=1.4, ms=8, label="64 ms 전에 60 mV 아래로 (값은 그보다 낮음)")
    a.set_xlabel("갓 저장한 셀의 읽기 '1' 마진 (mV)"); a.set_ylabel("64 ms 유지 후 읽기 '1' 마진 (mV)")
    a.set_title("(a) 처음 읽기가 좋다고 64 ms 뒤에도 좋은 것은 아니다"); a.legend(fontsize=8, loc="upper left")
    b.set_xscale("log"); b.axhline(60, color="k", ls="--", lw=0.8); b.axvline(64, color="k", ls="--", lw=0.8)
    b.axvspan(64, 3000, color=PASS, alpha=0.06)
    b.text(70, b.get_ylim()[0] + 0.3, "통과 영역", fontsize=8, color=PASS)
    b.set_xlabel("리텐션 시간 t_ret '1' (ms, 로그 축)"); b.set_ylabel("64 ms 유지 후 읽기 '1' 마진 (mV)")
    b.set_title("(b) 후보 전체: 통과 영역은 오른쪽 위"); b.set_xlim(0.08, 3000)
    save(fig, "fig09_p2_tradeoff.png")


# ---------------------------------------------------------------------------
def fig10_p2_cap():
    p = D / "verify_cap.json"
    if not p.exists():
        return
    d = json.load(open(p, encoding="utf-8"))
    mats = [m for m in ("SiO2", "Al2O3", "HfO2", "ZrO2") if m in d]
    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    x = np.arange(len(mats))
    cm = [d[m]["c_meas"] for m in mats]; cc = [d[m]["c_closed"] for m in mats]
    ax.bar(x - 0.2, cc, 0.4, color="#bbbbbb", label=r"해석식 C = $\varepsilon_r \varepsilon_0 \cdot 2h / t \times W$")
    ax.bar(x + 0.2, cm, 0.4, color="#f39c12", label="시뮬레이션: plate 전하의 dQ/dV")
    for xi, a_, b_ in zip(x, cc, cm):
        ax.text(xi + 0.2, b_ + 0.3, f"{b_:.2f}\n({(b_/a_-1)*100:+.2f}%)", ha="center", fontsize=8)
    lab = {"SiO2": "SiO$_2$", "Al2O3": "Al$_2$O$_3$", "HfO2": "HfO$_2$", "ZrO2": "ZrO$_2$"}
    ax.set_xticks(x); ax.set_xticklabels([f"{lab[m]}\n(εr {d[m]['eps']:g})" for m in mats])
    ax.set_ylabel("C_STORE (fF), W = 100 nm"); ax.set_ylim(0, max(cc) * 1.32)
    ax.set_title("C_STORE 추출 검증 · 과제 예시 구조 (h 1.0 µm, t 5 nm, 기둥 40 nm)")
    ax.legend(loc="upper left", fontsize=8)
    save(fig, "fig10_p2_cap.png")


if __name__ == "__main__":
    for f in (fig01_p1_structure, fig02_p1_idvg, fig03_p1_physics, fig04_p1_search, fig05_p1_ablation,
              fig06_p2_structure, fig07_p2_read, fig08_p2_retention, fig09_p2_tradeoff, fig10_p2_cap):
        try:
            f()
        except Exception as e:
            import traceback
            print("FAILED", f.__name__, type(e).__name__, e); traceback.print_exc()


# ---------------------------------------------------------------------------
def fig11_bodytap():
    """body 접점 위치 비교: 바닥 접점(탐색 구조) 대 p+ 탭(제출 구조)."""
    p2 = R / "bodytap" / "p2_G2_bodytap.json"
    if not p2.exists():
        print("skip fig11: no", p2); return
    d2 = json.load(open(p2, encoding="utf-8"))
    rb = json.load(open(D / "retention_G2.json", encoding="utf-8"))["1"]
    rt = d2["retention_bodytap"]["1"]
    fig = plt.figure(figsize=(13.6, 4.7))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 0.82], wspace=0.34)
    a, b, c = (fig.add_subplot(gs[0, k]) for k in range(3))

    # (a) 탭 구조의 실제 도핑 + p+ 경계(보라 등고선)
    xs, ys, ns = [], [], []
    with open(D / "struct_p2_tap.csv") as fh:
        for r in csv.DictReader(fh):
            if r["region"] == "bulk":
                xs.append(float(r["x_um"])); ys.append(float(r["y_um"])); ns.append(float(r["NetDoping"]))
    xs, ys, ns = map(np.asarray, (xs, ys, ns))
    v = np.sign(ns) * np.log10(np.maximum(np.abs(ns), 1e14)); v = np.where(np.abs(v) < 14.0, 0.0, v)
    tri = mtri.Triangulation(xs, ys)
    cs = a.tricontourf(tri, v, levels=np.linspace(-21, 21, 43), cmap="RdBu_r", extend="both")
    a.tricontour(tri, ns, levels=[0.0], colors="k", linewidths=0.9)
    a.tricontour(tri, ns, levels=[-1e19], colors="#6a1b9a", linewidths=1.3)
    a.set_aspect("equal"); a.set_xlim(-0.47, 0.82); a.set_ylim(1.02, -0.26); a.grid(False)
    a.set_xlabel("x (µm)"); a.set_ylabel("y (µm)")
    a.set_title("(a) p+ 탭 구조: 최종 셀의 실제 메시 순 도핑", pad=6)
    cb = fig.colorbar(cs, ax=a, shrink=0.8, ticks=[-20, -17, 0, 18, 20], pad=0.02)
    cb.ax.set_yticklabels([r"p $10^{20}$", r"p $10^{17}$", "0", r"n $10^{18}$", r"n $10^{20}$"], fontsize=7.5)
    for (x0, x1), lab, col in (((-0.45, -0.30), "body", "#6a1b9a"), ((0.0, 0.10), "source", "#333"),
                               ((0.55, 0.80), "drain", "#333")):
        a.plot([x0, x1], [-0.012, -0.012], color=col, lw=3.2, solid_capstyle="butt")
        a.text((x0 + x1) / 2, -0.03, lab, ha="center", va="bottom", fontsize=7.8, color=col)
    a.add_patch(mp.Rectangle((0.1505, -0.24), 0.049, 0.23, fc=HK, ec="k", lw=0.6))
    a.text(0.21, -0.20, "커패시터\n(위로 1.4 µm)", ha="left", va="center", fontsize=7.2)
    a.add_patch(mp.Rectangle((0.25, -0.05), 0.30, 0.04, fc=GM, ec="k", lw=0.6))
    a.text(0.40, -0.055, "게이트 (Mo)", ha="center", va="bottom", fontsize=7.4)
    ann(a, "p+ 탭 1e20 (폭 0.15 µm)\n보라 선 = 1e19 경계", (-0.37, 0.07), (-0.44, 0.40), fs=7.4, color="#6a1b9a")
    a.annotate("", xy=(-0.30, -0.075), xytext=(0.0, -0.075), arrowprops=dict(arrowstyle="<->", lw=0.9))
    a.text(-0.15, -0.09, "간격 0.30 µm\n(p 기판)", ha="center", va="bottom", fontsize=7.2)
    a.text(0.18, 0.96, "바닥: 접점 없음 (전류가 나가지 않는 면)", ha="center", va="bottom", fontsize=7.6,
           bbox=dict(fc="white", ec="none", alpha=0.8))

    # (b) '1' 유지 곡선
    tb, vb = np.asarray(rb["t"]) * 1e3, np.asarray(rb["v"])
    tt, vt = np.asarray(rt["t"]) * 1e3, np.asarray(rt["v"])
    b.plot(tb, vb, color="0.45", lw=1.8, ls="--", label=f"바닥 접점 (탐색 구조): t_ret {rb['t_ret']*1e3:.1f} ms")
    b.plot(tt, vt, color=C1, lw=1.8, label=f"p+ 탭 (제출 구조): t_ret {rt['t_ret']*1e3:.1f} ms")
    b.plot([tb[-1]], [vb[-1]], "o", color="0.45", ms=5); b.plot([tt[-1]], [vt[-1]], "o", color=C1, ms=5)
    vc = min(rb["vcrit"], rt["vcrit"])
    b.axhline(rb["vcrit"], color="0.3", lw=0.8, ls=":")
    b.text(95, vc - 0.004, f"V_crit {rb['vcrit']:.3f} V (두 구조 같음)", fontsize=7.6, va="top")
    b.axvline(64, color="k", lw=0.8, ls="-.")
    b.text(66, 1.76, "64 ms", fontsize=7.8, va="bottom")
    b.set_xlabel("유지 시간 (ms)"); b.set_ylabel("셀 전압 V_cell (V)")
    b.set_title("(b) '1' 유지 (125 °C, BL 0 V, V_B −0.5 V)")
    b.set_xlim(0, max(tb[-1], tt[-1]) * 1.05); b.set_ylim(vc - 0.03, 2.01)
    b.legend(loc="upper right", fontsize=7.6)

    # (c) 저장 노드 누설이 어느 단자로 나가나 (V_cell = 2 V)
    sb, st = d2.get("split_bottom", {}), d2.get("split_bodytap", {})
    def pick(s):
        k = min(s, key=lambda kk: abs(float(kk) - 2.0))
        return abs(s[k]["drain"]), abs(s[k]["body"])
    rows = [("바닥 접점\n(탐색 구조)", *pick(sb)), ("p+ 탭\n(제출 구조)", *pick(st))]
    x = np.arange(len(rows))
    dr = np.array([r[1] for r in rows]); bo = np.array([r[2] for r in rows])
    c.bar(x, dr, 0.55, color="#9ecae1", ec="k", lw=0.5, label="BL(드레인) 단자로")
    c.bar(x, bo, 0.55, bottom=dr, color="#fdae6b", ec="k", lw=0.5, label="body 단자로")
    for xi, d_, b_ in zip(x, dr, bo):
        c.text(xi, d_ + b_ + 0.006, f"합 {d_ + b_:.4f}\n(body {b_ / (d_ + b_) * 100:.1f}%)", ha="center",
               va="bottom", fontsize=7.6)
    c.set_xticks(x); c.set_xticklabels([r[0] for r in rows], fontsize=8)
    c.set_ylabel("저장 노드 누설 (pA/µm)"); c.set_ylim(0, max(dr + bo) * 1.3)
    c.set_title("(c) '1' 유지 누설이 나가는 단자 (V_cell 2 V)")
    c.legend(loc="lower center", fontsize=7.3, framealpha=0.9)
    save(fig, "fig11_p2_bodytap.png")
