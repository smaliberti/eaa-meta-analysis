"""
step7_meta_analysis.py
======================
Meta-analysis: Modifiable determinants of epigenetic age acceleration
Step 7: Random-effects model (DerSimonian-Laird) + Forest plots

STATISTICAL APPROACH:
- Estimator: DerSimonian-Laird (DL) tau²
- Weights: inverse-variance (wi = 1 / (sei² + tau²))
- Pooled estimate: weighted mean of yi
- Heterogeneity: Cochran Q, I², tau²
- Confidence intervals: 95% (z = 1.96)
- Analyses: (1) overall across all studies
             (2) subgroup by exposure category (≥3 studies)

OUTPUT FILES:
- step7_results.xlsx
- forest_plot_overall_part1.pdf, part2.pdf, ...
- forest_plot_[category].pdf
- forest_plot_summary_subgroups.pdf
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

warnings.filterwarnings('ignore')

# ── 0. Load data ────────────────────────────────────────────────────────────
INPUT = "STEP6_.xlsx"
OUTDIR = "step7_outputs"
os.makedirs(OUTDIR, exist_ok=True)

df = pd.read_excel(INPUT, sheet_name='pool_A_primary_beta')
df['yi'] = pd.to_numeric(df['yi'], errors='coerce')
df['sei'] = pd.to_numeric(df['sei'], errors='coerce')
df = df.dropna(subset=['yi', 'sei'])
df = df[df['sei'] > 0].reset_index(drop=True)

print(f"Pool A loaded: {len(df)} studies | yi range: [{df['yi'].min():.3f}, {df['yi'].max():.3f}]")

# ── 1. Exposure category harmonisation ──────────────────────────────────────
def harmonize(e):
    if pd.isna(e):
        return 'Other'
    e = str(e).lower()

    if 'smok' in e:
        return 'Smoking'

    if any(x in e for x in ['physical activity', 'exercise', 'sedent', 'walking', 'fitness', 'gardening']):
        return 'Physical activity'

    if any(x in e for x in [
        'diet', 'nutriti', 'food', 'vitamin', 'fiber', 'fatty', 'flavon',
        'iron', 'sugar', 'alcohol', 'drink', 'cannabis', 'marijuana'
    ]):
        return 'Diet & Nutrition / Alcohol'

    if 'sleep' in e:
        return 'Sleep'

    if any(x in e for x in [
        'stress', 'psycho', 'anxiety', 'depress', 'mental', 'adversit',
        'trauma', 'childhood', 'personality'
    ]):
        return 'Psychosocial stress'

    if any(x in e for x in [
        'pollut', 'air', 'environmental', 'toxin', 'chemical', 'pesticide',
        'toluene', 'greenness', 'temperature', 'climate', 'urban'
    ]):
        return 'Environmental exposure'

    if any(x in e for x in [
        'socioeconom', 'ses', 'income', 'education', 'occupat', 'neighborhood',
        'social', 'deprivat', 'sdh', 'mobility'
    ]):
        return 'Socioeconomic / Social'

    if any(x in e for x in [
        'bmi', 'obes', 'metabol', 'inflam', 'immune', 'biomark', 'insulin',
        'lipid', 'body', 'syndrome', 'sii'
    ]):
        return 'Metabolic / Inflammatory'

    if any(x in e for x in ['lifestyle', 'physiolog', 'supplement']):
        return 'Lifestyle (mixed)'

    return 'Other'

df['category'] = df['exposure_category'].apply(harmonize)

# ── 2. DerSimonian-Laird random-effects model ───────────────────────────────
def dl_meta(yi, sei):
    """
    DerSimonian-Laird random-effects meta-analysis.
    Returns dict with all statistics.
    """
    yi = np.array(yi, dtype=float)
    sei = np.array(sei, dtype=float)
    vi = sei ** 2
    k = len(yi)

    # Fixed-effects weights and estimate
    w_fe = 1.0 / vi
    theta_fe = np.sum(w_fe * yi) / np.sum(w_fe)

    # Cochran Q
    Q = np.sum(w_fe * (yi - theta_fe) ** 2)
    df_q = k - 1

    # tau² (DL estimator — floored at 0)
    C = np.sum(w_fe) - np.sum(w_fe ** 2) / np.sum(w_fe)
    tau2 = max(0.0, (Q - df_q) / C) if C > 0 else 0.0

    # Random-effects weights
    w_re = 1.0 / (vi + tau2)
    theta_re = np.sum(w_re * yi) / np.sum(w_re)
    se_re = np.sqrt(1.0 / np.sum(w_re))

    # 95% CI
    ci_lo = theta_re - 1.96 * se_re
    ci_hi = theta_re + 1.96 * se_re

    # I²
    I2 = max(0.0, (Q - df_q) / Q * 100) if Q > 0 else 0.0

    # p-values
    z = theta_re / se_re
    from scipy.stats import norm, chi2
    p_theta = 2 * (1 - norm.cdf(abs(z)))
    p_Q = 1 - chi2.cdf(Q, df_q)

    return {
        'k': k,
        'theta': theta_re,
        'se': se_re,
        'ci_lo': ci_lo,
        'ci_hi': ci_hi,
        'tau2': tau2,
        'tau': np.sqrt(tau2),
        'Q': Q,
        'df_Q': df_q,
        'p_Q': p_Q,
        'I2': I2,
        'z': z,
        'p_theta': p_theta,
        'w_re': w_re,
        'theta_fe': theta_fe,
    }

# ── 3. Graphic helpers ──────────────────────────────────────────────────────
COLORS = {
    'diamond':  '#185FA5',
    'square':   '#378ADD',
    'ci_line':  '#444441',
    'zeroline': '#E24B4A',
    'grid':     '#D3D1C7',
    'text':     '#2C2C2A',
    'header':   '#2E4C6D',
    'I2_low':   '#3B6D11',
    'I2_mid':   '#BA7517',
    'I2_high':  '#A32D2D',
    'nonsig':   '#888780',
}

def i2_color(i2):
    if i2 < 25:
        return COLORS['I2_low']
    if i2 < 75:
        return COLORS['I2_mid']
    return COLORS['I2_high']

def short_author(name, max_len=18):
    name = str(name).strip()
    return name if len(name) <= max_len else name[:max_len - 1] + "…"

def safe_year(x):
    return str(int(x)) if pd.notna(x) else ""

def robust_xlim(lo_values, hi_values, pooled_lo=None, pooled_hi=None):
    """
    Robust x-limits to avoid extreme outliers making the plot unreadable.
    """
    lo_values = np.array(lo_values, dtype=float)
    hi_values = np.array(hi_values, dtype=float)

    x_lo = np.percentile(lo_values, 5)
    x_hi = np.percentile(hi_values, 95)

    if pooled_lo is not None:
        x_lo = min(x_lo, pooled_lo)
    if pooled_hi is not None:
        x_hi = max(x_hi, pooled_hi)

    if np.isclose(x_lo, x_hi):
        x_lo -= 1.0
        x_hi += 1.0

    xpad = (x_hi - x_lo) * 0.10
    return x_lo - xpad, x_hi + xpad

# ── 4. Simplified readable forest plot ──────────────────────────────────────
def forest_plot(data, meta, title, filename, show_summary=True):
    """
    Readable forest plot:
    - only First Author + Year on the left
    - forest graphic on the right
    - pooled subtitle on top
    - no long text columns inside plot
    """
    data = data.copy().sort_values('yi').reset_index(drop=True)

    # ricalcolo meta per avere pesi allineati all'ordine del plot
    meta_plot = dl_meta(data['yi'], data['sei'])
    w_re = meta_plot['w_re']
    rel_w = w_re / w_re.max() if len(w_re) > 0 else np.ones(len(data))

    k = len(data)
    fig_h = max(8, 1.8 + k * 0.42)
    fig, ax = plt.subplots(figsize=(12, fig_h))
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')

    yi_arr = data['yi'].values
    sei_arr = data['sei'].values
    all_lo = yi_arr - 1.96 * sei_arr
    all_hi = yi_arr + 1.96 * sei_arr

    # per evitare che pochi outlier schiaccino tutto
    pooled_lo = meta['ci_lo'] if show_summary else None
    pooled_hi = meta['ci_hi'] if show_summary else None
    xmin, xmax = robust_xlim(all_lo, all_hi, pooled_lo=pooled_lo, pooled_hi=pooled_hi)

    # y positions
    row_step = 1.15
    y_studies = [row_step * (k - i) for i in range(k)]
    y_summary = -0.9
    y_header = row_step * k + 1.1

    ymin = y_summary - 0.8 if show_summary else -0.4
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, y_header + 0.8)

    # zero line
    ax.axvline(0, color=COLORS['zeroline'], lw=1.1, ls='--', alpha=0.7, zorder=1)

    # horizontal grid
    for y in y_studies:
        ax.axhline(y, color=COLORS['grid'], lw=0.4, alpha=0.5, zorder=0)

    # left text columns
    COL_AUTH = -0.24
    COL_YEAR = -0.08

    def ax_text(txt, ax_x, y_data, ha='left', fontsize=8, bold=False, color=COLORS['text']):
        ax.text(
            ax_x, y_data, txt,
            transform=ax.get_yaxis_transform(),
            ha=ha, va='center',
            fontsize=fontsize,
            fontweight='bold' if bold else 'normal',
            color=color
        )

    # header
    ax_text("First Author", COL_AUTH, y_header, bold=True, fontsize=9.5, color=COLORS['header'])
    ax_text("Year",         COL_YEAR, y_header, bold=True, fontsize=9.5, color=COLORS['header'])
    ax.axhline(y_header - 0.45, color=COLORS['header'], lw=0.9, alpha=0.6)

    # studies
    for i, (_, row) in enumerate(data.iterrows()):
        y = y_studies[i]
        yi_i = row['yi']
        sei_i = row['sei']
        lo_i = yi_i - 1.96 * sei_i
        hi_i = yi_i + 1.96 * sei_i

        lo_plot = max(lo_i, xmin)
        hi_plot = min(hi_i, xmax)

        # CI line
        ax.plot([lo_plot, hi_plot], [y, y], color=COLORS['ci_line'], lw=1.0, zorder=2)

        # truncated CI arrows
        if lo_i < xmin:
            ax.plot(xmin, y, marker='<', color=COLORS['ci_line'], markersize=5, zorder=2)
        if hi_i > xmax:
            ax.plot(xmax, y, marker='>', color=COLORS['ci_line'], markersize=5, zorder=2)

        # square marker
        marker_size = 28 + 80 * rel_w[i]
        ax.scatter(
            yi_i, y,
            s=marker_size,
            marker='s',
            color=COLORS['square'],
            edgecolor=COLORS['ci_line'],
            linewidth=0.6,
            zorder=3
        )

        auth = short_author(row.get('first_author', ''), 18)
        year = safe_year(row.get('year', np.nan))

        ax_text(auth, COL_AUTH, y, fontsize=8.2)
        ax_text(year, COL_YEAR, y, fontsize=8.2)

    # pooled diamond
    if show_summary:
        diam_h = 0.34
        diamond = plt.Polygon(
            [[meta['ci_lo'], y_summary],
             [meta['theta'], y_summary + diam_h],
             [meta['ci_hi'], y_summary],
             [meta['theta'], y_summary - diam_h]],
            closed=True,
            facecolor=COLORS['diamond'],
            edgecolor='white',
            lw=0.8,
            zorder=4
        )
        ax.add_patch(diamond)
        ax_text("Pooled effect", COL_AUTH, y_summary, bold=True, fontsize=8.8, color=COLORS['diamond'])

    # axes
    ax.set_xlabel("Effect size β (95% CI)", fontsize=9.5, color=COLORS['text'], labelpad=10)
    ax.set_yticks([])
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.spines['bottom'].set_color(COLORS['grid'])
    ax.tick_params(axis='x', labelsize=8.5, colors=COLORS['text'])

    # title + subtitle
    ax.set_title(title, fontsize=10.5, fontweight='bold', color=COLORS['header'], loc='left', pad=14)

    subtitle = (
        f"Pooled β = {meta['theta']:+.3f} "
        f"(95% CI {meta['ci_lo']:+.3f}, {meta['ci_hi']:+.3f})   |   "
        f"I² = {meta['I2']:.1f}%   |   τ² = {meta['tau2']:.4f}   |   k = {meta['k']}"
    )
    fig.text(
        0.125, 0.965, subtitle,
        ha='left', va='top',
        fontsize=8.5,
        color=i2_color(meta['I2'])
    )

    plt.subplots_adjust(left=0.30, right=0.97, top=0.93, bottom=0.08)

    out = os.path.join(OUTDIR, filename)
    plt.savefig(out, dpi=220, facecolor='white')
    plt.close()
    print(f"  Saved: {out}")
    return out

# ── 5. Paginated overall forest plot ────────────────────────────────────────
def forest_plot_paginated(data, meta, title_prefix, filename_prefix, chunk_size=25):
    """
    Split a large forest plot into multiple readable pages.
    Final page includes pooled summary diamond.
    """
    data = data.copy().sort_values('yi').reset_index(drop=True)
    n = len(data)
    files = []

    for part_idx, start in enumerate(range(0, n, chunk_size), start=1):
        end = min(start + chunk_size, n)
        chunk = data.iloc[start:end].copy().reset_index(drop=True)
        is_last = (end == n)

        title = f"{title_prefix}\nStudies {start + 1}–{end} of {n}"
        filename = f"{filename_prefix}_part{part_idx}.pdf"

        forest_plot(
            chunk,
            meta=meta,
            title=title,
            filename=filename,
            show_summary=is_last
        )
        files.append(filename)

    return files

# ── 6. Run overall model ────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("OVERALL MODEL — Pool A")
print("=" * 60)

meta_overall = dl_meta(df['yi'], df['sei'])
print(f"  k          = {meta_overall['k']}")
print(f"  Pooled β   = {meta_overall['theta']:+.4f}")
print(f"  95% CI     = [{meta_overall['ci_lo']:+.4f}, {meta_overall['ci_hi']:+.4f}]")
print(f"  SE         = {meta_overall['se']:.4f}")
print(f"  z          = {meta_overall['z']:.3f},  p = {meta_overall['p_theta']:.4f}")
print(f"  tau²       = {meta_overall['tau2']:.4f}")
print(f"  I²         = {meta_overall['I2']:.1f}%")
print(f"  Q({meta_overall['df_Q']})    = {meta_overall['Q']:.2f},  p_Q = {meta_overall['p_Q']:.4f}")

overall_files = forest_plot_paginated(
    df,
    meta_overall,
    title_prefix="Forest Plot — Overall random-effects model | Outcome: Epigenetic Age Acceleration (years)",
    filename_prefix="forest_plot_overall",
    chunk_size=25
)

# ── 7. Subgroup analyses by exposure category ───────────────────────────────
print("\n" + "=" * 60)
print("SUBGROUP ANALYSES BY EXPOSURE CATEGORY")
print("=" * 60)

cat_counts = df['category'].value_counts()
cats_poolable = cat_counts[cat_counts >= 3].index.tolist()
cats_narrative = cat_counts[cat_counts < 3].index.tolist()

print(f"\nCategories with ≥3 studies (poolable): {len(cats_poolable)}")
print(f"Categories with <3 studies (narrative only): {cats_narrative}")

subgroup_results = []

for cat in sorted(cats_poolable):
    sub = df[df['category'] == cat].copy().reset_index(drop=True)
    meta = dl_meta(sub['yi'], sub['sei'])

    print(f"\n  [{cat}]  k={meta['k']}")
    print(f"    β = {meta['theta']:+.4f}  [{meta['ci_lo']:+.4f}, {meta['ci_hi']:+.4f}]")
    print(f"    I² = {meta['I2']:.1f}%  |  Q({meta['df_Q']})={meta['Q']:.2f} p={meta['p_Q']:.3f}")

    fname = f"forest_plot_{cat.lower().replace(' ', '_').replace('/', '_')}.pdf"
    forest_plot(
        sub,
        meta,
        title=f"Forest Plot — {cat}\nOutcome: Epigenetic Age Acceleration (years)",
        filename=fname,
        show_summary=True
    )

    subgroup_results.append({
        'Exposure category':       cat,
        'k (studies)':             meta['k'],
        'Pooled β':                round(meta['theta'], 4),
        'SE':                      round(meta['se'], 4),
        '95% CI lower':            round(meta['ci_lo'], 4),
        '95% CI upper':            round(meta['ci_hi'], 4),
        'z':                       round(meta['z'], 3),
        'p (pooled β)':            round(meta['p_theta'], 4),
        'tau²':                    round(meta['tau2'], 4),
        'tau':                     round(meta['tau'], 4),
        'I² (%)':                  round(meta['I2'], 1),
        'Q':                       round(meta['Q'], 2),
        'df (Q)':                  meta['df_Q'],
        'p (Q heterogeneity)':     round(meta['p_Q'], 4),
        'Interpretation':          ('Significant' if meta['p_theta'] < 0.05 else 'Not significant'),
        'Heterogeneity level':     ('Low' if meta['I2'] < 25 else ('Moderate' if meta['I2'] < 75 else 'High')),
    })

# ── 8. Summary forest plot (one diamond per category) ───────────────────────
print("\nCreating summary forest plot...")

if len(subgroup_results) > 0:
    fig_h_summary = max(6, len(subgroup_results) * 0.75 + 3.5)
    fig, ax = plt.subplots(figsize=(15, fig_h_summary))
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')

    res_sorted = sorted(subgroup_results, key=lambda x: x['Pooled β'])
    k_s = len(res_sorted)

    xs_lo = [r['95% CI lower'] for r in res_sorted]
    xs_hi = [r['95% CI upper'] for r in res_sorted]
    xmin, xmax = robust_xlim(xs_lo, xs_hi)

    ax.set_xlim(xmin, xmax)

    row_step_s = 1.15
    y_vals = [row_step_s * (i + 1) for i in range(k_s)]
    y_h = row_step_s * (k_s + 1)
    ax.set_ylim(-0.5, y_h + 1.0)

    ax.axvline(0, color=COLORS['zeroline'], lw=1.2, ls='--', alpha=0.7)

    def ax_t(txt, ax_x, y_d, bold=False, fs=8.5, col=COLORS['text']):
        ax.text(
            ax_x, y_d,
            txt,
            transform=ax.get_yaxis_transform(),
            ha='left',
            va='center',
            fontsize=fs,
            fontweight='bold' if bold else 'normal',
            color=col
        )

    COL_CAT = -0.24
    COL_K   = -0.03
    COL_B   = 0.05
    COL_I2  = 0.31

    ax_t("Exposure category", COL_CAT, y_h, bold=True, fs=10, col=COLORS['header'])
    ax_t("k",                 COL_K,   y_h, bold=True, fs=10, col=COLORS['header'])
    ax_t("β (95% CI)",        COL_B,   y_h, bold=True, fs=10, col=COLORS['header'])
    ax_t("I²",                COL_I2,  y_h, bold=True, fs=10, col=COLORS['header'])
    ax.axhline(y_h - 0.35, color=COLORS['header'], lw=1.0, alpha=0.6)

    for i, r in enumerate(res_sorted):
        y = y_vals[i]
        lo = r['95% CI lower']
        hi = r['95% CI upper']
        th = r['Pooled β']

        ax.plot([lo, hi], [y, y], color=COLORS['ci_line'], lw=1.5, zorder=2)

        diam_h = 0.25
        col = COLORS['diamond'] if r['p (pooled β)'] < 0.05 else COLORS['nonsig']
        diamond = plt.Polygon(
            [[lo, y], [th, y + diam_h], [hi, y], [th, y - diam_h]],
            closed=True, facecolor=col, edgecolor='white', lw=0.6, zorder=4
        )
        ax.add_patch(diamond)

        i2c = i2_color(r['I² (%)'])
        ax_t(r['Exposure category'],               COL_CAT, y, fs=8.5)
        ax_t(str(r['k (studies)']),                COL_K,   y, fs=8.5)
        ax_t(f"{th:+.3f} ({lo:+.3f}, {hi:+.3f})",  COL_B,   y, fs=8.5)
        ax_t(f"{r['I² (%)']:.0f}%",                COL_I2,  y, fs=8.5, col=i2c)

    ax.set_xlabel("Pooled β (years of EAA per unit exposure)", fontsize=10,
                  color=COLORS['text'], labelpad=10)
    ax.set_yticks([])
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.spines['bottom'].set_color(COLORS['grid'])
    ax.tick_params(axis='x', labelsize=8.5)

    ax.set_title(
        "Summary forest plot — Subgroup analyses by exposure category\n"
        "Pool A: unstandardized β coefficients",
        fontsize=11, fontweight='bold', color=COLORS['header'], pad=10, loc='left'
    )

    leg = [
        mpatches.Patch(color=COLORS['diamond'], label='p < 0.05'),
        mpatches.Patch(color=COLORS['nonsig'], label='p ≥ 0.05'),
        Line2D([0], [0], color=COLORS['I2_low'],  lw=2, label='I² < 25% (low)'),
        Line2D([0], [0], color=COLORS['I2_mid'],  lw=2, label='25% ≤ I² < 75% (moderate)'),
        Line2D([0], [0], color=COLORS['I2_high'], lw=2, label='I² ≥ 75% (high)')
    ]
    ax.legend(handles=leg, fontsize=8, loc='lower right',
              frameon=True, framealpha=0.9, edgecolor=COLORS['grid'])

    plt.subplots_adjust(left=0.30, right=0.97, top=0.90, bottom=0.12)

    out_sum = os.path.join(OUTDIR, "forest_plot_summary_subgroups.pdf")
    plt.savefig(out_sum, dpi=220, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {out_sum}")

# ── 9. Save Excel results ────────────────────────────────────────────────────
print("\nSaving Excel results...")

overall_row = {
    'Exposure category':   f'OVERALL ({meta_overall["k"]} studies)',
    'k (studies)':         meta_overall['k'],
    'Pooled β':            round(meta_overall['theta'], 4),
    'SE':                  round(meta_overall['se'], 4),
    '95% CI lower':        round(meta_overall['ci_lo'], 4),
    '95% CI upper':        round(meta_overall['ci_hi'], 4),
    'z':                   round(meta_overall['z'], 3),
    'p (pooled β)':        round(meta_overall['p_theta'], 4),
    'tau²':                round(meta_overall['tau2'], 4),
    'tau':                 round(meta_overall['tau'], 4),
    'I² (%)':              round(meta_overall['I2'], 1),
    'Q':                   round(meta_overall['Q'], 2),
    'df (Q)':              meta_overall['df_Q'],
    'p (Q heterogeneity)': round(meta_overall['p_Q'], 4),
    'Interpretation':      ('Significant' if meta_overall['p_theta'] < 0.05 else 'Not significant'),
    'Heterogeneity level': ('Low' if meta_overall['I2'] < 25 else ('Moderate' if meta_overall['I2'] < 75 else 'High')),
}

results_df = pd.DataFrame([overall_row] + subgroup_results)

study_weights = df[['study_id', 'first_author', 'year', 'exposure_category', 'category', 'yi', 'sei']].copy()
study_weights['weight_re (%)'] = (meta_overall['w_re'] / meta_overall['w_re'].sum() * 100).round(2)
study_weights['vi'] = (df['sei'] ** 2).round(6)

out_xlsx = os.path.join(OUTDIR, "step7_results.xlsx")
with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
    results_df.to_excel(writer, sheet_name='meta_analysis_results', index=False)
    study_weights.to_excel(writer, sheet_name='study_weights', index=False)

    narrative_df = df[df['category'].isin(cats_narrative)][
        ['first_author', 'year', 'exposure_category', 'yi', 'sei']
    ].copy()
    narrative_df.to_excel(writer, sheet_name='narrative_only_lt3', index=False)

print(f"  Saved: {out_xlsx}")

# ── 10. Final report ─────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 7 COMPLETE — SUMMARY")
print("=" * 60)

print(f"\nOverall pooled β = {meta_overall['theta']:+.4f} "
      f"(95% CI: {meta_overall['ci_lo']:+.4f} to {meta_overall['ci_hi']:+.4f})")
print(f"I² = {meta_overall['I2']:.1f}%  |  tau² = {meta_overall['tau2']:.4f}")
print(f"Overall forest plot files: {overall_files}")

print(f"\nSubgroup results:")
for r in sorted(subgroup_results, key=lambda x: x['Pooled β'], reverse=True):
    sig = "**" if r['p (pooled β)'] < 0.05 else "  "
    print(f"  {sig} {r['Exposure category']:35s} β={r['Pooled β']:+.4f} "
          f"[{r['95% CI lower']:+.4f},{r['95% CI upper']:+.4f}] "
          f"I²={r['I² (%)']:.0f}% k={r['k (studies)']}")

print("\n** = statistically significant (p < 0.05)")
print(f"\nOutput files in: {OUTDIR}/")
