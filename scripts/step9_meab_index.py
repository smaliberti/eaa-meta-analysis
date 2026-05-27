"""
step9_meab_index.py
====================
Meta-analysis: Modifiable determinants of epigenetic age acceleration
Step 9: MEAB-Index (Modifiable Epigenetic Aging Burden Index)

METHODOLOGICAL FRAMEWORK:
--------------------------
The MEAB-Index is a novel composite index that quantifies the total burden
of biological aging attributable to modifiable determinants, as measured
by epigenetic clocks (EAA). It is constructed from the precision-weighted
pooled beta estimates obtained in Step 7 (DerSimonian-Laird random-effects
model, Pool A, n=60 studies).

THREE COMPLEMENTARY METRICS:
1. MEAB-composite: precision-weighted average across all 7 interpretable
   exposure categories. Represents the average expected years of EAA per
   unit of modifiable exposure, across the full spectrum of determinants.

2. MEAB-positive: precision-weighted average restricted to statistically
   significant (p<0.05) categories with positive beta. Represents the
   conservative evidence-based estimate of accelerating exposures.

3. Cumulative Preventable Burden (CPB): sum of statistically significant
   positive pooled betas. Represents the upper-bound estimate of years of
   biological aging potentially preventable if all significant risk factors
   were simultaneously eliminated. Stated assumption: independent effects.

BOOTSTRAP CI: 10,000 iterations (non-parametric, percentile method)
              to account for uncertainty in the input pooled estimates.

This index is novel: no equivalent composite measure exists in the
epigenetic aging literature.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import os, warnings
from scipy import stats
warnings.filterwarnings('ignore')
np.random.seed(42)

OUTDIR = "step9_outputs"
os.makedirs(OUTDIR, exist_ok=True)

# ── palette ──────────────────────────────────────────────────────────────────
C = dict(
    blue   = '#185FA5', light_blue = '#B5D4F4', blue50 = '#E6F1FB',
    teal   = '#0F6E56', teal50     = '#E1F5EE',
    amber  = '#BA7517', amber50    = '#FAEEDA',
    red    = '#A32D2D', red50      = '#FCEBEB',
    green  = '#3B6D11', green50    = '#EAF3DE',
    gray   = '#888780', gray50     = '#F1EFE8',
    header = '#2E4C6D', txt        = '#2C2C2A',
    grid   = '#D3D1C7',
    sig    = '#185FA5',
    ns     = '#888780',
)

# ── Input data from Step 7 ───────────────────────────────────────────────────
# All 8 pooled estimates from Pool A subgroup analyses
# 'Other' excluded from MEAB (heterogeneous residual category)
CATEGORIES = [
    # name,                        k,   beta,    se,      p
    ('Environmental exposure',     9,  +0.4657, 0.1526, 0.0023),
    ('Metabolic / Inflammatory',   6,  +0.9130, 0.2271, 0.0001),
    ('Psychosocial stress',        8,  +0.1877, 0.0593, 0.0016),
    ('Socioeconomic / Social',     4,  +0.1742, 0.1176, 0.1386),
    ('Sleep',                      3,  +0.0063, 0.0784, 0.9360),
    ('Diet & Nutrition',          12,  -0.0626, 0.1326, 0.6371),
    ('Physical activity',          7,  -0.2783, 0.3501, 0.4266),
]

df = pd.DataFrame(CATEGORIES,
                  columns=['category','k','beta','se','p'])
df['vi']     = df['se']**2
df['wi']     = 1.0 / df['vi']
df['ci_lo']  = df['beta'] - 1.96*df['se']
df['ci_hi']  = df['beta'] + 1.96*df['se']
df['sig']    = df['p'] < 0.05
df['direction'] = df['beta'].apply(
    lambda x: 'Accelerates EAA' if x > 0 else 'Decelerates EAA')

# ── MEAB calculations ─────────────────────────────────────────────────────────
def meab_composite(data):
    """Precision-weighted composite across all categories."""
    w  = data['wi'].values
    b  = data['beta'].values
    W  = w.sum()
    mu = (w * b).sum() / W
    se = np.sqrt(1.0 / W)
    return mu, se, mu-1.96*se, mu+1.96*se

def meab_positive(data):
    """Precision-weighted composite — significant positive categories only."""
    sub = data[data['sig'] & (data['beta'] > 0)]
    if len(sub) == 0:
        return np.nan, np.nan, np.nan, np.nan
    w  = sub['wi'].values
    b  = sub['beta'].values
    W  = w.sum()
    mu = (w * b).sum() / W
    se = np.sqrt(1.0 / W)
    return mu, se, mu-1.96*se, mu+1.96*se

def cpb(data):
    """Cumulative Preventable Burden: sum of significant positive betas."""
    sub = data[data['sig'] & (data['beta'] > 0)]
    total = sub['beta'].sum()
    # SE via error propagation (assuming independence)
    se_total = np.sqrt((sub['se']**2).sum())
    return total, se_total, total-1.96*se_total, total+1.96*se_total

# Point estimates
m_comp, se_comp, lo_comp, hi_comp = meab_composite(df)
m_pos,  se_pos,  lo_pos,  hi_pos  = meab_positive(df)
m_cpb,  se_cpb,  lo_cpb,  hi_cpb  = cpb(df)

print("=" * 60)
print("MEAB-Index — Point Estimates")
print("=" * 60)
print(f"\nMEAB-composite : {m_comp:+.4f}  SE={se_comp:.4f}"
      f"  95% CI [{lo_comp:+.4f}, {hi_comp:+.4f}]")
print(f"MEAB-positive  : {m_pos:+.4f}  SE={se_pos:.4f}"
      f"  95% CI [{lo_pos:+.4f}, {hi_pos:+.4f}]")
print(f"CPB            : {m_cpb:+.4f}  SE={se_cpb:.4f}"
      f"  95% CI [{lo_cpb:+.4f}, {hi_cpb:+.4f}]")

# ── Bootstrap CI (10,000 iterations) ─────────────────────────────────────────
N_BOOT = 10_000
boot_comp = np.zeros(N_BOOT)
boot_pos  = np.zeros(N_BOOT)
boot_cpb  = np.zeros(N_BOOT)

for i in range(N_BOOT):
    # Simulate each category's estimate ~ N(beta_k, se_k²)
    boot_betas = np.array([
        np.random.normal(row['beta'], row['se'])
        for _, row in df.iterrows()
    ])
    df_b = df.copy()
    df_b['beta'] = boot_betas
    df_b['sig']  = df['sig']  # keep original significance flags

    # composite
    W = df_b['wi'].sum()
    boot_comp[i] = (df_b['wi'] * df_b['beta']).sum() / W

    # positive
    sub_p = df_b[df_b['sig'] & (df_b['beta'] > 0)]
    if len(sub_p) > 0:
        Wp = sub_p['wi'].sum()
        boot_pos[i] = (sub_p['wi'] * sub_p['beta']).sum() / Wp
    else:
        boot_pos[i] = np.nan

    # cpb
    sub_c = df_b[df_b['sig'] & (df_b['beta'] > 0)]
    boot_cpb[i] = sub_c['beta'].sum()

# Percentile CI
bci_comp = np.percentile(boot_comp, [2.5, 97.5])
bci_pos  = np.nanpercentile(boot_pos, [2.5, 97.5])
bci_cpb  = np.percentile(boot_cpb, [2.5, 97.5])

print(f"\nBootstrap 95% CI (10,000 iterations):")
print(f"  MEAB-composite : [{bci_comp[0]:+.4f}, {bci_comp[1]:+.4f}]")
print(f"  MEAB-positive  : [{bci_pos[0]:+.4f},  {bci_pos[1]:+.4f}]")
print(f"  CPB            : [{bci_cpb[0]:+.4f},  {bci_cpb[1]:+.4f}]")

# Contribution to MEAB-composite
df['contrib_pct'] = df['wi'] / df['wi'].sum() * 100
df['weighted_contribution'] = df['wi'] / df['wi'].sum() * df['beta']

# ── Figure 1: MEAB Component Chart ───────────────────────────────────────────
print("\nBuilding Figure 1: MEAB component chart...")

fig = plt.figure(figsize=(14, 9), facecolor='white')
gs  = gridspec.GridSpec(1, 2, figure=fig,
                        width_ratios=[0.58, 0.42],
                        left=0.04, right=0.97,
                        top=0.90, bottom=0.10,
                        wspace=0.06)

ax_l = fig.add_subplot(gs[0])  # left: dot plot per category
ax_r = fig.add_subplot(gs[1])  # right: MEAB summary bars

# ── LEFT: category dot plot ───────────────────────────────────────────────────
ax_l.set_facecolor('white')
for sp in ['top','right']:
    ax_l.spines[sp].set_visible(False)
ax_l.spines['bottom'].set_color(C['grid'])
ax_l.spines['left'].set_color(C['grid'])

df_plot = df.sort_values('beta').reset_index(drop=True)
y_pos   = np.arange(len(df_plot))

for i, row in df_plot.iterrows():
    col  = C['sig'] if row['sig'] else C['ns']
    fcol = C['blue50'] if (row['sig'] and row['beta']>0) else \
           C['teal50'] if (row['sig'] and row['beta']<=0) else \
           C['gray50']

    # CI bar
    ax_l.plot([row['ci_lo'], row['ci_hi']], [i, i],
              color=col, lw=1.8, alpha=0.7, zorder=2)

    # dot sized by k
    ms = 40 + row['k']*12
    ax_l.scatter(row['beta'], i, s=ms, color=col,
                 edgecolors='white', lw=0.8, zorder=4)

    # weight % text
    ax_l.text(max(row['ci_hi'], row['beta'])+0.04, i,
              f"  {row['contrib_pct']:.1f}%",
              va='center', fontsize=8, color=C['txt'])

    # significance marker
    if row['sig']:
        ax_l.text(row['beta'], i+0.32, '★',
                  ha='center', fontsize=9, color=col)

ax_l.axvline(0,  color=C['red'],    lw=0.9, ls=':', alpha=0.7)
ax_l.axvline(m_comp, color=C['blue'], lw=1.2, ls='--', alpha=0.8,
             label=f'MEAB-composite = {m_comp:+.3f}')
ax_l.axvspan(bci_comp[0], bci_comp[1],
             alpha=0.07, color=C['blue'])

ax_l.set_yticks(y_pos)
ax_l.set_yticklabels(df_plot['category'], fontsize=9)
ax_l.tick_params(axis='x', labelsize=9)
ax_l.set_xlabel('\u03b2  (years of EAA per unit exposure)',
                fontsize=10, color=C['txt'], labelpad=8)
ax_l.legend(fontsize=8.5, frameon=True, framealpha=0.9,
            edgecolor=C['grid'], loc='lower right')

# dot size legend
for k_val, lbl in [(3,'k = 3'), (8,'k = 8'), (12,'k = 12')]:
    ax_l.scatter([], [], s=40+k_val*12, color=C['gray'],
                 label=lbl, edgecolors='white', lw=0.6)
ax_l.legend(fontsize=8, frameon=True, framealpha=0.9,
            edgecolor=C['grid'], loc='lower right', ncol=2)

# ── RIGHT: MEAB summary ───────────────────────────────────────────────────────
ax_r.set_facecolor('white')
for sp in ['top','right','left']:
    ax_r.spines[sp].set_visible(False)
ax_r.spines['bottom'].set_color(C['grid'])
ax_r.yaxis.set_visible(False)

metrics = [
    ('MEAB-composite\n(all 7 categories)',
     m_comp, bci_comp[0], bci_comp[1], C['blue'],   C['blue50']),
    ('MEAB-positive\n(3 significant)',
     m_pos,  bci_pos[0],  bci_pos[1],  C['teal'],   C['teal50']),
    ('Cumulative\nPreventable Burden',
     m_cpb,  bci_cpb[0],  bci_cpb[1],  C['amber'],  C['amber50']),
]

y_m = [2.2, 1.1, 0.0]

for (label, mu, lo, hi, col, fcol), y in zip(metrics, y_m):
    # bar from 0 to mu
    ax_r.barh(y, mu, height=0.55, left=0,
              color=fcol, edgecolor=col, linewidth=1.0, zorder=2)
    # bootstrap CI
    ax_r.plot([lo, hi], [y, y], color=col, lw=2.0, zorder=3,
              solid_capstyle='round')
    ax_r.scatter(mu, y, s=80, color=col, zorder=4,
                 edgecolors='white', lw=0.8)

    # value label
    ax_r.text(hi + 0.02, y,
              f" {mu:+.3f}\n [{lo:+.3f}, {hi:+.3f}]",
              va='center', fontsize=8.5, color=col)

    # metric label
    ax_r.text(-0.02, y, label,
              ha='right', va='center', fontsize=8.5,
              color=col, fontweight='bold')

ax_r.axvline(0, color=C['grid'], lw=0.8)
ax_r.set_xlim(-0.25, 2.4)
ax_r.set_ylim(-0.5, 3.0)
ax_r.tick_params(axis='x', labelsize=9)
ax_r.set_xlabel('Years of EAA', fontsize=10, color=C['txt'], labelpad=8)

# title
fig.suptitle(
    "MEAB-Index — Modifiable Epigenetic Aging Burden Index\n"
    "Left: pooled \u03b2 per exposure category (dot size \u221d k).  "
    "\u2605 = statistically significant (p\u00a0<\u00a00.05).  "
    "Right: three MEAB metrics with bootstrap 95% CI.",
    fontsize=10.5, fontweight='bold', color=C['header'],
    x=0.01, ha='left', y=0.97, va='top')

out1 = os.path.join(OUTDIR, "meab_index_chart.pdf")
fig.savefig(out1, dpi=180, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f"  Saved: {out1}")

# ── Figure 2: Bootstrap distributions ────────────────────────────────────────
print("Building Figure 2: Bootstrap distributions...")

fig2, axes = plt.subplots(1, 3, figsize=(13, 4.5), facecolor='white')
titles = ['MEAB-composite', 'MEAB-positive', 'Cumulative Preventable Burden']
boot_arrays = [boot_comp, boot_pos[~np.isnan(boot_pos)], boot_cpb]
vals  = [m_comp,  m_pos,  m_cpb]
bcis  = [bci_comp, bci_pos, bci_cpb]
fcolors = [C['blue50'], C['teal50'], C['amber50']]
scolors = [C['blue'],   C['teal'],   C['amber']]

for ax, title, arr, val, bci, fc, sc in zip(
        axes, titles, boot_arrays, vals, bcis, fcolors, scolors):
    ax.set_facecolor('white')
    for sp in ['top','right']:
        ax.spines[sp].set_visible(False)
    ax.spines['left'].set_color(C['grid'])
    ax.spines['bottom'].set_color(C['grid'])

    ax.hist(arr, bins=80, color=fc, edgecolor=sc,
            linewidth=0.3, alpha=0.9, density=True)

    # CI shading
    mask = (arr >= bci[0]) & (arr <= bci[1])
    if mask.sum() > 0:
        ax.hist(arr[mask], bins=80, color=sc, alpha=0.25,
                density=True)

    ax.axvline(val,    color=sc,       lw=1.8, ls='-',  alpha=0.9,
               label=f'Estimate = {val:+.3f}')
    ax.axvline(bci[0], color=sc,       lw=1.2, ls='--', alpha=0.7,
               label=f'95% CI [{bci[0]:+.3f}, {bci[1]:+.3f}]')
    ax.axvline(bci[1], color=sc,       lw=1.2, ls='--', alpha=0.7)
    ax.axvline(0,      color=C['red'], lw=0.8, ls=':',  alpha=0.6)

    ax.set_title(title, fontsize=9.5, fontweight='bold', color=C['header'], pad=6)
    ax.set_xlabel('Years of EAA', fontsize=9, color=C['txt'])
    ax.set_ylabel('Density' if ax == axes[0] else '', fontsize=9, color=C['txt'])
    ax.legend(fontsize=7.5, frameon=True, framealpha=0.9, edgecolor=C['grid'])
    ax.tick_params(labelsize=8.5)

fig2.suptitle(
    "Bootstrap distributions of MEAB-Index metrics  (n\u00a0=\u00a010,000 iterations)",
    fontsize=10.5, fontweight='bold', color=C['header'],
    x=0.01, ha='left')
fig2.tight_layout(rect=[0, 0, 1, 0.93])

out2 = os.path.join(OUTDIR, "meab_bootstrap_distributions.pdf")
fig2.savefig(out2, dpi=180, bbox_inches='tight', facecolor='white')
plt.close(fig2)
print(f"  Saved: {out2}")

# ── Figure 3: Sensitivity analysis of MEAB ───────────────────────────────────
print("Building Figure 3: Sensitivity analysis...")

# Remove one category at a time and recompute all 3 MEAB metrics
sens_rows = []
for i, row in df.iterrows():
    df_s = df.drop(i)
    mc, _, lc, hc = meab_composite(df_s)
    mp, _, lp, hp = meab_positive(df_s)
    mb, _, lb, hb = cpb(df_s)
    sens_rows.append(dict(
        removed=row['category'],
        meab_comp=mc, meab_comp_lo=lc, meab_comp_hi=hc,
        meab_pos=mp,  meab_pos_lo=lp,  meab_pos_hi=hp,
        cpb=mb,       cpb_lo=lb,       cpb_hi=hb,
    ))
sens_df = pd.DataFrame(sens_rows)

fig3, ax = plt.subplots(figsize=(11, 5.5), facecolor='white')
ax.set_facecolor('white')
for sp in ['top','right']:
    ax.spines[sp].set_visible(False)
ax.spines['bottom'].set_color(C['grid'])
ax.spines['left'].set_color(C['grid'])

x_pos = np.arange(len(sens_df))
w3 = 0.25

for j, (col_val, col_lo, col_hi, color, label, offset) in enumerate([
    ('meab_comp', 'meab_comp_lo', 'meab_comp_hi', C['blue'],  'MEAB-composite', -w3),
    ('meab_pos',  'meab_pos_lo',  'meab_pos_hi',  C['teal'],  'MEAB-positive',   0),
    ('cpb',       'cpb_lo',       'cpb_hi',        C['amber'], 'CPB',            +w3),
]):
    xp = x_pos + offset
    ax.bar(xp, sens_df[col_val], width=w3*0.85,
           color=color, alpha=0.75, edgecolor='white', lw=0.5, label=label)
    ax.errorbar(xp, sens_df[col_val],
                yerr=[sens_df[col_val]-sens_df[col_lo],
                      sens_df[col_hi]-sens_df[col_val]],
                fmt='none', color=color, elinewidth=1.2, capsize=3)

# reference lines
ax.axhline(m_comp, color=C['blue'],  lw=1.0, ls='--', alpha=0.5)
ax.axhline(m_pos,  color=C['teal'],  lw=1.0, ls='--', alpha=0.5)
ax.axhline(m_cpb,  color=C['amber'], lw=1.0, ls='--', alpha=0.5)
ax.axhline(0,      color=C['grid'],  lw=0.6)

ax.set_xticks(x_pos)
ax.set_xticklabels(sens_df['removed'], rotation=25, ha='right', fontsize=8.5)
ax.tick_params(axis='y', labelsize=9)
ax.set_ylabel('MEAB value (years of EAA)', fontsize=9.5, color=C['txt'])
ax.legend(fontsize=9, frameon=True, framealpha=0.9, edgecolor=C['grid'])
ax.set_title(
    "Sensitivity analysis: MEAB-Index when one category is removed\n"
    "Dashed lines = full-model values",
    fontsize=10, fontweight='bold', color=C['header'], pad=8, loc='left')

fig3.tight_layout()
out3 = os.path.join(OUTDIR, "meab_sensitivity_categories.pdf")
fig3.savefig(out3, dpi=180, bbox_inches='tight', facecolor='white')
plt.close(fig3)
print(f"  Saved: {out3}")

# ── Save Excel ────────────────────────────────────────────────────────────────
print("\nSaving Excel...")

# Category-level detail
df_out = df[['category','k','beta','se','ci_lo','ci_hi','p','sig',
             'vi','wi','contrib_pct']].copy()
df_out.columns = ['Exposure category','k','Pooled β','SE',
                  '95% CI lower','95% CI upper','p-value','Significant',
                  'Variance (vi)','Weight (wi)','Contribution to MEAB (%)']
df_out = df_out.round(4)

# MEAB summary
summary = pd.DataFrame([
    {'Metric':'MEAB-composite',
     'Value':round(m_comp,4), 'SE':round(se_comp,4),
     '95% CI lower (analytical)':round(lo_comp,4),
     '95% CI upper (analytical)':round(hi_comp,4),
     'Bootstrap CI lower':round(bci_comp[0],4),
     'Bootstrap CI upper':round(bci_comp[1],4),
     'N categories':7,
     'Definition':'Precision-weighted average across all 7 interpretable categories'},
    {'Metric':'MEAB-positive',
     'Value':round(m_pos,4), 'SE':round(se_pos,4),
     '95% CI lower (analytical)':round(lo_pos,4),
     '95% CI upper (analytical)':round(hi_pos,4),
     'Bootstrap CI lower':round(bci_pos[0],4),
     'Bootstrap CI upper':round(bci_pos[1],4),
     'N categories':3,
     'Definition':'Precision-weighted average of significant positive categories'},
    {'Metric':'Cumulative Preventable Burden (CPB)',
     'Value':round(m_cpb,4), 'SE':round(se_cpb,4),
     '95% CI lower (analytical)':round(lo_cpb,4),
     '95% CI upper (analytical)':round(hi_cpb,4),
     'Bootstrap CI lower':round(bci_cpb[0],4),
     'Bootstrap CI upper':round(bci_cpb[1],4),
     'N categories':3,
     'Definition':'Sum of significant positive pooled betas (assumes independent effects)'},
])

sens_df_out = sens_df.round(4)

out_xlsx = os.path.join(OUTDIR, "step9_meab_results.xlsx")
with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
    summary.to_excel(writer,      sheet_name='MEAB_index_summary',     index=False)
    df_out.to_excel(writer,        sheet_name='category_contributions',  index=False)
    sens_df_out.to_excel(writer,   sheet_name='sensitivity_by_category', index=False)

print(f"  Saved: {out_xlsx}")

# ── Final report ──────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("STEP 9 — MEAB-INDEX COMPLETE")
print("="*60)
print(f"""
MEAB-composite = {m_comp:+.3f} years
  Bootstrap 95% CI: [{bci_comp[0]:+.3f}, {bci_comp[1]:+.3f}]
  Interpretation: across all modifiable determinants studied,
  the average association with EAA is +{m_comp:.2f} years per unit exposure.

MEAB-positive  = {m_pos:+.3f} years
  Bootstrap 95% CI: [{bci_pos[0]:+.3f}, {bci_pos[1]:+.3f}]
  Interpretation: among proven accelerating exposures (environmental,
  metabolic, psychosocial), the weighted average effect is +{m_pos:.2f} years.

Cumulative Preventable Burden = {m_cpb:+.3f} years
  Bootstrap 95% CI: [{bci_cpb[0]:+.3f}, {bci_cpb[1]:+.3f}]
  Interpretation: if environmental, metabolic, and psychosocial risk
  factors were simultaneously eliminated, up to {m_cpb:.2f} years of
  biological aging could potentially be prevented.

Dominant contributor to MEAB: Psychosocial stress (44.0% of weight)
  — due to its high precision (SE = 0.059), not its largest beta.
""")
