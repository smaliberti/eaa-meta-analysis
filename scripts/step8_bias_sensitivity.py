"""
step8_bias_sensitivity.py
=========================
Meta-analysis: Modifiable determinants of epigenetic age acceleration
Step 8: Publication bias assessment + Sensitivity analysis

ANALYSES:
1. Funnel plot (overall + by exposure category)
2. Egger's test (linear regression of effect on precision)
3. Begg's test (rank correlation)
4. Leave-one-out sensitivity analysis (overall + by category)
5. Trim-and-fill correction (Duval & Tweedie)

All analyses on Pool A (n = 60, unstandardized beta coefficients)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os, warnings
from scipy import stats
warnings.filterwarnings('ignore')

OUTDIR = "step8_outputs"
os.makedirs(OUTDIR, exist_ok=True)

# ── palette ──────────────────────────────────────────────────────────────────
C = dict(
    blue   = '#185FA5', light_blue = '#B5D4F4',
    red    = '#E24B4A', orange     = '#BA7517',
    green  = '#3B6D11', gray       = '#888780',
    header = '#2E4C6D', txt        = '#2C2C2A',
    grid   = '#D3D1C7', row_alt    = '#F5F4F0',
    fill   = '#E24B4A',  # trim-and-fill imputed studies
)

# ── DL model ─────────────────────────────────────────────────────────────────
def dl_meta(yi, sei):
    from scipy.stats import norm, chi2 as chi2d
    yi, sei = np.array(yi, float), np.array(sei, float)
    vi = sei**2; k = len(yi)
    w_fe = 1/vi
    th_fe = np.sum(w_fe*yi)/np.sum(w_fe)
    Q = np.sum(w_fe*(yi-th_fe)**2)
    Cv = np.sum(w_fe) - np.sum(w_fe**2)/np.sum(w_fe)
    tau2 = max(0., (Q-(k-1))/Cv)
    w_re = 1/(vi+tau2)
    th = np.sum(w_re*yi)/np.sum(w_re)
    se = np.sqrt(1/np.sum(w_re))
    I2 = max(0., (Q-(k-1))/Q*100) if Q > 0 else 0.
    z = th/se
    return dict(k=k, theta=th, se=se,
                ci_lo=th-1.96*se, ci_hi=th+1.96*se,
                tau2=tau2, I2=I2, Q=Q, dfQ=k-1,
                pQ=1-chi2d.cdf(Q,k-1),
                p=2*(1-norm.cdf(abs(z))), w_re=w_re)

# ── Egger's test ──────────────────────────────────────────────────────────────
def egger_test(yi, sei):
    """
    Egger's regression test for funnel plot asymmetry.
    Regresses standardized effect (yi/sei) on precision (1/sei).
    Intercept ≠ 0 indicates asymmetry (potential publication bias).
    H0: intercept = 0 (no asymmetry)
    """
    yi, sei = np.array(yi, float), np.array(sei, float)
    precision = 1.0 / sei
    std_effect = yi / sei
    slope, intercept, r, p, se_intercept = stats.linregress(precision, std_effect)
    # t-test on intercept
    n = len(yi)
    t_stat = intercept / se_intercept
    p_intercept = 2 * stats.t.sf(abs(t_stat), df=n-2)
    return dict(intercept=intercept, slope=slope,
                se_intercept=se_intercept, t=t_stat,
                p=p_intercept, n=n,
                interpretation='Asymmetry detected (p < 0.10)' if p_intercept < 0.10
                               else 'No significant asymmetry (p ≥ 0.10)')

# ── Begg's test ───────────────────────────────────────────────────────────────
def begg_test(yi, sei):
    """
    Begg & Mazumdar rank correlation test.
    Tests whether standardized effect sizes are correlated with their variances.
    Uses Kendall's tau.
    """
    yi, sei = np.array(yi, float), np.array(sei, float)
    vi = sei**2
    # Adjusted effect sizes (subtract pooled estimate)
    m = dl_meta(yi, sei)
    adj_yi = yi - m['theta']
    tau, p = stats.kendalltau(vi, adj_yi)
    return dict(tau=tau, p=p, n=len(yi),
                interpretation='Asymmetry detected (p < 0.10)' if p < 0.10
                               else 'No significant asymmetry (p ≥ 0.10)')

# ── Trim and fill ─────────────────────────────────────────────────────────────
def trim_and_fill(yi, sei, side='left', max_iter=100):
    """
    Duval & Tweedie trim-and-fill method.
    Estimates number of missing studies and adjusted pooled estimate.
    side: 'left' = missing studies on left (negative) side
          'right' = missing studies on right side
    """
    yi, sei = np.array(yi, float), np.array(sei, float)
    n = len(yi)

    def pool_estimate(y, s):
        m = dl_meta(y, s)
        return m['theta']

    theta0 = pool_estimate(yi, sei)
    k0 = 0

    for _ in range(max_iter):
        # Center effects
        centered = yi - theta0
        if side == 'left':
            ranks = stats.rankdata(np.abs(centered))
            signs = np.sign(centered)
            # Count studies on the smaller side
            n_pos = np.sum(centered > 0)
            n_neg = np.sum(centered < 0)
            k_est = max(0, n_pos - n_neg)
        else:
            n_pos = np.sum(centered > 0)
            n_neg = np.sum(centered < 0)
            k_est = max(0, n_neg - n_pos)

        if k_est == k0:
            break
        k0 = k_est

        # Add imputed studies (mirror of the most extreme on dominant side)
        if k_est > 0:
            if side == 'left':
                extreme = np.sort(centered)[-k_est:]
                imputed_yi  = theta0 - extreme
            else:
                extreme = np.sort(centered)[:k_est]
                imputed_yi  = theta0 - extreme
            imputed_sei = sei[np.argsort(np.abs(centered))[-k_est:]]
            aug_yi  = np.concatenate([yi,  imputed_yi])
            aug_sei = np.concatenate([sei, imputed_sei])
            theta0  = pool_estimate(aug_yi, aug_sei)
        else:
            break

    # Final adjusted model
    if k0 > 0:
        if side == 'left':
            extreme = np.sort(yi - theta0)[-k0:]
            imp_yi  = theta0 - extreme
        else:
            extreme = np.sort(yi - theta0)[:k0]
            imp_yi  = theta0 - extreme
        imp_sei = sei[np.argsort(np.abs(yi - theta0))[-k0:]]
        aug_yi  = np.concatenate([yi,  imp_yi])
        aug_sei = np.concatenate([sei, imp_sei])
        adj_meta = dl_meta(aug_yi, aug_sei)
    else:
        adj_meta = dl_meta(yi, sei)
        imp_yi  = np.array([])
        imp_sei = np.array([])

    return dict(k_imputed=k0, imputed_yi=imp_yi, imputed_sei=imp_sei,
                adj_theta=adj_meta['theta'],
                adj_ci_lo=adj_meta['ci_lo'],
                adj_ci_hi=adj_meta['ci_hi'],
                adj_se=adj_meta['se'])

# ── Leave-one-out ─────────────────────────────────────────────────────────────
def leave_one_out(yi, sei, labels):
    yi, sei = np.array(yi, float), np.array(sei, float)
    results = []
    for i in range(len(yi)):
        mask = np.ones(len(yi), dtype=bool); mask[i] = False
        m = dl_meta(yi[mask], sei[mask])
        results.append(dict(
            removed=labels[i],
            theta=round(m['theta'], 4),
            ci_lo=round(m['ci_lo'], 4),
            ci_hi=round(m['ci_hi'], 4),
            I2=round(m['I2'], 1),
            p=round(m['p'], 4),
        ))
    return pd.DataFrame(results)

# ════════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ════════════════════════════════════════════════════════════════════════════
df = pd.read_excel('STEP6_.xlsx', sheet_name='pool_A_primary_beta')
df['yi']  = pd.to_numeric(df['yi'],  errors='coerce')
df['sei'] = pd.to_numeric(df['sei'], errors='coerce')
df = df[df['sei'] > 0].dropna(subset=['yi','sei']).reset_index(drop=True)

def harmonize(e):
    if pd.isna(e): return 'Other'
    e = str(e).lower()
    if 'smok' in e: return 'Smoking'
    if any(x in e for x in ['physical activity','exercise','sedent','walking','fitness','gardening']): return 'Physical activity'
    if any(x in e for x in ['diet','nutriti','food','vitamin','fiber','fatty','flavon','iron','sugar','alcohol','drink','cannabis','marijuana']): return 'Diet & Nutrition'
    if 'sleep' in e: return 'Sleep'
    if any(x in e for x in ['stress','psycho','anxiety','depress','mental','adversit','trauma','childhood','personality']): return 'Psychosocial stress'
    if any(x in e for x in ['pollut','air','environmental','toxin','chemical','pesticide','toluene','greenness','temperature','climate','urban']): return 'Environmental exposure'
    if any(x in e for x in ['socioeconom','ses','income','education','occupat','neighborhood','social','deprivat','sdh','mobility']): return 'Socioeconomic/Social'
    if any(x in e for x in ['bmi','obes','metabol','inflam','immune','biomark','insulin','lipid','body','syndrome','sii']): return 'Metabolic/Inflammatory'
    return 'Other'

df['category'] = df['exposure_category'].apply(harmonize)
yi_all  = df['yi'].values
sei_all = df['sei'].values
labels  = df['first_author'].values

meta_all = dl_meta(yi_all, sei_all)

print(f"Pool A loaded: {len(df)} studies")
print(f"Overall pooled β = {meta_all['theta']:+.4f} "
      f"(95% CI: {meta_all['ci_lo']:+.4f} to {meta_all['ci_hi']:+.4f})")

# ════════════════════════════════════════════════════════════════════════════
# 1. FUNNEL PLOT — OVERALL
# ════════════════════════════════════════════════════════════════════════════
print("\n--- Funnel plot (overall) ---")

def funnel_plot(yi, sei, meta, title, filename, imputed_yi=None, imputed_sei=None):
    fig, ax = plt.subplots(figsize=(9, 7), facecolor='white')
    ax.set_facecolor('white')
    for sp in ['top','right']:
        ax.spines[sp].set_visible(False)
    ax.spines['left'].set_color(C['grid'])
    ax.spines['bottom'].set_color(C['grid'])

    th   = meta['theta']
    tau2 = meta['tau2']
    se_max = sei.max() * 1.15

    # pseudo-CI boundaries (based on pooled + tau)
    se_range = np.linspace(0, se_max, 200)
    ci95_lo = th - 1.96 * np.sqrt(se_range**2 + tau2)
    ci95_hi = th + 1.96 * np.sqrt(se_range**2 + tau2)
    ci99_lo = th - 3.29 * np.sqrt(se_range**2 + tau2)
    ci99_hi = th + 3.29 * np.sqrt(se_range**2 + tau2)

    ax.fill_betweenx(se_range, ci99_lo, ci99_hi,
                     alpha=0.08, color=C['blue'], label='99% pseudo-CI')
    ax.fill_betweenx(se_range, ci95_lo, ci95_hi,
                     alpha=0.14, color=C['blue'], label='95% pseudo-CI')
    ax.axvline(th, color=C['blue'], lw=1.2, ls='--', alpha=0.8)
    ax.axvline(0,  color=C['red'],  lw=0.8, ls=':', alpha=0.6)

    # observed studies
    ax.scatter(yi, sei, s=55, color=C['blue'], edgecolors='white',
               linewidths=0.5, alpha=0.85, zorder=4, label='Observed studies')

    # imputed studies (trim-and-fill)
    if imputed_yi is not None and len(imputed_yi) > 0:
        ax.scatter(imputed_yi, imputed_sei, s=55, color=C['fill'],
                   edgecolors='white', linewidths=0.5, alpha=0.85,
                   marker='D', zorder=5, label=f'Imputed studies (k={len(imputed_yi)})')

    ax.set_xlabel('\u03b2  (years of EAA per unit exposure)',
                  fontsize=10, color=C['txt'], labelpad=8)
    ax.set_ylabel('Standard Error (SE)',
                  fontsize=10, color=C['txt'], labelpad=8)
    ax.invert_yaxis()   # SE=0 at top (most precise) — standard convention
    ax.set_ylim(se_max, -se_max*0.05)
    ax.tick_params(labelsize=9, colors=C['txt'])
    ax.legend(fontsize=8.5, frameon=True, framealpha=0.9,
              edgecolor=C['grid'], loc='upper right')

    ax.set_title(title, fontsize=10, fontweight='bold',
                 color=C['header'], pad=10, loc='left')

    het_txt = (f"Pooled \u03b2 = {th:+.3f}  \u00b7  "
               f"I\u00b2 = {meta['I2']:.1f}%  \u00b7  "
               f"k = {meta['k']} studies")
    ax.text(0.02, 0.02, het_txt, transform=ax.transAxes,
            fontsize=8, color=C['header'],
            bbox=dict(boxstyle='round,pad=0.25', facecolor='#EBF0F7',
                      edgecolor=C['blue'], lw=0.7))

    plt.tight_layout()
    out = os.path.join(OUTDIR, filename)
    fig.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {out}")

funnel_plot(yi_all, sei_all, meta_all,
            title=("Funnel Plot \u2014 Overall  (Pool A, n\u00a0=\u00a060)\n"
                   "Outcome: Epigenetic Age Acceleration (years)"),
            filename="funnel_plot_overall.pdf")

# ════════════════════════════════════════════════════════════════════════════
# 2. EGGER'S TEST
# ════════════════════════════════════════════════════════════════════════════
print("\n--- Egger's test ---")
egger = egger_test(yi_all, sei_all)
print(f"  Intercept = {egger['intercept']:+.4f} (SE = {egger['se_intercept']:.4f})")
print(f"  t = {egger['t']:.3f},  p = {egger['p']:.4f}")
print(f"  → {egger['interpretation']}")

# ════════════════════════════════════════════════════════════════════════════
# 3. BEGG'S TEST
# ════════════════════════════════════════════════════════════════════════════
print("\n--- Begg's test ---")
begg = begg_test(yi_all, sei_all)
print(f"  Kendall's tau = {begg['tau']:+.4f},  p = {begg['p']:.4f}")
print(f"  → {begg['interpretation']}")

# ════════════════════════════════════════════════════════════════════════════
# 4. TRIM AND FILL
# ════════════════════════════════════════════════════════════════════════════
print("\n--- Trim-and-fill (Duval & Tweedie) ---")
tf = trim_and_fill(yi_all, sei_all, side='left')
print(f"  Imputed studies (k) = {tf['k_imputed']}")
print(f"  Adjusted pooled β   = {tf['adj_theta']:+.4f} "
      f"(95% CI: {tf['adj_ci_lo']:+.4f} to {tf['adj_ci_hi']:+.4f})")
print(f"  Original pooled β   = {meta_all['theta']:+.4f}")

# Funnel plot with trim-and-fill
funnel_plot(yi_all, sei_all, meta_all,
            title=("Funnel Plot with Trim-and-Fill \u2014 Pool A (n\u00a0=\u00a060)\n"
                   "Red diamonds = imputed studies (Duval & Tweedie)"),
            filename="funnel_plot_trim_fill.pdf",
            imputed_yi  = tf['imputed_yi'],
            imputed_sei = tf['imputed_sei'])

# ════════════════════════════════════════════════════════════════════════════
# 5. LEAVE-ONE-OUT SENSITIVITY ANALYSIS
# ════════════════════════════════════════════════════════════════════════════
print("\n--- Leave-one-out sensitivity analysis ---")
loo_df = leave_one_out(yi_all, sei_all, labels)

# range of pooled estimates when removing each study
theta_range = loo_df['theta'].max() - loo_df['theta'].min()
most_influential = loo_df.loc[(loo_df['theta'] - meta_all['theta']).abs().idxmax(), 'removed']
print(f"  Range of pooled β when removing one study: "
      f"[{loo_df['theta'].min():+.4f}, {loo_df['theta'].max():+.4f}]")
print(f"  Most influential study: {most_influential}")
print(f"  Original pooled β: {meta_all['theta']:+.4f}")
print(f"  All LOO estimates remain {'positive' if (loo_df['theta'] > 0).all() else 'mixed sign'}")

# Leave-one-out plot
def loo_plot(loo_df, meta, filename):
    k = len(loo_df)
    fig_h = max(k * 0.32 + 2.5, 8)
    fig, ax = plt.subplots(figsize=(12, fig_h), facecolor='white')
    ax.set_facecolor('white')
    for sp in ['top','right','left']:
        ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_color(C['grid'])

    # sort by theta for readability
    loo_s = loo_df.sort_values('theta').reset_index(drop=True)
    y_pos = np.arange(len(loo_s))

    # CI bars
    for i, row in loo_s.iterrows():
        col = C['blue'] if row['p'] < 0.05 else C['gray']
        ax.plot([row['ci_lo'], row['ci_hi']], [i, i],
                color=col, lw=1.2, alpha=0.6, zorder=2)
        ax.scatter(row['theta'], i, s=30, color=col,
                   edgecolors='white', lw=0.4, zorder=3)

    # original pooled estimate
    ax.axvline(meta['theta'], color=C['blue'], lw=1.5, ls='--',
               alpha=0.9, label=f"Original pooled \u03b2 = {meta['theta']:+.3f}")
    ax.axvline(0, color=C['red'], lw=0.8, ls=':', alpha=0.5)

    # shaded CI band of original estimate
    ax.axvspan(meta['ci_lo'], meta['ci_hi'],
               alpha=0.08, color=C['blue'], label='Original 95% CI')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(loo_s['removed'], fontsize=7)
    ax.tick_params(axis='x', labelsize=8.5)
    ax.set_xlabel('Pooled \u03b2 when study removed  (years of EAA)',
                  fontsize=9.5, color=C['txt'], labelpad=8)
    ax.legend(fontsize=8.5, frameon=True, framealpha=0.9,
              edgecolor=C['grid'], loc='lower right')
    ax.set_title("Leave-One-Out Sensitivity Analysis \u2014 Pool A (n\u00a0=\u00a060)\n"
                 "Each point shows the pooled \u03b2 when that study is excluded",
                 fontsize=10, fontweight='bold', color=C['header'],
                 pad=10, loc='left')

    plt.tight_layout()
    out = os.path.join(OUTDIR, filename)
    fig.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {out}")

loo_plot(loo_df, meta_all, "leave_one_out_overall.pdf")

# ════════════════════════════════════════════════════════════════════════════
# 6. FUNNEL PLOTS BY CATEGORY (categories with k >= 5)
# ════════════════════════════════════════════════════════════════════════════
print("\n--- Funnel plots by category ---")
cat_counts = df['category'].value_counts()
cats_funnel = cat_counts[cat_counts >= 5].index.tolist()

cat_egger = []
for c in sorted(cats_funnel):
    sub  = df[df['category']==c].reset_index(drop=True)
    meta = dl_meta(sub['yi'], sub['sei'])
    eg   = egger_test(sub['yi'], sub['sei'])
    cat_egger.append(dict(
        category=c, k=meta['k'],
        beta=round(meta['theta'],3),
        egger_intercept=round(eg['intercept'],4),
        egger_p=round(eg['p'],4),
        interpretation=eg['interpretation']
    ))
    fname = f"funnel_plot_{c.lower().replace(' ','_').replace('/','_').replace('&','and')}.pdf"
    funnel_plot(sub['yi'].values, sub['sei'].values, meta,
                title=f"Funnel Plot \u2014 {c}  (k\u00a0=\u00a0{meta['k']})",
                filename=fname)
    print(f"  {c}: Egger p = {eg['p']:.4f} — {eg['interpretation']}")

# ════════════════════════════════════════════════════════════════════════════
# 7. SAVE EXCEL RESULTS
# ════════════════════════════════════════════════════════════════════════════
print("\nSaving Excel results...")

bias_summary = pd.DataFrame([{
    'Analysis':          'Egger\'s test (overall)',
    'Statistic':         f"t = {egger['t']:.3f}",
    'Intercept/Tau':     round(egger['intercept'], 4),
    'SE':                round(egger['se_intercept'], 4),
    'p-value':           round(egger['p'], 4),
    'Interpretation':    egger['interpretation'],
    'Studies (k)':       egger['n'],
}, {
    'Analysis':          "Begg's test (overall)",
    'Statistic':         f"Kendall tau = {begg['tau']:.4f}",
    'Intercept/Tau':     round(begg['tau'], 4),
    'SE':                '',
    'p-value':           round(begg['p'], 4),
    'Interpretation':    begg['interpretation'],
    'Studies (k)':       begg['n'],
}, {
    'Analysis':          'Trim-and-fill (Duval & Tweedie)',
    'Statistic':         f"k imputed = {tf['k_imputed']}",
    'Intercept/Tau':     '',
    'SE':                '',
    'p-value':           '',
    'Interpretation':    (f"Adjusted β = {tf['adj_theta']:+.4f} "
                          f"[{tf['adj_ci_lo']:+.4f}, {tf['adj_ci_hi']:+.4f}]"),
    'Studies (k)':       60 + tf['k_imputed'],
}])

tf_comparison = pd.DataFrame([{
    'Model':     'Original (k = 60)',
    'Pooled β':  round(meta_all['theta'], 4),
    '95% CI lower': round(meta_all['ci_lo'], 4),
    '95% CI upper': round(meta_all['ci_hi'], 4),
    'I² (%)':    round(meta_all['I2'], 1),
}, {
    'Model':     f"Trim-and-fill adjusted (k = {60+tf['k_imputed']})",
    'Pooled β':  round(tf['adj_theta'], 4),
    '95% CI lower': round(tf['adj_ci_lo'], 4),
    '95% CI upper': round(tf['adj_ci_hi'], 4),
    'I² (%)':    '',
}])

cat_egger_df = pd.DataFrame(cat_egger)

out_xlsx = os.path.join(OUTDIR, "step8_results.xlsx")
with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
    bias_summary.to_excel(writer,    sheet_name='bias_tests_summary',    index=False)
    tf_comparison.to_excel(writer,   sheet_name='trim_fill_comparison',  index=False)
    loo_df.to_excel(writer,          sheet_name='leave_one_out_overall', index=False)
    cat_egger_df.to_excel(writer,    sheet_name='egger_by_category',     index=False)

print(f"  Saved: {out_xlsx}")

# ════════════════════════════════════════════════════════════════════════════
# 8. FINAL REPORT
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("STEP 8 COMPLETE — SUMMARY")
print("="*60)
print(f"\nEgger's test:    intercept = {egger['intercept']:+.4f}, "
      f"p = {egger['p']:.4f} → {egger['interpretation']}")
print(f"Begg's test:     tau = {begg['tau']:+.4f}, "
      f"p = {begg['p']:.4f} → {begg['interpretation']}")
print(f"Trim-and-fill:   k imputed = {tf['k_imputed']}, "
      f"adjusted β = {tf['adj_theta']:+.4f} "
      f"[{tf['adj_ci_lo']:+.4f}, {tf['adj_ci_hi']:+.4f}]")
print(f"\nLeave-one-out:   β range = "
      f"[{loo_df['theta'].min():+.4f}, {loo_df['theta'].max():+.4f}]")
print(f"  All LOO estimates positive: {(loo_df['theta'] > 0).all()}")
print(f"  All LOO estimates significant (p<0.05): {(loo_df['p'] < 0.05).all()}")
print(f"\nOutput files in: {OUTDIR}/")
