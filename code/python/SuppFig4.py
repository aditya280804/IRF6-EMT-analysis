import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------- Paths ----------------
base_dir = "./"
out_dir = os.path.join(base_dir, "figures/Suppfig")
os.makedirs(out_dir, exist_ok=True)

variance_file = os.path.join(base_dir, "data/pc1_variances.tsv")
loadings_file = os.path.join(base_dir, "data/pc1_loadings_all_networks.tsv")

# ---------------- Load data ----------------
variance_df = pd.read_csv(variance_file, sep='\t')
loadings_df = pd.read_csv(loadings_file, sep='\t', index_col=0)

teamscore_wt = 0.65

# ================== FIG S1A ==================
plt.figure(figsize=(5.3, 3))
sns.histplot(variance_df["PC1_Variance"], bins=20, element="step")

plt.axvline(x=teamscore_wt, color='red')
plt.text(
    teamscore_wt - 0.02,
    0.16 * plt.ylim()[1],
    "WT PC1 Variance (65%)",
    fontsize=9,
    color="red",
    rotation=90
)

plt.xlabel("PC1 Variance Explained", fontsize=10)
plt.ylabel("Frequency", fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "SuppFig_PC1_variance.png"), dpi=300)
plt.close()

# ================== FIG S1B ==================
mean_loadings = loadings_df.mean(axis=0)
std_loadings = loadings_df.std(axis=0)

plt.figure(figsize=(5.3, 3))
plt.bar(
    mean_loadings.index,
    mean_loadings.values,
    yerr=std_loadings.values,
    capsize=4,
    color="skyblue",
    edgecolor="black"
)

plt.xticks(fontsize=10, rotation=45)
plt.yticks(fontsize=10)
plt.ylabel("PC1 Loading Coefficient", fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "SuppFig_PC1_loadings.png"), dpi=300)
plt.close()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# -----------------------------
# Load your data
# -----------------------------
pc1_df = pd.read_csv(f"{base_dir}/data/pc1_variances.tsv", sep="\t")
team_df = pd.read_csv(f"{base_dir}/data/team_strength_random_networks.tsv", sep="\t")
df = pd.merge(pc1_df, team_df, on="network_id")

# -----------------------------
# Hill function with baseline
# -----------------------------
def hill_baseline(x, y0, Vmax, K, n):
    return y0 + Vmax * (x**n) / (K**n + x**n)

# Data
x = df["team_strength"].values
y = df["PC1_Variance"].values

# Initial guesses
p0 = [
    min(y),           # y0 (baseline)
    max(y)-min(y),    # Vmax
    np.median(x),     # K
    2.0               # n
]

# Optional bounds to prevent extreme values
bounds = (
    [0, 0, 0, 0.1],   # lower bounds: y0, Vmax, K, n
    [1.0, 1.0, 1.0, 5.0]  # upper bounds
)

params, _ = curve_fit(hill_baseline, x, y, p0=p0, bounds=bounds, maxfev=10000)
y0, Vmax, K, n = params

# Print fitted equation
print(f"Fitted Hill function with baseline:")
print(f"y = {y0:.2f} + {Vmax:.2f} * (x^{n:.2f} / (x^{n:.2f} + {K:.2f}^{n:.2f}))")

# -----------------------------
# Plot
# -----------------------------
x_fit = np.linspace(x.min(), x.max(), 300)
y_fit = hill_baseline(x_fit, y0, Vmax, K, n)
# Construct the equation string with actual values
eq_text = rf"$y = {y0:.2f} + {Vmax:.2f} \frac{{x^{{{n:.2f}}}}}{{x^{{{n:.2f}}} + {K:.2f}^{{{n:.2f}}}}}$"



plt.figure(figsize=(5.3, 3))
plt.scatter(x, y, s=15, alpha=0.5, edgecolors="none")
plt.plot(x_fit, y_fit, color="black", linewidth=1.5)

# WT point
plt.scatter(0.46, 0.65, color="red", s=50, label="WT network")
plt.annotate(
    "WT Network",
    xy=(0.46, 0.65),
    xytext=(0.95, 0.5),
    textcoords='axes fraction',
    arrowprops=dict(facecolor='red', arrowstyle='->', lw=1.5),
    fontsize=9, ha='right', va='bottom', color='red'
)

# Add the equation at bottom right
plt.text(
    0.95, 0.05,
    eq_text,
    transform=plt.gca().transAxes,
    ha='right',
    va='bottom',
    fontsize=9,
    color='black'
)

plt.xlabel("Team Strength", fontsize=10)
plt.ylabel("PC1 Variance Explained", fontsize=10)
plt.tight_layout()
plt.savefig(f"{base_dir}/figures/Suppfig/SuppFig_PC1_vs_TeamStrength_Hill_baseline.png", dpi=300)
plt.close()





