import re
import numpy as np
import rasterio
from scipy.optimize import curve_fit
from scipy.spatial.distance import euclidean
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# Input files

ref_path = r"path/to/reference_timeseries.tif"
cluster_stack_path = r"path/to/cluster_timeseries.tif"
cluster_map_path = r"path/to/cluster_map.tif"


# Parameters
MIN_PIXELS_PER_CLUSTER = 50
MIN_DOY = 20
INITIAL_PARAMETERS = [0.5, 0.5, 120, 250, 0]
OUTPUT_TABLE = "double_logistic_cluster_ranking.csv"
OUTPUT_FIGURE = "Phenology_Matching.png"

# Double-logistic model

def double_logistic(t, a, b, t1, t2, c):
    return (1 / (1 + np.exp(-a * (t - t1)))) - (1 / (1 + np.exp(-b * (t - t2)))) + c


# ============================================================
# HELPERS
# ============================================================

def get_doy(src):
    dates = [re.search(r"\d{8}", d).group() for d in src.descriptions]
    dates = [datetime.strptime(d, "%Y%m%d").timetuple().tm_yday for d in dates]
    return np.array(dates)


def mean_curve(stack):
    bands = stack.shape[0]
    X = stack.reshape(bands, -1).T
    X = X[np.all(np.isfinite(X), axis=1)]
    return np.nanmean(X, axis=0)


# ============================================================
# Load reference
# ============================================================

print("Loading reference...")

with rasterio.open(ref_path) as src:
    ref_stack = src.read()
    ref_doy = get_doy(src)

ref_curve = mean_curve(ref_stack)

# normalize
ref_curve = (ref_curve - np.mean(ref_curve)) / np.std(ref_curve)

# restrict to green-up period (important)
if MIN_DOY is not None:
    mask = ref_doy >= MIN_DOY
    ref_doy = ref_doy[mask]
    ref_curve = ref_curve[mask]

# initial guess
p0 = INITIAL_PARAMETERS


# fit reference model
ref_params, _ = curve_fit(double_logistic, ref_doy, ref_curve, p0=p0, maxfev=10000)

print("Reference parameters:", ref_params)


# ============================================================
# READ CLUSTER STACK
# ============================================================

print("Loading cluster stack...")

with rasterio.open(cluster_stack_path) as src:
    stack = src.read()
    dates = get_doy(src)

bands, rows, cols = stack.shape

X = stack.reshape(bands, -1).T
valid = np.all(np.isfinite(X), axis=1)
X = X[valid]

# normalize pixel-wise
X = (X - X.mean(axis=1, keepdims=True)) / (X.std(axis=1, keepdims=True) + 1e-6)

# load cluster map
with rasterio.open(cluster_map_path) as src:
    cluster_map = src.read(1).flatten()[valid]

clusters = np.unique(cluster_map)
clusters = clusters[clusters > 0]

print(f"Clusters found: {len(clusters)}")

# ============================================================
# FIT EACH CLUSTER
# ============================================================

results = []

for cid in clusters:

    pixels = X[cluster_map == cid]

    if len(pixels) < MIN_PIXELS_PER_CLUSTER:
        continue

    curve = np.mean(pixels, axis=0)

    try:
        params, _ = curve_fit(
            double_logistic,
            dates,
            curve,
            p0=p0,
            maxfev=10000
        )

        # distance in parameter space
        param_dist = euclidean(params, ref_params)

        # curve similarity
        fitted_curve = double_logistic(dates, *params)
        ref_fit = double_logistic(dates, *ref_params)

        corr = np.corrcoef(fitted_curve, ref_fit)[0, 1]

        rmse = np.sqrt(np.mean((fitted_curve - ref_fit) ** 2))

        results.append([
            cid,
            len(pixels),
            corr,
            rmse,
            param_dist,
            *params
        ])

    except (RuntimeError, ValueError):
        continue


# ============================================================
# RESULTS TABLE
# ============================================================

df = pd.DataFrame(results, columns=[
    "Cluster", "Pixels", "Corr", "RMSE", "ParamDist",
    "a", "b", "t1", "t2", "c"
])

df = df.sort_values(["Corr", "ParamDist"], ascending=[False, True])

print(df)

df.to_csv(OUTPUT_TABLE, index=False)

# ============================================================
# PLOT TOP 3 CLUSTERS
# ============================================================

top = df.head(3)

plt.figure(figsize=(10,6))

t = np.linspace(min(dates), max(dates), 200)

# reference
plt.plot(t, double_logistic(t, *ref_params), 'k', linewidth=3, label="Reference")

for _, row in top.iterrows():
    plt.plot(
        t,
        double_logistic(t, row.a, row.b, row.t1, row.t2, row.c),
        label=f"Cluster {int(row.Cluster)}"
    )

plt.legend()
plt.title("Double Logistic Phenology Matching")
plt.xlabel("DOY")
plt.ylabel("Normalized NDVI")
plt.tight_layout()
plt.savefig(OUTPUT_FIGURE, dpi=300)
plt.show()



print("Done.")
