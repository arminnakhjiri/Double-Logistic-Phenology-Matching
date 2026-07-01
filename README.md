# Double-Logistic Phenology Matching

This repository provides a **Python workflow** for identifying raster clusters whose seasonal dynamics most closely resemble a reference phenological signature using **double-logistic curve fitting**.

Instead of directly comparing observed time series, the workflow fits a double-logistic model to both the reference and each cluster, allowing phenological similarity to be assessed using both fitted model parameters and reconstructed seasonal curves.

---

## 📌 Overview

- **Input 1:** Reference multi-band time-series GeoTIFF
- **Input 2:** Clustered multi-band time-series GeoTIFF
- **Input 3:** Cluster map
- **Method:** Double-logistic curve fitting
- **Output:** Ranked cluster similarity table and fitted phenology curves

---

## 🧰 Python Dependencies

Install the required packages:

```bash
pip install numpy scipy rasterio pandas matplotlib
```

---

## 🧪 Methodology

The workflow consists of the following steps:

1. Load the reference and clustered raster time series.
2. Extract acquisition dates from raster band descriptions.
3. Compute the mean temporal profile of the reference pixels.
4. Fit a double-logistic model to the reference curve.
5. Compute the mean temporal profile of every cluster.
6. Fit a double-logistic model to each cluster.
7. Compare fitted models using:
   - Pearson correlation
   - Root Mean Square Error (RMSE)
   - Euclidean distance between model parameters
8. Rank clusters according to their phenological similarity with the reference.
9. Export the ranked results and visualize the best-matching phenological curves.

---

## ⚙️ Adjustable Parameters

Several parameters can be modified at the beginning of the script:

```python
MIN_PIXELS_PER_CLUSTER = 50
MIN_DOY = 20
INITIAL_PARAMETERS = [0.5, 0.5, 120, 250, 0]
```

These allow the workflow to be adapted to different vegetation types, growing seasons, and temporal datasets.

---

## 📁 Outputs

The script produces:

- **double_logistic_cluster_ranking.csv** — Ranked similarity table for all fitted clusters
- **Phenology_Matching.png** — Comparison of the fitted reference curve and the highest-ranked clusters

---

## 🌍 Applications

Typical applications include:

- Crop identification
- Crop type mapping
- Phenological similarity analysis
- Vegetation monitoring
- Seasonal dynamics comparison
- Remote sensing time-series analysis
- Agricultural monitoring

---

## 📄 License

MIT License

---

## Author

**Armin Nakhjiri**

Remote Sensing Scientist

✉️ Nakhjiri.Armin@gmail.com

---

*Simple tools for efficient geospatial data processing.*
