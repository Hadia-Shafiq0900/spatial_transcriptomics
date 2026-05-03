# Spatial Transcriptomics Analysis using Scanpy & Squidpy

## Overview

This project performs spatial transcriptomics analysis on a 10x Genomics Visium
mouse brain dataset using two powerful Python libraries — **Scanpy** and **Squidpy**.
Spatial transcriptomics allows us to measure gene expression while preserving the
spatial location of cells within a tissue, giving us biological context that
traditional single-cell RNA-seq cannot provide.

---

## Purpose

The goal of this assignment is to:
- Learn how to load and preprocess spatial transcriptomics data
- Perform quality control, normalization, and dimensionality reduction
- Cluster cells and identify spatial gene expression patterns
- Visualize results both in reduced dimensional space (UMAP) and on actual tissue images

---

## Dataset

- **Source:** 10x Genomics Visium — Mouse Brain Section
- **Loaded via:** `scanpy.datasets.visium_sge()`
- **Format:** AnnData object (`.h5ad`) containing:
  - Gene expression matrix
  - High-resolution tissue image
  - Spatial spot coordinates

---

## Libraries Used

| Library | Purpose |
|---|---|
| `scanpy` | Single-cell and spatial data analysis |
| `squidpy` | Spatial statistics and neighborhood analysis |
| `matplotlib` | Plotting and visualization |
| `seaborn` | Statistical visualizations |
| `igraph` + `leidenalg` | Graph-based clustering |

---

## Installation

```bash
pip install scanpy squidpy matplotlib seaborn igraph leidenalg
```

---

## Step-by-Step Workflow

### Step 1 — Load the Data
```python
import scanpy as sc
adata = sc.datasets.visium_sge(sample_id="V1_Mouse_Brain_Sagittal_Posterior")
adata.var_names_make_unique()
```
The dataset is loaded as an **AnnData** object. This is the standard data structure
used in single-cell and spatial analysis. It contains the gene expression matrix,
tissue image, and spot coordinates all in one object.

---

### Step 2 — Quality Control
```python
sc.pp.calculate_qc_metrics(adata, inplace=True)
```
We calculate quality control metrics like:
- **Total counts per spot** — spots with very low counts may be empty or damaged
- **Number of genes detected** — spots with very few genes are likely low quality

**Violin plots** are used to visualize the distribution of these metrics across all spots.

---

### Step 3 — Normalization & Log Transformation
```python
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
```
- **Normalization:** Scales each spot to have the same total count, removing
  sequencing depth bias
- **Log transformation:** Compresses the data range and makes gene expression
  distributions more normal, which improves downstream analysis

---

### Step 4 — Feature Selection
```python
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
```
We select the **top 2000 most variable genes** — these carry the most biological
signal and help us distinguish different cell types and regions.

---

### Step 5 — Dimensionality Reduction (PCA)
```python
sc.pp.pca(adata)
```
**Principal Component Analysis (PCA)** reduces thousands of gene dimensions down
to 50 principal components that capture the most variance in the data. This makes
computation faster and removes noise.
<img width="660" height="1079" alt="img 5" src="https://github.com/user-attachments/assets/82338325-1c8a-4e9d-a7bd-9870ab03416f" />

---

### Step 6 — Neighborhood Graph
```python
sc.pp.neighbors(adata)
```
A **k-nearest neighbor graph** is built in PCA space. Each spot is connected to
its most similar spots. This graph is the foundation for both clustering and UMAP.

---

### Step 7 — UMAP Visualization
```python
sc.tl.umap(adata)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts"])
```
**UMAP** projects the data into 2D for visualization. Spots that are
transcriptionally similar appear close together on the UMAP plot. This helps
us visually identify distinct cell populations.

<img width="661" height="921" alt="img6" src="https://github.com/user-attachments/assets/a5837d54-945e-4a1a-b0fc-f6c7fdfb12b2" />

---

### Step 8 — Leiden Clustering
```python
sc.tl.leiden(adata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)
```
**Leiden clustering** groups spots into clusters based on the neighborhood graph.
Each cluster ideally represents a distinct cell type or brain region. Results are
stored in `adata.obs["clusters"]`.

<img width="2294" height="1091" alt="img7" src="https://github.com/user-attachments/assets/37e9e41c-0391-4b41-81d3-fee4f51b9b56" />

---

### Step 9 — Spatial Visualization
```python
sc.pl.spatial(adata, img_key="hires", color=["clusters"])
```
Clusters are overlaid on the **actual tissue image**. This is the key advantage
of spatial transcriptomics — we can see exactly where each cell type is located
in the brain tissue.

<img width="2239" height="1091" alt="img8" src="https://github.com/user-attachments/assets/4ea12e7f-6575-4449-8766-dda28147c917" />

---

### Step 10 — Gene Expression in Space
```python
sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)
```
Individual gene expression is plotted on the tissue image:
- **COL1A2** — a marker for connective tissue / fibroblasts
- **SYPL1** — associated with synaptic vesicles in neurons

This shows us exactly which spatial regions of the brain express these genes.

<img width="1147" height="1079" alt="img9" src="https://github.com/user-attachments/assets/e7e20594-296c-4d05-bbde-6da2efb4cd74" />


---

## Results & Graphs Explained

| Plot | What it Shows |
|---|---|
| Violin plot (QC) | Distribution of total counts and genes per spot |
| Spatial QC plot | Which spots on the tissue have high/low counts |
| UMAP (QC metrics) | How quality metrics spread across transcriptomic space |
| UMAP (clusters) | Distinct transcriptional groups in 2D |
| Spatial cluster map | Cluster locations overlaid on brain tissue image |
| Spatial gene plots | COL1A2 and SYPL1 expression on tissue |

---

## Key Concepts

- **Spot:** A physical location on the Visium slide capturing ~10 cells
- **AnnData:** The data object storing expression matrix + metadata + images
- **UMAP:** Non-linear dimensionality reduction for 2D visualization
- **Leiden:** Graph-based community detection algorithm for clustering
- **Spatial transcriptomics:** Measuring gene expression with spatial coordinates preserved

---

## References

- [Scanpy Documentation](https://scanpy.readthedocs.io/)
- [Squidpy Documentation](https://squidpy.readthedocs.io/)
- [10x Genomics Visium](https://www.10xgenomics.com/spatial-transcriptomics)
- Tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html
