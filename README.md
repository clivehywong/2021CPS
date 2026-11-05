# Fronto-Cerebellar Connectivity Mediating Cognitive Processing Speed

[![Journal](https://img.shields.io/badge/Journal-NeuroImage%20(2021)-blue.svg)](https://doi.org/10.1016/j.neuroimage.2020.117556)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-green.svg)](https://www.r-project.org/)

Official data repository and analysis codebase for the empirical research paper:  
> **"Fronto-Cerebellar Connectivity Mediating Cognitive Processing Speed"**  
> Published in *NeuroImage*, Volume 226, 117556 (2021).  
> **Authors:** Clive H. Y. Wong, Jiao Liu, Tatia M. C. Lee, Jing Tao, Alex W. K. Wong, Bolton K. H. Chau, Lidian Chen, Chetwyn C. H. Chan  
> **DOI:** [10.1016/j.neuroimage.2020.117556](https://doi.org/10.1016/j.neuroimage.2020.117556)

---

## 📌 Executive Summary

This repository provides preprocessed functional MRI (fMRI) datasets, region of interest (ROI) masks, and full R analysis scripts to reproduce the findings of our study on individual differences in **Cognitive Processing Speed (CS)**. 

By employing a multi-modal (auditory & visual) cross-compatibility Arrow Task, first-order vector autoregressive (VAR(1)) modeling, and Mixed-Effect LASSO regression, this project explores the interregional effective connectivities underlying human processing efficiency.

### Key Scientific Findings:
- **Network Architecture:** Identified a 6-path effective connectivity network (3 long-range, 3 short-range) that reliably predicts individual cognitive processing speed.
- **Dual Mechanism:** Demonstrates that cognitive speed is driven by the dynamic interplay between **cerebellar-guided automaticity** (short-range) and **frontal top-down attentional control** (long-range).
- **Methodological Advance:** Multivariate vector autoregression (VAR(1)) significantly outperformed traditional bivariate functional connectivity measures (e.g., Pearson's correlation, gPPI) in predictive performance.

---

## 📁 Repository Layout

```text
2021CPS/
├── Connectivity/        # R scripts for calculating 6 connectivity metrics (gPPI, VAR(1), VAR(2), Pearson, etc.)
├── MixedEffectLasso/    # R scripts for feature selection & Mixed-Effect LASSO modeling
├── data/                # Preprocessed fMRI time-series, ANTs transformation matrices, and watershed ROI masks
├── LICENSE              # GNU General Public License v3.0
└── README.md            # Project documentation
