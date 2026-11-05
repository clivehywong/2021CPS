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

## 📁 Repository Structure & Key Scripts

```text
2021CPS/
├── Connectivity/
│   ├── funcConn.R        # R script for computing multivariate VAR(1), VAR(2), Pearson, and partial correlations
│   └── PPI.R             # R script for generalized psychophysiological interaction (gPPI) estimation
├── MixedEffectLasso/
│   └── MMLasso.R         # R script for Mixed-Effect LASSO variable selection and model fitting
├── data/                 # ICA-cleaned fMRI time-series, ANTs matrices, and ROI masks
├── .gitignore
├── LICENSE               # GNU General Public License v3.0
└── README.md             # Main project documentation
```

---

## 🌐 Core Findings: Optimal VAR(1) Model

Among the 6 connectivity measures, the **first-order vector autoregressive model (VAR(1))** combined with `MMLasso.R` yielded the best predictive fit for Cognitive Speed (CS). The optimal model identified **6 significant effective connectivity paths**:

### Long-range Paths (Cognitive Control):

* **LMFC → MCV6** (Left Medial Frontal to Cerebellar Vermis VI)
* **RIFJ → RCH6** (Right Inferior Frontal Junction to Right Cerebellum VI)
* **MCV6 → RMFC** (Cerebellar Vermis VI to Right Medial Frontal)

### Short-range Paths (Rule-based S-R Processes):

* **LFEF → RFEF** (Bilateral Frontal Eye-Fields)
* **LCH6 → RCH6** (Bilateral Cerebellar Lobule VI)
* **MCV6 → LTHAL** (Cerebellar Vermis to Left Thalamus)

---

## 📊 Dataset & Region of Interest (ROI) Overview

The analysis incorporates 13 predefined ROIs derived from publicly available processing speed activation maps (NeuroVault:857):

| Label | Anatomical Substrate | Voxels | MNI Coordinates (X, Y, Z) |
| --- | --- | --- | --- |
| **LMFC / RMFC** | Left / Right Medial Frontal Cortex | 155 / 155 | (-5.8, 51.1, 8.3) / (5.2, 49.4, 13.1) |
| **LFEF / RFEF** | Left / Right Frontal Eye-Field | 150 / 148 | (-28.8, -4.4, 51.2) / (35.4, 50.3, -2.5) |
| **LIFJ / RIFJ** | Left / Right Inferior Frontal Junction | 153 / 150 | (-44.8, 2.7, 34.2) / (45.0, 6.7, 32.7) |
| **LIPS / RIPS** | Left / Right Intraparietal Sulcus | 155 / 154 | (-30.8, -50.5, 45.0) / (31.5, -49.5, 45.0) |
| **LTHAL / RTHAL** | Left / Right Thalamus | 153 / 153 | (-11.6, -19.6, 9.2) / (11.9, -17.1, 10.5) |
| **LCH6 / RCH6** | Left / Right Cerebellum Lobule VI | 154 / 150 | (-30.5, -59.7, -25.5) / (24.5, -53.0, -22.1) |
| **MCV6** | Medial Cerebellar Vermis VI | 153 | (2.0, -68.6, -20.1) |

---

## 🛠️ Installation & Environment Setup

### Prerequisites

* **R Version:** `>= 4.0.0`
* **Required R Packages:** `lme4`, `glmmLasso`, `emmeans`, `retime`, `vars`, `ppcor`, `tseries`, `MuMIn`

### Quick Setup

# Clone the repository
git clone [https://github.com/clivehywong/2021CPS.git](https://github.com/clivehywong/2021CPS.git)
cd 2021CPS

# Install required R dependencies
Rscript -e 'install.packages(c("lme4", "glmmLasso", "emmeans", "vars", "ppcor", "tseries", "MuMIn"))'

---

## 🚀 How to Reproduce Analysis

1. **Connectivity Estimation:**
Run the scripts in `Connectivity/` to extract windowed time-series from the raw datasets in `data/` and compute the connectivity matrices (VAR(1), VAR(2), gPPI, and correlation metrics).
2. **Mixed-Effect LASSO Regression:**
Execute the R scripts in `MixedEffectLasso/` to run cross-validated variable selection and reproduce the final mixed-effect models predicting Cognitive Speed.

---


---

## 🖼️ Figures & Model Architecture

### Fig. 1: Experimental Task & Workflow
> Schematic diagram describing the adapted Arrow Task in 3 conditions (compatible, incompatible, and control) crossed with 2 modalities (visual and audial).
![Fig 1: Task Design](https://ars.els-cdn.com/content/image/1-s2.0-S1053811920310417-gr1.jpg)

### Fig. 2: Regions of Interest (ROIs) & Time-Series Processing
> (A) Workflow for extracting ROIs. (B) Selected 13 ROI masks derived from processing speed activation maps. (C) Extraction of task-specific windowed BOLD time-series.
![Fig 2: ROI Extraction](https://ars.els-cdn.com/content/image/1-s2.0-S1053811920310417-gr2.jpg)

### Fig. 3: Model Estimation & Feature Selection Procedure
> **(A)** Overview of model estimation: Task connectivity matrices are preselected via one-sample *t*-tests (*p* < 0.05), followed by Mixed-Effect LASSO for feature selection, and final linear mixed-model fitting.  
> **(B)** Differentiation of positive vs. negative interregional connectivity.  
> **(C)** Alignment of connectivity direction and regression slope for speed-connectivity correlates.
![Fig 3: Model Estimation Procedure]
<img width="652" height="524" alt="1-s2 0-S1053811920310417-gr3" src="https://github.com/user-attachments/assets/d4012d6c-5878-45a1-9898-2455fd5c9a6d" />


### Fig. 4: Significant Fronto-Cerebellar Effective Connectivity Network
> Predictive VAR(1) connectivity network for Cognitive Speed (CS) identified by Mixed-Effect LASSO.
![Fig 4: Connectivity Predictors](https://ars.els-cdn.com/content/image/1-s2.0-S1053811920310417-gr4.jpg)

---

## 💾 Data and Code Availability

- **fMRI Datasets & Transformations:** Raw/ICA-cleaned functional time-series, ANTs transformation matrices, and watershed ROI masks are fully hosted and available in the [`data/`](data/) directory of this repository.
- **ROI Activation Maps:** Public activation maps used to define the 13 ROIs were retrieved from [NeuroVault (Collection 857)](https://identifiers.org/neurovault.collection:857).
- **Execution Scripts:** All analysis scripts for connectivity estimation (`PPI.R`, `funcConn.R`) and mixed-effect modeling (`MMLasso.R`) are included in `Connectivity/` and `MixedEffectLasso/`.

---

## 📖 Citation

If you use this codebase or dataset in your academic research, please cite our paper:

@article{Wong2021FrontoCerebellar,
  title={Fronto-Cerebellar Connectivity Mediating Cognitive Processing Speed},
  author={Wong, Clive H. Y. and Liu, Jiao and Lee, Tatia M. C. and Tao, Jing and Wong, Alex W. K. and Chau, Bolton K. H. and Chen, Lidian and Chan, Chetwyn C. H.},
  journal={NeuroImage},
  volume={226},
  pages={117556},
  year={2021},
  doi={10.1016/j.neuroimage.2020.117556}
}

---

## 📬 Contact & License

* **Corresponding Authors:** Clive H. Y. Wong ([hyewong@gmail.com](https://www.google.com/search?q=mailto%3Ahyewong%40gmail.com)), Chetwyn C. H. Chan?
* **License:** Distributed under the GNU General Public License v3.0 (GPL-3.0)

---

## 👏 Acknowledgments

This study was supported by:
- The 12th Five-Year Plan project of the Ministry of Science and Technology of the People's Republic of China (Grant No. **2013BAI10B01**)
- Fujian Rehabilitation Tech Co-Innovation of China
- The University of Hong Kong May Endowed Professorship in Neuropsychology
- Science and Technology Program of Guangdong (Grant No. **2018B030334001**)

*Special thanks to the University Research Facility in Behavioral and Systems Neuroscience, The Hong Kong Polytechnic University, for providing technical support and scanning equipment throughout the study.*
