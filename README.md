# Functional Singular Value Decomposition (FSVD)

This README accompanies the paper "Functional Singular Value Decomposition" by Jianbin Tan and Anru Zhang, detailing the associated data and code.

## 1. Data
### Abstract

This repository contains both simulated and actual datasets utilized in our study. The dynamic COVID-19 dataset is publicly available on the [COVID-19 Data Repository by CSSE at Johns Hopkins University](https://github.com/CSSEGISandData/COVID-19). The longitudinal electronic health records are sourced from the [MIMIC-IV Database](https://physionet.org/content/mimiciv/3.0/).

### Availability
All datasets necessary to replicate our findings are provided.

### Data Dictionary
Included within the "Data" directory are the raw data files for COVID-19 dynamic trajectories and longitudinal electronic health records, along with additional requisite information.

## 2. Code
### Abstract
Our proposed FSVD methodology is applied to tasks including optimal dimension reduction, clustering, factor modeling, and functional completion. We benchmark our method against several prominent techniques:
- [Functional Principal Component Analysis (FPCA)](https://cran.r-project.org/web/packages/fdapace/)
- [Identification of substructures in longitudinal data through functional clustering](https://cran.r-project.org/web/packages/fdapace/)
- [Clustering for sparsely sampled functional data](https://www.tandfonline.com/doi/abs/10.1198/016214503000189)
- [Factor models for high-dimensional time series](https://cran.r-project.org/web/packages/HDTSA/index.html)
- Singular value decomposition (SVD)
- Techniques for matrix completion [available here](https://cran.r-project.org/web/packages/filling/index.html)
- [K-Nearest Neighbors (K-NN)](https://www.jmlr.org/papers/v18/17-073.html)
- Smoothing spline techniques

### Reproducibility
- The results in Section 5 are reproducible by running "Simulation_result.R".
- Analysis of the COVID-19 dynamic data (Section 6.1) is performed using "Data_analysis_COVID19.R".
- Analysis of longitudinal electronic health records (Section 6.2) is conducted using "Data_analysis_EHR.R".