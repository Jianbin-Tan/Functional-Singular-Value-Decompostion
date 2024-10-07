# Functional Singular Value Decomposition (FSVD)

This repository contains the data and code scripts for the paper "Functional Singular Value Decomposition" by Jianbin Tan, Pixu Shi, and Anru R. Zhang. This paper is accessible at [Arxiv](https://arxiv.org/abs/2410.03619).

## 1. Data

This repository contains both simulated and real datasets utilized in our study. All datasets necessary to replicate our findings are provided. Specifically: 
- The dynamic COVID-19 dataset is publicly available at the [COVID-19 Data Repository by CSSE at Johns Hopkins University](https://github.com/CSSEGISandData/COVID-19). 
- The longitudinal electronic health records are sourced from the [MIMIC-IV Database](https://physionet.org/content/mimiciv/3.0/).

### Data Dictionary
The "Data" directory contains the raw data files, including:
- COVID-19 dynamic trajectories
- Longitudinal EHR data
Additional information and documentation are also provided for these datasets.

## 2. Code
### Overview
The code in this repository demonstrates the application of our proposed FSVD methodology to various tasks, such as optimal dimension reduction, clustering, factor modeling, and functional completion. We benchmark our method against several prominent approaches:
- [Functional principal component analysis (FPCA)](https://cran.r-project.org/web/packages/fdapace/)
- [Identification of substructures in longitudinal data through functional clustering](https://cran.r-project.org/web/packages/fdapace/)
- [Clustering for sparsely sampled functional data](https://www.tandfonline.com/doi/abs/10.1198/016214503000189)
- [Factor models for high-dimensional time series](https://cran.r-project.org/web/packages/HDTSA/index.html)
- [Matrix completion](https://cran.r-project.org/web/packages/filling/index.html)
- [Predictive methods to missing data imputation](https://www.jmlr.org/papers/v18/17-073.html)
- Singular value decomposition (SVD)
- Smoothing spline

### Reproducibility
- **Simulation Results**: The results presented in Section 5 can be reproduced by running the script `Simulation_result.R`.
- **COVID-19 Dynamic Data Analysis**: The analysis in Section 6.1 can be performed using `Data_analysis_COVID19.R`.
- **Longitudinal EHR Analysis**: The analysis in Section 6.2 can be conducted using `Data_analysis_EHR.R`.
