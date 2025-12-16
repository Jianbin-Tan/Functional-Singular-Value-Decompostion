# Functional Singular Value Decomposition (FSVD)

This repository provides the data and code associated with the paper "Functional Singular Value Decomposition" by Jianbin Tan, Pixu Shi, and Anru Zhang. The paper is available on [arXiv](https://arxiv.org/abs/2410.03619), and an accompanying R package, **FSVD**, can be found at [https://github.com/Tan-jianbin/FSVD](https://github.com/Tan-jianbin/FSVD).

## 1. Data

This repository contains both real datasets utilized in our study. Specifically: 
- The dynamic COVID-19 dataset was downloaded from the COVID-19 Data Repository by CSSE at Johns Hopkins University (https://github.com/CSSEGISandData/COVID-19).
Specifically, we used the raw files: [time_series_covid19_confirmed_global.csv](https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv) and [UID_ISO_FIPS_LookUp_Table.csv](https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv).

- The longitudinal electronic health records (EHR) were obtained from the [MIMIC-IV Database on PhysioNet](https://physionet.org/content/mimiciv/3.0/).
Specifically, our analysis uses the file: `lab_drg_870_872_Nov_iv.csv`. After obtaining credentialed access and downloading these files, we extracted the EHR data from one patient saved as `dat_ehr.rda`.

### Processing

- **COVID-19 dataset:** We load the two raw files from the JHU CSSE repository, select a fixed set of countries/regions, normalize their regional epidemic counts by the corresponding populations, and apply a log10 transformation. We then create the analysis inputs as (i) the transformed observations and (ii) their corresponding time indices, with time rescaled to [0,1]. The preprocessing code is implemented in `Data_analysis_COVID19.rda`.

- **EHR dataset (MIMIC-IV):** We load the extracted EHR files `lab_drg_870_872_Nov_iv.csv`. For each patient, we subset the records by `SUBJECT_ID` and split the data by clinical feature (`FEATURE_NAME`). Within each feature, we rescale the observed times to [0,1] over that patient’s time window, remove missing measurements, and snap times to a common grid. When multiple measurements fall on the same grid time, we aggregate them by taking the mean value. Finally, we retain only features with sufficient information (at least 5 observed non-zero values) for downstream analysis. The preprocessing code is implemented in `Data_analysis_EHR_whole.rda`.

## 2. Code
### Overview
The code in this repository demonstrates the application of our proposed FSVD methodology to various tasks, such as optimal dimension reduction, functional clustering, functional regression, factor model, and data completion. We benchmark our method against several prominent approaches:
- [Functional principal component analysis (FPCA)](https://cran.r-project.org/web/packages/fdapace/)
- [Identification of substructures in longitudinal data through functional clustering](https://cran.r-project.org/web/packages/fdapace/)
- [Clustering for sparsely sampled functional data](https://www.tandfonline.com/doi/abs/10.1198/016214503000189)
- [Functional linear regression analysis for longitudinal data](https://projecteuclid.org/journals/annals-of-statistics/volume-33/issue-6/Functional-linear-regression-analysis-for-longitudinal-data/10.1214/009053605000000660.full)
- [Penalized functional regression](https://www.tandfonline.com/doi/abs/10.1198/jcgs.2010.10007?casa_token=2eQCx5RtRgYAAAAA:Wrh0wz1Qs2MemK1Q4ysRWUQ1uiop5I4lUOuLnZSlbJQZF4Fqc72Nggw3Cb-lSvxmjUE-MghFS6cHYA)
- [Factor models for high-dimensional time series](https://cran.r-project.org/web/packages/HDTSA/index.html)
- [Matrix completion](https://cran.r-project.org/web/packages/filling/index.html)
- [Variational autoencoders](https://doi.org/10.1016/j.patcog.2020.107501)
- [Predictive methods to missing data imputation](https://www.jmlr.org/papers/v18/17-073.html)
- Singular value decomposition (SVD)
- Smoothing spline

### Reproducibility
- **Simulation Results**: The results presented in Section 7 can be reproduced by running the script `Simulation_result.R`.
- **COVID-19 Dynamic Data Analysis**: The analysis in Section 8 for the COVID-19 data can be performed using `Data_analysis_COVID19.R`.
- **Longitudinal EHR Analysis**: The analysis in Section 8 for the EHR data can be conducted using `Data_analysis_EHR.R`.
