# Functional Singular Value Decomposition (FSVD)

This repository provides the data and code associated with the paper "Functional Singular Value Decomposition" by Jianbin Tan, Pixu Shi, and Anru Zhang. The paper is available on [arXiv](https://arxiv.org/abs/2410.03619), and an accompanying R package, **FSVD**, can be found at [https://github.com/Tan-jianbin/FSVD](https://github.com/Tan-jianbin/FSVD).

## 1. Data

This repository contains both simulated and real datasets utilized in our study. Specifically: 
- The dynamic COVID-19 dataset was downloaded from the COVID-19 Data Repository by CSSE at Johns Hopkins University (https://github.com/CSSEGISandData/COVID-19).
Specifically, we used the raw files: [time_series_covid19_confirmed_global.csv](https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv) and [UID_ISO_FIPS_LookUp_Table.csv](https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv).

- The longitudinal electronic health records were obtained from the MIMIC-IV Database on PhysioNet (https://physionet.org/content/mimiciv/3.0/).
Specifically, our analysis uses the extracted files death_drg_870_872_Nov_iv.csv and lab_drg_870_872_Nov_iv.csv. After obtaining credentialed access and downloading the official MIMIC-IV release from PhysioNet, we extracted and preprocessed the relevant tables and saved the processed dataset as `dat_ehr.rda`.

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
