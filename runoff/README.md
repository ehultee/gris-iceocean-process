# runoff
The notebooks and code within this folder run the processing to create runoff datasets. This workflow performs several steps:
1. Aggregate runoff datasets (baseline and CMIP) by drainage basin
2. Bias correct the CMIP runoff dataset to a baseline dataset

# Process
To run the runoff workflow, follow these steps:
1. Gather input datasets:
    - Baseline runoff dataset: the current notebook is configured to use RACMO
    - CMIP runoff dataset
3. Run `Step1-Process-RACMO.ipynb`
4. Run `Step2-CProcess-CESM_WACCM.ipynb`
5. Run `Step3-Bias-Correct-CESM_WACCM.ipynb`
6. Run `Step4-Final-Adjustments-CESM_WACCM.ipynb`
