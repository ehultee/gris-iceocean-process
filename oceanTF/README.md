# oceanTF
The notebooks and code within this folder run the processing to create ocean thermal forcing datasets.

# Process
To run the oceanTF workflow, follow these steps:
1. Gather input datasets:
    - Baseline ocean datasets (temperature and salinity): the current notebook is configured to use Hadley EN4
    - CMIP ocean datasets: thetao (temperature) and so (salinity)
3. Run `Step1-TF_process-EN4.ipynb`
4. Run `Step2-Compute-CESM_WACCM_TF.ipynb`
5. Run `Step3-QDM_CESM_WACCM-CCR.ipynb`
6. Run `Step4-To_ISMIP_Grid.ipynb`
