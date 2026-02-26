# oceanTF
The notebooks and code within this folder run the processing to create ocean thermal forcing datasets.  The scripts were initially written by Lizz Ultee (elizabeth.ultee@nasa.gov), with contributions from Donald Slater and Denis Felikson.

# Process
To run the oceanTF workflow, follow these steps:
1. Gather input datasets:
    - Baseline ocean datasets (temperature and salinity): the current notebook is configured to use Hadley EN4
    - CMIP ocean datasets: thetao (temperature) and so (salinity)
2. Run `Step1-TF_process-EN4.ipynb`
3. Run `Step2-Compute-CESM_WACCM_TF.ipynb`
4. Run `Step3-QDM_CESM_WACCM-CCR.ipynb`
5. Run `Step4-To_ISMIP_Grid.ipynb`

Note that you will only run Step1 once.  If processing a new GCM, you will need to re-run steps 2-4, but you will re-use the same climate baseline data (EN4, already processed in Step1 the first time you ran it).