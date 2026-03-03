# Subglacial runoff/discharge
The notebooks and code within this folder run the processing to create the subglacial discharge/runoff datasets for the Greenland ocean forcing in ISMIP7. The scripts were initially written by Donald Slater (donald.slater@ed.ac.uk) and subsequently cleaned up by Denis Felikson and Donald Slater.

## Required inputs
1. The ISMIP7 runoff basins definition `ismip7_runoff_basins_nearestbelowsl.nc`.
    - Included in this repo but also on Globus at /ISMIP6/ISMIP7_Prep/CMIP6_test_protocol/Tools/ismip7-gris-ocean-forcing/ismip7_runoff_basins_nearestbelowsl.nc
2. EITHER monthly RACMO surface runoff for the period over which we do the bias correction (currently 1985-2014), OR the aggregated subglacial discharge on the basis of RACMO (see further explanation below on whether you need to run Step 1 of the processing).
    - The aggregated subglacial discharge file is on Globus at /ISMIP6/ISMIP7_Prep/CMIP6_test_protocol/Tools/ismip7-gris-ocean-forcing/Q_RACMO_1958_2024.nc
3. Monthly surface runoff for the CMIP model we are processing.
    - E.g., for CESM2-WACCM these are on Globus at /ISMIP6/ISMIP7_Prep/CMIP6_test_protocol/CESM2-WACCM/ssp585/PROCESSEDmon/runoff/SDBN1/ etc
4. You will need to adjust the directory and file path definitions in the script for the system you are running on.

## Steps
1. Run `Step1-Process-RACMO.ipynb` if necessary
    - Sums past RACMO surface runoff over the ISMIP7 runoff basis to obtain subglacial runoff/discharge for the past decades. Outputs the aggregated subglacial discharge which is subsequently used for bias correction. The aggregated subglacial discharge can be found on Globus at /ISMIP6/ISMIP7_Prep/CMIP6_test_protocol/Tools/ismip7-gris-ocean-forcing/Q_RACMO_1958_2024.nc and this step doesn't need to be run if you already have that file.
2. Run `Step2-Process-CMIP.ipynb`
    - Aggregates CMIP-derived surface runoff over the ISMIP7 runoff basins to obtain CMIP-derived subglacial runoff/discharge.
3. Run `Step3-Bias-Correct-CMIP.ipynb`
    - Bias-corrects the CMIP-derived subglacial runoff using the RACMO-derived subglacial runoff.
4. Run `Step4-Final-Adjustments-CMIP.ipynb`
    - Produces final files with the precise formatting required.

## Possible improvements
There is plenty of scope for making the code more robust and efficient. It would be very possible to combine steps 2-4 into one script so that definitions (e.g., of directories) are only made once, and potentially avoiding writing out and reading intermediate steps, which is currently time consuming.

## Likely problems
1. The workflow has only been run on CESM2-WACCM output. Applying to a different model may lead to file path issues (due to the different file string), or spatial or time grid issues if they are formatted differently from CESM2-WACCM.
2. The definition of directories in the present code is not efficient and could lead to not finding the right files if applied to a different model or run on a different system.