# covidm-mtd-Omi

Data and `R` code to support Barnard *et al.*, Modelling the medium-term dynamics of SARS-CoV-2 transmission in England in the Omicron era, [*Nature Communications*](https://www.nature.com/articles/s41467-022-32404-y) (2022). 

The original pre-print submitted in November 2021 and later versions can be found [here](https://www.medrxiv.org/content/10.1101/2021.11.22.21266584)

[![DOI](https://zenodo.org/badge/511480146.svg)](https://zenodo.org/badge/latestdoi/511480146)

## Study description

England has experienced a heavy burden of COVID-19, with multiple waves of SARS-CoV-2 transmission since early 2020 and high infection levels following the emergence and spread of Omicron variants since late 2021. In response to rising Omicron cases, booster vaccinations were accelerated and offered to all adults in England. Using a model fitted to more than 2 years of epidemiological data, we project potential dynamics of SARS-CoV-2 infections, hospital admissions and deaths in England to December 2022. We consider key uncertainties including future behavioural change and waning immunity and assess the effectiveness of booster vaccinations in mitigating SARS-CoV-2 disease burden between October 2021 and December 2022. If no new variants emerge, SARS-CoV-2 transmission is expected to decline, with low levels remaining in the coming months. The extent to which projected SARS-CoV-2 transmission resurges later in 2022 depends largely on assumptions around waning immunity and to some extent, behaviour, and seasonality.

## Main repository files

A description of how these results were generated using the files in this repository follows below.

### Data processing

The following `R` scripts are used to process data prior to model fitting:

#### Data input 1 (e.g. to produce `processed-data-20220506122858.qs`)

1. `build_sitrep_data.R`: `R` script to process NHS England sitrep data on COVID-19 deaths, hospital admissions and hospital and ICU bed occupancy (sensitive data provided as part of SPI-M membership, not provided in this repository)
2. `build_virusprev_data.R`: `R` script to process SARS-CoV-2 PCR prevalence data for NHS England regions (data publicly available from the [Office for National Statistics Coronavirus Infection Survey](https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata))
3. `build_seroprev_data.R`: `R` script to process SARS-CoV-2 seroprevalence data for NHS England regions (data publicly available from three sources: the [Office for National Statistics Coronavirus Infection Survey](https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk), the [UK Biobank](https://www.ukbiobank.ac.uk/learn-more-about-uk-biobank/covid-19-hub) and Imperial College London's [REACT-2 study](https://www.imperial.ac.uk/medicine/research-and-impact/groups/react-study/))
4. `build_delta_data.R`: `R` script to process genomic sequencing data for the Delta B.1.617.2 variant of concern in NHS England regions (sensitive data provided as part of SPI-M membership, not provided in this repository)
5. `build_sgtf_data.R`: `R` script to process S gene target failure data used as a proxy for the proportion of infections attributable to the Alpha B.1.1.7 and Omicron B.1.1.529 variants of concern in NHS England regions (sensitive data provided as part of SPI-M membership, not provided in this repository)
6. `build_processed_data_for_fitting.R` `R` script which takes data outputs from the previous 5 scripts and builds a single `.qs` file for model fitting beginning `processed-data-`

#### Data input 2 (e.g. to produce `schedule3-MTPs-20220506121302.rds`)

 1. `build_schoolattendance.R`: `R` script to process data on school attendance in England (data publicly available from the UK government's [Department for Education statistics](https://explore-education-statistics.service.gov.uk/find-statistics/attendance-in-education-and-early-years-settings-during-the-coronavirus-covid-19-outbreak))
 2. `build_schedule_google_mobility.R`: `R` script to combine school attendance data from the script above with publicly available [Google COVID-19 Community Mobility Reports](https://www.google.com/covid19/mobility/) data to build a single `.RDS` schedule file for model fitting, beginning 'schedule3-MTPs-'
 
#### Data input 3 (e.g. to produce `vax-covidm20220505205235.rds`)
 
 1. `build_vax_data.R`: `R` script to process COVID-19 immunisations data and to generate vaccination schedules recording the number of individuals receiving vaccinations by vaccine product, date, age and NHS England region
 

### Model fitting (two-stage)

#### MCMC fitting

We use `fit-omi.R` to generate MCMC fits. We call this script from the command line with the following arguments:

`Rscript fit-omi.R FIT_TYPE POP_SET REP_START REP_END WANE_YN VAC_EFF V3_SEVERITY DFR SEAS_YN SEAS_AMP START_FILE BFOLD OMI_PROTECTION OMI_SEV OMI_CRIT BA2_RELU GAM_DISP`

For example, we use the following inputs for model fits in the manuscript:

`Rscript fit-omi.R relu all 1 5 yeswane central 2.0 1.0 seaslate 0.1 relu_yeswane_sev2.0_22021104 2.5 0.551 0.5 0.5 1.5 0.3`

#### Particle filter fitting
    
We use `pfilter.R` to generate particle-filtered fits from MCMC fits. We call this script from the command line with the following arguments:

`Rscript pfilter.R FIT_FILE`

For example, we use the following inputs for model fits in the manuscript:

`Rscript pfilter.R relu_yeswane_sev2.0_22050605`

### Generating mobility schedules (for projections)

We use code within `paper_mobility.R` to generate mobility scenarios to use as inputs for the projections

### Generating sample fits (for projections)

We use `make_projections.R` to generate model projections from model fits. We call this script from the command line with the following arguments:

`Rscript make_projections.R`

### Producing Figures & Tables

We use the following scripts to produce Figure files contained in the manuscript:

1. `vaccine_plotting.R` to generate processed files comparing vaccine schedules (for Figure S7)
2. `paper_figs_updated.R` (Figures 2, 3, S1A, S1B, S1C, S2, S3, S4A, S4B, S5-S8)
3. `plot_story_fig_updated.R` (Figure 4)
4. `plot_summary_fig.R` (Figure 5)
5. `tableS1.R` (Table S1)

Figure 1 is generated in these [Google slides](https://docs.google.com/presentation/d/1eaaWGhDxt3lVmOtWg7yNLzNU9bs9Vne_pLf6tTHlqpY/edit?usp=sharing)

### Other scripts & subfolders

Some scripts in the main repository are not directly used but are sourced from within other scripts:

* `booster_schedule.R` is sourced from `fit-omi.R`, `pfilter.R` and `make_projections.R` and is used to generate first booster vaccine timings for each 5-year age group
* `check_fit.R` is sourced from `fit-omi.R`, `pfilter.R` and `make_projections.R` and is used to check the model fit
* `commit.R` is sourced from `build_sitrep_data.R` and is used to compare new data builds to older data builds
* `cpp_funcs.R` is sourced from `fit-omi.R`, `pfilter.R` and `make_projections.R` and generates C++ functions
* `khoury_lookup.R` is sourced from `cpp_funcs.R` and `params.R` and is used to convert vaccine effectiveness assumptions assuming changes in neutralisation fold
* `params.R` is sourced from `fit-omi.R`, `pfilter.R` and `make_projections.R` and sets up model parameters
* `processes.R` is sourced from `fit-omi.R`, `pfilter.R` and `make_projections.R` and sets up burden processes
* `spim_output.R` is sourced from `build_processed_data_for_fitting.R`, `fit-omi.R`, `pfilter.R`, `paper_figs_updated.R`, `plot_story_fig_updated.R` and `make_projections.R` and generates sample model fits and projections
* `vax_analysis.R` was used to generate summary statistics based on the rollout of COVID-19 vaccines in 2021 (e.g. the typical delays between first and second doses)
* `vax_funcs.R` is sourced from `build_vax_data.R` and `vax_analysis.R` and contains functions for processing vaccinations data

Other subfolders:

* `./covidm_for_fitting/` contains the covidm model itself (we are using version 3, contained in `covidm_for_fitting/model_v3/`)
* `./scripts/` contains `fitmodel.sh` which is a shell script used to run the MCMC fitting (i.e. `fit-omi.R`) from the command line
* `./data/` contains data files required to run the model (e.g. information on demographics and NHS England regions) and data required to generate model assumptions (e.g. on booster vaccination uptake and on school attendance)
* `./fits/` is where model fits are saved and read from
* `./fitting_data/` usually contains all of the data required to run the model fitting (this data is not publicly available and thus is not provided in this repository)
* `./output/` is where model output files get saved

### Other details & contact

The analysis was performed using `R` version 4.0.3 (2020-10-10).

For any issues with the code please contact [Rosanna Barnard](https://www.lshtm.ac.uk/aboutus/people/barnard.rosanna).
