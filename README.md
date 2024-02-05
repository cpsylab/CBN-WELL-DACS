# Code for CANBIND-WELLNESS Depression Anxiety Coupling Strength Prediction Study

## Paper 

Nunes A, Pavlova P, Cunningham JEA, Nuñez JJ, Quilty LC, Foster JA, Harkness KL, Ho K, Lam RW, Li QS, Milev R, Rotzinger S, Soares CN, Taylor VH, Turecki G, Kennedy S, Frey BN, Rudzicz F, Uher R. Depression-Anxiety Coupling Strength as a Predictor of Relapse in Major Depressive Disorder: A CAN-BIND Wellness Monitoring Study Report

## Requirements 

### Software 
- Julia 1.9.2. See `Project.toml` for package requirements. 
- R version 4.2.2. 
    - Packages: 
		- `dplyr`
		- `EValue`
		- `ggplot2`
		- `lmerTest`
		- `perm`
		- `reshape2`
		- `rms`
		- `sjPlot`
		- `survival`
		- `survminer`
		- `tableone`

### Data requirements 
The `data/` directory requires the following files: 
	- `data/20210706_CBN_Wellness_Relapser.xlsx`
	- `data/Wellness.DEID.Subject.List.xlsx`
	- `data/ae.sas7bdat`
	- `data/ce.sas7bdat`
	- `data/cm.sas7bdat`
	- `data/dm.sas7bdat`
	- `data/fa.sas7bdat`
	- `data/ho.sas7bdat`
	- `data/mh.sas7bdat`
	- `data/qs.sas7bdat`
	- `data/sc.sas7bdat`
	- `data/sv.sas7bdat`
	- `Wellness MDE FINAL 20231101.xlsx`

## Running

- Run `main.sh` to create necessary folders in your directory. You will need to have access to the CAN-BIND data to then populate the `data` folder with the necessary files above. 
- Run `julia --project dataio.jl` in the terminal, which will parse all of the data files for the analyses. 
- R Code: 
	- Note: you can run these via the terminal, but some of the formatted tables (generated via `tab_model`) may not be generated. You may need to run them in interactive mode. 
	- `analysis.R`: Runs the main analysis (Relapse-C), as well sensitivity analyses evaluating the time horizon at which DACS predicts Relapse-C, and the number of QIDS/GAD7 measurements required. 
	- `analysis-initial4weeks.R`: Runs the analyis predicting Relapse-C using only the first 4 weeks of available QIDS/GAD7 measurements
	- `analysis-4weekspre.R`: Runs the analyis predicting Relapse-C using only the final 4 weeks of available QIDS/GAD7 measurements
	- `analysis-relapseredef.R`: Runs the analysis predicting Relapse-I
	- `analysis-clearrelapse.R`: Runs the analysis that excludes patients without clear relapses (generates Table S8)