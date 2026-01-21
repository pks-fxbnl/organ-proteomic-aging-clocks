Demo code for proteomic aging model (LASSO + bootstrap ensemble)

Summary
This package contains a fully runnable demo script designed to illustrate the workflow used in our manuscript for developing proteomic aging models (conventional age and organ-specific age). To comply with journal policies and protect participant privacy, the demo uses synthetic proteomics data with the same dimensional structure as the real dataset but no real protein names or biological information.

The code demonstrates:
	1	Generation of demo proteomics datasets (training: n=98 healthy; testing: n=1850 HIV)
	2	Training a chronological age predictor using all 1254 proteins
	3	Training an organ-specific age model using a randomly selected subset of 30 proteins
	4	Ensemble prediction using a bootstrap-based lambda aggregation strategy (500 bootstrap models)
	5	Output of age predictions for the test dataset
No real biological data are included.


File list
README.txt                     <-- this file
agingclock_demo_code.R         <-- runnable demo script


Requirements
	•	R version 4.4.2
	•	glmnet 4.1.8


To install dependencies:

install.packages("glmnet")


The script will:
	1	Generate a synthetic proteomic dataset (1254 features, standardized)
	2	Train a LASSO model using bootstrap sampling and lambda aggregation
	3	Output predicted proteomic ages for the test dataset

Once finished, the script returns the following R objects:
Object                      Description
conventional_age_demo       Proteomic age predictions using all 1254 proteins
organ_specific_age_demo     Proteomic age predictions using 30 organ-specific proteins (randomly selected for demo)
conv_res                    Full model training results for conventional model: lambdas, coefficients, loss table, bootstrap predictions
organ_res                   Full model training results for organ-specific model
