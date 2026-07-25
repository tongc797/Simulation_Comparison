This repository contains all code to replicate simulation results presented in "*Evaluating the Performance of Structural Equation Models: Process-Based vs. Covariance-Based Simulation*"

Scripts are organized in three folders:<br>

**1. Simulate Data & Fit Model:** contains code to simulate data using process-based or covariance-based simulation, and code for using the simulated data to fit the CLPM or the CoT model.<br>
* **CLPM:** this folder contains code for the CLPM<br>
**CLPM_FT_generate_data.R:** generate data using process-based simulation<br>
				**CLPM_MVN_generate_data.R:** generate data using covariance-based simulation<br>
				**CLPM_FT_FitModel.R:** fit the CLPM using data generated from process-based simulation<br>
				**CLPM_MVN_FitModel.R:** fit the CLPM using data generated from covariance-based simulation<br>
* **CoT:** this folder contains code for the CoT model<br>
				**CoT_GE_generate_data.R:** generate data using process-based simulation, needs to be run with GeneEvolve scripts, use this script to replace the GE-75.R script in the GeneEvolve repository<br>
				**CoT_GE_FitModel.R:** fit the CoT model using data generated from process-based simulation<br>
				**CoT_MVN_generate_data_FitModel.R:** this code does two things: first generates data using covariance-based simulation and then fits the CoT model using the data generated

**2. True Values:** contains code to obtain true values of parameters in different models (CLPM, CoT).<br>
* **CLPM_W18_TrueValues.R:** get true values of the variances and covariance of X and Y at Wave 18 for the CLPM<br>
* **CoT_GE_TrueValues.R:** get empirical observed values of model parameters from GeneEvolve simulated datasets to be used as "true values"<br>
* **CoT_nonlinear constraints.R:** get true values of parameters in the CoT model that involve nonlinear constraints
		
**3. Obtain results:** contains code for obtaining various results<br>
* **Calculate bias precision power from results.R:** get bias, precision, and statistical power of all parameters from model fitting results<br>
* **Test significance of bias difference.R:** test whether differences in biases obtained from the two simulation approaches are statistically significant<br>
* **density plots of bias.R:** create density plots of bias (Fig. 2 & Fig. S1-4)
