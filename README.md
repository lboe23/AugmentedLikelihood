# AugmentedLikelihood
Code for simulation studies in paper, "An Augmented Likelihood Approach for the Discrete Proportional Hazards Model Using Auxiliary and Validated Outcome Data – with Application to HCHS/SOL Study", Boe and Shaw

Code files include:

(1) R code that runs simulations for data from a simple random sample to assess the performance of the proposed method compared to the the standard interval-censored approach that does not incorporate auxiliary data: Proposed_FinalSims_RandomSample.R

(2)(a) R code modified from Baldoni et al. 2021 that simulates a superpopulation and draws samples from it using a stratified three-stage sampling scheme. Our code combines the files generateTargetPopulation.R, generateTargetPopulationData.R, and sampleTargetPopulation.R  from https://github.com/plbaldoni/HCHSsim/tree/master/R into one code file and modifies the code in order to simulate our covariate and event times of interest and adjust the sample size. Our resulting code file is called GenerateSampleComplexSurveyData.R, which produces 1000 simulated data sets to be used for the simulation in part (b).
(b) R code that runs simulations for complex survey design data to assess the performance of the proposed method compared to the the standard interval-censored approach that does not incorporate auxiliary data using the samples generated in (2)(a): Proposed_FinalSims_ComplexSurvey.R

(3) An R script, PROPOSED_AUGMENTEDLIKELIHOOD_FUNCTIONS.R, containing two functions, "log_like_proposed," which calculates the log-likelihood for the proposed method and "gradient_proposed" which calculates the gradient/estimating function for the proposed method. Both functions ask you to specify function purpose which should either be "SUM" or "INDIVIDUAL." SUM returns the sum of the log-likelihood or gradient contributions, respectively, while INDIVIDUAL returns a vector/matrix of each person's contributions.

(4) An RCPP file creating the function "cmat" which is used in file (2). The function "cmat" is inspired by the function "dmat" described below.

(5) A simulated dataset representing hypothetical data with similar features to the HCHS/SOL cohort that will be used to illustrate our method: SampleData.RData

(6) R code that performs a sample data analysis using the simulated data from (5) above: Sample_Data_Analysis_Final.R. This code illustrates the analysis methods we used in our manuscript to compute hazard ratio (HR) and 95% confidence interval (CI) estimates of a simulated health outcome based on (a) the proposed method and (b) the standard, no auxiliary data method. In this applied analysis, we show how to use the proposed method with complex survey data and additionally apply regression calibration to adjust for covariate error. For variance estimation in the presence of regression calibration and a complex survey design, we will use the parametric multiple imputation (MI) procedure proposed by Baldoni et al. 2021. The code used for this MI-based variance estimation procedure is adapted from the following R code: https://github.com/plbaldoni/HCHSsim/blob/master/R/runMI.R. 

(7) An R data file containing the output obtained from applying this parametric MI procedure for variance estimation: SampleAnalysis_MIVarianceResults.RData. We provide this output so that the user does not have to wait to run the code for the MI variance estimation procedure and simply has access to the results. 

Note: Files (1)-(2)(a)(b) rely on Files (3) and (4) in order to calculate the log likelihood. File (3) relies on functions "dmat" and "getrids" from the following link https://github.com/XiangdongGu/icensmis/blob/master/src/dataproc.cpp in order to run. These Rcpp functions were developed by Gu, Ma, and Balasubramanian (2015) for the following paper: "Semiparametric time to event models in the presence of error-prone, self-reported outcomes—With application to the women’s health initiative."
