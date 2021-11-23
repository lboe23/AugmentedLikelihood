# AugmentedLikelihood
Code for simulation studies in paper, "An Augmented Likelihood Approach for the Discrete Proportional Hazards Model Using Auxiliary and Validated Outcome Data – with Application to HCHS/SOL Study", Boe and Shaw

Code files include:

(1) R code that runs simulations to assess the performance of the proposed method compared to the the standard interval-censored approach that does not incorporate auxiliary data : Proposed_FinalSims_RandomSample.R

(2) An R script, PROPOSED_AUGMENTEDLIKELIHOOD_FUNCTIONS.R, containing two functions, "log_like_proposed," which calculates the log-likelihood for the proposed method and "gradient_proposed" which calculates the gradient/estimating function for the proposed method. Both functions ask you to specify function purpose which should either be "SUM" or "INDIVIDUAL." SUM returns the sum of the log-likelihood or gradient contributions, respectively, while INDIVIDUAL returns a vector/matrix of each person's contributions.

(3) An RCPP file creating the function "cmat" which is used in file (2). The function "cmat" is inspired by the function "dmat" described below.

File (1) relies on Files (2) and (3) in order to calculate the log likelihood. File (2) relies on functions "dmat" and "getrids" from the following link https://github.com/XiangdongGu/icensmis/blob/master/src/dataproc.cpp in order to run. These Rcpp functions were developed by Gu, Ma, and Balasubramanian (2015) for the following paper: "Semiparametric time to event models in the presence of error-prone, self-reported outcomes—With application to the women’s health initiative."
