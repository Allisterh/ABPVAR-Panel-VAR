# ABPVAR
Bayesian mixed-frequency panel VAR

This contains all the R codes used for implementing the methods proposed in the paper 'BAYESIAN GROUP-SHRINKAGE BASED ESTIMATION FOR PANEL
 VECTOR AUTOREGRESSIVE MODELS WITH MIXED FREQUENCY DATA' by Chakraborty, Khare and Michailidis.

Details on the key files:

functions.R

Contains all the R functions related to the estimation and forecasting by the proposed ABPVAR model.

truth_gen_sim_gr1.R and truth_gen_sim_gr2.R

Generates the true VAR transition matrices for simulations using data generated according to Grouping 1 and Grouping 2 of the ABPVAR model, respectively.

sim_main_gr1.R and sim_main_gr2.R

Performs estimation and prediction by ABPVAR model with Grouping 1 and 2 respectively, used in simulations studies.

data_prep_states.R

Code for preparing data on multiple US states used in our macroeconomic data analysis using identical frequency.

forecast_states_gr1.R,forecast_states_gr2.R

Code for forecasting analysis by ABPVAR model with Grouping 1 and 2 respectively, using multiple US states data .

data_prep_EU.R

Code for preparing data on European countries used in our macroeconomic data analysis using mixed frequency.

forecast_EU_gr1.R,forecast_EU_gr2.R

Code for forecasting analysis by ABPVAR model with Grouping 1 and 2 respectively, using data on European countries.

