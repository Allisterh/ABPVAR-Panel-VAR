# ABPVAR
Bayesian mixed-frequency panel VAR

This contains all the R codes used for implementing the methods proposed in the paper 'BAYESIAN GROUP-SHRINKAGE BASED ESTIMATION FOR PANEL
 VECTOR AUTOREGRESSIVE MODELS WITH MIXED FREQUENCY DATA' by Chakraborty, Khare and Michailidis.

Details on the key files:
1.functions.R: Contains all the R functions related to estimation and forecasting by the proposed ABPVAR model with mixed frequency data.

2.function_states.R: Contains all the R functions related to estimation and forecasting by the proposed ABPVAR model in an identical frequency setting.

3.sim_gr1.R: Generates true VAR transition matrices, simulates data according to ABPVAR Grouping 1 model and performs estimation and forecasting on the simulated data using ABPVAR model with Grouping 1.

4.sim_gr2.R: Generates true VAR transition matrices, simulates data according to ABPVAR Grouping 2 model and performs estimation and forecasting on the simulated data using ABPVAR model with Grouping 2.

5.data_US_states_gr1.R and data_US_states_gr2.R: Prepares data for analysis of macroeconomic data for multiple US states used in the paper, for application of the ABPVAR model with Grouping 1 and 2 respectively.

6.forecast_states_gr1.R and forecast_states_gr2.R: Performs forecasting exercise using ABPVAR model with Grouping 1 and 2 respectively, based on datasets prepared by the R scripts 'data_US_states_gr1.R' and 'data_US_states_gr2.R'.

7.data_US_states_gr1_2nd_vintage.R and data_US_states_gr2_2nd_vintage.R: Prepares data for analysis of macroeconomic data for multiple US states used in the paper, using the second vintage of the most recent data, for application of the ABPVAR model with Grouping 1 and 2 respectively.

8.data_EU_gr1.R and data_EU_gr2.R: Prepares mixed frequency data for analysis of macroeconomic data for all sets of European countries used in the paper, for application of the ABPVAR model with Grouping 1 and 2 respectively.

9.forecast_EU_g1.R and forecast_EU_g2.R: Performs forecasting exercise using mixed frequency ABPVAR model with Grouping 1 and 2 respectively, based on a dataset prepared by the R scripts data_EU_gr1.R and data_EU_gr2.R.

10.data_EU_gr1_2nd_vintage.R and data_EU_gr2_2nd_vintage.R: Prepares mixed frequency data for analysis of macroeconomic data for a set of European countries used in the paper, using the second vintage of the most recent data, for application of the ABPVAR model with Grouping 1 and 2 respectively. The forecasting exercise using ABPVAR model with Grouping 1 and 2 can be performed on these datasets using the R scripts 'forecast_EU_g1.R' and 'forecast_EU_g2.R' respectively, in a similar way.

11.sim_missp.R: Generates true VAR transition matrices, simulates data and performs analysis using ABPVAR model to assess the effect of misspecification on its forecasting performance, as described in the supplement.

12.Sim_Theta_mat.R: Generates true VAR transition matrices, simulates data and performs analysis using ABPVAR model to evaluate the role of the similarity matrix, as described in the paper.
