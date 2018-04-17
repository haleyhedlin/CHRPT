# CHRPT
This repository contains the code to develop and validate the risk prediction tool described in the manuscript 'Development and Validation of a Comprehensive Health Risk Prediction Tool for Post-Menopausal Women' with the hopes that others will apply these methods to new cohorts. The code used to create the web app using Shiny is also provided.

model_build.R contains the R code for model building and internal validation code for the event-first model

model_build_ee.R contains the R code for model building and internal validation code for the event-ever model

predict_ext.R contains the R code for obtaining predictions from model estimates based on the WHI models when applied to new cohorts

val.csv gives the data structure for the data used in predict_ext.R

validation_data_dictionary.xls is a data dictionary for val.csv

app.R contains the code to build the web app in Shiny
