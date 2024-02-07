# Sample-Size-Calculations-for-Computer-Experiments
Code for our 2018 Statistica Sinica paper "Computer Experiments: Prediction Accuracy, Sample Size and Model Complexity Revisited".


This repository contains the following:

**********************************************************************************************
** File Piston_RAUV_Simulation.R
**********************************************************************************************

	- R code to reproduce the results of section 6 in our paper ("An Illustration").

	- Fits model to a size 1000 design to estimate the variance parameter sigma^2.
	
	- Reproduces the simulation study, summarized in Table 1.

	- Reproduces Figure 6 and Figure 7.


**********************************************************************************************
** Sub-directory Sample_Size_Shiny
**********************************************************************************************

	- An R-Shiny application for sample size calculations

	- Could be used to recreate Figures 2, 3, 4 and 8

	- provides the results documented in Example 3.1.1, Example 3.1.2, Section 3.3, 
	  the middle column of Table 1, Table 2 and the 3rd paragraph in the Discussion 
	  section.

	- Can be run using an R Studio client, by either executing the "server.R" or the 
	  "ui.R" file.

	- Make sure to save the directory locally as is, including the .rmd and the .html 
          files.

**********************************************************************************************
**********************************************************************************************
