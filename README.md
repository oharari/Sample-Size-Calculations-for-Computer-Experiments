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

# About this application

The aim of this application is to provide sample size calculations for the Gaussian process model. For more information read [the manuscript](https://b38efb82-a-62cb3a1a-s-sites.googlegroups.com/site/ofirhararishomepage/Sample_Size_Revisited.pdf?attachauth=ANoY7cocRxLf7yLkbT3uauaHg-e7eHSgED_RRyIGZ74c9y9ub5xYK7_jckKQXYft63DWVOiPcSCDLEIA0RPVUx6ItLhSV9Zo3hio_4EiFy6ayOVN9WfnHU8QeCuEjg814Hz3BI0YW8jnuhZFRLqKpIElpREIIXktVtqqG0ECflxi2VYWMV5SkBI83fInw0lMhs-Op-X7gsuOZccTd8eQD0sR-dOKVGvOrani8ZI04LVnsuDKIW5bxMI%3D&attredirects=0). 

The focus of these calculations is the **_Average Unexplained Variation_** of a predictor $\widehat{y}\left(\boldsymbol{x}\right)$ given a sample $\mathcal{D}=\left\lbrace\boldsymbol{x}_1\ldots,\boldsymbol{x}_n\right\rbrace\subset [0,1]^d$, defined by 

```math
\mathrm{AUV}\left(\widehat{y};\mathcal{D}\right) = \int_{[0,1]^d}\dfrac{\mathbb{V}\mathrm{ar}\left\{y\left(\boldsymbol{x}\right)\big|\mathcal{D}\right\}}{\mathbb{V}\mathrm{ar}\left\{y\left(\boldsymbol{x}\right)\right\}}\hspace{.1em}\mathrm{d}\boldsymbol{x}.
```

We assume that the observed data follow a zero mean, stationary Gaussian process with correlation function $R\left(\boldsymbol{x},\boldsymbol{x}'\right)$. Micchelli & Wahba (1979) and Harari & Steinberg (2014) have shown that if we use the Mercer expansion
```math
R\left(\boldsymbol{x},\boldsymbol{x}'\right) = \sum\limits_{k=1}^{\infty}\lambda_k\varphi_k\left(\boldsymbol{x}\right)\varphi_k\left(\boldsymbol{x}'\right),
```
a lower bound for the AUV is given by
```math
\hspace{20em}\mathrm{AUV}\left(\widehat{y};\mathcal{D}\right)\geq\sum\limits_{k\geq n+1}\lambda_k.\hspace{20em}(1)
```
The program solves (numerically) the Fredholm integral equation of the second kind to derive the Mercer expansion of $R\left(\boldsymbol{x},\boldsymbol{x}'\right)$, and consequently (square root of) the lower bound $(1)$.

The goal of the "Sample Paths" tab is to assist the researcher with the choice of correlation function (and hyperparameters) for his model, by plotting one dimensional realizations.

The families of correlation functions supported by this application are --

1. **The Power Exponential correlation function:**
```math
R\left(\boldsymbol{x},\boldsymbol{x}'\right) = 
\exp\left\{-\sum\limits_{i=1}^d\dfrac{\left|x_i-x_i'\right|^{p_i}}{\theta_i}\right\}
```
where $\theta_i$ and $p_i$ are the correlation length and rougness parameters in the $i$th direction, respectively.

2. **The Gaussian correlation function:**
a special case of the Power Exponential correlation function for $p_1=\cdots=p_d=2$.

3. **The (product) Matern correlation function:**
```math
R\left(\boldsymbol{x},\boldsymbol{x}'\right) = \prod\limits_{i=1}^d\dfrac{1}{\Gamma(\nu)2^{\nu-1}}\left(\dfrac{2\sqrt{\nu}\left|x_i-x_i'\right|}{\theta_i}\right)^{\nu}\mathcal{K}_{\nu}\left(\dfrac{2\sqrt{\nu}\left|x_i-x_i'\right|}{\theta_i}\right)
```
where $\theta_i$ is the correlation length parameter along the $i$th direction, $\nu$ is the smoothness parameter and $\mathcal{K}_{\nu}\left(\cdot\right)$ is the modified
Bessel function of order $\nu$.

