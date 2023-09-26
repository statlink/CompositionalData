---
name: CompositionalData
topic: Compositional Data Analysis
maintainer: Michail Tsagris
email: mtsagris@uoc.gr
version: 2023-09-26
source: https://github.com/cran-task-views/CompositionalData/
---

*Compositional data are positive multivariate data where the sum of the values of each vector sums to the same constant, typically taken to be unity for convenience purposes.*

This task view is a collection of packages intended to supply R code to deal with 
Compositional data analysis. The packages can be roughly structured into the following topics
(although several of them have functionalities which go across these
categories):


### General purpose packages

-   `r pkg("Compositional")`: The package performs regression, classification, contour plots, hypothesis testing and fitting of distributions with or without transformations. The packages further includes functions for percentages (or proportions). 

-   `r pkg("compositions")`: Provides functions for the consistent analysis of compositional data (e.g. portions of substances) and positive numbers (e.g. concentrations).

-   `r pkg("coda.base")`: A minimum set of functions to perform compositional data analysis using the log-ratio approach introduced by Aitchison (1982). Main functions have been implemented in C++ for better performance.

-   `r pkg("easyCODA")`: Univariate and multivariate methods for compositional data analysis, based on log-ratios. Accent is given to simple pairwise log-ratios. Selection can be made of log-ratios that account for a maximum percentage of log-ratio variance. Various multivariate analyses of log-ratios are included in the package. 

-   `r pkg("ToolsForCoDa")`: Provides functions for multivariate analysis with compositional data. Includes a function for doing compositional canonical correlation analysis. This analysis requires two data matrices of compositions, which can be adequately transformed and used as entries in a specialized program for canonical correlation analysis, that is able to deal with singular covariance matrices. A function for log-ratio principal component analysis with condition number computations has been added to the package.


### Robust methods for compositional data analysis

-   `r pkg("complmrob")`: Robust regression methods for compositional data. The distribution of the estimates can be approximated with various bootstrap methods. These bootstrap methods are available for the compositional as well as for standard robust regression estimates. This allows for direct comparison between them.

-   `r pkg("robCompositions")`: Methods for analysis of compositional data including robust methods, imputation of missing values, methods to replace rounded zeros, count zeros, methods to deal with essential zeros, (robust) outlier detection for compositional data, (robust) principal component analysis for compositional data, (robust) factor analysis for compositional data, (robust) discriminant analysis for compositional data (Fisher rule), robust regression with compositional predictors, functional data analysis and p-splines, contingency and compositional tables and (robust) Anderson-Darling normality tests for compositional data as well as popular log-ratio transformations. In addition, visualisation and diagnostic tools are implemented as well as high and low-level plot functions for the ternary diagram.

-   `r pkg("robregcc")`: The package implements the algorithm estimating the parameters of the robust regression model with compositional covariates. The model simultaneously treats outliers and provides reliable parameter estimates. 


### Regression modelling

-   `r pkg("codalm")`: Implements the expectation-maximization algorithm for transformation-free linear regression for compositional outcomes and predictors.

-   `r pkg("codaredistlm")`: Provided data containing an outcome variable, compositional variables and additional covariates (optional); linearly regress the outcome variable on an isometric log ratio transformation of the linearly dependent compositional variables. The package provides predictions (with confidence intervals) in the change in the outcome/response variable based on the multiple linear regression model and evenly spaced reallocations of the compositional values.

-   `r pkg("Compack")`: Regression methodologies with compositional covariates, including (1) sparse log-contrast regression with compositional covariates, and (2) sparse log-contrast regression with functional compositional predictors.

-   `r pkg("FLORAL")`: Log-ratio lasso regression for continuous, binary, and survival outcomes with compositional features.

-   `r pkg("ocomposition")`: Regression model where the response variable is a rank-indexed compositional vector (non-negative values that sum up to one and are ordered from the largest to the smallest). Parameters are estimated in the Bayesian framework using MCMC methods.


### Bioiformatics related packages

-   `r pkg("ampir")`: A toolkit to predict antimicrobial peptides from protein sequences on a genome-wide scale. It incorporates two support vector machine models ("precursor" and "mature") trained on publicly available antimicrobial peptide data using calculated physico-chemical and compositional sequence properties. In order to support genome-wide analyses, these models are designed to accept any type of protein as input and calculation of compositional properties has been optimised for high-throughput use. 

-   `r pkg("ArArRedux")`: Processes noble gas mass spectrometer data to determine the isotopic composition of argon (comprised of Ar36, Ar37, Ar38, Ar39 and Ar40) released from neutron-irradiated potassium-bearing minerals. Then uses these compositions to calculate precise and accurate geochronological ages for multiple samples as well as the covariances between them. Error propagation is done in matrix form, which jointly treats all samples and all isotopes simultaneously at every step of the data reduction process. Includes methods for regression of the time-resolved mass spectrometer signals to t=0 ('time zero') for both single- and multi-collector instruments, blank correction, mass fractionation correction, detector intercalibration, decay corrections, interference corrections, interpolation of the irradiation parameter between neutron fluence monitors, and (weighted mean) age calculation. All operations are performed on the logs of the ratios between the different argon isotopes so as to properly treat them as 'compositional data'.

-   `r pkg("BRACoD.R")`: The goal of this method is to identify associations between bacteria and an environmental variable in 16S or other compositional data. The environmental variable is any variable which is measure for each microbiome sample, for example, a butyrate measurement paired with every sample in the data. Microbiome data is compositional, meaning that the total abundance of each sample sums to 1, and this introduces severe statistical distortions. This method takes a Bayesian approach to correcting for these statistical distortions, in which the total abundance is treated as an unknown variable. This package runs the python implementation using reticulate.

-   `r pkg("coda4microbiome")`: Functions for microbiome data analysis that take into account its compositional nature. Performs variable selection through penalized regression for both, cross-sectional and longitudinal studies, and for binary and continuous outcomes. 

-   `r pkg("codacore")`: In the context of high-throughput genetic data, the package identifies a set of sparse biomarkers that are predictive of a response variable of interest. More generally, the package can be applied to any regression problem where the independent variable is compositional, to derive a set of scale-invariant log-ratios that are maximally associated to a dependent variable.

-   `r pkg("MicrobiomeStat")`: A suite of methods for powerful and robust microbiome data analysis addressing zero-inflation, phylogenetic structure and compositional effects. The methods can be applied to the analysis of other (high-dimensional) compositional data arising from sequencing experiments.

-   `r pkg("seq2R")`: This package is useful for loading '.fasta' or '.gbk' files, and for retrieving sequences from 'GenBank' dataset. This package allows to detect differences or asymmetries based on nucleotide composition by using local linear kernel smoothers. Also, it is possible to draw inference about critical points (i. e. maximum or minimum points) related with the derivative curves. Additionally, bootstrap methods have been used for estimating confidence intervals and speed computational techniques (binning techniques) have been implemented.

-   `r pkg("VDAP")`: The package analyzes Peptide Array Data and characterizes peptide sequence space. It aAllows for high level visualization of global signal, Quality control based on replicate correlation and/or relative Kd, calculation of peptide Length/Charge/Kd parameters, Hits selection based on RFU Signal, and amino acid composition/basic motif recognition with RFU signal weighting. Basic signal trends can be used to generate peptides that follow the observed compositional trends.

-   `r pkg("ZIBseq")`: The package detects abundance differences across clinical conditions. Besides, it takes the sparse nature of metagenomic data into account and handles compositional data efficiently.


### Other packages

-   `r pkg("alc")`: A set of tests for compositional pathologies. Tests for coherence of correlations, compositional dominance of distance with, compositional perturbation invariance and singularity of the covariation matrix.

-   `r pkg("ccmm")`: The package estimates the direct and indirect (mediation) effects of treatment on the outcome when intermediate variables (mediators) are compositional and high-dimensional.

-   `r pkg("CIM")`: The package produces statistical indicators of the impact of migration on the socio-demographic composition of an area. Three measures can be used: ratios, percentages and the Duncan index of dissimilarity. The input data files are assumed to be in an origin-destination matrix format, with each cell representing a flow count between an origin and a destination area. Columns are expected to represent origins, and rows are expected to represent destinations. The first row and column are assumed to contain labels for each area. 

-   `r pkg("CMMs")`: A compositional mediation model for continuous outcome and binary outcomes to deal with mediators that are compositional data.

-   `r pkg("countprop")`: The package calculates metrics of proportionality using the logit-normal multinomial model. It can also provide empirical and plugin estimates of these metrics. 

-   `r pkg("FlexDir")`: The package provides tools to work with the Flexible Dirichlet distribution. The main features are an E-M algorithm for computing the maximum likelihood estimate of the parameter vector and a function based on conditional bootstrap to estimate its asymptotic variance-covariance matrix. It contains also functions to plot graphs, to generate random observations and to handle compositional data.

-   `r pkg("gmGeostats")`: The package offers support for geostatistical analysis of multivariate data, in particular data with restrictions, e.g. positive amounts, compositions, distributional data, microstructural data, etc. It includes descriptive analysis and modelling for such data, both from a two-point Gaussian perspective and multipoint perspective.

-   `r pkg("lba")`: Latent budget analysis is a method for the analysis of a two-way contingency table with an exploratory variable and a response variable. It is specially designed for compositional data.

-   `r pkg("lnmCluster")`: An implementation of logistic normal multinomial (LNM) clustering. It is an extension of the LNM mixture model, and is designed for clustering compositional data. The package includes 3 extended models: LNM Factor Analyzer, LNM Bicluster Mixture Model and Penalized LNM Factor Analyzer. 

-   `r pkg("multilevelcoda")`: Implements Bayesian multilevel modelling for compositional data in a multilevel framework. Compute multilevel compositional data and isometric log-ratio at between and within-person levels, fit Bayesian multilevel models for compositional predictors and outcomes, and run post-hoc analyses such as isotemporal substitution models.

-   `r pkg("SARP.compo")`: Provides a set of functions to interpret changes in compositional data based on a network representation of all pairwise ratio comparisons: computation of all pairwise ratios, construction of a p-value matrix of all pairwise tests of these ratios between conditions, conversion of this matrix to a network.

-   `r pkg("zCompositions")`: Principled methods for the imputation of zeros, left-censored and missing data in compositional data sets

-   `r pkg("zetadiv")`: Functions to compute compositional turnover using zeta-diversity, the number of species shared by multiple assemblages. The package includes functions to compute zeta-diversity for a specific number of assemblages and to compute zeta-diversity for a range of numbers of assemblages. It also includes functions to explain how zeta-diversity varies with distance and with differences in environmental variables between assemblages, using generalised linear models, linear models with negative constraints, generalised additive models,shape constrained additive models, and I-splines.


### Links
-   [gR initiative homepage and mailing list](http://www.R-project.org/gR/)
-   [Bioconductor](http://www.Bioconductor.org/)

### Other resources

-   [CRAN Task View: Cluster](https://cran.r-project.org/view=Cluster)
-   [CRAN Task View: MachineLearning](https://cran.r-project.org/view=MachineLearning)
-   [CRAN Task View: Robust](https://cran.r-project.org/view=Robust)
-   [CRAN Task View: SpatioTemporal](https://cran.r-project.org/view=SpatioTemporal)

