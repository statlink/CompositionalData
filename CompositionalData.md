---
name: CompositionalData
topic: Compositional Data Analysis
maintainer: Michail Tsagris and Patrice Kiener
email: mtsagris@uoc.gr
version: 2023-09-26
source: https://github.com/cran-task-views/CompositionalData/
---

*Compositional data are positive multivariate data where the sum of the values of each vector sums to the same constant, typically taken to be unity for convenience purposes. In this case, their sample space is the standard simplex. The most popular approach is to use the logarithm transformation applied to ratios of the variables, initially suggested by [Aitchison  (1982)](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x). However this approach has drawbacks and for this many alternative transformations have been developed throughout the years. Additionally, the literature numbers many different methodologies and approaches, some of which are publicly available in the packages listed below.*

This task view is a collection of packages intended to supply R code to deal with 
Compositional data analysis. The packages can be roughly structured into the following topics
(although several of them have functionalities which go across these
categories):


### General purpose packages

-   `r pkg("Compositional", priority = "core")`: The package performs regression, classification, contour plots, hypothesis testing and fitting of distributions with or without transformations. Further, machine learning techniques are included, see for instance the CTV view("MachineLearning"). The packages further includes functions for percentages (or proportions). Notable references include [Tsagris, Preston and Wood (2011)](https://doi.org/10.48550/arXiv.1106.1451), [Tsagris and Stewart (2018)](https://doi.org/10.1134/S1995080218030198), [Alenazi (2023)](https://doi.org/10.1080/03610926.2021.2014890) and [Tsagris, Alenazi and Stewart (2023)](https://doi.org/10.1007/s11222-023-10277-5).

-   `r pkg("compositions", priority = "core")`: Provides functions for the consistent analysis of compositional data (e.g. portions of substances) and positive numbers (e.g. concentrations). References include the book by [Aitchison (1986)](https://www.amazon.com/Statistical-Analysis-Compositional-Data/dp/1930665784) and [van den Boogaart. and Tolosana-Delgado (2008)](https://doi.org/10.1016/j.cageo.2006.11.017).

-   `r pkg("easyCODA", priority = "core")`: Univariate and multivariate methods for compositional data analysis, based on log-ratios. Accent is given to simple pairwise log-ratios. Selection can be made of log-ratios that account for a maximum percentage of log-ratio variance. Various multivariate analyses of log-ratios are included in the package. The package implements the approach in the book by [Greenacre (2018)](https://doi.org/10.1201/9780429455537).

-   `r pkg("zCompositions", priority = "core")`: Principled methods for the imputation of zeros, left-censored and missing data in compositional data sets as descibed in [Palarea-Albaladejo and Martin-Fernandez (2015)](hhttps://doi.org/10.1016/j.chemolab.2015.02.019).


### Robust methods 

-   `r pkg("complmrob")`: Robust regression methods for compositional data. The distribution of the estimates can be approximated with various bootstrap methods. These bootstrap methods are available for the compositional as well as for standard robust regression estimates. This allows for direct comparison between them. See [Hron, Filzmoser and Thompson (2012)](https://doi.org/10.1080/02664763.2011.644268). 

-   `r pkg("robCompositions", priority = "core")`: Methods for analysis of compositional data including robust methods, imputation of missing values, methods to replace rounded zeros, count zeros, methods to deal with essential zeros, (robust) outlier detection for compositional data, (robust) principal component analysis for compositional data, (robust) factor analysis for compositional data, (robust) discriminant analysis for compositional data (Fisher rule), robust regression with compositional predictors, functional data analysis and p-splines, contingency and compositional tables and (robust) Anderson-Darling normality tests for compositional data as well as popular log-ratio transformations. In addition, visualisation and diagnostic tools are implemented as well as high and low-level plot functions for the ternary diagram. References include [Filzmoser, Hron and Templ (2018) ](https://doi.org/10.1007/978-3-319-96422-5) and [Hron et al. (2016)](https://doi.org/10.1016/j.csda.2015.07.007).

-   `r pkg("robregcc")`: The package implements an algorithm estimating the parameters of the robust regression model with compositional covariates. The model simultaneously treats outliers and provides reliable parameter estimates. The algorithm is described in [Mishra and Mueller (2019)](https://doi.org/10.48550/arXiv.1909.04990).

- See also the CTV view("Robust") for more robust techniques.  

### Regression modelling

-   `r pkg("codalm")`: Implements the expectation-maximization algorithm for transformation-free linear regression for compositional outcomes and predictors of [Fiksel, Zeger and Datta (2021)](https://doi.org/10.1111/biom.13465).

-   `r pkg("codaredistlm")`: Provided data containing an outcome variable, compositional variables and additional covariates (optional); linearly regress the outcome variable on an isometric log ratio transformation of the linearly dependent compositional variables. The package provides predictions (with confidence intervals) in the change in the outcome/response variable based on the multiple linear regression model and evenly spaced reallocations of the compositional values. References include [Dumuid et al. (2017a)](https://doi.org/10.1177/0962280217710835) and [Dumuid et al. (2017b)](https://doi.org/10.1177/0962280217737805).

-   `r pkg("Compack")`: Regression methodologies with compositional covariates, including (1) sparse log-contrast regression with compositional covariates by [Lin et al. (2014)](https://doi.org/10.1093/biomet/asu031), and (2) sparse log-contrast regression with functional compositional predictors by [Sun et al. (2020)]( https://doi.org/10.1214/20-AOAS1357).

-   `r pkg("DirichletReg")`: The package implements Dirichlet regression models. 

-   `r pkg("FLORAL")`: Log-ratio lasso regression for continuous, binary, and survival outcomes with compositional features by [Fei et al. (2023)](https://doi.org/10.1101/2023.05.02.538599).

-   `r pkg("ocomposition")`: Regression model where the response variable is a rank-indexed compositional vector (non-negative values that sum up to one and are ordered from the largest to the smallest). Parameters are estimated in the Bayesian framework using MCMC methods. The relevant paper is [Rozenas and  Alvarez (2012)](https://doi.org/10.1093/pan/mpr041). See also the CTV view("Bayesian") for more packages dedicated to Bayesian statistics.  

-   `r pkg("zoid")`: Fits Dirichlet regression and zero-and-one inflated Dirichlet regression with Bayesian methods implemented in Stan. These models are sometimes referred to as trinomial mixture models; covariates and overdispersion can optionally be included. See also the CTV view("Bayesian") for more packages dedicated to Bayesian statistics.  



### Other packages

-   `r pkg("aIc")`: Tests for coherence of correlations as suggested by [Erb et al. (2020)](https://doi.org/10.1016/j.acags.2020.100026), compositional dominance of distance, compositional perturbation invariance as suggested by [Aitchison (1992)](https://doi.org/10.1007/BF00891269) and singularity of the covariation matrix. 

-   `r pkg("ccmm")`: The package estimates the direct and indirect (mediation) effects of treatment on the outcome when intermediate variables (mediators) are compositional and high-dimensional. For more information see [Sohn and Li (2019)]( https://doi.org/10.1214/18-AOAS1210).

-   `r pkg("CIM")`: The package produces statistical indicators of the impact of migration on the socio-demographic composition of an area. Three measures can be used: ratios, percentages and the Duncan index of dissimilarity. The input data files are assumed to be in an origin-destination matrix format, with each cell representing a flow count between an origin and a destination area. Columns are expected to represent origins, and rows are expected to represent destinations. The first row and column are assumed to contain labels for each area. See [Rodriguez-Vignoli and Rowe (2018)](https://doi.org/10.1080/00324728.2017.1416155) for technical details. 

-   `r pkg("CMMs")`: A compositional mediation model for continuous outcome and binary outcomes to deal with mediators that are compositional data. For more information see [Lin et al. (2022)](https://doi.org/10.1016/j.jad.2021.12.019).

-   `r pkg("coda.base")`: A minimum set of functions to perform compositional data analysis using the log-ratio approach introduced by [Aitchison (1982)](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x). Main functions have been implemented in C++ for better performance.

-   `r pkg("countprop")`: The package calculates metrics of proportionality using the logit-normal multinomial model. It can also provide empirical and plugin estimates of these metrics. 

-   `r pkg("FlexDir")`: The package provides tools to work with the Flexible Dirichlet distribution by [Migliorati, Ongaro and Monti (2017)](https://doi.org/10.1007/s11222-016-9665-y). The main features are an E-M algorithm for computing the maximum likelihood estimate of the parameter vector and a function based on conditional bootstrap to estimate its asymptotic variance-covariance matrix. It contains also functions to plot graphs, to generate random observations and to handle compositional data.

-   `r pkg("gmGeostats")`: The package offers support for geostatistical analysis of multivariate data, in particular data with restrictions, e.g. positive amounts, compositions, distributional data, microstructural data, etc. It includes descriptive analysis and modelling for such data, both from a two-point Gaussian perspective and multipoint perspective.  The methods mainly follow [Tolosana-Delgado, Mueller and van den Boogaart (2018)](https://doi.org/10.1007/s11004-018-9769-3). See also the CTV view("SpatioTemporal").

-   `r pkg("lba")`: Latent budget analysis is a method for the analysis of a two-way contingency table with an exploratory variable and a response variable. It is specially designed for compositional data. The package is based on the PhD thesis of an der Ark (1999) at the University of Utrecht titled "Contributions to Latent Budget Analysis, a tool for the analysis of compositional data".

-   `r pkg("multilevelcoda")`: Implements Bayesian multilevel modelling for compositional data in a multilevel framework. Compute multilevel compositional data and isometric log-ratio at between and within-person levels, fit Bayesian multilevel models for compositional predictors and outcomes, and run post-hoc analyses such as isotemporal substitution models. See also the CTV view("Bayesian") for more packages dedicated to Bayesian statistics.  

-   `r pkg("NBDdirichlet")`: Implements NBD-Dirichlet Model of Consumer Buying Behavior for Marketing
Research. The relevant paper is [Goodhardt, Ehrenberg and Chatfield (1984)](https://doi.org/10.2307/2981696).

-   `r pkg("rBeta2009")`: The packate contains a function for generating Dirichlet random vectors using a C translation. The relevant paper is [Hung, Balakrishnan and Cheng (2011)](https://doi.org/10.1080/00949650903409999).

-   `r pkg("SARP.compo")`: Provides a set of functions to interpret changes in compositional data based on a network representation of all pairwise ratio comparisons: computation of all pairwise ratios, construction of a p-value matrix of all pairwise tests of these ratios between conditions, conversion of this matrix to a network.

-   `r pkg("ToolsForCoDa")`: Provides functions for multivariate analysis with compositional data. Includes a function for doing compositional canonical correlation analysis. This analysis requires two data matrices of compositions, which can be adequately transformed and used as entries in a specialized program for canonical correlation analysis, that is able to deal with singular covariance matrices. A function for log-ratio principal component analysis with condition number computations has been added to the package. The methodology is described in [Graffelman et al. (2017)](https://doi.org/10.1101/144584).

-   `r pkg("zetadiv")`: Functions to compute compositional turnover using zeta-diversity by [Hui and  McGeoch (2014)](https://doi.org/10.1086/678125), the number of species shared by multiple assemblages. The package includes functions to compute zeta-diversity for a specific number of assemblages and to compute zeta-diversity for a range of numbers of assemblages. It also includes functions to explain how zeta-diversity varies with distance and with differences in environmental variables between assemblages, using generalised linear models, linear models with negative constraints, generalised additive models,shape constrained additive models, and I-splines.


### Ternary diagrams

-   `r pkg("Ternary", priority = "core")`: Plots ternary diagrams (simplex plots/Gibbs triangles) and Holdridge life zone plots (see [Holdridge (1947)](https://www.science.org/doi/10.1126/science.105.2727.367)) using the standard graphics functions. An alternative to 'ggtern', which uses the 'ggplot2' family of plotting functions. Includes a 'Shiny' user interface for point-and-click ternary plotting. 

-   `r pkg("ggtern", priority = "core")`: The package extends the functionality of 'ggplot2', providing the capability to plot ternary diagrams for (subset of) the 'ggplot2' geometries. Additionally, 'ggtern' has implemented several new geometries which are unavailable to the standard 'ggplot2' release. For further examples and documentation, please proceed to the 'ggtern' website.

-   `r pkg("isopleuros", priority = "core")`: Ternary plots made simple. This package allows to create ternary plots using 'graphics'. It provides functions to display the data in the ternary space, to add or tune graphical elements and to display statistical summaries. It also includes common ternary diagrams which are useful for the archaeologist (e.g. soil texture charts, ceramic phase diagram).

-   `r pkg("provenance")`: The package includes tools to plot compositional and count data on ternary diagrams and point-counting data on radial plots, to calculate the sample size required for specified levels of statistical precision, and to assess the effects of hydraulic sorting on detrital compositions. Includes an intuitive query-based user interface for users who are not proficient in R.



### Bioinformatics/ecology related packages

-   `r pkg("ampir")`: A toolkit to predict antimicrobial peptides from protein sequences on a genome-wide scale. It incorporates two support vector machine models ("precursor" and "mature") trained on publicly available antimicrobial peptide data using calculated physico-chemical and compositional sequence properties. In order to support genome-wide analyses, these models are designed to accept any type of protein as input and calculation of compositional properties has been optimised for high-throughput use. For more information see [Fingerhut et al. (2020)](hhttps://doi.org/10.1093/bioinformatics/btaa653).

-   `r pkg("ArArRedux")`: Processes noble gas mass spectrometer data to determine the isotopic composition of argon (comprised of Ar36, Ar37, Ar38, Ar39 and Ar40) released from neutron-irradiated potassium-bearing minerals. Then uses these compositions to calculate precise and accurate geochronological ages for multiple samples as well as the covariances between them. Error propagation is done in matrix form, which jointly treats all samples and all isotopes simultaneously at every step of the data reduction process. Includes methods for regression of the time-resolved mass spectrometer signals to t=0 ('time zero') for both single- and multi-collector instruments, blank correction, mass fractionation correction, detector intercalibration, decay corrections, interference corrections, interpolation of the irradiation parameter between neutron fluence monitors, and (weighted mean) age calculation. All operations are performed on the logs of the ratios between the different argon isotopes so as to properly treat them as 'compositional data'.

-   `bioc("banocc")`: This is a Bioconductor package designed for compositional data, where each sample sums to one. It infers the approximate covariance of the unconstrained data using a Bayesian model coded with `r pkg("ArArRedux")`. It provides as output the "stanfit" object as well as posterior median and credible interval estimates for each correlation element. See also the CTV view("Bayesian") for more packages dedicated to Bayesian statistics.  

-   `r pkg("BRACoD.R")`: The goal of this method is to identify associations between bacteria and an environmental variable in 16S or other compositional data. The environmental variable is any variable which is measure for each microbiome sample, for example, a butyrate measurement paired with every sample in the data. Microbiome data is compositional, meaning that the total abundance of each sample sums to 1, and this introduces severe statistical distortions. This method takes a Bayesian approach to correcting for these statistical distortions, in which the total abundance is treated as an unknown variable. This package runs the python implementation using reticulate.
See also the CTV view("Bayesian") for more packages dedicated to Bayesian statistics.  

-   `r pkg("coda4microbiome")`: Functions for microbiome data analysis that take into account its compositional nature. Performs variable selection through penalized regression for both, cross-sectional and longitudinal studies, and for binary and continuous outcomes. For more information see [Bokulich et al. (2016)](https://doi.org/10.1126/scitranslmed.aad7121).

-   `r pkg("codacore")`: In the context of high-throughput genetic data, the package identifies a set of sparse biomarkers that are predictive of a response variable of interest [Gordon-Rodriguez, Quinn and Cunningham (2022)](https://academic.oup.com/bioinformatics/article/38/1/157/6366546?login=false). More generally, the package can be applied to any regression problem where the independent variable is compositional, to derive a set of scale-invariant log-ratios that are maximally associated to a dependent variable.

-   `bioc("banocc")`: This is a Bioconductor package and contains the explorative ordination method that combines quasi-likelihood estimation, compositional regression models and latent variable models for integrative visualization of several omics datasets. Both unconstrained and constrained integration are available. The results are shown as interpretable, compositional multiplots.

-   `r pkg("lnmCluster")`: An implementation of logistic normal multinomial (LNM) clustering. It is an extension of the LNM mixture model proposed by [Fang and Subedi (2020)](https://doi.org/10.1038/s41598-023-41318-8), and is designed for clustering compositional data. The package includes 3 extended models: LNM Factor Analyzer, LNM Bicluster Mixture Model and Penalized LNM Factor Analyzer. Details for model assumptions and interpretation can be found in the papers: [Tu and Subedi (2021)](https://doi.org/10.48550/arXiv.2101.01871) and [Tu and Subedi (2022)](https://doi.org/10.1002/sam.11555). See also the CTV view("Cluster").

-   `r pkg("MicrobiomeStat")`: A suite of methods for powerful and robust microbiome data analysis addressing zero-inflation, phylogenetic structure and compositional effects [Zhou et al. (2022)](https://doi.org/10.1186/s13059-022-02655-5). The methods can be applied to the analysis of other (high-dimensional) compositional data arising from sequencing experiments.

-   `r pkg("QFASA")`: Accurate estimates of the diets of predators are required in many areas of ecology, but for many species current methods are imprecise, limited to the last meal, and often biased. The diversity of fatty acids and their patterns in organisms, coupled with the narrow limitations on their biosynthesis, properties of digestion in monogastric animals, and the prevalence of large storage reservoirs of lipid in many predators, led to the development of quantitative fatty acid signature analysis (QFASA) to study predator diets. Some relevant papers are [Aitchison (1992)](https://doi.org/10.1007/BF00891269), [Stewart, Iverson and Field (2014)](https://doi.org/10.1007/s10651-014-0280-9) and [Stewart (2017)](https://doi.org/10.1080/02664763.2016.1193846).

-   `bioc("scDDboost")`: This is a Bioconductor package to analyze changes in the distribution of single-cell expression data between two experimental conditions. Compared to other methods that assess differential expression, scDDboost benefits uniquely from information conveyed by the clustering of cells into cellular subtypes. Through a novel empirical Bayesian formulation it calculates gene-specific posterior probabilities that the marginal expression distribution is the same (or different) between the two conditions. The implementation in scDDboost treats gene-level expression data within each condition as a mixture of negative binomial distributions.
See also the CTV view("Bayesian") for more packages dedicated to Bayesian statistics.  

-   `r pkg("VDAP")`: The package analyzes Peptide Array Data and characterizes peptide sequence space. It aAllows for high level visualization of global signal, Quality control based on replicate correlation and/or relative Kd, calculation of peptide Length/Charge/Kd parameters, Hits selection based on RFU Signal, and amino acid composition/basic motif recognition with RFU signal weighting. Basic signal trends can be used to generate peptides that follow the observed compositional trends.

-   `r pkg("ZIBseq")`: The package detects abundance differences across clinical conditions [Xiaoling, Gang, and Zhenqiu (2016)](https://doi.org/10.1089/cmb.2015.0157). Besides, it takes the sparse nature of metagenomic data into account and handles compositional data efficiently.


### Links

-   [CRAN Task View: Bayesian](https://cran.r-project.org/view=Bayesian)
-   [CRAN Task View: Cluster](https://cran.r-project.org/view=Cluster)
-   [CRAN Task View: MachineLearning](https://cran.r-project.org/view=MachineLearning)
-   [CRAN Task View: Robust](https://cran.r-project.org/view=Robust)
-   [CRAN Task View: SpatioTemporal](https://cran.r-project.org/view=SpatioTemporal)

