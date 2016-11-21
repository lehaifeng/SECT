# The Smooth Euler Characteristic Transform

Radiogenomics, or imaging genomics, focuses on understanding the relationship between clinical imaging and functional genomic variation. In [Crawford et al. (2016)](http://arxiv.org/), we explore the use of a novel statistic, the smooth Euler characteristic transform (SECT), as an automated procedure to extract geometric or topological statistics from tumor images. The argument for using the SECT for radiogenomics is that more precise geometric information is captured in robust ways when using topological summaries of a shape. The SECT is a variation of the persistent homology transform (PHT), which was initially introduced by [Turner et al. (2014)](https://arxiv.org/abs/1206.2790) to represent shapes and measure distances between CT scans of primate bones as an application in geometric morphometrics. The advantage of the SECT over the PHT is that the statistical summary computed by the SECT is a collection of smooth vectors, while the summary generated by the PHT is a collection of persistent homology diagrams ([Edelsbrunner et al., 2000](https://users.cs.duke.edu/~edels/Papers/2002-J-04-TopologicalPersistence.pdf)) and thus has a complicated representation and geometry ([Turner et al., 2014](https://arxiv.org/abs/1206.2790)). Unlike the SECT, therefore, the PHT is difficult to incorporate into standard statistical models such a linear mixed models. In this work, we examine in detail how SECT features applied to tumor images from The Cancer Genome Atlas (TCGA) compare to gene expression data, as well as classical volumetric and morphometric features, in predicting two clinical outcomes: disease free survival (DFS) and overall survival (OS). We show that the SECT features outperform gene expression, volumetric features, and morphometric features in predicting DFS.

The SECT for 3D images is implemented as a set of MATLAB routines, which can be carried out within both MATLAB and R environments. The Bayesian linear mixed models (LMMs) we used to incorporate shape statistics in predictive analyses is carried out within the R Environment. 

### The MATLAB Environment



### The R Environment
R is a widely used, free, and open source software environment for statistical computing and graphics. The most recent version of R can be downloaded from the 
[Comprehensive R Archive Network (CRAN)](http://cran.r-project.org/)
CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions that are likely sufficient for many users' needs.  Users can also install R from source code;  however, this may require a significant amount of effort.  For specific details on how to compile, install, and manage R and R-packages, refer to the manual [R Installation and Administration](http://cran.r-project.org/doc/manuals/r-release/R-admin.html).


### R Packages Required for SECT and the Bayesian LMMs
The SECT tutorial requires the installation of the following R libraries:

[doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)

[R.matlab](https://cran.r-project.org/web/packages/R.matlab/index.html)

[BGLR](https://cran.r-project.org/web/packages/BGLR/index.html)

The easiest method to install these packages is with the following example command entered in an R shell:

install.packages("BGLR", dependecies = TRUE)

Alternatively, one can also [install R packages from the command line]
                             (http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

### Segmented TCIA Magnetic Resonance Images


### Tutorial for Running SECT
The tutorial provided here is based on 

### Questions and Feedback
For questions or concerns with the SECT functions, please contact
[Lorin Crawford](mailto:lac55@stat.duke.edu) or [Anthea Monod](mailto:rr2579@cumc.columbia.edu).

We appreciate any feedback you may have with our repository and instructions.
