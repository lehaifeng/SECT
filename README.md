# The Smooth Euler Characteristic Transform

Radiogenomics, or imaging genomics, focuses on understanding the relationship between clinical imaging and functional genomic variation. In [Crawford et al. (2016)](http://arxiv.org/), we explore the use of a novel statistic, the smooth Euler characteristic transform (SECT), as an automated procedure to extract geometric or topological statistics from tumor images. The argument for using the SECT for radiogenomics is that more precise geometric information is captured in robust ways when using topological summaries of a shape. The SECT is a variation of the persistent homology transform (PHT), which was initially introduced by [Turner et al. (2014)](https://arxiv.org/abs/1206.2790) to represent shapes and measure distances between CT scans of primate bones as an application in geometric morphometrics. The advantage of the SECT over the PHT is that the statistical summary computed by the SECT is a collection of smooth vectors, while the summary generated by the PHT is a collection of persistent homology diagrams ([Edelsbrunner et al., 2000](https://users.cs.duke.edu/~edels/Papers/2002-J-04-TopologicalPersistence.pdf)) and thus has a complicated representation and geometry ([Turner et al., 2014](https://arxiv.org/abs/1206.2790)). Therefore, unlike the SECT, the PHT is difficult to incorporate into standard statistical models such a linear mixed models. In this work, we examine in detail how SECT features applied to tumor images from The Cancer Genome Atlas (TCGA) compare to gene expression data, as well as classical volumetric and morphometric features, in predicting two clinical outcomes: disease free survival (DFS) and overall survival (OS). We show that the SECT features outperform gene expression, volumetric features, and morphometric features in predicting DFS.

The SECT for 3D images is implemented as a set of MATLAB routines, which can be carried out within both MATLAB and R environments. The Bayesian linear mixed models (LMMs) we used to incorporate shape statistics in predictive analyses is carried out within the R Environment. 

### The MATLAB Environment
MATLAB is a multi-paradigm numerical computing environment and a proprietary programming language developed by [MathWorks](https://www.mathworks.com/index-c.html). MATLAB allows matrix manipulations, plotting of functions and data, implementation of algorithms, creation of user interfaces, and interfacing with programs written in other languages. For more on licensing options, please visit [here](https://www.mathworks.com/campaigns/products/ppc/google/matlab-toolbox-price-request.html?form_seq=reg).

### The R Environment
R is a widely used, free, and open source software environment for statistical computing and graphics. The most recent version of R can be downloaded from the 
[Comprehensive R Archive Network (CRAN)](http://cran.r-project.org/)
CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions that are likely sufficient for many users' needs.  Users can also install R from source code;  however, this may require a significant amount of effort.  For specific details on how to compile, install, and manage R and R-packages, refer to the manual [R Installation and Administration](http://cran.r-project.org/doc/manuals/r-release/R-admin.html).

### R Packages Required for SECT and the Bayesian LMMs
The statistical implementation of the SECT topological summaries using Bayesian LMMs requires the installation of the following R libraries:

[doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)

[R.matlab](https://cran.r-project.org/web/packages/R.matlab/index.html)

[BGLR](https://cran.r-project.org/web/packages/BGLR/index.html)

[diggitdata](https://www.bioconductor.org/packages/devel/data/experiment/html/diggitdata.html)

The easiest method to install these packages is with the following example command entered in an R shell:

    install.packages("BGLR", dependecies = TRUE)

Alternatively, one can also [install R packages from the command line]
                             (http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

### Segmented TCIA Magnetic Resonance Images (MRIs)
MRIs of primary GBM tumors were collected from 92 patients archived by the [The Cancer Imaging Archive (TCIA)] (https://wiki.cancerimagingarchive.net/display/Public/TCGA-GBM), which is a publicly accessible data repository of medical images of cancer patients with matched data in The Cancer Genome Atlas (TCGA) — a collection of a variety of genomic and clinical data for 33 types of cancer. The 92 patients were selected based on two sets of criteria, namely, that they had post-contrast T1 axial MRIs taken at the time of their diagnosis, and that they had available matching (mRNA) gene expression data and clinical correlates on [cBioPortal](http://www.cbioportal.org).

We segmented the TCIA MRI images using a computer-assisted segmentation program to extract tumor lesions from the surrounding brain tissue, which first converts MRI images to a grayscale, and then thresholds to generate binary images. Morphological segmentation is then applied to delineate connected components. More specifically, the program selects contours corresponding to enhanced tumor lesions, which are lighter than healthy brain tissue. As previously noted, necrosis presents as dark regions nested within the indicated lesion. Examples of the raw image obtained from TCIA, and the final segmented result, is given in [Crawford et al. (2016)](http://arxiv.org/) under Figure 7(a) and Figure 7(b), respectively. All segmented TCIA images used in our study can be found in a zipped file in the Data folder.

### Tutorial for Running SECT
The tutorial for computing the SECT topological summaries of the 92 segmented TCIA MRIs is provided here in the MATLAB code folder. Note that the current version of the code only takes .png image files. 

### Tutorial for Running Bayesian LMMs
The tutorial for running a predictive analysis, similar to the one presented in [Crawford et al. (2016)](http://arxiv.org/), can be found in the R Code folder. This code looks at a subset of the 92 TCGA patients, and corresponds to those included in the "diggitdata" Bioconductor R package. This script serves as a simple proof of concept. In order to obtain the mRNA gene expression measurements for all 92 patients, please visit [cBioPortal](http://www.cbioportal.org) or [The Genomic Data Commons Data Portal](https://gdc-portal.nci.nih.gov).

### Questions and Feedback
For questions or concerns with the SECT functions, please contact
[Lorin Crawford](mailto:lac55@stat.duke.edu) or [Anthea Monod](mailto:rr2579@cumc.columbia.edu).

We appreciate any feedback you may have with our repository and instructions.
