# RBiomirGS
Simple to use package for miRNA gene set analysis, with miRNA-mRNA interaction searching and human mRNA orthologs conversion functionalites.

Installation:

  - Install devtools
  
        install.packages("devtools")
    
  - Install bioconductor
  
        source("https://bioconductor.org/biocLite.R")
      
        biocLite()
    
  - Install the package
  
        devtools::install_github("jzhangc/git_RBiomiRGS/RBiomiRGS", repos = BiocInstaller::biocinstallRepos())   

Update log

    0.2.4 (May.17.2017)
      - miRNA and mRNA score output file name adjusted
      - Bug fixes

    0.2.3
      - miRNA:mRNA interaction weight functionality fixed
      - Parameter optimization for rbiomirgs_logistic() adjusted
      - Other bug fixes
    
    
    0.2.2 
      - arugment adjusted for the volcano plot function
    
    0.2.1 
      - mRNA score output added

    0.2.0
      - Package name changed to RBiomirGS
      - Plot functions added: rbiomirgs_volcano(), rbiomirgs_histogram()
      - Adjustments on the variable names for logitstic regression resutls
      - mRNA scoring modified with reversed sign, so that positive value means actiation and negative means inhibition. Same goes with the parameters obtained from logistic regression
      - Function names are now all lower cases
      - Bug fixes
    
    0.1.8
      - Preparation for the plotting module
      - Description file updated
      - zzz.R file added
      - Bug fixes

    0.1.6 - 0.1.7
      - Bug fixes

    0.1.4 - 0.1.5 
      - logistic regression-based miRNA gene set analysis function rbiomirGS_logistic() added, featuring multiple parameter optimization algorithms as well as parallel computing
      - many "quality of life" features added for both rbiomirGS_logistic() and rbiomirGS_mrnalist() function
      - codes optmized to reduce redundancy 
      - Bug fixes
    
    0.1.3
      - rbiomirGS_mrnalist() now is able to connect to ensembl database to convert mmu/rno entrez gene IDs to the hsa ortholog entrez IDs
      - rbiomirGS_gmt() function added to load gmt gene set files from ensembl databases
      - Bug fixes

    0.1.2
      - In addition to the environment, rbiomirGS_mrnalist() now outputs csv files of the mRNA resutls, as well as a hsa entrez list to the environment for modelling use

    0.1.1
      - miRNA to mRNA function rbiomirGS_mrnalist() added.

    0.1.0 
      - Initial release
