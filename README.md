# RBiomirGS
Simple to use package for miRNA gene set analysis, with miRNA-mRNA interaction searching and human mRNA orthologs conversion functionalities.


To cite in publication

    Zhang J, Storey KB. 2018. RBiomirGS: An all-in-one miRNA gene set analysis solution featuring target mRNA mapping and expression profile integration. PeerJ. 6: e4262.
  


Installation:

  - Install devtools
  
        install.packages("devtools")
    
  - Install bioconductor
  
        source("https://bioconductor.org/biocLite.R")
      
        biocLite()
    
  - Install stable release
  
        devtools::install_github("jzhangc/git_RBiomiRGS/RBiomirGS", repos = BiocInstaller::biocinstallRepos())
  
  - Install development build
  
        devtools::install_github("jzhangc/git_RBiomiRGS/RBiomirGS", repos = BiocInstaller::biocinstallRepos(), ref = "beta")

Update log

    0.2.9 (Feb.22.2018)
      - ratioFC arugment added to rbiomirgs_logistic()
      - The default value for var_mirnaFC of rbiomirgs_logistic() changed to "logFC"

    0.2.7 - 0.2.8
      - Citation information added
      - Bug fixes
    
    0.2.6
      - rbiomirgs_histogram() changed to rbiomirgs_bar()
      - Bug fixes
    
    0.2.5
      - p value thresholding added for rbiomirgs_histogram()
      - Bug fixes

    0.2.4
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
