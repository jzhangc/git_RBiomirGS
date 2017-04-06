# RBioMiRGS
Simple to use package for both microarray and RNAseq data analysis (limma based)

Installation:

  - Install devtools
  
        install.packages("devtools")
    
  - Install bioconductor
  
        source("https://bioconductor.org/biocLite.R")
      
        biocLite()
    
  - Install the package
  
        devtools::install_github("jzhangc/git_RBiomiRGS/RBiomiRGS", repos = BiocInstaller::biocinstallRepos())   

Update log

    0.1.4 - 0.1.5 (April.6.2017)
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
