# RBiomirGS
Simple to use package for miRNA gene set analysis, with miRNA-mRNA interaction searching and human mRNA orthologs conversion functionalities.


To cite in publication

    Zhang J, Storey KB. 2018. RBiomirGS: An all-in-one miRNA gene set analysis solution featuring target mRNA mapping and expression profile integration. PeerJ. 6: e4262.
  


Installation:

  - Install devtools
  
        install.packages("devtools")
    
  - Install bioconductor
  
        if (!requireNamespace("BiocManager"))
            install.packages("BiocManager")
            
        BiocManager::install()
    
  - Install stable release
  
        devtools::install_github("jzhangc/git_RBiomiRGS/RBiomirGS", repos = BiocManager::repositories())
  
  - Install development build
  
        devtools::install_github("jzhangc/git_RBiomiRGS/RBiomirGS", repos = BiocManager::repositories(), ref = "beta")

Update log
    
    0.2.20 (January.3.2024)
      - The hsa_entrez list now is also a "mir_entrez_list" class
      - rbiomirgs_logistic and rbiomirgs_logisticV2 updated with argument check
    
    0.2.19 (September.6.2023)
      - A bug fixed for rbiomirgs_volcano() where the function will crash when no FDR significance was found while FDR correction is active

    0.2.18 (March.25.2023)
      - A bug fixed for rbiomirgs_logistic() where the "mir_entrez_list" check fails to give proper error messages

    0.2.17 (November.16.2022)
      - rbiomirgs_logistic() updated with error message ("incompabile") if to use a "mir_entrez_list" class

    0.2.16 (March.21.2022)
      - It is now possible to custom file name for rbiomirgs_volcano() and rbiomirgs_bar()
      - A bug fixed for rbiomirgs_volcano() where label for target passing alpha wasn't properly displayed when fdr=TRUE
    
    0.2.15 (December.20.2021)
      - The entrez list derived from mirnascan() now has both "list" and "mir_entrez_list" classes
      - A new rbiomirgs_logisticV2() function added that always checks if the "mrnalist" is a "mir_entrez_list" class
        - The original rbiomirgs_logistic() function still available for compatibility purposes
        - The new "mir_entrez_list" classes works with both rbiomirgs_logistic and rbiomirgs_logisticV2() functions
        - The old entrez list (i.e. only a "list" class) will NOT work with the rbiomirgs_logisticV2() function
        
    0.2.14 (April.16.2021)
      - p_line_offset added to the rbiomirgs_volcano() function. See manual for details
    
    0.2.13 (Mar.25.2021)
      - rbiomirgs_logistic() updated to handle empty miRNA ID entries
        - Small updates to rbiomirgs_logistic() message delivery
    
    0.2.12 (Aug.12.2019)
      - The default server url updated to "http://multimir.org/cgi-bin/multimir_univ.pl"
      - rbiomirgs_mrnascan() updated with updated biomaRt entrezgene identifier
      - Small fixes
    
    0.2.10
      - Additional argument check added to rbiomirgs_mrnascan()
      - New bioconductor installation instructions added

    0.2.9
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
