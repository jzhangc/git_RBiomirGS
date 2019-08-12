#' @title rbiomirgs_mrnascan
#'
#' @description Obtain target mRNA list for miRNAs of interest. Results can either be predicted or validated miRNA-mRNA interactions. The function uses multiple databsaes hosted at \code{multimir.ucdenver.edu/}.This function needs a internet connection.
#' @param mir Input miRNAs vector. Can't be a matrix or data frame.
#' @param sp Species. Options are \code{"hsa"} (default), \code{"mmu"} and \code{"rno"}.
#' @param addhsaEntrez When \code{sp = "mmu"} or \code{sp = "rno"}, users can set this argument to \code{TRUE} so that a new list containing hsa ortholog entrez ID will be exported to the environment. The function connects to the up-to-date \code{ensembl} databases.
#' @param queryType Type of reuslts. Options are \code{"validated"} and \code{"predicted"}.
#' @param predictPrecentage Set only if \code{queryType = "predicted"}. The percentage of the top scored mRNA predictions to return. Default is \code{5}.
#' @param url The database host server: \code{"http://multimir.org/cgi-bin/multimir_univ.pl"}
#' @param parallelComputing If to use parallel computing or not. Default is \code{FALSE}.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @return The function returns a \code{list} object with mRNA targets for the query miRNAs, as well as detailed resutls to the working directory as \code{csv} files. When \code{sp = "mmu"} or \code{sp = "rno"}, and \code{addhsaEntrez = TRUE}, a new list containing hsa ortholog entrez ID will be exported to the environment.
#' @import doParallel
#' @import foreach
#' @import biomaRt
#' @importFrom XML readHTMLTable
#' @importFrom RCurl postForm
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @examples
#' \dontrun{
#' rbiomirgs_mrnascan(objTitle = "mmu_mirna", mir = c("mmu-miR-26a-5p", "mmu-miR-100-5p"), sp = "mmu", addhsaEntrez = TRUE, queryType = "predicted", predictPercentage = 10)
#' }
#' @export
rbiomirgs_mrnascan <- function(objTitle = "miRNA", mir =  NULL, sp = "hsa", addhsaEntrez = FALSE,
                               queryType = NULL, predictPercentage = 10,
                               url = "http://multimir.org/cgi-bin/multimir_univ.pl",
                               parallelComputing = FALSE, clusterType = "PSOCK"){

  #### check the arguments
  if (is.null(mir)) stop("Please set the input miRNA(s). Either single targets or a vector of mutliple ones can be used.")
  if (!is.null(dim(mir))) stop("The input argument mir needs to be a vector. ")

  if (!sp %in% c("hsa", "mmu", "rno")){ # check species
    stop(cat("Only human, mouse or rat are supported for now. Please choose either \"hsa\" (default), \"mmu\", or \"rno\" for species."))
  }
  if (sp == "hsa" && addhsaEntrez) {
    cat("Argument addhsaEntrez automatically set to FALSE when sp = \"hsa\".\n")
    addhsaEntrez = FALSE
  }

  #### set up the search input
  if (is.null(queryType)){
    stop(cat("Please set the queryType argument. Options are \"validated\" and \"predicted\"."))
  } else if (queryType != "validated" & queryType != "predicted"){
    stop(cat("queryType only takes \"validated\" or \"predicted\". Check the spell."))
  } else if (queryType == "validated"){
    db <- c("mirecords", "mirtarbase", "tarbase")
  } else if (queryType == "predicted"){
    if (sp == "rno"){
      db <- c("elmmo", "microcosm", "miranda", "mirdb")  # no "diana_microt", "pictar", "pita", "targetscan" for rno
    } else {
      db <- c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb",
              "pictar", "pita", "targetscan")
    }
  }

  #### set up a tmpfunc for query
  # i - target miRNA
  # j - single database
  # mode - queryType ("validated", "predicted")
  # percenetage - predictPercentage
  tmpfunc_q <- function(i, j, mode = NULL, percentage = NULL){
    target.table <- "target"
    mirna.table <- "mirna"

    tmpmirna <- i
    tmpmirna <- paste(tmpmirna, collapse = "','")
    tmpmirna <- paste("('", tmpmirna, "')", sep = "")

    ## set up query syntax
    if (mode == "validated"){ # validated databases
      q <- paste("SELECT m.mature_mirna_id,", "t.target_symbol, t.target_entrez, t.target_ensembl,",
                 "i.experiment, i.pubmed_id FROM",
                 mirna.table,
                 "AS m INNER JOIN",
                 j,
                 "AS i INNER JOIN",
                 target.table,
                 "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid AND", "i.target_uid=t.target_uid) WHERE",
                 sep = " ") # set databases for query
      q <- paste(q,
                 "(m.mature_mirna_acc IN",
                 tmpmirna,
                 "OR m.mature_mirna_id IN",
                 tmpmirna,
                 ")",
                 sep = " ") # query for miRNA
      q <- paste(q,
                 " AND m.org = '",
                 sp,
                 "' AND t.org = '",
                 sp,
                 "'",
                 sep = "") # set up species

    } else { # predicted databases
      # from multiMiR pacakge, to get the score from the database
      tmpfunc_cutoff <- function(cutoff.file = "http://multimir.ucdenver.edu/multimir_cutoffs.rda"){
        multimir_cutoffs <- NULL
        url.file <- url(cutoff.file)
        on.exit(close(url.file))
        load(url.file)
        return(multimir_cutoffs)
      }

      # initial syntax
      q <- paste("SELECT m.mature_mirna_acc, m.mature_mirna_id,",
                 "t.target_symbol, t.target_entrez, t.target_ensembl FROM",
                 mirna.table, "AS m INNER JOIN", j, "AS i INNER JOIN",
                 target.table, "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid",
                 "AND i.target_uid=t.target_uid) WHERE",
                 sep = " ")

      if (j == "diana_microt") {
        q <- sub(" FROM ", ", i.miTG_score AS score FROM ", q)
      } else if (j == "elmmo") {
        q <- sub(" FROM ", ", i.p AS score FROM ", q)
      } else if (j %in% c("microcosm", "mirdb", "pictar")) {
        q <- sub(" FROM ", ", i.score FROM ", q)
      } else if (j == "miranda") {
        q <- sub(" FROM ", ", i.mirsvr_score AS score FROM ", q)
      } else if (j == "pita") {
        q <- sub(" FROM ", ", i.ddG AS score FROM ", q)
      } else if (j == "targetscan") {
        q <- sub(" FROM ", ", i.context_plus_score AS score FROM ", q)
      }
      q <- paste(q, "(m.mature_mirna_acc IN", tmpmirna,
                 "OR m.mature_mirna_id IN", tmpmirna, ")",
                 sep = " ")
      q <- paste(q, " AND m.org = '", sp, "' AND t.org = '", sp, "'",
                 sep = "")  #add organism to the query

      # add dataset-specific cutoff to the query
      nm <- paste(j, sp, sep = ".")
      tmp_cutoff <- tmpfunc_cutoff()
      cutoff <- tmp_cutoff[[nm]][[paste(percentage, "%", sep = "")]]

      if (j == "diana_microt") {
        q <- paste(q, "AND i.miTG_score >=", cutoff,
                   "ORDER BY i.miTG_score DESC", sep = " ")
      } else if (j == "elmmo") {
        q <- paste(q, "AND i.p >=", cutoff, "ORDER BY i.p DESC",
                   sep = " ")
      } else if (j %in% c("microcosm", "mirdb", "pictar")) {
        q <- paste(q, "AND i.score >=", cutoff,
                   "ORDER BY i.score DESC", sep = " ")
      } else if (j == "miranda") {
        q <- paste(q, "AND i.mirsvr_score <=", cutoff,
                   "ORDER BY i.mirsvr_score", sep = " ")
      } else if (j == "pita") {
        q <- paste(q, "AND i.ddG <=", cutoff, "ORDER BY i.ddG",
                   sep = " ")
      } else if (j == "targetscan") {
        # q <- paste(q, 'AND i.site_type == 3', sep=' ')
        q <- paste(q, "AND i.context_plus_score <=", cutoff,
                   "ORDER BY i.context_plus_score",
                   sep = " ")
      }
    }

    ## query
    tmprslt <- postForm(url, query = q, .cgifields = c("query")) # submit query and fetch the resutls in a
    tmprslt <- readHTMLTable(tmprslt)

    ## parse
    tmpout <- NULL
    l <- length(tmprslt)

    if (l == 2) {
      tmpout <- tmprslt[[2]]
    } else if (l == 1) {
      warning(paste("No records returned for miRNA: ", i, " in database: ", j, ".", sep = ""))
    } else if (l == 0) {
      cat(paste("Request to multiMiR web server failed. check the ",
                "your query, or",
                "use the multiMiR web server:",
                "http://multimir.ucdenver.edu is temporarily down.\n"))
    }

    if (!is.null(tmpout)){
      tmpout <- data.frame(database = j, tmpout, stringsAsFactors = FALSE)
    }

    return(tmpout)
  }

  #### set up a temfunc for extract miRNA-mRNA entrez info for modelling from the output list
  # x - the dataframe from the output list
  # y - the output element (entrez ID) of the list
  tmpfunc_lst <- function(x){
    if (!is.null(x)){
      y <- as.character(unique(x$target_entrez))
      y[y == ""] <- NA
      y <- y[!is.na(y)]
      y <- unique(y)
    } else {
      y <- NA
    }
    return(y)
  }

  #### query and output
  out <- vector(mode = "list", length = length(mir)) # output list
  names(out) <- mir
  out_entrez <- vector(mode = "list", length = length(mir)) # output entrez list for modelling
  names(out_entrez) <- mir

  #### (optional) convert mmu/rno entrez ID to hsa entrez ID
  if (addhsaEntrez){
    # starting message
    cat(paste("Obtaining hsa orthologs information from ensembl databases for ", sp, "... May be slow depending on internet connectivity...", sep = ""))

    # set the target species
    if (sp == "mmu"){
      martsp <- "mmusculus"
    } else if (sp == "rno"){
      martsp <- "rnorvegicus"
    }

    # extract hsa ortholog information
    martsp_ensembl <- useMart("ensembl", dataset = paste0(martsp, "_gene_ensembl"))
    attr <- c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene")
    martsp_hsa_orth <- getBM(attr, filters = "with_hsapiens_homolog", values = TRUE,
                             mart = martsp_ensembl)
    names(martsp_hsa_orth)[names(martsp_hsa_orth) == "ensembl_gene_id"] <- paste0(sp, "_ensembl_gene_id") # generalized term for change column names

    # extract hsa entrezgene ID
    hsa_ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # establish the human set
    attr_hsa <- c("ensembl_gene_id", "entrezgene")
    hsa_entrez <- getBM(attr_hsa, filters = "", values = TRUE,
                        mart = hsa_ensembl)

    # merge the two dataframes
    martsp_hsa_orth_entrez <- merge(martsp_hsa_orth, hsa_entrez,
                                    by.x = "hsapiens_homolog_ensembl_gene", by.y = "ensembl_gene_id",
                                    all.x = TRUE)
    names(martsp_hsa_orth_entrez)[names(martsp_hsa_orth_entrez) == "entrezgene"] <- "hsa_entrezgene"
    martsp_hsa_orth_entrez <- martsp_hsa_orth_entrez[!duplicated(martsp_hsa_orth_entrez[, paste0(sp, "_ensembl_gene_id")]), ]

    # tmpfunc
    # d - input mmu/rno mRNA results dataframe
    # h - output vector containing hsa entrez
    tmpfunc_hsa_entrez <- function(d){
      if (!is.null(d)){
        h <- as.character(martsp_hsa_orth_entrez[martsp_hsa_orth_entrez[, paste0(sp, "_ensembl_gene_id")] %in% d$target_ensembl, "hsa_entrezgene"])
        h[h == ""] <- NA
        h <- h[!is.na(h)]
        h <- unique(h)
      } else {
        h <- NA
      }
      return(h)
    }

    # ending message
    cat("done!\n")

    # output list
    out_hsa_entrez <- vector(mode = "list", length = length(out_entrez))
    names(out_hsa_entrez) <- names(out_entrez)
  }

  if (!parallelComputing){
    ## populate output list and output
    out[] <- lapply(mir, function(m){
      tmpout <- foreach(n = db, .combine = rbind, .packages = c("RCurl", "XML")) %do% {
        cat(paste("searching ", n, " for ", m, " ...", sep = ""))
        tmp <- tmpfunc_q(i = m, j = n, mode = queryType, percentage = predictPercentage)
        cat("done!\n")
        return(tmp)}
    })

    # write mRNA results into files
    foreach(x = 1: length(out)) %do% {
      write.csv(out[[x]], file = paste(names(out)[x], "_mRNA.csv", sep = ""),  na = "NA", row.names = FALSE)
    }

    ## populate the entrez list
    out_entrez[] <- foreach(o = 1:length(out)) %do% tmpfunc_lst(out[[o]])

    ## populate the mmu/rno to hsa entrez ID list
    if (addhsaEntrez){
      out_hsa_entrez[] <- foreach(p = 1:length(out)) %do% tmpfunc_hsa_entrez(out[[p]])
    }

  } else { # parallel computing
    # set up cpu core number
    n_cores <- detectCores() - 1

    if (clusterType == "PSOCK"){ # for all OS systems
      ## set up cpu cluster for PSOCK
      cl <- makeCluster(n_cores, type = clusterType, outfile = "")
      registerDoParallel(cl) # part of doParallel package
      on.exit(stopCluster(cl)) # close connect when exiting the function

      ## populate the output list
      out[] <- foreach(m = mir, .packages = "foreach") %dopar% {
        tmpout <- foreach(n = db, .combine = rbind, .packages = c("RCurl", "XML")) %do% {
          cat(paste("searching ", n, " for ", m, " ...", sep = ""))
          tmp <- tmpfunc_q(i = m, j = n, mode = queryType, percentage = predictPercentage)
          cat("done!\n")
          return(tmp)}
      }

      # write mRNA results into files
      foreach(x = 1: length(out)) %dopar% {
        write.csv(out[[x]], file = paste(names(out)[x], "_mRNA.csv", sep = ""),  na = "NA", row.names = FALSE)
      }

      ## populate the entrez list
      out_entrez[] <- foreach(o = 1:length(out)) %dopar% tmpfunc_lst(out[[o]])

      ## populate the mmu/rno to hsa entrez ID list
      if (addhsaEntrez){
        out_hsa_entrez[] <- foreach(p = 1:length(out)) %dopar% tmpfunc_hsa_entrez(out[[p]])
      }

    } else if (clusterType == "FORK"){ # macOS and Unix-like systmes only
      ## message
      cat(paste("searching ", db, " ...", sep = ""))

      ## use mclapply from parallel pacakge for the FORK method to populate the ouput list
      out[] <- mclapply(mir, FUN = function(m){
        tmp <- foreach(n = db, .combine = rbind, .packages = c("RCurl", "XML")) %do% tmpfunc_q(m, n, mode = queryType, percentage = predictPercentage)
        return(tmp)
      }, mc.cores = n_cores, mc.preschedule = FALSE)

      ## message
      cat("done!\n")

      # write mRNA results into files
      mclapply(1:length(out), FUN = function(x){
        write.csv(out[[x]], file = paste(names(out)[x], "_mRNA.csv", sep = ""), na = "NA", row.names = FALSE)
      }, mc.cores = n_cores, mc.preschedule = FALSE)

      ## use mclapply to populate the entrez list
      out_entrez[] <- mclapply(mir, FUN = function(m){
        tmp <- foreach(o = 1:length(out)) %do% tmpfunc_lst(out[[o]])
        return(tmp)
      }, mc.cores = n_cores, mc.preschedule = FALSE)

      ## use mclapply to populate the hsa entrez ID list
      if (addhsaEntrez){
        out_hsa_entrez[] <- mclapply(mir, FUN = function(m){
          tmp <- foreach(p = 1:length(out)) %do% tmpfunc_hsa_entrez(out[[p]])
          return(tmp)
        }, mc.cores = n_cores, mc.preschedule = FALSE)
      }
    }
  }

  #### message
  if (addhsaEntrez){
    message(paste("...all done! And entrez ID for hsa orthologs added for ", sp, ".", sep = ""))
  } else {
    message("...all done!")
  }

  #### output
  ## the output list to the environment
  assign(paste(objTitle, "_mrna_list", sep = ""), out, envir = .GlobalEnv)
  ## the entrez list to the environment
  assign(paste(objTitle, "_mrna_entrez_list", sep = ""), out_entrez, envir = .GlobalEnv)
  ## the hsa entrez list to the environment
  if (addhsaEntrez){
    assign(paste(objTitle, "_mrna_hsa_entrez_list", sep = ""), out_hsa_entrez, envir = .GlobalEnv)
  }
}
