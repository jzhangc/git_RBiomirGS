.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Written by Jing Zhang, Ph.D. Please direct questions to jzhangcad@gmail.com.")
  suppressPackageStartupMessages(require(pathview))
  return(TRUE)
}

#' @title rbiomirGS_mrnalist
#'
#' @description Obtain target mRNA list for miRNAs of interest. Resuls can be either predicted or validated. The function uses multiple databsaes hosted at \code{multimir.ucdenver.edu/}.
#' @param mir Input miRNAs vector.
#' @param sp Species. Options are \code{"hsa"} (default), \code{"mmu"} and \code{"rno"}.
#' @param queryType Type of reuslts. Options are \code{"validated"} and \code{"predicted"}.
#' @param predictPrecentage Set only if \code{queryType = "predicted"}. The percentage of the top scored mRNA predictions to return. Default is \code{5}.
#' @param url The database host server: "http://multimir.ucdenver.edu/cgi-bin/multimir.pl"
#' @param parallelComputing If to use parallel computing or not. Default is \code{FALSE}.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @return The function returns a \code{list} object with mRNA targets for the query miRNAs.
#' @import doParallel
#' @import foreach
#' @importFrom XML readHTMLTable
#' @importFrom RCurl postForm
#' @importFrom parallel detectCores makeCluster stopCluster
#' @examples
#' \dontrun{
#' mRNA <- rbiomirGS_mrnalist(mir = c("hsa-miR-26a-5p", "hsa-miR-100-5p"), sp = "hsa", queryType = "predicted", predictPercentage = 10)
#' }
#' @export
rbiomirGS_mrnalist <- function(mir =  NULL, sp = "hsa",
                               queryType = NULL, predictPercentage = 5,
                               url = "http://multimir.ucdenver.edu/cgi-bin/multimir.pl",
                               parallelComputing = FALSE, clusterType = "PSOCK"){

  #### check the arguments
  if (is.null(mir)){
    stop(cat("Please set the input miRNA(s). Either single targets or a vector of mutliple ones can be used."))
  }


  if (sp != "hsa" & sp != "mmu" & sp != "rno"){ # check species
    stop(cat("Only human, mouse or rat are supported for now. Please choose either \"hsa\" (default), \"mmu\", or \"rno\" for species."))
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

  #### set up a tmpfunc for query, modified from multiMiR pacakge.
  # i - target miRNA
  # j - single database
  # mode - queryType ("validated", "predicted")
  # percenetage - predictPercentage
  tmpfunc <- function(i, j, mode = NULL, percentage = NULL){

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
                 mirna.table, "AS m INNER JOIN", table, "AS i INNER JOIN",
                 target.table, "AS t ON (m.mature_mirna_uid=i.mature_mirna_uid",
                 "AND i.target_uid=t.target_uid) WHERE",
                 sep = " ")

      if (table == "diana_microt") {
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

      if (!is.na(score.cutoff)) {
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
      cat(paste("Request to multiMiR web server faile. check the ",
                "your query, or",
                "use the multiMiR web server:",
                "http://multimir.ucdenver.edu is temporarily down.\n"))
    }

    if (!is.null(tmpout)){
      tmpout <- data.frame(database = j, tmpout, stringsAsFactors = FALSE)
    }

    return(tmpout)
  }

  #### query
  out <- vector(mode = "list", length = length(mir))
  names(out) <- mir

  if (!parallelComputing){

    out[] <- lapply(mir, function(m){
      tmp <- foreach(n = db, .combine = rbind, .packages = c("RCurl", "XML")) %do% tmpfunc(i = m, j = n, mode = queryType, percentage = predictPercentage)
      return(tmp)
    })

  } else { # parallel computing
    # set up cpu core number
    n_cores <- detectCores() - 1


    if (clusterType == "PSOCK"){ # for all OS systems
      # set up cpu cluster for PSOCK
      cl <- makeCluster(n_cores, type = clusterType, outfile = "")
      registerDoParallel(cl) # part of doParallel package
      on.exit(stopCluster(cl)) # close connect when exiting the function

      out[] <- foreach(m = mir, .packages = "foreach") %dopar% {
        tmpout <- foreach(n = db, .combine = rbind, .packages = c("RCurl", "XML")) %do% tmpfunc(m, n, mode = queryType, percentage = predictPercentage)
      }

    } else if (clusterType == "FORK"){ # macOS and Unix-like systmes only
      # use mclapply from parallel pacakge for the FORK method
      out[] <- mclapply(mir, FUN = function(m){
        tmp <- foreach(n = db, .combine = rbind, .packages = c("RCurl", "XML")) %do% tmpfunc(m, n, mode = queryType, percentage = predictPercentage)
        return(tmp)
      }, mc.cores = n_cores, mc.preschedule = FALSE)

    }
  }
  return(out)
}
