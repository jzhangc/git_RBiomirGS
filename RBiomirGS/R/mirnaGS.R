#' @title rbiomirgs_gmt
#'
#' @description load \code{gmt} gene set files as lists
#' @param file Input \code{gmt} file.
#' @return The function loads the \code{gmt} gene set database file and returns a \code{list} object.
#' @examples
#' \dontrun{
#' geneset <- rbiomirGS_gmt(file = "kegg.gmt")
#' }
#' @export
rbiomirgs_gmt <- function(file){
  # open connection
  gmt <- file(file)

  # check if the file can be read
  filecheck <- try(suppressWarnings(open(gmt)), silent = TRUE)
  if (class(filecheck) == "try-error") {
    stop("Bad gmt file.")
  } else {
    tmpfile <- scan(gmt, what = "", quiet = T, sep = "\n")
    tmplist <- strsplit(tmpfile, "\t")
    names(tmplist) <- sapply(tmplist, '[[', 1) # extract the gene set name as the names()
    outlist <- lapply(tmplist, '[', -c(1, 2)) # remove the annotation info
  }

  # close connection
  close(gmt)

  # output
  message(paste(length(outlist), " gene sets sucessfully loaded from gmz file.", sep = ""))
  return(outlist)
}

#' @title rbiomirgs_logistic
#'
#' @description Logistic regression-based gene set analysis using measured miRNA p value and fold change,  with the option of custom-setting parameter optimization algorithms.
#' @param objTitle Output GS result \code{csv} file name prefix. Default is \code{mirna_mrna}.
#' @param mirnascoreTitle Output miRNA score \code{csv} file name. Default is \code{mirnascore}.
#' @param mrnascoreTitle Output mRNA score \code{csv} file name. Default is \code{mrnascore}.
#' @param defile The input \code{csv} file containing miRNA list and DE resutls.
#' @param mirna_DE DE list of miRNAs of interest. This can be a \code{data.frame}, \code{matrix} or \code{list} object.
#' @param var_mirnaName Variable name for miRNA names in the DE list. Default is \code{"miRNA"}.
#' @param var_mirnaFC Variable name for miRNA fold change (or log transformed FC) in the DE list. Default is \code{"logFC"}.
#' @param ratioFC Whether the FC provided is a ratio value. Default is \code{FALSE}.
#' @param var_mirnaP Variable name for miRNA p value in the DE list. Default is \code{"p.value"}. Note that the value will be -log10 transformed before calculating the miRNA score.
#' @param mrnalist List containing the mRNA targets for the miRNAs of interest. This is a \code{list} object and can be obtained from \code{\link{rbiomirgs_mrnascan}} function.
#' @param mrna_Weight A vector weight for the miRNA-mRNA interaction. Default is \code{NULL}.
#' @param gs_file Input \code{gmt} for gene set, and can be obtained from \code{ensembl} databases.
#' @param optim_method The parameter optimization method for the logistic regression model. Options are \code{"L-BFGS-B"}, \code{"BFGS"}, and \code{"IWLS"}. Default is \code{"IWLS"}.
#' @param p.adj P value adjustment methods. Options are \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"}. Default is \code{"fdr"}.
#' @param ... Additional arguments for \code{optim} function.
#' @param parallelComputing If to use parallel computing or not. Default is \code{FALSE}.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Whether to display
#' @return A \code{data.frame} object containing GS results to to the environment, as well as a \code{csv} to the working directory. A \code{txt} file containing the iteration information will be generated if \code{opti_method = "L-BFGS-B" or "BFGS"}.
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @examples
#' \dontrun{
#' rbiomirgs_logistic(objTitle = "mirna_mrna",
#'                    mirna_DE = tstdfm2, var_mirnaName = "miRNA", var_mirnaFC = "logFC", var_mirnaP = "pvalue",
#'                    mrnalist = hsa_mrna_entrez_list_woNA, mrna_Weight = NULL,
#'                    gs_file = "~/OneDrive/my papers/my papers/potential_DRDC_paper 2 (diving)/data/kegg.v5.2.entrez.gmt",
#'                    optim_method = "L-BFGS-B", p.adj = "fdr",
#'                    parallelComputing = FALSE, clusterType = "PSOCK")

#' }
#' @export
rbiomirgs_logistic <- function(objTitle = "mirna_mrna",
                               mirnascoreTitle = "mirnascore",
                               mrnascoreTitle = "mrnascore",
                               defile = NULL,
                               mirna_DE = NULL, var_mirnaName = "miRNA", var_mirnaFC = "logFC", ratioFC = FALSE,
                               var_mirnaP = "p.value",
                               mrnalist = NULL, mrna_Weight = NULL,
                               gs_file = NULL,
                               optim_method = "IWLS", p.adj = "fdr",
                               ...,
                               parallelComputing = FALSE, clusterType = "PSOCK",
                               verbose = TRUE){
  #### check arguments
  if (is.null(mirna_DE)){
    stop("Please set the input object.")
  }

  for (var in c(var_mirnaName, var_mirnaFC, var_mirnaP)) {
    if (!var %in% names(mirna_DE)) stop(paste0(var, "not found in the input mirna_DE"))
  }

  if (any(class(mrnalist) %in% "mir_entrez_list")) {
    stop("mrnalist is a mir_entrez_list object, please use rbiomirgs_logisticV2().")
  }

  if (is.null(mrnalist) & class(mrnalist) != "list"){
    stop("Please set the proper mRNA target list. Currently, only list is supported.")
  }

  if (!optim_method %in% c("BFGS", "L-BFGS-B", "IWLS")){
    stop("Please set the proper optimization method. Options are \"L-BFGS-B\" (default), \"BFGS\" and \"IWLS\" ")
  }

  #### calculate the miRNA score
  if (is.null(defile)){
    if (class(mirna_DE) == "data.frame"){
      mirna.DE <- mirna_DE
      if (ratioFC){
        mirna.score <- sign(log2(mirna.DE[, var_mirnaFC])) * (-log10(mirna.DE[, var_mirnaP]))
      } else {
        mirna.score <- sign(mirna.DE[, var_mirnaFC]) * (-log10(mirna.DE[, var_mirnaP]))
      }
      names(mirna.score) <- mirna.DE[, var_mirnaName]
    } else if (class(mirna_DE == "matrix")){
      mirna.DE <- mirna_DE
      if (ratioFC){
        mirna.score <- sign(log2(as.numeric(mirna.DE[, var_mirnaFC]))) * (-log10(as.numeric(mirna.DE[, var_mirnaP])))
      } else {
        mirna.score <- sign(as.numeric(mirna.DE[, var_mirnaFC])) * (-log10(as.numeric(mirna.DE[, var_mirnaP])))
      }
      names(mirna.score) <- mirna.DE[, var_mirnaName]
    } else if (class(mirna_DE) == "list"){
      mirna.DE <- mirna_DE
      if (ratioFC){
        mirna.score <- sign(log2(mirna.DE[[var_mirnaFC]])) * (-log10(mirna.DE[[var_mirnaP]]))
      } else {
        mirna.score <- sign(mirna.DE[[var_mirnaFC]]) * (-log10(mirna.DE[[var_mirnaP]]))
      }
      names(mirna.score) <- mirna.DE[[var_mirnaName]]
    } else {
      stop("Currently, the input only supports dataframe, list or matrix")
    }
  } else {
    mirna.DE <- read.csv(file = defile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    mirna.score <- sign(mirna.DE[, var_mirnaFC]) * (-log10(mirna.DE[, var_mirnaP]))
    names(mirna.score) <- mirna.DE[, var_mirnaName]
  }

  #### prepare the mRNA list
  if (length(mrnalist[names(mrnalist) == ""]) > 0) {
    mrnalist <- mrnalist[names(mrnalist) != ""]  # remove empty miRNA symbols
    if (verbose) cat("Empty miRNA ID entires removed. \n")
  }
  if (length(mrnalist[!is.na(mrnalist)]) > 0) {
    mrna.raw.list <- mrnalist[!is.na(mrnalist)] # remove the NAs
    if (verbose) cat("miRNAs without mRNA targets removed. \n")
  }
  mrna <- sort(unique(unlist(mrna.raw.list))) # all mRNAs identified as miRNA targets

  mirna.in.DE <- as.character(mirna.DE[, var_mirnaName]) # extract miRNAs present in the DE results
  if (length(mirna.in.DE[mirna.in.DE == ""]) > 0) {
    mirna.in.DE <- mirna.in.DE[mirna.in.DE != ""]  # remove any empty entries
    if (verbose) cat("mRNAs without an miRNA ID removed. \n\n")
  }

  mirna.in.mrna <- names(mrna.raw.list) # extract miRNAs with mRNA target information
  mirna.working <- intersect(mirna.in.DE, mirna.in.mrna) # extract miRNAs from DE resutls with mRNA target info

  # message to show all the missing miRNAs (the ones without mRNA information)

  #### calculate mRNA score
  ## give product scaler for the scores: 1L for the miRNAs that are the upstream miRNA to the mRNA tagets (rows) (1L: TRUE, 0L: FALSE)
  mat <- foreach(i = mirna.working, .combine = cbind) %do% {tmp <- as.numeric(mrna %in% mrna.raw.list[[i]])} # logic to numeric
  rownames(mat) <- mrna
  colnames(mat) <- mirna.working

  ## calculate
  # the step below uses linear algebra, which is the fastest. other options are t(t(mat) * mirna.score[mirna.working]),
  # and sweep(mat, MARGIN = 2, mirna.score[mirna.working], "*") (slowest)
  mat <- mat %*% diag(mirna.score[mirna.working]) # index[mirnas] is to obtain the score for only the miRNAs with BOTH DE and downstream mRNA information
  colnames(mat) <- mirna.working

  if (is.null(mrna_Weight)){
    mrna.score <- -rowSums(mat) # reversed sign from miRNA to mRNA. Positive number means activation on mRNA and GS from this point on.
  } else {
    w <- as.matrix(mrna_Weight)
    if (identical(dim(mat), dim(w))){
      mat_w <- mat * w
      colnames(mat_w) <- mirna.working
      rownames(mat_w) <- mrna
      mrna.score <- -rowSums(mat_w) # reversed sign from miRNA to mRNA. Positive number means activation on mRNA and GS from this point on.
    } else {
      stop(cat("The miRNA:mRNA interaction weight matrix doesn't match the dimension of the miRNA:mRNA score matrix. Please check."))
    }
  }

  names(mrna.score) <- rownames(mat)
  mrna.score <- as.matrix(mrna.score)

  #### set up GS
  GS <- rbiomirgs_gmt(file = gs_file)
  gs_names <- names(GS) # extract GS names

  #### logistic modelling
  ## tmpfuncs
  # cost function
  cost <- function(vTh, mX, vY){
    vH <- 1 / (1 + exp(-mX %*% vTh)) # logit
    cst <- 1 / nrow(mX) * sum(-vY * log(vH) - (1 - vY) * log(1 - vH))
    return(cst)
  }

  # set up tmpfunc for not IWLS methods
  if (optim_method == "IWLS"){
    tmpfunc_IWLS <- function(i){

      ## setup y
      y <- as.numeric(mrna %in% GS[[i]])

      # message
      if (verbose) cat(paste("Assessing gene set: ", i, "...", sep = ""))

      model <- glm.fit(x = x, y = y, family = quasibinomial())
      model.smry <- summary.glm(model)

      # message
      if (verbose) cat("done!\n")

      ## gather results
      genes.tested <- sum(y)
      coef <- model.smry$coefficients[, 1] # extract coef for loss calculation
      loss <- cost(vTh = coef, mX = x, vY = y)

      ## output
      tmp <- c(model$converged, loss, genes.tested, model.smry$coefficients[2, ])
    }

  } else { # set up tmpfuncs for optm methods
    # logit (vectorized function, x is a matrix aand vTh is a vector)
    vlogit <- function(mX, vTh){
      vH <- 1 / (1 + exp(-mX %*% vTh))
      return(vH)
    }

    # log likelihood
    # vY -  1 and 0s (label)
    vloglike <- function(mX, vY, vTh){
      vH <- vlogit(mX, vTh) # get the estimated probablity
      logvH <- log(vH) # log probability
      vY0 <- 1 - vY # 1 - y
      logvH0 <- log(1 - vH) # log 1 - P
      vY_t <- t(vY) # transpose. Check if this is really needed
      vY0_t <- t(vY0) # transpose Check if this is really needed

      # below is negative log likelihood, or the loss function for logit.
      # the goal is to minimize this (or maximize the positive likelihood)
      logl <- -sum(vY_t %*% logvH + vY0_t %*% logvH0)
      return(logl)
    }

    if (optim_method %in% c("BFGS","L-BFGS-B")){
      vGr <- function(mX, vY, vTh){
        grad <- vTh * 0
        vH <- vlogit(mX, vTh) # get the estimated probablity
        for (i in 1:ncol(mX)){ # for the variables only
          grad[i] <- sum(mX[, i] * (vY - vH))
        }
        return(-grad)
      }
    } else {
      vGr <- NULL
    }

    tmpfunc_optim <- function(i, altX, ...){
      # set up vY
      y <- as.numeric(mrna %in% GS[[i]])

      # initialize parameters
      startvalue.model <- lm(y ~ altX) # note that since we use lm function, we don't need the x0 = 1 term here. So we use altX.
      startv <- startvalue.model$coefficients

      # message
      if (verbose) cat(paste("Assessing gene set: ", i, "...", sep = ""))

      # modelling/optimization
      sink(file = paste(objTitle, "_", optim_method, "_log.txt", sep = ""), append = TRUE) # dump iteration messages
      if (verbose) cat(paste(i, " \n", sep = ""))
      vTh.optim <- optim(par = startv, fn = vloglike, gr = vGr, mX = x, vY = y,
                         method = optim_method,
                         control = list(trace = TRUE, REPORT = 1), hessian = TRUE, ...)
      if (verbose) cat("\n\n") # add a new line between gene sets
      sink() # close dump

      # message
      if (verbose) cat("done!\n")

      # gather the resutls
      genes.tested <- sum(y)
      coef <- vTh.optim$par
      covmat <- solve(vTh.optim$hessian)
      stderr <- sqrt(diag(covmat))
      loss <- cost(vTh = coef, mX = x, vY = y)
      z <- coef / stderr
      pvalue <- 2*(1 - pnorm(abs(z)))
      model.sum <- cbind(genes.tested, coef, stderr, loss, z, pvalue)
      tmp <- model.sum[2 , ]
    }
  }

  ## logistic regression
  # set up variables
  x <- cbind(rep(1, length(mrna)), mrna.score) # set up vairable
  altX <- x[, 2] # this for automatically setting the initial values for the parameter

  ## modelling/optimization
  # set up result matrix
  if (optim_method == "IWLS"){
    res <- matrix(NA, nrow = length(gs_names), ncol = 7)
  } else {
    res <- matrix(NA, nrow = length(gs_names), ncol = 6)
  }

  if (!parallelComputing){
    if (optim_method == "IWLS"){
      res[] <- foreach(i = gs_names, .combine = rbind) %do% tmpfunc_IWLS(i)
    } else {
      res[] <- foreach(i = gs_names, .combine = rbind) %do% tmpfunc_optim(i, altX = altX)
    }
  } else {
    # set up cpu core number
    n_cores <- detectCores() - 1

    if (clusterType == "PSOCK"){ # for all OS systems
      ## set up cpu cluster for PSOCK
      cl <- makeCluster(n_cores, type = clusterType, outfile = "")
      registerDoParallel(cl) # part of doParallel package
      on.exit(stopCluster(cl)) # close connect when exiting the function

      ## modelling
      if (optim_method == "IWLS"){
        res[] <- foreach(i = gs_names, .combine = rbind) %dopar% tmpfunc_IWLS(i)
      } else {
        res[] <- foreach(i = gs_names, .combine = rbind) %dopar% tmpfunc_optim(i, altX = altX)
      }
    } else { # macOS and Unix-like systems only
      # message
      if (verbose) cat("Assessing gene sets...")

      if (optim_method == "IWLS"){
        res <- as.data.frame(do.call(rbind, mclapply(gs_names, FUN = tmpfunc_IWLS, mc.cores = n_cores, mc.preschedule = FALSE)))
      } else {
        res <- as.data.frame(do.call(rbind, mclapply(gs_names, FUN = tmpfunc_optim, altX = altX, mc.cores = n_cores, mc.preschedule = FALSE)))
      }
      # message
      if (verbose) cat("done!\n")
    }
  }

  #### add matrix information and export
  rownames(res) <- gs_names

  if (optim_method == "IWLS"){
    colnames(res) <- c("converged", "loss", "gene.tested", "coef", "std.err", "t.value", "p.value")
  } else {
    colnames(res) <- c("gene.tested", "coef", "std.err", "loss", "z.score", "p.value") # see if to add "converged"
  }

  ## p value adjustment
  out <- data.frame(GS = rownames(res), res, adj.p.val = p.adjust(res[, "p.value"], method = p.adj))

  if (optim_method == "IWLS") {
    out$converged <- factor(out$converged, levels = c(1, 0), labels = c("Y", "N"))
  }

  ## output
  mrna.score.out <- data.frame(EntrezID = rownames(mrna.score), S_mrna = mrna.score)
  mirna.score.out <-data.frame(miRNA = names(mirna.score), S_mirna = mirna.score)

  write.csv(out, file = paste(objTitle, "_GS.csv", sep = ""), na = "NA", row.names = FALSE)
  write.csv(mrna.score.out, file = paste(mrnascoreTitle, ".csv", sep = ""), na = "NA", row.names = FALSE)
  write.csv(mirna.score.out, file = paste(mirnascoreTitle, ".csv", sep = ""), na = "NA", row.names = FALSE)
  assign(paste(objTitle, "_GS", sep = ""), out, envir = .GlobalEnv)
  assign(mrnascoreTitle, mrna.score.out, envir = .GlobalEnv)
  assign(mirnascoreTitle, mirna.score.out, envir = .GlobalEnv)
}


#' @title rbiomirgs_logisticV2
#'
#' @description The V2 version fo the \code{\link{rbiomirgs_logistic}} function that takes the \code{mir_entrez_list} class for \code{mrnalist}.
#' @param objTitle Output GS result \code{csv} file name prefix. Default is \code{mirna_mrna}.
#' @param mirnascoreTitle Output miRNA score \code{csv} file name. Default is \code{mirnascore}.
#' @param mrnascoreTitle Output mRNA score \code{csv} file name. Default is \code{mrnascore}.
#' @param defile The input \code{csv} file containing miRNA list and DE resutls.
#' @param mirna_DE DE list of miRNAs of interest. This can be a \code{data.frame}, \code{matrix} or \code{list} object.
#' @param var_mirnaName Variable name for miRNA names in the DE list. Default is \code{"miRNA"}.
#' @param var_mirnaFC Variable name for miRNA fold change (or log transformed FC) in the DE list. Default is \code{"logFC"}.
#' @param ratioFC Whether the FC provided is a ratio value. Default is \code{FALSE}.
#' @param var_mirnaP Variable name for miRNA p value in the DE list. Default is \code{"p.value"}. Note that the value will be -log10 transformed before calculating the miRNA score.
#' @param mrnalist A \code{mir_entrez_list} class containing the entrez ID of the mRNA targets for the miRNAs of interest. This is a \code{list} object and can be obtained from \code{\link{rbiomirgs_mrnascan}} function.
#' @param mrna_Weight A vector weight for the miRNA-mRNA interaction. Default is \code{NULL}.
#' @param gs_file Input \code{gmt} for gene set, and can be obtained from \code{ensembl} databases.
#' @param optim_method The parameter optimization method for the logistic regression model. Options are \code{"L-BFGS-B"}, \code{"BFGS"}, and \code{"IWLS"}. Default is \code{"IWLS"}.
#' @param p.adj P value adjustment methods. Options are \code{"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"}. Default is \code{"fdr"}.
#' @param ... Additional arguments for \code{optim} function.
#' @param parallelComputing If to use parallel computing or not. Default is \code{FALSE}.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Whether to display
#' @return A \code{data.frame} object containing GS results to to the environment, as well as a \code{csv} to the working directory. A \code{txt} file containing the iteration information will be generated if \code{opti_method = "L-BFGS-B" or "BFGS"}.
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @examples
#' \dontrun{
#' rbiomirgs_logistic(objTitle = "mirna_mrna",
#'                    mirna_DE = tstdfm2, var_mirnaName = "miRNA", var_mirnaFC = "logFC", var_mirnaP = "pvalue",
#'                    mrnalist = hsa_mrna_entrez_list_woNA, mrna_Weight = NULL,
#'                    gs_file = "~/OneDrive/my papers/my papers/potential_DRDC_paper 2 (diving)/data/kegg.v5.2.entrez.gmt",
#'                    optim_method = "L-BFGS-B", p.adj = "fdr",
#'                    parallelComputing = FALSE, clusterType = "PSOCK")

#' }
#' @export
rbiomirgs_logisticV2 <- function(objTitle = "mirna_mrna",
                               mirnascoreTitle = "mirnascore",
                               mrnascoreTitle = "mrnascore",
                               defile = NULL,
                               mirna_DE = NULL, var_mirnaName = "miRNA", var_mirnaFC = "logFC", ratioFC = FALSE,
                               var_mirnaP = "p.value",
                               mrnalist = NULL, mrna_Weight = NULL,
                               gs_file = NULL,
                               optim_method = "IWLS", p.adj = "fdr",
                               ...,
                               parallelComputing = FALSE, clusterType = "PSOCK",
                               verbose = TRUE){
  #### check arguments
  if (is.null(mirna_DE)){
    stop("Please set the input object.")
  }

  for (var in c(var_mirnaName, var_mirnaFC, var_mirnaP)) {
    if (!var %in% names(mirna_DE)) stop(paste0(var, "not found in the input mirna_DE"))
  }

  # if (is.null(mrnalist) & class(mrnalist) != "list"){
  #   stop("Please set the proper mRNA target list. Currently, only list is supported.")
  # }
  if (!any(class(mrnalist) %in% "mir_entrez_list")) stop("The mrnalist must be a \"mir_entrez_list\" class.")

  if (!optim_method %in% c("BFGS", "L-BFGS-B", "IWLS")){
    stop("Please set the proper optimization method. Options are \"L-BFGS-B\" (default), \"BFGS\" and \"IWLS\" ")
  }

  #### calculate the miRNA score
  if (is.null(defile)){
    if (class(mirna_DE) == "data.frame"){
      mirna.DE <- mirna_DE
      if (ratioFC){
        mirna.score <- sign(log2(mirna.DE[, var_mirnaFC])) * (-log10(mirna.DE[, var_mirnaP]))
      } else {
        mirna.score <- sign(mirna.DE[, var_mirnaFC]) * (-log10(mirna.DE[, var_mirnaP]))
      }
      names(mirna.score) <- mirna.DE[, var_mirnaName]
    } else if (class(mirna_DE == "matrix")){
      mirna.DE <- mirna_DE
      if (ratioFC){
        mirna.score <- sign(log2(as.numeric(mirna.DE[, var_mirnaFC]))) * (-log10(as.numeric(mirna.DE[, var_mirnaP])))
      } else {
        mirna.score <- sign(as.numeric(mirna.DE[, var_mirnaFC])) * (-log10(as.numeric(mirna.DE[, var_mirnaP])))
      }
      names(mirna.score) <- mirna.DE[, var_mirnaName]
    } else if (class(mirna_DE) == "list"){
      mirna.DE <- mirna_DE
      if (ratioFC){
        mirna.score <- sign(log2(mirna.DE[[var_mirnaFC]])) * (-log10(mirna.DE[[var_mirnaP]]))
      } else {
        mirna.score <- sign(mirna.DE[[var_mirnaFC]]) * (-log10(mirna.DE[[var_mirnaP]]))
      }
      names(mirna.score) <- mirna.DE[[var_mirnaName]]
    } else {
      stop("Currently, the input only supports dataframe, list or matrix")
    }
  } else {
    mirna.DE <- read.csv(file = defile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    mirna.score <- sign(mirna.DE[, var_mirnaFC]) * (-log10(mirna.DE[, var_mirnaP]))
    names(mirna.score) <- mirna.DE[, var_mirnaName]
  }

  #### prepare the mRNA list
  if (length(mrnalist[names(mrnalist) == ""]) > 0) {
    mrnalist <- mrnalist[names(mrnalist) != ""]  # remove empty miRNA symbols
    if (verbose) cat("Empty miRNA ID entires removed. \n")
  }
  if (length(mrnalist[!is.na(mrnalist)]) > 0) {
    mrna.raw.list <- mrnalist[!is.na(mrnalist)] # remove the NAs
    if (verbose) cat("miRNAs without mRNA targets removed. \n")
  }
  mrna <- sort(unique(unlist(mrna.raw.list))) # all mRNAs identified as miRNA targets

  mirna.in.DE <- as.character(mirna.DE[, var_mirnaName]) # extract miRNAs present in the DE results
  if (length(mirna.in.DE[mirna.in.DE == ""]) > 0) {
    mirna.in.DE <- mirna.in.DE[mirna.in.DE != ""]  # remove any empty entries
    if (verbose) cat("mRNAs without an miRNA ID removed. \n\n")
  }

  mirna.in.mrna <- names(mrna.raw.list) # extract miRNAs with mRNA target information
  mirna.working <- intersect(mirna.in.DE, mirna.in.mrna) # extract miRNAs from DE resutls with mRNA target info

  # message to show all the missing miRNAs (the ones without mRNA information)

  #### calculate mRNA score
  ## give product scaler for the scores: 1L for the miRNAs that are the upstream miRNA to the mRNA tagets (rows) (1L: TRUE, 0L: FALSE)
  mat <- foreach(i = mirna.working, .combine = cbind) %do% {tmp <- as.numeric(mrna %in% mrna.raw.list[[i]])} # logic to numeric
  rownames(mat) <- mrna
  colnames(mat) <- mirna.working

  ## calculate
  # the step below uses linear algebra, which is the fastest. other options are t(t(mat) * mirna.score[mirna.working]),
  # and sweep(mat, MARGIN = 2, mirna.score[mirna.working], "*") (slowest)
  mat <- mat %*% diag(mirna.score[mirna.working]) # index[mirnas] is to obtain the score for only the miRNAs with BOTH DE and downstream mRNA information
  colnames(mat) <- mirna.working

  if (is.null(mrna_Weight)){
    mrna.score <- -rowSums(mat) # reversed sign from miRNA to mRNA. Positive number means activation on mRNA and GS from this point on.
  } else {
    w <- as.matrix(mrna_Weight)
    if (identical(dim(mat), dim(w))){
      mat_w <- mat * w
      colnames(mat_w) <- mirna.working
      rownames(mat_w) <- mrna
      mrna.score <- -rowSums(mat_w) # reversed sign from miRNA to mRNA. Positive number means activation on mRNA and GS from this point on.
    } else {
      stop(cat("The miRNA:mRNA interaction weight matrix doesn't match the dimension of the miRNA:mRNA score matrix. Please check."))
    }
  }

  names(mrna.score) <- rownames(mat)
  mrna.score <- as.matrix(mrna.score)

  #### set up GS
  GS <- rbiomirgs_gmt(file = gs_file)
  gs_names <- names(GS) # extract GS names

  #### logistic modelling
  ## tmpfuncs
  # cost function
  cost <- function(vTh, mX, vY){
    vH <- 1 / (1 + exp(-mX %*% vTh)) # logit
    cst <- 1 / nrow(mX) * sum(-vY * log(vH) - (1 - vY) * log(1 - vH))
    return(cst)
  }

  # set up tmpfunc for not IWLS methods
  if (optim_method == "IWLS"){
    tmpfunc_IWLS <- function(i){

      ## setup y
      y <- as.numeric(mrna %in% GS[[i]])

      # message
      if (verbose) cat(paste("Assessing gene set: ", i, "...", sep = ""))

      model <- glm.fit(x = x, y = y, family = quasibinomial())
      model.smry <- summary.glm(model)

      # message
      if (verbose) cat("done!\n")

      ## gather results
      genes.tested <- sum(y)
      coef <- model.smry$coefficients[, 1] # extract coef for loss calculation
      loss <- cost(vTh = coef, mX = x, vY = y)

      ## output
      tmp <- c(model$converged, loss, genes.tested, model.smry$coefficients[2, ])
    }

  } else { # set up tmpfuncs for optm methods
    # logit (vectorized function, x is a matrix aand vTh is a vector)
    vlogit <- function(mX, vTh){
      vH <- 1 / (1 + exp(-mX %*% vTh))
      return(vH)
    }

    # log likelihood
    # vY -  1 and 0s (label)
    vloglike <- function(mX, vY, vTh){
      vH <- vlogit(mX, vTh) # get the estimated probablity
      logvH <- log(vH) # log probability
      vY0 <- 1 - vY # 1 - y
      logvH0 <- log(1 - vH) # log 1 - P
      vY_t <- t(vY) # transpose. Check if this is really needed
      vY0_t <- t(vY0) # transpose Check if this is really needed

      # below is negative log likelihood, or the loss function for logit.
      # the goal is to minimize this (or maximize the positive likelihood)
      logl <- -sum(vY_t %*% logvH + vY0_t %*% logvH0)
      return(logl)
    }

    if (optim_method %in% c("BFGS","L-BFGS-B")){
      vGr <- function(mX, vY, vTh){
        grad <- vTh * 0
        vH <- vlogit(mX, vTh) # get the estimated probablity
        for (i in 1:ncol(mX)){ # for the variables only
          grad[i] <- sum(mX[, i] * (vY - vH))
        }
        return(-grad)
      }
    } else {
      vGr <- NULL
    }

    tmpfunc_optim <- function(i, altX, ...){
      # set up vY
      y <- as.numeric(mrna %in% GS[[i]])

      # initialize parameters
      startvalue.model <- lm(y ~ altX) # note that since we use lm function, we don't need the x0 = 1 term here. So we use altX.
      startv <- startvalue.model$coefficients

      # message
      if (verbose) cat(paste("Assessing gene set: ", i, "...", sep = ""))

      # modelling/optimization
      sink(file = paste(objTitle, "_", optim_method, "_log.txt", sep = ""), append = TRUE) # dump iteration messages
      if (verbose) cat(paste(i, " \n", sep = ""))
      vTh.optim <- optim(par = startv, fn = vloglike, gr = vGr, mX = x, vY = y,
                         method = optim_method,
                         control = list(trace = TRUE, REPORT = 1), hessian = TRUE, ...)
      if (verbose) cat("\n\n") # add a new line between gene sets
      sink() # close dump

      # message
      if (verbose) cat("done!\n")

      # gather the resutls
      genes.tested <- sum(y)
      coef <- vTh.optim$par
      covmat <- solve(vTh.optim$hessian)
      stderr <- sqrt(diag(covmat))
      loss <- cost(vTh = coef, mX = x, vY = y)
      z <- coef / stderr
      pvalue <- 2*(1 - pnorm(abs(z)))
      model.sum <- cbind(genes.tested, coef, stderr, loss, z, pvalue)
      tmp <- model.sum[2 , ]
    }
  }

  ## logistic regression
  # set up variables
  x <- cbind(rep(1, length(mrna)), mrna.score) # set up vairable
  altX <- x[, 2] # this for automatically setting the initial values for the parameter

  ## modelling/optimization
  # set up result matrix
  if (optim_method == "IWLS"){
    res <- matrix(NA, nrow = length(gs_names), ncol = 7)
  } else {
    res <- matrix(NA, nrow = length(gs_names), ncol = 6)
  }

  if (!parallelComputing){
    if (optim_method == "IWLS"){
      res[] <- foreach(i = gs_names, .combine = rbind) %do% tmpfunc_IWLS(i)
    } else {
      res[] <- foreach(i = gs_names, .combine = rbind) %do% tmpfunc_optim(i, altX = altX)
    }
  } else {
    # set up cpu core number
    n_cores <- detectCores() - 1

    if (clusterType == "PSOCK"){ # for all OS systems
      ## set up cpu cluster for PSOCK
      cl <- makeCluster(n_cores, type = clusterType, outfile = "")
      registerDoParallel(cl) # part of doParallel package
      on.exit(stopCluster(cl)) # close connect when exiting the function

      ## modelling
      if (optim_method == "IWLS"){
        res[] <- foreach(i = gs_names, .combine = rbind) %dopar% tmpfunc_IWLS(i)
      } else {
        res[] <- foreach(i = gs_names, .combine = rbind) %dopar% tmpfunc_optim(i, altX = altX)
      }
    } else { # macOS and Unix-like systems only
      # message
      if (verbose) cat("Assessing gene sets...")

      if (optim_method == "IWLS"){
        res <- as.data.frame(do.call(rbind, mclapply(gs_names, FUN = tmpfunc_IWLS, mc.cores = n_cores, mc.preschedule = FALSE)))
      } else {
        res <- as.data.frame(do.call(rbind, mclapply(gs_names, FUN = tmpfunc_optim, altX = altX, mc.cores = n_cores, mc.preschedule = FALSE)))
      }
      # message
      if (verbose) cat("done!\n")
    }
  }

  #### add matrix information and export
  rownames(res) <- gs_names

  if (optim_method == "IWLS"){
    colnames(res) <- c("converged", "loss", "gene.tested", "coef", "std.err", "t.value", "p.value")
  } else {
    colnames(res) <- c("gene.tested", "coef", "std.err", "loss", "z.score", "p.value") # see if to add "converged"
  }

  ## p value adjustment
  out <- data.frame(GS = rownames(res), res, adj.p.val = p.adjust(res[, "p.value"], method = p.adj))

  if (optim_method == "IWLS") {
    out$converged <- factor(out$converged, levels = c(1, 0), labels = c("Y", "N"))
  }

  ## output
  mrna.score.out <- data.frame(EntrezID = rownames(mrna.score), S_mrna = mrna.score)
  mirna.score.out <-data.frame(miRNA = names(mirna.score), S_mirna = mirna.score)

  write.csv(out, file = paste(objTitle, "_GS.csv", sep = ""), na = "NA", row.names = FALSE)
  write.csv(mrna.score.out, file = paste(mrnascoreTitle, ".csv", sep = ""), na = "NA", row.names = FALSE)
  write.csv(mirna.score.out, file = paste(mirnascoreTitle, ".csv", sep = ""), na = "NA", row.names = FALSE)
  assign(paste(objTitle, "_GS", sep = ""), out, envir = .GlobalEnv)
  assign(mrnascoreTitle, mrna.score.out, envir = .GlobalEnv)
  assign(mirnascoreTitle, mirna.score.out, envir = .GlobalEnv)
}
