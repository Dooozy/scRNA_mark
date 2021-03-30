
###### Taken from utilities_internal.R ######
SetIfNull <- function(x, default) {
  if(is.null(x = x)){
    return(default)
  } else {
    return(x)
  }
}
SetCalcParams <- function(object, calculation, time = TRUE, ...) {
  object@calc.params[calculation] <- list(...)
  object@calc.params[[calculation]]$object <- NULL
  object@calc.params[[calculation]]$object2 <- NULL
  if(time) {
    object@calc.params[[calculation]]$time <- Sys.time()
  }
  return(object)
}
#############################################

# RegressOutResid_2 was changed to use doMC instead of doSNOW or parallelisation since it was throwing errors
# RegressOutResid is located at preprocessing_internal.R
RegressOutResid_2 <- function(
  object,
  vars.to.regress,
  genes.regress = NULL,
  model.use = 'linear',
  use.umi = FALSE,
  display.progress = TRUE,
  do.par = FALSE,
  num.cores = 1
) {
  possible.models <- c("linear", "poisson", "negbinom")
  if (!model.use %in% possible.models) {
    stop(
      paste0(
        model.use,
        " is not a valid model. Please use one the following: ",
        paste0(possible.models, collapse = ", "),
        "."
      )
    )
  }
  genes.regress <- SetIfNull(x = genes.regress, default = rownames(x = object@data))
  genes.regress <- intersect(x = genes.regress, y = rownames(x = object@data))
  latent.data <- FetchData(object = object, vars.all = vars.to.regress)
  bin.size <- ifelse(test = model.use == 'negbinom', yes = 5, no = 100)
  bin.ind <- ceiling(x = 1:length(x = genes.regress) / bin.size)
  max.bin <- max(bin.ind)
  if (display.progress) {
    message(paste("Regressing out:", paste(vars.to.regress, collapse = ", ")))
    pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
  }
  data.resid <- c()
  data.use <- object@data[genes.regress, , drop = FALSE];
  if (model.use != "linear") {
    use.umi <- TRUE
  }
  if (use.umi) {
    data.use <- object@raw.data[genes.regress, object@cell.names, drop = FALSE]
  }
  # input checking for parallel options
  if (do.par) {
    if (num.cores == 1) {
      num.cores <- detectCores() / 2
    } else if (num.cores > detectCores()) {
      num.cores <- detectCores() - 1
      warning(paste0("num.cores set greater than number of available cores(", detectCores(), "). Setting num.cores to ", num.cores, "."))
    }
  } else if (num.cores != 1) {
    num.cores <- 1
    warning("For parallel processing, please set do.par to TRUE.")
  }
  # cl <- parallel::makeCluster(num.cores)#, outfile = "")
  # # using doSNOW library because it supports progress bar update
  # registerDoSNOW(cl)
  opts <- list()
  if (display.progress) {
    # define progress bar function
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    time_elapsed <- Sys.time()
  }
  reg.mat.colnames <- c(colnames(x = latent.data), "GENE")
  fmla_str = paste0("GENE ", " ~ ", paste(vars.to.regress, collapse = "+"))
  if (model.use == "linear") {
    # In this code, we'll repeatedly regress different Y against the same X
    # (latent.data) in order to calculate residuals.  Rather that repeatedly
    # call lm to do this, we'll avoid recalculating the QR decomposition for the
    # latent.data matrix each time by reusing it after calculating it once and
    # storing it in a fastResiduals function.
    regression.mat <- cbind(latent.data, data.use[1,])
    colnames(regression.mat) <- reg.mat.colnames
    qr = lm(as.formula(fmla_str), data = regression.mat, qr = TRUE)$qr
    rm(regression.mat)
  }
  data.resid <- foreach(i = 1:max.bin, .combine = "c", .options.snow = opts) %do% {
    genes.bin.regress <- rownames(x = data.use)[bin.ind == i]
    gene.expr <- as.matrix(x = data.use[genes.bin.regress, , drop = FALSE])
    empty_char = character(length = dim(gene.expr)[1]) # Empty vector to reuse
    new.data <- sapply(
      X = genes.bin.regress,
      FUN = function(x) {
        # Fast path for std. linear models
        if(model.use=="linear") {
          resid <- qr.resid(qr, gene.expr[x,])
        } else {
          regression.mat <- cbind(latent.data, gene.expr[x,])
          colnames(x = regression.mat) <- reg.mat.colnames
          fmla = as.formula(fmla_str)
          resid <- switch(
            EXPR = model.use,
            'poisson' = residuals(
              object = glm(
                formula = fmla,
                data = regression.mat,
                family = "poisson"
              ),
              type = 'pearson'
            ),
            'negbinom' = NBResiduals(
              fmla = fmla,
              regression.mat = regression.mat,
              gene = x,
              return.mode = TRUE
            )
          )
        }
        if (!is.list(x = resid)) {
          resid <- list('resid' = resid, 'mode' = empty_char)
        }
        return(resid)
      }
    )
    new.data.resid <- new.data[seq.int(from = 1, to = length(x = new.data), by = 2)]
    new.data.resid = matrix(unlist(new.data.resid), nrow = length(new.data.resid[[1]]))
    colnames(x = new.data.resid) <- genes.bin.regress
    new.data.mode <- unlist(x = new.data[seq.int(from = 2, to = length(x = new.data), by = 2)])
    names(x = new.data.mode) <- genes.bin.regress
    new.data <- list('resid' = new.data.resid, 'mode' = new.data.mode)
    return(new.data)
  }
  if (display.progress) {
    time_elapsed <- Sys.time() - time_elapsed
    cat(paste("\nTime Elapsed: ",time_elapsed, units(time_elapsed)))
    close(pb)
  }
  # stopCluster(cl)
  modes <- unlist(x = data.resid[seq.int(from = 2, to = length(x = data.resid), by = 2)])
  modes <- modes[modes == 'scale']
  names(x = modes) <- gsub(
    pattern = 'mode.',
    replacement = '',
    x = names(x = modes),
    fixed = TRUE
  )
  data.resid <- data.resid[seq.int(from = 1, to = length(x = data.resid), by = 2)]
  data.resid <- as.matrix(x = as.data.frame(x = data.resid))
  data.resid <- t(x = data.resid)
  if (length(x = modes)) {
    message(
      "The following genes failed with glm.nb, and fell back to scale(log(y+1))\n\t",
      paste(names(x = modes), collapse = ', ')
    )
  }
  rownames(x = data.resid) <- genes.regress
  suppressWarnings(expr = gc(verbose = FALSE))
  if (use.umi) {
    data.resid <- log1p(
      x = sweep(
        x = data.resid,
        MARGIN = 1,
        STATS = apply(X = data.resid, MARGIN = 1, FUN = min),
        FUN = "-"
      )
    )
  }
  return(data.resid)
}

###### From RcppExports.R ######
FastRowScale <- function(mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
  .Call('_Seurat_FastRowScale', PACKAGE = 'Seurat', mat, scale, center, scale_max, display_progress)
}
FastSparseRowScale <- function(mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
  .Call('_Seurat_FastSparseRowScale', PACKAGE = 'Seurat', mat, scale, center, scale_max, display_progress)
}
#################################

# ScaleData_2 will only call the changed function RegressOutResid_2, nothing else was changed
# ScaleData is at preprocessing.R
ScaleData_2 <- function(
  object,
  genes.use = NULL,
  data.use = NULL,
  vars.to.regress,
  model.use = 'linear',
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  display.progress = TRUE,
  assay.type = "RNA",
  do.cpp = TRUE,
  check.for.norm = TRUE,
  do.par = FALSE,
  num.cores = 1
) {
  data.use <- SetIfNull(
    x = data.use,
    default = GetAssayData(
      object = object,
      assay.type = assay.type,
      slot = "data"
    )
  )

  if (check.for.norm) {
    if (!("NormalizeData" %in% names(object@calc.params))) {
      cat("NormalizeData has not been run, therefore ScaleData is running on non-normalized values. Recommended workflow is to run NormalizeData first.\n")
    }
    if (is.null(object@calc.params$NormalizeData$normalization.method)) {
      cat("ScaleData is running on non-normalized values. Recommended workflow is to run NormalizeData first.\n")
    }
  }
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  genes.use <- as.vector(
    x = intersect(
      x = genes.use,
      y = rownames(x = data.use)
    )
  )
  data.use <- data.use[genes.use, ]
  if (!missing(x = vars.to.regress) && !is.null(x = vars.to.regress)) {
    data.use <- RegressOutResid_2(
      object = object,
      vars.to.regress = vars.to.regress,
      genes.regress = genes.use,
      use.umi = use.umi,
      model.use = model.use,
      display.progress = display.progress,
      do.par = do.par,
      num.cores = num.cores
    )
    if (model.use != "linear") {
      use.umi <- TRUE
    }
    if (use.umi && missing(scale.max)) {
      scale.max <- 50
    }
  }
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("ScaleData"))]
  parameters.to.store$data.use <- NULL

  object <- SetCalcParams(
    object = object,
    calculation = "ScaleData",
    ... = parameters.to.store
  )
  if (!do.cpp) {
    return(ScaleDataR(
      object = object,
      data.use = data.use,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = scale.max,
      genes.use = genes.use
    ))
  }
  scaled.data <- matrix(
    data = NA,
    nrow = length(x = genes.use),
    ncol = ncol(x = object@data
    )
  )
  rownames(scaled.data) <- genes.use
  if (length(object@cell.names) <= min.cells.to.block) {
    block.size <- length(genes.use)
  }
  gc()
  colnames(scaled.data) <- colnames(object@data)
  max.block <- ceiling(x = length(x = genes.use) / block.size)
  gc()
  if (display.progress) {
    message("Scaling data matrix")
    pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
  }
  for (i in 1:max.block) {
    my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
    my.inds <- my.inds[my.inds <= length(x = genes.use)]
    if (class(x = data.use) == "dgCMatrix" | class(x = data.use) == "dgTMatrix") {
      data.scale <- FastSparseRowScale(
        mat = data.use[genes.use[my.inds], , drop = F],
        scale = do.scale,
        center = do.center,
        scale_max = scale.max,
        display_progress = FALSE
      )
    } else {
      data.scale <- FastRowScale(
        mat = as.matrix(x = data.use[genes.use[my.inds], , drop = F]),
        scale = do.scale,
        center = do.center,
        scale_max = scale.max,
        display_progress = FALSE
      )
    }
    dimnames(x = data.scale) <- dimnames(x = data.use[genes.use[my.inds], ])
    scaled.data[genes.use[my.inds], ] <- data.scale
    rm(data.scale)
    gc()
    if (display.progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (display.progress) {
    close(pb)
  }
  object <- SetAssayData(
    object = object,
    assay.type = assay.type,
    slot = 'scale.data',
    new.data = scaled.data
  )
  gc()
  object@scale.data[is.na(object@scale.data)] <- 0
  return(object)
}
