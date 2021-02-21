# jonashaslbeck@gmail.com; December 3, 2020

# -----------------------------------------------------
# -------- DataGen from FMs with varying loadings -----
# -----------------------------------------------------

f_datagen2 <- function(nfactor, # number of factors
                       items_pf, # items per factor
                       psi, # correlation matrix between factors
                       n, # sample size
                       ra = c(0.2, 1), # range of factor loadings
                       cats = NULL) # no of response categories

{

  ## Generate model specification

  # Draw factor loadings
  l_loadings <- list()
  for(i in 1:nfactor) l_loadings[[i]] <- round(runif(items_pf, ra[1], ra[2]), 2)

  # Each item loads on one factor
  spec_items <- rep(NA, nfactor)
  for(i in 1:nfactor) {
    base <- (0:(nfactor-1))[i] * items_pf
    spec_items[i] <- paste0("F", i, " =~ ", paste0(l_loadings[[i]], "*x", (base+1):(base + items_pf), collapse = " + "))
  }
  spec_items <- paste(spec_items, collapse = " \n ")
  full_model <- spec_items

  # Create factor correlations
  if(nfactor>1) {
    nfac_combn <- combn(1:nfactor, 2)
    n_combn <- ncol(nfac_combn)
    spec_psi <- rep(NA, n_combn)
    for(i in 1:n_combn) spec_psi[i] <- paste0("F", nfac_combn[1, i], " ~~", psi," * F", nfac_combn[2, i], " \n")
    spec_psi <- paste(spec_psi, collapse = "")

    full_model <- paste(spec_items, "\n", spec_psi)
  }

  # Sample data
  data <- simulateData(full_model,
                       sample.nobs = n)

  # Collapse into ordinal categories
  if(!is.null(cats)) {

    p <- ncol(data)

    data2 <- matrix(NA, n, p)
    for(i in 1:p) {
      data2[, i] <- as.numeric(cut(data[, i], breaks=cats, labels=1:cats))
    }
    data <- data2

  }


  # Create output list
  outlist <- list("l_loadings" = l_loadings,
                  "model_spec" = full_model,
                  "data" = data)

  # Return
  return(outlist)

} # eoF


