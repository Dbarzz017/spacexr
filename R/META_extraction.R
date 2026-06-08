#estimates extraction: 
# first is per-rep extraction, then wrapper. 
extract_estimates <- function(RCTD, rep_id) {
  #start actual function code here
  gene_fits <- RCTD@de_results$gene_fits #this has a majority of the useful information we need to extract
  n_params <- dim(gene_fits$all_vals)[[2]] # this covers only X2, which is fine since any X1
  cell_types <- dimnames(gene_fits$all_vals)[[3]]
  param_names <- colnames(RCTD@internal_vars_de$X2) #covers all the parameters/covariates in our cell-type specific design matrix, not just params_to_test

  rows_to_add <- list() # empty list of rows to rbind into estimate_df
  row_i <- 1 # count/indexing variable 

  for (cell_type in cell_types) {
      ct_index <- match(cell_type, cell_types) # get the numerical index

      # this only keeps genes which converged for our particular cell_type 
      for (gene in rownames(gene_fits$con_mat)[gene_fits$con_mat[, cell_type]]) 
          {
          # next we loop through the different covariates in the matrix. Note that we start at 2 since the first columm is conventionally the intercept, which we do not want in our estimates matrix for meta analysis, as that adds computation for nothing. IF NEEDED, can add it later. 

          for (cov in 2:n_params) {

              #get the name of the covariate
              effect_id <- param_names[cov]

              #get the de estimate
              estimate <- gene_fits$all_vals[gene, cov, cell_type]
              
              #This is exactly how its acccessed in de_pop as well
              se_index <- (ct_index - 1) * n_params + cov

              se <- gene_fits$s_mat[gene, se_index]

              #add new entry in list
              converged <- isTRUE(gene_fits$con_mat[gene, cell_type]) &&
                  isTRUE(gene_fits$con_all[gene, se_index])

              if (!converged || is.na(estimate) || is.na(se) || !is.finite(estimate) || !is.finite(se) || se <= 0) {
                  next
              }

              rows_to_add[[row_i]] <- data.frame(
                  replicate_id = rep_id,   # passed down from function which wraps this (for reps)
                  gene = gene,
                  cell_type = cell_type,
                  effect_id = effect_id,
                  param_index = cov,
                  estimate = estimate,
                  se = se,
                  stringsAsFactors = FALSE
              )

              row_i <- row_i + 1
          }
      }
  }
  #in empty case
  if (length(rows_to_add) == 0) {
      return(data.frame(
          replicate_id = character(),
          gene = character(),
          cell_type = character(),
          effect_id = character(),
          param_index = integer(),
          estimate = numeric(),
          se = numeric(),
          stringsAsFactors = FALSE
      ))
      print("estimates matrix was empty")
  }

  estimates_df <- do.call(rbind, rows_to_add) #rbind all the rows (which are secretly dataframes) in our list - this is causing speed issues, need to optimize this in the second pass through
  rownames(estimates_df) <- NULL

  return(estimates_df)
}

#wrapper
extract_estimates_reps <- function(RCTD_reps) {

    if (length(RCTD_reps@RCTD.reps) == 0) {
    stop("RCTD_reps contains no replicate RCTD objects.")
    }

    estimates_df <- data.frame()
    #loop through replicates, and just keep merging them to estimate_df
    for (rep in seq_along(RCTD_reps@RCTD.reps)) {
        new_df <- extract_estimates(RCTD_reps@RCTD.reps[[rep]], names(RCTD_reps@RCTD.reps)[[rep]])
        estimates_df <- rbind(estimates_df, new_df)
    }
    rownames(estimates_df) <- NULL
    return(estimates_df)
}


#metadata extraction next
extract_metadata <- function(replicate_id, sample_id, subject_id, group_id, batch_id, moderator_df, rep_id_col) {

  #First step should be to validate rep_id's 
  if (anyNA(replicate_id)) {
      stop("replicate_ids cannot have any NA")
  }
  if (anyDuplicated(replicate_id)) {
      stop("replicate_ids cannot have any duplicated elements")
  }
  if (length(replicate_id) == 0 || is.null(replicate_id)) {
      stop("replicate_ids cannot be NULL or length 0")
  }
  if (any(replicate_id == "")) {
      stop("replicate_ids cannot be empty strings")
  }

  #can't think of any more right now, make sure to add more as see fit

  #initialize our dataframe with replicate_ids as the first column: 
  metadata_df <- data.frame(replicate_id = replicate_id, stringsAsFactors = FALSE)


  #validate_metadata_col is a helper function defined below. 
  metadata_df$sample_id <- validate_metadata_col(sample_id, replicate_id, "sample_id")
  metadata_df$subject_id <- validate_metadata_col(subject_id, replicate_id, "subject_id")
  metadata_df$group_id <- validate_metadata_col(group_id, replicate_id, "group_id")
  metadata_df$batch_id <- validate_metadata_col(batch_id, replicate_id, "batch_id")

  #Now we have to deal with moderators - for simplicity, assume design matrix is already properly formatted, since its a design matrix for moderators it should already be aligned to rep_ids
  moderator_df <- validate_moderators(moderator_df, replicate_id, rep_id_col = rep_id_col)

  if(!is.null(moderator_df)) {
      #add the moderator df on 
      metadata_df <- merge.data.frame(metadata_df, moderator_df, by = "replicate_id", all = FALSE, sort = FALSE)

      # make sure no rows were lost, this should be extra redundancy assuming validate_moderators did its job
      if (nrow(metadata_df) != length(replicate_id)) {
      stop("Merging moderators changed the number of metadata rows.")
      }
      # just for style, ensuring ordering is kept the same and we don't need any rownames now since our metadata is formatted properly.
      metadata_df <- metadata_df[match(replicate_id, metadata_df$replicate_id), , drop = FALSE]



  }
    rownames(metadata_df) <- NULL
    return(metadata_df)
}

#effects dictionary
make_effect_dictionary <- function(RCTD) {
    #extract the X2 matrix which contains cell-type specific covariates. This essentially contains all the data we wish to reorganize and document. 
    X2 <- RCTD@internal_vars_de$X2
    params_to_test <- RCTD@internal_vars_de$params_to_test

    effect_ids <- colnames(X2)
    #in case no colnames were designated in our design matrix to begin with, enforce a standardized naming routine
    if (is.null(effect_ids)) {
        effect_ids <- paste0("X2_param_", seq_len(ncol(X2)))
    }

    effect_dictionary <- data.frame(
        effect_id = effect_ids,
        param_index = seq_len(ncol(X2)),
        source_matrix = "X2",
        tested_by_cside = seq_len(ncol(X2)) %in% params_to_test,
        is_intercept = effect_ids == "intercept",
        stringsAsFactors = FALSE
    )

    return(effect_dictionary)
}

#then we wrap effects_dict
extract_effects_dictionary <- function(RCTD_reps) {
    if (length(RCTD_reps@RCTD.reps) == 0) {
    stop("RCTD_reps contains no replicate RCTD objects.")
    }
    # use the first rep to create the effect dictionary
    effects_dict <- make_effect_dictionary(RCTD_reps@RCTD.reps[[1]])
    #assume an invariant: all reps should have identical column names for their X2 matrices, as each replicate should be getting ran with an identical within-sample design matrix (e.g. explanatory_variable)
    for (rep in seq_along(RCTD_reps@RCTD.reps)) {
        new_dict <- make_effect_dictionary(RCTD_reps@RCTD.reps[[rep]])
        if (!identical(effects_dict$effect_id, new_dict$effect_id)) {
            stop("effect_id values do not match across RCTD replicates.")
        }
        if (!identical(effects_dict$param_index, new_dict$param_index)) {
            stop("param_index values do not match across RCTD replicates.")
        }
        if (!identical(effects_dict$tested_by_cside, new_dict$tested_by_cside)) {
            stop("tested_by_cside values do not match across RCTD replicates.")
        }
    }

    rownames(effects_dict) <- NULL
    return(effects_dict)
}