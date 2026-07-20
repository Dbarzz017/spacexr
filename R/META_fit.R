#main fitting function
fit_metaset <- function(metaset, mode = "uni", submode = "random", formula = ~1, random = NULL, analysis_id = NULL, overwrite = FALSE, min_reps = 2L, se_max = 4, ct_prop_min = NULL, test = 'z', dfs = 'residual') {
    # first we must do all sorts of preprocessing/preparation of the data with regards to universal things. That is, focusing on aspects of the function which must be applied to all data regardless of mode. 
    # So far, merging the estimates and metadata, isolating all unique pairs of ct/gene/effect, and filtering based off how many reps per pair (assume k reps per pair at minimum, k >= 3). Incorporate other filters later: 

    merged_df <- merge(metaset@estimates, metaset@metadata, by = "replicate_id", all = FALSE, sort = FALSE)

    # Normalize submode first.
    if (is.null(submode) || submode == "") {
        submode <- "random"
    }

    active_formula <- if (mode == "uni" && submode %in% c("fixed", "random")) {
        ~1
    } else {
        formula
    }

    # validate test and df parameters
    if (!test %in% c("z", "t", "knha", "adhoc")) {
        stop("test must be one of: 'z', 't', 'knha', or 'adhoc'.")
    }

    if (mode == "mv") {
        if (!test %in% c("z", "t")) {
            stop("For mode = 'mv', test must be either 'z' or 't'.")
        }

        if (!dfs %in% c("residual", "contain")) {
            stop("For mode = 'mv', dfs must be 'residual' or 'contain'.")
        }
    }

    # Validate formula variables before model.matrix().
    formula_vars <- all.vars(active_formula)
    missing_formula_vars <- setdiff(formula_vars, colnames(merged_df))

    if (length(missing_formula_vars) > 0) {
        stop(
            "Formula variables missing from merged data: ",
            paste(missing_formula_vars, collapse = ", ")
        )
    }

    # Apply row-level filters.
    merged_df <- merged_df[
        is.finite(merged_df$estimate) &
        is.finite(merged_df$se) &
        merged_df$se > 0 &
        merged_df$se < se_max,
        ,
        drop = FALSE
    ]

    if (!is.null(ct_prop_min)) {
        if (
            length(ct_prop_min) != 1L ||
            !is.numeric(ct_prop_min) ||
            !is.finite(ct_prop_min) ||
            ct_prop_min < 0 ||
            ct_prop_min > 1
        ) {
            stop("ct_prop_min must be NULL or one numeric value between 0 and 1.")
        }

        if (!"ct_prop" %in% colnames(merged_df)) {
            stop("ct_prop_min was provided, but estimates does not contain ct_prop.")
        }

        merged_df <- merged_df[
            is.finite(merged_df$ct_prop) &
            merged_df$ct_prop >= ct_prop_min,
            ,
            drop = FALSE
        ]
    }

    if (length(formula_vars) > 0) {
        merged_df <- merged_df[
            complete.cases(merged_df[, formula_vars, drop = FALSE]),
            ,
            drop = FALSE
        ]
    }

    if (nrow(merged_df) == 0) {
        stop("No estimate rows remain after preprocessing.")
    }

    # Global design check.
    X <- model.matrix(active_formula, data = merged_df)

    if (qr(X)$rank < ncol(X)) {
        stop("Design matrix is rank deficient.")
    }

    required_min_reps <- max(
        min_reps,
        qr(X)$rank + 1L
    )




    # Construct pairs after all row filtering.
    pairs <- interaction(
        merged_df$gene,
        merged_df$cell_type,
        merged_df$effect_id,
        drop = TRUE
    )

    pairs_df <- split(merged_df, pairs)

    pair_rep_counts <- vapply(
        pairs_df,
        function(pair) length(unique(pair$replicate_id)),
        integer(1)
    )

    pairs_df <- pairs_df[pair_rep_counts >= required_min_reps]

    if (length(pairs_df) == 0) {
        stop("No pairs have enough replicates after preprocessing.")
    }

    message("Number of pairs with sufficient reps: ", length(pairs_df))

    #include more preprocessing steps below this line: 

    #Now we work on the casework necessary: 
    #first up is a normal meta analysis (the easy part). 
    # we have two subcases for a normal meta analysis (so far): either you perform a fixed-effects meta analysis, or random effects. We include both, and they will be sub methods. If no sub method is chosen, assume random effects (for now) 
    if (mode == "uni") {
        #run whatever needs to be done for both fixed/random in the specifically rma.uni case 

        if (submode == "fixed") {
            #logic for fixed effect meta analysis goes here
            res_list <- lapply(pairs_df, pair_rma_calc_fe, min_reps = min_reps, se_max = se_max)

        } else if (submode == "random") {
            # logic for random-effects meta analysis goes here 
            res_list <- lapply(pairs_df, pair_rma_calc_re, min_reps = min_reps, se_max = se_max, test = test)
        } else if (submode == "mixed") {
            # logic for mixed-effects meta regression
            # need to isolate design matrix.  To do so, assume FOR NOW that our metadata_df is hardcoded and that we know that our moderator data starts at column 12 until the end in pairs_df. This should probably be fixed because it seems very hardcoded. 
            res_list <- lapply(pairs_df, pair_rma_calc_me, formula = active_formula, min_reps = min_reps, se_max = se_max, test = test)

        } else {
            # throw an error
            stop("when mode is uni, submode must be either 'random', 'fixed', 'mixed' or empty.")
        }
    # next is the meta regression case. The submodes here are still up for debate. Current idea is to have a couple of the heavy hitting  common experiemnt schema (such as samples -> subjects as a hierarchy) and then to include a custom mode which can be modular and designed by user. Thinking some sort of graph to represent data could be an interesting path to take this, however not sure yet/implementation details for later on. if no submode is elected, raise a flag and stop.
    } else if (mode == "mv") { 
        # run whatever needs to be done for all possible submodes of meta regression - likely checks for if parameters empty/nonempty

        if (is.null(random)) {
            stop("random must be supplied for mode = 'mv'. If no random structure, then use mode = 'uni'.")
        }

        res_list <- lapply(pairs_df, pair_rma_calc_mv, formula = active_formula, random = random, min_reps = min_reps, se_max = se_max, test = test, dfs = dfs)

        #for now, assume user passes through their nesting structure for the variance-covariance matrix with a value in random 

    } else {
        stop("mode must be either 'uni' or 'mv'.")
    }

    #clean up results, maybe rbind or something similar, organize everything, adjust p-vals, etc... not sure what needs to be done here but this is the last step so it'll be done last. 

    # check to see if analysis_id is null, at which point use our automatic formula: 
    if (is.null(analysis_id)) {
    formula_id <- paste(deparse(active_formula), collapse = "")
    formula_id <- gsub("[^A-Za-z0-9_]+", "_", formula_id)
    formula_id <- gsub("^_+|_+$", "", formula_id)

    random_id <- if (is.null(random)) {
        "no_random"
    } else {
        paste(deparse(random), collapse = "")
    }
    random_id <- gsub("[^A-Za-z0-9_]+", "_", random_id)
    random_id <- gsub("^_+|_+$", "", random_id)

    analysis_id <- paste(
        "meta",
        mode,
        submode,
        formula_id,
        random_id,
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        sep = "_"
    )
    }

    if (analysis_id %in% names(metaset@results) && !overwrite) {
        stop("analysis_id already exists in metaset@results. Set overwrite = TRUE to replace it.")
    }

    # summary of the analysis (not data, but the analysis in general) 
    analysis_info <- list(
        analysis_id = analysis_id,
        mode = mode,
        submode = submode,
        formula = paste(deparse(active_formula), collapse = ""),
        random = if (is.null(random)) NA_character_ else paste(deparse(random), collapse = ""),
        min_reps = min_reps,
        se_max = se_max,
        ct_prop_min = ct_prop_min,
        n_pairs = length(res_list),
        n_converged = sum(vapply(res_list, function(x) isTRUE(x$converged), logical(1))),
        n_failed = sum(!vapply(res_list, function(x) isTRUE(x$converged), logical(1))),
        call = match.call()
    )

    # summary of the main results ordered by significance
    results_summary <- do.call(rbind, lapply(res_list, function(rec) {
        # in case pair didn't converge or is null
        if (!isTRUE(rec$converged) || is.null(rec$fit)) {
            return(data.frame(
                gene = rec$gene,
                cell_type = rec$cell_type,
                effect_id = rec$effect_id,
                term = NA_character_,
                estimate = NA_real_,
                se = NA_real_,
                zval = NA_real_,
                p_val = NA_real_,
                ci_low = NA_real_,
                ci_high = NA_real_,
                n_reps = rec$n_reps,
                tau2 = NA_real_,
                I2 = NA_real_,
                QE = NA_real_,
                QEp = NA_real_,
                converged = FALSE,
                message = rec$message,
                stringsAsFactors = FALSE
            ))
        }
        # otherwise populate: this is all simple parsing
        fit <- rec$fit
        terms <- rownames(fit$b)
        if (is.null(terms)) {
            terms <- paste0("term_", seq_along(as.numeric(fit$b)))
        }
        
        data.frame(
            gene = rec$gene,
            cell_type = rec$cell_type,
            effect_id = rec$effect_id,
            term = terms,
            estimate = as.numeric(fit$b),
            se = as.numeric(fit$se),
            zval = as.numeric(fit$zval),
            p_val = as.numeric(fit$pval),
            ci_low = as.numeric(fit$ci.lb),
            ci_high = as.numeric(fit$ci.ub),
            n_reps = rec$n_reps,
            tau2 = if (!is.null(fit$tau2)) fit$tau2 else NA_real_,
            I2 = if (!is.null(fit$I2)) fit$I2 else NA_real_,
            QE = if (!is.null(fit$QE)) fit$QE else NA_real_,
            QEp = if (!is.null(fit$QEp)) fit$QEp else NA_real_,
            converged = TRUE,
            message = NA_character_,
            stringsAsFactors = FALSE
        )
    }))

    # make qvals and sort by significance
    rownames(results_summary) <- NULL
    results_summary$q_val <- NA_real_
    finite <- is.finite(results_summary$p_val)
    if (any(finite)) {
        q_groups <- interaction(
            results_summary$cell_type,
            results_summary$effect_id,
            results_summary$term,
            drop = TRUE
        )
        results_summary$q_val[finite] <- unsplit(
            lapply(
                split(results_summary$p_val[finite], q_groups[finite]),
                p.adjust,
                method = "BH"
            ),
            q_groups[finite]
        )
    }
    q_order <- ifelse(is.na(results_summary$q_val), Inf, results_summary$q_val)
    results_summary <- results_summary[
        order(q_order, -abs(results_summary$estimate)),
        ,
        drop = FALSE
    ]
    rownames(results_summary) <- NULL
    #lastly, return a list at entry "analysis_id" in teh results slot of our metaset. 
    metaset@results[[analysis_id]] <- list(
        analysis_info = analysis_info,
        results_summary = results_summary,
        fit_records = res_list
    )

    return(metaset) 
}