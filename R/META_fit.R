#main fitting function
fit_metaset <- function(metaset, mode = "uni", submode = "random", formula = ~1, random = NULL, analysis_id = NULL, overwrite = FALSE, min_reps = 2L, se_max = 4, ct_prop_min = NULL) {
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
            res_list <- lapply(pairs_df, pair_rma_calc_re, min_reps = min_reps, se_max = se_max)
        } else if (submode == "mixed") {
            # logic for mixed-effects meta regression
            # need to isolate design matrix.  To do so, assume FOR NOW that our metadata_df is hardcoded and that we know that our moderator data starts at column 12 until the end in pairs_df. This should probably be fixed because it seems very hardcoded. 
            res_list <- lapply(pairs_df, pair_rma_calc_me, formula = active_formula, min_reps = min_reps, se_max = se_max)

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

        res_list <- lapply(pairs_df, pair_rma_calc_mv, formula = active_formula, random = random, min_reps = min_reps, se_max = se_max)

        #for now, assume user passes through their nesting structure for the variance-covariance matrix with a value in random 

    } else {
        stop("mode must be either 'uni' or 'mv'.")
    }

    #clean up results, maybe rbind or something similar, organize everything, adjust p-vals, etc... not sure what needs to be done here but this is the last step so it'll be done last. 

    # check to see if analysis_id is null, at which point use our automatic formula: 
    if (is.null(analysis_id)) {
    formula_id <- paste(deparse(formula), collapse = "")
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

    #lastly, return a list at entry "analysis_id" in teh results slot of our metaset. 
    metaset@results[[analysis_id]] <- list(
    analysis_id = analysis_id,
    call = match.call(),
    mode = mode,
    submode = submode,
    formula = formula,
    random = random,
    n_pairs = length(res_list),
    fit_records = res_list
    )

    return(metaset) 
}