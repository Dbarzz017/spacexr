#pair functions for all modes of main function: 
#fixed effects
pair_rma_calc_fe <- function(pair, min_reps = 2L, se_max = 4) {
    # good drafted result schematic. fit will be the index at which you obtain an entire rma object. This is heavy, but good for now. 
    out <- list(
        gene = pair$gene[[1]],
        cell_type = pair$cell_type[[1]],
        effect_id = pair$effect_id[[1]],
        n_reps = length(unique(pair$replicate_id)), 
        formula = ~1,
        random = NULL,
        fit = NULL,
        converged = FALSE,
        message = NA_character_
    )

    res <- tryCatch(
        {
            validate_pair(
                pair,
                min_reps = min_reps,
                se_max = se_max
            )

            metafor::rma.uni(
                yi = pair$estimate,
                sei = pair$se,
                method = "FE",
                data = pair
            )
        },
        error = function(e) e,
        warning = function(w) w
    )

    if (inherits(res, "error") || inherits(res, "warning")) {
        out$message <- conditionMessage(res)
        return(out)
    }

    out$fit <- res
    out$converged <- TRUE
    return(out)
}

#random effects
pair_rma_calc_re <- function(pair, min_reps = 2L, se_max = 4, test = 'z') {
    # good drafted result schematic. fit will be the index at which you obtain an entire rma object. This is heavy, but good for now. 
    out <- list(
        gene = pair$gene[[1]],
        cell_type = pair$cell_type[[1]],
        effect_id = pair$effect_id[[1]],
        n_reps = length(unique(pair$replicate_id)), 
        formula = ~1,
        random = NULL,
        fit = NULL,
        converged = FALSE,
        message = NA_character_
    )

    res <- tryCatch(
        {
            validate_pair(
                pair,
                min_reps = min_reps,
                se_max = se_max
            )

            metafor::rma.uni(
                yi = pair$estimate,
                sei = pair$se,
                method = "REML",
                data = pair,
                test = test
            )
        },
        error = function(e) e,
        warning = function(w) w
    )

    if (inherits(res, "error") || inherits(res, "warning")) {
        out$message <- conditionMessage(res)
        return(out)
    }

    out$fit <- res
    out$converged <- TRUE
    return(out)
}

#mixed effects
pair_rma_calc_me <- function(pair, formula = ~1, min_reps = 2L, se_max = 4, test = 'z') {
    # good drafted result schematic. fit will be the index at which you obtain an entire rma object. This is heavy, but good for now. 
    out <- list(
        gene = pair$gene[[1]],
        cell_type = pair$cell_type[[1]],
        effect_id = pair$effect_id[[1]],
        n_reps = length(unique(pair$replicate_id)), 
        formula = formula,
        random = NULL,
        fit = NULL,
        converged = FALSE,
        message = NA_character_
    )

    res <- tryCatch(
        {
            validate_pair(
                pair,
                formula = formula,
                min_reps = min_reps,
                se_max = se_max
            )

            metafor::rma.uni(
                yi = pair$estimate,
                sei = pair$se,
                mods = formula, 
                method = "REML",
                data = pair,
                test = test
            )
        },
        error = function(e) e,
        warning = function(w) w
    )

    if (inherits(res, "error") || inherits(res, "warning")) {
        out$message <- conditionMessage(res)
        return(out)
    }

    out$fit <- res
    out$converged <- TRUE
    return(out)
}


#multilevel
pair_rma_calc_mv <- function(pair, formula = ~1, random, min_reps = 2L, se_max = 4, test = 'z', dfs = 'residual') {
    #check random is nonempty
    if (missing(random) || is.null(random)) {
        stop("random argument must be nonempty and not null. If no random structure, then use rma.uni mode")
    }
    # good drafted result schematic. fit will be the index at which you obtain an entire rma object. This is heavy, but good for now. 
    out <- list(
        gene = pair$gene[[1]],
        cell_type = pair$cell_type[[1]],
        effect_id = pair$effect_id[[1]],
        n_reps = length(unique(pair$replicate_id)), 
        formula = formula,
        random = random,
        fit = NULL,
        converged = FALSE,
        message = NA_character_
    )

    res <- tryCatch(
        {
            validate_pair(
                pair,
                formula = formula,
                random = random,
                min_reps = min_reps,
                se_max = se_max
            )

            metafor::rma.mv(
                yi = pair$estimate,
                V = pair$se^2,
                mods = formula,
                random = random,
                data = pair,
                test = test, 
                dfs = dfs
            )
        },
            error = function(e) e,
            warning = function(w) w
    )

    if (inherits(res, "error") || inherits(res, "warning")) {
        out$message <- conditionMessage(res)
        return(out)
    }

    out$fit <- res
    out$converged <- TRUE
    return(out)
}

