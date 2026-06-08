#metadata column validators
validate_metadata_col <- function(x, rep_ids, col_name) {
    #idea is that this function makes sure everything with the column makes sense. We need the following to hold: 
    #length of x must be either length(rep_ids) or 1 (then we just assume its the same for all reps)
    #if names(x) = rep_ids, we must align from that.Otherwise, assume pairwise or replicated. 

    if (length(x) == 1) {
        return(rep(x, length(rep_ids)))
    } else if (length(x) == length(rep_ids)) {
        #check for names
        if (!is.null(names(x))) {

            if (anyDuplicated(names(x))) {
                stop(col_name, " has duplicated names.")
            }
            missing_ids <- setdiff(rep_ids, names(x))

            if(length(missing_ids) > 0) {
                stop(
                    col_name,
                    " is missing values for replicate_id: ",
                    paste(missing_ids, collapse = ", ")
                )
            }
            

            # otherwise, we have a valid list of names
            return(unname(x[rep_ids]))

        }
        # no list of names, match pairwise
        return(x)

    }
    # x doesn't have a permitted length
    stop(col_name,"must have length 1, or same length as rep_ids")
}

#also need to validate the moderators dataframe (design matrix) 
#user must input moderators to be either a dataframe, or NULL

validate_moderators <- function(moderators, rep_ids, rep_id_col = NULL) {
  if (is.null(moderators)) {
    return(NULL)
  } else if (!is.data.frame(moderators)) {
    stop("moderators must be NULL or a data.frame.")
  }

  moderators <- as.data.frame(moderators, stringsAsFactors = FALSE)
  reserved_cols <- c(
    "replicate_id",
    "sample_id",
    "subject_id",
    "group_id",
    "batch_id"
    )
    reserved_overlap <- intersect(
    setdiff(colnames(moderators), c("replicate_id", rep_id_col)),
    reserved_cols
    )


    if (length(reserved_overlap) > 0) {
        stop(
            "Moderator columns cannot use reserved names: ",
            paste(reserved_overlap, collapse = ", ")
        )
    }


  #first check to see if rep_ids is a column 
  if (!is.null(rep_id_col)) {
    if (!rep_id_col %in% colnames(moderators)) {
        #then the specified column, which is not null, is invalid
        stop(rep_id_col, " is not a valid colname")
    }

    ids <- as.character(moderators[[rep_id_col]])
    if (anyNA(ids)) {
    stop("moderators[[rep_id_col]] cannot have any NA")
    }
    if (anyDuplicated(ids)) {
        stop("moderators[[rep_id_col]] cannot have any duplicated elements")
    }

    missing_ids <- setdiff(rep_ids, ids)
    extra_ids <- setdiff(ids, rep_ids)

    if(length(missing_ids) > 0 || length(extra_ids) > 0) {
        stop(
            "moderators replicate IDs do not match rep_ids. Missing: ",
            paste(missing_ids, collapse = ", "),
            ". Extra: ",
            paste(extra_ids, collapse = ", ")
        )
    }


    #our column is valid so now we match the order of rows in moderators to our list of replicates so that when we join them later they are aligned. 

    moderators <- moderators[match(rep_ids, ids), , drop = FALSE]
    #make sure there is no preexisting column named replicate_id
    if (rep_id_col != "replicate_id" && "replicate_id" %in% colnames(moderators)) {
    stop(
        "moderators already has a replicate_id column, but rep_id_col = ",
        rep_id_col,
        ". Remove one of these columns or set rep_id_col = 'replicate_id'."
    )
    }
    colnames(moderators)[colnames(moderators) == rep_id_col] <- "replicate_id"
    #make sure we put replicate_id as our first column
    moderators <- moderators[, c("replicate_id", setdiff(colnames(moderators), "replicate_id")), drop = FALSE]
    rownames(moderators) <- NULL
    return(moderators)


  }
  #pretty much exactly the same thing but with the rownames since no colname was specified
  rn <- rownames(moderators)

  if (
    !is.null(rn) &&
      !anyNA(rn) &&
      !any(rn == "") &&
      !anyDuplicated(rn) &&
      setequal(rn, rep_ids)
  ) {
    moderators <- moderators[rep_ids, , drop = FALSE]
    #make sure we aren't overwriting an existing replicate_id column
    if ("replicate_id" %in% colnames(moderators)) {
      stop("moderators already has a replicate_id column. Either remove it or pass rep_id_col.")
    }

    moderators$replicate_id <- rep_ids

    moderators <- moderators[ , c("replicate_id", setdiff(colnames(moderators), "replicate_id")), drop = FALSE]

    rownames(moderators) <- NULL
    return(moderators)
  }
  #neither were valid, so raise error
  stop(
    "Could not align moderators to rep_ids. ",
    "Provide rep_id_col or set rownames(moderators) to rep_ids."
  )

}