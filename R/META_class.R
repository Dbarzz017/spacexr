#metaset class definition
setClass("metaset", slots = c(
    estimates = "data.frame",
    metadata = "data.frame",
    effects_dict = "data.frame",
    results = "list"
))


#constructor 
create.metaset <- function(estimates_df, metadata_df, effects_dict) {
    new(
        "metaset",
        estimates = estimates_df,
        metadata = metadata_df,
        effects_dict = effects_dict,
        results = list()
    )
}