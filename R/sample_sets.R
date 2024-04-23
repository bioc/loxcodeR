# library(plyr)
#' Summary table (v3)
#'
#' @param lox loxcode_experiment object
#' @param sample_set sample set
#' @return table of sample information
#' @rdname summary_table
#' @export
#' @examples
#' # Example usage:
#' # Assuming 'experiment' is a loxcode_experiment object and 'sample_set' is the sample set to summarize
#' # Generate a summary table for the specified sample set in the experiment
#' summary_table <- summary_table(experiment, sample_set)
#' summary_table
setGeneric("summary_table", function(lox, sample_set = "all_samples") {
    standardGeneric("summary_table")
})

#' @rdname summary_table
setMethod("summary_table", "loxcode_experiment", function(lox, sample_set =
                                                              "all_samples") {
    sample_name = NULL
    # initialize variables
    samples = lox@samples
    counts = lox@count_matrixes[[sample_set]]
    aliases = lox@alias[[sample_set]]
    all_codes = lox@code_sets[["all_codes"]]
    invalid_codes = lox@code_sets[["invalid_codes"]]
    meta = lox@meta

    # subset data
    samples = samples[names(samples) %in% names(counts)]
    metadata = subset(meta, sample_name %in% names(counts))

    # add data to data table
    table = data.frame(matrix(ncol = 2, nrow = length(samples)))
    names(table) = c("sample_name", "Alias")
    table$sample_name = names(counts)
    table$Alias = aliases$alias
    table = merge(
        table,
        metadata,
        all = TRUE,
        by = c("sample_name"),
        sort = FALSE
    )
    table$`Barcode Count` = sapply(samples, function(x)
        nrow(x@decode@data))
    table$`Invalid Count` = sapply(samples, function(x)
        nrow(subset(x@decode@data, is_valid == FALSE)))
    table$`Number of Reads` = sapply(samples, function(x)
        x@decode_stats$tot_reads)
    table$`Max Complexity` = sapply(samples, function(x)
        max(x@decode@data$dist_orig, na.rm = TRUE))
    table$`Consensus Filtered` = sapply(samples, function(x)
        x@decode_stats$consensus_filtered)
    table$`Percent Filtered` = sapply(samples, function(x)
        round(
            100 * x@decode_stats$consensus_filtered / x@decode_stats$tot_reads,
            2
        ))
    return (table)
})

#' fill set aliases
#' helper function complete the aliases for the specified sample set
#'
#' @param lox loxcode_experiment object
#' @param count_matrix count matrix of loxcode object
#' @return loxcode object with filled alias
fillSetAliases <- function(lox, count_matrix) {
    if (is.null(lox) || is.null(count_matrix) ||
        !count_matrix %in% names(lox@count_matrixes)) {
        return ()
    }

    aliases = data.frame(
        sample_name = names(lox@count_matrixes[[count_matrix]]),
        alias = "",
        stringsAsFactors = FALSE
    )
    aliases$alias = c(sapply(seq_len(nrow(aliases)), function(x)
        paste("Sample", x)))

    lox@alias[[count_matrix]] = aliases
    return (lox)
}

#' rename sample set
#'
#' @param x loxcode_experiment object
#' @param s sample set
#' @param n new code set name
#' @return new loxcode_experiment object
#' @rdname rename_sampleset
#' @export
#' @examples
#' # Example usage:
#' # Assuming 'experiment' is a loxcode_experiment object and 'old_set_name' is the old sample set name to be renamed
#' # Rename the sample set named 'old_set_name' to 'new_set_name'
#' new_experiment <- rename_sample_set(x = experiment, s = "old_set_name", n = "new_set_name")
#' new_experiment
setGeneric("rename_sampleset", function(x, s, n) {
    standardGeneric("rename_sampleset")
})

#' @rdname rename_sampleset
setMethod("rename_sampleset", "loxcode_experiment", function(x, s, n) {
    if (s %in% c("all_samples")) {
        return(x)
    }

    temp = x@count_matrixes[[s]]
    x@count_matrixes[[s]] <- NULL
    x@count_matrixes[[n]] = temp
    temp = x@alias[[s]]
    x@alias[[s]] <- NULL
    x@alias[[n]] = temp
    return(x)
})

#' Get metadata of a merged experiment
#'
#' @param x loxcode_experiment object
#' @param s set of merged samples
#' @return a data frame of meta data
#' @rdname get_collapsed_meta
#' @export
#' @examples
#' # Load necessary libraries and data
#' library(LoxCodeR2024)
#' data <- read.csv("your_data.csv")
#'
#' # Create a loxcode experiment object
#' lox_experiment <- create_lox_experiment(data)
#'
#' # Merge samples
#' merged_samples <- merge_samples(lox_experiment, sample_set = "Sample_Set")
#'
#' # Get collapsed metadata
#' meta <- get_collapsed_meta(lox_experiment, merged_samples)

setGeneric("get_collapsed_meta", function(x, s) {
    standardGeneric("get_collapsed_meta")
})

#' @rdname get_collapsed_meta
setMethod("get_collapsed_meta", "loxcode_experiment", function(x, s) {
    sample_name = NULL
    counts = x@count_matrixes[[s]]
    sample_names = names(counts)

    # Divide samples into collapsed vs non-collapsed
    non_col_names = intersect(names(x@samples), names(counts))
    col_names = setdiff(sample_names, non_col_names)

    # metadata of collapsed samples
    meta = sapply(x@meta, unique)
    df1 = data.frame(matrix(ncol = length(meta), nrow = length(col_names)))
    names(df1) = names(meta)
    df1$sample_name = col_names
    if (length(col_names)) {
        for (i in seq_along(col_names)) {
            components = unlist(strsplit(col_names[i], split = "__"))
            if (length(components) == 1 &
                grepl("NA", components[[1]])) {
                # no metadata to insert
                break
            }
            else {
                used_cols = c()
                for (metadata in components) {
                    possible_cols = names(meta)[sapply(meta, function(x)
                        metadata %in% x)]
                    possible_cols = setdiff(possible_cols,
                                            intersect(possible_cols, used_cols))
                    df1[i, possible_cols[[1]]] = metadata
                    used_cols = c(used_cols, possible_cols[1])
                }
            }
        }
    }

    # metadata of non_collapsed samples
    df2 = subset(x@meta, sample_name %in% non_col_names)


    df = rbind.fill(df2, df1)
    df = df[!is.na(df$sample_name),]
    row.names(df) = df$sample_name
    df$sample_name = NULL

    return(df)
})

#' Collapse selected samples (v2)
#'
#' @param lox loxcode_experiment object
#' @param count_matrix sample set/ count_matrix
#' @param code_set code set
#' @param index indices of samples to collapse
#' @param name new sample name
#' @param union boolean: use union or intersection of samples
#' @param average boolean: use average of sum of counts
#' @return updated loxcode_experiment object
#' @rdname collapse_selection
#' @export
#' @examples
#' # Load required packages
#' library(LoxCodeR2024)
#'
#' # Example usage
#' # Assuming lox, count_matrix, code_set, index, name, union, and average are defined and have the required structure
#' lox <- readRDS("~/Desktop/LoxCodeR2024/LoxcodeR_app/Week2.rds")
#' collapse_selection(lox, count_matrix = "all_samples", code_set = "all_codes", index = c(1, 3, 5), name = NULL, union = TRUE, average = FALSE)

setGeneric("collapse_selection", function(lox,
                                          count_matrix = "all_samples",
                                          code_set = "all_codes",
                                          index,
                                          name = NULL,
                                          union = TRUE,
                                          average = FALSE) {
    standardGeneric("collapse_selection")
})

#' @rdname collapse_selection
setMethod("collapse_selection", "loxcode_experiment", function(lox,
                                                               count_matrix = "all_samples",
                                                               code_set = "all_codes",
                                                               index,
                                                               name = NULL,
                                                               union = TRUE,
                                                               average = FALSE) {
    counts = lox@count_matrixes[[count_matrix]]
    samples_list = names(counts)[index]

    # make new loxcode_sample
    new_sample = merge_sample_list(lox, samples_list, union, average)
    name = switch(is.null(name), paste0(samples_list, collapse = "_"), name)
    new_sample@name = name
    lox@samples[[name]] = new_sample

    # add to metadata
    metadata = cbind(data.frame(sample_name = name), new_sample@meta)
    lox@meta = rbind.fill(lox@meta, metadata)

    # add to current count_matrix
    if (length(index) == 1) {
        counts[[name]] = counts[[samples_list]]
    } else if (union & !average) {
        counts[[name]] = rowSums(counts[, samples_list])
    } else if (union & average) {
        counts[[name]] = rowSums(counts[, samples_list]) / length(samples_list)
    } else if (!union & !average) {
        counts[[name]] = (rowSums(counts[, samples_list]) * matrixStats::rowProds(as.matrix(counts)[, samples_list] >
                                                                                      0))
    } else if (!union & average) {
        counts[[name]] = rowSums(counts[, samples_list]) * matrixStats::rowProds(as.matrix(counts)[, samples_list] >
                                                                                     0) / length(samples_list)
    }
    lox@count_matrixes[[count_matrix]] = counts

    # add to aliases
    aliases = lox@alias[[count_matrix]]
    row = data.frame("sample_name" = name,
                     "alias" = paste("Sample", nrow(aliases) + 1))
    lox@alias[[count_matrix]] = rbind(aliases, row)

    return (lox)
})

#' Collapses the count matrix into samples with the sample metadata (v2)
#'
#' @param lox loxcode_experiment object
#' @param count_matrix current count_matrix
#' @param collapse column names of metadata on which to collapse
#' @param name name of new count_matrix
#' @param union boolean, True if barcodes in either should be counted
#' @param average boolean, True if counts should be averaged instead of summed
#' @return loxcode_experiment object with new collapsed samples
#' @rdname collapse
#' @export
#' @examples
#' # Load required packages
#' library(LoxCodeR2024)
#'
#' # Example usage
#' # Assuming lox, count_matrix, collapse, name, union, and average are defined and have the required structure
#' lox <- readRDS("~/Desktop/LoxCodeR2024/LoxcodeR_app/Week2.rds")
#' count_matrix="all_samples"
#' name="new"
#' collapse="Organ"
#' collapse(lox, count_matrix, collapse, name = "new_count_matrix", union = TRUE, average = FALSE)
setGeneric("collapse", function(lox,
                                count_matrix,
                                collapse,
                                name,
                                union = TRUE,
                                average = FALSE) {
    standardGeneric("collapse")
})

#' @rdname collapse
setMethod("collapse", "loxcode_experiment", function(lox,
                                                     count_matrix,
                                                     collapse,
                                                     name,
                                                     union = TRUE,
                                                     average = FALSE) {
    counts = lox@count_matrixes[[count_matrix]]
    meta = subset(lox@meta, sample_name %in% names(counts))
    params = match(collapse, names(meta))

    # group samples based on which share metadata
    track = list()
    for (i in seq_len(nrow(meta))) {
        sample = meta$sample_name[i]
        metadata = paste0(meta[i, params], collapse = "_")
        if (metadata %in% names(track))
            track[[metadata]] = c(track[[metadata]], sample)
        else
            track[[metadata]] = c(sample)
    }
    # remove metadata with NA
    track[grepl("NA", names(track))] <- NULL
    # early return if no samples to merge
    if (length(track) == 0) {
        return (lox)
    }

    # merge samples that share the same meta data
    new_samples = lapply(track, function(samples)
        merge_sample_list(lox, samples, union, average))
    for (i in seq_along(new_samples)) {
        sample_name = names(track)[i]
        new_samples[[i]]@name = sample_name
        metadata = cbind(data.frame(sample_name = sample_name),
                         new_samples[[i]]@meta)
        lox@meta = rbind.fill(lox@meta, metadata)
    }
    lox@samples = c(lox@samples, new_samples)

    # modify count_matrixes
    codes = lox@code_sets[["all_codes"]]$code
    new_count_matrix = data.frame(matrix(nrow = length(codes), ncol = length(new_samples)))
    row.names(new_count_matrix) = codes
    names(new_count_matrix) = names(track)
    for (i in seq_along(track)) {
        if (length(track[[i]]) == 1)
            new_count_matrix[[i]] = counts[[track[[i]][1]]]
        else if (union &
                 !average)
            new_count_matrix[[i]] = rowSums(counts[, track[[i]]])
        else if (union &
                 average)
            new_count_matrix[[i]] = rowSums(counts[, track[[i]]]) / length(track[[i]])
        else if (!union &
                 !average)
            new_count_matrix[[i]] = (rowSums(counts[, track[[i]]]) * matrixStats::rowProds(as.matrix(counts)[, track[[i]]] >
                                                                                               0))
        else if (!union &
                 average)
            new_count_matrix[[i]] = rowSums(counts[, track[[i]]]) * matrixStats::rowProds(as.matrix(counts)[, track[[i]]] >
                                                                                              0) / length(track[[i]])
    }
    lox@count_matrixes[[name]] = new_count_matrix

    # fill aliases
    lox = fillSetAliases(lox, name)

    return(lox)
})

#' Merge list of samples
#' Merges a list of samples into one loxcode_sample object
#'
#' @param lox loxcode_experiment object
#' @param samples list of sample names to be merged
#' @param union True if barcodes in either samples should be counted
#' @param average True if counts should be averaged instead of summed
#' @return merged loxcode_sample object
#' @rdname merge_sample_list
#' @export
#' @examples
#' # Example usage:
#' # Assuming 'experiment' is a loxcode_experiment object and 'sample_list' is a list of sample names to be merged
#' # Merge the samples in 'sample_list' into a single loxcode_sample object, summing their counts
#' merged_sample <- merge_sample_list(lox = experiment, samples = sample_list, union = TRUE, average = FALSE)
#' merged_sample
setGeneric("merge_sample_list", function(lox,
                                         samples,
                                         union = TRUE,
                                         average = FALSE) {
    standardGeneric("merge_sample_list")
})

#' @rdname merge_sample_list
setMethod("merge_sample_list", "loxcode_experiment", function(lox,
                                                              samples,
                                                              union = TRUE,
                                                              average = FALSE) {
    # initialize and declare variables
    new = new("loxcode_sample")
    FIXED_PARAMS = c("code", "size", "is_valid", "id", "dist_orig")
    samples_list = lox@samples[samples]

    # fill decode@data slot
    if (length(samples_list) == 1) {
        merged_data = samples_list[[1]]@decode@data
        counts = merged_data$count
        firstreads = merged_data$firstread
    }
    else {
        data_list = c()
        data_list = lapply(samples_list, function(x)
            c(data_list, x@decode@data))
        merged_data = Reduce(function(x, y)
            merge(x, y, all = union, by = FIXED_PARAMS),
            data_list)
        counts = rowSums(merged_data[, grepl("count", names(merged_data))], na.rm = TRUE)
        if (average) {
            counts = counts / length(samples)
        }
        firstreads = rowSums(merged_data[, grepl("firstread", names(merged_data))], na.rm = TRUE) / length(samples)
    }
    new@decode@data = data.frame(
        count = counts,
        firstread = firstreads,
        code = merged_data$code,
        size = merged_data$size,
        is_valid = merged_data$is_valid,
        id = merged_data$id,
        dist_orig = merged_data$dist_orig
    )

    new@decode@saturation = integer(0)
    new@decode@read_ids = list()

    # fill in the metadata slot
    meta_list = c()
    meta_list = lapply(samples_list, function(x)
        c(meta_list, x@meta))
    meta = Reduce(function(x, y)
        merge(x, y, all = TRUE), meta_list)
    meta[] <-
        data.frame(lapply(meta, as.character), stringsAsFactors = FALSE)
    if (length(samples_list) == 1)
        new@meta = data.frame(meta)
    else
        new@meta = meta[lapply(meta, function(x)
            length(unique(x))) == 1][1,]

    # files slot
    new@files = list()
    new@files = lapply(samples_list, function(x)
        list(new@files, x@files))

    # decode stats slot
    for (stat in names(samples_list[[1]]@decode_stats)) {
        values = lapply(samples_list, function(x)
            x@decode_stats[[stat]])
        new@decode_stats[[stat]] = sum(unlist(values))
    }

    return (new)
})

#' Creates a new count_matrix with dataframe of samples
#'
#' @param lox loxcode_experiment object
#' @param count_matrix name of existing count_matrix to choose from
#' @param indices indices to include
#' @param name name of the new count_matrix
#' @return new loxcode_experiment object with a new count_matrix
#' @rdname make_count_matrix
#' @export
#' @examples
#' # Example usage:
#' # Assuming 'experiment' is a loxcode_experiment object
#' # Create a new count_matrix named 'new_count_matrix' by selecting indices 1, 3, and 5 from an existing count_matrix 'old_count_matrix'
#' new_experiment <- make_count_matrix(experiment, count_matrix = "old_count_matrix", indices = c(1, 3, 5), name = "new_count_matrix")
#' new_experiment

setGeneric("make_count_matrix", function(lox,
                                         count_matrix = "all_samples",
                                         indices,
                                         name = "test") {
    standardGeneric("make_count_matrix")
})

#' @rdname make_count_matrix
setMethod("make_count_matrix", "loxcode_experiment", function(lox,
                                                              count_matrix = "all_samples",
                                                              indices,
                                                              name = "test") {
    # handle null values
    if (is.null(lox) || is.null(indices)) {
        return ()
    }


    indices = sort(indices)
    new_matrix = lox@count_matrixes[[count_matrix]][indices]
    lox@count_matrixes[[name]] = new_matrix
    lox = fillSetAliases(lox, name)
    return(lox)
})
