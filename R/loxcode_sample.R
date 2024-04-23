#' Load distance maps for distance-to-origin
#'
#' @param path path to directory containing maps named 0, 1, 2, 3, 4, corresponding to the size index of
#' 13, 9, 7, 5, 3 element cassettes respectively
#' @return distance map files
#' @export
#' @examples
#' # Example usage:
#' path <- "path/to/your/directory"
#' distance_maps <- load_distance_maps(path)
#' distance_maps
load_origin_distmaps <- function(path) {
    load_origin_files_wrapper(c(
        paste(path, 0, sep = '/'),
        paste(path, 1, sep = '/'),
        paste(path, 2, sep = '/'),
        paste(path, 3, sep = '/'),
        paste(path, 4, sep = '/')
    ))
}

#' Load distance maps for pairwise distances
#'
#' @param path path to directory containing the maps 0 (for 13-element distances) and 1 (for 9-element distances)
#' @return distance maps files for pairwise distance
#' @export
#' @examples
#' # Example usage:
#' path <- "path/to/your/directory"
#' pairwise_distance_maps <- load_distance_maps(path)
#' pairwise_distance_maps
load_pair_distmaps <- function(path) {
    load_pair_files_wrapper(c(paste(path, 0, sep = '/'),
                              paste(path, 1, sep = '/')))
}

#' Load probability tables (Markov chain)
#'
#' @param path path to directory containing the Markov chain tables
#' @return Probability table
#' @export
#' @examples
#' # Example usage:
#' path <- "path/to/markov_chain_tables_directory"
#' probability_table <- load_probability_tables(path)
#' probability_table
load_prob_files <- function(path) {
    size_idx <- 0:4
    rec_dist <- seq_len(15)
    l <-
        lapply(size_idx, function(s)
            sapply(rec_dist, function(r)
                paste(
                    path , paste0('size_', s, '_rec', r), sep = '/'
                )))
    load_prob_files_wrapper(l)
}

#' S4 class to contain output of decode()
#'
#' @slot data contains raw output of decode() as a data.frame
#' @slot read_ids list of FASTQ
#' @slot saturation list of colected files
#' @return decoded output
#' @export
#' @examples
#' decode_data <- data.frame(
#'   barcode = c("ATCG", "CGTA", "GATC"),
#'   count = c(10, 20, 15)
#' )
#' read_ids <- list(FASTQ1 = c("read1", "read2"), FASTQ2 = c("read3", "read4"))
#' saturation <- list("file1.txt", "file2.txt")
#'
#' # Create an instance of decode_output class
#' decode_output_instance <- new("decode_output",
#'                                data = decode_data,
#'                                read_ids = read_ids,
#'                                saturation = saturation)
#'
#' # Print the decode_output object
#' print(decode_output_instance)
decode_output <- setClass (
    "decode_output",

    # Defining slot type
    representation (
        data = "data.frame",
        saturation = "vector",
        read_ids = "list"
    ),

    # Initializing slots
    prototype = list(
        saturation = c(),
        data = data.frame(),
        read_ids = list()
    )
)

remove_existing <- function(x, n) {
    if (!any(is.na(match(n, names(x))))) {
        x <- x[,-match(n, names(x))]
    }
    return(x)
}

#' S4 class to represent loxcode experimental sample data
#'
#' @slot decode A data.frame to contain raw decode data from LoxCodeR2024::decode()
#' @slot meta A data.frame for user-defined sample metadata
#' @slot name Contains name of the experiment
#' @slot files A vector containing file names of loxcoder experiment
#' @slot decode_stats A list containing statistics of decoded data
#' @slot consensus_filtered_data A vector containing filtered data
#' @return Loxcode experiment sample data
#' @export
#' @examples
#' # Define sample metadata
#' sample_meta <- data.frame(
#'   sample_name = "Sample1",
#'   meta_value1 = 123,
#'   meta_value2 = "abc"
#' )
#'
#' # Create a LoxcodeExperimentSample object
#' sample_data <- LoxcodeExperimentSample(
#'   decode = data.frame(...),  # Raw decode data
#'   meta = sample_meta,        # Sample metadata
#'   name = "Sample1",          # Name of the experiment
#'   files = c("file1.fastq", "file2.fastq"),  # File names
#'   decode_stats = list(...),  # Decode statistics
#'   consensus_filtered_data = c(...)  # Filtered data
#' )
loxcode_sample <- setClass (
    "loxcode_sample",

    # Defining slot type
    representation (
        decode = "decode_output",
        name = "character",
        meta = "data.frame",
        files = "vector",
        decode_stats = "list",
        consensus_filtered_data = "vector"
    ),
    # Initializing slots
    prototype = list(
        decode = new("decode_output"),
        name = '',
        meta = data.frame(),
        files = c("", ""),
        decode_stats = list(),
        consensus_filtered_data = c("")
    )
)

#' Get number of rows in loxcode data
#'
#' @param x Loxcode object
#' @return length of rows in loxcode
#' @export
setMethod("length", "loxcode_sample", function(x)
    nrow(x@decode@data))

setMethod("nrow", "loxcode_sample", function(x)
    length(x))

#' Get loxcode sample name (description)
#'
#' @rdname name
#' @export
#' @examples
#' # Example usage:
#' # Assuming 'experiment' is a loxcode_sample object
#' # Retrieve the name of the experiment
#' experiment_name <- name(experiment)
#' experiment_name
setGeneric("name", function(x) {
    standardGeneric("name")
})

#' @rdname name
setMethod("name", "loxcode_sample", function(x)
    x@name)

# #' Set loxcode sample name (description)
# #'
# #' @param x sample to be named
# #' @param v name of sample
# #' @export
# setGeneric("name<-", function(x, v){standardGeneric("name<-")})
#
# setMethod("name<-", "loxcode_sample", function(x, v){
#   x@name <- v
# })

#' Add cassette validation column to decoded cassette data
#'
#' @param x loxcode cassette object
#' @return Cassette with validation column
#' @rdname validate
#' @export
#' @examples
#' # Example usage:
#' # Assuming 'cassette_data' is a loxcode cassette object
#' # Add a validation column to the decoded cassette data
#' validated_cassette <- validate(cassette_data)
#' validated_cassette

setGeneric("validate", function(x) {
    standardGeneric("validate")
})

#' @rdname validate
setMethod("validate", "loxcode_sample", function(x) {
    x@decode@data <- remove_existing(x@decode@data, 'is_valid')
    x@decode@data <-
        cbind(x@decode@data, data.frame(is_valid = is_valid(x@decode@data$code)))
    return(x)
})

#' Impute missing code in 13-element cassettes
#'
#' For 13-element cassettes that are missing a single element, the
#' missing element is imputed to minimise the resulting dist_orig.
#' @param x loxcode object
#' @return cassette with imputed missing code
#' @rdname impute
#' @export
#' @examples
#' # Example usage:
#' # Assuming 'lox_data' is a loxcode object containing 13-element cassettes
#' # Impute the missing code in each 13-element cassette
#' imputed_cassettes <- impute_missing_code_in_13_element_cassettes(lox_data)
#' imputed_cassettes

setGeneric("impute", function(x) {
    standardGeneric("impute")
})

#' @rdname impute
setMethod("impute", "loxcode_sample", function(x) {
    x@decode@data$code <-
        impute_13(x@decode@data$code, x@decode@data$size)
    return(x)
})

#' Get cassette IDs
#'
#' Appends a column of packed cassette IDs, or -1 if it cannot be packed
#' @param x loxcode data object
#' @return cassette Id
#' @rdname makeid
#' @export
#' @examples
#' # Example usage:
#' # Assuming 'lox_data' is a loxcode data object
#' # Get cassette IDs and append them as a new column
#' cassette_ids <- get_cassette_IDs(lox_data)
#' cassette_ids
setGeneric("makeid", function(x) {
    standardGeneric("makeid")
})


#' @rdname makeid
setMethod("makeid", "loxcode_sample", function(x) {
    x@decode@data <- remove_existing(x@decode@data, 'id')
    x@decode@data <-
        cbind(x@decode@data, data.frame(id = pack(
            x@decode@data$code, x@decode@data$is_valid
        )))
    return(x)
})

#' Fetch distances from origin (dist_orig)
#'
#' Appends a column of dist_orig values for valid cassettes
#' @param x loxcode object
#' @return Distance from origin of valid cassette
#' @rdname get_origin_dist
#' @export
#' @examples
#' # Load necessary libraries and data
#' library(LoxCodeR2024)
#' data <- read.csv("your_data.csv")
#'
#' # Create a loxcode object
#' lox_obj <- create_lox_object(data)
#'
#' # Fetch distances from origin
#' dist_orig_values <- get_origin_dist(lox_obj)
setGeneric("get_origin_dist", function(x) {
    standardGeneric("get_origin_dist")
})

#' @rdname get_origin_dist
setMethod("get_origin_dist", "loxcode_sample", function(x) {
    print("G")
    x@decode@data <- remove_existing(x@decode@data, 'dist_orig')
    print("H")
    print(x@files)
    x@decode@data <-
        cbind(x@decode@data,
              data.frame(dist_orig = retrieve_dist_origin(
                  x@decode@data$id, x@decode@data$size
              )))
    print("I")
    return(x)
})

#' Get recombination distance distribution
#'
#' Returns the proportion of valid cassettes for each valid recombination distance (dist_orig)
#' @param x Loxcode object
#' @param size size of the loxcode
#' @return distribution if distance of recombination
#' @rdname get_rec_prob
#' @export
#' @examples
#' # Load necessary libraries and data
#' library(LoxCodeR2024)
#' data <- read.csv("your_data.csv")
#'
#' # Assuming you have a Loxcode object called lox
#' size <- 10  # Set the size of the Loxcode
#' distribution <- get_rec_prob(lox, size)
setGeneric("get_rec_prob", function(x, size) {
    standardGeneric("get_rec_prob")
})

#' @rdname get_rec_prob
setMethod("get_rec_prob", "loxcode_sample", function(x, size) {
    r <- data.frame(table(valid(x)$dist_orig))
    names(r) <- c('rec', 'prob')
    r$prob <- r$prob / sum(r$prob)
    r$rec <- as.numeric(r$rec)
    return(r)
})

#' Get ensemble generation probability
#'
#' Retrieve ensemble probabilities as a weighted linear combination of Markov probabilities
#' where distances (dist_orig) are weighted using the sample distribution of recombination distances.
#'
#' Results are appended to readout data as a 'prob' column.
#'
#' @param x loxcode_sample object
#' @return ensemble genration probability
#' @rdname retrieve_prob_ensemble
#' @export
#' @examples
#' # Example usage:
#' # Assuming 'sample' is a loxcode_sample object
#' # Retrieve ensemble generation probability for the sample
#' ensemble_probability <- get_ensemble_generation_probability(sample)
#' ensemble_probability
setGeneric("retrieve_prob_ensemble", function(x) {
    standardGeneric("retrieve_prob_ensemble")
})

#' @rdname retrieve_prob_ensemble
setMethod("retrieve_prob_ensemble", "loxcode_sample", function(x) {
    sizes <- unique(LoxCodeR2024::valid(x)$size)
    x@decode@data$prob <- NA
    for (i in sizes) {
        # stratify by size
        r <- get_rec_prob(x, i)
        mask <-
            (x@decode@data$is_valid == TRUE & x@decode@data$size == i)
        probs <- rep(0, sum(mask))
        for (j in seq_len(nrow(r))) {
            print(paste('Weight: ', r$prob[j]))
            probs <-
                probs + r$prob[j] * LoxCodeR2024::retrieve_prob(x@decode@data[mask,]$id,
                                                                x@decode@data[mask,]$size,
                                                                rep(r$rec[j], sum(mask)))
        }
        x@decode@data[mask,]$prob <- probs
    }
    return(x)
})

#' Access decoded cassette data, valid cassettes only
#'
#' @param x loxcode_sample object
#' @return valid cassettes
#' @rdname valid
#' @export
#' @examples
#' # Example usage:
#' # Assuming 'sample' is a loxcode_sample object
#' # Access the decoded cassette data for valid cassettes only
#' valid_cassettes <- valid(sample)
#' valid_cassettes
setGeneric("valid", function(x) {
    standardGeneric("valid")
})

#' @rdname valid
setMethod("valid", "loxcode_sample", function(x) {
    v = x@decode@data

    return(v[v$is_valid == TRUE,])
})

#' #' Access decoded cassette data
#' #'
#' #' @param x loxcode_sample object
#' #' @rdname data
#' #' @export
#' setGeneric("data", function(x){ standardGeneric("data") })
#'
#' #' @rdname data
#' setMethod("data", "loxcode_sample", function(x){
#'   return(x@decode@data)
#' })
