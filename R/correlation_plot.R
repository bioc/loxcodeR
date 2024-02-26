#' correlation plot
#'
#' Function used to create a correlation plot
#'
#' @param lox loxcode_experiment object
#' @param count_matrix A matrix of the selected samples
#' @param code_set A matrix of the selected code sets
#' @param parameter metadata type to display
#' @param data metadata to display
#' @param method_ correlation method (p)
#' @export

setGeneric("correlation_plot", function(lox,
                                        count_matrix = "all_samples",
                                        code_set = "all_codes",
                                        # parameter,
                                        # data,
                                        method_ = "pearson") {
  standardGeneric("correlation_plot")
})

setMethod("correlation_plot", "loxcode_experiment",
          function(lox,
                   count_matrix = "all_samples",
                   code_set = "all_codes",
                   # parameter,
                   # data,
                   method_ = "pearson") {

            if (is.null(lox) || is.null(count_matrix) || is.null(method_) ||
                is.null(code_set))
           # || is.null(parameter) || is.null(data) || data=="")
                {

              return ()
            }

            # initialize data
            counts = lox@count_matrixes[[count_matrix]]
            codes = lox@code_sets[[code_set]]
            # samp_table = lox@meta
            # params = getMetaParameters(lox)
            #
            # # select required meta data
            # primary_data = samp_table[[parameter]][match(names(counts), samp_table$sample_name)]
            # other_parameters = params[!params %in% c(parameter)]
            #
            # # update sample names in count_matrix
            # names = paste(primary_data)
            # for (param in other_parameters) {
            #   row = samp_table[[param]][match(names(counts), samp_table$sample_name)]
            #   names = paste(names, row)
            # }
            # names(counts) = names

            # subset from count_matrix relevant data
            counts[counts < 10] = 0
            counts = counts[rowSums(counts) > 0, ]
            #counts = counts[, primary_data %in% c(data)]
            counts = counts[order(row.names(counts)), order(names(counts))]

            # subset based on code complexity and flip distance
            #codes = codes[codes$dist_orig > 4 & codes$flp_dist > 1, ]
            counts = counts[row.names(counts) %in% codes$code, ]


            # shared = (sum(rowSums(counts[, 1:2]) > 0 &
                            # rowSums(counts[, 3:4]) > 0) / sum(rowSums(counts[, 1:2]) > 0 |
                            #                                     rowSums(counts[, 3:4]) > 0))

            # find correlations
            #counts = as.matrix(counts)
            COR = as.data.frame(cor(log(counts + 1), method= method_))
            if (nrow(COR) == 0) {
              return ()
            }

            # plot the data
            heatmap(
              as.matrix(COR),
              #Rowv = NA,
              #Colv = NA,
              # main = paste0(parameter, " #", data, " barcode shared P5/P6: ", 100 * (round(shared, 2)), "%"),
              #main = paste0(parameter, " #", data),
              scale = "none"
            )

})

# getMetaParameters <- function(lox) {
#   params = names(lox@meta)[!names(lox@meta) %in% c("sample_name")]
#   params = params[!grepl("path", params)]
#   return(params)
# }
