# #library('scatterpie')
# library('ggplot2')
# library('tidyr')
# library('ggbeeswarm')
# library('plotly')
# library('plyr')
# library('pals')
# library('reshape')
# library('base')
# library(comprehenr)





#' Size Plot
#'
#' Plots distribution of loxcode sizes found in sample
#'
#' @param lox loxcode experiment object of current experiment
#' @param sample loxcode sample object
#' @param count_matrix A matrix of the selected samples
#' @param code_set A matrix of the selected code sets
#' @param labels A string containing sample name or alias
#' @return a plot of the distribution of code sizes in the sample
#' @import tidyr
#' @import ggbeeswarm
#' @import dplyr
#' @import pals
#' @import comprehenr
#' @rdname size_plot
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' # Example usage
#' # Assuming lox, sample, count_matrix, code_set, and labels are defined
#' size_plot(lox, sample="Sample 1", count_matrix="all_samples", code_set="all_codes", labels = "alias")
setGeneric("size_plot", function(lox,
                                 sample,
                                 count_matrix = "all_samples",
                                 code_set = "all_codes",
                                 labels) {
    standardGeneric("size_plot")
})

#' @rdname size_plot
setMethod("size_plot", "loxcode_experiment", function(lox,
                                                      sample,
                                                      count_matrix = "all_samples",
                                                      code_set = "all_codes",
                                                      labels) {
    code <- NULL
    size <- NULL
    # If labels is equal to "alias", certain operations are performed,
    # including filtering the counts data based on the specified sample,
    # subset operations on codes and data,
    # and creating a ggplot object g1 with stacked bar plots.
    if (labels == "alias") {
        #print("hi1")
        # which(grepl("Sample 1", Week2@alias$all_samples))
        #temp = which(grepl(count_matrix, lox@alias$all_samples))
        counts <- lox@count_matrixes[[count_matrix]]
        #print("hi2")
        #print(counts)
        counts$codes <- row.names(counts)
        #print("hi3")
        # Week2@alias$all_samples[[1]][[2]]
        #View(lox@alias$all_samples)
        #Week2@alias$all_samples[Week2@alias$all_samples$alias == "Sample 1",]$sample_name
        tempp <- lox@alias$all_samples[lox@alias$all_samples$alias == sample, ]$sample_name
        #print(tempp)
        #print("hi4")
        #fin = lox@alias$count_matrix[[1]][[tempp]]
        #print("hi5")
        #print(fin)
        counts <- subset(counts, counts[, tempp] > 0)
        #print("hi6")
        codes <- lox@code_sets[[code_set]]
        #print("hi7")
        #print(codes)
        data <- subset(codes, code %in% counts$codes)
        #print("hi8")
        #print(data)
        #title = lables
        g1 <- ggplot(data) + geom_bar(aes(as.factor(size), fill = as.factor(is_valid)),
                                     position = "stack",
                                     stat = "count") +
            xlab("size") + ylab("diversity")


        #print(g1)
        return(g1)
    }
    # If labels is not equal to "alias", different operations are performed,
    # including filtering the counts data based on the specified sample,
    # subset operations on codes and data,
    # and creating a ggplot object g with stacked bar plots.
    else{
        counts <- lox@count_matrixes[[count_matrix]]
        counts$codes <- row.names(counts)
        counts <- subset(counts, counts[, sample] > 0)
        codes <- lox@code_sets[[code_set]]
        data <- subset(codes, code %in% counts$codes)
        #title = labels
    }

    title <- switch(labels,
                   "sample" = sample,
                   "alias" = get_alias(lox, count_matrix, sample))

    g <- ggplot(data) + geom_bar(aes(as.factor(size), fill = as.factor(is_valid)),
                                position = "stack",
                                stat = "count") +
        xlab("size") + ylab("diversity") +
        ggtitle(title)

    return(g)
})

# #' Size Plot (v2)
# #'
# #' Plots distribution of loxcode sizes found in sample. Works for loxcode_samples without loxcode_experiment
# #'
# #' @param sample loxcode sample object
# #' @return a plot of the distribution of code sizes in the sample
# #' @rdname size_plot_sample
# #' @export
# #' @examples
# #' # Load required packages
# #' library(loxcodeR)
# #'
# #' # Example usage
# #' # Assuming sample is a loxcode sample object
# #' size_plot_sample(sample)
# setGeneric("size_plot_sample", function(sample) {
#     standardGeneric("size_plot_sample")
# })
#
# #' @rdname size_plot_sample
# setMethod("size_plot_sample", "loxcode_sample", function(sample) {
#     size = NULL
#     dist_orig = NULL
#     g = ggplot(sample@decode@data) +
#         geom_bar(aes(as.factor(size), fill = as.factor(dist_orig)),
#                  position = "stack",
#                  stat = "count") +
#         xlab("size") + ylab("diversity") +
#         ggtitle(sample@name) +
#         labs(fill = "dist_orig")
#     return(g)
# })

#' Dist_orig Plot
#'
#' Plots distribution of distance from origin found in sample
#'
#' @param lox loxcode experiment object of current experiment
#' @param sample loxcode sample object
#' @param count_matrix A matrix of the selected samples
#' @param code_set A matrix of the selected code sets
#' @param labels A string containing sample name or alias
#' @return a plot of the distribution of distance from origin in the sample
#' @rdname dist_orig_plot
#' @import ggplot2
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 ggtitle
#' @importFrom scatterpie geom_scatterpie
#' @export
#' @examples
#' # Load a sample loxcode experiment object and sample
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#'
#' # Plot the distribution of distance from origin
#' dist_orig_plot(lox, sample="Sample 1", count_matrix= "all_samples", code_set= "all_codes", labels = "alias")

setGeneric("dist_orig_plot", function(lox,
                                      sample,
                                      count_matrix = "all_samples",
                                      code_set = "all_codes",
                                      labels = "alias") {
    standardGeneric("dist_orig_plot")
})

#' @rdname dist_orig_plot
setMethod("dist_orig_plot", "loxcode_experiment", function(lox,
                                                           sample,
                                                           count_matrix = "all_samples",
                                                           code_set = "all_codes",
                                                           labels = "alias") {
    code <- NULL
    dist_orig <- NULL
    # x = loxcode_sample
    # u <- valid(x)
    # #u <- u[u$size == size, ]
    # fill_scale <- scale_fill_manual(breaks = 0:15, values = rep('blue', 16)) # use twice since gradient is hard to see otherwise
    # g <- ggplot(data = u) + facet_wrap(~size) + geom_bar(aes(x = dist_orig, fill = factor(dist_orig)), show.legend = F) + scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +
    #   xlab("Distance from origin") + ylab("Diversity") + ggtitle(loxcodeR::name(x))

    # If labels is equal to "alias", certain operations are performed,
    # including filtering the counts data based on the specified sample,
    # subset operations on codes and data,
    # and creating a ggplot object g with faceted bar plots based on the "size" variable.
    if (labels == "alias") {
        counts <- lox@count_matrixes[[count_matrix]]
        counts$codes <- row.names(counts)
        tempp <- lox@alias$all_samples[lox@alias$all_samples$alias == sample, ]$sample_name
        counts <- subset(counts, counts[, tempp] > 0)
        codes <- lox@code_sets[[code_set]]
        codes[is.na(codes)] <- 0
        data <- subset(codes, code %in% counts$codes)
        g <- ggplot(data) + facet_wrap( ~ size, ncol = 2) + geom_bar(aes(dist_orig, fill = factor(dist_orig)), show.legend = TRUE) +
            scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +
            xlab("Distance from origin") + ylab("diversity")

        return(g)
    }

    # If labels is not equal to "alias", different operations are performed,
    # including filtering the counts data based on the specified sample,
    # subset operations on codes and data,
    # and creating a ggplot object g with faceted bar plots based on the "size" variable.
    else{
        counts <- lox@count_matrixes[[count_matrix]]
        counts$codes <- row.names(counts)
        counts <- subset(counts, counts[, sample] > 0)
        codes <- lox@code_sets[[code_set]]
        codes[is.na(codes)] <- 0
        data <- subset(codes, code %in% counts$codes)
    }


    title <- switch(labels,
                   "sample" = sample,
                   "alias" = get_alias(lox, count_matrix, sample))

    g <- ggplot(data) + facet_wrap( ~ size, ncol = 2) + geom_bar(aes(dist_orig, fill = factor(dist_orig)), show.legend = FALSE) +
        scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +
        xlab("Distance from origin") + ylab("diversity")
    ggtitle(title)

    return(g)
})






#' Barcode union
#'
#' Function to perform union of barcodes based on size or complexity
#'
#' @param rep1 loxcode object containing data of a sample
#' @param rep2 loxcode object containing data of a sample
#' @param type string containing "size" or "complexity"
#' @param range A list containg the start and end of size or dist_origin to perform union on
#' @return A unique list of the barcodes after performing union operation
#' @export
#' @importFrom dplyr filter
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assume rep1, rep2, type, and range are defined
#' rep1 <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")@samples$N712_N501
#' rep2 <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")@samples$N712_N506
#' barcode_union(rep1, rep2, type = "size", range = list(start = 10, end = 100))
barcode_union <- function(rep1, rep2, type = "size", range) {
    size <- NULL
    dist_orig <- NULL
    # if-else condition based on the value of the type argument. If type is equal to "size",
    # the function performs the following steps:
    # 1. It filters the rep1 object's decode@data based on the size values within the specified range.
    # 2. It filters the rep2 object's decode@data based on the size values within the specified range.
    # 3. It extracts the unique codes from both filtered results using the unique function.
    # 4. It combines the unique codes from both rep1 and rep2 into a single vector x.
    if (type == "size") {
        x <- unique(c(
            filter(rep1@decode@data, size >= range[1], size <= range[2])$code,
            filter(rep2@decode@data, size >= range[1], size <= range[2])$code
        ))
    }
    # If type is equal to "complexity", the function performs similar steps as above,
    # but filters the dist_orig values instead of the size values.
    else if (type == "complexity") {
        x <- unique(c(
            filter(
                rep1@decode@data,
                dist_orig >= range[1],
                dist_orig <= range[2]
            )$code,
            filter(
                rep2@decode@data,
                dist_orig >= range[1],
                dist_orig <= range[2]
            )$code
        ))
    }
    return (x)
}

#' Barcode Stats
#'
#' Number of times a barcode is seen in a sample
#'
#' @param union_bc Barcodes data to check in a sample
#' @param rep loxcode object containing data of a sample
#' @return A object containing the barcode and data on the matched barcode in the sample
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming union_bc and rep are defined
#' rep1 <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")@samples$N712_N501
#' rep2 <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")@samples$N712_N506
#' bc_union <- barcode_union(rep1, rep2, type="size", range=c(1, 13))
#' get_barcode_stats_rep(bc_union, rep1)
get_barcode_stats_rep <- function(union_bc, rep) {
    # the match function is used to find the indices of the elements in union_bc that match the
    # codes in rep@decode@data$code. These indices are stored in the index variable.
    index <- match(union_bc, rep@decode@data$code)
    # u <- loxcodeR::valid(rep)[index, ]
    # The function then selects the corresponding rows from rep@decode@data using the indices in
    # index and assigns the result to the variable u. This step retrieves the rows in
    # rep@decode@data that match the codes in union_bc
    u <- rep@decode@data[index,]
    # Next, any missing values in the count column of u are replaced with zero using the
    # assignment operator <- and the logical condition is.na(u$count)
    u$count[is.na(u$count)] <- 0
    # Finally, the code column in u is updated to have the values from union_bc for all rows,
    # and the modified u data frame is returned.
    u$code <- union_bc
    return(u)
}



#' Pair Comparison Plot 2
#'
#' Function to create a pair comparision plot by passing code set and sample set. This function can use custom sample and code set for plot creation.
#'
#' @param lox Loxcode experiment object
#' @param sampleset A matrix of the selected samples
#' @param codeset A matrix of the selected code sets
#' @param s1 First Loxcode sample
#' @param s2 Second Loxcode sample
#' @param colorBy Color points by size or dist_orig or firstread
#' @param labels Sample name or alias as labels
#' @param sizeRange Range of size available
#' @param dist_origRange Range of distance from origin available
#' @param firstreadRange Range of first read available
#' @return A plot containing pairwise comparison of the samples
#' @rdname pair_comparison_plot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming lox, sampleset, codeset, s1, and s2 are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' #s1 = switch("N709_N501","sample" = lox@samples[["N709_N501"]],"alias" = lox@samples[["N709_N501"]])
#' #s2 = switch("N709_N502","sample" = lox@samples[["N709_N502"]],"alias" = lox@samples[["N709_N502"]])
#' #pair_comparison_plot2(lox, sampleset = "all_samples", codeset = "all_codes", s1, s2, colorBy = "size", labels = "alias", sizeRange = NULL, dist_origRange = NULL, firstreadRange = NULL)
setGeneric("pair_comparison_plot2", function(lox,
                                             sampleset = "all_samples",
                                             codeset = "all_codes",
                                             s1,
                                             s2,
                                             colorBy = "size",
                                             labels = "alias",
                                             sizeRange = NULL,
                                             dist_origRange = NULL,
                                             firstreadRange = NULL) {
    standardGeneric("pair_comparison_plot2")
})

#' @rdname pair_comparison_plot2
setMethod("pair_comparison_plot2", "loxcode_experiment", function(lox,
                                                                  sampleset = "all_samples",
                                                                  codeset = "all_codes",
                                                                  s1,
                                                                  s2,
                                                                  colorBy = "size",
                                                                  labels = "alias",
                                                                  sizeRange = NULL,
                                                                  dist_origRange = NULL,
                                                                  firstreadRange = NULL) {
    size <- NULL
    dist_orig <- NULL
    firstread <- NULL
    # Labels
    name1 <- switch(labels,
                   "alias" = get_alias(lox, sampleset, s1@name),
                   "sample" = s1@name)
    name2 <- switch(labels,
                   "alias" = get_alias(lox, sampleset, s2@name),
                   "sample" = s2@name)

    # Comparison Table
    comp_table <- get_comparison_table2(s1, s2)


    # Subset Comparison Table
    comp_table <-
        comp_table[comp_table$code %in% lox@code_sets[[codeset]]$code,]
    if (!is.null(sizeRange)) {
        comp_table <-
            subset(comp_table, size >= sizeRange[1] & size <= sizeRange[2])
    }
    if (!is.null(dist_origRange)) {
        comp_table <-
            subset(comp_table,
                   dist_orig >= dist_origRange[1] & dist_orig <= dist_origRange[2])
    }
    if (!is.null(firstreadRange)) {
        comp_table <-
            subset(comp_table,
                   firstread >= firstreadRange[1] & firstread <= firstreadRange[2])
    }

    nonzero_mask <-
        (comp_table$s1_count > 0 & comp_table$s2_count > 0)
    s1_zero_mask <- (comp_table$s1_count == 0)
    s2_zero_mask <- (comp_table$s2_count == 0)

    theme_update(plot.title = element_text(hjust = 0.5))
    color <- switch(
        colorBy,
        "size" = as.factor(comp_table$size[nonzero_mask]),
        "dist_orig" = as.factor(comp_table$dist_orig[nonzero_mask]),
        "firstread" = log10(comp_table$firstread[nonzero_mask])
    )


    #Main plot
    #original
    # g = ggplot() +
    #   geom_point(aes(y = comp_table$log_s1[nonzero_mask] ,
    #                  x = comp_table$log_s2[nonzero_mask],
    #                  color = color)) +
    #   geom_abline(alpha = 0.75) +
    #   ylab(paste0('log10(1 + reads in ', name1, ')')) +
    #   xlab(paste0('log10(1 + reads in ', name2, ')')) +
    #   ggtitle(paste(name1, 'vs', name2))


    g <- ggplot() +
        #s1_count[nonzero_mask]
        geom_point(aes(
            y = comp_table$s1_count[nonzero_mask] ,
            x = comp_table$s2_count[nonzero_mask],
            color = color
        )) +
        geom_abline(alpha = 0.75) +
        ylab(paste0(name1)) +
        xlab(paste0(name2)) +
        ggtitle(paste(name1, 'vs', name2)) +
        scale_y_log10() + scale_x_log10()


    if (colorBy == "firstread") {
        g <- g + scale_color_continuous(name = colorBy)
    } else {
        g <- g + scale_color_discrete(name = colorBy)
    }

    # Barcodes only in one sample
    if (name1 != name2) {
        color1 <- switch(
            colorBy,
            "size" = as.factor(comp_table$size[s1_zero_mask]),
            "dist_orig" = as.factor(comp_table$dist_orig[s1_zero_mask]),
            "firstread" = log10(comp_table$firstread[s1_zero_mask])
        )
        color2 <- switch(
            colorBy,
            "size" = as.factor(comp_table$size[s2_zero_mask]),
            "dist_orig" = as.factor(comp_table$dist_orig[s2_zero_mask]),
            "firstread" = log10(comp_table$firstread[s2_zero_mask])
        )

        # g = g +
        #   ggbeeswarm::geom_quasirandom(aes(y = rep(-0.2, sum(s1_zero_mask)),
        #                             x = log10(1 + comp_table$s2_count[s1_zero_mask]),
        #                             color = color1), groupOnX = FALSE, width = 0.2) +
        #   ggbeeswarm::geom_quasirandom(aes(y = log10(1 + comp_table$s1_count[s2_zero_mask]),
        #                        x = rep(-0.2, sum(s2_zero_mask)), color = color2),
        #                    groupOnX = TRUE, width = 0.2)

        g <- g +
            ggbeeswarm::geom_quasirandom(
                aes(
                    y = 0.5,
                    x = comp_table$s2_count[s1_zero_mask],
                    color = color1
                ),
                groupOnX = FALSE,
                width = 0.2
            ) +
            ggbeeswarm::geom_quasirandom(
                aes(
                    y = comp_table$s1_count[s2_zero_mask],
                    x = 0.5,
                    color = color2
                ),
                groupOnX = TRUE,
                width = 0.2
            ) + scale_y_log10() + scale_x_log10()





        #
        #       ggbeeswarm::geom_quasirandom(aes(0.5, y = rep(-0.2, sum(s1_zero_mask)),
        #                            #x = log10(1 + comp_table$s2_count[s1_zero_mask]),
        #                            x = comp_table$s2_count[s1_zero_mask],
        #                            color = color1), groupOnX = FALSE, width = 0.2) +
        #       ggbeeswarm::geom_quasirandom(aes(0.5, #y = log10(1 + comp_table$s1_count[s2_zero_mask]),
        #                            y = comp_table$s1_count[s2_zero_mask],
        #                            x = rep(-0.2, sum(s2_zero_mask)), color = color2),
        #                        groupOnX = TRUE, width = 0.2)+
        #theme_minimal() +

    }
    #g <- ggplotly(g)
    return(g)
})

#' Pair Comparison Plot all
#'
#' @param lox Loxcode_experiment object
#' @param min_reads Code set
#' @param matrix First Loxcode_sample
#' @param codeset Second Loxcode_sample
#' @return A pair comparision plot comparing all samples agaist each other
#' @rdname pair_comparison_plot_all
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming lox, min_reads, matrix, and codeset are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' #pair_comparison_plot_all(lox,matrix = lox@count_matrixes$`all_samples_ST114_Bcells ST114_2 NN167_Week2`,codeset = lox@code_sets$invalid_codes_ST111)

setGeneric("pair_comparison_plot_all", function(lox,
                                                min_reads = 10,
                                                matrix = "all_samples",
                                                codeset = "all_codes")
{
    standardGeneric("pair_comparison_plot_all")
})

#' @rdname pair_comparison_plot_all
setMethod("pair_comparison_plot_all", "loxcode_experiment",
          function(lox,
                   min_reads = 10,
                   matrix = "all_samples",
                   codeset = "all_codes") {
              .xvalue <- NULL
              . <- NULL
              .yvalue <- NULL

              # The select variable is assigned the value of the matrix argument.
              # This represents the count matrix to be used for the plot.
              select <- matrix
              # The select matrix is modified such that any entry less than min_reads is set to 0.
              select[select < min_reads] <- 0

              # Rows in the select matrix with zero sums are removed using subset and rowSums.
              select <- subset(select, rowSums(select) > 0)
              #g = list()
              #marrangeGrob(react$pairs, nrow=2, ncol=2)

              # The gatherpairs function from the tidyverse package is used to convert the select matrix
              # into a long format suitable for plotting.
              select %>% gatherpairs(names(select)) %>%  {
                  #ggplot
                  ggplot(
                      subset(., .xvalue > 0 & .yvalue > 0),
                      aes(
                          x = .xvalue,
                          y = .yvalue,
                          size = .xvalue + .yvalue,
                          color = abs((.xvalue - .yvalue) / (.xvalue +
                                                                 .yvalue))
                      )
                  ) +
                      geom_point(alpha = 0.8) +
                      ggbeeswarm::geom_quasirandom(
                          aes(0.5, .yvalue , size = .yvalue, name = ""),
                          subset(., .xvalue == 0 &
                                     .yvalue > 0),
                          alpha = 0.8
                      ) +
                      ggbeeswarm::geom_quasirandom(
                          aes(.xvalue, 0.5, size = .xvalue),
                          subset(., .xvalue > 0 &
                                     .yvalue == 0),
                          groupOnX = FALSE,
                          alpha = 0.8
                      ) +
                      scale_size_area(max_size = 6) +
                      geom_smooth(method = 'lm') +
                      facet_wrap(.xkey ~ .ykey,
                                 ncol = length(unique(.$.ykey)),
                                 scales = 'free') +
                      scale_color_gradientn(colours = pals::brewer.purples(10),
                                            name = "") +
                      theme_minimal() + scale_y_log10(name = "") + scale_x_log10(name = "")
              }
              #ggplotly(g)
              # return(g)
          })

#' Gather Pairs
#'
#' Function used to
#'
#' @param data Data required for function
#' @param ... variable as quosures
#' @param na.rm default as FALSE for remove na values
#' @param convert default as False
#' @param factor_key default as false
#' @param xkey .xkey
#' @param xvalue .xvalue
#' @param ykey .ykey
#' @param yvalue .yvalue
#' @return A gathered pair of key pair values
#' @importFrom rlang quos
#' @importFrom rlang enquo
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr select
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 scale_fill_discrete
#' @importFrom ggplot2 scale_size_area
#' @importFrom ggplot2 scale_color_gradient
gatherpairs <- function(data,
                        ...,
                        xkey = '.xkey',
                        xvalue = '.xvalue',
                        ykey = '.ykey',
                        yvalue = '.yvalue',
                        na.rm = FALSE,
                        convert = FALSE,
                        factor_key = FALSE) {
    . <- NULL
    size <- NULL
    # Check if 'value' argument is passed in '...'
    if ("value" %in% names(list(...))) {
        # If 'value' is passed, capture its value and remove it from '...'
        value <- list(...)$value
        ...$value <- NULL
    } else {
        # If 'value' is not passed, set it to NULL
        value <- NULL
    }

    # The quos function is used to capture the variables specified in ... as quosures.
    vars <-
        quos(...)
    xkey <- enquo(xkey)
    xvalue <- enquo(xvalue)
    ykey <- enquo(ykey)

    yvalue <- enquo(yvalue)


    # The data frame data is piped into a series of transformations using the %>% operator.
    # The gather function is used to gather the variables specified in ... into key-value pairs,
    # with the key column assigned to xkey and the value column assigned to xvalue. The na.rm,
    # convert, and factor_key arguments are passed through to the gather function.
    # The result of the previous step is combined with the original data frame using cbind,
    # and the select function is used to select the variables specified in ....
    data %>% {
        cbind(
            gather(
                .,
                key = !!xkey,
                value = !!xvalue,
                !!!vars,
                na.rm = na.rm,
                convert = convert,
                factor_key = factor_key
            ),
            select(.,!!!vars)
        )
        # Another call to gather is made to gather the variables specified in ... into key-value
        # pairs, with the key column assigned to ykey and the value column assigned to yvalue.
        # The na.rm, convert, and factor_key arguments are again passed through to the gather function.
    } %>% gather(
        .,
        key = !!ykey,
        value = !!yvalue,
        !!!vars,
        na.rm = na.rm,
        convert = convert,
        factor_key = factor_key
    )
}

#' Comparison Table
#'
#' Function used to compare two samples based on barcode counts
#'
#' @param s1 loxcode object containing data of 1st sample
#' @param s2 loxcode object containing data of 2nd sample
#' @return A data frame of comparision of two samples
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' s1=lox@samples$N712_N501
#' s2=lox@samples$N712_N506
#' # Example usage
#' # Assuming s1 and s2 are defined
#' get_comparison_table2(s1, s2)
get_comparison_table2 <- function(s1, s2) {
    # table of barcodes in both samples
    # The counts data frame is created by merging the decode data of s1 and s2 based on the
    # common barcode ("code") using the merge function. The all=TRUE argument ensures that
    # all barcodes from both samples are included in the result. The columns related to
    # "dist_orig," "size," and "firstread" are extracted and stored in the counts data frame.
    counts <- merge(s1@decode@data,
                   s2@decode@data,
                   by = "code",
                   all = TRUE)
    counts$dist_orig <- rowMeans(counts[, grepl("dist_orig", names(counts))], na.rm = TRUE)
    counts$size <- rowMeans(counts[, grepl("size", names(counts))], na.rm = TRUE)
    counts$firstread <- rowMeans(counts[, grepl("firstread", names(counts))], na.rm = TRUE)
    # The counts data frame is modified by replacing any missing values with 0 using the is.na
    # and assignment operations.
    counts$count.x[is.na(counts$count.x)] <- 0
    counts$count.y[is.na(counts$count.y)] <- 0

    # The table data frame is created with columns for "code," "size," "dist_orig," "firstread,"
    # "s1_count," and "s2_count" using the data from the counts data frame.
    table <- data.frame(
        code = counts$code,
        size = counts$size,
        dist_orig = counts$dist_orig,
        firstread = counts$firstread,
        s1_count = counts$count.x,
        s2_count = counts$count.y,
        stringsAsFactors = FALSE
    )

    # scale by total number of reads



    # A mask called nonzero_mask is created to identify rows where both s1_count and s2_count are
    # greater than zero. Separate masks called s1_zero_mask and s2_zero_mask are created to
    # identify rows where either s1_count or s2_count is zero.
    nonzero_mask <- (table$s1_count > 0 & table$s2_count > 0)
    s1_zero_mask <- (table$s1_count == 0)
    s2_zero_mask <- (table$s2_count == 0)

    # log of counts
    # The log_s1 and log_s2 columns in the table data frame are populated with the corresponding
    # barcode counts from s1_count and s2_count for the rows where both counts are greater than
    # zero. These columns represent the log-transformed counts.
    table$log_s1[nonzero_mask] <- table$s1_count[nonzero_mask]
    table$log_s2[nonzero_mask] <- table$s2_count[nonzero_mask]

    # table$log_s1[nonzero_mask] = log10(1 + table$s1_count[nonzero_mask])
    # table$log_s2[nonzero_mask] = log10(1 + table$s2_count[nonzero_mask])

    return(table)
}

# #' Show recombination distance distribution by size as beeswarm plot
# #'
# #' Produces a beeswarm plot where for each cassette size, individual distinct barcodes are shown as
# #' points. Size and color of points correspond to the read count of each barcode.
# #' @param x loxcode_sample object
# #' @param count_threshold counts threshold, barcodes with a count number exceeding this threshold are ignored. This is in order
# #' to avoid having very large points resulting from barcodes with disproportionate read counts.
# #' @return A beeswarm plot of the recombinaton distance
# #' @rdname dist_count_beeswarm_plot
# #' @export
# #' @examples
# #' # Load required packages
# #' library(loxcodeR)
# #'
# #' # Example usage
# #' # Assume x and count_threshold are defined
# #' x <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
# #' dist_count_beeswarm_plot(x, count_threshold = 100)
# #'
# setGeneric("dist_count_beeswarm_plot", function(x, count_threshold) {
#     standardGeneric("dist_count_beeswarm_plot")
# })
#
# #' @rdname dist_count_beeswarm_plot
# setMethod("dist_count_beeswarm_plot", "loxcode_sample", function(x, count_threshold) {
#     count = NULL
#     dist_orig = NULL
#     size = NULL
#     # The function first applies the valid function from the loxcodeR package to the input
#     # sample x, and filters the results to include only barcodes with a "count" value below the
#     # specified "count_threshold". The filtered data is assigned to the variable y.
#     loxcodeR::valid(x) %>% filter(count < count_threshold) -> y
#     # The ggplot function is used to initialize the plotting object, with the data set to y.
#     # The ggbeeswarm::geom_quasirandom function is called to create a beeswarm plot.
#     g <- ggplot(data = y) +
#         ggbeeswarm::geom_quasirandom(
#             width = 0.9,
#             groupOnX = FALSE,
#             aes(
#                 x = dist_orig,
#                 y = size,
#                 size = count,
#                 color = count
#             ),
#             alpha = 0.2
#         ) +
#         scale_size_area("num_reads") +
#         scale_x_continuous("distance from origin", breaks = 0:15) +
#         scale_y_continuous("size", breaks = c(3, 5, 7, 9, 13)) +
#         scale_color_gradient(low = 'blue', high = 'red') +
#         ggtitle(loxcodeR::name(x))
#     return(g)
# })

#' Heatmap Plot
#'
#' Generate heatmap from count_matrix and also perform clustering while plotting
#'
#' @param loxcode_experiment loxcode_experiment object
#' @param count_matrix A matrix of the selected samples
#' @param code_set A matrix of the selected code sets
#' @param style string containing "ggplot" or "plotly"
#' @param labels A string containing the name by which to sample
#' @param clustering A string containing parameter by which clustering is to be performed
#' @param min_reads Integer containing value for min number of reads
#' @param max_repeats Integer containing value for max number of repeats
#' @param min_repeats Integer containing value for min number of repeats
#' @param split_by1 String containing parameter to perform first split by
#' @param split_by2 String containing parameter to perform second split by
#' @param agglomeration setting default as complete. Used for clustering
#' @return Heat map plot of the selected samples. Also produced with splits and clusters in some cases.
#' @rdname heatmap_plot
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_radius
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_step
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 scale_color_continuous
#' @importFrom ggplot2 scale_color_discrete
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom stats na.omit
#' @importFrom grid unit
#' @importFrom ggplot2 sec_axis
#' @importFrom ggplot2 scale_x_log10
#' @importFrom ggplot2 scale_y_log10
#' @importFrom ggplot2 theme_update
#' @importFrom plyr revalue
#' @importFrom plyr rbind.fill
#' @importFrom ggplot2 dup_axis
#' @importFrom methods new
#' @importFrom stats aggregate
#' @importFrom stats alias
#' @importFrom utils str
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming loxcode_experiment, count_matrix, and code_set are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' heatmap_plot(lox)
setGeneric("heatmap_plot", function(loxcode_experiment,
                                    count_matrix = "all_samples",
                                    code_set = "all_codes",
                                    style = "ggplot",
                                    labels = "sample",
                                    clustering = "none",
                                    agglomeration = "complete",
                                    min_reads = 0,
                                    max_repeats = 100,
                                    min_repeats = 1,
                                    split_by1 = "none",
                                    split_by2 = "none") {
    standardGeneric("heatmap_plot")
})

#' @rdname heatmap_plot
setMethod("heatmap_plot", "loxcode_experiment",
          function(loxcode_experiment,
                   count_matrix = "all_samples",
                   code_set = "all_codes",
                   style = "ggplot",
                   labels = "sample",
                   clustering = "none",
                   agglomeration = "complete",
                   min_reads = 0,
                   max_repeats = 100,
                   min_repeats = 1,
                   split_by1,
                   split_by2) {
              pos <- NULL
              code <- NULL
              variable <- NULL
              value <- NULL
              # xtracting the necessary parameters from the loxcode_experiment object and assigning them to local variables.
              x <- loxcode_experiment

              #clustering
              if (clustering != "none") {
                  cluster_row <- switch(
                      clustering,
                      "none" = FALSE,
                      "row" = TRUE,
                      "col" = FALSE,
                      "both" = TRUE
                  )
                  cluster_col <- switch(
                      clustering,
                      "none" = FALSE,
                      "row" = FALSE,
                      "col" = TRUE,
                      "both" = TRUE
                  )
              }
              else{
                  cluster_row <- FALSE
                  cluster_col <- FALSE
              }

              ##select only codes that are present in code_set
              # The count matrix is obtained from the loxcode_experiment object using the specified
              # count_matrix and code_set. Barcodes that do not meet the filtering criteria
              # (e.g., minimum reads, maximum repeats, minimum repeats) are removed from the matrix.
              m <- x@count_matrixes[[count_matrix]]
              m <- m[row.names(m) %in% x@code_sets[[code_set]]$code, ]
              m[m < min_reads] <- 0
              m <- m[rowSums(m) > 0, ]
              m <- m[rowSums(m > 0) < max_repeats, ]
              m <- m[rowSums(m > 0) >= min_repeats, ]
              m <- m[, colSums(m) > 0]
              #m=sweep(m,2,colSums(m),`/`)
              #print(m)


              if (cluster_row == TRUE) {
                  row.order <-
                      hclust(dist(m), method = agglomeration)$order
                  m <- m[row.order,]
              }
              if (cluster_col == TRUE) {
                  col.order <-
                      hclust(dist(t(m)), method = agglomeration)$order
                  m <- m[, col.order]
              }

              # If clustering is specified (not "none"), hierarchical clustering is performed on the rows
              # and/or columns of the count matrix based on the specified agglomeration method
              # (agglomeration). The rows and/or columns of the count matrix are reordered accordingly.
              if (cluster_row == FALSE &
                  cluster_col == FALSE & agglomeration == "binary")
              {
                  binary <- array()
                  for (i in seq_len(nrow(m)))
                      binary[i] <- paste0(as.numeric((m > 0)[i, ]), collapse = "")
                  m <- m[order(binary, rowSums(m)), ]
              }


              m$pos <- seq_len(nrow(m))
              m$code <- row.names(m)
              mm <- reshape::melt(m, id.vars = c("pos", "code"))
              m <- subset(m, select = -c(pos, code))


              mm[mm$value == 0, ]$value <- NA
              if (split_by1 != "none") {
                  new <- mm$variable
                  selectfrom <- x@meta[x@meta$sample_name %in% new, ]

                  mm[, split_by1] <- NA
                  for (i in unique(new)) {
                      m1 <- x@meta[x@meta$sample_name == i, ][, split_by1]
                      mm[mm$variable == i, ][, split_by1] <- m1
                      #print(mouse)
                  }
              }
              if (split_by2 != "none") {
                  new <- mm$variable
                  selectfrom <- x@meta[x@meta$sample_name %in% new, ]

                  mm[, split_by2] <- NA
                  for (i in unique(new)) {
                      m1 <- x@meta[x@meta$sample_name == i, ][, split_by2]
                      mm[mm$variable == i, ][, split_by2] <- m1
                      #print(mouse)
                  }
              }
              #print(mouse)
              labs <-
                  switch(
                      labels,
                      "sample" = names(m),
                      "alias" = lapply(names(m), get_alias, lox = x, set = count_matrix)
                  )
              #ggplot(mm)+
              #print(mm)
              g <- c()
              if (split_by2 != "none" &&
                  split_by1 != "none" && split_by2 != split_by1) {
                  g <- ggplot(mm, aes(x = split_by1, y = split_by2)) +
                      geom_tile(aes(
                          pos,
                          as.numeric(variable),
                          fill = log10(value + 1),
                          group = code
                      )) +
                      scale_fill_gradient2(
                          midpoint = 0,
                          high = "darkred",
                          mid = "lightgray",
                          name = "log10(counts+1)"
                      ) +
                      theme_minimal() +
                      xlab("") +
                      ylab("samples") + xlab("loxcodes") +
                      scale_x_continuous(expand = c(0, 0)) +
                      scale_y_continuous(
                          breaks = seq_along(labs),
                          labels = labs,
                          sec.axis = dup_axis(labels = colSums(m >
                                                                   0), name = "# codes")
                      ) +
                      theme(legend.position = "top") +
                      facet_wrap(c(vars(mm[, split_by1]), vars(mm[, split_by2])), scales = "free", ncol =
                                     2, )
              } else if (split_by1 != "none") {
                  #    for(i in unique(mm$mouse)){
                  #      g = g + ggplot(mm[mm$mouse==i,])+
                  #        geom_tile(aes(pos,as.numeric(variable),fill=log10(value+1),group=code))+
                  #        scale_fill_gradient2(midpoint = 0,high = "darkred",mid = "lightgray",name="log10(counts+1)")+
                  #        theme_minimal()+
                  #        xlab("")+
                  #        ylab("samples")+xlab("loxcodes")+
                  #        scale_x_continuous(expand = c(0, 0))+
                  #        scale_y_continuous(breaks=seq_along(labs),labels=labs,
                  #                           sec.axis = dup_axis(labels=colSums(m>0),name="# codes"))+
                  #        theme(legend.position="top")
                  #    }
                  g <- ggplot(mm) +
                      geom_tile(aes(
                          pos,
                          as.numeric(variable),
                          fill = log10(value + 1),
                          group = code
                      )) +
                      scale_fill_gradient2(
                          midpoint = 0,
                          high = "darkred",
                          mid = "lightgray",
                          name = "log10(counts+1)"
                      ) +
                      theme_minimal() +
                      xlab("") +
                      ylab("samples") + xlab("loxcodes") +
                      scale_x_continuous(expand = c(0, 0)) +
                      scale_y_continuous(
                          breaks = seq_along(labs),
                          labels = labs,
                          sec.axis = dup_axis(labels = colSums(m >
                                                                   0), name = "# codes")
                      ) +
                      theme(legend.position = "top") +
                      facet_wrap(vars(mm[, split_by1]), scales = "free", ncol =
                                     2, )
              } else if (split_by2 != "none") {
                  g <- ggplot(mm) +
                      geom_tile(aes(
                          pos,
                          as.numeric(variable),
                          fill = log10(value + 1),
                          group = code
                      )) +
                      scale_fill_gradient2(
                          midpoint = 0,
                          high = "darkred",
                          mid = "lightgray",
                          name = "log10(counts+1)"
                      ) +
                      theme_minimal() +
                      xlab("") +
                      ylab("samples") + xlab("loxcodes") +
                      scale_x_continuous(expand = c(0, 0)) +
                      scale_y_continuous(
                          breaks = seq_along(labs),
                          labels = labs,
                          sec.axis = dup_axis(labels = colSums(m >
                                                                   0), name = "# codes")
                      ) +
                      theme(legend.position = "top") +
                      facet_wrap(vars(mm[, split_by2]), scales = "free", ncol =
                                     2, )
              } else {
                  g <- ggplot(mm) +
                      geom_tile(aes(
                          pos,
                          as.numeric(variable),
                          fill = log10(value + 1),
                          group = code
                      )) +
                      scale_fill_gradient2(
                          midpoint = 0,
                          high = "darkred",
                          mid = "lightgray",
                          name = "log10(counts+1)"
                      ) +
                      theme_minimal() +
                      xlab("") +
                      ylab("samples") + xlab("loxcodes") +
                      scale_x_continuous(expand = c(0, 0)) +
                      scale_y_continuous(
                          breaks = seq_along(labs),
                          labels = labs,
                          sec.axis = dup_axis(labels = colSums(m >
                                                                   0), name = "# codes")
                      ) +
                      theme(legend.position = "top")
              }


              return(g)
          })


#' Generate bubble from count_matrix
#'
#' Generate bubble plot from count_matrix and also perform clustering while plotting
#'
#' @param loxcode_experiment loxcode_experiment object
#' @param count_matrix A matrix of the selected samples
#' @param code_set A matrix of the selected code sets
#' @param style string containing "ggplot" or "plotly"
#' @param labels A string containing the name by which to sample
#' @param clustering A string containing parameter by which clustering is to be performed
#' @param agglomeration Selection of the collection
#' @param min_reads Integer containing value for min number of reads
#' @param max_repeats Integer containing value for max number of repeats
#' @param min_repeats Integer containing value for min number of repeats
#' @param split_by1 String containing parameter to perform first split by
#' @param split_by2 String containing parameter to perform second split by
#' @return Bubble plot of the selected samples. Also produced with splits and clusters in some cases.
#' @rdname bubble_plot
#' @importFrom stats hclust
#' @importFrom stats dist
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_radius
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 facet_grid
#' @importFrom dplyr vars
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assume loxcode_experiment, count_matrix, and code_set are defined
#' loxcode_experiment <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' bubble_plot(loxcode_experiment,count_matrix="all_samples",code_set="all_codes",style="ggplot", labels="sample", clustering="none",agglomeration="complete",min_reads=0,max_repeats=100,min_repeats=1,split_by1="none",split_by2="none")

setGeneric("bubble_plot", function(loxcode_experiment,
                                   count_matrix = "all_samples",
                                   code_set = "all_codes",
                                   style = "ggplot",
                                   labels = "sample",
                                   clustering = "none",
                                   agglomeration = "complete",
                                   min_reads = 0,
                                   max_repeats = 100,
                                   min_repeats = 1,
                                   split_by1 = "none",
                                   split_by2 = "none") {
    standardGeneric("bubble_plot")
})
#' @rdname bubble_plot
setMethod("bubble_plot", "loxcode_experiment", function(loxcode_experiment,
                                                        count_matrix = "all_samples",
                                                        code_set = "all_codes",
                                                        style = "ggplot",
                                                        labels = "sample",
                                                        clustering = "none",
                                                        agglomeration = "complete",
                                                        min_reads = 0,
                                                        max_repeats = 100,
                                                        min_repeats = 1,
                                                        split_by1,
                                                        split_by2) {
    x <- loxcode_experiment

    variable <- NULL
    index <- NULL
    code <- NULL
    value <- NULL
    #clustering
    cluster_row <- switch(
        clustering,
        "none" = FALSE,
        "row" = TRUE,
        "col" = FALSE,
        "both" = TRUE
    )
    cluster_col <- switch(
        clustering,
        "none" = FALSE,
        "row" = FALSE,
        "col" = TRUE,
        "both" = TRUE
    )

    ##select only codes that are present in code_set
    m <- x@count_matrixes[[count_matrix]]
    m <- m[row.names(m) %in% x@code_sets[[code_set]]$code, ]
    m[m < min_reads] <- 0

    m <- m[rowSums(m) > 0, ]
    m <- m[rowSums(m > 0) < max_repeats, ]
    m <- m[rowSums(m > 0) >= min_repeats, ]
    m <- m[, colSums(m) > 0]
    m <- sweep(m, 2, colSums(m), `/`)

    if (cluster_row == TRUE) {
        row.order <-
            hclust(dist(m), method = agglomeration)$order
        m <- m[row.order,]
    }
    if (cluster_col == TRUE) {
        col.order <-
            hclust(dist(t(m)), method = agglomeration)$order
        m <- m[, col.order]
    }

    if (cluster_row == FALSE &
        cluster_col == FALSE & agglomeration == "binary")
    {
        binary <- array()
        for (i in seq_len(nrow(m)))
            binary[i] <- paste0(as.numeric((m > 0)[i, ]), collapse = "")
        m <- m[order(binary, rowSums(m)), ]
        print(table(binary))
    }

    m$index <- seq_len(nrow(m))
    m$code <- row.names(m)
    mm <- reshape::melt(m, id.var = c("index", "code"))

    if (split_by1 != "none") {
        new <- mm$variable
        selectfrom <- x@meta[x@meta$sample_name %in% new, ]

        mm[, split_by1] <- NA
        for (i in unique(new)) {
            m1 <- x@meta[x@meta$sample_name == i, ][, split_by1]
            mm[mm$variable == i, ][, split_by1] <- m1
            #print(mouse)
        }
    }
    if (split_by2 != "none") {
        new <- mm$variable
        selectfrom <- x@meta[x@meta$sample_name %in% new, ]

        mm[, split_by2] <- NA
        for (i in unique(new)) {
            m1 <- x@meta[x@meta$sample_name == i, ][, split_by2]
            mm[mm$variable == i, ][, split_by2] <- m1
            #print(mouse)
        }
    }
    labs <-
        switch(
            labels,
            "sample" = names(m),
            "alias" = lapply(names(m), get_alias, lox = x, set = count_matrix)
        )

    if (split_by2 != "none" &&
        split_by1 != "none" && split_by2 != split_by1) {
        ggplot() +
            geom_hline(
                aes(yintercept = as.vector(seq_len(10) %o% 10 ^ (0:-7))),
                color = "lightgray",
                size = 0.1,
                alpha = 0.3
            ) +
            scale_color_manual(
                values = unname(base::sample(
                    pals::brewer.dark2(nrow(m)), replace = TRUE
                )),
                breaks = rownames(m),
                guide = FALSE
            ) +
            geom_point(
                aes(
                    factor(variable),
                    index,
                    color = factor(code),
                    size = log(value)
                ),
                data = mm,
                alpha = 0.8
            ) +
            theme(legend.position = "none") +
            scale_radius(range = c(1, 3)) +
            theme_minimal() +
            theme(plot.title = element_text(size = 7),
                  legend.position = "none") +
            xlab("samples") + ylab("") +
            scale_x_discrete(labels = labs) +
            theme(
                panel.grid.major.x = element_blank() ,
                panel.grid.minor.y = element_blank() ,
                axis.text = element_text(size = 10)
                #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
            ) + coord_flip() +
            facet_grid(vars(mm[, split_by1]), vars(mm[, split_by2]), scales = "free")
    } else if (split_by1 != "none") {
        ggplot() +
            geom_hline(
                aes(yintercept = as.vector(seq_len(10) %o% 10 ^ (0:-7))),
                color = "lightgray",
                size = 0.1,
                alpha = 0.3
            ) +
            scale_color_manual(
                values = unname(base::sample(
                    pals::brewer.dark2(nrow(m)), replace = TRUE
                )),
                breaks = rownames(m),
                guide = FALSE
            ) +
            geom_point(
                aes(
                    factor(variable),
                    index,
                    color = factor(code),
                    size = log(value)
                ),
                data = mm,
                alpha = 0.8
            ) +
            theme(legend.position = "none") +
            scale_radius(range = c(1, 3)) +
            theme_minimal() +
            theme(plot.title = element_text(size = 7),
                  legend.position = "none") +
            xlab("samples") + ylab("") +
            scale_x_discrete(labels = labs) +
            theme(
                panel.grid.major.x = element_blank() ,
                panel.grid.minor.y = element_blank() ,
                axis.text = element_text(size = 10)
                #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
            ) + coord_flip() +
            facet_grid(vars(mm[, split_by1]))
    } else {
        ggplot() +
            geom_hline(
                aes(yintercept = as.vector(seq_len(10) %o% 10 ^ (0:-7))),
                color = "lightgray",
                size = 0.1,
                alpha = 0.3
            ) +
            scale_color_manual(
                values = unname(base::sample(
                    pals::brewer.dark2(nrow(m)), replace = TRUE
                )),
                breaks = rownames(m),
                guide = FALSE
            ) +
            geom_point(
                aes(
                    factor(variable),
                    index,
                    color = factor(code),
                    size = log(value)
                ),
                data = mm,
                alpha = 0.8
            ) +
            theme(legend.position = "none") +
            scale_radius(range = c(1, 3)) +
            theme_minimal() +
            theme(plot.title = element_text(size = 7),
                  legend.position = "none") +
            xlab("samples") + ylab("loxcodes") +
            scale_x_discrete(labels = labs) +
            theme(
                panel.grid.major.x = element_blank() ,
                panel.grid.minor.y = element_blank() ,
                axis.text = element_text(size = 10)
                #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
            ) + coord_flip()
    }

})



#' Read stats plot
#'
#' Generate size, complexity, both, ratio plot for the samples
#'
#' @param x loxcode experiment object
#' @param count_matrix A matrix of the selected samples
#' @param code_set A matrix of the selected code sets
#' @param plot specify which plot (size/complexity/both/ratio)
#' @param fill specify whether reads should be normalized in boxplots
#' @param labels specify if graph should be labelled by sample name or alias
#' @return A plot of the selected type explaing details of the samples
#' @importFrom ggplot2 margin
#' @rdname readstats_plot_old
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming x, count_matrix, code_set, plot, fill, and labels are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' readstats_plot_old(lox)
setGeneric("readstats_plot_old", function(x,
                                          count_matrix = "all_samples",
                                          code_set = "all_codes",
                                          plot = "size",
                                          fill = TRUE,
                                          labels = "sample") {
    standardGeneric("readstats_plot_old")
})

#' @rdname readstats_plot_old
setMethod("readstats_plot_old", "loxcode_experiment", function(x,
                                                               count_matrix = "all_samples",
                                                               code_set = "all_codes",
                                                               plot = "size",
                                                               fill = TRUE,
                                                               labels = "sample") {
    size <- NULL
    count <- NULL
    dist_orig <- NULL
    reads <- NULL
    ratio <- NULL
    pos <- "stack"
    if (fill == TRUE)
        pos <- "fill"
    P <- ggplot()

    if (plot == "size") {
        #iterating over samples and sizes and count # reads per condition
        d1 <- data.frame()

        for (i in x@samples) {
            if (i@name %in% names(x@count_matrixes[[count_matrix]])) {
                dd <- i@decode@data

                for (j in c(seq_len(13))) {
                    ddd <- data.frame(
                        sample = name(i),
                        size = j,
                        count = sum(subset(dd, size == j)$count)
                    )
                    d1 <- plyr::rbind.fill(d1, ddd)

                    print(ddd)
                }
            }
        }
        d1 <- subset(d1, count > 0)
        P <- (
            ggplot(d1) + geom_bar(
                aes(
                    fill = factor(size),
                    y = count,
                    x = sample
                ),
                position = pos,
                stat = "identity"
            ) +
                scale_fill_discrete(name = "size") +
                theme_minimal()
        )
    }

    if (plot == "complexity") {
        #iterating over samples and complexity and count # reads per condition
        d2 <- data.frame()

        for (i in x@samples) {
            if (i@name %in% names(x@count_matrixes[[count_matrix]])) {
                dd <- i@decode@data

                for (j in c(0:15)) {
                    ddd <- data.frame(
                        sample = name(i),
                        dist_orig = j,
                        count = sum(subset(dd, dist_orig == j)$count)
                    )
                    d2 <- plyr::rbind.fill(d2, ddd)

                }
            }
        }
        d2 <- subset(d2, count > 0)
        P <- (
            ggplot(d2) + geom_bar(
                aes(
                    fill = factor(dist_orig),
                    y = count,
                    x = sample
                ),
                position = pos,
                stat = "identity"
            ) +
                scale_fill_discrete(name = "complexity") +
                theme_minimal()
        )
    }

    ## plotting reads for both complexity and size
    if (plot == "both") {
        #iterating over sizes and complexities
        d3 <- data.frame()

        for (i in x@samples) {
            if (i@name %in% names(x@count_matrixes[[count_matrix]])) {
            dd <- i@decode@data

            for (j in c(0:15))
                for(k in c(seq_len(13))) {
                    ddd <- data.frame(
                        sample = name(i),
                        dist_orig = j,
                        size = k,
                        reads = sum(subset(
                            dd, dist_orig == j & size == k
                        )$count)
                    )
                    d3 <- plyr::rbind.fill(d3, ddd)

                }
            }
        }

        d3 <- subset(d3, reads > 0)
        P <- (
            ggplot(d3) + facet_wrap( ~ sample, ncol = 2,scales="free") + geom_point(aes(
                factor(size),
                factor(dist_orig),
                size = reads,
                color = (reads)
            )) +
                scale_fill_discrete(name = "complexity")
            + xlab("size") + ylab("complexity") +
                theme(
                    plot.margin = margin(
                        t = 0,
                        # Top margin
                        r = 5,
                        # Right margin
                        b = 5,
                        # Bottom margin
                        l = 5
                    ),
                    panel.spacing.y = unit(5, "pt"),
                ) + xlab("size") + ylab("complexity") +
                scale_size_area() + scale_color_gradient(low = "blue", high = "red")
        )
    }


    ## plotting ratio barcodes/reads
    if (plot == "ratio") {
        #iterating over sizes and complexities
        d4 <- data.frame()

        for (i in x@samples) {
            if (i@name %in% names(x@count_matrixes[[count_matrix]])) {
                dd <- i@decode@data

                for (j in c(0:15))
                    for(k in c(seq_len(13))) {
                        ddd <- data.frame(
                            sample = name(i),
                            dist_orig = j,
                            size = k,
                            ratio = nrow(subset(
                                dd, dist_orig == j &
                                    size == k
                            )) / sum(subset(
                                dd, dist_orig == j & size == k
                            )$count)
                        )
                        d4 <- plyr::rbind.fill(d4, ddd)

                    }
            }
        }
        #print("-------------------------------here-------------------------------")
        d4 <- subset(d4, ratio > 0)
        P <- (
            ggplot(d4) + facet_wrap( ~ sample, ncol = 2, scales = "free") + geom_point(aes(
                factor(size),
                factor(dist_orig),
                size = ratio,
                color = (ratio)
            )) +
                scale_fill_discrete(name = "complexity") +
                theme(
                    plot.margin = margin(
                        t = 0,
                        # Top margin
                        r = 5,
                        # Right margin
                        b = 5,
                        # Bottom margin
                        l = 5
                    ),
                    panel.spacing.y = unit(5, "pt"),
                ) + xlab("size") + ylab("complexity") +
                scale_size_area() + scale_color_gradient(low = "blue", high = "red")
        )
    }
    P
})

# ggplot()+geom_histogram(aes(NN167@samples[["N701-N505"]]@decode@data$size), stat="count")

#' readstats plot
#'
#' Generate size, complexity, both, ratio plot for the samples
#'
#' @param loxcode_experiment loxcode_experiment object
#' @param count_matrix A matrix of the selected samples
#' @param code_set A matrix of the selected code sets
#' @param labels labels to use in the plot
#' @param plot specify which plot (size/complexity/both/ratio)
#' @param fill specify whether reads should be normalized in boxplots
#' @return A plot of the selected type explaing details of the samples
#' @rdname readstats_plot2
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming loxcode_experiment, count_matrix, code_set, plot, fill, and labels are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' readstats_plot2(lox)
setGeneric("readstats_plot2", function(loxcode_experiment,
                                       count_matrix = "all_samples",
                                       code_set = "all_codes",
                                       plot = "size",
                                       fill = TRUE,
                                       labels = "alias") {
    standardGeneric("readstats_plot2")
})

#' @rdname readstats_plot2
setMethod("readstats_plot2", "loxcode_experiment", function(loxcode_experiment,
                                                            count_matrix = "all_samples",
                                                            code_set = "all_codes",
                                                            plot = "size",
                                                            fill = TRUE,
                                                            labels = "sample") {
    count <- NULL
    dist_orig <- NULL
    ratio <- NULL
    x <- loxcode_experiment
    # graph config
    pos <- "stack"
    if (fill == TRUE) {
        pos <- "fill"
    }

    # initialize
    p <- ggplot()
    counts <- x@count_matrixes[[count_matrix]]
    codes <- x@code_sets[[code_set]]
    samples <- names(counts)

    ### SIZE PLOT
    if (plot == "size") {
        sizeBySample <- data.frame()
        # create size count table
        for (i in seq_along(samples)) {
            for (j in c(seq_len(13))) {
                row <- data.frame(
                    sample = samples[[i]],
                    alias = get_alias(x, count_matrix,
                                      samples[[i]]),
                    size = j,
                    count = 0
                )
                sizeBySample <- plyr::rbind.fill(sizeBySample, row)
            }
        }

        # fill size count table
        for (i in seq_len(nrow(counts))) {
            row <- counts[i,]
            code <- rownames(row)
            if (code %in% codes$code) {
                size <- length(gregexpr(" ", code)[[1]]) + 1 # one more than number of spaces
                for (j in names(row)) {
                    curr <- row[[j]]
                    index <- which(sizeBySample$sample == j &
                                      sizeBySample$size == size)
                    prev <- sizeBySample$count[[index]]
                    sizeBySample$count[[index]] <- curr + prev
                }
            }
        }

        sizeBySample <- subset(sizeBySample, count > 0)
        p <- ggplot(sizeBySample) +
            geom_bar(
                aes(
                    fill = factor(size),
                    y = count,
                    x = switch(labels, "alias" = alias, "sample" = sample)
                ),
                position = pos,
                stat = "identity"
            ) +
            scale_fill_discrete(name = "size") +
            labs(x = "sample")
        saveRDS(p, "lastplot.rds")
    }

    ### COMPLEXITY PLOT
    if (plot == "complexity") {
        # create complexity count table
        complexityBySample <- data.frame()
        for (i in seq_along(samples)) {
            for (j in c(0:15)) {
                row <- data.frame(
                    sample = samples[[i]],
                    alias = get_alias(x,
                                      count_matrix,
                                      samples[[i]]),
                    dist_orig = j,
                    count = 0
                )
                complexityBySample <- plyr::rbind.fill(complexityBySample, row)
            }
        }

        # fill complexity count table
        for (i in seq_len(nrow(counts))) {
            row <- counts[i,]
            code <- rownames(row)
            if (code %in% codes$code) {
                dist <- codes$dist_orig[[match(code, codes$code)]]
                if (is.na(dist)) {
                    dist <- 0
                }
                for (j in names(row)) {
                    curr <- row[[j]]
                    index <- which(
                        complexityBySample$sample == j & complexityBySample$dist_orig == dist
                    )
                    prev <- complexityBySample$count[[index]]
                    complexityBySample$count[[index]] <- curr + prev
                }
            }
        }

        complexityBySample <- subset(complexityBySample, count > 0)
        p <- ggplot(complexityBySample) +
            geom_bar(
                aes(
                    fill = factor(dist_orig),
                    y = count,
                    x = switch(labels, "alias" = alias, "sample" = sample)
                ),
                position = pos,
                stat = "identity"
            ) +
            scale_fill_discrete(name = "complexity") +
            labs(x = "sample") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    }

    ### SIZE COMPLEXITY PLOT
    if (plot == "both") {
        # create counts table
        bothBySample <- data.frame()
        for (i in seq_along(samples)) {
            for (j in 0:15)
                for (k in seq_len(13)) {
                    row <- data.frame(
                        sample = samples[[i]],
                        alias = get_alias(x, count_matrix, samples[[i]]),
                        dist_orig = j,
                        size = k,
                        count = 0
                    )
                    bothBySample <- plyr::rbind.fill(bothBySample, row)
                }
        }

        # fill counts table
        for (i in seq_len(nrow(counts))) {
            row <- counts[i,]
            code <- rownames(row)
            if (code %in% codes$code) {
                size <- length(gregexpr(" ", code)[[1]]) + 1 # one more than number of spaces
                dist <- codes$dist_orig[[match(code, codes$code)]]
                if (is.na(dist)) {
                    dist <- 0
                }
                for (j in names(row)) {
                    curr <- row[[j]]
                    index <- which(
                        bothBySample$sample == j &
                            bothBySample$dist_orig == dist & bothBySample$size == size
                    )
                    prev <- bothBySample$count[[index]]
                    bothBySample$count[[index]] <- curr + prev
                }
            }
        }

        bothBySample <- subset(bothBySample, count > 0)

        p <- (
            ggplot(bothBySample) + facet_wrap( ~ switch(
                labels, "alias" = alias, "sample" = sample
            )) + geom_point(aes(
                factor(size),
                factor(dist_orig),
                size = count,
                color = count
            )) +
                scale_fill_discrete(name = "complexity") +
                theme_minimal() + xlab("size") + ylab("complexity") +
                scale_size_area() + scale_color_gradient(low = "blue", high =
                                                             "red")
        )
    }

    ### RATIO PLOT
    if (plot == "ratio") {
        # create ratio table
        ratioBySample <- data.frame()
        for (i in seq_along(samples)) {
            for (j in 0:15)
                for (k in seq_len(13)) {
                    row <- data.frame(
                        sample = samples[[i]],
                        alias = get_alias(x, count_matrix, samples[[i]]),
                        dist_orig = j,
                        size = k,
                        ratio = 0,
                        ncodes = 0,
                        count = 0
                    )
                    ratioBySample <- plyr::rbind.fill(ratioBySample, row)
                }
        }
        # fill ratio table
        for (i in seq_len(nrow(counts))) {
            row <- counts[i,]
            code <- rownames(row)
            if (code %in% codes$code) {
                size <- length(gregexpr(" ", code)[[1]]) + 1 # one more than number of spaces
                dist <- codes$dist_orig[[match(code, codes$code)]]
                if (is.na(dist)) {
                    dist <- 0
                }
                for (j in names(row)) {
                    curr <- row[[j]]
                    index <- which(
                        ratioBySample$sample == j &
                            ratioBySample$dist_orig == dist & ratioBySample$size == size
                    )
                    prev <- ratioBySample$count[[index]]
                    ratioBySample$count[[index]] <- curr + prev
                    ratioBySample$ncodes[[index]] <- ratioBySample$ncodes[[index]] + 1
                    ratioBySample$ratio[[index]] <- ratioBySample$ncodes[[index]] / ratioBySample$count[[index]]
                }
            }
        }

        ratioBySample$count <- NULL
        ratioBySample$ncodes <- NULL
        ratioBySample <- subset(ratioBySample, ratio > 0)
        ratioBySample <- subset(ratioBySample, ratio != Inf)
        ratioBySample
        print("******************here**************************")
        p <- (
            ggplot(ratioBySample) + facet_wrap( ~ switch(
                labels, "alias" = alias, "sample" = sample
            ), ncol = 3) + geom_point(aes(
                factor(size),
                factor(dist_orig),
                size = ratio,
                color = (ratio)
            )) +
                scale_fill_discrete(name = "complexity") +
                theme_minimal() + xlab("size") + ylab("complexity") +
                scale_size_area() + scale_color_gradient(low = "blue", high = "red")
        )

    }


    return (p)
})

# #' sample comparison pie
# #'
# #' Pie comparion plot of the samples
# #'
# #' @param x loxcode_experiment object
# #' @param count_matrix A matrix of the selected samples
# #' @param code_set A matrix of the selected code sets
# #' @param scale scale the size of pies
# #' @param labels by sample or alias
# #' @return A plot containg comparision of the samples in the form of a pie
# #' @export
# setGeneric("sample_comparison_pie", function(x,count_matrix="all_samples",code_set="all_codes", scale=1, labels="sample") {standardGeneric("sample_comparison_pie")})
#
# setMethod("sample_comparison_pie", "loxcode_experiment", function(x,count_matrix="all_samples",code_set="all_codes", scale=1, labels="sample"){
#
#   M=x@count_matrixes[[count_matrix]]
#   M=M[row.names(M)%in%x@code_sets[[code_set]]$code,]
#
#   m=data.frame();
#
#   for(i in 1:ncol(M))for(j in 1:ncol(M)){
#     both=sum(M[,i]>0 & M[,j]>0,na.rm=TRUE)
#     in.1=sum(M[,i]>0 & M[,j]==0,na.rm=TRUE)
#     in.2=sum(M[,j]>0 & M[,i]==0,na.rm=TRUE)
#     d=data.frame(sample1=i, sample2=j,both=both,in.1=in.1,in.2=in.2,stringsAsFactors = FALSE)
#     m=plyr::rbind.fill(m,d)
#   }
#   library('scatterpie')
#   sample_names = names(x@count_matrixes[[count_matrix]])
#   labels = switch(labels, "sample"=sample_names, "alias"=lapply(sample_names, get_alias, lox=x, set=count_matrix))
#   print(unique(m$sample1))
#   print(unique(m$sample2))
#   g=ggplot(m)+
#     scatterpie::geom_scatterpie(aes(x=sample1,y=sample2,r=(both+in.1+in.2)/(1200/scale)),data=m, cols=c(3:5),size=0.1)+
#     scale_x_continuous(breaks=unique(m$sample1),labels=labels)+
#     scale_y_continuous(breaks=unique(m$sample2),labels=labels)+
#     theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1))+
#     scale_fill_discrete(name="")
#   print('hello')
#   return(g)
#
# })

#' Plotting the saturation of samples
#'
#' @param loxcode_experiment loxcode_experiment object
#' @param loxcode_sample A matrix of the selected samples
#' @param code_set A matrix of the selected code sets
#' @return Saturation plot of the loxcode sample
#' @rdname saturation_plot
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming loxcode_experiment, loxcode_sample, and code_set are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' saturation_plot(lox,"N712_N501")
setGeneric("saturation_plot", function(loxcode_experiment,
                                       loxcode_sample,
                                       code_set = "all_codes") {
    standardGeneric("saturation_plot")
})

#' @rdname saturation_plot
setMethod("saturation_plot", "loxcode_experiment", function(loxcode_experiment,
                                                            loxcode_sample,
                                                            code_set = "all_codes") {
    code <- NULL
    firstread <- NULL
    index <- NULL
    dist_orig <- NULL
    x <- loxcode_experiment
    i <- x@samples[[loxcode_sample]]
    cs <- i@decode@data

    cs <- subset(cs, code %in% x@code_sets[[code_set]]$code)
    cs <- cs[order(cs$firstread), ]
    cs$index <- c(seq_len(nrow(cs)))
    new <- data.frame(
        firstread = cs$firstread,
        index = cs$index,
        dist_orig = cs$dist_orig
    )

    g <-
        ggplot(new) + geom_step(aes(firstread, index, color = factor(dist_orig))) +
        geom_point(aes(firstread, index, color = factor(dist_orig)))
    # g <- ggplot(new) + geom_line(aes(firstread1,index1,color=dist_orig1)) + geom_line(aes(firstread2,index2,color=dist_orig2))+scale_colour_manual()
    return (g)
})

#' Multiple lined saturation plot
#'
#' @param loxcode_experiment loxcode_experiment object
#' @param loxcode_samples A matrix of the selected samples
#' @param codesets A matrix of the selected code sets
#' @param colorby colour by dist_orig or sample
#' @return Saturation plot of the loxcode sample
#' @rdname saturation_multi
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming loxcode_experiment, loxcode_samples, and codesets are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' saturation_multi(lox,"N712_N501","all_codes")

setGeneric("saturation_multi", function(loxcode_experiment,
                                        loxcode_samples,
                                        codesets,
                                        colorby = "dist_orig") {
    standardGeneric("saturation_multi")
})

#' @rdname saturation_multi
setMethod("saturation_multi", "loxcode_experiment", function(loxcode_experiment,
                                                             loxcode_samples,
                                                             codesets,
                                                             colorby = "dist_orig") {
    firstread <- NULL
    index <- NULL
    dist_orig <- NULL
    x <- loxcode_experiment
    data <- data.frame()

    # create the data.frame
    for (i in seq_along(loxcode_samples)) {
        s <- loxcode_samples[[i]]
        c <- codesets[[i]]
        sample <- x@samples[[s]]
        counts <- sample@decode@data
        counts$dist_orig[is.na(counts$dist_orig)] <- 0
        counts <- subset(counts, counts$code %in% x@code_sets[[codesets[[i]]]]$code)
        counts <- counts[order(counts$firstread),]
        counts$index <- c(seq_len(nrow(counts)))
        one_sample_data <- data.frame(
            firstread = counts$firstread,
            index = counts$index,
            dist_orig = counts$dist_orig,
            sample = s
        )
        colour <- switch(colorby,
                        "dist_orig" = one_sample_data$dist_orig,
                        "sample" = one_sample_data$sample)

        if (i == 1) {
            g <- ggplot(one_sample_data, aes(firstread, index)) +
                geom_line(data = one_sample_data) +
                geom_point(data = one_sample_data, aes(firstread, index,
                                                       colour = switch(
                                                           colorby,
                                                           "dist_orig" = factor(dist_orig),
                                                           "sample" = sample
                                                       ))) +
                scale_color_discrete(name = "")

        }
        else {
            g <- g +
                geom_line(data = one_sample_data) +
                geom_point(data = one_sample_data,
                           aes(firstread, index, colour = switch(
                               colorby,
                               "dist_orig" = factor(dist_orig),
                               "sample" = sample
                           )))
        }
    }
    return(g)
})




#' Sample Table
#'
#' Tabulates the barcode counts, total number of reads, maximum complexity and consensus filtered for each sample in an experiment
#'
#' @param x loxcode_experiment object
#' @param c count_matrix of samples to display
#' @return A table of the samples selected
#' @rdname sample_table
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming x and c are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' sample_table(lox,"all_samples")
setGeneric("sample_table", function(x, c) {
    standardGeneric("sample_table")
})

#' @rdname sample_table
setMethod("sample_table", "loxcode_experiment", function(x, c) {
    d <- data.frame()
    for (i in seq_along(x@samples)) {
        curr_sample <- x@samples[[i]]
        if (name(curr_sample) %in% names(x@count_matrixes[[c]])) {
            temp <- data.frame(
                "Sample_Name" = name(curr_sample),
                curr_sample@meta,
                stringsAsFactors = FALSE
            )
            d <- plyr::rbind.fill(d, temp)
        }
    }
    return(d)
})

strCap <- function(y) {
    c <- strsplit(y, " ")[[1]]
    paste(toupper(substring(c, 1, 1)),
          substring(c, 2),
          sep = "",
          collapse = " ")
}

capitalize <- function(w) {
    cap <- c()
    for (word in w) {
        cap <- append(cap, strCap(word))
    }
    return (cap)
}

#' Summary Table
#'
#' Tabulates the barcode counts, total number of reads, maximum complexity and consensus filtered for each sample in an experiment
#'
#' @param x loxcode_experiment object
#' @param c count_matrix of samples to display
#' @return A table of the selected samples
#' @rdname summary_table1
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming x and c are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' summary_table1(lox,"all_samples")
setGeneric("summary_table1", function(x, c = "all_samples") {
    standardGeneric("summary_table1")
})

#' @rdname summary_table1
setMethod("summary_table1", "loxcode_experiment", function(x, c = "all_samples") {
    d <- data.frame()

    if (!length(x@count_matrixes[[c]])) {
        return()
    }

    # check if samples have been collapsed
    collapsed <- TRUE
    if (length(intersect(names(x@samples), names(x@count_matrixes[[c]]))) != 0) {
        collapsed <- FALSE
    }

    # extract data from loxcode experiment
    if (collapsed == FALSE) {
        for (i in seq_along(x@samples)) {
            curr_sample <- x@samples[[i]]
            metadata <- curr_sample@meta
            names(metadata) <- capitalize(names(metadata))
            if (name(curr_sample) %in% names(x@count_matrixes[[c]])) {
                temp <- data.frame(
                    "Sample_Name" = name(curr_sample),
                    "Alias" = x@alias[[c]]$alias[[match(name(curr_sample), x@alias[[c]]$sample_name)]],
                    metadata,
                    "Barcode_Count" = length(curr_sample@decode@data$count),
                    "Invalid_Count" = length(curr_sample@decode@data$count) - sum(curr_sample@decode@data$is_valid, na.rm =
                                                                                      FALSE),
                    "Number_of_Reads" = curr_sample@decode_stats$tot_reads,
                    "Max_Complexity" = max(na.omit(
                        curr_sample@decode@data$dist_orig
                    )),
                    "Consensus_Filtered" = curr_sample@decode_stats$consensus_filtered,
                    "Percent_Filtered" = round(
                        100 * curr_sample@decode_stats$consensus_filtered / curr_sample@decode_stats$tot_reads,
                        2
                    ),
                    stringsAsFactors = FALSE
                )
                d <- plyr::rbind.fill(d, temp)
            }
        }
    }
    else {
        counts <- x@count_matrixes[[c]]
        all_codes <- x@code_sets[["all_codes"]]
        invalid_codes <- x@code_sets[["invalid_codes"]]
        meta <- get_collapsed_meta(x, c)
        for (i in seq_len(ncol(counts))) {
            curr <- names(counts)[i]
            barcodes <- row.names(counts[which (counts[[curr]] > 0),])
            metadata <- meta[curr,]
            names(metadata) <- capitalize(names(meta))
            temp <- data.frame(
                "Sample_Name" = curr,
                "Alias" = x@alias[[c]]$alias[[match(curr, x@alias[[c]]$sample_name)]],
                metadata,
                "Barcode_Count" = sum(counts[, i] > 0),
                "Invalid_Count" = sum(barcodes %in% row.names(invalid_codes)),
                "Number_of_Reads" = round(sum(counts[, i]), 0),
                "Max_Complexity" = max(all_codes[which (all_codes$code %in% barcodes),]$dist_orig, na.rm =
                                           TRUE),
                stringsAsFactors = FALSE
            )
            d <- plyr::rbind.fill(d, temp)
        }
        names(d)[names(d) == "metadata"] <- names(meta)
    }


    return(d)
})

#' Summary table (v2)
#' @param lox loxcode_experiment object
#' @param s sample set
#' @return table of sample information
#' @rdname summary_table2
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming lox and s are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' summary_table2(lox, "all_samples")
setGeneric("summary_table2", function(lox, s = "all_samples") {
    standardGeneric("summary_table2")
})

#' @rdname summary_table2
setMethod("summary_table2", "loxcode_experiment", function(lox, s = "all_samples") {
    counts <- lox@count_matrixes[[s]]
    aliases <- lox@alias[[s]]
    all_codes <- lox@code_sets[["all_codes"]]
    invalid_codes <- lox@code_sets[["invalid_codes"]]
    meta <- get_collapsed_meta(lox, s)
    d <- data.frame()
    for (i in seq_len(ncol(counts))) {
        curr <- names(counts)[i]
        barcodes <- row.names(counts[which (counts[[curr]] > 0),])
        metadata <- meta[curr,]

        if (curr %in% names(lox@samples)) {
            stats <- lox@samples[[curr]]@decode_stats
            nreads <- stats$tot_reads
            cfiltered <- stats$consensus_filtered
            pfiltered <- round(100 * cfiltered / nreads, 2)
        } else {
            nreads <- NA
            cfiltered <- round(sum(counts[, i]), 0)
            pfiltered <- NA
        }

        temp <- data.frame(
            "Sample_Name" = curr,
            "Alias" = aliases$alias[match(curr, aliases$sample_name)],
            metadata,
            "Barcode_Count" = sum(counts[, i] > 0),
            "Invalid_Count" = sum(barcodes %in% row.names(invalid_codes)),
            "Number_of_Reads" = nreads,
            "Max_Complexity" = max(all_codes[which (all_codes$code %in% barcodes),]$dist_orig, na.rm =
                                       TRUE),
            "Consensus_Filtered" = cfiltered,
            "Percent_Filtered" = pfiltered,
            stringsAsFactors = FALSE
        )
        d <- plyr::rbind.fill(d, temp)
    }
    d$Max_Complexity <- replace(d$Max_Complexity, is.infinite(d$Max_Complexity), 0)
    return(d)
})

#' codeset table
#'
#' Tabulates the codes in a codeset
#'
#' @param x loxcode_experiment object
#' @param n name of codeset to display
#' @return dataframe of codeset information
#' @rdname codeset_table
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assume x and n are defined
#' x <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' n="all_codes"
#' codeset_table(x, n)
#'
setGeneric("codeset_table", function(x, n) {
    standardGeneric("codeset_table")
})

#' @rdname codeset_table
setMethod("codeset_table", "loxcode_experiment", function(x, n) {
    if ((!is.null(x@code_sets$all_codes$size))) {
        x@code_sets$all_codes$size <- as.integer(x@code_sets$all_codes$size)
    }

    if ((!is.null(x@code_sets$all_codes$dist_orig))) {
        x@code_sets$all_codes$dist_orig <-
            as.integer(x@code_sets$all_codes$dist_orig)
    }

    if ((!is.null(x@code_sets$all_codes$is_valid))) {
        x@code_sets$all_codes$is_valid <-
            as.integer(x@code_sets$all_codes$is_valid)

        x@code_sets$all_codes$is_valid <- factor(x@code_sets$all_codes$is_valid,
                                                labels = c(unique(x@code_sets$all_codes$is_valid)))


    }

    table <- x@code_sets[[n]]
    # if (!is.null(d@code_sets$all_codes$Size)){
    #   as.integer(d@code_sets$all_codes$Size)
    # }
    names(table) <- capitalize(names(table))
    return(table)
})

#' Experiments Table
#'
#' Presents a list of loxcode experiments in table form
#'
#' @param exp list of loxcode experiment objects
#' @return a data frame containing the information of loxcode experiment
#' @rdname exp_table
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assume exp is a list of loxcode experiment objects
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' exp_table(c(lox))
#'
setGeneric("exp_table", function(exp) {
    standardGeneric("exp_table")
})

#' @rdname exp_table
setMethod("exp_table", "list", function(exp) {
    d <- data.frame()

    if (length(exp) > 0) {
        for (i in seq_along(exp)) {
            directories <- ""
            curr <- exp[[i]]
            for (dir in curr@dir) {
                directories <- paste(directories, dir)
            }
            temp <- data.frame(
                "Experiment_Name" = curr@name,
                "Directory" = directories,
                "Samples" = length(curr@samples),
                "Code_Sets" = length(curr@code_sets),
                "Sample_Sets" = length(curr@count_matrixes),
                stringsAsFactors = FALSE
            )
            d <- plyr::rbind.fill(d, temp)
        }
    }
    return(d)
})

#' Code frequency table
#'
#' Creates a data.frame of the frequency for each size and complexity combination
#'
#' @param x loxcode experiment object
#' @param s A matrix of the selected samples
#' @param c A matrix of the selected code sets
#' @return data frame of code frequencies
#' @rdname code_freq_table
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assume x, s, and c are defined
#' x <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' code_freq_table(x, s="all_samples", c="all_codes")
#'
setGeneric("code_freq_table", function(x, s = "all_samples", c = "all_codes") {
    standardGeneric("code_freq_table")
})

#' @rdname code_freq_table
setMethod("code_freq_table", "loxcode_experiment", function(x, s = "all_samples", c =
                                                                "all_codes") {
    size <- NULL
    dist_orig <- NULL
    codeset <- x@code_sets[[c]]
    counts <- x@count_matrixes[[s]]
    Y <- data.frame()

    if (!length(counts)) {
        return()
    } # prevents LoxCodeR from crashing
    #if (has_warning(max(codeset$dist_orig, na.rm=TRUE))) {print("crash prevented"); return() }

    # create a column in code sets recording frequency of barcodes
    index <- match(row.names(counts), codeset$code)
    index <- subset(index,!is.na(index))
    codeset$rep[index] <- rowSums(counts > 0)

    # create a data.frame with repetition frequencies for each size and complexity
    for (i in c(seq_len(max(codeset$size, na.rm = TRUE))))
        for (j in c(seq_len(max(codeset$dist_orig, na.rm = TRUE)))) {
            y <- data.frame(size = i, dist_orig = j)
            z <- as.data.frame(table(subset(codeset, size == i &
                                               dist_orig == j)$rep), stringsAsFactors = FALSE)
            if (nrow(z) == 0)
                next

            n <- z$Var1

            z <- as.data.frame(z[, -1]) ##to-do: find a way to order count entries
            y <- cbind(y, sum(z), t(z))
            names(y) <- c("size", "dist_orig", "radius", n)
            Y <- plyr::rbind.fill(Y, y)

        }
    Y[is.na(Y)] <- 0

    return(Y)
})

#' Filtered Code freq table
#'
#' Data frame of filtered codes specified by user parameters
#'
#' @param x loxcode_experiment object
#' @param s independent sample set
#' @param c A matrix of the selected code sets
#' @param t tolerance threshold for repeated proportions
#' @param m maximum repeats tolerated
#' @return scatterpie of filtered codes
#' @rdname filtered_codes_table
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assume x, s, c, t, and m are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' filtered_codes_table(lox, s = "all_samples", c = "all_codes", t = 5, m = 3)
setGeneric("filtered_codes_table", function(x,
                                            s = "all_samples",
                                            c = "all_codes",
                                            t = 0.05,
                                            m = 2) {
    standardGeneric("filtered_codes_table")
})

#' @rdname filtered_codes_table
setMethod("filtered_codes_table", "loxcode_experiment", function(x,
                                                                 s = "all_samples",
                                                                 c = "all_codes",
                                                                 t = 0.05,
                                                                 m = 2) {
    # data-frame of repetition frequencies
    Y <- code_freq_table(x, s, c)

    if (!length(x@code_sets[[c]]) |
        !length(x@count_matrixes[[s]])) {
        print("errorrrr")
        return()
    }

    YY <- data.frame()
    total <- max(as.numeric(names(Y[, !names(Y) %in% c("size", "dist_orig", "radius")])))
    for (i in seq_len(nrow(Y))) {
        row <- Y[i, ]
        proportion <- sum(row[as.character(seq(m, total))] / sum(row[as.character(seq(1, total))]))
        if (proportion * 100 <= t) {
            YY <- plyr::rbind.fill(YY, row)
        }
    }
    return(YY)
})

#' Code frequency pie
#'
#' Display the frequency proportions of codes of every permutation of size and complexity
#'
#' @param x loxcode_experiment object
#' @param s sample set of independent samples
#' @param c A matrix of the selected code sets
#' @return scatterpie plot of the codes in every permutation of size and complexity
#' @rdname code_frequency_pie
#' @importFrom scatterpie geom_scatterpie
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme_bw
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#' library(scatterpie)
#'
#' # Example usage
#' # Assuming x, s, and c are defined and have the required structure
#' x <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' code_frequency_pie(x, s="all_samples", c="all_codes")
setGeneric("code_frequency_pie", function(x, s = "all_samples", c = "all_codes") {
    standardGeneric("code_frequency_pie")
})

#' @rdname code_frequency_pie
setMethod("code_frequency_pie", "loxcode_experiment", function(x, s = "all_samples", c =
                                                                   "all_codes") {
    size <- NULL
    dist_orig <- NULL
    radius <- NULL
    Y <- code_freq_table(x, s, c)
    size <- Y$size
    dist_orig <- Y$dist_orig
    radius <- Y$radius
    # scatterpie plot
    g <- ggplot(Y) +
        scatterpie::geom_scatterpie(
            aes(
                x = size,
                y = dist_orig,
                r = log10(radius) / 10
            ),
            data = Y,
            cols = c(4:ncol(Y)),
            size = 0.1
        ) +
        scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13)) + scale_y_continuous(breaks =
                                                                                       c(seq_len(10))) +
        theme_bw() + scale_fill_discrete(name = "Repeats")

    return(g)
})

#' Filtered code frequency pie
#'
#' Filters out codes that have an overlap between independent samples
#'
#' @param x loxcode_experiment object
#' @param s independent sample set
#' @param c code set
#' @param t tolerance threshold for repeated proportions
#' @param m maximum repeats tolerated
#' @return scatter pie of filtered codes
#' @rdname filtered_codes_pie
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme_bw
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assume x, s, c, t, and m are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' filtered_codes_pie(lox, s = "all_samples", c = "all_codes", t = 5, m = 3)
setGeneric("filtered_codes_pie", function(x,
                                          s = "all_samples",
                                          c = "all_codes",
                                          t = 0.05,
                                          m = 2) {
    standardGeneric("filtered_codes_pie")
})

#' @rdname filtered_codes_pie
setMethod("filtered_codes_pie", "loxcode_experiment", function(x,
                                                               s = "all_samples",
                                                               c = "all_codes",
                                                               t = 0.05,
                                                               m = 2) {
    Y <- code_freq_table(x, s, c)
    YY <- filtered_codes_table(x, s, c, t, m)

    size <- YY$size
    dist_orig <- YY$dist_orig
    radius <- YY$radius
    g <- ggplot(YY) +
        scatterpie::geom_scatterpie(
            aes(
                x = size,
                y = dist_orig,
                r = log10(radius) / 5
            ),
            data = YY,
            cols = c(4:ncol(YY)),
            size = 0.1
        ) +
        scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13)) + scale_y_continuous(breaks =
                                                                                       c(seq_len(max(YY$dist_orig)) + 1)) +
        theme_bw() + scale_fill_discrete(name = "Repeats")

    return(g)
})

#' Get sample alias
#'
#' @param lox loxcode_experiment object
#' @param set sample set name
#' @param name sample_name
#' @return sample alias
#' @rdname get_alias
#' @export
#' @examples
#' # Load necessary libraries and data
#' library(loxcodeR)
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#'
#' # Get sample alias
#' alias <- get_alias(lox,"N712_N501","all_sample")
setGeneric("get_alias", function(lox, set, name) {
    standardGeneric("get_alias")
})

#' @rdname get_alias
setMethod("get_alias", "loxcode_experiment", function(lox, set, name) {
    aliases <- lox@alias[[set]]
    index <- match(name, aliases$sample_name)
    alias <- NA
    if (!is.na(index)) {
        alias <- aliases$alias[[index]]
    }
    return (alias)
})

#' Get sample name from alias
#'
#' @param lox loxcode_experiment object
#' @param set sample set name
#' @param name alias
#' @return sample alias
#' @rdname get_samplename
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming lox, set, and name are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' get_samplename(lox, set = "all_samples", name = "Sample 1")
setGeneric("get_samplename", function(lox, set, name) {
    standardGeneric("get_samplename")
})

#' @rdname get_samplename
setMethod("get_samplename", "loxcode_experiment", function(lox, set, name) {
    aliases <- lox@alias[[set]]
    index <- match(name, aliases$alias)
    samplename <- NA
    if (!is.na(index)) {
        samplename <- aliases$sample_name[[index]]
    }
    return (samplename)
})

# #' Density distributions
# #'
# #' Plot the density of the distribution of dist_orig or sizes in simulated samples
# #'
# #' @param ref reference sample_distribution
# #' @param npois poisson mean
# #' @param dist_type sample distribution from "sample" or "poisson"
# #' @param samples list of simulated samples
# #' @param plot_type density plot of "size" or "dist_orig"
# #' @return density plot of the distribution of dist_orig or sizes
# #' @rdname density_plot
# #' @importFrom ggplot2 geom_density
# #' @importFrom ggplot2 ggtitle
# #' @export
# #' @examples
# #' # Load required packages
# #' library(loxcodeR)
# #' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
# #' # Example usage
# #' # Assume ref, npois, dist_type, samples, and plot_type are defined
# #' density_distributions(ref, npois, dist_type, samples, plot_type = "size")
# setGeneric("density_plot", function(ref,
#                                     npois,
#                                     dist_type = "sample",
#                                     samples = list(),
#                                     plot_type = "dist_orig") {
#     standardGeneric("density_plot")
# })
#
# #' @rdname density_plot
# setMethod("density_plot", "loxcode_sample", function(ref,
#                                                      npois,
#                                                      dist_type = "sample",
#                                                      samples = list(),
#                                                      plot_type = "dist_orig") {
#     size = NULL
#     dist_orig = NULL
#     data = getData(ref, plot_type)
#     for (sample in samples) {
#         data = plyr::rbind.fill(data, getData(sample[[1]], plot_type))
#     }
#
#     if (dist_type == "sample") {
#         g = ggplot(data, aes(
#             x = switch(plot_type, "size" = size, "dist_orig" = dist_orig),
#             fill = source,
#             colour = source
#         )) +
#             xlab(plot_type) +
#             geom_density(alpha = 0.3) +
#             ggtitle(paste0(plot_type, " density plot"))
#     }
#     else if (dist_type == "poisson") {
#         g = ggplot()
#     }
#
#     return (g)
# })

#' Get sample data
#'
#' Get a list of sizes from the sample
#'
#' @param sample loxcode sample object
#' @param plot_type "size" or "dist_orig"
#' @return data frame with sizes in a sample
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming sample is defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' getData(lox@samples$N712_N501)
getData <- function(sample, plot_type = "size") {
    size <- NULL
    dist_orig <- NULL
    data <- sample@decode@data
    if (plot_type == "size") {
        d <- data.frame("size" = unique(data$size), "count" = NA)
        d$count <- vapply(d$size, function(x)
            sum(subset(data, size == x)$count), FUN.VALUE = numeric(1))
        v <- data.frame("size" = as.vector(rep(d$size, d$count)))
    }
    if (plot_type == "dist_orig") {
        d <- data.frame("dist_orig" = unique(data$dist_orig),
                        "count" = NA)
        d$count <- vapply(d$dist_orig, function(x)
            sum(subset(data, dist_orig == x)$count), FUN.VALUE = numeric(1))
        v <- data.frame("dist_orig" = as.vector(rep(d$dist_orig, d$count)))
    }
    v$source <- sample@name
    return(v)
}

## plot number of read

#' Read Plot
#'
#' plots the number of reads
#'
#' @param lox loxcode_experiment object
#' @param code_set code set in selected samples
#' @param count_matrix A matrix of the selected samples
#' @param labels specify if graph should be labelled by sample name or alias
#' @return A plot of the number of reads
#' @rdname read_plot
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming lox, code_set, and count_matrix are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' read_plot(lox)

setGeneric("read_plot",
           function(lox,
                    count_matrix = "all_samples",
                    code_set = "all_codes",
                    labels = "alias") {
               standardGeneric("read_plot")
           })

#' @rdname read_plot
setMethod("read_plot",
          "loxcode_experiment",
          function(lox,
                   count_matrix = "all_samples",
                   code_set = "all_codes",
                   labels = "sample") {
              # initialize variables
              counts <- lox@count_matrixes[[count_matrix]]
              codes <- lox@code_sets[[code_set]]

              #Crerate table
              samp_name <- names(counts)

              sampleAlias <- comprehenr::to_vec(for (x in samp_name)
                  get_alias(lox, count_matrix, x))

              num_read <- colSums(counts)
              read_table <- data.frame(num_read)

              # generate plot
              p <- ggplot(data = read_table,
                          aes(#x = samp_name,
                              x = switch(
                                  labels, "alias" = sampleAlias, "sample" = samp_name
                              ),
                              y = num_read)) +
                  geom_bar(stat = "identity", fill = "#66B2FF") +
                  labs(x = "Sample", y = "Num of Reads") +
                  scale_y_continuous(expand = c(0, 0)) +
                  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
              #fig <- ggplotly(p) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

              return (p)
          })


#' Generate bar code table
#'
#'A function used to generate table for the barcode details.
#'
#' @param lox loxcode_experiment object
#' @param code_set code set of the selected sample
#' @param count_matrix A matrix of the selected samples
#' @param labels specify if graph should be labelled by sample name or alias
#' @return A table containing the data about the bar codes
#' @rdname barcode_table
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assume lox, code_set, and count_matrix are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' barcode_table(lox, code_set= "all_codes", count_matrix= "all_samples", labels = "alias")

setGeneric("barcode_table",
           function(lox,
                    count_matrix = "all_samples",
                    code_set = "all_codes",
                    labels = "alias") {
               standardGeneric("barcode_table")
           })

#' @rdname barcode_table
setMethod("barcode_table",
          "loxcode_experiment",
          function(lox,
                   count_matrix = "all_samples",
                   code_set = "all_codes",
                   labels = "sample") {
              # initialize variables
              counts <- lox@count_matrixes[[count_matrix]]

              #Crerate table
              samp_name <- names(counts)

              alias <- comprehenr::to_vec(for (x in samp_name)
                  get_alias(lox, count_matrix, x))

              num_read <- colSums(counts)
              num_barcode <- colSums(counts > 0)
              barcode_table <- data.frame(alias, num_read, num_barcode)


              return (barcode_table)
          })


# #' correlation plot
# #'
# #' Function used to create a correlation plot
# #'
# #' @param lox loxcode_experiment object
# #' @param count_matrix A matrix of the selected samples
# #' @param code_set A matrix of the selected code sets
# #' @param method_ correlation method (p)
# #' @return ggplot object of plot
# #' @rdname correlation_plot
# #' @importFrom stats cor
# #' @importFrom stats heatmap
# #' @export
#
# setGeneric("correlation_plot", function(lox,
#                                         count_matrix = "all_samples",
#                                         code_set = "all_codes",
#                                         # parameter,
#                                         # data,
#                                         method_ = "pearson") {
#     standardGeneric("correlation_plot")
# })
#
# #' @rdname correlation_plot
# setMethod("correlation_plot", "loxcode_experiment",
#           function(lox,
#                    count_matrix = "all_samples",
#                    code_set = "all_codes",
#                    # parameter,
#                    # data,
#                    method_ = "pearson") {
#               if (is.null(lox) || is.null(count_matrix) || is.null(method_) ||
#                   is.null(code_set))
#                   # || is.null(parameter) || is.null(data) || data=="")
#               {
#                   return ()
#               }
#
#               # initialize data
#               counts = lox@count_matrixes[[count_matrix]]
#               codes = lox@code_sets[[code_set]]
#               # samp_table = lox@meta
#               # params = getMetaParameters(lox)
#               #
#               # # select required meta data
#               # primary_data = samp_table[[parameter]][match(names(counts), samp_table$sample_name)]
#               # other_parameters = params[!params %in% c(parameter)]
#               #
#               # # update sample names in count_matrix
#               # names = paste(primary_data)
#               # for (param in other_parameters) {
#               #   row = samp_table[[param]][match(names(counts), samp_table$sample_name)]
#               #   names = paste(names, row)
#               # }
#               # names(counts) = names
#
#               # subset from count_matrix relevant data
#               counts[counts < 10] = 0
#               counts = counts[rowSums(counts) > 0,]
#               #counts = counts[, primary_data %in% c(data)]
#               counts = counts[order(row.names(counts)), order(names(counts))]
#
#               # subset based on code complexity and flip distance
#               #codes = codes[codes$dist_orig > 4 & codes$flp_dist > 1, ]
#               counts = counts[row.names(counts) %in% codes$code,]
#
#
#               # shared = (sum(rowSums(counts[, 1:2]) > 0 &
#               # rowSums(counts[, 3:4]) > 0) / sum(rowSums(counts[, 1:2]) > 0 |
#               #                                     rowSums(counts[, 3:4]) > 0))
#
#               # find correlations
#               #counts = as.matrix(counts)
#               COR = as.data.frame(cor(log(counts + 1), method = method_))
#               if (nrow(COR) == 0) {
#                   return ()
#               }
#
#               # plot the data
#               heatmap(as.matrix(COR),
#                       #Rowv = NA,
#                       #Colv = NA,
#                       # main = paste0(parameter, " #", data, " barcode shared P5/P6: ", 100 * (round(shared, 2)), "%"),
#                       #main = paste0(parameter, " #", data),
#                       scale = "none")
#
#           })

# getMetaParameters <- function(lox) {
#   params = names(lox@meta)[!names(lox@meta) %in% c("sample_name")]
#   params = params[!grepl("path", params)]
#   return(params)
# }


#' readstats plot (v3)
#'
#' Function to create size, complexity, both or ratio plot
#'
#' @param lox loxcode_experiment object
#' @param count_matrix Count matrix of loxcode
#' @param code_set Code set of loxcode
#' @param plot specify which plot (size/complexity/both/ratio)
#' @param fill specify whether reads should be normalized in boxplots
#' @param labels specify if graph should be labelled by sample name or alias
#' @return Plot based on the type of plot selected
#' @rdname readstats_plot
#' @export
#' @examples
#' # Load required packages
#' library(loxcodeR)
#'
#' # Example usage
#' # Assuming lox, count_matrix, code_set, plot, and fill are defined
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/data-2024-05-20.rds")
#' readstats_plot(lox)

setGeneric("readstats_plot",
           function(lox,
                    count_matrix = "all_samples",
                    code_set = "all_codes",
                    plot = "size",
                    fill = TRUE,
                    labels = "alias") {
               standardGeneric("readstats_plot")
           })

#' @rdname readstats_plot
setMethod("readstats_plot",
          "loxcode_experiment",
          function(lox,
                   count_matrix = "all_samples",
                   code_set = "all_codes",
                   plot = "size",
                   fill = TRUE,
                   labels = "sample") {
              # initialize variables
              counts <- lox@count_matrixes[[count_matrix]]
              codes <- lox@code_sets[[code_set]]
              pos <- ifelse(fill, "fill", "stack")

              # generate plot
              if (plot == "ratio") {
                  p <- get_ratio_plot(lox, counts, codes, count_matrix, labels)
              } else {
                  p <- stats_bar_plot(lox, counts, codes, count_matrix, plot, pos, labels)
              }
              return (p)
          })

#' stats bar plot
#' helper function for readstats_plot to generate size and complexity plots
#'
#' @param counts count_matrix selected
#' @param codes code_set selected
#' @param count_matrix name of count_matrix selected
#' @param plot plot type (size/complexity)
#' @param pos fill type
#' @param lox loxcode object
#' @param labels specify if samples should be labelled by sample_name or alias
#' @return ggplot object of stats bar plot

stats_bar_plot <-
    function(lox,
             counts,
             codes,
             count_matrix,
             plot,
             pos,
             labels) {
        # initiliaze variables
        data <- switch(
            plot,
            "size" = "size",
            "complexity" = "dist_orig",
            "both" = c("size", "dist_orig")
        )
        samples <- names(counts)

        # create counts table
        countsBySample <- getCountsTable(lox, samples, codes, count_matrix, counts, data)

        # create plot
        if (plot %in% c("size", "complexity")) {
            p <- getCountsPlot(countsBySample, labels, pos, data)
        } else {
            p <- getBothPlot(countsBySample, labels)
        }

        return(p)
    }

#' get counts table
#' helper function to generate table of counts for each sample by selected data
#'
#' @param samples list of samples by name
#' @param codes code_set selected
#' @param count_matrix name of count_matrix selected
#' @param counts count_matrix selected
#' @param lox loxcode object
#' @param data data to count (size/complexity)
#' @return Table of counts for each sample

getCountsTable <-
    function(lox,
             samples,
             codes,
             count_matrix,
             counts,
             data) {
        countsBySample <- data.frame()
        for (i in seq_along(samples)) {
            sampleName <- samples[i]
            sampleAlias <- get_alias(lox, count_matrix, sampleName)
            indices <- codes$code %in% row.names(counts)[counts[[sampleName]] >
                                                            0]
            countsTable <- data.frame(table(codes[data][indices, ]))
            if (nrow(countsTable) != 0) {
                row <- data.frame("sample" = sampleName,
                                 "alias" = sampleAlias,
                                 countsTable)
                countsBySample <- plyr::rbind.fill(countsBySample, row)
            }
        }
        return (countsBySample)
    }

#' get ratio plot
#' helper function to generate plot of ratios
#'
#' @param lox loxcode_experiment
#' @param labels list of labels to samples
#' @param codes code_set selected
#' @param count_matrix name of count_matrix selected
#' @param counts count_matrix selected
#' @return ggplot object of plot
#' @rdname get_ratio_plot

setGeneric("get_ratio_plot",
           function(lox,
                    counts,
                    codes,
                    count_matrix,
                    labels) {
               standardGeneric("get_ratio_plot")
           })

#' @rdname get_ratio_plot
setMethod("get_ratio_plot",
          "loxcode_experiment",
          function(lox,
                   counts,
                   codes,
                   count_matrix,
                   labels) {
              # initialize
              samples <- names(counts)
              data <- c("size", "dist_orig")

              # get table of ratios by size and complexity
              ratioBySample <- data.frame()
              for (i in seq_along(samples)) {
                  sampleName <- samples[i]
                  sampleAlias <- get_alias(lox, count_matrix, sampleName)
                  codes <- lox@samples[[sampleName]]@decode@data
                  indices <- codes$code %in% row.names(counts)[counts[[sampleName]] >
                                                                  0]
                  subsetCodes <- codes[indices, ]
                  countsBySample <- data.frame(table(subsetCodes[data]))
                  if (nrow(countsBySample) != 0) {
                      totalsBySample <- aggregate(
                          subsetCodes$count,
                          by = list(Category = subsetCodes$size, subsetCodes$dist_orig),
                          FUN = sum
                      )
                      names(totalsBySample) <- c("size", "dist_orig", "total")
                      ratioTable <- merge(countsBySample, totalsBySample, by =
                                             data)
                      ratioTable$Freq <- ratioTable$Freq / ratioTable$total
                      row <- data.frame("sample" = sampleName,
                                       "alias" = sampleAlias,
                                       ratioTable)
                      ratioBySample <- plyr::rbind.fill(ratioBySample, row)
                  }
              }

              # generate plot
              p <- getBothPlot(ratioBySample, labels)
              return (p)
          })


#' get counts plot
#' helper function to generate a bar graph of sizes/complexities for each sample
#'
#' @param countsTable table of counts by sample and size/complexity
#' @param labels label plot by sample_name or alias
#' @param pos fill or stack
#' @param data data to plot (size/complexity)
#' @return ggplot object of plot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text

getCountsPlot <- function(countsTable, labels, pos, data) {
    Var1 <- NULL
    Freq <- NULL
    p <- ggplot(countsTable) +
        geom_bar(aes(
            fill = factor(Var1),
            y = Freq,
            x = switch(labels, "alias" = alias, "sample" = sample)
        ),
        position = pos,
        stat = "identity") +
        scale_fill_discrete(name = data) +
        labs(x = "Sample", y = "Counts") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    return(p)
}

#' get both plot
#' helper function to generate a bar graph of sizes/complexities for each sample
#'
#' @param countsTable table of counts by sample and size/complexity
#' @param labels label plot by sample_name or alias
#' @return ggplot object of plot

getBothPlot <- function(countsTable, labels) {
    size <- NULL
    dist_orig <- NULL
    Freq <- NULL

    p <- ggplot(countsTable) +
        facet_wrap(~ switch(labels, "alias" = alias, "sample" = sample)) +
        geom_point(aes(
            factor(size),
            factor(dist_orig),
            size = Freq,
            color = Freq
        )) +
        scale_fill_discrete(name = "complexity") +
        theme_minimal() + xlab("size") + ylab("complexity") +
        scale_size_area() + scale_color_gradient(low = "blue", high = "red")

    return (p)
}
