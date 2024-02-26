#' readstats plot (v3)
#'
#' Function to create size, complexity, both or ratio plot
#'
#' @param lox loxcode_experiment object
#' @param count_matrix
#' @param code_set
#' @param plot specify which plot (size/complexity/both/ratio)
#' @param fill specify whether reads should be normalized in boxplots
#' @param labels specify if graph should be labelled by sample name or alias
#' @return Plot based on the type of plot selected
#' @export

setGeneric("readstats_plot",
           function(lox,
                    count_matrix = "all_samples",
                    code_set = "all_codes",
                    plot = "size",
                    fill = TRUE,
                    labels = "alias") {
             standardGeneric("readstats_plot")
           })

setMethod("readstats_plot",
          "loxcode_experiment",
          function(lox,
                   count_matrix = "all_samples",
                   code_set = "all_codes",
                   plot = "size",
                   fill = TRUE,
                   labels = "sample") {

            # initialize variables
            counts = lox@count_matrixes[[count_matrix]]
            codes = lox@code_sets[[code_set]]
            pos = ifelse(fill, "fill", "stack")

            # generate plot
            if (plot == "ratio") {
              p = get_ratio_plot(lox, counts, codes, count_matrix, labels)
            } else {
              p = stats_bar_plot(lox, counts, codes, count_matrix, plot, pos, labels)
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
#' @param labels specify if samples should be labelled by sample_name or alias

stats_bar_plot <- function(lox, counts, codes, count_matrix, plot, pos, labels) {

  # initiliaze variables
  data = switch(plot,
                "size" = "size",
                "complexity" = "dist_orig",
                "both" = c("size", "dist_orig"))
  samples = names(counts)

  # create counts table
  countsBySample = getCountsTable(lox, samples, codes, count_matrix, counts, data)

  # create plot
  if (plot %in% c("size", "complexity")) {
    p = getCountsPlot(countsBySample, labels, pos, data)
  } else {
    p = getBothPlot(countsBySample, labels)
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
#' @param data data to count (size/complexity)

getCountsTable <- function(lox, samples, codes, count_matrix, counts, data) {

  countsBySample = data.frame()
  for (i in 1:length(samples)) {
    sampleName = samples[i]
    sampleAlias = get_alias(lox, count_matrix, sampleName)
    indices = codes$code%in%row.names(counts)[counts[[sampleName]]>0]
    countsTable = data.frame(table(codes[data][indices,]))
    if (nrow(countsTable) != 0) {
      row = data.frame("sample" = sampleName,
                       "alias" = sampleAlias,
                       countsTable)
      countsBySample = plyr::rbind.fill(countsBySample, row)
    }
  }
  return (countsBySample)
}

#' get ratio plot
#' helper function to generate plot of ratios
#'
#' @param lox loxcode_experiment
#' @param samples list of samples by name
#' @param codes code_set selected
#' @param count_matrix name of count_matrix selected
#' @param counts count_matrix selected

setGeneric("get_ratio_plot",
           function(lox,
                    counts,
                    codes,
                    count_matrix,
                    labels) {
             standardGeneric("get_ratio_plot")
           })

setMethod("get_ratio_plot",
          "loxcode_experiment",
          function(lox,
                   counts,
                   codes,
                   count_matrix,
                   labels) {

            # initialize
            samples = names(counts)
            data = c("size", "dist_orig")

            # get table of ratios by size and complexity
            ratioBySample = data.frame()
            for (i in 1:length(samples)) {
              sampleName = samples[i]
              sampleAlias = get_alias(lox, count_matrix, sampleName)
              codes = lox@samples[[sampleName]]@decode@data
              indices = codes$code%in%row.names(counts)[counts[[sampleName]]>0]
              subsetCodes = codes[indices,]
              countsBySample = data.frame(table(subsetCodes[data]))
              if (nrow(countsBySample) != 0) {
                totalsBySample = aggregate(
                  subsetCodes$count,
                  by=list(Category=subsetCodes$size, subsetCodes$dist_orig),
                  FUN=sum)
                names(totalsBySample) = c("size", "dist_orig", "total")
                ratioTable = merge(countsBySample, totalsBySample, by=data)
                ratioTable$Freq = ratioTable$Freq / ratioTable$total
                row = data.frame("sample" = sampleName,
                                 "alias" = sampleAlias,
                                 ratioTable)
                ratioBySample = plyr::rbind.fill(ratioBySample, row)
              }
            }

            # generate plot
            p = getBothPlot(ratioBySample, labels)
            return (p)
          })


#' get counts plot
#' helper function to generate a bar graph of sizes/complexities for each sample
#'
#' @param countsTable table of counts by sample and size/complexity
#' @param labels label plot by sample_name or alias
#' @param pos fill or stack
#' @param data data to plot (size/complexity)

getCountsPlot <- function(countsTable, labels, pos, data) {

  p = ggplot(countsTable) +
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
#' @param plot plot type (both/ratio)

getBothPlot <- function(countsTable, labels) {

  p = ggplot(countsTable) +
    facet_wrap( ~ switch(labels, "alias" = alias, "sample" = sample)) +
    geom_point(aes(factor(size), factor(dist_orig), size = Freq, color = Freq)) +
    scale_fill_discrete(name = "complexity") +
    theme_minimal() + xlab("size") + ylab("complexity") +
    scale_size_area() + scale_color_gradient(low = "blue", high = "red")

  return (p)
}
