library('plyr')
library('base')
#' Lox barcode casette object
#'
#' @slot lox_sites vector of cre-binding sites
#' @slot code_elements vector of distinguishable code elements

lox_casette <- setClass(
  "lox_casette",

  representation(
    lox_sites = "vector",
    code_elements = "vector"
  ),

  prototype = list(
    lox_sites = c(1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1),
    code_elements = c(1:13)
  )
)

#' Perform an excision on the sequence
#'
#' @param lox lox_casette object
#' @param a lox_site position 1
#' @param b lox_site position 2
#' @return an excised lox_casette object
#' @export

setGeneric("excise", function(lox, a, b) {standardGeneric("excise")})

setMethod("excise", "lox_casette", function(lox, a, b) {
  lox_sites = lox@lox_sites
  code_elements = lox@code_elements

  site_start = lox_sites[1:a]
  site_end = c()
  if (b < length(lox_sites)) site_end = lox_sites[(b+1):length(lox_sites)]
  lox@lox_sites = c(site_start, site_end)

  code_start = c()
  if (a > 1) code_start = code_elements[1:(a-1)]
  code_end = c()
  if (b <= length(code_elements)) code_end = code_elements[b:length(code_elements)]
  lox@code_elements = c(code_start, code_end)

  return(lox)
})

#' Perform an inversion on the sequence
#'
#' @param lox vector of lox sites
#' @param a site position 1
#' @param b site position 2
#' @return a new vector of lox sites
#' @export

setGeneric("invert", function(lox, a, b) {standardGeneric("invert")})

setMethod("invert", "lox_casette", function(lox, a, b) {
  lox_sites = lox@lox_sites
  code_elements = lox@code_elements

  code_start = c()
  if (a > 1) code_start = code_elements[1:(a-1)]
  code_mid = rev(code_elements[a:(b-1)])
  code_end = c()
  if (b <= length(code_elements)) code_end = code_elements[b:length(code_elements)]
  for (i in 1:length(code_mid)) {
    code_mid[i] = -code_mid[i]
  }
  lox@code_elements = c(code_start, code_mid, code_end)

  return(lox)
})

#' Simulate the lox recombination
#'
#' @param lox lox_casette object
#' @param site1 of casette
#' @param site2 of casette
#' @return lox_casette after random recombination
#' @export

setGeneric("simulate", function(lox, site1=NULL, site2=NULL) {standardGeneric("simulate")})

setMethod("simulate", "lox_casette", function(lox, site1=NULL, site2=NULL){
  lox_sites = lox@lox_sites
  code_elements = lox@code_elements
  size = length(lox_sites)

  # randomly choose two lox sites
  if (is.null(site1) || is.null(site2)) {
    sites = sort(base::sample(2:(size-1), 2, replace=FALSE))
    site1 = sites[[1]]
    site2 = sites[[2]]
  }

  # simulate
  if (sign(lox_sites[[site1]]) == sign(lox_sites[[site2]])) {
    lox = excise(lox, site1, site2)
  }
  else {
    lox = invert(lox, site1, site2)
  }

  return(lox)
})

#' Generate all recombination events
#'
#' @param lox lox_casette object
#' @return all recombination events
#' @export
setGeneric("simulate_all_init", function(lox) {standardGeneric("simulate_all_init")})

setMethod("simulate_all_init", "lox_casette", function(lox) {
  recombinations = c()
  for (site1 in 1:13) {
    for (site2 in (site1+1):14) {
      recombinations = c(recombinations, simulate(lox, site1, site2))
    }
  }
  return (recombinations)
})

#' Simulate n lox recombination events
#'
#' @param lox vector of lox sites
#' @param n number of recombination events to run
#' @return vector of lox sites after n recombinations
#' @export

setGeneric("simulate_runs", function(lox, n) {standardGeneric("simulate_runs")})

setMethod("simulate_runs", "lox_casette", function(lox, n) {
  for (i in 1:n) {
    lox = simulate(lox)
  }
  return(lox)
})

#' Simulate recombination of multiple casettes
#'
#' @param lox original casette
#' @param c number of casettes
#' @param n number of recombination events per casette
#' @return list of casettes after cre lox
#' @export

setGeneric("simulate_reads", function(lox, c, n) {standardGeneric("simulate_reads")})

setMethod("simulate_reads", "lox_casette", function(lox, c, n) {
  casettes = list(simulate_runs(lox, n))
  for (i in 1:(c-1)) {
    new_lox = simulate_runs(lox, n)
    casettes = list.append(new_lox, casettes)
  }
  return(casettes)
})

#' Simulate a loxcode_sample
#'
#' @param lox original casette
#' @param c number of casettes/ number of cells in which a barcode is made
#' @param n mean number of recombination events per casette
#' @param ref reference sample from which distribution of dist_orig is obtained
#' @param name sample_name
#' @return a loxcode_sample with simulated data
#' @export

setGeneric("simulate_sample", function(lox=NULL, c=3000, n=3, ref=NULL, name="simulation") {standardGeneric("simulate_sample")})

setMethod("simulate_sample", "lox_casette", function(lox=NULL, c=3000, n=3, ref = NULL, name="simulation") {

  if (is.null(lox)) lox = new("lox_casette")
  sample = new("loxcode_sample")
  decode = new("decode_output")

  # decode the simulated casettes
  if (!is.null(ref)) {
    dist_orig_dist = ref@decode@data$dist_orig[!is.na(ref@decode@data$dist_orig)]
    dist_orig_dist = getData(ref, "dist_orig")$dist_orig
    dist_orig_sample = base::sample(dist_orig_dist, c, replace=TRUE)
  } else {
    if (is.null(n)) n = 3
    dist_orig_sample = stats::rpois(c, n)
  }

  casettes = lapply(dist_orig_sample, function(x) simulate_runs(lox, x))
  data = data.frame()

  # unique_casettes = unique(casettes)
  # data = data.frame(matrix(ncol = 7, nrow = length(unique_casettes)), stringsAsFactors = FALSE)
  # code_elements = lapply(unique_casettes, function(x) x@code_elements)
  # names(data) = c("count", "firstread", "code", "size", "is_valid", "id", "dist_orig")
  # data$code = lapply(code_elements, function(x) paste0(x, collapse=" "))
  # data$size = lapply(code_elements, function(x) length(x))
  # data$is_valid = lapply(data$code, is_valid)
  # data$id = pack(unlist(data$code), unlist(data$is_valid))
  # data$dist_orig = dist_orig_sample
  # data$dist_orig[!unlist(data$is_valid)] = NA

  for (i in 1:length(casettes)) {
    code = casettes[[i]]@code_elements
    codeAsString = paste0(code, collapse=" ")
    size = length(code)
    is_valid = is_valid(codeAsString)
    id = pack(codeAsString, is_valid)
    dist_orig = dist_orig_sample[i]
    if (!is_valid) dist_orig = NA

    if (!(codeAsString %in% data$code)) {
      row = data.frame("count" = 1,
                       "firstread" = i,
                       "code" = codeAsString,
                       "size" = size,
                       "is_valid" = is_valid,
                       "id" = id,
                       "dist_orig" = dist_orig,
                       stringsAsFactors = FALSE)
      data = plyr::rbind.fill(data, row)
    }
    else {
      index = match(codeAsString, data$code)
      data$count[[index]] = data$count[[index]]+1
    }
  }

  decode@data = data
  decode@saturation = integer(0)
  sample@decode = decode

  sample@name = name
  sample@meta = data.frame()
  sample@decode_stats = list(tot_reads=c,
                             missing_start=0,
                             multi_start=0,
                             missing_end=0,
                             multi_end=0,
                             consensus_filtered=0)
  return(sample)
})

#' Simulate Loxcode Samples
#'
#' @param lox original casette
#' @param nsamples number of loxcode_samples
#' @param ncodes number of casettes/ number of cells in which a barcode is made
#' @param dist_type model from existing sample or poisson distribution
#' @param npois mean number of recombination events per casette
#' @param ref reference sample from which distribution of dist_orig is obtained
#' @param name sample_name
#' @return a list of loxcode_sample with simulated data
#' @export

setGeneric("simulate_nsamples", function(lox, nsamples=10, ncodes=3000, dist_type="Poisson", npois=NULL, ref=NULL, name="simulation") {standardGeneric("simulate_nsamples")})

setMethod("simulate_nsamples", "lox_casette", function(lox, nsamples=10, ncodes=300, dist_type="Poisson", npois=NULL, ref=NULL, name="simulation") {
  samples_list = c()
  samples_list = lapply(c(1:nsamples), function(x) c(samples_list, simulate_sample(lox, c=ncodes, n=npois, ref=ref, name=paste0(name, "_", x))))
  names(samples_list) = paste0(name, "_", c(1:nsamples))
  return(samples_list)
})
