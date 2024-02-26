#########0#########1#########2#########3#########4#########5#########6#########7
#/Users/zoeeeeeewu/Desktop/lox/
# Import libraries
library(shiny)
library(plotly)
library(DT)
library(loxcoder)
library(shinydashboard)
library(rlist)
library(shinyFiles)
library(shinyalert)
library(knitr)
library(rmarkdown)
library(yaml)
library(gridExtra)
library(testit)
library(shinyjs)
library(stringr)
library(tinytex)
library(pals)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(flexdashboard)
#########0#########1#########2#########3#########4#########5#########6#########7

# Source modules
source('../R/plots.R')
#/Users/zoeeeeeewu/Desktop/lox/LoxCodeR2022-master/LoxcodeR_app

#########0#########1#########2#########3#########4#########5#########6#########7
# upload size
options(shiny.maxRequestSize=50*1024^2)
#########0#########1#########2#########3#########4#########5#########6#########7

# Initialize dist maps
load_origin_distmaps('/media/arnav/Storage/origin/')
# load_pair_distmaps('/wehisan/general/user_managed/grpu_naik.s_2/TW/maps/origin/')

#########0#########1#########2#########3#########4#########5#########6#########7
# Loxcode experiment object

lox <- readRDS("Week2.rds")

#lox <- readRDS("test_data.rds")

#########0#########1#########2#########3#########4#########5#########6#########7
exp = list(lox)
d <- summary_table(lox, "all_samples")

chart_choices = c("Statistics Plots", "Heatmap", "Saturation Plot",
                  "Pair Comparison Plot")
codeset_selectionID = c("codeset_stats", "view_codeset",
                        "codeset_stats", "codeset_heat", "codeset_sat",
                        "filter_code_name", "filter_codeset", "codeset_pair")
sample_selectionID = c("view_sample", "matrix_stats", "matrix_heat",
                       "filter_name", "independent_samples", "sampleset_sat",
                       "sampleset_pair", "sim_sampleset")
samplesID = c("size_sample", "complexity_sample", "sample_sat", "sim_sample")

### ACTIVITY LOG
startSession = "Session started."
refreshSession = "Session refreshed."

### VARIABLES
react <- reactiveValues(curr=lox, name=lox@name, exp=exp, samp=d, curr_pair=NULL, pairs=list())
global <- reactiveValues()
params <- reactiveValues(functions=list(), types=list(),
                         inputs=list(), annotations=list(),
                         #notes = list(),
                         loxcodes=list()
)
########nSure
#plots = list())
##############
logs <- reactiveValues(activity=list(startSession), timestamps=list(paste(Sys.time())))
sat <- reactiveValues(samples=list(), codesets=list(), plots=list(), curr_plot = NULL)
sims <- reactiveValues(samples=list())
simulationValues = reactiveValues(samples=list(), sim_params=list())



### UPDATE FUNCTIONS
#' Update codeset based on selection
#'
#'\code{updateCodesetSelection} updates the codeset based on user selection in the app
#'
#' @param session Session object of the shiny server
#' @param selectionID The lost of IDs of all the input objects in codesets
#' @param selected A string with name of the selected codeset
#'
#' @return A session object of selected value of codeset by client
#' @export
#'
#' @examples
#' updateCodesetSelection(session,c("codeset_stats","view_codeset","codeset_stats"),"test")
#' Runs inside shiny app
updateCodesetSelection <- function(session, selectionID, selected) {
  print(session)
  print(selectionID)
  print(selected)
  for (ID in selectionID){
    updateSelectInput(session, ID, choices=names(react$curr@code_sets), selected=selected)
  }
}

#' Update Sample set based on selection
#'
#' @param session Session object of the shiny server
#' @param selectionID The lost of IDs of all the input objects in codesets
#' @param selected A string with name of the selected codeset
#'
#' @return A session object of selected value of codeset by client
#' @export
#'
#' @examples
#' updateCodesetSelection(session,c("view_sample","matrix_stats","matrix_heat"),"test")
#' Runs inside shiny app
updateSampleSelection <- function(session, selectionID, selected) {
  print(session)
  print(selectionID)
  print(selected)
  for (ID in selectionID){
    updateSelectInput(session, ID, choices=names(react$curr@count_matrixes), selected=selected)
  }
}

updateSamples <- function(session) {
  choices = names(react$curr@samples)
  updateSelectInput(session, "sample1_pair", "Sample 1:", choices=choices)
  updateSelectInput(session, "sample2_pair", "Sample 2:", choices=choices)
  for (ID in samplesID) {
    updateSelectInput(session, ID, "", choices=choices)
  }
}

updateCurrentExp <- function(session, curr, exp) {
  index = match(curr@name, loxcoder::exp_table(exp)$Experiment_Name)
  exp = list.remove(exp, index)
  exp = list.append(exp, curr)
  return(exp)
}

validateFastq <- function(session, samplesheet, files) {
  files = sort(files[grepl(".fastq$", files)])
  R1 = sort(files[grepl("R1_001.", files)])
  R2 = sort(files[grepl("R2_001.", files)])

  if ("sample" %in% names(samplesheet)){
    sample_names = sort(samplesheet$sample)
  } else {
    showNotification("Sample sheet is invalid. Missing `sample` column.")
    return(FALSE)
  }

  # find the fastq files for each sample
  for (i in sample_names) {
    if (sum(grepl(i,R1))!=1 || sum(grepl(i,R2))!=1) {
      showNotification(paste("Could not find *.fastq file for sample", i, "in fastq directory."))
      return(FALSE)
    }
  }

  # checks if there are two runs each in fastq directory
  if (length(R1) != length(R2)) {
    for (s in R1){
      if ((gsub("_R1_001", "_R2_001", s) %in% R2) == FALSE) {
        showNotification(paste(s, "is missing `R2_001` run in fastq directory."))
      }
    }
    for (s in R2){
      if ((gsub("_R2_001", "_R1_001", s) %in% R1) == FALSE) {
        showNotification(paste(s, "is missing `R1_001` run in fastq directory."))
      }
    }
  }

  return (TRUE)
}

# Function to call in place of dropdownMenu
dropdownMenuCustom <- function (..., type = c("messages", "notifications", "tasks"),
                                badgeStatus, icon = NULL, .list = NULL, customSentence = customSentence)
{
  type <- match.arg(type)
  if (!is.null(badgeStatus)) shinydashboard:::validateStatus(badgeStatus)
  items <- c(list(...), .list)
  lapply(items, shinydashboard:::tagAssert, type = "li")
  dropdownClass <- paste0("dropdown ", type, "-menu")
  if (is.null(badgeStatus)) {badge <- NULL}
  else {badge <- span(class = paste0("label label-", badgeStatus), numItems)}
  tags$li(
    class = dropdownClass,
    a(
      href = "#",
      class = "dropdown-toggle",
      `data-toggle` = "dropdown",
      icon,
      badge
    ),
    tags$ul(
      class = "dropdown-menu",
      tags$li(
        class = "header",
        customSentence(numItems, type)
      ),
      tags$li(
        tags$ul(class = "menu", items)
      )
    )
  )
}

customSentence <- function(numItems, type) {
  paste("Current Loxcode Experiment")
}

shinyalertAnnotate <- function(session, callbackR) { # shinymodel
  shinyalert(
    "Plot Added",
    text = "To add comments: \n Go to the 'Report section' and dobble click the 'Annotation' cell",
    type = "success"
  )
  callbackR()
}

shinyalertDescribe <- function(session, type, callbackR) {
  shinyalert(
    title = "Describe",
    text = paste("Write a short description that describes how you filtered the", switch(type, "sample"="sample", "code"="code"), "set: "),
    closeOnClickOutside = TRUE,
    showConfirmButton = TRUE,
    closeOnEsc = TRUE,
    showCancelButton = TRUE,
    type = "input",
    inputType = "text",
    callbackR = callbackR
  )}

shinyalertName <- function(session, callbackR) {
  shinyalert(
    title = "Name",
    text = "Name the merged experiment: ",
    type = "input",
    inputType = "text", callbackR = callbackR
  )
}

### INCLUDE IN REPORT
include <- function(value, plot, type, input) {
  params$functions = list.append(params$functions, plot)
  params$types = list.append(params$types, type)
  params$inputs = list.append(params$inputs, input)
  params$annotations = list.append(params$annotations, value)
  #params$notes = list.append(params$notes, comment)

  params$loxcodes = list.append(params$loxcodes, react$curr)

  ########nSure
  #params$plots = list.append(params$plots, graph)
  ##################
  addToLog(session, logReport(session, react$curr, type, input))
}

### ACTIVITY LOG FUNCTIONS
# Add item to activity log
addToLog <- function(session, item) {
  logs$activity = list.append(logs$activity, item)
  logs$timestamps = list.append(logs$timestamps, paste(Sys.time()))
}

# log import

#' Upload the logs
#'
#' @param session
#' @param lox
#' @param method
#'
#' @return
#' @export
#'
#' @examples
logUpload <- function(session, lox, method) {
  item = paste("Uploaded loxcode_experiment ", lox@name, "by",
               switch(method, "rds"="'RDS object upload'.", "fastq"="'fastq files and directory upload'."))
  return (item)
}

# log merge experiments
logMerge <- function(session, experiments, lox) {
  names = ""
  for (i in 1:length(experiments)) {
    if (i != (length(experiments))) { names = paste0(names, experiments[[i]]@name, sep = ", ") }
    else { names = paste0(names, experiments[[i]]@name, sep = " ")}
  }
  item = paste0("Merged Loxcode experiments: ", names, ". Created: ", lox@name, ". ")
  return (item)
}

# log created sample sets or code sets
logCreate <- function(session, lox, set_name, type, method, description="") {
  item = paste("Created new ", type, " set in ", lox@name, " (", set_name, ") by ",
               switch(method, "selection"="'Create from Selection'.", "all"="'Create from All'."), description)
  return (item)
}

# log rename or delete sample sets or code sets
logUpdate <- function(session, lox, set_name, type, method, new_name=NULL) {
  item = paste(switch(method, "rename"="Renamed", 'delete'="Deleted"),
               type, " set in ", lox@name, " (", set_name, ") ",
               switch(method, "rename"=paste("to '", new_name, "'"), "delete"=""), ".")
  return (item)
}

# log collapse samples
logCollapse <- function(session, lox, new_set, union, average, params) {
  type <- function() {
    if (union & average) { type = paste("(union and average)") }
    else if (union & !average) { type = paste("(union and sum)") }
    else if (!union & average) { type = paste("(intersection and average)")}
    else { type = paste("(intersection and sum)")}
    return (type)
  }
  item = paste("Collapsed", type(), "of samples in", lox@name, "on parameters:",
               paste(params, collapse=", "), ". Created sample set '", new_set, "'.")
  return (item)
}

# log filter codes
logFilter <- function(session, lox, sample, code, max_reps, tolerance, new_name) {
  item = paste("Filtered codes in", lox@name, "on parameters: Sample_set =", sample, ", Code_set =", code,
               ", Max_reps = ", max_reps, ", Tolerance_level =", tolerance, ". Created code set '", new_name, "'.")
  return (item)
}

# log add to report
logReport <- function(session, lox, plot, parameters) {
  item = paste("Added plot to report.", str_to_title(plot), "of", lox@name)
  if (plot!="pair_comparison_plot") {
    item = paste(item, "with parameters: ", paramsAsText(parameters))
  }
  item = paste(item, ".")
  return (item)
}

# log download report
logDownloadReport <- function(session, file, type) {
  item = paste("Downloaded", type, "file.", file, sep=" ")
  return (item)
}

### converts the parameters into text
paramsAsText <- function(params) {
  parametersAsText = list()
  for (i in 1:length(params)) {
    n = names(params)[[i]]
    p = params[[i]]
    if (is(p, "loxcode_experiment") | is(p, "loxcode_sample")) {
      parametersAsText = list.append(parametersAsText, paste(n, "=", p@name))
    }
    else if (is.character(p)) {
      parametersAsText = list.append(parametersAsText, paste(n, "=", p))
    }
    else if (is.numeric(p)) {
      parametersAsText = list.append(parametersAsText, paste(n, "=", paste(p, ",")))
    }
    else if (rapportools::is.boolean(p)) {
      parametersAsText = list.append(parametersAsText, paste(n, "=", p))
    }
  }
  return(paste(parametersAsText, collapse=", "))
}

function(input, output, session) {

  # current loxcode_experiment object
  output$curr_lox = renderMenu({
    dropdownMenuCustom(
      type = "messages",
      icon = icon("bookmark"),
      badgeStatus = NULL,
      customSentence = customSentence,
      messageItem(from=react$curr@name, message="", icon=icon("dna"), href=NULL)
    )})


  ### IMPORT
  # upload loxcode_experiment object
  observeEvent(
    input$submit_rds, {
      if (is.null(input$rds_file)){
        showNotification("Please specify a file path.")
        return
      } else {
        if (grepl(".rds$", input$rds_file[[length(input$rds_file)]])){

          obj = readRDS(file=input$rds_file$datapath)
          if (is(obj, "loxcode_experiment")){
            react$curr <- obj

            # change double to int
            #View(Week2@code_sets$all_codes$size)
            if (!is.null(react$curr@code_sets$all_codes$size)){
              as.integer(react$curr@code_sets$all_codes$size)
            }


            react$samp <- sample_table(react$curr, "all_samples")
            react$exp <- rlist::list.append(react$exp, react$curr)
            updateSampleSelection(session, sample_selectionID, NULL)
            updateCodesetSelection(session, codeset_selectionID, NULL)
            updateSamples(session)
            addToLog(session, logUpload(session, react$curr, "rds"))
          } else {
            showModal(modalDialog(
              title = "Oops!",
              "Object uploaded was not of class <loxcode_experiment>",
              footer = modalButton("Ok")
            ))
          }
        } else {
          plotOutput("distPlot")

          showModal(modalDialog(
            title = "Oops!",
            "File uploaded was not an R object (*.rds).",
            footer = modalButton("Ok")
          ))
        }
      }
    }
  )

  observeEvent(
    input$submit_fastq, {
      files = list.files(input$dir_input)
      print("10")
      if (is.null(input$samplesheet)){
        showNotification("Please specify a file path.")
      } else {
        if (grepl(".xls$", input$samplesheet$datapath) | grepl(".xlsx$", input$samplesheet$datapath)){ # validate file extension
          samplesheet = read_excel(input$samplesheet$datapath)
          print(input$samplesheet$datapath)


          if (validateFastq(session, samplesheet, files)) { # validate file contents
            newlox <- load_from_xlsx(
              name = input$name_exp,
              s=input$samplesheet$datapath,
              dir=input$dir_input,
              load = TRUE,
              full = FALSE)
            print("11")
            react$curr <- newlox
            react$samp <- sample_table(react$curr, "all_samples")
            react$exp <- rlist::list.append(react$exp, react$curr)
            updateSampleSelection(session, sample_selectionID, NULL)
            updateCodesetSelection(session, codeset_selectionID, NULL)
            updateSamples(session)
            addToLog(session, logUpload(session, react$curr, "fastq"))
          } else {
            print("12")
            showModal(modalDialog(
              title = "Oops!",
              "Invalid files uploaded.",
              footer = modalButton("Ok")
            ))

          }

        } else {
          print("13")
          showModal(modalDialog(
            title = "Oops!",
            "File uploaded was not an excel file (*.xls or *.xlsx).",
            footer = modalButton("Ok")
          ))

        }
      }
    }
  )

  #table of loxcode_experiment objects
  output$experiments_table = renderDataTable({datatable(
    exp_table(react$exp),
    rownames = FALSE,
    class = "cell-border stripe",
    filter = 'top',
    selection = 'multiple'
  )})

  observeEvent(
    input$select_exp, {
      if (is.null(input$experiments_table_rows_selected) |
          length(input$experiments_table_rows_selected) > 1) {

        showModal(modalDialog(
          title = "Oops!",
          "Please select one experiment!",
          footer = modalButton("Ok")
        ))
        return ()
      }
      else if (length(input$experiments_table_rows_selected) == 1) {
        react$curr = react$exp[[input$experiments_table_rows_selected]]
        showNotification(paste(react$curr@name, " selected."))
      }
    }
  )

  observeEvent(
    input$merge_exp, {
      if (length(input$experiments_table_rows_selected) != 2) {

        showModal(modalDialog(
          title = "Oops!",
          "Please select two experiments to merge!",
          footer = modalButton("Ok")
        ))

        return()
      }
      else {
        shinyalertName(session, mergeExperiments)
      }
    })

  mergeExperiments <- function(value) {
    index = input$experiments_table_rows_selected[1:2]
    experiments = react$exp[index]
    showNotification("Merging experiments...")
    react$curr = merge_experiments(experiments[[1]], experiments[[2]], name = value)
    showNotification("Experiments merged!")
    react$samp <- sample_table(react$curr, "all_samples")
    react$exp <- rlist::list.append(react$exp, react$curr)
    updateSampleSelection(session, sample_selectionID, NULL)
    updateCodesetSelection(session, codeset_selectionID, NULL)
    updateSamples(session)
    addToLog(session, logMerge(session, experiments, react$curr))
  }

  # samplesheet view
  output$samplesheet = renderDataTable({
    d = data.frame()
    if (!is.null(input$samplesheet)) {
      d = read_excel(input$samplesheet$datapath)
      d$Status = ""
    }
    else
    {
      d = react$curr@samp_table
    }
    datatable(d,options = list(
      # dom = 't',
      scrollX = TRUE,
      # scrollY = TRUE,
      # scroller = TRUE,
      fixedColumns = list(leftColumns = 2)
    ))
  })

  output$save_exp = downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".rds", sep="")
    },
    content = function(file) {
      saveRDS(react$curr, file)
    }
  )

  observeEvent(
    input$del_exp, {
      if (!is.null(input$experiments_table_rows_selected)){
        react$exp = list.remove(react$exp, input$experiments_table_rows_selected)
      }
    }
  )

  ### Research details

  observeEvent(
    input$detail_submit,{
      rmd_content=""
      print("***")
      print(input$name)
      print(input$comment)
      print("***")

      rmd_title=paste("### Researscher Name :","",sep="")
      rmd_content=paste(rmd_content,rmd_title,sep = "\n")
      for(ii in 1:length(str_split(input$name,"\n", simplify = TRUE))){
        rmd_content=paste(rmd_content,str_split(input$name,"\n", simplify = TRUE)[[ii]],sep = "\n")
        rmd_content=paste(rmd_content,"",sep = "\n")
      }
      rmd_title=paste("### Experiment Name :","",sep="")
      rmd_content=paste(rmd_content,rmd_title,sep = "\n")
      for(ii in 1:length(str_split(input$expname,"\n", simplify = TRUE))){
        rmd_content=paste(rmd_content,str_split(input$expname,"\n", simplify = TRUE)[[ii]],sep = "\n")
        rmd_content=paste(rmd_content,"",sep = "\n")
      }
      rmd_title=paste("### Comments :","",sep="")
      rmd_content=paste(rmd_content,rmd_title,sep = "\n")
      for(ii in 1:length(str_split(input$comment,"\n", simplify = TRUE))){
        rmd_content=paste(rmd_content,str_split(input$comment,"\n", simplify = TRUE)[[ii]],sep = "\n")
        print(str_split(input$comment,"\n", simplify = TRUE))
        rmd_content=paste(rmd_content,"",sep = "\n")
      }
      writeLines(rmd_content, "./chapter5.Rmd")
    }
  )



  ### SIMULATIONS
  output$simulated_samples = renderDataTable({datatable(
    data.frame(simulated_samples = names(sims$samples)),
    options = list(dom = 't')
  )})

  observeEvent(
    input$simulate_sample, {
      if (input$dist_orig_type == "Poisson") {
        ref = NULL
        npois = input$poisson_mean
      } else {
        ref = react$curr@samples[[input$sim_sample]]
        npois = NULL
      }
      sims$samples = loxcoder::simulate_nsamples(
        lox = new("lox_casette"),
        nsamples = input$nsamples,
        ncodes = input$ncodes,
        name = input$sim_names,
        ref = ref,
        npois = npois
      )
    }
  )

  output$dist_orig_sim = renderPlotly({
    if (length(sims$samples) == 0 |
        is.null(input$simulated_samples_rows_selected)) { return() }
    else {
      index = input$simulated_samples_rows_selected
      index = index[length(index)]
      sample = sims$samples[[index]][[1]]
      dist_orig_plot_sample(sample)
    }
  })

  output$size_sim = renderPlotly({
    if (length(sims$samples) == 0 |
        is.null(input$simulated_samples_rows_selected))
    { return() }
    else {
      index = input$simulated_samples_rows_selected
      index = index[length(index)]
      sample = sims$samples[[index]][[1]]
      size_plot_sample(sample)
    }
  })

  output$dist_orig_density = renderPlotly({
    if (length(sims$samples) == 0 |
        is.null(input$simulated_samples_rows_selected))
    { return() }
    else {
      index = input$simulated_samples_rows_selected
      samples = sims$samples[index]
      density_plot(
        ref = react$curr@samples[[input$sim_sample]],
        npois = input$poisson_mean,
        dist_type = "sample",
        samples = samples,
        plot_type = "dist_orig"
      )
    }
  })

  output$size_density = renderPlotly({
    if (length(sims$samples) == 0 | is.null(input$simulated_samples_rows_selected)) { return() }
    else {
      index = input$simulated_samples_rows_selected
      samples = sims$samples[index]
      density_plot(
        ref = react$curr@samples[[input$sim_sample]],
        npois = input$poisson_mean,
        dist_type = "sample",
        samples = samples,
        plot_type = "size"
      )
    }
  })

  observe({
    if (is.null(react$curr)) { return() }
    else {
      aliases = react$curr@alias[[input$sim_sampleset]]
      if (is.null(aliases)) { return() }
      else {
        updateSelectInput(session, "sim_alias",
                          choices=react$curr@alias[[input$sim_sampleset]]$alias)
        updateSelectInput(session, "sim_sample",
                          choices=react$curr@alias[[input$sim_sampleset]]$sample_name)
      }
    }
  })

  # coordinate sample names and aliases
  # observe(
  #   if (is.null(react$curr)) { return() }
  #   else{
  #     aliases = react$curr@alias[[input$sim_sampleset]]
  #
  #     if (is.null(aliases) | input$sim_sample=="" | input$sim_alias=="") { return () }
  #
  #     else if (input$dist_orig_type=='Sample Alias') {
  #       selected_sample = get_samplename(react$curr, input$sim_sampleset, input$sim_alias)
  #       updateSelectInput(session, "sim_sample", "Samples:",
  #                         choices = aliases$sample_name,
  #                         selected = selected_sample)
  #     }
  #
  #     else if (input$dist_orig_type=='Sample Name') {
  #       selected_sample = get_alias(react$curr, input$sim_sampleset, input$sim_sample)
  #       updateSelectInput(session, "sim_alias", "Samples:",
  #                         choices = aliases$alias,
  #                         selected = selected_sample)
  #     }
  #   }
  # )

  ### CREATE CODESET
  # d <- summary_table(react$curr, input$view_sample)
  # react$samp <- d

  output$codeset_table = renderDataTable({

    d=codeset_table(react$curr, input$view_codeset)
    if (!is.null(d$Size)) {
      d$Size = factor(d$Size)
      #           labels = c(unique(d$Size)))
    }


    datatable(
      d,
      rownames = FALSE,
      class = "cell-border stripe",
      filter = 'top',
      options = list(
        scrollX=T,
        fixedColumns=list(leftColumns=2)
      ),
      extensions = c('FixedColumns','Scroller')
    )} %>% formatStyle(columns=c(seq(2, ncol(react$samp))), 'text-align'='center'))

  output$selected_codes = renderPrint({
    s = input$codeset_table_rows_selected

    d <- codeset_table(react$curr, input$view_codeset)
    if (length(s)) {
      if (length(s)==1) { cat(length(s),'Code Selected:\n\n') }
      else { cat(length(s),'Codes Selected:\n\n') }
      cat(d$Code[s], sep = ', ')
    }
  })




  observeEvent(
    input$delete_codeset, {
      react$curr <- delete_codeset(react$curr, input$view_codeset)
      addToLog(session, logUpdate(session, react$curr, input$view_codeset, "code", "delete"))
      updateCodesetSelection(session, codeset_selectionID, input$name_codeset)
      react$exp = updateCurrentExp(session, react$curr, react$exp)
    }
  )

  observeEvent(
    input$create_codeset, {
      selection = input$codeset_table_rows_selected
      if (is.null(selection)) {
        showModal(modalDialog(
          title = "Oops!",
          "Please select at least one code set",
          footer = modalButton("Ok")
        ))
        return() }

      if (input$name_codeset == "") {
        showModal(modalDialog(
          title = "Oops!",
          "Please name the code set",
          footer = modalButton("Ok")
        ))
        return() }

      # print(input$name_codeset)


      # if (length(selection)) {
      react$curr = make_codeset_index(react$curr, c=input$view_codeset,
                                      I=selection, n=input$name_codeset)
      addToLog(session, logCreate(session, react$curr,
                                  input$name_codeset,
                                  "code", "selection"))
      updateCodesetSelection(session, codeset_selectionID, input$name_codeset)
      updateTextInput(session, "name_codeset", label="Name of new codeset:",
                      placeholder="Codeset Name", value="")
      react$exp = updateCurrentExp(session, react$curr, react$exp)
      #}
    }
  )

  observeEvent(
    input$create_all_codeset, {
      if (input$name_codeset == ""){
        showModal(modalDialog(
          title = "Oops!",
          "Please name the code set",
          footer = modalButton("Ok")
        ))
        return()
      }
      react$curr = make_codeset_index(react$curr,
                                      c=input$view_codeset,
                                      I=input$codeset_table_rows_all,
                                      n=input$name_codeset)
      global$name = input$name_codeset

      #shinyalertDescribe(session, "code", logCodeFilter)
      updateCodesetSelection(session, codeset_selectionID, input$name_codeset)
      updateTextInput(session, "name_codeset", label="Name of new codeset:", placeholder="Codeset Name", value="")
      react$exp = updateCurrentExp(session, react$curr, react$exp)
    }
  )

  logCodeFilter <- function(value) {
    addToLog(session, logCreate(session, react$curr, global$name, "code", "all", value))
  }

  includeRatio <- function(value) {
    include(value=value,
            plot=readstats_plot,
            type="ratio_plot",
            input=list(loxcode_experiment=react$curr,
                       count_matrix=input$matrix_stats,
                       code_set=input$codeset_stats, plot="ratio")
    )}

  observeEvent(
    input$rename_codeset, {
      react$curr = rename_codeset(react$curr, c=input$view_codeset, n=input$name_codeset)
      addToLog(session, logUpdate(session, react$curr, input$view_codeset,
                                  "code", "rename", input$name_codeset))
      updateCodesetSelection(session, codeset_selectionID, input$name_codeset)
      updateTextInput(session, "name_codeset", label="Name of new codeset:",
                      placeholder="Codeset Name", value="")
      react$exp = updateCurrentExp(session, react$curr, react$exp)
    }
  )

  ### CREATE SAMPLE SET
  output$sample_table = renderDataTable({
    d <- summary_table(react$curr, input$view_sample)
    react$samp <- d

    if (!is.null(d$Population)) {
      d$Population = factor(d$Population)
      #     labels = c(unique(d$Population)))
    }
    #
    #
    # if (!is.null(d$population)) {
    #   d$population = factor(d$population,
    #                         labels = c(unique(d$population)))
    # }
    #
    if (!is.null(d$Barcode.induction.age)) {
      d$Barcode.induction.age = factor(d$Barcode.induction.age)
      #             labels = c(unique(d$Barcode.induction.age)))
    }

    if (!is.null(d$timepoint)) {
      d$timepoint = factor(d$timepoint)
      #         labels = c(unique(d$timepoint)))
    }

    if (!is.null(d$experiment)) {
      d$experiment = factor(d$experiment)
      #          labels = c(unique(d$experiment)))
    }

    if (!is.null(d$Experiment)) {
      d$Experiment = factor(d$Experimentt)
      #          labels = c(unique(d$Experiment)))
    }

    if (!is.null(d$mouse)) {
      d$mouse = factor(d$mouse)
      #           labels = c(unique(d$mouse)))
    }

    if (!is.null(d$Mouse)) {
      d$Mouse = factor(d$Mouse)
      #        labels = c(unique(d$Mouse)))
    }

    if (!is.null(d$PCR.replicate)) {
      d$PCR.replicate = factor(d$PCR.replicate)
      #         labels = c(unique(d$PCR.replicate)))
    }
    #
    #     if (!is.null(d$population)) {
    #       d$population = factor(react$curr@samp_table$population,
    #                             labels = c(unique(react$curr@samp_table$population)))
    #     }
    #
    #     if (!is.null(d$Barcode.induction.age)) {
    #       d$Barcode.induction.age = factor(react$curr@samp_table$Barcode.induction.age,
    #                                        labels = c(unique(react$curr@samp_table$Barcode.induction.age)))
    #     }
    #
    #     if (!is.null(d$timepoint)) {
    #       d$timepoint = factor(react$curr@samp_table$timepoint,
    #                            labels = c(unique(react$curr@samp_table$timepoint)))
    #     }
    #
    #     if (!is.null(d$experiment)) {
    #       d$experiment = factor(react$curr@samp_table$experiment,
    #                            labels = c(unique(react$curr@samp_table$experiment)))
    #     }
    #
    #     if (!is.null(d$Experiment)) {
    #       d$Experiment = factor(react$curr@samp_table$Experiment,
    #                             labels = c(unique(react$curr@samp_table$Experiment)))
    #     }
    #
    #     if (!is.null(d$mouse)) {
    #       d$mouse = factor(react$curr@samp_table$mouse,
    #                             labels = c(unique(react$curr@samp_table$mouse)))
    #     }
    #
    #     if (!is.null(d$Mouse)) {
    #       d$Mouse = factor(react$curr@samp_table$mouse,
    #                        labels = c(unique(react$curr@samp_table$Mouse)))
    #     }
    datatable(
      d,
      filter = 'top',
      rownames = FALSE,
      class = "cell-border stripe",
      editable = list(target="cell", disable=list(columns=c(0, seq(2, ncol(react$samp))))),
      selection = list(target = 'row'),
      options = list(
        # dom = 't',
        scrollX = TRUE,
        # scrollY = TRUE,
        # scroller = TRUE,
        fixedColumns = list(leftColumns = 2)
      ),
      extensions = c('FixedColumns','Scroller')
    )} %>% formatStyle(columns=c(seq(3, ncol(react$samp))), 'text-align'='center'))

  observeEvent(
    input$view_sample, {
      d <- summary_table(react$curr, input$view_sample)
      react$samp <- d
    })

  output$selected_samples = renderPrint({
    s = input$sample_table_rows_selected
    d <- summary_table(react$curr, input$view_sample)
    if (length(s)) {
      if (length(s)==1) { cat(length(s),'Sample Selected:\n\n') }
      else { cat(length(s),'Samples Selected:\n\n') }
      cat(d$Sample_Name[s], sep = ', ')
    }
  })

  # renaming samples
  proxy = dataTableProxy("sample_table")
  observeEvent(
    input$sample_table_cell_edit, {
      d = react$samp
      info = input$sample_table_cell_edit
      i = info$row
      j = info$col + 1  # column index offset by 1
      v = info$value
      d[i, j] <<- coerceValue(v, d[i, j])
      sample = d[i,j-1]
      replaceData(proxy, d, resetPaging=FALSE, rownames=FALSE)
      react$curr = new_alias(react$curr, input$view_sample, sample, v)
      updateSelectInput(session, "sample1_pair", "Sample 1:",
                        choices=names(react$curr@samples),
                        selected=input$sample1_pair)
      updateSelectInput(session, "sample2_pair", "Sample 2:",
                        choices=names(react$curr@samples),
                        selected=input$sample1_pair)
      react$exp = updateCurrentExp(session, react$curr, react$exp)
    })

  # create sample set from selection
  observeEvent(
    input$create_sample, {
      selection = input$sample_table_rows_selected
      print(input$sample_table_rows_selected)

      #print(is.null(input$view_sample))
      # print(is.null(input$name_sample))
      # print(length(input$name_sample))
      # print((input$name_sample) == "")
      #
      # print(input$name_sample)


      # name_sample = input$name_sample
      if (is.null(selection)) {
        showModal(modalDialog(
          title = "Oops!",
          "Please select at least one sample",
          footer = modalButton("Ok")
        ))
        return() }

      if ((input$name_sample) == "") {
        showModal(modalDialog(
          title = "Oops!",
          "Please name the sample",
          footer = modalButton("Ok")
        ))
        return() }

      # if (length(selection)) {
      react$curr = make_count_matrix(react$curr,
                                     c=input$view_sample,
                                     i=selection, n=input$name_sample)



      addToLog(session, logCreate(session, react$curr, input$name_sample,
                                  "sample set", "selection"))
      updateSampleSelection(session, sample_selectionID, input$name_sample)
      updateTextInput(session, "name_sample",
                      label="Name of new collection of samples:",
                      placeholder="Sample Collection Name", value="")
      react$exp = updateCurrentExp(session, react$curr, react$exp)
      #}
    }
  )

  # create sample set from all
  observeEvent(
    input$create_all_sample, {
      if (input$name_sample == ""){
        showModal(modalDialog(
          title = "Oops!",
          "Please name the sample set",
          footer = modalButton("Ok")
        ))
        return()
      }
      react$curr = make_count_matrix(react$curr, c=input$view_sample,
                                     i=input$sample_table_rows_all,
                                     n=input$name_sample)
      global$name = input$name_sample
      #shinyalertDescribe(session, "sample", logSampleFilter)
      updateSampleSelection(session, sample_selectionID, input$name_sample)
      updateTextInput(session, "name_sample",
                      label="Name of new collection of samples:",
                      placeholder="Sample Collection Name", value="")
      react$exp = updateCurrentExp(session, react$curr, react$exp)
    }
  )

  logSampleFilter <- function(value) {
    addToLog(session, logCreate(session, react$curr, global$name, "sample set", "all", value))
  }

  # delete a sample set
  observeEvent(
    input$delete_sample, {
      react$curr <- delete_count_matrix(react$curr, input$view_sample)
      addToLog(session, logUpdate(session, react$curr, input$view_sample, "sample", "delete"))
      updateSampleSelection(session, sample_selectionID, "all_samples")
      react$exp = updateCurrentExp(session, react$curr, react$exp)
    }
  )

  # save data for console migration
  observeEvent(
    input$save_data,{
      saveRDS(react$curr, "temp-data.rds")
      saveRDS(ggplotly(),"lastplot.rds")
    }
  )
  # load data from console
  observeEvent(
    input$reload_data,{
      #react$curr <- readRDS("temp-data.rds")
      obj = readRDS("temp-data.rds")
      if (is(obj, "loxcode_experiment")){
        react$curr <- obj

        # change double to int
        #View(Week2@code_sets$all_codes$size)
        if (!is.null(react$curr@code_sets$all_codes$size)){
          as.integer(react$curr@code_sets$all_codes$size)
        }


        react$samp <- sample_table(react$curr, "all_samples")
        react$exp <- rlist::list.append(react$exp, react$curr)
        updateSampleSelection(session, sample_selectionID, NULL)
        updateCodesetSelection(session, codeset_selectionID, NULL)
        updateSamples(session)
      }

      print("loaded")
    }
  )

  # rename sample set
  observeEvent(
    input$rename_sample, {
      if (input$view_sample == "all_samples") { showNotification("Sample set 'all_samples' cannot be renamed.") }
      else {
        react$curr = rename_sampleset(react$curr, input$view_sample, input$name_sample)
        addToLog(session, logUpdate(session, react$curr, input$view_sample, "sample", "rename", input$name_sample))
        updateSampleSelection(session, sample_selectionID, input$name_sample)
        updateTextInput(session, "name_sample", label="Name of new collection of samples:", placeholder="Sample Collection Name", value="")
        react$exp = updateCurrentExp(session, react$curr, react$exp)
      }
    }
  )

  # generate aliases
  observeEvent(
    input$generate_alias,
    if (!length(input$alias_parameters)) { return() }
    else {
      react$curr = generate_alias(react$curr, input$view_sample, input$alias_parameters)
    }
  )

  observe({
    updateCheckboxGroupInput(
      session,
      "alias_parameters",
      "Choose alias parameters:",
      choices=names(get_collapsed_meta(react$curr, input$view_sample)))
  })

  # collapse samples
  observeEvent(
    input$collapse_samples, {
      print(input$collapse_parameters)
      if (!length(input$collapse_parameters)) { return() }
      else {
        react$curr <- collapse(react$curr, input$view_sample, input$collapse_parameters, input$collapse_name, input$collapse_union, input$collapse_average)
        addToLog(session, logCollapse(session, react$curr, input$collapse_name, input$collapse_union, input$collapse_average, input$collapse_parameters))
        updateSampleSelection(session, sample_selectionID, input$collapse_name)
        updateCheckboxInput(session, "collapse_union", "Union", value=NULL)
        updateCheckboxInput(session, "collapse_average", "Average", value=NULL)
        updateTextInput(session, "collapse_name", label="Name of new sample set:", placeholder="Sample Set Name", value="")
        updateCheckboxGroupInput(session, "collapse_parameters", "Choose parameters to collapse:", choices=names(get_collapsed_meta(react$curr, input$view_sample)), selected=NULL)
        react$exp = updateCurrentExp(session, react$curr, react$exp)
      }
    }
  )

  observeEvent(
    input$collapse_selection, {
      if (length(input$sample_table_rows_selected) < 1) { return () }
      else {
        react$curr = collapse_selection(lox=react$curr,
                                        count_matrix=input$view_sample,
                                        index=input$sample_table_rows_selected,
                                        union=input$collapse_union,
                                        average=input$collapse_average)
      }
    }
  )

  # observe({
  #   updateCheckboxGroupInput(
  #     session,
  #     "collapse_parameters",
  #     "Choose parameters to collapse:",
  #     choices=names(get_collapsed_meta(react$curr, input$view_sample)))
  # })

  ### FILTER CODES
  output$unfiltered_codes = renderPlot(
    code_frequency_pie(react$curr, input$independent_samples, input$filter_codeset)
  )

  output$filtered_codes = renderPlot(
    filtered_codes_pie(react$curr, input$independent_samples, input$filter_codeset, input$filter_tolerance, input$filter_reps)
  )

  observe({
    if (is.null(react$curr)) {
      showModal(modalDialog(
        title = "Oops!",
        "The dataset is empty",
        footer = modalButton("Ok")
      ))
      return() }


    else {
      updateSelectInput(session, "filter_codeset",
                        choices = names(react$curr@code_sets)[!names(react$curr@code_sets) == "invalid_codes"])
    }
  })

  # observe({
  #   if (is.null(react$curr@code_sets[[input$filter_codeset]])) { return () }
  #   else {
  #     Y = code_freq_table(react$curr, input$independent_samples, input$filter_codeset)
  #     total = max(as.numeric(names(Y[,!names(Y)%in%c("size", "dist_orig", "radius")])))
  #     updateSliderInput(session, "filter_reps", label="Maximum allowed code repetitions", min=2, max=total, value=input$filter_reps, step=1)
  #   }
  # })

  # filter codes
  observeEvent(
    input$create_filtered, {

      if (input$filter_code_name == "") {
        showModal(modalDialog(
          title = "Oops!",
          "Please name the filtered code set",
          footer = modalButton("Ok")
        ))
        return() }

      params = list(react$curr,
                    input$independent_samples,
                    input$filter_codeset, input$filter_tolerance,
                    input$filter_reps, input$filter_code_name)
      react$curr = do.call(make_filtered_codeset, params)
      addToLog(session, do.call(logFilter, list.append(session, params)))
      react$exp = updateCurrentExp(session, react$curr, react$exp)
      updateCodesetSelection(session, codeset_selectionID, input$filter_code_name)
      updateSampleSelection(session, sample_selectionID, input$independent_samples)
    }
  )

  ## STATISTICS PLOT
  # bar code table
  output$barcode_table = renderDataTable({
    temp_table = do.call(barcode_table, list(react$curr, count_matrix=input$matrix_stats,
                                             labels=input$labels_stats))
  })



  observeEvent(
    input$includeBartable, {
      shinyalertAnnotate(session, includeBartable)
    })

  includeBartable <- function(value="value") {
    include(value=value,
            plot=barcode_table,
            type="barcode_table",
            input=list(react$curr, count_matrix=input$matrix_stats,
                       labels=input$labels_stats)
            #comment = comment
    )}





  # Size plot
  output$size_plot = renderPlotly({
    do.call(readstats_plot, list(react$curr, count_matrix=input$matrix_stats,
                                 code_set=input$codeset_stats, plot="size",
                                 fill=input$fill_size, labels=input$labels_stats))

  })

  observeEvent(
    input$includeSize, {
      includeSize()
    }, ignoreInit = FALSE, ignoreNULL = FALSE)

  includeSize <- function(value="value") {
    include(value=value,
            plot=readstats_plot,
            type="size_plot",
            input=list(react$curr, count_matrix=input$matrix_stats,
                       code_set=input$view_codeset1, plot="size", fill=input$fill_size,
                       labels=input$labels_stats)
            #comment = comment
    )}

  # Complexity plot
  output$complexity_plot = renderPlotly({
    #graph <-
    readstats_plot(react$curr, count_matrix=input$matrix_stats,
                   code_set=input$codeset_stats, plot="complexity",
                   fill=input$fill_complexity,
                   labels=input$labels_stats)
    #+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5)))
  })

  observeEvent(
    input$includeComplexity, {
      includeComplexity()
    }, ignoreInit = FALSE, ignoreNULL = FALSE)

  includeComplexity <- function(value="value") {
    include(value=value,
            plot=readstats_plot,
            type="complexity_plot",
            input=list(react$curr, count_matrix=input$matrix_stats,
                       code_set=input$codeset_stats, plot="complexity",
                       fill=input$fill_complexity,
                       labels=input$labels_stats)
    )}

  # Ratio plot
  output$ratio_plot = renderPlotly({
    #graph <-
    readstats_plot_old(react$curr, count_matrix=input$matrix_stats,
                       code_set=input$codeset_stats, plot="ratio",
                       labels=input$labels_stats)
    #+ theme(axis.text.y = element_text(angle = 90, vjust = 1))
    #+ scale_y_discrete(guide = guide_axis(angle = 90)))
  })

  observeEvent(
    input$includeRatio, {
      includeRatio()
    }, ignoreInit = FALSE, ignoreNULL = FALSE)

  includeRatio <- function(value="value") {
    print(input$matrix_stats)
    include(value=value,
            plot=readstats_plot_old,
            type="ratio_plot",
            input=list(react$curr, count_matrix=input$matrix_stats,
                       code_set=input$codeset_stats, plot="ratio",
                       labels=input$labels_stats)
    )}

  # Both plot
  output$both_plot = renderPlotly({
    #graph <-
    readstats_plot_old(react$curr, count_matrix=input$matrix_stats,
                       code_set=input$codeset_stats, plot="both",
                       labels=input$labels_stats)

  })

  observeEvent(
    input$includeBoth, {
      includeBoth()
    }, ignoreInit = FALSE, ignoreNULL = FALSE)

  includeBoth <- function(value="value") {
    include(value=value,
            plot=readstats_plot_old,
            type="both_plot",
            input=list(react$curr, count_matrix=input$matrix_stats,
                       code_set=input$codeset_stats, plot="both",
                       labels=input$labels_stats)
    )}

  # sample size plot

  output$sample_size_plot = renderPlotly({
    #size_plot_sample(input$size_sample)
    # print(input$size_sample)
    # print(input$matrix_stats)
    #print(input$size_alias)


    size_plot(lox = react$curr,
              sample = switch(input$labels_stats, "sample" = input$size_sample,
                              "alias" = input$size_alias),
              count_matrix = input$matrix_stats,
              code_set = input$codeset_stats,
              labels = input$labels_stats)
    #print(input$labels_stats)

  })

  observeEvent(
    input$includeSampleSize, {
      includeSampleSize()
    }, ignoreInit = FALSE, ignoreNULL = FALSE)

  includeSampleSize <- function(value="value") {
    if (switch(input$labels_stats, "sample" = input$size_sample,
               "alias" = input$size_alias)==""){
      include(value=value,
              plot=size_plot,
              type="sample_size_plot",
              input=list(lox = react$curr,
                         sample = lox@samp_table$sample[1],
                         count_matrix = input$matrix_stats,
                         code_set = input$codeset_stats,
                         labels = input$labels_stats))

    } else{
      include(value=value,
              plot=size_plot,
              type="sample_size_plot",
              input=list(lox = react$curr,
                         sample = switch(input$labels_stats, "sample" = input$size_sample,
                                         "alias" = input$size_alias),
                         count_matrix = input$matrix_stats,
                         code_set = input$codeset_stats,
                         labels = input$labels_stats))
    }
  }

  # sample complexity plot
  output$sample_complexity_plot = renderPlotly({
    dist_orig_plot(lox = react$curr,
                   sample = switch(input$labels_stats, "sample" = input$complexity_sample,
                                   "alias" = input$complexity_alias),
                   count_matrix = input$matrix_stats,
                   code_set = input$codeset_stats,
                   labels = input$labels_stats)

    # dist_orig_plot(lox = react$curr,
    #                sample = input$complexity_sample,
    #                count_matrix = input$matrix_stats,
    #                code_set = input$codeset_stats,
    #                labels = input$labels_stats)

  })

  observeEvent(
    input$includeSampleComplexity, {
      includeSampleComplexity()
    }, ignoreInit = FALSE, ignoreNULL = FALSE
  )

  includeSampleComplexity <- function(value="value") {
    if (switch(input$labels_stats, "sample" = input$size_sample,
               "alias" = input$size_alias)==""){
      include(value=value,
              plot=dist_orig_plot,
              type="sample_complexity",
              input=list(lox = react$curr,
                         sample = lox@samp_table$sample[1],
                         count_matrix = input$matrix_stats,
                         code_set = input$codeset_stats,
                         labels = input$labels_stats))

    } else{
      include(value=value,
              plot=dist_orig_plot,
              type="sample_complexity",
              input=list(lox = react$curr,
                         sample = input$complexity_sample,
                         count_matrix = input$matrix_stats,
                         code_set = input$codeset_stats,
                         labels = input$labels_stats))
    }
  }

  # selection by sample name or alias
  observe(
    if (is.null(react$curr)) { return() }
    else{
      aliases = react$curr@alias[[input$matrix_stats]]
      if (is.null(aliases)) { return() }
      else {
        updateSelectInput(session, "size_sample", "Samples", choices=aliases$sample_name)
        updateSelectInput(session, "complexity_sample", "Samples", choices=aliases$sample_name)
        updateSelectInput(session, "size_alias", "Samples", choices=aliases$alias)
        updateSelectInput(session, "complexity_alias", "Samples", choices=aliases$alias)
      }
    }
  )

  # coordinate sample names and aliases
  # observe(
  #   if (is.null(react$curr)) { return() }
  #   else{
  #     aliases = react$curr@alias[[input$matrix_stats]]
  #
  #     if (is.null(aliases) | input$size_sample=="" | input$size_alias==""| input$complexity_sample=="" | input$complexity_alias=="") { return () }
  #
  #     else if (input$labels_stats=='alias') {
  #       selected_size_sample = get_samplename(react$curr, input$matrix_stats, input$size_alias)
  #       updateSelectInput(session, "size_sample", "Samples",
  #                         choices = aliases$sample_name,
  #                         selected = selected_size_sample)
  #       selected_complexity_sample = get_samplename(react$curr, input$matrix_stats, input$complexity_alias)
  #       updateSelectInput(session, "complexity_sample", "Samples",
  #                         choices = aliases$sample_name,
  #                         selected = selected_complexity_sample)
  #     }
  #
  #     else if (input$labels_stats=='sample') {
  #       selected_size_sample = get_alias(react$curr, input$matrix_stats, input$size_sample)
  #       updateSelectInput(session, "size_alias", "Samples",
  #                         choices = aliases$alias,
  #                         selected = selected_size_sample)
  #       selected_complexity_sample = get_alias(react$curr, input$matrix_stats, input$complexity_sample)
  #       updateSelectInput(session, "complexity_alias", "Samples",
  #                         choices = aliases$alias,
  #                         selected = selected_complexity_sample)
  #     }
  #   }
  # )



  ### HEATMAP PLOT
  # selection parameters for the correlation plot
  observe({
    if (is.null(react$curr)) {
      return ()
    }
    else {
      updateSelectInput(session,
                        "correlation_method",
                        label="Method",
                        choices = c("pearson","kendall", "spearman"),
                        selected = input$correlation_method)

      # params = getMetaParameters(react$curr)
      # updateSelectInput(session,
      #                   "correlation_parameter",
      #                   label="Parameter",
      #                   choices = params,
      #                   selected = input$correlation_parameter)

      if (!is.null(input$correlation_parameter)) {
        data = unique(react$curr@meta[[input$correlation_parameter]])
        updateSelectInput(session,
                          "correlation_data",
                          label="Data",
                          choices=data,
                          selected = input$correlation_data)
      }
    }
  })

  # correlation plot output
  output$correlation_plot = renderPlot({
    correlation_plot(lox = react$curr,
                     count_matrix = input$matrix_heat,
                     code_set = input$codeset_heat,
                     # parameter = input$correlation_parameter,
                     # data = input$correlation_data,
                     method_ = input$correlation_method)
  })



  output$heatmap_ggplot = renderPlot({
    if(input$type_heat=="ggplot"){
      #print("hi")
      heatmap_plot(react$curr,
                   count_matrix=input$matrix_heat,
                   code_set=input$codeset_heat,
                   style=input$type_heat, labels=input$labels_heat,
                   clustering=input$clustering,
                   agglomeration=input$agglomeration,
                   min_reads=input$min_reads,
                   max_repeats=input$max_repeats,
                   min_repeats=input$min_repeats,
                   split_by1=input$split_by1,
                   split_by2=input$split_by2)
      #+ scale_y_discrete())
    }
  })

  # output$heatmap_plotly = renderPlotly({
  #   if(input$type_heat=="plotly"){
  #   heatmap_plot(react$curr, count_matrix=input$matrix_heat,
  #                code_set=input$codeset_heat, style=input$type_heat,
  #                labels=input$labels_heat, clustering=input$clustering,
  #                agglomeration=input$agglomeration,
  #                min_reads=input$min_reads,max_repeats=input$max_repeats,
  #                min_repeats=input$min_repeats)
  #     #+ scale_y_discrete()
  #   }
  # })

  output$bubble_ggplot = renderPlot({
    if(input$type_heat=="ggplot"){
      bubble_plot(react$curr, count_matrix=input$matrix_heat,
                  code_set=input$codeset_heat, style=input$type_heat,
                  labels=input$labels_heat, clustering=input$clustering,
                  agglomeration=input$agglomeration,min_reads=input$min_reads,
                  max_repeats=input$max_repeats,min_repeats=input$min_repeats,split_by1=input$split_by1,split_by2=input$split_by2)
    }
  })
  ################## bug############fix#############
  output$bubble_plotly = renderPlotly({
    if(input$type_heat=="plotly"){bubble_plot(react$curr, count_matrix=input$matrix_heat,
                                              code_set=input$codeset_heat, style=input$type_heat,
                                              labels=input$labels_heat,
                                              clustering=input$clustering,
                                              agglomeration=input$agglomeration,
                                              min_reads=input$min_reads,
                                              max_repeats=input$max_repeats,
                                              min_repeats=input$min_repeats,split_by1=input$split_by1,split_by2=input$split_by2)
    }
  })

  observeEvent(
    input$includeHeatmap, {
      shinyalertAnnotate(session, includeHeatmap)
    })

  includeHeatmap <- function(value="value") {
    print('added')
    include(value=value,
            plot=heatmap_plot,
            type="heatmap_plot",
            input=list(loxcode_experiment=react$curr,
                       count_matrix=input$matrix_heat,
                       code_set=input$codeset_heat,
                       style=input$type_heat, labels=input$labels_heat,
                       clustering=input$clustering,
                       agglomeration=input$agglomeration,
                       min_reads=input$min_reads,
                       max_repeats=input$max_repeats,
                       min_repeats=input$min_repeats,
                       split_by1=input$split_by1,
                       split_by2=input$split_by2)
    )
    print('added')
  }

  observeEvent(
    input$includeBubble, {
      shinyalertAnnotate(session, includeBubble)
    })

  includeBubble <- function(value="value") {
    print(value)
    include(value=value,
            plot=bubble_plot,
            type="bubble_plot",
            input=list(loxcode_experiment=react$curr, count_matrix=input$matrix_heat, code_set=input$codeset_heat, style=input$type_heat, labels=input$labels_heat, clustering=input$clustering,agglomeration=input$agglomeration,min_reads=input$min_reads,max_repeats=input$max_repeats,min_repeats=input$min_repeats,split_by1=input$split_by1,split_by2=input$split_by2)
    )}

  output$sample_comparison_pie = renderPlotly({
    sample_comparison_pie(react$curr, input$matrix_heat,
                          input$codeset_heat,
                          scale=as.numeric(input$scale_pie),
                          labels=input$labels_heat)
  })

  ### SATURATION PLOT
  output$saturation = renderPlotly({

    # worked:
    # sapply(sat$plots, grid.arrange)
    # ggplotly(sat$curr_plot)

    #sapply(sat$plots, grid.arrange)
    #ggplotly(sat$curr_plot)
    s=subplot(sat$plots,
              nrows = ceiling(length(sat$plots)/2),
              shareX = TRUE, shareY = TRUE,
              titleX = FALSE, titleY = FALSE)
    print('hello')
    ggplotly(s,height=100)

    # plots <- lapply(vars, function(var) {
    #   plot_ly(economics, x = ~date, y = as.formula(paste0("~", var))) %>%
    #     add_lines(name = var)
    # })


  })

  observeEvent(
    input$includeSaturation, {
      includeSaturation()
    }, ignoreInit = FALSE, ignoreNULL = FALSE)

  includeSaturation <- function(value="value") {
    if (switch(input$labels_stats, "sample" = input$size_sample,
               "alias" = input$size_alias)==""){
      include(value=value,
              plot=saturation_plot,
              type="saturation_plot",
              input=list(loxcode_experiment=react$curr,
                         loxcode_sample=lox@samp_table$sample[1],
                         code_set=input$view_codeset1))

    } else{
      include(value=value,
              plot=saturation_plot,
              type="saturation_plot",
              input=list(loxcode_experiment=react$curr,
                         loxcode_sample=input$sample_sat,
                         code_set=input$view_codeset1)
      )}
  }

  # add plot
  observeEvent(
    input$add_sat, {
      sat$curr_plot <- saturation_multi(react$curr,
                                        input$sample_sat,
                                        input$codeset_sat)
      print(countsat)


      sat$samples = list.append(sat$samples, input$sample_sat)
      sat$codesets = list.append(sat$codesets, input$codeset_sat)
      sat$plots = list.append(sat$plots, sat$curr_plot)
    })

  # remove last plot
  observeEvent(
    input$remove_sat, {
      n = length(sat$plots)
      if (n==0) { return() }
      else if (n==1) {
        sat$samples <- list()
        sat$codesets <- list()
        sat$plots <- list()
      }
      else {
        sat$samples[[n]] <- NULL
        sat$codesets[[n]] <- NULL
        sat$plots[[n]] <- NULL
      }
    })

  # clear plot
  observeEvent(
    input$clear_sat, {
      sat$samples <- list()
      sat$codesets <- list()
      sat$plots <- list()
    })

  # switch between sample names and aliases
  observe(
    if (is.null(react$curr)) { return() }
    else{
      aliases = react$curr@alias[[input$sampleset_sat]]
      if (is.null(aliases)) { return() }
      else {
        updateSelectInput(session, "sample_sat", "Samples", choices=aliases$sample_name)
        updateSelectInput(session, "alias_sat", "Samples", choices=aliases$alias)
      }
    }
  )

  #Num of reads
  output$read_plot = renderPlotly({

    do.call(read_plot, list(react$curr, count_matrix=input$matrix_stats,
                            code_set=input$codeset_stats,
                            labels=input$labels_stats))

  })

  observeEvent(
    input$includeRead, {
      includeRead()
    }, ignoreInit = FALSE, ignoreNULL = FALSE)

  includeRead <- function(value="value") {
    include(value=value,
            plot=read_plot,
            type="read_plot",
            input=list(react$curr, count_matrix=input$matrix_stats,
                       code_set=input$codeset_stats,
                       labels=input$labels_stats)

    )}


  ########### Here
  # coordinate sample names and aliases
  observe(
    if (is.null(react$curr)) { return() }
    else{
      aliases = react$curr@alias[[input$sampleset_sat]]

      if (is.null(aliases) | input$sample_sat=="" | input$alias_sat=="") { return () }

      else if (input$name_sat=='alias') {
        selected_sample = get_samplename(react$curr, input$sampleset_sat, input$alias_sat)
        updateSelectInput(session, "sample_sat", "Samples",
                          choices = aliases$sample_name,
                          selected = selected_sample)
      }

      else if (input$name_sat=='sample') {
        selected_sample = get_alias(react$curr, input$sampleset_sat, input$sample_sat)
        updateSelectInput(session, "alias_sat", "Samples",
                          choices = aliases$alias,
                          selected = selected_sample)
      }
    }
  )

  ############# here
  ## PAIR COMPARISON PLOT V1

  # plot page can extend
  # ggplot
  output$pair_ggplot = renderPlot({
    range = switch(input$colour_pair, "size"=input$size_slider_pair,
                   "complexity"=input$complexity_slider_pair)

    do.call(grid.arrange, react$pairs)

  })


  # plotly
  output$pair_plotly = renderPlotly({
    range = switch(input$colour_pair, "size"=input$size_slider_pair,
                   "complexity"=input$complexity_slider_pair)


    out12=sapply(react$pairs, "[[",1)
    ggplotly(react$curr_pair)
  })

  # output$pair_plotly = renderPlotly({
  #   range = switch(input$colour_pair, "size"=input$size_slider_pair,
  #                  "complexity"=input$complexity_slider_pair)
  #   react$curr_pair <- pair_comparison_plot2(
  #     lox = react$curr,
  #     s1 = react$curr@samples[[input$sample1_pair]],
  #     s2 = react$curr@samples[[input$sample2_pair]],
  #     sampleset = input$sampleset_pair,
  #     codeset = input$codeset_pair,
  #     colorBy = input$colour_pair,
  #     sizeRange = input$size_slider_pair,
  #     dist_origRange = input$complexity_slider_pair,
  #     firstreadRange = input$firstread_slider_pair
  #   )
  #   do.call(grid.arrange, list.append(react$pairs, react$curr_pair))
  # })

  # add new plot
  observeEvent(
    input$add_pair, {
      if (input$name_pair == "alias"){
        f1 = which((react$curr@alias[["all_samples"]]$alias) == input$alias1_pair)
        f2 = which((react$curr@alias[["all_samples"]]$alias) == input$alias2_pair)
        sample1_pair <- react$curr@alias[["all_samples"]]$sample_name[f1]
        sample2_pair <- react$curr@alias[["all_samples"]]$sample_name[f2]
      }
      react$curr_pair <- pair_comparison_plot2(
        lox = react$curr,
        s1 = switch(input$name_pair,
                    "sample" = react$curr@samples[[input$sample1_pair]],
                    "alias" = react$curr@samples[[sample1_pair]]),
        s2 = switch(input$name_pair,
                    "sample" = react$curr@samples[[input$sample2_pair]],
                    "alias" = react$curr@samples[[sample2_pair]]),
        sampleset = input$sampleset_pair,
        codeset = input$codeset_pair,
        colorBy = input$colour_pair,
        labels = input$name_pair,
        sizeRange = input$size_slider_pair,
        dist_origRange = input$complexity_slider_pair,
        firstreadRange = input$firstread_slider_pair
      )
      # View(react$curr@samples[[input$sample1_pair]])
      # View(react$curr@samples[[input$sample2_pair]])
      # View(react$curr@samples[[input$alias2_pair]])
      # print(input$sample1_pair)
      # print(input$sample2_pair)
      # "alias1_pair"
      # View(react$curr@samples[[input$alias2_pair]])
      # View(react$curr)
      # print(input$alias2_pair)
      # print(react$curr@alias[["all_samples"]]$alias == sample)

      #print(get_samplename(input$alias1_pair))
      # f = which((react$curr@alias[["all_samples"]]$alias) == input$alias1_pair)
      # print(input$alias1_pair)
      # print(f)
      # print(react$curr@alias[["all_samples"]]$sample_name[f])




      react$pairs = list.append(react$pairs, react$curr_pair) }
  )

  #add all comb V1 (replace plots)
  observeEvent(
    input$add_all, {
      react$curr_pair <- pair_comparison_plot_all(
        lox = react$curr,
        matrix = react$curr@count_matrixes[[input$sampleset_pair]],
        codeset=input$codeset_pair
      )
      react$pairs <- list()
      react$pairs = list.append(react$pairs, react$curr_pair)

    })



  ### PAIR COMPARISON PLOT V2
  # output$pair_ggplot = renderPlot({
  #   range = switch(input$colour_pair, "size"=input$size_slider_pair,
  #                  "complexity"=input$complexity_slider_pair)
  #
  #   react$curr_pair <- pair_comparison_plot2(
  #     lox = react$curr,
  #     s1 = react$curr@samples[[input$sample1_pair]],
  #     s2 = react$curr@samples[[input$sample2_pair]],
  #     sampleset = input$sampleset_pair,
  #     codeset = input$codeset_pair,
  #     colorBy = input$colour_pair,
  #     labels = input$name_pair,
  #     sizeRange = input$size_slider_pair,
  #     dist_origRange = input$complexity_slider_pair,
  #     firstreadRange = input$firstread_slider_pair
  #   )
  #
  #   do.call(grid.arrange, list.append(react$pairs, react$curr_pair))
  # })
  #
  #
  # # # plotly
  # # output$pair_plotly = renderPlotly({
  # #   range = switch(input$colour_pair, "size"=input$size_slider_pair,
  # #                  "complexity"=input$complexity_slider_pair)
  # #   react$curr_pair <- pair_comparison_plot2(
  # #     lox = react$curr,
  # #     s1 = react$curr@samples[[input$sample1_pair]],
  # #     s2 = react$curr@samples[[input$sample2_pair]],
  # #     sampleset = input$sampleset_pair,
  # #     codeset = input$codeset_pair,
  # #     colorBy = input$colour_pair,
  # #     sizeRange = input$size_slider_pair,
  # #     dist_origRange = input$complexity_slider_pair,
  # #     firstreadRange = input$firstread_slider_pair
  # #   )
  # #   do.call(grid.arrange, list.append(react$pairs, react$curr_pair))
  # # })
  #
  # # add new plot
  # observeEvent(
  #   input$add_pair, { react$pairs = list.append(react$pairs, react$curr_pair) }
  # )
  #
  # # add all comb V2 (preserve all plots)
  # observeEvent(
  #   input$add_all, {
  #     react$curr_pair <- pair_comparison_plot_all( #renderPlot({
  #         lox = react$curr,
  #         matrix = react$curr@count_matrixes[[input$sampleset_pair]],
  #         codeset=input$codeset_pair
  #         )
  #     #})
  #
  #     react$pairs = list.append(react$pairs, react$curr_pair)
  # print(react$pairs)
  # print(class(react$pairs))
  #})

  ######################################################################################

  # remove previous pair plot
  observeEvent(
    input$remove_pair, {
      n = length(react$pairs)
      if (n==0) { return() }
      else if (n==1) { react$pairs <- list() }
      else { react$pairs[[n]] <- NULL }
    }
  )

  # clear all pair plots
  observeEvent(
    input$clear_pair, {
      react$pairs <- list()}
  )

  # add to report
  observeEvent(
    input$includePair, {
      shinyalertAnnotate(session, includePair)
    })

  includePair <- function(value="value") {
    include(value=value,
            type = "pair_comparison_plot",
            input = list.append(react$pairs, react$curr_pair),
            plot= grid.arrange
            #plot= marrangeGrob(react$pairs, nrow=2, ncol=2)
    )}


  #switch between sample names and aliases
  observe(
    if (is.null(react$curr)) { return() }
    else{
      aliases = react$curr@alias[[input$sampleset_pair]]
      if (is.null(aliases)) { return() }
      else {
        updateSelectInput(session, "sample1_pair", "Samples", choices=aliases$sample_name)
        updateSelectInput(session, "sample2_pair", "Samples", choices=aliases$sample_name)

        updateSelectInput(session, "alias1_pair", "Samples", choices=aliases$alias)
        updateSelectInput(session, "alias2_pair", "Samples", choices=aliases$alias)
      }
    }
  )

  # #coordinate sample names and aliases
  # observe(
  #   if (is.null(react$curr)) { return() }
  #   else{
  #     aliases = react$curr@alias[[input$sampleset_pair]]
  #
  #     #if (is.null(aliases) | input$sample1_pair=="" | input$sample2_pair=="" | input$alias1_pair=="" | input$alias2_pair=="") { return () }
  #     if (is.null(aliases)) { return() }
  #
  #     else if (input$name_pair=='sample') {
  #       selected_sample1 = get_samplename(react$curr, input$sampleset_pair, input$alias1_pair)
  #       updateSelectInput(session, "sample1_pair", "Samples",
  #                         choices = aliases$sample_name,
  #                         selected = selected_sample1)
  #       selected_sample2 = get_samplename(react$curr, input$sampleset_pair, input$alias2_pair)
  #       updateSelectInput(session, "sample2_pair", "Samples",
  #                         choices = aliases$sample_name,
  #                         selected = selected_sample2)
  #     }
  #
  #     else if (input$name_pair=='alias') {
  #       selected_sample1 = get_alias(react$curr, input$sampleset_pair, input$sample1_pair)
  #       updateSelectInput(session, "alias1_pair", "Samples",
  #                         choices = aliases$alias,
  #                         selected = selected_sample1)
  #       selected_sample2 = get_alias(react$curr, input$sampleset_pair, input$sample2_pair)
  #       updateSelectInput(session, "alias2_pair", "Samples",
  #                         choices = aliases$alias,
  #                         selected = selected_sample2)
  #     }
  #   }
  # )

  observe({
    # updates the slider based on the distance range of the samples selected
    samples = react$curr@samples
    if (is.null(samples) | is.null(samples[[input$sample1_pair]]) |
        is.null(samples[[input$sample2_pair]])) { return() }
    else {
      updateRange(samples, "complexity_slider_pair", "dist_orig")
      updateRange(samples, "size_slider_pair", "size")
      updateRange(samples, "firstread_slider_pair", "firstread")
    }
  })

  updateRange <- function(samples, slider, type) {
    min_one <- min(na.omit(samples[[input$sample1_pair]]@decode@data[[type]]))
    min_two <- min(na.omit(samples[[input$sample2_pair]]@decode@data[[type]]))
    max_one <- max(na.omit(samples[[input$sample1_pair]]@decode@data[[type]]))
    max_two <- max(na.omit(samples[[input$sample2_pair]]@decode@data[[type]]))
    newmin <- min(min_one, min_two)
    newmax <- max(max_one, max_two)
    updateSliderInput(session, slider, value = c(newmin,newmax), min=newmin, max=newmax)
  }

  ### DOWNLOAD REPORT

  output$components_table = renderDataTable(server=FALSE,{
    datatable(components_table(params),
              class = "cell-border stripe",
              editable = list(target="cell", disable=list(columns=c(0, 1, 2))),
              colnames = c(ID = 1),
              extensions = 'RowReorder',
              options = list(rowReorder=TRUE, order = list(c(0 , 'asc'))),
              callback=JS("// pass on data to R
                          table.on('row-reorder', function(e, details, changes) {
                            Shiny.onInputChange('components_table_row_reorder', JSON.stringify(details));
                          });")
    )})

  components_table <- function(params) {
    d = data.frame()
    if (length(params$functions) == 0) {
      return(data.frame())
    }
    for (i in 1:length(params$functions)) {
      row = data.frame("Plot_type" = params$types[[i]],
                       "Experiment" = params$loxcodes[[i]]@name,
                       "Annotation" = params$annotations[[i]],
                       "Sample name" = names(params$inputs[[i]][[1]]@count_matrixes)[length(names(params$inputs[[i]][[1]]@count_matrixes))],
                       #"Notes" = params$notes[[i]],
                       stringsAsFactors = FALSE)
      d = plyr::rbind.fill(d, row)
    }
    return(d)
  }

  # edit annotations
  comp_proxy = dataTableProxy("components_table")
  observeEvent(
    input$components_table_cell_edit, {
      d = components_table(params)
      info = input$components_table_cell_edit
      i = info$row
      j = info$col
      v = info$value
      d[i, j] <<- coerceValue(v, d[i, j])
      replaceData(proxy, d, resetPaging=FALSE, rownames=FALSE)
      params$annotations[[i]] <- v
    })

  # reorder the rows
  observeEvent(input$components_table_row_reorder, {
    info <- input$components_table_row_reorder
    # error checking
    if(is.null(info) | class(info) != 'character') { return() }

    info <- read_yaml(text=info)
    saveRDS(info,"row.rds")

    if(length(info) == 0) { return() }
    reorder(info)
  })
  print("*********************")
  print(params)
  print("*********************")
  reorder <- function(info) {
    temp = list(functions=list(), types=list(), inputs=list(), annotations=list())
    #, notes = list())
    temp$functions <- params$functions
    temp$types <- params$types
    temp$inputs <- params$inputs
    temp$annotations <- params$annotations
    #temp$notes <- params$notes

    temp$loxcodes <- params$loxcodes
    #####uSure
    #temp$plots <- params$plots
    #####uSure
    for (i in 1:length(info)) {
      curr=info[[i]]
      temp$functions[[curr$newPosition + 1]] = params$functions[[curr$oldPosition + 1]]
      temp$types[[curr$newPosition + 1]] = params$types[[curr$oldPosition + 1]]
      temp$inputs[[curr$newPosition + 1]] = params$inputs[[curr$oldPosition + 1]]
      temp$annotations[[curr$newPosition + 1]] = params$annotations[[curr$oldPosition + 1]]
      #temp$notes[[curr$newPosition + 1]] = params$notes[[curr$oldPosition + 1]]

      temp$loxcodes[[curr$newPosition + 1]] = params$loxcodes[[curr$oldPosition + 1]]
      #####uSure
      #temp$plots[[curr$newPosition + 1]] = params$plots[[curr$oldPosition + 1]]
      ##############
    }
    params$functions = temp$functions
    params$types = temp$types
    params$inputs = temp$inputs
    params$annotations = temp$annotations
    #params$notes = temp$notes

    params$loxcodes = temp$loxcodes
    #####uSure
    #params$plots = temp$plots
    ##############
  }

  # remove a component
  observeEvent(
    input$remove_component, {
      rows = input$components_table_rows_selected
      print(rows)
      if (is.null(rows)) {
        showModal(modalDialog(
          title = "Oops!",
          "Please select at least one plot to remove",
          footer = modalButton("Ok")
        ))
        return() }
      params$functions[rows] <- NULL
      params$types[rows] <- NULL
      params$inputs[rows] <- NULL
      params$annotations[rows] <- NULL
      #params$notes[rows] <- NULL

      params$loxcodes[rows] <- NULL
      #####uSure
      #params$plots[rows] <- NULL
      ##############
    }
  )
  output$experiments_table1 = renderDataTable({datatable(
    exp_table(react$exp),
    rownames = FALSE,
    class = "cell-border stripe",
    filter = 'top',
    selection = 'multiple'
  )})

  observeEvent(
    input$select_exp1, {
      if (is.null(input$experiments_table1_rows_selected) |
          length(input$experiments_table1_rows_selected) > 1) {
        print(length(input$experiments_table1_rows_selected))
        showModal(modalDialog(
          title = "Oops!",
          "Please select one experiment!",
          footer = modalButton("Ok")
        ))

        return ()
      }
      else if (length(input$experiments_table1_rows_selected) == 1) {
        react$curr = react$exp[[input$experiments_table1_rows_selected]]
        showNotification(paste(react$curr@name, " selected."))
      }
    }
  )


  ##Size Plot box
  observeEvent(
    input$includeSize1, {
      if(input$includeSize1==TRUE){
        includeSize()
      }
      else{
        for (i in 1:length(params$functions)){
          if (params$types[i] == "size_plot"){
            params$functions[i] <- NULL
            params$types[i] <- NULL
            params$inputs[i] <- NULL
            params$annotations[i] <- NULL
            #params$notes[i] <- NULL

            params$loxcodes[i] <- NULL
            #####uSure
            #params$plots[i] <- NULL
            ############
          }
        }
      }
    }, ignoreInit = TRUE)
  ##Complexity plot box
  observeEvent(
    input$includeComplexity1, {
      if(input$includeComplexity1==TRUE){
        includeComplexity()
      }
      else{
        for (i in 1:length(params$functions)){
          if (params$types[i] == "complexity_plot"){
            params$functions[i] <- NULL
            params$types[i] <- NULL
            params$inputs[i] <- NULL
            params$annotations[i] <- NULL
            #params$notes[i] <- NULL

            params$loxcodes[i] <- NULL
          }
        }
      }
    }, ignoreInit = TRUE)

  ##Ratio plot box
  observeEvent(
    input$includeRatio1, {
      if(input$includeRatio1==TRUE){
        includeRatio()
      }
      else{
        for (i in 1:length(params$functions)){
          if (params$types[i] == "ratio_plot"){
            params$functions[i] <- NULL
            params$types[i] <- NULL
            params$inputs[i] <- NULL
            params$annotations[i] <- NULL
            #params$notes[i] <- NULL

            params$loxcodes[i] <- NULL
          }
        }
      }
    }, ignoreInit = TRUE)

  ##Both plot box
  observeEvent(
    input$includeBoth1, {
      if(input$includeBoth1==TRUE){
        includeBoth()
      }
      else{
        for (i in 1:length(params$functions)){
          if (params$types[i] == "both_plot"){
            params$functions[i] <- NULL
            params$types[i] <- NULL
            params$inputs[i] <- NULL
            params$annotations[i] <- NULL
            #params$notes[i] <- NULL

            params$loxcodes[i] <- NULL
          }
        }
      }
    }, ignoreInit = TRUE)

  ##Sample size plot box
  observeEvent(
    input$includeSampleSize1, {
      if(input$includeSampleSize1==TRUE){
        includeSampleSize()
      }
      else{
        for (i in 1:length(params$functions)){
          if (params$types[i] == "sample_size_plot"){
            params$functions[i] <- NULL
            params$types[i] <- NULL
            params$inputs[i] <- NULL
            params$annotations[i] <- NULL
            #params$notes[i] <- NULL

            params$loxcodes[i] <- NULL
          }
        }
      }
    }, ignoreInit = TRUE)

  ##Sample Compx plot box
  observeEvent(
    input$includeSampleComplexity1, {
      if(input$includeSampleComplexity1==TRUE){
        includeSampleComplexity()
      }
      else{
        for (i in 1:length(params$functions)){
          if (params$types[i] == "sample_complexity"){
            params$functions[i] <- NULL
            params$types[i] <- NULL
            params$inputs[i] <- NULL
            params$annotations[i] <- NULL
            #params$notes[i] <- NULL

            params$loxcodes[i] <- NULL
          }
        }
      }
    }, ignoreInit = TRUE)

  ##Heat map plot box
  observeEvent(
    input$includeHeatmap1, {
      if(input$includeHeatmap1==TRUE){
        includeHeatmap()
      }
      else{
        for (i in 1:length(params$functions)){
          if (params$types[i] == "heatmap_plot"){
            params$functions[i] <- NULL
            params$types[i] <- NULL
            params$inputs[i] <- NULL
            params$annotations[i] <- NULL
            #params$notes[i] <- NULL

            params$loxcodes[i] <- NULL
          }
        }
      }
    }, ignoreInit = TRUE)

  ##bubble plot box
  observeEvent(
    input$includeBubble1, {
      if(input$includeBubble1==TRUE){
        includeBubble()
      }
      else{
        for (i in 1:length(params$functions)){
          if (params$types[i] == "bubble_plot"){
            params$functions[i] <- NULL
            params$types[i] <- NULL
            params$inputs[i] <- NULL
            params$annotations[i] <- NULL
            #params$notes[i] <- NULL

            params$loxcodes[i] <- NULL
          }
        }
      }
    }, ignoreInit = TRUE)

  ## Saturation plot box
  observeEvent(
    input$includeSaturation1, {
      if(input$includeSaturation1==TRUE){
        includeSaturation()
      }
      else{
        for (i in 1:length(params$functions)){
          if (params$types[i] == "saturation_plot"){
            params$functions[i] <- NULL
            params$types[i] <- NULL
            params$inputs[i] <- NULL
            params$annotations[i] <- NULL
            #params$notes[i] <- NULL

            params$loxcodes[i] <- NULL
          }
        }
      }
    }, ignoreInit = TRUE)





  # download report
  output$downloadReport <- downloadHandler(
    filename = function() {
      file = paste('my-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
      addToLog(session, logDownloadReport(session, file, input$format))


      l1temp=c()
      for (i in 7:length(params$types)){
        print(params$inputs[[i]][[2]])
        if (!(params$inputs[[i]][[2]] %in% l1temp)){
          l1temp=c(l1temp,params$inputs[[i]][[2]])
          #cat(paste("### text"))
        }
      }
      tnum=2
      rmd_content=""
      print("######################")
      print(l1temp)
      print("######################")
      if(length(l1temp)>2){
      for (i in l1temp[3:length(l1temp)]){
        rmd_title=paste("### ",i,sep="")
        rmd_content=paste(rmd_content,rmd_title,sep = "\n")
        rmd_content=paste(rmd_content,"```{r results='asis'}
plt <- htmltools::tagList()
if (length(params$types)>5){
for (i in 8:length(params$types)){
  if(length(list1)>0){
  if (params$inputs[[i]][[2]]==list1[1]){
if(params$types[[i]]=='ratio_plot' || params$types[[i]]=='both_plot'){
g <- do.call(params$functions[[i]], params$inputs[[i]])
plt[[i]]=ggplotly(g,height = 6000) }else if(params$types[[i]]=='sample_complexity'){
  g <- do.call(params$functions[[i]], params$inputs[[i]])
plt[[i]]=ggplotly(g) }else{
 g <- do.call(params$functions[[i]], params$inputs[[i]])
plt[[i]]=ggplotly(g)
#params$inputs[[i]]

cat('\n')
}
}
}
}
  if(length(list1)<2){list1=c()}
   else{list1=list1[2:length(list1)]}
}
```

```{r results='asis',autodep=TRUE}
plt
```",sep="\n")}
}
writeLines(rmd_content, "./chapter1.Rmd")


l1temp=c()
l2temp=c()
for (i in 8:length(params$types)){
  print(params$inputs[[i]][[1]]@code_sets[["test"]])
  if (length(names(params$inputs[[i]][[1]]@code_sets))>21){
    for (j in 22:length(names(params$inputs[[i]][[1]]@code_sets))){
      if(!(names(params$inputs[[i]][[1]]@code_sets)[j] %in% l1temp))
        l1temp=c(names(params$inputs[[i]][[1]]@code_sets)[22],l1temp)
      l2temp=c(i,l2temp)
    }
    #cat(paste("### text"))
  }
}
rmd_content=""
print("***")
print(l1temp)
print(l2temp)
print("***")
if (length(l2temp)>0){
  for (i in 1:length(l2temp)){
    rmd_title=paste("### ",l1temp[i],sep="")
    rmd_content=paste(rmd_content,rmd_title,sep = "\n")
    rmd_content=paste(rmd_content,"",sep = "\n")
    rmd_content=paste(rmd_content,"```{r }
datatable(params$inputs[[",l2temp[i],"]][[1]]@code_sets[['",l1temp[i],"']], options = list(pageLength = 5, scrollY = '2000px'))
```",sep="")}}
writeLines(rmd_content, "./chapter3.Rmd")




l1temp=c()
l2temp=c()
for (i in 8:length(params$types)){
  print(params$inputs[[i]][[2]])
  if (!(params$inputs[[i]][[2]] %in% l1temp)){
    if (!(params$inputs[[i]][[2]] %in% c("all_samples"))){
      l1temp=c(params$inputs[[i]][[2]],l1temp)
      l2temp=c(i,l2temp)
    }
    #cat(paste("### text"))
  }
}
rmd_content=""
print("***")
print(l1temp)
print(l2temp)
print("***")
if(length(l2temp)>0){
  for (i in 1:length(l2temp)){
    rmd_title=paste("### ",l1temp[i],sep="")
    rmd_content=paste(rmd_content,rmd_title,sep = "\n")
    rmd_content=paste(rmd_content,"",sep = "\n")
    rmd_content=paste(rmd_content,"```{r }
datatable(params$inputs[[",l2temp[i],"]][[1]]@samp_table[match(names(params$inputs[[",l2temp[i],"]][[1]]@count_matrixes[[names(params$inputs[[",l2temp[i],"]][[1]]@count_matrixes)[length(names(params$inputs[[",l2temp[i],"]][[1]]@count_matrixes))]]]),params$inputs[[",l2temp[i],"]][[1]]@samp_table[[1]]),], options = list(pageLength = 5, scrollY = '2000px'))
```",sep="")}
}
writeLines(rmd_content, "./chapter4.Rmd")



return(file)
    },

content = function(file) {
  src <- switch(input$format,
                PDF = normalizePath('report2.Rmd'),
                HTML = normalizePath('repotemp.Rmd'),
                Word = normalizePath('report2.Rmd'))
  # temporarily switch to the temp dir, in case you do not have write
  # permission to the current working directory
  owd <- setwd(tempdir())
  on.exit(setwd(owd))
  file.copy(src, switch(input$format,
                        PDF = 'report2.Rmd',
                        HTML = 'repotemp.Rmd',
                        Word = 'report2.Rmd'),
            overwrite = TRUE)

  # set up parameters
  params <- list(
    format = input$format,
    functions = params$functions,
    types = params$types,
    inputs = params$inputs,
    annotations = params$annotations,
    #notes = params$notes,

    loxcodes = params$loxcodes
  )

  out <- render(
    switch(input$format,
           PDF = 'report2.Rmd',
           HTML = 'repotemp.Rmd',
           Word = 'report2.Rmd'),
    switch(input$format, PDF = pdf_document(),
           HTML = flex_dashboard(),
           Word = word_document()),
    params=params,
    envir = new.env(parent = globalenv()))
  file.rename(out, file)
}
  )

  #   content = function(file) {
  #     params <- list(
  #       format = input$format,
  #       functions = params$functions,
  #       types = params$types,
  #       inputs = params$inputs,
  #       annotations = params$annotations,
  #       loxcodes = params$loxcodes
  #     )
  #
  #     if (params$format == HTML){
  #       #src <- normalizePath('report2.Rmd')
  #       src <- normalizePath('repo.Rmd')
  #       # temporarily switch to the temp dir, in case you do not have write
  #       # permission to the current working directory
  #       owd <- setwd(tempdir())
  #       on.exit(setwd(owd))
  #       #file.copy(src, 'report2.Rmd', overwrite = TRUE)
  #       file.copy(src, 'repo.Rmd', overwrite = TRUE)
  #
  #       out <- render(#'report2.Rmd',
  #         'repo.Rmd',
  #         flex_dashboard(),
  #         params=params,
  #         envir = new.env(parent = globalenv()))
  #       file.rename(out, file)}
  #
  #     else{
  #       src <- normalizePath('report2.Rmd')
  #       # temporarily switch to the temp dir, in case you do not have write
  #       # permission to the current working directory
  #       owd <- setwd(tempdir())
  #       on.exit(setwd(owd))
  #       file.copy(src, 'report2.Rmd', overwrite = TRUE)
  #
  #     out <- render('report2.Rmd',
  #                   switch(input$format,
  #                          PDF = pdf_document(),
  #                          Word = word_document()),
  #                   params=params,
  #                   envir = new.env(parent = globalenv()))
  #
  #     }
  # )


  ##############worked
  # content = function(file) {
  #   src <- normalizePath('report.Rmd')
  #   # temporarily switch to the temp dir, in case you do not have write
  #   # permission to the current working directory
  #   owd <- setwd(tempdir())
  #   on.exit(setwd(owd))
  #   file.copy(src, 'report.Rmd', overwrite = TRUE)
  #
  #   # set up parameters
  #   params <- list(
  #     format = input$format,
  #     functions = params$functions,
  #     types = params$types,
  #     inputs = params$inputs,
  #     annotations = params$annotations,
  #     #notes = params$notes,
  #
  #     loxcodes = params$loxcodes
  #   )
  #
  #   out <- render(
  #     'report.Rmd',
  #     switch(input$format, PDF = pdf_document(),
  #            #HTML = html_document(),
  #            HTML = includeMarkdown(knitr::knit('report.Rmd')),
  #            Word = word_document()),
  #     params=params,
  #     envir = new.env(parent = globalenv()))
  #   file.rename(out, file)
  # }
  # )


  ### ACTIVITY LOG
  output$log_table = renderDataTable({datatable(
    log_table(logs),
    rownames = FALSE,
    class = "cell-border stripe",
    selection = 'none',
    options = list(dom='t')
  )})

  log_table <- function(timestamps, activity) {
    d = data.frame()
    for (i in 1:length(logs$timestamps)) {
      row = data.frame("Time" = logs$timestamps[[i]],
                       "Activity" = logs$activity[[i]],
                       stringsAsFactors=FALSE)
      d = plyr::rbind.fill(d, row)
    }
    return(d)
  }


  # restart the session
  observeEvent(
    input$restart, {
      js$reset()
      addToLog(session, refreshSession)
    })

  # download log
  output$downloadLog = downloadHandler(
    filename="temp.html",
    content=NULL
  )
}


