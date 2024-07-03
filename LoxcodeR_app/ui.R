

#########0#########1#########2#########3#########4#########5#########6#########7

# Import Libraries
library(shiny)
library(shinyDirectoryInput)
library(plotly)
library(DT)
library(loxcodeR)
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
library(tinytex)
library(pals)
library(ggplot2)
library(plyr)
library(dplyr)
#########0#########1#########2#########3#########4#########5#########6#########7

# Source modules


#########0#########1#########2#########3#########4#########5#########6#########7

### INITIALIZE
# load_origin_distmaps('/wehisan/general/user_managed/grpu_naik.s_2/TW/maps/origin/')
# load_pair_distmaps('/wehisan/general/user_managed/grpu_naik.s_2/TW/maps/origin/')

datards <- system.file("extdata","data-2024-05-20.rds",package="loxcodeR")
lox <- readRDS(datards)
exp = list(lox)





### CONSTANTS
chart_choices = c("Statistics Plots", "Heatmap", "Saturation Plot", "Pair Comparison Plot")

header <- dashboardHeader(title = "LoxCodeR",
                          dropdownMenuOutput("curr_lox"))

#########0#########1#########2#########3#########4#########5#########6#########7

# Side bar menu
sidebar <- dashboardSidebar(

    tabPanel("View code set", tabName = "codeset-create", "codes-filter",
             selectInput("view_codeset","filter_codeset", label="Choose code set:",
                         choices=names(lox@code_sets)),
             actionButton("delete_codeset", "Delete Code Set")
    ),

    tabPanel(
        "View sample set", tabName = "sample-create",
        selectInput("view_sample", "independent_samples", label="Choose sample set:",
                    choices=names(lox@count_matrixes)),
        actionButton("delete_sample", "Delete Sample Set")),


    # box(
    #   width = 6,
    #   title = "View codes",
    #   status = NULL,
    #   color = NULL,
    #   solidHeader = TRUE,
    #   collapsible = TRUE,
    #   selectInput("independent_samples",
    #               label="Sample Set", choices=names(lox@count_matrixes)),
    #   selectInput("filter_codeset", label="Code Set", choices=c()),
    # ),

    sidebarMenu(
        menuItem("Import", tabName="import", icon=icon("upload")),
        menuItem(
            "Create", tabName="create", icon=icon("folder-plus"),
            menuSubItem("Create Code Sets", tabName="codeset-create"),
            menuSubItem("Create Sample Sets", tabName="sample-create"),
            menuSubItem("Filter Codes", tabName="codes-filter")
        ),
        menuItem(
            "Plots", tabName="plots", icon=icon("chart-bar"),
            menuSubItem("Statistics Plots", tabName="stats-plot"),
            menuSubItem("Heat Map", tabName="heatmap-plot"),
            menuSubItem("Saturation Plot", "saturation-plot"),
            menuSubItem("Pair Comparison Plot", "pair-plot")
        ),
        #menuItem("Simulate", tabName="simulate", icon=icon("vials")),
        menuItem("Report", tabName="report", icon=icon("file-alt")),
        menuItem("Log", tabName="log", icon=icon("history")),
        tabPanel(
            "save_data", tabName = "save-data",
            actionButton("save_data", "Export")),
        tabPanel(
            "Reload Data", tabName = "reload-data",
            actionButton("reload_data", "Import")),
        tabPanel("Refresh Graph",tabName = "refresh-graph",actionButton("refresh", "Refresh Graph"),)

        # menuItem(
        #   "View code set", tabName = "codeset-create",
        #   selectInput("view_codeset", label="Choose code set:",
        #               choices=names(lox@code_sets)),
        #   actionButton("delete_codeset", "Delete Code Set"),
        #   startExpanded = TRUE, collapsed = FALSE, collapsible = FALSE
        # ),
        #
        # menuItem(
        #   "View sample set", tabName = "sample-create",
        #   selectInput("view_sample", label="Choose sample set:",
        #               choices=names(lox@count_matrixes)),
        #   actionButton("delete_sample", "Delete Sample Set"),
        #   startExpanded = TRUE
        # )

    )


    # sidebarUserPanel("User Name",
    #                  subtitle = a(href = "#", icon("circle", class = "text-success"), "Online"),
    #                  # Image file should be in www/ subdir
    #                  image = "userimage.png"
    # ),


)

# tabName = "sample-create",
# fluidRow(
#   box(
#     width = 6,
#     title = "View a sample set",
#     status = NULL,
#     color = NULL,
#     solidHeader = TRUE,
#     collapsible = TRUE,
#     selectInput("view_sample", label="Choose sample set:",
#                 choices=names(lox@count_matrixes)),
#     actionButton("delete_sample", "Delete Sample Set")

# tabName = "codeset-create",
# fluidRow(
#   box(
#     title = "View a code set",
#     status = NULL,
#     color = NULL,
#     solidHeader = TRUE,
#     collapsible = TRUE,
#     selectInput("view_codeset", label="Choose code set:", choices=names(lox@code_sets)),
#     actionButton("delete_codeset", "Delete Code Set")
#   )



#########0#########1#########2#########3#########4#########5#########6#########7

# Simulation tab
simulationTab <- tabItem(
    tabName = "simulate",
    fluidRow(
        box(
            title = "Simulation Parameters",
            status = NULL,
            color = NULL,
            solidHeader = TRUE,
            collapsible = TRUE,
            numericInput(
                "nsamples",
                label = "Number of samples to generate:",
                10,
                min = 1,
                max = 50
            ),
            numericInput(
                "ncodes",
                label = "Number of cells with barcodes",
                3000,
                min = 100,
                max = 10000
            ),
            textInput("sim_names", label = "Name of simulated samples:", value =
                          "simulation"),
            actionButton("simulate_sample", label = "Simulate")

        ),

        box(
            title = "Sampling Distribution",
            status = NULL,
            color = NULL,
            solidHeader = TRUE,
            collapsible = TRUE,
            selectInput("dist_orig_type", label="Sample or Poisson distribution", choices=c("Sample Name", "Sample Alias", "Poisson")),
            conditionalPanel(
                condition = "input.dist_orig_type == 'Sample Name' | input.dist_orig_type == 'Sample Alias'",
                selectInput("sim_sampleset", label="Sample sets:", choices=names(lox@count_matrixes)),
            ),
            conditionalPanel(
                condition = "input.dist_orig_type == 'Sample Name'",
                selectInput("sim_sample", label="Samples:", choices=names(lox@samples)),
            ),
            conditionalPanel(
                condition = "input.dist_orig_type == 'Sample Alias'",
                selectInput("sim_alias", label="Samples:", choices=c())
            ),
            conditionalPanel(
                condition = "input.dist_orig_type == 'Poisson'",
                numericInput("poisson_mean", label="Poisson distribution mean:", 3, min=1, max=14),
            )
        ),
    ),
    fluidRow(
        box(
            width = 4,
            title = "Simulated Samples",
            status = "danger",
            collapsible = TRUE,
            collapsed = TRUE,
            dataTableOutput("simulated_samples")
        )
    )
)

#########0#########1#########2#########3#########4#########5#########6#########7

# Dashboard
body <- dashboardBody(
    tabItems(
        tabItem(
            tabName = "import",
            box(
                width = 6,
                title = "Method 1: Upload a loxcode experiment object",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                fileInput("rds_file", "Choose an rds file:", multiple=TRUE, accept=c("rds")),
                actionButton("submit_rds", "Upload")
            ),
            box(
                width = 6,
                title = "Method 2: Upload samplesheet and fastq directory",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                textInput("name_exp", "Name of the loxcode experiment:", placeholder="Experiment Name"),
                fileInput("samplesheet", "Choose an xlsx file of samples", multiple=FALSE, accept=c("xlsx")),
                directoryInput('directory', label = 'Choose a fastq directory', value = '~'),
                #textInput("dir_input", "Choose a fastq directory"),
                actionButton("submit_fastq", "Upload")
            ),
            box(
                width = 12,
                title = "Loxcode Experiments",
                status = "danger",
                colour = NULL,
                solidHeader = FALSE,
                collapsible = FALSE,
                dataTableOutput("experiments_table"),
                actionButton("select_exp", "Select"),
                actionButton("merge_exp", "Merge selected"),
                actionButton("del_exp", "Delete"),
                actionButton("refresh", "Refresh Graph"),
                downloadButton("save_exp", "Download Current")
            ),
            box(
                width = 12,
                title = "Sample Sheet",
                status = "danger",
                collapsible = TRUE,
                collapsed = TRUE,
                wellPanel(DTOutput("samplesheet"))
            )
        ),

        simulationTab,

        tabItem(
            tabName = "codeset-create",
            fluidRow(
                # box(
                #   title = "View a code set",
                #   status = NULL,
                #   color = NULL,
                #   solidHeader = TRUE,
                #   collapsible = TRUE,
                #   selectInput("view_codeset", label="Choose code set:", choices=names(lox@code_sets)),
                #   actionButton("delete_codeset", "Delete Code Set")
                # ),

                #   box(width = 6,
                #   title = "View code sets",
                #   status = NULL,
                #   color = NULL,
                #   solidHeader = TRUE,
                #   collapsible = TRUE,
                #   selectInput("view_codeset", label="Choose a code set:", choices=c(),
                #   actionButton("delete_codeset", "Delete Code Set")
                # ),

                box(
                    width = 12,
                    title = "Edit code sets",
                    status = NULL,
                    color = "red",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    textInput("name_codeset", label="Name of new codeset:", placeholder="Code Set Name"),
                    fluidRow(
                        column(4, actionButton("create_codeset", "Create from Selection")),
                        column(4, actionButton("create_all_codeset", "Create from All")),
                        column(4, actionButton("rename_codeset", "Rename Current"))
                    )
                ),
                box(
                    width = 12,
                    title = "Code Set Table",
                    status = "danger",
                    solidHeader = FALSE,
                    wellPanel(DTOutput("codeset_table")),
                    verbatimTextOutput("selected_codes")
                )
            )
        ),

        tabItem(
            tabName = "sample-create",
            fluidRow(
                # box(
                #   width = 6,
                #   title = "View smaples",
                #   status = NULL,
                #   color = NULL,
                #   solidHeader = TRUE,
                #   collapsible = TRUE,
                #   selectInput("view_sample", label="Choose a sample:", choices=names(lox@count_matrixes)),
                #   actionButton("delete_sample", "Delete Sample Set")
                #   ),

                # box(
                #   width = 6,
                #   title = "View a sample set",
                #   status = NULL,
                #   color = NULL,
                #   solidHeader = TRUE,
                #   collapsible = TRUE,
                #   selectInput("view_sample", label="Choose sample set:",
                #               choices=names(lox@count_matrixes)),
                #   actionButton("delete_sample", "Delete Sample Set")
                # ),
                box(
                    width = 12,
                    title = "Edit sample sets",
                    status = NULL,
                    color = "red",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    textInput("name_sample", label="Name of new sample set:",
                              placeholder="Sample Set Name"),
                    fluidRow(
                        column(4, actionButton("create_sample", "Create from Selection")),
                        column(4, actionButton("create_all_sample", "Create from All")),
                        column(4, actionButton("rename_sample", "Rename Current"))
                    )
                ),
                box(
                    width = 12,
                    title = "Sample Set Table",
                    status = "danger",
                    solidHeader = FALSE,
                    wellPanel(dataTableOutput("sample_table")),
                    verbatimTextOutput("selected_samples"),
                ),
                box(
                    width = 6,
                    title = "Generate aliases",
                    status = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    collapsed = TRUE,
                    checkboxGroupInput("alias_parameters",
                                       "Choose alias parameters:", choices=c("")),
                    actionButton("generate_alias", "Generate")
                ),
                box(
                    width = 6,
                    title = "Collapse samples",
                    status = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    collapsed = TRUE,
                    textInput("collapse_name", label="Name of new sample set:",
                              placeholder="Sample Set Name"),
                    tags$strong("Collapse type: "),
                    fluidRow(
                        column(3, checkboxInput("collapse_union", "Union")),
                        column(3, checkboxInput("collapse_average", "Average"))
                    ),
                    checkboxGroupInput("collapse_parameters",
                                       "Choose parameters to collapse:", choices=c("Population","PCR.replicate","sample_name","mouse","Organ","tamoxifen","DOB","DOD","Barcode.induction.age","timepoint","Cre.type","fatq._path","experiment")),
                    actionButton("collapse_samples", "Collapse"),
                    actionButton("collapse_selection", "Collapse selected")
                )
            )
        ),

        tabItem(
            tabName = "codes-filter",
            box(
                width = 6,
                title = "View codes",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                selectInput("independent_samples",
                            label="Sample Set", choices=names(lox@count_matrixes)),
                selectInput("filter_codeset", label="Code Set", choices=c()),
            ),
            box(
                width = 6,
                title = "Filter parameters",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                column(6, sliderInput("filter_reps", label="Maximum allowed code repetitions", min=1, max=10, value=1)),
                column(6, sliderInput("filter_tolerance", label="Tolerance Level (%)", min=0.1, max=100, value=5, step=0.1)),
                textInput("filter_code_name", label="Name of new filtered code set:", placeholder="Code Set Name"),
                actionButton("create_filtered", "Create Filtered Code Set")
            ),
            box(
                width = 6,
                status = "danger",
                title = "Unfiltered codes",
                collapsible = TRUE,
                collapsed = TRUE,
                plotOutput("unfiltered_codes"),
            ),
            box(
                width = 6,
                status = "danger",
                title = "Filtered codes",
                collapsible = TRUE,
                collapsed = TRUE,
                plotOutput("filtered_codes"),
            ),
        ),

        tabItem(
            tabName = "stats-plot",
            fluidRow(
                box(
                    title = "View a sample set",
                    width = 4,
                    status = NULL,
                    color = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    selectInput("matrix_stats", "Sample:", choices = names(lox@count_matrixes))
                ),
                box(
                    title = "View a code set",
                    width = 4,
                    status = NULL,
                    color = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    selectInput("codeset_stats", "Codes:", choices = names(lox@code_sets))
                ),
                box(
                    title = "Configure the plot",
                    width = 4,
                    status = NULL,
                    color = NULL,
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    selectInput("labels_stats", "Plot Labels:",
                                choices = c("sample", "alias"), selected = "sample")
                ),
                box(
                    width = 12,
                    title = "Bar Code Table",
                    status = "danger",
                    solidHeader = FALSE,
                    wellPanel(dataTableOutput("barcode_table")),
                    actionButton("includeBartable", "Add to report")
                ),

                box(
                    width = 12,
                    title = "Statistics Plots",
                    status = "danger",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    tabsetPanel(
                        type = "tabs",
                        tabPanel(
                            "Reads",
                            plotlyOutput("read_plot"),
                            actionButton("includeRead", "Add to report")),
                        tabPanel(
                            "Size",
                            plotlyOutput("size_plot"),
                            checkboxInput("fill_size", "Fill"),
                            actionButton("includeSize", "Add to report")),
                        # actionButton("includeSizeNote", "Add Notes")),
                        tabPanel(
                            "Complexity",
                            plotlyOutput("complexity_plot"),
                            checkboxInput("fill_complexity", "Fill"),
                            actionButton("includeComplexity", "Add to report")),

                        tabPanel("Ratio", plotlyOutput("ratio_plot",height="1000",inline = FALSE),
                                 actionButton("includeRatio", "Add to report")),

                        tabPanel("Both", plotlyOutput("both_plot",height="1000",inline = FALSE),
                                 actionButton("includeBoth", "Add to report"))

                    )
                ),



                box(
                    width = 6,
                    title = "Sample Size Plot",
                    status = "danger",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    conditionalPanel(
                        condition = "input.labels_stats=='sample'",
                        selectInput("size_sample", "", choices=c()),
                    ),
                    conditionalPanel(
                        condition = "input.labels_stats=='alias'",
                        selectInput("size_alias", "", choices=c()),
                    ),
                    plotlyOutput("sample_size_plot"),
                    actionButton("includeSampleSize", "Add to report")
                    #actionButton("includeSampleSizeNote", "Add Notes")
                ),
                box(
                    width = 6,
                    title = "Sample Complexity Plot",
                    status = "danger",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    conditionalPanel(
                        condition = "input.labels_stats=='sample'",
                        selectInput("complexity_sample", "", choices=c()),
                    ),
                    conditionalPanel(
                        condition = "input.labels_stats=='alias'",
                        selectInput("complexity_alias", "", choices=c()),
                    ),
                    plotlyOutput("sample_complexity_plot"),
                    actionButton("includeSampleComplexity", "Add to report")
                    #actionButton("includeSampleComplexityNote", "Add Notes")
                )
            )
        ),

        tabItem(
            tabName = "heatmap-plot",
            box(
                width = 6,
                title = "View a sample set",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                selectInput("matrix_heat", "Sample:", choices = names(lox@count_matrixes))
            ),
            box(
                width = 6,
                title = "View a code set",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                selectInput("codeset_heat", "Codes:", choices = names(lox@code_sets))
            ),
            # box(
            #   width = 6,
            #   title = "Distance Plots",
            #   status = "danger",
            #   color = NULL,
            #   collapsible = TRUE,
            #   collapsed = TRUE,
            #   tabsetPanel(
            #     tabPanel(
            #       "Comparison Pie",
            #       plotlyOutput("sample_comparison_pie"),
            #       textInput("scale_pie", "Scale factor:", value=1)
            #     ),
            #     tabPanel(
            #       "Correlation Plot",
            #       selectInput("correlation_method", label="",
            #                   choices=c('pearson', 'kendall', 'spearman')),
            #       #selectInput("correlation_parameter", label="", choices=c()),
            #       # conditionalPanel(
            #       #   condition = "!is.null(input.correlation_parameter)",
            #       #   selectInput("correlation_data", label="", choices=c())
            #       # ),
            #       plotOutput("correlation_plot")
            #     )
            #   ),
            #
            # ),
            box(
                width = 6,
                title = "Configure the heatmap plot",
                status = "danger",
                color = NULL,
                collapsible = TRUE,
                collapsed = TRUE,
                selectInput("type_heat", "Plot Types:", choices = c("ggplot", "plotly")),
                selectInput("labels_heat", "Plot Labels:", choices=c("sample", "alias"), selected = "alias"),
                selectInput("split_by1", "Split By1:", choices=c("mouse","Organ","Population","PCR.repicate", "none"), selected = "none"),
                selectInput("split_by2", "Split By2:", choices=c("mouse","Organ","Population","PCR.repicate", "none"), selected = "none"),
                selectInput("clustering", "Clustering", choices=c("none", "row", "col", "both")),
                selectInput("agglomeration", "Agglomeration", choices=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid","binary"),selected="complete"),
                sliderInput("min_reads", label="Minimum # reads", min=1, max=100, value=1),
                sliderInput("max_repeats", label="Maximum # repeats", min=1, max=100, value=100),
                sliderInput("min_repeats", label="Minimum # repeats", min=1, max=100, value=1)
            ),
            box(
                width = 12,
                title = "Heat Map",
                status = "danger",
                color = NULL,
                collapsible = TRUE,
                collapsed = TRUE,
                conditionalPanel(

                    condition = "input.type_heat=='plotly'",
                    plotlyOutput("heatmap_plotly",height = 1600)
                ),
                conditionalPanel(

                    condition = "input.type_heat=='ggplot'",
                    plotOutput("heatmap_ggplot",height = 1600)
                ),
                actionButton("includeHeatmap", "Add to report")
                #actionButton("includeHeatmapNote", "Add Notes")
            ),
            box(
                width = 12,
                title = "Bubble Map",
                status = "danger",
                color = NULL,
                collapsible = TRUE,
                collapsed = TRUE,
                conditionalPanel(
                    condition = "input.type_heat=='ggplot'",
                    plotOutput("bubble_ggplot",height = 1600)
                ),
                conditionalPanel(
                    condition = "input.type_heat=='plotly'",
                    plotlyOutput("bubble_plotly",height = 1600)
                ),
                actionButton("includeBubble", "Add to report")
                #actionButton("includeBubbleNote", "Add Notes")
            ),
        ),
        tabItem(
            tabName = "saturation-plot",
            box(
                width = 6,
                title = "View a sample set",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                selectInput("sampleset_sat", "Sample sets:", choices = names(lox@count_matrixes)),
                selectInput("name_sat", "Plot labels", choices = c("sample", "alias"), selected = "alias")
            ),
            box(
                width = 6,
                title = "View plot",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                conditionalPanel(
                    condition = "input.name_sat=='sample'",
                    selectInput("sample_sat", "Samples:", choices = names(lox@samples))
                ),
                conditionalPanel(
                    condition = "input.name_sat=='alias'",
                    selectInput("alias_sat", "Samples:", choices = lox@alias[["all_samples"]]$alias)
                ),
                selectInput("codeset_sat", "Codes:", choices = names(lox@code_sets)),
                actionButton("add_sat", "Add new"),
                actionButton("remove_sat", "Remove"),
                actionButton("clear_sat", "Clear")
            ),
            box(
                width = 12,
                title = "Saturation Plot",
                status = "danger",
                collapsible = TRUE,
                collapsed = TRUE,
                plotlyOutput("saturation"),
                actionButton("includeSaturation", "Add to report"),
                #actionButton("includeSaturationNote", "Add Notes")
            )
        ),

        tabItem(
            tabName = "pair-plot",
            box(
                width = 6,
                title = "Choose your samples",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                selectInput("sampleset_pair", "Sample sets:", choices = names(lox@alias)),
                selectInput("codeset_pair", "Code Sets:", choices = names(lox@code_sets)),
                selectInput("name_pair", "Plot labels", choices = c("sample", "alias"), selected = "sample"),
                conditionalPanel(
                    condition = "input.name_pair=='sample'",
                    selectInput("sample1_pair", "Samples:", choices = c()),
                    selectInput("sample2_pair", "Samples:", choices = c())
                ),
                conditionalPanel(
                    condition = "input.name_pair=='alias'",
                    # selectInput("sample1_pair", "Samples:", choices = c()),
                    # selectInput("sample2_pair", "Samples:", choices = c())
                    selectInput("alias1_pair", "Samples:", choices = c()),
                    selectInput("alias2_pair", "Samples:", choices = c())
                ),
                actionButton("add_pair", "Add new plot"),
                actionButton("remove_pair", "Remove plot"),
                actionButton("clear_pair", "Clear plots"),
                actionButton("add_all", "Add all combinations")
            ),
            box(
                width = 6,
                title = "Configure the plot",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                selectInput("type_pair", "Plot type", choices = c("ggplot", "plotly")),
                selectInput("colour_pair", "Colour by:", choices = c("size", "dist_orig", "firstread")),
                sliderInput("complexity_slider_pair", "Filter Distance Origin Range:", min = 0, max = 0, value = c(0,0)),
                sliderInput("size_slider_pair", "Filter Size Range:", min = 0, max = 0, value = c(0,0)),
                sliderInput("firstread_slider_pair", "Filter Firstreads Range", min = 0, max = 0, value = c(0,0))
            ),
            box(
                width = 12,
                #height = 10,
                #height = 1000,
                title = "Pair Comparison Plot",
                status = "danger",
                collapsible = TRUE,
                collapsed = TRUE,
                conditionalPanel(
                    condition = "input.type_pair=='ggplot'",
                    plotOutput("pair_ggplot"),
                ),
                conditionalPanel(
                    condition = "input.type_pair=='plotly'",
                    plotlyOutput("pair_plotly"),
                ),
                #actionButton("includePair", "Add to report")
                #actionButton("includePairNote", "Add Notes")
            )
        ),
        tabItem(
            tabName = "report",
            box(
                width = 12,
                title = "Generate report",
                solidHeader = TRUE,
                radioButtons("format", "Document format", c("PDF", "HTML", "Word"), inline = TRUE),
                downloadButton("downloadReport")
            ),
            box(
                width = 12,
                title = "Experiment Details",
                textAreaInput("name","Researcher Name :",value ="Your Name"),
                textAreaInput("expname","Experiment Name :",value ="Experiment 1"),
                textAreaInput("comment","Additional Comments :",value ="Comment"),
                actionButton("detail_submit", "Submit")
            ),
            box(
                width = 12,
                status = "danger",
                title = "Report Components",
                mainPanel(code("Double click 'Annotation' section to edit comment"),
                          p("\n")),
                dataTableOutput("components_table"),
                #actionButton("components_table_cell_edit", "edit"),
                actionButton("remove_component", "Remove")
            ),
            box(
                width = 12,
                title = "Select Loxcode Experiments",
                status = "danger",
                colour = NULL,
                solidHeader = FALSE,
                collapsible = FALSE,
                dataTableOutput("experiments_table1"),
                actionButton("select_exp1", "Select")
            ),
            box(
                title = "Select a code set",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                selectInput("view_codeset1", label="Choose code set:", choices=names(lox@code_sets))
            ),
            box(
                width = 6,
                title = "View a sample set",
                status = NULL,
                color = NULL,
                solidHeader = TRUE,
                collapsible = TRUE,
                selectInput("view_sample1", label="Choose sample set:",
                            choices=names(lox@count_matrixes)),
            ),
            box(
                width = 12,
                status = "danger",
                title = "Report Components",
                checkboxInput("includeSize1", label="Size Plot", value = FALSE, width = 220),
                checkboxInput("includeComplexity1", label="Complexity Plot", value = FALSE, width = 220),
                checkboxInput("includeRatio1", label="Ratio Plot", value = FALSE, width = 220),
                checkboxInput("includeBoth1", label="Both Plot", value = FALSE, width = 220),
                checkboxInput("includeSampleSize1", label="Sample Size Plot", value = FALSE, width = 220),
                checkboxInput("includeSampleComplexity1", label="Sample Complexity Plot", value = FALSE, width = 220),
                checkboxInput("includeHeatmap1", label="Heatmap Plot", value = FALSE, width = 220),
                checkboxInput("includeBubble1", label="Bubble Plot", value = FALSE, width = 220),
                checkboxInput("includeSaturation1", label="Saturation Plot", value = FALSE, width = 220)
                #actionButton("includeSize1", "Add Size Plot")
            )
        ),
        tabItem(
            tabName = "log",
            box(
                width = 12,
                title = "Actions",
                solidHeader = TRUE,
                actionButton("restart", "Restart", icon=icon("refresh")),
                bookmarkButton(),
                downloadButton("downloadLog"),
            ),
            box(
                width = 12,
                status = "danger",
                title = "Activity Log",
                collapsible = TRUE,
                collapsed = TRUE,
                dataTableOutput("log_table")
            )
        )
    )
)

# Define the js method that resets the page
jsResetCode <- "shinyjs.reset1 = function() {history.go(0)}"
enableBookmarking(store = "server")

ui <- function(request) {
    fluidPage(
        useShinyalert(),
        useShinyjs(),
        basicPage(),
        extendShinyjs(text=jsResetCode, functions=c("reset1")),
        dashboardPage(header, sidebar, body, skin = "red")
    )
}
