# LoxcodeR
###### Web application to visualise and analyse different aspects of a Loxcode.
###### Manipulate and extract meaningful insights from experiment data.
###### Pdf and Web page based report creation to present data analysed in the application.

## Motivation
Main objective of the app is to perform data manipulation and data visualization in a more efficient and user friendly environment. The app can perform some general visualization with just the click of a few buttons. Minimal to no coding expertise required to use the app. The app is an end to end solution. It can load up loxcode raw data and convert it into a rds data format. It can also export all the visualization created by the user to a report format for easy presentation.

## Technology Used
- R programming
- C programming
- RMD
- RDS Dataset
- Shiny Front-end

## Installation
- For Ubuntu 20.04.02 or debian user download dependencies using sudo apt-get install from terminal:
	- `sudo apt-get install libcurl4-openssl-dev`
	- `sudo apt-get install libssl-dev`
	- `sudo apt-get install libxml2-dev`
	- `sudo apt-get install libgit2-dev`
	- `sudo aptitude install libgdal-dev`
	- `sudo apt-get install texlive-latex-extra`
	- `sudo apt-get install texlive-fonts-extra`
- For MacOS does not need to do the above installation as it already had those functions.
- open the file loxcoder.Rproj
- open the file server.R from LoxcodeR_app folder
- go to build (next to Environment and History), click "install and Restart" button to install the loxcoder package
- Install the following packages:
    - install.packages(c('devtools','pals','viridis','ggplot2','tidyr','ggbeeswarm','dygraphs','shiny','dplyr','htmlwidgets','digest','bit','flexdashboard','ggthemes','highcharter','scatterD3','comprehenr','VennDiagram','scatterpie','tidyverse','reshape','VennDiagram','scatterpie','tidyverse','reshape'))
    - tinytex::install_tinytex()
    - devtools::install_github("rstudio/crosstalk")
    - devtools::install_github("jcheng5/d3scatter")
    - devtools::install_github("hadley/shinySignals")
- Click on RunApp or type runApp('LoxcodeR_app') into the console to run the app

## Data
App can take raw data in the format 
Raw data is converted to rds format which is easier to store and manipulate by the app.

## Data Manipulation
All data manipulation is done on the rds data object. Code set and the sample set can be both manipulated to select only the desired subset from the experiment. 
For sample set, collapse on particular attribute feature is also available. 
For code set, we can filter it based on code repetition and Tolerance level as well.

## Analysis
We have many different types of plot that can be made for the experiment:
- Statistic Plot
    - Read Plot
        This plot shows the number of read for each sample in the experiment.<br/>
![Read Plot](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/ReadPlot.png)
    - Size Plot
        This plot compares the count of different loxcode size for each sample in the experiment.<br/>
![Size Plot](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/SizePlot.png)
    - Complexity Plot
        This plot compares the count of distance from origin for each sample in the experiment.<br/>
![Complexity Plot](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/ComplexityPlot.png)        
    - Ratio Plot
        This plots the ratio of complexity to size for each sample in the experiment.<br/>
![Ratio Plot](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/RatioPlot.png)        
    - Both Plot
        This plot compares the complexity and size based on the number of reads for each sample in the experiment.<br/>
      
- Sample Size Plot
    This plots diversity of valid and invalid loxcode for each size in 1 sample.<br/>
  ![Sample Size Plot](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/SampleSizePlot.png)
- Sample Complexity Plot 
    This plots the diversity of complexity versus the distance from origin for each size loxcode in 1 sample.<br/>
  ![Sample Complexity Plot](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/SampleComplexityPlot.png)
- Heat Map
    This plots a heatmap of each sample in the experiment for the log of loxcode count. This plot can be split by 2 factorts by configuring the heatmap.<br/>
![Heat Map](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/HeatMap.png)
- Bubble Map
    This plots a bubbple plot of each sample in the experiment for the log of loxcode count. This plot can be split by 2 factorts by configuring the heatmap.<br/>
![Bubble Plot](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/BubblePlot.png)
- Saturation Plot
    (Need to add label)<br/>
![Saturation Plot](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/SaturationPlot.png)
- Pair Comparision Plot
    It compares 2 samples based on the distribution size or distance from origin or firstread.<br/>
![PairComparision Plot](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/PairComparisionPlot.png)
![Pair Comparision Plot 2](https://github.com/tomsergeweber/LoxCodeR2022/blob/master/Docs/PairComparisionPlot2.png)
## Future Work
