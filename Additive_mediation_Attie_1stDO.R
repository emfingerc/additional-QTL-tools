#' Version 0.1 of the additive mediation shiny tool that allows mediation for
#' using intermediate and bmediatR packages
#' 
#' intermediate: https://github.com/byandell/intermediate
#' bmediatR: https://github.com/wesleycrouse/bmediatR
#' 
#' In addition to the relevant package citations, parts of the datatable code 
#' followed formats shown on 
#' https://clarewest.github.io/blog/post/making-tables-shiny/
#' 
#' Also note that the select covariates sub-menu includes operators so that the 
#' final matrix can be generated using as.formula() and model.matrix()
#' E.g. when one selects "Sex" "+" and "Generation" using the menu, in that order,
#' the final matrix is generated using the formula '~Sex+Generation'


# load libraries=============================================================================
require("rstudioapi")
require("bmediatR")
require("dplyr")
require("stringr")
require("tidyverse")
require("BiocManager")
require("intermediate")
require("ggplot2")
require("cowplot")
require("qtl2")
require("grid")
require("ggrepel")
require("gridGraphics")
require("ggpubr")
require("shiny")
require("shinyFiles")
require("bslib")
require("spsComps")
require("DT")
require("shinyjs")

# set shiny options==========================================================================
options(shiny.maxRequestSize = 20000*1024^2)  # Increase to 20GB needed for Genoprobs, etc

# load data==================================================================================
# datasets
# note that this uses a number of datafiles that for the moment I've stored locally. This will
# change in the future
# these are the RDS objects formatted for the QTL viewer (Churchill group format) 
liver_RNA <- readRDS("~/Liver.RDS")
islet_RNA <- readRDS("~/Islet.RDS")
heart_RNA <- readRDS("~Heart.RDS")
skel_muscle_RNA <- readRDS("~/SkeletalMuscle.RDS")
clinical <- readRDS("~/Clinical_Phenotypes_V11.RDS")
adipose_RNA <- readRDS("~/Adipose.RDS")
islet_proteins <- readRDS("~dataset.islet.proteins.RDS")
ex_vivo <- readRDS("~Ex_vivo.RDS")

# files for QTL scans
allele_probabilities <- readRDS("~/attie_DO500_genoprobs_v5.rds")
kinship <- readRDS("~/Calculated_K.rds")
map <- readRDS("~/grid_pmap.rds")
markers <- read.csv("~/attie_v10_markers")
chrom_sep <- read.csv("~/chromosomal_sep.csv")

# files to mediate with
# these are derived from the QTL2 RDS objects 
mediation_phenotypes <- readRDS("~/Phenotypes_V11_combined.rds")
mediation_tissue_specific <- readRDS("~/Combined_v11_RNA_objects.RDS")
mediator_proteins <- readRDS("~/New_islet_proteomics.rds")

# diagrams for the mediation models
diagrams <- readRDS("~/corrected_mediation_diagrams_v2.rds")

# listing objects for the shiny app to call
mediator_list_object <- c("mediator_proteins", "mediation_tissue_specific", "mediation_phenotypes")
names(mediator_list_object) <- c("Proteins","RNA","Phenotypes")
data_list <- c("liver_RNA", "islet_RNA", "heart_RNA", "skel_muscle_RNA", "adipose_RNA", "clinical", "ex_vivo")
names(data_list) <- c("Liver RNA", "Islet RNA", "Heart RNA", "Skeletal muscle RNA", "Adipose RNA", "Clinical traits", "Ex vivo")
categor_covar <- c("Sex","Wave","Batch","GenLit","Diet","DietName", "BirthDate", "Generation")
current_env <- environment()
completed_mediation <- c()
covariates <- "~"

# set parameters ===================================================================================
#DO_set
DO_set <- "first"
# select RNA type
RNA_type <- "genes"
# peak type
peak_type <- "additive"
# ln_prior_c
ln_prior_c <- "complete"
#set cores
ncors = 1

#set functions======================================================================================
# the mediation function itself is a file
source("~/Additive_mediator_script_for_shiny.R")

# reformat the phenotypes
pheno_reformat <- function(selected_dataset){
  # get the data to be mediated
  if (selected_dataset$datatype == "pheno"){
    pheno <- data.frame(selected_dataset$data)
  }
  if (selected_dataset$datatype != "pheno"){
    pheno <- data.frame(selected_dataset$data$rz)
  }
  pheno$Mouse <- rownames(pheno)
  # note- the RDS objects consider the covariate factors separately so these must be 
  # gathered from the sample information in the dataset
  ms_info <- data.frame(selected_dataset$annot.samples)
  colnames(ms_info)[1]<-"Mouse"
  rownames(ms_info)<-ms_info$Mouse
  pheno <- ms_info %>% inner_join(., pheno, by=c("Mouse"="Mouse"))
  rownames(pheno)<-pheno$Mouse
  return(pheno)
}

# generate QTL peak lists
make_QTL_list <- function(selected_dataset, chosen_trait){
  # populate lists
  QTL_list <- data.frame(selected_dataset$lod.peaks[["additive"]])
  # make sure the columns are named appropriately
  if (selected_dataset$datatype == "mrna"){
    QTL_list$lodcolumn <- QTL_list[,grep(pattern="^gene",x=colnames(QTL_list))]
  }
  if (selected_dataset$datatype == "pheno"){
    QTL_list$lodcolumn <- QTL_list[grep(pattern="^data",x=colnames(QTL_list))]
  }
  if (selected_dataset$datatype == "protein"){
    QTL_list$lodcolumn <- QTL_list[grep(pattern="^protein",x=colnames(QTL_list))]
  }
  # select the phenotype to scan within the peaks
  pheno_list <- data.frame(unique(QTL_list$lodcolumn))
  colnames(pheno_list)<-"Phenotype"
  # option selected is 'chosen_type'
  pheno_list <- subset(pheno_list, Phenotype == chosen_trait)
  QTL_list <- subset(QTL_list, lodcolumn == chosen_trait)
  # add the necessary information on the peak location
  if (DO_set == "first"){
    QTL_list <- markers[c("marker.id","chr","pos")] %>% inner_join(., QTL_list, by=c("marker.id"="marker.id"))
  }
  if (DO_set == "second"){
    markers <- markers %>% mutate(pos = bp_mm10/(10^6))
    QTL_list <- markers[c("marker.id","chr","pos")] %>% inner_join(., QTL_list, by=c("marker.id"="marker.id"))
  }
  return(QTL_list)
}

# generate the mediation list object
make_mediator_list <- function(mediator_object_choice, mediator_compartment){
  chosen_set <- .GlobalEnv[[mediator_list_object[[mediator_object_choice]]]]
  if (mediator_object_choice=="RNA"){
    compartments <- names(chosen_set$RNA_values)
    major_object <- "RNA_values"
    mediator_type <- "RNA"
    major_info_object <- "RNA_info"
    overall_mediator_type <- "RNA/Protein"
    object_name <- mediator_compartment
    tissue_compartment <- mediator_compartment
    info_object <- mediator_compartment
  }
  if (mediator_object_choice=="Proteins"){
    compartments <- names(chosen_set$Protein_values)
    major_object <- "Protein_values"
    mediator_type <- "Protein"
    major_info_object <- "Protein_info"
    overall_mediator_type <- "RNA/Protein"
    object_name <- mediator_compartment
    tissue_compartment <- mediator_compartment
    info_object <- mediator_compartment
  }
  if (mediator_object_choice=="Phenotypes"){
    compartments <- names(chosen_set$Data)
    major_object <- "Data"
    mediator_type <- "Phenotypes"
    major_info_object <- "Info"
    overall_mediator_type <- "Phenotypes"
    object_name <- mediator_compartment
    tissue_compartment <- mediator_compartment
    info_object <- mediator_compartment
  }
  mediator_list <- data.frame(matrix(nrow=1, ncol=6))
  colnames(mediator_list) <- c("Major_object","Object.names", "Tissue.compartment", "Type", "Major_info_object", "Info_object")
  mediator_list$Major_object[1] <- major_object
  mediator_list$Object.names[1] <- object_name
  mediator_list$Tissue.compartment[1] <- tissue_compartment
  mediator_list$Type[1] <- mediator_type
  mediator_list$Major_info_object[1] <- major_info_object
  mediator_list$Info_object[1] <- info_object
  return(mediator_list)
}

# set UI======================================================================================================
ui <- page_sidebar(
  # sets the page
  # sets the title
  titlePanel("Additive mediation tool"),
  
  # sets sidebar
  sidebar = sidebar(
    id = "side_panel",
    helpText(
      "Select your data to mediate, what to mediate with, and options."
    ),
    
                        selectInput(
                          "selected_dataset",
                          label = "Choose a dataset to display",
                          choices = names(data_list),
                          selected = NULL
                        ),
                        sliderInput(
                          "LOD_thr",
                          label = "LOD threshold for evaluation",
                          min = 4,
                          max = 20,
                          value = 6,
                          round = TRUE
                        ),
                          selectizeInput(
                          inputId = "mediator_object_choice",
                          label = "Choose the dataset to mediate with (RNA, protein, or other phenotypes.",
                          choices = names(mediator_list_object),
                          multiple = FALSE,
                          options = list(
                            placeholder = 'Search...'
                          )
                        ),
    selectizeInput(
      "mediator_compartment",
      label = "Which compartment/dataset (e.g. Liver)? Note that you can only use one at a time. This means must delete whatever is selected by default if it isn't what you want before you type your choice.",
      choices = NULL,
      multiple = TRUE,
      options = list(
        placeholder = 'Search...'
      )
    ),  
    selectizeInput(
      inputId = "which_covar",
      label = "Choose covariates. Note that if you want multiple, also choose the operators (e.g. +, *, etc) between them. You cannot choose operators by themselves or multiple traits without operators between them.",
      choices = NULL,
      multiple = TRUE,
      options = list(
        placeholder = 'Search...'
      )),
   
    selectizeInput(
      inputId = "which_trait",
      label = "Choose the trait to mediate",
      choices = NULL,
      multiple = TRUE,
      options = list(
        placeholder = 'Search...'
      )),
    
    selectizeInput(
      inputId = "which_peak",
      label = "Choose the peak to mediate. Select it by its marker id from the table",
      choices = NULL,
      multiple = TRUE,
      options = list(
        placeholder = 'Search...'
      )),
    sliderInput(
      "window",
      label = "What window for analysis (Mbp)? 0 means all genome",
      min = 0,
      max = 150,
      value = 4,
      round = TRUE
    ),
    sliderInput(
      "lod_drop_thr",
      label = "What fractional LOD drop is considered of interest?",
      min = 0.3,
      max = 1,
      value = 0.4,
      round = TRUE
    ),
    sliderInput(
      "probability_effect_thr",
      label = "What probability is considered of interest?",
      min = 0.4,
      max = 1,
      value = 0.4,
      round = TRUE
    ),
    pre(id = "console"),
    actionButton("run", "Run"),
    actionButton("reset_input", "Reset inputs")
  ),
  # sets main panel
  mainPanel(
    # set scrollbar
    div(style = "overflow-y: scroll;"),
    position = "right",
    #set starting layout
    card(
      card_header("Selected trait & covariates used"),
      verbatimTextOutput("selected_trait"),
      verbatimTextOutput("covariates"),
      verbatimTextOutput("covariate_set")
      
    ),
    card(
      card_header("Peaks"),
      DT::dataTableOutput("peaks")
    ),
    card(
      card_header("QTL rescan"),
      plotOutput("QTL_plot")
    ),
    card(
      card_header("LOD drop"),
      plotOutput("LOD_drop_plot")
    ),
    card(
      card_header("Complete mediation"),
     plotOutput("complete_mediation")
    ),
    card(
      card_header("Colocal effects"),
      plotOutput("colocal")
    ),
    card(
      card_header("Results"),
      DT::dataTableOutput("mediation_table")
    ),
    card(
      card_header("Console"),
      "Activity console",
      shinyjs::useShinyjs(),
      textOutput("text")
    )  
  ),
)

# set server==================================================================================================
server <- function(input, output, session) {
 
    
 # monitor the console
 #  observeEvent(input$run, {
  #  withCallingHandlers({
  #   shinyjs::html("text", "")
  #  shinyCatch({message("a message")}, prefix = '')
  #  },
  #  message = function(m) {
  #    shinyjs::html(id = "text", html = m$message, add = TRUE)
  # })
  # })
  
  # create the mediator list
  observeEvent(input$selected_dataset, {
    # update expression
    req(input$selected_dataset)
    req(input$mediator_object_choice)
    updateSelectizeInput(session,
                         "mediator_object_choice", 
                         choices = names(mediator_list_object),
                         options = list(maxItems = 1,
                                        maxOptions = 5),
                         server = TRUE
    )
  })
  
  output$selected_dataset <- renderText({
    paste0(input$selected_dataset)
  })
  
  # create the compartment list
  observeEvent(input$mediator_object_choice, {
    # update expression
    req(input$mediator_object_choice)
    chosen_set <- .GlobalEnv[[mediator_list_object[[input$mediator_object_choice]]]]
    if (input$mediator_object_choice=="RNA"){
      compartments <- names(chosen_set$RNA_values)
    }
    if (input$mediator_object_choice=="Proteins"){
      compartments <- names(chosen_set$Protein_values)
    }
    if (input$mediator_object_choice=="Phenotypes"){
      compartments <- names(chosen_set$Data)
    }
    updateSelectizeInput(session,
                      "mediator_compartment", 
                      choices = compartments,
                      options = list(maxItems = 1,
                                     maxOptions = 5),
                      server = TRUE
    )
  })
  
  output$mediator_object_choice <- renderText({
    paste0(input$mediator_object_choice)
  })
  
  # the reformatting of the selected traits
  pheno <- reactive({
    chosen_data <- .GlobalEnv[[data_list[[input$selected_dataset]]]]
    phenotypes <- pheno_reformat(chosen_data)
    return(phenotypes)
  })
  
  # select the phenotypes
   myCovar <- reactive({
     covariates <- ""
     covars <- input$which_covar
     if (is.null(covars)) {
       return("")
     }
       covariates <- paste0(covariates,input$which_covar)
       assign("covar_statement",covariates, envir = .GlobalEnv)
   })
   
  # update the selection for the traits of interest
     observeEvent(input$selected_dataset, {
       # update expression
         req(input$selected_dataset)
         updateSelectizeInput(session,
                            "which_covar", 
                            choices = c(colnames(pheno()),"+","*","-","/"),
                            options = list(maxItems = 5,
                                           maxOptions = 5),
                            server = TRUE
                            )
     })
     
     # initial message
     output$covariates <- renderText({
       paste0(myCovar())
     })
     output$covariate_set <- renderPrint({
       paste0(myCovar())
     })
     
     # create the covariate matrix
     observeEvent(input$run,{
       covariates <- as.formula(paste0("~",paste(covar_statement,collapse = "")))
       phenotypes <- pheno()
       if (any(covar_statement %in% c("+","*","-","/"))){
         covar_statement <- covar_statement[-which((covar_statement %in% c("+","*","-","/"))==TRUE)]
       }
       if (any(covar_statement %in% categor_covar)){
         for (i in 1:nrow(phenotypes)){
           phenotypes[i,covar_statement[which((covar_statement %in% categor_covar)==TRUE)]] <- as.character(phenotypes[i,covar_statement[which((covar_statement %in% categor_covar)==TRUE)]])
         }
       }
       covar_matrix <- model.matrix(covariates, data=phenotypes)[,-1]
       #covar_matrix <- data.frame(covar_matrix)
       #colnames(covar_matrix) <- covar_statement
       #covar_matrix <- as.matrix(covar_matrix)
       if(length(covar_statement)==1){
         covar_matrix <- data.frame(t(t(covar_matrix)))
         #colnames(covar_matrix) <- covar_statement
         covar_matrix <- as.matrix(covar_matrix)
       }
       assign("covar_matrix", covar_matrix, envir = .GlobalEnv)
       assign("covariate_formula", covariates, envir = .GlobalEnv)
       assign("length_covar_resp", length(covar_statement), envir = .GlobalEnv)
       assign("phenotype_obj", phenotypes, envir = .GlobalEnv)
     })
  
    # create the trait list
     observeEvent(input$selected_dataset, {
       # update expression
       req(input$selected_dataset)
       updateSelectizeInput(session,
                            "which_trait", 
                            choices = colnames(pheno()),
                            options = list(maxItems = 1,
                                           maxOptions = 5),
                            server = TRUE
       )
     })
     
     output$selected_trait <- renderText({
       paste0(input$which_trait)
     })
   
     # populate a list of peaks
     output$peaks <- renderDT({DT::datatable(
                 make_QTL_list(.GlobalEnv[[data_list[[input$selected_dataset]]]], input$which_trait),
                 options = list(paging = TRUE,    ## paginate the output
                                pageLength = 5,  ## number of rows to output for each page
                                scrollX = TRUE,   ## enable scrolling on X axis
                                scrollY = TRUE,   ## enable scrolling on Y axis
                                autoWidth = TRUE, ## use smart column width handling
                                server = TRUE,   ## use client-side processing
                                dom = 'Bfrtip',
                                buttons = c('csv', 'excel'),
                                columnDefs = list(list(targets = '_all', className = 'dt-center'),
                                                  list(targets = c(0, 8, 9), visible = FALSE))
                 ),
                 extensions = 'Buttons',
                 selection = 'single', ## enable selection of a single row
                 filter = 'bottom',              ## include column filters at the bottom
                 rownames = TRUE                ## do show row numbers/names
     )
     })
     
     # update the selection for the peaks of interest
     observeEvent(input$which_trait, {
       # update expression
       req(input$selected_dataset)
       updateSelectizeInput(session,
                            "which_peak", 
                            choices = make_QTL_list(.GlobalEnv[[data_list[[input$selected_dataset]]]], input$which_trait)$marker.id,
                            options = list(maxItems = 1,
                                           maxOptions = 5),
                            server = TRUE
       )
     })
     output$selected_peak <- renderText({
       paste0(input$which_peak)
     })
     
     
     
     # run the mediation analysis
     observeEvent(input$run, {
       # check availability of data
       req(input$selected_dataset)
       req(input$which_trait)
       req(input$which_peak)
       # process the data
       pheno_list <- data.frame(matrix(nrow=1, ncol=1))
       colnames(pheno_list) <- "Phenotype"
       pheno_list$Phenotype[1]<- input$which_trait        
       chosen_mediators <- .GlobalEnv[[mediator_list_object[[input$mediator_object_choice]]]]
       mediator_list <- make_mediator_list(input$mediator_object_choice, input$mediator_compartment)
       if (input$mediator_object_choice=="RNA"){
         mediator_type <- "RNA/Protein"
       }
       if (input$mediator_object_choice=="Proteins"){
         mediator_type <- "RNA/Protein"
       }
       if (input$mediator_object_choice=="Phenotypes"){
         mediator_type <- "Phenotypes"
       }
       if (input$window == 0){
         usewindow <- FALSE
       }
       if (input$window > 0){
         usewindow <- TRUE
         window <- as.numeric(input$window)
       }
       phenotypes <- pheno()
       if (any(covar_statement) %in% categor_covar){
       for (i in 1:nrow(phenotypes)){
           phenotypes[i,covar_statement[which((covar_statement %in% categor_covar)==TRUE)]] <- as.character(phenotypes[i,covar_statement[which((covar_statement %in% categor_covar)==TRUE)]])
       }
       }
       QTL_list <- make_QTL_list(.GlobalEnv[[data_list[[input$selected_dataset]]]], input$which_trait)
       QTL_list <- subset(QTL_list, marker.id == input$which_peak)
       completed_mediation <- Emfinger_mediation(
         wr_dir = getwd(),
         K = kinship,
         genoprobs = allele_probabilities,
         map = map,
         pheno = phenotypes,
         covar = covar_matrix,
         test_these = pheno_list,
         mediator = chosen_mediators,
         cc_colors = "new",
         QTL_list = QTL_list,
         mediator_list = mediator_list,
         plot_mediation = TRUE,
         save_mediation = FALSE,
         show_fits = FALSE,
         lod_drop_thr = input$lod_drop_thr,
         ln_prior_c = ln_prior_c,
         diagrams = diagrams,
         save_rescan=FALSE,
         probability_effect_thr = input$probability_effect_thr,
         save_plot_obj =FALSE,
         covar_pheno = NULL,
         covar_mediator = NULL,
         weights_all = NULL,
         weights_pheno = NULL,
         weights_mediator = NULL,
         pval = TRUE,
         pval_model_type = "mediation",
         LOD_thr = input$LOD_thr,
         mediator_type = mediator_type,
         markers = markers,
         chr_breaks = chrom_sep,
         DO_set = DO_set,
         pdf_name_prompt = FALSE,
         ncors = ncors,
         RNA_type = RNA_type,
         usewindow = usewindow,
         window = window,
         save_PDF_file = FALSE,
         perm_p = FALSE,
         combine_graphs = FALSE
       )
       assign("completed_mediation", completed_mediation, envir = .GlobalEnv)
       
       # populate a list of peaks
       output$mediation_table <- renderDT({DT::datatable(
         completed_mediation$tabular_results,
         options = list(paging = TRUE,    ## paginate the output
                        pageLength = 5,  ## number of rows to output for each page
                        scrollX = TRUE,   ## enable scrolling on X axis
                        scrollY = TRUE,   ## enable scrolling on Y axis
                        autoWidth = TRUE, ## use smart column width handling
                        server = TRUE,   ## use client-side processing
                        dom = 'Bfrtip',
                        buttons = c('csv', 'excel'),
                        columnDefs = list(list(targets = '_all', className = 'dt-center'),
                                          list(targets = c(0, 8, 9), visible = FALSE))
         ),
         extensions = 'Buttons',
         selection = 'single', ## enable selection of a single row
         filter = 'bottom',              ## include column filters at the bottom
         rownames = FALSE                ## don't show row numbers/names
       )
       })
       
       #QTL plot
       output$QTL_plot <- renderPlot({completed_mediation$QTL_plot})
       
       #LOD drop plot
       output$LOD_drop_plot <- renderPlot({completed_mediation$LOD_drop_plot})
       
       #complete mediation plot
       output$complete_mediation <- renderPlot({completed_mediation$bmediatR_plots$complete})
       
       #colocal plot
       output$colocal <- renderPlot({completed_mediation$bmediatR_plots$colocal})
     })
     
     # reset events
     observeEvent(input$reset_input, {
       shinyjs::reset("side_panel")
       #QTL plot
       output$QTL_plot <- NULL
       
       #LOD drop plot
       output$LOD_drop_plot <- NULL
       
       #complete mediation plot
       output$complete_mediation <- NULL
       
       #colocal plot
       output$colocal <- NULL
       
       #mediation table
       output$mediation_table <- NULL
     })
}


# set app ====================================================================================================
shiny::shinyApp(ui = ui, server = server)
