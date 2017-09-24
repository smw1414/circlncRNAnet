#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# circRNA 
library(shiny)
library(shinyjs)
library(ggplot2) #devtools::install_github('hadley/ggplot2')
library(googleVis)
library(data.table)
require(visNetwork)
library(magrittr)
library(DT)
library(plotly)
data=fread("../output/linc_coexp_pairs.txt",header=T,sep="\t",na.strings = c(".",NA))
data=data[cor!=""]
selected_cols= c("co_exp_gene","cor","cor_p","distance","co_exp_gene_chr","co_exp_gene_start","co_exp_gene_end")
qsum=fread("../output/q_cor_summary.txt",header=T,sep="\t")
colnames(qsum)=c("Query","cor > 0.5 (% of coexpressed genes)","pvalue < 0.05 (% of coexpressed genes)")
jobid= gsub("/data/lancer_jobs/|/coexpc","",getwd())

run_circos="./circos_zoom.r"
run_gene_enrichment <- "./gene_enrichment.r"
run_scatterplot <- "./scatterplot.r"
run_heatmap <- "./heatmap.R"
run_triple_network_circRNA_RBP_2step <- "./triple_network_circRNA_RBP_2step.R"
run_triple_network_circRNA_sponge_2step <- "./triple_network_circRNA_sponge_2step.R"
run_triple_network_lncRNA_RBP_2step <- "./triple_network_lncRNA_RBP_2step.R"
run_triple_network_lncRNA_sponge_2step <- "./triple_network_lncRNA_sponge_2step.R"
run_randomization <- "./randomization.R"

#data=data[,.(gene,log2FoldChange,pvalue,Gene_full_name,Gene_old_names,Gene_other_names,Coexpressed_gene="Link")][order(-abs(log2FoldChange))]
# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
   useShinyjs(),
   # Application title
   #titlePanel(h4("coexpressed genes")),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        #verbatimTextOutput('out1'),
        selectInput("qgene","Query gene",unique(data$circRNA),selectize = TRUE,selected =  NULL),hr(),
        helpText("Co-expressed gene"),
      	sliderInput("corr", "Pearson's correlation coefficient (default: |cor| > 0.5)",
                    min = -1, max = 1, value = c(-0.5, 0.5),step = 0.05),
        checkboxInput('exclude', 'Exclude range', TRUE),hr(),
        hr(),
        tags$p(paste("Number of co-expressed gene(s):")),
        verbatimTextOutput("number_coexp"),
       # numericInput("corr", "Pearson's correlation coefficient > 0.5 (default):", 0.5, step = 0.1, max = 1),
       # numericInput("pvalue", "P-value < 0.05 (default):", 0.05,step = 0.01,max=1,min=0),hr(),
       
#####################################################
        tabsetPanel(
          id = 'summary',
          tabPanel('Summary',
          #      column(12, br(),htmlOutput("qsumchart"))
		column(12, br(),plotlyOutput("qsumchart"))
          ),
          tabPanel('BoxPlot',
                   column(12,br(),h5("Overall Transcriptome Profile"),imageOutput("boxplot"))
          ),
          tabPanel('ScatterPlot',
                   column(12,br(),br(),
                          uiOutput("co_exp_gene"),
                           imageOutput("scatterplot")
                          #plotOutput("scatterplot")
                          )
          ),
		tabPanel('Randomized correlation',
		         value = "rand_panel",
		         column(
		           12,
		           br(),
		           hidden(div(
		             id = "loading-rand",
		             img(src = "../../../images/loading-1.gif")
		           )),
		           strong("Obs: Observation, Rand: Randomization"),
		           plotlyOutput('randomization')#, height = "300px")
		         ))
        )
#####################################################

      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         tabsetPanel(
           id = 'coexp',
           tabPanel(
             'Co-expressed Genes',
             fluidRow(column(
               12,
               br(),
               checkboxInput("selcols", tags$b("Modify columns:")),
               conditionalPanel(
                 condition = "input.selcols == true",
                 checkboxGroupInput('sel_cols', 'Columns to Display:',
                                    colnames(data), selected = selected_cols)
               ),
               downloadButton('data_f_dl',"CSV"),
               DT::dataTableOutput('tbl')
             )), 
             fluidRow(
               actionButton(
                 "cs",
                 "Previous Page",
                 style = "color: #fff; background-color: #337ab7; border-color: #2e6da4",
                 icon = icon("chevron-left", lib = "glyphicon"),
                 onclick = "window.close()"
               ),
               actionButton(
                 "en",
                 "GO&KEGG Enrichment",
                 icon = icon("th"),
                 style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
               ),
               p(),
               uiOutput("tb")
             )
             ),
           tabPanel(title = 'Co-expressed Genes (Heatmap)',
                    value = "heatmap_panel",
                    column(
                      12,
                      br(),
                      hidden(div(
                        id = "loading-heatmap",
                        img(src = "../../../images/loading-1.gif")
                      )),
                      #plotOutput("heatmap")#, width = "100%", height = "100%")
                      strong("Limited to most correlated 500 co-expressed genes"),
                      plotlyOutput('heatmap',height = "1500px")
                    )),
           tabPanel(
             title = 'Coordinates of Co-expressed Genes (Circos)',
             value = "circos",
             column(
               12,
               br(),
               numericInput(
                 "top",
                 "# Top ranked co-expressed genes (default):",
                 25,
                 step = 1,
                 max = 100,
                 min = 1
               ),
               hidden(div(
                 id = "loading-circos",
                 img(src = "../../../images/loading-1.gif")
               )),
               imageOutput("circos")
             )
           ),
           # output$network_rbp level1
           # output$circ_rbp_tbl table
           # output$RBP print
           # output$network_circrbpg level2
           tabPanel(title = 'circRNA-RBP (Networks)',
                    value = "crbp_network",
                    br(),
                    strong("Click Apply to start analysis: "),
                    actionButton("apply_value_rbp", "Apply"),hr(),
                    fluidRow(
                      column(
                        12,
                       # br(),
                        hidden(div(
                          id = "loading-crbp",
                          img(src = "../../../images/loading-1.gif")
                        )),
                        strong("circRNA-RBP"),
                        downloadButton('rbp_tbl_dl',"CSV"),
                        DT::dataTableOutput('rbp_tbl'),
                        br(),
                        strong("Co-expressed gene(s)-RBP"),
                        downloadButton('rbp_coexp_tbl_dl',"CSV"),
                        DT::dataTableOutput('rbp_coexp_tbl'),
                        verbatimTextOutput('RBP'),
                        br()
                      )
                    ),
                    fluidRow(column(6,
                                    visNetworkOutput("network_rbp"),
                                    hidden(div(
                                      id = "loading-crbpg",
                                      img(src = "../../../images/loading-1.gif")
                                    ))),
                             column(6,
                                    visNetworkOutput("network_rbpg")))),
         
           # output$network_mir level1 
           # output$sponge_tbl table
           # output$sponge print
           # output$network_mirg level2
           tabPanel(title = 'circRNA-miRNA (Networks)',
                    value = "cmir_network",
                    br(),
                    strong("Click Apply to start analysis: "),
                    actionButton("apply_value_sponge", "Apply"),hr(),
                    fluidRow(
                      column(
                        12,
                        # br(),
                        hidden(div(
                          id = "loading-cmir",
                          img(src = "../../../images/loading-1.gif")
                        )),
                        strong("circRNA-miRNA"),
                        downloadButton('sponge_tbl_dl',"CSV"),
                        DT::dataTableOutput('sponge_tbl'),
                        br(),
                        strong("Co-expressed gene(s)-miRNA"),
                        downloadButton('sponge_coexp_tbl_dl',"CSV"),
                        DT::dataTableOutput('sponge_coexp_tbl'),
                        verbatimTextOutput('sponge'),
                        br()
                      )
                    ),
                    fluidRow(column(6,
                                    visNetworkOutput("network_mir"),
                                    hidden(div(
                                      id = "loading-cmirg",
                                      img(src = "../../../images/loading-1.gif")
                                    ))),
                             column(6,
                                    visNetworkOutput("network_mirg"))))),
                  br()#,
##########################################
## actionButtons
##########################################
#         # actionButton("cs", "Previous Page",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("chevron-left",lib = "glyphicon"),onclick = paste0("window.open('http://120.126.1.61/shiny/circlnc/jobs/",jobid,"/degc/','_self')") ),  
#          actionButton("cs", "Previous Page",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("chevron-left",lib = "glyphicon"),onclick = "window.close()" ),
#          actionButton("en", "GO&KEGG Enrichment",icon = icon("th"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
# 
#       	 actionButton("hm", "Heatmap",icon = icon("th"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
#          actionButton("circos", "Circos",icon = icon("cog"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
#       	 actionButton("crbp", "circRNA-RBP Networks",icon = icon("asterisk"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
#       	 actionButton("cmir", "circRNA-miRNA Networks",icon = icon("asterisk"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
# 	       p(),
#          uiOutput("tb")

##########################################         
      )
   )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output,session) {
      data_f <- reactive({
         if(!input$exclude){
           coexp_tbl = data[(cor < input$corr[2] & cor > input$corr[1]) & circRNA==input$qgene,]
          # coexp_tbl = data[(cor < input$corr[2] & cor > input$corr[1]) & cor_p < input$pvalue & circRNA==input$qgene,] # with pvalue filter
      	 }else{
    	     coexp_tbl = data[(cor > input$corr[2] | cor< input$corr[1])  & circRNA==input$qgene,]
    	    # coexp_tbl = data[(cor > input$corr[2] | cor< input$corr[1]) & cor_p < input$pvalue & circRNA==input$qgene,]
    	   }
           list(coexp_tbl=unique(coexp_tbl))
      }) 

      n_coexp_gene <- reactive({
        if(!input$exclude){
          coexp_tbl = data[(cor < input$corr[2] & cor > input$corr[1]) & circRNA==input$qgene,]
        }else{
          coexp_tbl = data[(cor > input$corr[2] | cor< input$corr[1])  & circRNA==input$qgene,]
        }
        nrow(coexp_tbl)
      }) 
      
      cor_exclude <- reactive({
        message(input$corr[1])
        if(!input$exclude){
          "in"
          #list(range="in")
        }else{
          "ex"
          #list(range="ex")
        }
      }) 
      
   observeEvent(input$hm, {
                updateTabsetPanel(session, "coexp",
                        selected = "heatmap_panel"
                )
        })

   observeEvent(input$crbp, {
                updateTabsetPanel(session, "coexp",
                        selected = "crbp_network"
                )
        })

   observeEvent(input$cmir, {
                updateTabsetPanel(session, "coexp",
                        selected = "cmir_network"
                )
        })

   observeEvent(input$circos, {
                updateTabsetPanel(session, "coexp",
                        selected = "circos"
                )
        })

   output$qsumchart <- renderPlotly({
                m <- list(
                        l = 100,
                        r = 10,
                        b = 100,
                        t = 50,
                        pad = 4
                )
                plot_ly(qsum,y= ~Query, x = ~round(100*`cor > 0.5 (% of coexpressed genes)`), type = 'bar',orientation = 'h', name= 'Correlation coefficient > 0.5') %>% add_trace(x = ~round(100*`pvalue < 0.05 (% of coexpressed genes)`), name = 'P-value < 0.05') %>% layout(legend = list(x = 0, y = 1.2),autosize=F, width=350, height=500,title='' ,yaxis=list(title=''),xaxis=list(title='% of Co-expressed Transcriptome',ticksuffix = "%"),barmode='group',margin=m)

   })

#   output$qsumchart <- renderGvis({
#      gvisBarChart(qsum,options=list(title="% of Co-expressed Transcriptome", hAxis="{format:'#,###%'}",width=300,height=400,legend='bottom'))
#   })
   
   output$boxplot = renderImage({
     list(src = "../output/boxplot.png",
          contentType = 'image/png',
          width = 500,
          height = 500,
          alt = "This is alternate text")
   }, deleteFile = FALSE)
 
   output$co_exp_gene <- renderUI({
          selectInput("coxgene","Coexpressed gene",unique(data_f()$coexp_tbl$co_exp_gene),selectize = TRUE,selected =  unique(data_f()$coexp_tbl$co_exp_gene)[1])
   })

# output$scatterplot2<-renderPlot({
#      message(input$coxgene)
#      if(!is.null(input$coxgene)){
#        rds<-paste0("../output/",input$qgene,"_",input$coxgene,"_scatterplot.rds")
#        if(!file.exists(rds)){
#          message("if")
#          system(paste0("cd ..;",run_scatterplot," -q ",input$qgene," -c ",input$coxgene))
#        }
#        message(" read rds")
#        p<-readRDS(rds)
#        print(p)
#      } else { NULL }
# 
#    })
   
   output$scatterplot <- renderImage({
     if(!file.exists(paste0("../output/",input$qgene,"_",input$coxgene,"_scatterplot.png"))){
       system(paste0("cd ..;",run_scatterplot," -q ",input$qgene," -c ",input$coxgene))
     }
     
     list(src = paste0("../output/",input$qgene,"_",input$coxgene,"_scatterplot.png"),
          contentType = 'image/png',
                width = '100%',
                height = '100%',
          alt = "Generating Scatter Plot...")
   }, deleteFile = FALSE)
   
   output$randomization <- renderPlotly({
     
     rds<-paste0("../output/random_",input$qgene,".rds")
     message(rds)
     if(!file.exists(rds)){
       show("loading-rand")
       disable_act_but()
       system(paste0("cd ..;",run_randomization," -g ",input$qgene))
       enable_act_but()
       hide("loading-rand")
     }else{
       hide("loading-rand")
     }
     p<-readRDS(rds)
     ggplotly(p)
   })
   
   output$number_coexp <- renderPrint({
     cat(n_coexp_gene())
   })     
   
   output$tbl = DT::renderDataTable(
     data_f()$coexp_tbl[, input$sel_cols, with = F],
     extensions = 'Buttons',
     class = 'compact',
     filter = 'top',
     caption = '',
     options = list(
       # dom = 'Blfrtip',
       # buttons = list(list(extend='copy'),
       #                list(extend='csv',filename=paste("coexp_gene",input$corr[2],input$corr[1],sep="_")),
       #                list(extend='excel',filename=paste("coexp_gene",input$corr[2],input$corr[1],sep="_"))
       # ),
       pageLength = 5,
       lengthMenu = c(5, 10, 15, 20),
       autoWidth = TRUE
     ),
     rownames = FALSE
   )
   
   output$data_f_dl <- downloadHandler(
     filename = function() {paste0(paste("circRNA_co-expressed_gene",input$qgene,cor_exclude(),input$corr[2],input$corr[1],sep="_"),".csv") },
     content = function(file) {
       write.csv(data_f()$coexp_tbl, file, row.names = F)
     }
   )
   
#####################################################

##### pathway #####
   pathway_act<- eventReactive(input$apply_kegg, {# pathway_network_link_bp_AC021218.2_ex_0.5-0.5.txt
     # enrichment_res_kegg_AC021218.2_ex_0.5-0.5_.png
     if(!file.exists(paste("../output/enrichment_res_kegg",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"))){
       show("loading-pathway")
       disable_act_but()
       system(paste0("cd ..;",run_gene_enrichment," -p ",input$corr[2]," -n ",input$corr[1]," -r ",cor_exclude()," -g ",input$qgene," -e kegg"))
       # system(paste0("cd ..;",run_gene_enrichment," -p 0.5 -n -0.5 -r ex -g FIRRE -e kegg"))
       enable_act_but()
       hide("loading-pathway")
     }else{
       hide("loading-pathway")
     }
     
     list(src = paste("../output/enrichment_res_kegg",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"),
          contentType = 'image/png',
          width = "100%",
          height = "100%",
          alt = "No significant enriched pathway found! Please change the correlation range.")})
   output$pathway <- renderImage({
     pathway_act()
   }, deleteFile = FALSE)
   
   pathwayVis_act<- eventReactive(input$apply_kegg, {
     if(!file.exists(paste0("../output/pathway_network_link_kegg","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"))){
       show("loading-pathwayVis")
     }else{
       link_df=fread(paste0("../output/pathway_network_link_kegg","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       node_df=fread(paste0("../output/pathway_network_node_kegg","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       hide("loading-pathwayVis")
     }
     visNetwork(node_df, link_df, height = "1000px", width = "100%") %>%
       visOptions(selectedBy = "group",
                  highlightNearest = TRUE,
                  nodesIdSelection = TRUE) %>%
       visPhysics(maxVelocity=1,stabilization = T) %>%
       visLayout(randomSeed = 123) 
     #  %>%
     #  visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
     
   })
   
   output$pathwayVis <- renderVisNetwork({
     pathwayVis_act() 
   })
   
   pathwayhm_act<- eventReactive(input$apply_hm, {
      
     if(!file.exists(paste("../output/enrichment_res_hm",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"))){
       show("loading-pathwayhm")
       disable_act_but()
       system(paste0("cd ..;",run_gene_enrichment," -p ",input$corr[2]," -n ",input$corr[1]," -r ",cor_exclude()," -g ",input$qgene," -e hm"))
       enable_act_but()
       hide("loading-pathwayhm")
     }else{
       hide("loading-pathwayhm")
     }
     
     list(src = paste("../output/enrichment_res_hm",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"),
          contentType = 'image/png',
          width = "100%",
          height = "100%",
          alt = "No significant enriched pathway found! Please change the correlation range.")
   })
   output$pathwayhm <- renderImage({
     pathwayhm_act()
   }, deleteFile = FALSE)
   pathwayhmVis_act<- eventReactive(input$apply_hm, {
     if(!file.exists(paste0("../output/pathway_network_link_hm","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"))){
       show("loading-pathwayhmVis")
     }else{
       link_df=fread(paste0("../output/pathway_network_link_hm","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       node_df=fread(paste0("../output/pathway_network_node_hm","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       hide("loading-pathwayhmVis")
     }
     visNetwork(node_df, link_df, height = "1000px", width = "100%") %>%
       visOptions(selectedBy = "group",
                  highlightNearest = TRUE,
                  nodesIdSelection = TRUE) %>%
       visPhysics(maxVelocity=1,stabilization = T) %>%
       visLayout(randomSeed = 123) 
     #  %>%
     #  visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
   })
   output$pathwayhmVis <- renderVisNetwork({
     pathwayhmVis_act()
   })
   
   #### bp
   bp_act<- eventReactive(input$apply_bp, {
     if(!file.exists(paste("../output/enrichment_res_bp",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"))){
       show("loading-bp")
       disable_act_but()
       system(paste0("cd ..;",run_gene_enrichment," -p ",input$corr[2]," -n ",input$corr[1]," -r ",cor_exclude()," -g ",input$qgene," -e bp"))
       enable_act_but()
       hide("loading-bp")
     }else{
       hide("loading-bp")
     }
     
     list(src = paste("../output/enrichment_res_bp",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"),
          contentType = 'image/png',
          width = "100%",
          height = "100%",
          alt = "No significant enriched biological process found! Please change the correlation range.")
   })
   output$bp <- renderImage({
     bp_act() 
   }, deleteFile = FALSE)
   
   bpVis_act<- eventReactive(input$apply_bp,{
     if(!file.exists(paste0("../output/pathway_network_link_bp","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"))){
       show("loading-bpVis")
     }else{
       link_df=fread(paste0("../output/pathway_network_link_bp","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       node_df=fread(paste0("../output/pathway_network_node_bp","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       hide("loading-bpVis")
     }
     visNetwork(node_df, link_df, height = "1000px", width = "100%") %>%
       visOptions(selectedBy = "group",
                  highlightNearest = TRUE,
                  nodesIdSelection = TRUE) %>%
       visPhysics(maxVelocity=1,stabilization = T) %>%
       visLayout(randomSeed = 123) 
     #  %>%
     #  visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
   })
   output$bpVis<-renderVisNetwork({
     bpVis_act()
   })
   
   #### mf 
   mf_act<- eventReactive(input$apply_mf, {
     if(!file.exists(paste("../output/enrichment_res_mf",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"))){
       show("loading-mf")
       disable_act_but()
       system(paste0("cd ..;",run_gene_enrichment," -p ",input$corr[2]," -n ",input$corr[1]," -r ",cor_exclude()," -g ",input$qgene," -e mf"))
       enable_act_but()
       hide("loading-mf")
     }else{
       hide("loading-mf")
     }
     
     list(src = paste("../output/enrichment_res_mf",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"),
          contentType = 'image/png',
          width = "100%",
          height = "100%",
          alt = "No significant enriched molecular function found! Please change the correlation range.")
   })
   output$mf <- renderImage({
     mf_act()
   }, deleteFile = FALSE)
   
   
   mfVis_act<- eventReactive(input$apply_mf,{
     if(!file.exists(paste0("../output/pathway_network_link_mf","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"))){
       show("loading-mfVis")
     }else{
       link_df=fread(paste0("../output/pathway_network_link_mf","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       node_df=fread(paste0("../output/pathway_network_node_mf","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       hide("loading-mfVis")
     }
     visNetwork(node_df, link_df, height = "1000px", width = "100%") %>%
       visOptions(selectedBy = "group",
                  highlightNearest = TRUE,
                  nodesIdSelection = TRUE) %>%
       visPhysics(maxVelocity=1,stabilization = T) %>%
       visLayout(randomSeed = 123) 
     #  %>%
     #  visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
   })
   output$mfVis<-renderVisNetwork({
     mfVis_act()
   })
   
   
   #### cc
   cc_act<- eventReactive(input$apply_cc, {
     if(!file.exists(paste("../output/enrichment_res_cc",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"))){
       show("loading-cc")
       disable_act_but()
       system(paste0("cd ..;",run_gene_enrichment," -p ",input$corr[2]," -n ",input$corr[1]," -r ",cor_exclude()," -g ",input$qgene," -e cc"))
       enable_act_but()
       hide("loading-cc")
     }else{
       hide("loading-cc")
     }
     
     list(src = paste("../output/enrichment_res_cc",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"),
          contentType = 'image/png',
          width = "100%",
          height = "100%",
          alt = "No significant enriched cellular component found! Please change the correlation range.")
   })
   output$cc <- renderImage({
     cc_act()
   }, deleteFile = FALSE)
   
   ccVis_act<- eventReactive(input$apply_cc,{
     if(!file.exists(paste0("../output/pathway_network_link_cc","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"))){
       show("loading-ccVis")
     }else{
       link_df=fread(paste0("../output/pathway_network_link_cc","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       node_df=fread(paste0("../output/pathway_network_node_cc","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       hide("loading-ccVis")
     }
     visNetwork(node_df, link_df, height = "1000px", width = "100%") %>%
       visOptions(selectedBy = "group",
                  highlightNearest = TRUE,
                  nodesIdSelection = TRUE) %>%
       visPhysics(maxVelocity=1,stabilization = T) %>%
       visLayout(randomSeed = 123) 
     #  %>%
     #  visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
   })
   output$ccVis<-renderVisNetwork({
     ccVis_act()
   })
   
   
   ##### msigdb tf
   tf_msigdb_act<- eventReactive(input$apply_tf_msigdb, {
     if(!file.exists(paste("../output/enrichment_res_tf",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"))){
       show("loading-tf_msigdb")
       disable_act_but()
       system(paste0("cd ..;",run_gene_enrichment," -p ",input$corr[2]," -n ",input$corr[1]," -r ",cor_exclude()," -g ",input$qgene," -e tf"))
       enable_act_but()
       hide("loading-tf_msigdb")
     }else{
       hide("loading-tf_msigdb")
     }
     
     list(src = paste("../output/enrichment_res_tf",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"),
          contentType = 'image/png',
          width = "100%",
          height = "100%",
          alt = "No significant enriched cellular component found! Please change the correlation range.")
   })
   output$tf_msigdb <- renderImage({
     tf_msigdb_act()
   }, deleteFile = FALSE)
   
   tf_msigdbVis_act<- eventReactive(input$apply_tf_msigdb, {
     if(!file.exists(paste0("../output/pathway_network_link_tf","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"))){
       show("loading-tf_msigdbVis")
     }else{
       link_df=fread(paste0("../output/pathway_network_link_tf","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       node_df=fread(paste0("../output/pathway_network_node_tf","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       hide("loading-tf_msigdbVis")
     }
     visNetwork(node_df, link_df, height = "1000px", width = "100%") %>%
       visOptions(selectedBy = "group",
                  highlightNearest = TRUE,
                  nodesIdSelection = TRUE) %>%
       visPhysics(maxVelocity=1,stabilization = T) %>%
       visLayout(randomSeed = 123) 
     #  %>%
     #  visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
   })
   output$tf_msigdbVis <- renderVisNetwork({
     tf_msigdbVis_act()
   })
   
   #### encode tf
   tf_encode_act<- eventReactive(input$apply_tf_encode, {
     if(!file.exists(paste("../output/enrichment_res_encodetf",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"))){
       show("loading-tf_encode")
       disable_act_but()
       system(paste0("cd ..;",run_gene_enrichment," -p ",input$corr[2]," -n ",input$corr[1]," -r ",cor_exclude()," -g ",input$qgene," -e encodetf"))
       enable_act_but()
       hide("loading-tf_encode")
     }else{
       hide("loading-tf_encode")
     }
     
     list(src = paste("../output/enrichment_res_encodetf",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".png",sep="_"),
          contentType = 'image/png',
          width = "100%",
          height = "100%",
          alt = "No significant enriched cellular component found! Please change the correlation range.")
   })
   output$tf_encode <- renderImage({
     tf_encode_act()
   }, deleteFile = FALSE)
   
   tf_encodeVis_act<- eventReactive(input$apply_tf_encode, {
     if(!file.exists(paste0("../output/pathway_network_link_encodetf","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"))){
       show("loading-tf_encodeVis")
     }else{
       link_df=fread(paste0("../output/pathway_network_link_encodetf","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       node_df=fread(paste0("../output/pathway_network_node_encodetf","_",input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],".txt"),header=T,sep="\t")
       hide("loading-tf_encodeVis")
     }
     visNetwork(node_df, link_df, height = "1000px", width = "100%") %>%
       visOptions(selectedBy = "group",
                  highlightNearest = TRUE,
                  nodesIdSelection = TRUE) %>%
       visPhysics(maxVelocity=1,stabilization = T) %>%
       visLayout(randomSeed = 123) 
     #  %>%
     #  visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
   })
   output$tf_encodeVis <- renderVisNetwork({
     tf_encodeVis_act()
   })
   
   

##### heatmap #####

# output$heatmap2 <- renderImage({
#      if(!file.exists(paste("../output/heatmap",input$qgene,input$corr[2],input$corr[1],cor_exclude(),".png",sep="_"))){
#        show("loading-heatmap")
#        disable_act_but()
#        system(paste0("cd ..;",run_heatmap," -p ",input$corr[2]," -n ",input$corr[1]," -r ",cor_exclude()," -g ",input$qgene))
#        enable_act_but()
#        hide("loading-heatmap")
#      }else{
#        hide("loading-heatmap")
#      }
# 
#      list(src = paste("../output/heatmap",input$qgene,input$corr[2],input$corr[1],cor_exclude(),".png",sep="_"),
#           contentType = 'image/png',
#                 width = '60%',
#                 height = '60%',
#           alt = "No heatmap generated! Please change a different correlation range.")
#    }, deleteFile = FALSE)


output$heatmap <- renderPlotly({
        
    rds<-paste("../output/heatmap",input$qgene,input$corr[2],input$corr[1],cor_exclude(),".rds",sep="_")
    message(rds)
     if(!file.exists(rds)){
       show("loading-heatmap")
       disable_act_but()
       system(paste0("cd ..;",run_heatmap," -p ",input$corr[2]," -n ",input$corr[1]," -r ",cor_exclude()," -g ",input$qgene))
       enable_act_but()
       hide("loading-heatmap")
     }else{
       hide("loading-heatmap")
     }
       p<-readRDS(rds)
       ggplotly(p)
   })
   
   
#####

#----------------------------------------------------------------------------
## circRNA-RBP Networks
## triple-network sponge level1
# 
# loading-crbp
# loading-crbpg
# loading-cmir
# loading-cmirg
# disable_act_but()
network_rbp_check<-reactive({
  parameters1<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  out_txt_nodel1<<-paste0("../output/trinetwork_RBP_node_level1_",parameters1,".txt")
  out_txt_linkl1<<-paste0("../output/trinetwork_RBP_link_level1_",parameters1,".txt")
  
  if(!file.exists(out_txt_nodel1)){
    show("loading-crbp")
    disable_act_but()
    system(paste0("cd ..;",run_triple_network_circRNA_RBP_2step," -q ",input$qgene," -r ",cor_exclude()," -d 2"," -p ",input$corr[2]," -n ",input$corr[1]))
    enable_act_but()
    hide("loading-crbp")
  }else{
    hide("loading-crbp")
  }
})
     
network_rbp_act <- eventReactive(input$apply_value_rbp, {
    network_rbp_check()
    testnode_level1=read.table(out_txt_nodel1,header=T,sep="\t")
    testlink_level1=read.table(out_txt_linkl1,header=T,sep="\t")
    visNetwork(testnode_level1, testlink_level1, height = "1000px", width = "100%") %>%
      visOptions(selectedBy = "group",
                 highlightNearest = TRUE,
                 nodesIdSelection = TRUE) %>%
      visPhysics(maxVelocity=1,stabilization = T) %>%
      visLayout(randomSeed = 123) #%>%
    # visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
  })
output$network_rbp <- renderVisNetwork({
  network_rbp_act()
})

rbp_tbl_act <- eventReactive(input$apply_value_rbp, {
  network_rbp_check()
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","RBP","_",parameters)
  lnc_rbp_table=fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
  # lnc_rbp_table=fread("../output/RBP_co_expressed_merged.txt",header=T,sep="\t")
  
  if ("Error" %in% colnames(lnc_rbp_table)) {
    datatable(
      lnc_rbp_table,
      selection = list(
        mode = "single",
        target = "row",
        selected = c(1)
      ),
      extensions = 'Buttons',
      class = 'compact',
      filter = 'top',
      caption = '',
      options = list(
        # dom = 'Blfrtip',
        # buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
        pageLength = 5,
        lengthMenu = c(5, 10, 15, 20),
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  } else {
    datatable(
      lnc_rbp_table[query_symbol %in% input$qgene, ][order(-support_sources_count)],
      selection = list(
        mode = "single",
        target = "row",
        selected = c(1)
      ),
      extensions = 'Buttons',
      class = 'compact',
      filter = 'top',
      caption = '',
      options = list(
        # dom = 'Blfrtip',
        # buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
        pageLength = 5,
        lengthMenu = c(5, 10, 15, 20),
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  }
 })
output$rbp_tbl = DT::renderDataTable({
  rbp_tbl_act ()
}, server=F)

rbp_tbl_dl <- reactive({
  network_rbp_check()
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","RBP","_",parameters)
  fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
})
output$rbp_tbl_dl <- downloadHandler(
  filename = function() {paste0(paste("circRNA_RBP",input$qgene,cor_exclude(),input$corr[2],input$corr[1],sep="_"),".csv") },
  content = function(file) {
    write.csv(rbp_tbl_dl(), file, row.names = F)
  }
)

rbp_coexp_tbl_act<-eventReactive(input$apply_value_rbp, {
  network_rbp_check()
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","RBP","_",parameters)
  lnc_rbp_table=fread(paste0("../",base,"_co_expressed_along.txt"),header=T,sep="\t")
  # lnc_rbp_table=fread("../output/RBP_co_expressed_merged.txt",header=T,sep="\t")
  
  if ("Error" %in% colnames(lnc_rbp_table)) {
    datatable(
      lnc_rbp_table,
      selection = list(
        mode = "single",
        target = "row",
        selected = c(1)
      ),
      extensions = 'Buttons',
      class = 'compact',
      filter = 'top',
      caption = '',
      options = list(
        # dom = 'Blfrtip',
        # buttons = list(list(extend='copy'),
        #                list(extend='csv',filename=paste("lncRNA_RBP_coexp",input$qgene,input$corr[2],input$corr[1],sep="_")),
        #                list(extend='excel',filename=paste("lncRNA_RBP_coexp",input$qgene,input$corr[2],input$corr[1],sep="_"))
        # ),
        pageLength = 5,
        lengthMenu = c(5, 10, 15, 20),
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  } else {
    datatable(
      lnc_rbp_table[order(-support_sources_count)],
      selection = list(
        mode = "single",
        target = "row",
        selected = c(1)
      ),
      extensions = 'Buttons',
      class = 'compact',
      filter = 'top',
      caption = '',
      options = list(
        # dom = 'Blfrtip',
        # buttons = list(list(extend='copy'),
        #                list(extend='csv',filename=paste("lncRNA_RBP_coexp",input$qgene,input$corr[2],input$corr[1],sep="_")),
        #                list(extend='excel',filename=paste("lncRNA_RBP_coexp",input$qgene,input$corr[2],input$corr[1],sep="_"))
        # ),
        pageLength = 5,
        lengthMenu = c(5, 10, 15, 20),
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  }
})
output$rbp_coexp_tbl<-DT::renderDataTable({
  rbp_coexp_tbl_act()
}, server = T)

rbp_coexp_tbl_dl <- reactive({
  network_rbp_check()
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","RBP","_",parameters)
  fread(paste0("../",base,"_co_expressed_along.txt"),header=T,sep="\t")
})
output$rbp_coexp_tbl_dl <- downloadHandler(
  filename = function() {paste0(paste("circRNA_RBP_coexp",input$qgene,cor_exclude(),input$corr[2],input$corr[1],sep="_"),".csv") },
  content = function(file) {
    write.csv(rbp_coexp_tbl_dl(), file, row.names = F)
  }
)

RBP_act<-eventReactive(input$apply_value_rbp,{
  network_rbp_check()
  s=input$rbp_tbl_rows_selected
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","RBP","_",parameters)
  lnc_rbp_table=fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
  
  if ("Error" %in% colnames(lnc_rbp_table)){
    cat("No RBP available")
  } else { 
    #lnc_rbp_table=fread("../output/RBP_co_expressed_merged.txt",header=T,sep="\t")
    lnc_rbp_table_f=lnc_rbp_table[query_symbol%in%input$qgene,][order(-support_sources_count)]
    cat('Network: \n')
    g=lnc_rbp_table_f$RBP[s]
    cat(input$qgene,g,"co-expressed gene(s)", sep = '<>')
  }
})
output$RBP = renderPrint({
  RBP_act()
})


# triple-network level 2
network_rbpg_act <-eventReactive(input$apply_value_rbp,{
  network_rbp_check()
  s=input$rbp_tbl_rows_selected
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","RBP","_",parameters)
  lnc_rbp_table=fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
  # lnc_rbp_table=fread("../output/RBP_co_expressed_merged.txt",header=T,sep="\t")

  if(is.null(s) & !("Error" %in% colnames(lnc_rbp_table))){
   s<-1
   lnc_rbp_table_f=lnc_rbp_table[query_symbol%in%input$qgene,][order(-support_sources_count)]
   cat('Network: \n')
   g=lnc_rbp_table_f$RBP[s]
  }
  
  if ("Error" %in% colnames(lnc_rbp_table)){
    g<-NULL
  } else {
    lnc_rbp_table_f=lnc_rbp_table[query_symbol%in%input$qgene,][order(-support_sources_count)]
    cat('Network: \n')
    g=lnc_rbp_table_f$RBP[s]
  } 
  
  parameters2<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2","_",g)
  #    parameters2<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  out_txt_nodel2=paste0("../output/trinetwork_RBP_node_level2_",parameters2,".txt")
  out_txt_linkl2=paste0("../output/trinetwork_RBP_link_level2_",parameters2,".txt")
  
  if(!file.exists(out_txt_nodel2) ){
    show("loading-crbpg")
    disable_act_but()
    system(paste0("cd ..;",run_triple_network_circRNA_RBP_2step," -q ",input$qgene," -b ",g," -r ",cor_exclude()," -d 2"," -p ",input$corr[2]," -n ",input$corr[1]))
    #       system(paste0("cd ..;",run_triple_network_circRNA_RBP_2step," -q ",input$qgene," -r ex -d 2"," -p ",input$corr[2]," -n ",input$corr[1]))
    enable_act_but()
     hide("loading-crbpg")
  }else{
    hide("loading-crbpg")
  }
  
  testnode_level2=read.table(out_txt_nodel2,header=T,sep="\t")
  testlink_level2=read.table(out_txt_linkl2,header=T,sep="\t")
  visNetwork(testnode_level2, testlink_level2, height = "1000px", width = "100%") %>%
    visOptions(selectedBy = "group",
               highlightNearest = TRUE,
               nodesIdSelection = TRUE) %>%
    visPhysics(maxVelocity=1,stabilization = T) %>%
    visLayout(randomSeed = 123) #%>%
  # visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
})
output$network_rbpg <- renderVisNetwork({
  network_rbpg_act()
})



# output$network_rbp <- renderVisNetwork({
#   parameters1<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
#   out_txt_nodel1=paste0("../output/trinetwork_RBP_node_level1_",parameters1,".txt")
#   out_txt_linkl1=paste0("../output/trinetwork_RBP_link_level1_",parameters1,".txt")
# 
#   if(!file.exists(out_txt_nodel1)){
#     show("loading-lncrbp")
#     system(paste0("cd ..;",run_triple_network_circRNA_RBP_2step," -q ",input$qgene," -r ex -d 2"," -p ",input$corr[2]," -n ",input$corr[1]))
#     hide("loading-lncrbp")
#   }else{
#     hide("loading-lncrbp")
#   }
#   testnode_level1=read.table(out_txt_nodel1,header=T,sep="\t")
#   testlink_level1=read.table(out_txt_linkl1,header=T,sep="\t")
#   visNetwork(testnode_level1, testlink_level1, height = "1000px", width = "100%") %>%
#     visOptions(selectedBy = "group",
#                highlightNearest = TRUE,
#                nodesIdSelection = TRUE) %>%
#     visPhysics(maxVelocity=1,stabilization = T) %>%
#     visLayout(randomSeed = 123) #%>%
#   # visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
#   })

# output$circ_rbp_tbl = DT::renderDataTable({
#   
#   parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
#   base<-paste0("output/trinetwork_tab_","RBP","_",parameters)
#   lnc_rbp_table=fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
#   # lnc_rbp_table=fread("../output/RBP_co_expressed_merged.txt",header=T,sep="\t")
#   
#   if ("Error" %in% colnames(lnc_rbp_table)){
#     datatable(lnc_rbp_table,selection=list(mode="single",target="row",selected=c(1)),extensions='Buttons',class = 'compact',filter='top',caption = '',options = list(dom='Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),pageLength = 5,lengthMenu = c(5, 10, 15, 20), autoWidth = TRUE), rownames= FALSE)
#   } else{
#     datatable(lnc_rbp_table[query_symbol%in%input$qgene,][order(-support_sources_count)],selection=list(mode="single",target="row",selected=c(1)),extensions='Buttons',class = 'compact',filter='top',caption = '',options = list(dom='Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),pageLength = 5,lengthMenu = c(5, 10, 15, 20), autoWidth = TRUE), rownames= FALSE)
#   }
#   
# })

# output$RBP = renderPrint({
#   s=input$circ_rbp_tbl_rows_selected
#   
#   parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
#   base<-paste0("output/trinetwork_tab_","RBP","_",parameters)
#   lnc_rbp_table=fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
#   
#   if ("Error" %in% colnames(lnc_rbp_table)){
#     cat("No RBP available")
#   } else { 
#     #lnc_rbp_table=fread("../output/RBP_co_expressed_merged.txt",header=T,sep="\t")
#     lnc_rbp_table_f=lnc_rbp_table[query_symbol%in%input$qgene,][order(-support_sources_count)]
#     cat('Network: \n')
#     g=lnc_rbp_table_f$RBP[s]
#     cat(input$qgene,g,"co-expressed gene(s)", sep = '<>')
#   }
# })

# triple-network level 2
# output$network_circrbpg <- renderVisNetwork({
#   s=input$circ_rbp_tbl_rows_selected
#   parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
#   base<-paste0("output/trinetwork_tab_","RBP","_",parameters)
#   lnc_rbp_table=fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
#   # lnc_rbp_table=fread("../output/RBP_co_expressed_merged.txt",header=T,sep="\t")
#   if ("Error" %in% colnames(lnc_rbp_table)){
#     g<-NULL
#   } else {
#     lnc_rbp_table_f=lnc_rbp_table[query_symbol%in%input$qgene,][order(-support_sources_count)]
#     cat('Network: \n')
#     g=lnc_rbp_table_f$RBP[s]
#   } 
#   
#   parameters2<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2","_",g)
#   #    parameters2<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
#   out_txt_nodel2=paste0("../output/trinetwork_RBP_node_level2_",parameters2,".txt")
#   out_txt_linkl2=paste0("../output/trinetwork_RBP_link_level2_",parameters2,".txt")
#   
#   if(!file.exists(out_txt_nodel2)){
#     show("loading-lncrbpg")
#     system(paste0("cd ..;",run_triple_network_circRNA_RBP_2step," -q ",input$qgene," -b ",g," -r ex -d 2"," -p ",input$corr[2]," -n ",input$corr[1]))
#     #       system(paste0("cd ..;",run_triple_network_circRNA_RBP_2step," -q ",input$qgene," -r ex -d 2"," -p ",input$corr[2]," -n ",input$corr[1]))
#     hide("loading-lncrbpg")
#   }else{
#     hide("loading-lncrbpg")
#   }
#   testnode_level2=read.table(out_txt_nodel2,header=T,sep="\t")
#   testlink_level2=read.table(out_txt_linkl2,header=T,sep="\t")
#   visNetwork(testnode_level2, testlink_level2, height = "1000px", width = "100%") %>%
#     visOptions(selectedBy = "group",
#                highlightNearest = TRUE,
#                nodesIdSelection = TRUE) %>%
#     visPhysics(maxVelocity=1,stabilization = T) %>%
#     visLayout(randomSeed = 123) #%>%
#   # visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
# })
# 

#--------------------------------------------------------------------------------------------
# triple-network level 2

#----------------------------------------------------------------------------
## circRNA-miRNA Networks
## triple-network sponge level1
network_mir_check<-reactive({
  parameters1<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  out_txt_nodel1<<-paste0("../output/trinetwork_sponge_node_level1_",parameters1,".txt")
  out_txt_linkl1<<-paste0("../output/trinetwork_sponge_link_level1_",parameters1,".txt")
  
  if(!file.exists(out_txt_nodel1)){
    show("loading-cmir")
    disable_act_but()
    system(paste0("cd ..;",run_triple_network_circRNA_sponge_2step," -q ",input$qgene," -r ",cor_exclude()," -d 2"," -p ",input$corr[2]," -n ",input$corr[1]))
    enable_act_but()
    hide("loading-cmir")
  }else{
    hide("loading-cmir")
  }
})

network_mir_act <- eventReactive(input$apply_value_sponge, {
  network_mir_check()
  testnode_level1=read.table(out_txt_nodel1,header=T,sep="\t")
  testlink_level1=read.table(out_txt_linkl1,header=T,sep="\t")
  visNetwork(testnode_level1, testlink_level1, height = "1000px", width = "100%") %>%
    visOptions(selectedBy = "group",
               highlightNearest = TRUE,
               nodesIdSelection = TRUE) %>%
    visPhysics(maxVelocity=1,stabilization = T) %>%
    visLayout(randomSeed = 123) #%>%
  #visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
})
output$network_mir <- renderVisNetwork({
  network_mir_act()
})

sponge_tbl_act<-eventReactive(input$apply_value_sponge, {
  network_mir_check()
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","sponge","_",parameters)
  lnc_sponge_table=fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
  # lnc_sponge_table=fread("../output/sponge_co_expressed_merged.txt",header=T,sep="\t")
  if ("Error" %in% colnames(lnc_sponge_table)){
    datatable(
      lnc_sponge_table,
      selection = list(
        mode = "single",
        target = "row",
        selected = c(1)
      ),
      extensions = 'Buttons',
      class = 'compact',
      filter = 'top',
      caption = '',
      options = list(
        # dom = 'Blfrtip',
        # buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
        pageLength = 5,
        lengthMenu = c(5, 10, 15, 20),
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  } else{
    datatable(
      lnc_sponge_table[query_symbol %in% input$qgene, ][order(-support_sources_count)],
      selection = list(
        mode = "single",
        target = "row",
        selected = c(1)
      ),
      extensions = 'Buttons',
      class = 'compact',
      filter = 'top',
      caption = '',
      options = list(
        # dom = 'Blfrtip',
        # buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
        pageLength = 5,
        lengthMenu = c(5, 10, 15, 20),
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  }
})
output$sponge_tbl = DT::renderDataTable({
  sponge_tbl_act()
}, server=T)

sponge_tbl_dl<-reactive({
  network_mir_check()
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","sponge","_",parameters)
  fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
})
output$sponge_tbl_dl <- downloadHandler(
  filename = function() {paste0(paste("circRNA_miRNA",input$qgene,cor_exclude(),input$corr[2],input$corr[1],sep="_"),".csv") },
  content = function(file) {
    write.csv(sponge_tbl_dl(), file, row.names = F)
  }
)

sponge_coexp_tbl_act<-eventReactive(input$apply_value_sponge, {
  network_mir_check()
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","sponge","_",parameters)
  lnc_sponge_table=fread(paste0("../",base,"_co_expressed_along.txt"),header=T,sep="\t")
  # lnc_sponge_table=fread("../output/sponge_co_expressed_merged.txt",header=T,sep="\t")
  if ("Error" %in% colnames(lnc_sponge_table)) {
    datatable(
      lnc_sponge_table,
      selection = list(
        mode = "single",
        target = "row",
        selected = c(1)
      ),
      extensions = 'Buttons',
      class = 'compact',
      filter = 'top',
      caption = '',
      options = list(
        # dom = 'Blfrtip',
        # buttons = list(list(extend='copy'),
        #                list(extend='csv',filename=paste("lncRNA_RBP_coexp",input$qgene,input$corr[2],input$corr[1],sep="_")),
        #                list(extend='excel',filename=paste("lncRNA_RBP_coexp",input$qgene,input$corr[2],input$corr[1],sep="_"))
        # ),
        pageLength = 5,
        lengthMenu = c(5, 10, 15, 20),
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  } else {
    datatable(
      lnc_sponge_table[order(-support_sources_count)],
      selection = list(
        mode = "single",
        target = "row",
        selected = c(1)
      ),
      extensions = 'Buttons',
      class = 'compact',
      filter = 'top',
      caption = '',
      options = list(
        # dom = 'Blfrtip',
        # buttons = list(list(extend='copy'),
        #                list(extend='csv',filename=paste("lncRNA_RBP_coexp",input$qgene,input$corr[2],input$corr[1],sep="_")),
        #                list(extend='excel',filename=paste("lncRNA_RBP_coexp",input$qgene,input$corr[2],input$corr[1],sep="_"))
        # ),
        pageLength = 5,
        lengthMenu = c(5, 10, 15, 20),
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  }
})
output$sponge_coexp_tbl = DT::renderDataTable({
  sponge_coexp_tbl_act()
}, server=T)

sponge_coexp_tbl_dl<-reactive({
  network_mir_check()
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","sponge","_",parameters)
  fread(paste0("../",base,"_co_expressed_along.txt"),header=T,sep="\t")
})
output$sponge_coexp_tbl_dl <- downloadHandler(
  filename = function() {paste0(paste("circRNA_miRNA_coexp",input$qgene,cor_exclude(),input$corr[2],input$corr[1],sep="_"),".csv") },
  content = function(file) {
    write.csv(sponge_coexp_tbl_dl(), file, row.names = F)
  }
)

sponge_act<-eventReactive(input$apply_value_sponge, {
  network_mir_check()
  s=input$sponge_tbl_rows_selected
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","sponge","_",parameters)
  lnc_sponge_table=fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
  # lnc_sponge_table=fread("../output/sponge_co_expressed_merged.txt",header=T,sep="\t")
  
  if ("Error" %in% colnames(lnc_sponge_table)){
    cat("No miRNA targeting site available")
  } else {
    lnc_sponge_table_f=lnc_sponge_table[query_symbol%in%input$qgene,][order(-support_sources_count)]
    cat('Network: \n')
    g=lnc_sponge_table_f$miRNA[s]
    cat(input$qgene,g,"co-expressed gene(s)", sep = '<>')
  }
})
output$sponge = renderPrint({
  sponge_act()
})

# triple-network level 2
network_mirg_act<-eventReactive(input$apply_value_sponge,{
  network_mir_check()
  s=input$sponge_tbl_rows_selected
  parameters<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2")
  base<-paste0("output/trinetwork_tab_","sponge","_",parameters)
  lnc_sponge_table=fread(paste0("../",base,"_query_along.txt"),header=T,sep="\t")
  # lnc_sponge_table=fread("../output/sponge_co_expressed_merged.txt",header=T,sep="\t")
  
  if(is.null(s) & !("Error" %in% colnames(lnc_sponge_table))){
    s<-1
    lnc_sponge_table_f=lnc_sponge_table[query_symbol%in%input$qgene,][order(-support_sources_count)]
    cat('Network: \n')
    g=lnc_sponge_table_f$miRNA[s]
  }
  
  if ("Error" %in% colnames(lnc_sponge_table)){
    g<-NULL
  } else {
    lnc_sponge_table_f=lnc_sponge_table[query_symbol%in%input$qgene,][order(-support_sources_count)]
    cat('Network: \n')
    g=lnc_sponge_table_f$miRNA[s]
  }
  
  parameters2<-paste0(input$qgene,"_",cor_exclude(),"_",input$corr[2],input$corr[1],"_","2","_",g)
  out_txt_nodel2=paste0("../output/trinetwork_sponge_node_level2_",parameters2,".txt")
  out_txt_linkl2=paste0("../output/trinetwork_sponge_link_level2_",parameters2,".txt")
  
  if(!file.exists(out_txt_nodel2)){
    show("loading-cmirg")
    disable_act_but()
    system(paste0("cd ..;",run_triple_network_circRNA_sponge_2step," -q ",input$qgene," -b ",g," -r ",cor_exclude()," -d 2"," -p ",input$corr[2]," -n ",input$corr[1]))
    enable_act_but()
    hide("loading-cmirg")
  }else{
    hide("loading-cmirg")
  }
  testnode_level2=read.table(out_txt_nodel2,header=T,sep="\t")
  testlink_level2=read.table(out_txt_linkl2,header=T,sep="\t")
  visNetwork(testnode_level2, testlink_level2, height = "1000px", width = "100%") %>%
    visOptions(selectedBy = "group",
               highlightNearest = TRUE,
               nodesIdSelection = TRUE) %>%
    visPhysics(maxVelocity=1,stabilization = T) %>%
    visLayout(randomSeed = 123) #%>%
  #visEdges(smooth =  list(enabled = TRUE, type = 'continuous',roundness=0),physics=F)
})
output$network_mirg <- renderVisNetwork({
  network_mirg_act()
})
disable_act_but<-function()({
  message("disabling but")
  shinyjs::disable("apply_value_rbp")
  disable("apply_value_sponge")
  shinyjs::disable("apply_kegg")
  shinyjs::disable("apply_hm")
  shinyjs::disable("apply_bp")
  shinyjs::disable("apply_mf")
  shinyjs::disable("apply_cc")
  shinyjs::disable("data_f_dl")
  shinyjs::disable("rbp_tbl_dl")
  shinyjs::disable("rbp_coexp_tbl_dl")
  shinyjs::disable("sponge_tbl_dl")
  shinyjs::disable("sponge_coexp_tbl_dl")
  shinyjs::disable("cs")
  shinyjs::disable("en")
  shinyjs::disable("corr")
  shinyjs::disable("qgene")
  shinyjs::disable("exclude")
  shinyjs::disable("pvalue")
  shinyjs::disable("coxgene")
  shinyjs::disable("top")
  shinyjs::disable("apply_tf_msigdb")
  shinyjs::disable("apply_tf_encode")

  
})

enable_act_but<-function()({
  message("enabling but")
  shinyjs::enable("apply_value_rbp")
  shinyjs::enable("apply_value_sponge")
  shinyjs::enable("apply_kegg")
  shinyjs::enable("apply_hm")
  shinyjs::enable("apply_bp")
  shinyjs::enable("apply_mf")
  shinyjs::enable("apply_cc")
  shinyjs::enable("data_f_dl")
  shinyjs::enable("rbp_tbl_dl")
  shinyjs::enable("rbp_coexp_tbl_dl")
  shinyjs::enable("sponge_tbl_dl")
  shinyjs::enable("sponge_coexp_tbl_dl")
  shinyjs::enable("cs")
  shinyjs::enable("en")
  shinyjs::enable("corr")
  shinyjs::enable("qgene")
  shinyjs::enable("exclude")
  shinyjs::enable("pvalue")
  shinyjs::enable("coxgene")
  shinyjs::enable("top")
  shinyjs::enable("apply_tf_msigdb")
  shinyjs::enable("apply_tf_encode")
  
})


# circos
output$circos <- renderImage({
  if(!file.exists(paste0("../output/",input$qgene,"_",input$top,".png"))){
    show("loading-circos")
    disable_act_but()
    system(paste0("cd ..;",run_circos," -f output -q ",input$qgene," -t ",input$top))
    enable_act_but()
    hide("loading-circos")
  }else{
    hide("loading-circos")
  }
  
  list(src = paste0("../output/",input$qgene,"_",input$top,".png"),
       contentType = 'image/png',
       width = 800,
       height = 800,
       alt = "Generating Circos plot for Coordinates of Co-expressed Genes...")
}, deleteFile = FALSE)

#####################################################  12.16
# enrichment 
## table
# 
# output$kegg_tbl = DT::renderDataTable({
# #enrichment_res_kegg_AC021218.2_0.5-0.5_ex_.txt
#        # paste("../output/enrichment_res_kegg",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_")
#         kegg_tbl=fread(paste("../output/enrichment_res_kegg",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
# 
#         datatable(kegg_tbl,selection=list(mode="single",target="row",selected=c(1)),extensions='Buttons',class = 'compact',filter='top',caption = '',options = list(dom='Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),pageLength = 5,lengthMenu = c(5, 10, 15, 20), autoWidth = TRUE), rownames= FALSE)
# })
# 
# output$kegghm_tbl = DT::renderDataTable({
#         kegghm_tbl=fread(paste("../output/enrichment_res_hm",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
#         datatable(kegghm_tbl,selection=list(mode="single",target="row",selected=c(1)),extensions='Buttons',class = 'compact',filter='top',caption = '',options = list(dom='Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),pageLength = 5,lengthMenu = c(5, 10, 15, 20), autoWidth = TRUE), rownames= FALSE)
#    })
# 
# output$gobp_tbl = DT::renderDataTable({
#         gobp_tbl=fread(paste("../output/enrichment_res_bp",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
#         datatable(gobp_tbl,selection=list(mode="single",target="row",selected=c(1)),extensions='Buttons',class = 'compact',filter='top',caption = '',options = list(dom='Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),pageLength = 5,lengthMenu = c(5, 10, 15, 20), autoWidth = TRUE), rownames= FALSE)
#    })
# 
# output$gomf_tbl = DT::renderDataTable({
#         gomf_tbl=fread(paste("../output/enrichment_res_mf",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
#         datatable(gomf_tbl,selection=list(mode="single",target="row",selected=c(1)),extensions='Buttons',class = 'compact',filter='top',caption = '',options = list(dom='Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),pageLength = 5,lengthMenu = c(5, 10, 15, 20), autoWidth = TRUE), rownames= FALSE)
#    })
# 
# output$gocc_tbl = DT::renderDataTable({
#         gocc_tbl=fread(paste("../output/enrichment_res_cc",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
#         datatable(gocc_tbl,selection=list(mode="single",target="row",selected=c(1)),extensions='Buttons',class = 'compact',filter='top',caption = '',options = list(dom='Blfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),pageLength = 5,lengthMenu = c(5, 10, 15, 20), autoWidth = TRUE), rownames= FALSE)
#    })
# 
# 
# 
# 
# output$tb <- renderUI({
#   if (input$en == 0)
#     return()
#   else
#     tabsetPanel(
#       tabPanel(
#         strong("Enriched Pathway (KEGG)"),
#         fluidRow(
#           column(
#             6,
#             br(),
#             strong("KEGG PATHWAY"),
#             div(id = "loading-pathway",
#                 img(src = "../../../images/loading-1.gif")),
#             imageOutput("pathway")
#           ),
#           column(
#             6,
#             br(),
#             strong("NETWORK"),
#             div(id = "loading-pathwayVis",
#                 img(src = "../../../images/loading-1.gif")),
#             visNetworkOutput("pathwayVis")
#           )
#         ),
#         fluidRow(column(
#           12, br(), strong("TABLE"),
#           DT::dataTableOutput('kegg_tbl')
#         ))
#       ),
#       tabPanel(
#         strong("Enriched Pathway (MSigDB)"),
#         fluidRow(
#           column(
#             6,
#             strong("HALLMARK PATHWAY (MSigDB)"),
#             div(id = "loading-pathwayhm",
#                 img(src = "../../../images/loading-1.gif")),
#             imageOutput("pathwayhm")
#           ),
#           column(
#             6,
#             strong("NETWORK"),
#             div(id = "loading-pathwayhmVis",
#                 img(src = "../../../images/loading-1.gif")),
#             visNetworkOutput("pathwayhmVis")
#           )
#         ),
#         fluidRow(column(
#           12, strong("TABLE"),
#           DT::dataTableOutput('kegghm_tbl')
#         ))
#       ),
#       tabPanel(
#         strong("Enriched Biological Process (GO)"),
#         div(id = "loading-bp",
#             img(src = "../../../images/loading-1.gif")),
#         imageOutput("bp"),
#         fluidRow(column(
#           12, strong("TABLE"),
#           DT::dataTableOutput('gobp_tbl')
#         ))
#       ),
#       tabPanel(
#         strong("Enriched Molecular Function (GO)"),
#         div(id = "loading-mf",
#             img(src = "../../../images/loading-1.gif")),
#         imageOutput("mf"),
#         fluidRow(column(
#           12, strong("TABLE"),
#           DT::dataTableOutput('gomf_tbl')
#         ))
#       ),
#       tabPanel(
#         strong("Enriched Cellular Component (GO)"),
#         div(id = "loading-cc",
#             img(src = "../../../images/loading-1.gif")),
#         imageOutput("cc"),
#         fluidRow(column(
#           12, strong("TABLE"),
#           DT::dataTableOutput('gocc_tbl')
#         ))
#       )
#     )
# })

kegg_tbl_act<-eventReactive(input$apply_kegg,{
  #enrichment_res_kegg_AC021218.2_0.5-0.5_ex_.txt
  # paste("../output/enrichment_res_kegg",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_")
  kegg_tbl=fread(paste("../output/enrichment_res_kegg",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
  datatable(
    kegg_tbl,
    selection = list(
      mode = "single",
      target = "row",
      selected = c(1)
    ),
    extensions = 'Buttons',
    class = 'compact',
    filter = 'top',
    caption = '',
    options = list(
      lengthMenu=c(5,10,25,50,100),
      dom = 'lBfrtip',
      buttons = list(#list(extend='copy'),
        list(extend='csv',filename=paste("kegg",input$qgene,input$corr[2],input$corr[1],sep="_"))#,
        #list(extend='excel',filename=paste("kegg",input$qgene,input$corr[2],input$corr[1],sep="_"))
      ),
      pageLength = 5,
      autoWidth = TRUE, 
      scrollX=TRUE
    ),
    rownames = FALSE
  )
})
output$kegg_tbl = DT::renderDataTable({
  kegg_tbl_act()}, server = F)

kegghm_tbl_act<-eventReactive(input$apply_hm,{
  kegghm_tbl=fread(paste("../output/enrichment_res_hm",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
  datatable(
    kegghm_tbl,
    selection = list(
      mode = "single",
      target = "row",
      selected = c(1)
    ),
    extensions = 'Buttons',
    class = 'compact',
    filter = 'top',
    caption = '',
    options = list(
      dom = 'lBfrtip',
      buttons = list(#list(extend='copy'),
        list(extend='csv',filename=paste("hm",input$qgene,input$corr[2],input$corr[1],sep="_"))#,
        #list(extend='excel',filename=paste("hm",input$qgene,input$corr[2],input$corr[1],sep="_"))
      ),
      pageLength = 5,
      lengthMenu = c(5, 10, 15, 20),
      autoWidth = TRUE,scrollX=TRUE
    ),
    rownames = FALSE
  )
})
output$kegghm_tbl = DT::renderDataTable({
  kegghm_tbl_act()
}, server = F)

gobp_tbl_act<-eventReactive(input$apply_bp,{
  gobp_tbl=fread(paste("../output/enrichment_res_bp",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
  datatable(
    gobp_tbl,
    selection = list(
      mode = "single",
      target = "row",
      selected = c(1)
    ),
    extensions = 'Buttons',
    class = 'compact',
    filter = 'top',
    caption = '',
    options = list(
      dom = 'lBfrtip',
      buttons = list(#list(extend='copy'),
        list(extend='csv',filename=paste("go_bp",input$qgene,input$corr[2],input$corr[1],sep="_"))#,
        #list(extend='excel',filename=paste("go_bp",input$qgene,input$corr[2],input$corr[1],sep="_"))
      ),
      pageLength = 5,
      lengthMenu = c(5, 10, 15, 20),
      autoWidth = TRUE, scrollX=TRUE
    ),
    rownames = FALSE
  )
})
output$gobp_tbl = DT::renderDataTable({
  gobp_tbl_act()
}, server = F)

gomf_tbl_act<-eventReactive(input$apply_mf,{
  gomf_tbl=fread(paste("../output/enrichment_res_mf",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
  datatable(
    gomf_tbl,
    selection = list(
      mode = "single",
      target = "row",
      selected = c(1)
    ),
    extensions = 'Buttons',
    class = 'compact',
    filter = 'top',
    caption = '',
    options = list(
      dom = 'lBfrtip',
      buttons = list(#list(extend='copy'),
        list(extend='csv',filename=paste("go_mf",input$qgene,input$corr[2],input$corr[1],sep="_"))#,
        #list(extend='excel',filename=paste("go_mf",input$qgene,input$corr[2],input$corr[1],sep="_"))
      ),
      pageLength = 5,
      lengthMenu = c(5, 10, 15, 20),
      autoWidth = TRUE, scrollX=TRUE
    ),
    rownames = FALSE
  )
})
output$gomf_tbl = DT::renderDataTable({
  gomf_tbl_act()
}, server = F)

gocc_tbl_act<-eventReactive(input$apply_cc,{
  gocc_tbl=fread(paste("../output/enrichment_res_cc",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
  datatable(
    gocc_tbl,
    selection = list(
      mode = "single",
      target = "row",
      selected = c(1)
    ),
    extensions = 'Buttons',
    class = 'compact',
    filter = 'top',
    caption = '',
    options = list(
      dom = 'lBfrtip',
      buttons = list(#list(extend='copy'),
        list(extend='csv',filename=paste("go_cc",input$qgene,input$corr[2],input$corr[1],sep="_"))
        # list(extend='excel',filename=paste("go_cc",input$qgene,input$corr[2],input$corr[1],sep="_"))
      ),
      pageLength = 5,
      lengthMenu = c(5, 10, 15, 20),
      autoWidth = TRUE, scrollX=TRUE
    ),
    rownames = FALSE
  )
})
output$gocc_tbl = DT::renderDataTable({
  gocc_tbl_act()
}, server = F)

tf_msigdb_tbl_act<-eventReactive(input$apply_tf_msigdb,{
  tbl=fread(paste("../output/enrichment_res_tf",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
  datatable(
    tbl,
    selection = list(
      mode = "single",
      target = "row",
      selected = c(1)
    ),
    extensions = 'Buttons',
    class = 'compact',
    filter = 'top',
    caption = '',
    options = list(
      dom = 'lBfrtip',
      buttons = list(#list(extend='copy'),
        list(extend='csv',filename=paste("go_tf",input$qgene,input$corr[2],input$corr[1],sep="_"))
        # list(extend='excel',filename=paste("go_cc",input$qgene,input$corr[2],input$corr[1],sep="_"))
      ),
      pageLength = 5,
      lengthMenu = c(5, 10, 15, 20),
      autoWidth = TRUE, scrollX=TRUE
    ),
    rownames = FALSE
  )
})
output$tf_msigdb_tbl = DT::renderDataTable({
  tf_msigdb_tbl_act()
}, server = F)


tf_encode_tbl_act<-eventReactive(input$apply_tf_encode,{
  tbl=fread(paste("../output/enrichment_res_encodetf",input$qgene,cor_exclude(),paste0(input$corr[2],input$corr[1]),".txt",sep="_"),header=T,sep="\t")
  datatable(
    tbl,
    selection = list(
      mode = "single",
      target = "row",
      selected = c(1)
    ),
    extensions = 'Buttons',
    class = 'compact',
    filter = 'top',
    caption = '',
    options = list(
      dom = 'lBfrtip',
      buttons = list(#list(extend='copy'),
        list(extend='csv',filename=paste("go_encodetf",input$qgene,input$corr[2],input$corr[1],sep="_"))
        # list(extend='excel',filename=paste("go_cc",input$qgene,input$corr[2],input$corr[1],sep="_"))
      ),
      pageLength = 5,
      lengthMenu = c(5, 10, 15, 20),
      autoWidth = TRUE, scrollX=TRUE
    ),
    rownames = FALSE
  )
})
output$tf_encode_tbl = DT::renderDataTable({
  tf_encode_tbl_act()
}, server = F)

########################################################
output$tb <- renderUI({
  if (input$en == 0)
    return()
  else
    tabsetPanel(
      tabPanel(
        strong("Enriched Pathway (KEGG)"),
        strong("Click Apply to start analysis: "),
        actionButton("apply_kegg", "Apply"),
        fluidRow(
          column(
            6,
            strong("KEGG PATHWAY"),
            hidden(div(id = "loading-pathway",
                       img(src = "../../../images/loading-1.gif"))),
            imageOutput("pathway")#,
          ),
          column(
            6,
            strong("NETWORK"),
            hidden( div(id = "loading-pathwayVis",
                        img(src = "../../../images/loading-1.gif"))),
            visNetworkOutput("pathwayVis")
          )
        ),
        #fluidRow
        fluidRow(column(
          12, strong("TABLE") ,
          DT::dataTableOutput('kegg_tbl')
        ))
      ),
      tabPanel(
        strong("Enriched Pathway (MSigDB)"),
        strong("Click Apply to start analysis: "),
        actionButton("apply_hm", "Apply"),
        fluidRow(
          column(
            6,
            strong("HALLMARK PATHWAY (MSigDB)"),
            hidden(div(id = "loading-pathwayhm",
                       img(src = "../../../images/loading-1.gif"))),
            imageOutput("pathwayhm")
          ),
          column(
            6,
            strong("NETWORK"),
            hidden(div(id = "loading-pathwayhmVis",
                       img(src = "../../../images/loading-1.gif"))),
            visNetworkOutput("pathwayhmVis")
          )
        ),
        fluidRow(column(
          12, strong("TABLE"),
          DT::dataTableOutput('kegghm_tbl')
        ))
      ),
      tabPanel(
        strong("Enriched Biological Process (GO)"),
        strong("Click Apply to start analysis: "),
        actionButton("apply_bp", "Apply"),
        fluidRow(column(
          6,
          strong("BIOLOGICAL PROCESS"),
          hidden(div(id = "loading-bp",
                     img(src = "../../../images/loading-1.gif"))),
          imageOutput("bp")
        ),
        column(
          6,
          strong("NETWORK"),
          hidden( div(id = "loading-bpVis",
                      img(src = "../../../images/loading-1.gif"))),
          visNetworkOutput("bpVis")
        )),
        fluidRow(column(
          12, strong("TABLE"),
          DT::dataTableOutput('gobp_tbl')
        ))
      ),
      tabPanel(
        strong("Enriched Molecular Function (GO)"),
        strong("Click Apply to start analysis: "),
        actionButton("apply_mf", "Apply"),
        fluidRow(column(
          6,
          strong("MOLECULAR FUNCTION"),
          hidden(div(id = "loading-mf",
                     img(src = "../../../images/loading-1.gif"))),
          imageOutput("mf")
        ),
        column(
          6,
          strong("NETWORK"),
          hidden( div(id = "loading-mfVis",
                      img(src = "../../../images/loading-1.gif"))),
          visNetworkOutput("mfVis")
        )),
        fluidRow(column(
          12, strong("TABLE"),
          DT::dataTableOutput('gomf_tbl')
        ))
      ),
      tabPanel(
        strong("Enriched Cellular Component (GO)"),
        strong("Click Apply to start analysis: "),
        actionButton("apply_cc", "Apply"),
        fluidRow(column(
          6,
          strong("CELLULAR COMPONENT"),
          hidden(div(id = "loading-cc",
                     img(src = "../../../images/loading-1.gif"))),
          imageOutput("cc")
        ),
        column(
          6,
          strong("NETWORK"),
          hidden( div(id = "loading-ccVis",
                      img(src = "../../../images/loading-1.gif"))),
          visNetworkOutput("ccVis")
        )),
        fluidRow(column(
          12, strong("TABLE"),
          DT::dataTableOutput('gocc_tbl')
        ))
        
      ),
      #### MSigDB 
      tabPanel(
        strong("Transcription factor (MSigdb)"),
        strong("Click Apply to start analysis: "),
        actionButton("apply_tf_msigdb", "Apply"),
        fluidRow(
          column(
            6,
            strong("Transcription factor (MSigdb)"),
            hidden(div(id = "loading-tf_msigdb",
                       img(src = "../../../images/loading-1.gif"))),
            imageOutput("tf_msigdb")#,
          ),
          column(
            6,
            strong("NETWORK"),
            hidden( div(id = "loading-tf_msigdbVis",
                        img(src = "../../../images/loading-1.gif"))),
            visNetworkOutput("tf_msigdbVis")
          )
        ),
        #fluidRow
        fluidRow(column(
          12, strong("TABLE") ,
          DT::dataTableOutput('tf_msigdb_tbl')
        ))
      ),      
      ####ENCODE TF
      tabPanel(
        strong("Transcription factor (ENCODE)"),
        strong("Click Apply to start analysis: "),
        actionButton("apply_tf_encode", "Apply"),
        fluidRow(
          column(
            6,
            strong("Transcription factor (ENCODE)"),
            hidden(div(id = "loading-tf_encode",
                       img(src = "../../../images/loading-1.gif"))),
            imageOutput("tf_encode")#,
          ),
          column(
            6,
            strong("NETWORK"),
            hidden( div(id = "loading-tf_encodeVis",
                        img(src = "../../../images/loading-1.gif"))),
            visNetworkOutput("tf_encodeVis")
          )
        ),
        #fluidRow
        fluidRow(column(
          12, strong("TABLE") ,
          DT::dataTableOutput('tf_encode_tbl')
        ))
      )
    )
})


})
   
# Run the application 
shinyApp(ui = ui, server = server)

