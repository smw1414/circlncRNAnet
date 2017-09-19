#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(ggplot2)
library(plotly)
library(data.table)
library(DT)
#library(plotly)
data=fread("../output/DEGs_lncRNA_ann.txt",header=T,sep="\t",na.strings = c(".",NA))
data=data[,.(gene,log2FoldChange,padj,Gene_full_name,Gene_old_names,Gene_other_names)][order(-abs(log2FoldChange))]

#rbp=fread("../output/RBP_merged.txt",header=T,sep="\t")

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
   useShinyjs(), 
   # Application title
#   titlePanel(h4("")),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        sliderInput("lfc", "Log2 Fold Change:",
                    min = round(min(data$log2FoldChange)), max = round(max(data$log2FoldChange)), value = c(-1, 1)),
        checkboxInput('exclude', 'Exclude range', TRUE),hr(),
        numericInput("padj", "p-adjusted:", 0.05,step = 0.01),hr(),
        tabsetPanel(
          id = 'pca',
          tabPanel('2D-PCA',
                   column(12,tags$br(),imageOutput("pca2d")
                   )
          ),
          tabPanel('3D-PCA',
                   column(12,tags$br(),plotlyOutput("pca3d")
                   )
          )
        )
  #      plotOutput("distPlot")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         #plotOutput("distPlot"),br(),
        br(),
        tabsetPanel(id = 'lncRNAs',
                    tabPanel(
                      'Differentially Expressed lncRNAs',
                      column(
                        12,
                        br(),
                        checkboxInput("selcols", tags$b("Modify Columns:")),
                        conditionalPanel(
                          condition = "input.selcols == true",
                          checkboxGroupInput(
                            'sel_cols',
                            'Columns to Display:',
                            colnames(data),
                            selected = colnames(data)
                          )
                        ),
                        actionButton('clear', 'Clear selection'),
                        downloadButton('tbl_dl', "CSV"),
                        DT::dataTableOutput('tbl'),
                        hr(),
                        verbatimTextOutput('gene_s'),
                        htmlOutput("gene")
                      )
                    ))
         
      )
   )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, session) {
    
      data_f <- reactive({
           if(!input$exclude){
              data_tmp = data[log2FoldChange>input$lfc[1]&log2FoldChange<input$lfc[2],]
              data_tmp = data_tmp[padj<input$padj,]
           }else{
              data_tmp = data[!(log2FoldChange>input$lfc[1]&log2FoldChange<input$lfc[2]),]
              data_tmp = data_tmp[padj<input$padj,]
           }

           list(data_tmp=data_tmp)
      })

      output$pca2d <- renderImage({
        list(src = "../output/pca.png",
             contentType = 'image/png',
             width = "85%",
             height = "100%",
             alt = "loading 2D-PCA")
      }, deleteFile = FALSE)
      
      
      output$pca3d <- renderPlotly({
        pca=read.table("../output/pca.txt",header=T,sep="\t")
        plot_ly(
          pca,
         # alpha = 0.5,
          x =  ~ PC1,
          y = ~ PC2,
          z = ~ PC3,
          color = ~ group,
          text = ~ sample,
          # size = 10000,
          colors = c('#0C4B8E', '#BF382A'),
          opacity = 0.7
        ) %>% layout(title = '<br>PCA')
      })  
      
      
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      a <- data$log2FoldChange
      min=trunc(min(data$log2FoldChange))-1
      max=trunc(max(data$log2FoldChange))+1
      a1=melt(table(cut(a,breaks=min:max,ordered=TRUE)))
      a2<-data.frame(sapply(a1,function(x) gsub("\\(|\\]","",gsub("\\,","~",x))))
      colnames(a2)<-c("log2Fold_Change","# Differentially expressed genes")
      a2$`# Differentially expressed genes`=as.numeric(levels(droplevels(a2$`# Differentially expressed genes`)))
      p=ggplot(a2,aes(x=factor(log2Fold_Change),y=`# Differentially expressed genes`,fill=factor(log2Fold_Change)))+geom_bar(stat="identity")+coord_flip()+theme_bw()
      p
   })
  # https://yihui.shinyapps.io/DT-proxy/
  # reset tab
   output$tbl = DT::renderDataTable(
     data_f()$data_tmp[, input$sel_cols, with = F],
     selection = list(target = "row", selected = c(1)),
     server = TRUE,
     extensions = 'Buttons',
     class = 'compact',
     filter = 'top',
     caption = '',
     options = list(
       # dom = 'Bfrtip',
       # buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
       pageLength = 5,
       lengthMenu = c(5, 10, 15, 20),
       autoWidth = TRUE
     ),
     rownames = FALSE
   )
   
   output$tbl_dl <- downloadHandler(
     filename = function() {paste0(paste("lncRNA_DEGs",input$lfc[2],input$lfc[1],input$padj,sep="_"),".csv") },
     content = function(file) {
       write.csv(data_f()$data_tmp, file, row.names = F)
     }
   )
   
   proxy_tbl = dataTableProxy('tbl')

   observeEvent(input$clear, {
     proxy_tbl %>% selectRows(NULL)
   })
   
#   output$rbp_tbl = DT::renderDataTable(rbp,extensions='Buttons',class = 'compact',filter='top',caption = '',options = list(dom='Bfrtip',buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),pageLength = 5,lengthMenu = c(5, 10, 15, 20), autoWidth = TRUE), rownames= FALSE)

   output$gene_s = renderPrint({
        s=input$tbl_rows_selected
        cat('These genes were selected (max:10):\n\n')
        if (length(s)<10){
                g=data_f()$data_tmp$gene[s][1:length(s)]

        }else{
                g=data_f()$data_tmp$gene[s][1:10]

        }
        cat(g, sep = ',')
   })


   output$gene = renderUI({
        s=input$tbl_rows_selected
        jobid= gsub("/data/lancer_jobs/|/deg","",getwd())
        #cat('These genes were selected (max:10):\n\n')
	if (length(s)<10){
        	g = data_f()$data_tmp$gene[s][1:length(s)]
		
	}else{
		g = data_f()$data_tmp$gene[s][1:10]
               
	}
        # cat(g, sep = ',')

actionButton("cs", "Calculate co-expressed genes",style="color: #fff; background-color: #337ab7; border-color: #2e6da4",icon = icon("th"),onclick = paste0("window.open('",jobid,"&genelist=",paste0(g,collapse=","),"','_blank')"))
   })   

#   output$com_Ray2013 <- renderImage({
#        s=input$tbl_rows_selected
#        #cat('These genes were selected (max:10):\n\n')
#        if (length(s)<10){
#                g = data_f()$data_tmp$gene[s][1:length(s)]
#
#        }else{
#                g = data_f()$data_tmp$gene[s][1:10]

#        }       
#	list(src = paste0("../output/com_Ray2013_",paste(g,collapse="_"),".png"),
#                       contentType = 'image/png',
#                       width = 800,
#                       height = 300,
#                       alt = "This is alternate text")
#        },deleteFile = FALSE)


})
   
# Run the application 
shinyApp(ui = ui, server = server)

