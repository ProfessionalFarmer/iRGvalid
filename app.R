library(shiny)
library(dplyr)
library(data.table)

# Shiny UI -------
ui <- fluidPage(
    theme = "bootstrap.css",
    includeCSS("www/styles.css"),
    
    navbarPage(
        "iRGvalid: A Robust In Silico Method for Optimal Reference Gene Validation",
        id = "main_navbar",
        
        tabPanel(
            "IRGvalid Result",
            sidebarLayout(
                sidebarPanel(
                    width = 4,
                    # Input: Slider for the number of bins ----
                    selectInput("IRGvalid.study.select", h3("Select result from our manuscript"), 
                                choices = list(`TCGA`=list("COAD", "LUAD", "BRCA"), `GSE102349`=list("NPC") ), 
                                selected = 1),
                    sliderInput("IRGvalid.study.panel.size", h3("Select panel size"),
                                min = 1, max = 10, value = c(5, 6)),
                    sliderInput("IRGvalid.study.mincor", h3("Choose minimum correlation cutoff"),
                                min = 0, max = 1,
                                value = c(0.95,1), step = 0.01),
                ),
                mainPanel(
                    width = 8,
                    p("This page only shows the correlations between target genes (HER2, MMP11, MYBL2 for TCGA-BRCA, NOTCH2, BRCA1, PDC for TCGA-COAD, HLA-A, ERBB3, HIF1A for TCGA-LUAD, and HLA-A, FNDC3B, ANXA1 for GSE102349-NPC) and reference panel"),
                    br(),
                    # Output: Verbatim text for data summary ----
                    #verbatimTextOutput("iRGvalid.description"),
                    
                    tableOutput("iRGvalid.result")
                )
            )
            
          
            
        ),
        tabPanel(
            "User customed",
            sidebarLayout(
                sidebarPanel(
                    selectInput("uc.project", h3("Select result from our manuscript"), 
                                choices = c("TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-STAD", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC"), 
                                selected = 1),
                    textInput("uc.target", h3("Target gene to normalize"), "HLA-A"),
                    textAreaInput("uc.panel",h3("Expected reference gene panel"), "CIAO1\nCNBP\nHEY1\nUBC", height = "350px"),
                    actionButton("uc.showexpression", h4("Show expression")),
                    actionButton("uc.analysis", h4("Perform analysis"))
                ),
                mainPanel(
                    h3("Expression levels in log2(TPM)"),
                    htmlOutput("uc.description"),
                    dataTableOutput("uc.expression.table"),
                    
                    # Output: Plot
                    column(12,plotOutput("uc.before.correlation")),
                    
                    # column(12,plotOutput("uc.after.correlation")),
                    dataTableOutput("uc.beforeAndAfter.expression"),
                    
                    h3("All combination result"),
                    dataTableOutput("uc.allcombination.res")
                )
            )
            
        )
    )
)


load("./all.combination.rdata")
source("./functions.R")


# Shiny Server Side -------
server <- function(input, output, session) {
    
    ####### iRGvalid result
    # Return the requested dataset ----
    datasetInput <- reactive({
        switch(input$IRGvalid.study.select,
               "COAD" = all.projects.comb[,c(1,2,4,13,14)],
               "LUAD" = all.projects.comb[,c(1,2,5,8,9)],
               "BRCA" = all.projects.comb[,c(1,2,7,11,12)],
               "NPC" = all.projects.comb[,c(1,2,7,11,12)])
    })
    
    # output$iRGvalid.description = renderPrint({ input$IRGvalid.study.panel.size[1] })

    
    # Show 
    output$iRGvalid.result <- renderTable({
        datasetInput() %>% 
            filter(`Panel size` >= input$IRGvalid.study.panel.size[1] & `Panel size` <= input$IRGvalid.study.panel.size[2] ) %>%
            mutate(`Mean correlation` = rowMeans(dplyr::select(., 3:length(.)))  ) %>%
            filter(`Mean correlation` >= input$IRGvalid.study.mincor[1] & `Mean correlation` <= input$IRGvalid.study.mincor[2] )
    })
    
    
    
    ############### User customed
    
    output$uc.description = renderPrint({ c(input$uc.project, input$uc.target, input$uc.panel) })
    

    observeEvent(input$uc.showexpression, {
        output$uc.expression.table <- renderDataTable({
            getGeneExpressionDf(input$uc.project, input$uc.target, input$uc.panel) %>% 
                tibble::add_column(Sample=rownames(.), .before = 1)
        })
    })
                     
    
    
    observeEvent(input$uc.analysis, {
        target.gene = input$uc.target
        referenece.genes = unlist( strsplit(input$uc.panel,"\n") )
        
        # obtain gene expressino data
        gene.exp.df <- getGeneExpressionDf(input$uc.project, target.gene, referenece.genes)
        
        # correlation analysis
        library(corrplot)
        correlation.before <- round(cor(gene.exp.df[,-c(ncol(gene.exp.df))], method = c("pearson")  ), 3) # spearman pearson
        
        output$uc.before.correlation <- renderPlot( corrplot(correlation.before, type = "upper", addCoef.col = "black", order = "hclust", title = paste(input$uc.target," - before normalization"), mar=c(0,0,1,0))   )
        output$uc.expression.table <- NULL
        
        # normalizing
        norm.res <- normalize(gene.exp.df, target.gene, referenece.genes)
        
        output$uc.beforeAndAfter.expression = renderDataTable(norm.res$target.df)
        
        summary.word <- paste0(
            "Final correlation (Rt value) between target gene and referenece panel is <font color=\"red\"><h3>", 
            norm.res$correlation, "</h3></font>\n--------",
            "<br>1. The correlations between target gene and reference genes are shown",
            "<br>2, The expression levels before and after are shown as below",
            "<br>3, We also show the result of all the combination"
        )
        
        output$uc.description = renderText({ summary.word  })
        
        
        output$uc.allcombination.res = renderDataTable(normalize.all.combination(gene.exp.df, target.gene, referenece.genes))
        
        
    })
    
    
}

# Run the application
shinyApp(ui = ui, server = server)


