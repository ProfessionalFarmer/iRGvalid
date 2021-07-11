if(!require(shiny)){BiocManager::install("shiny")}
if(!require(dplyr)){BiocManager::install("dplyr")}
if(!require(data.table)){BiocManager::install("data.table")}
if(!require(ggpubr)){BiocManager::install("ggpubr")}

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
                                choices = list(`TCGA`=list("COAD", "LUAD", "BRCA"), `GSE102349`=list("NPC"), `ALL (COAD, LUAD, BRCA and NPC)`=list("ALL") ), 
                                selected = 1),
                    sliderInput("IRGvalid.study.panel.size", h3("Select panel size"),
                                min = 1, max = 10, value = c(5, 6)),
                    sliderInput("IRGvalid.study.mincor", h3("Range of correlation (Rt) value"),
                                min = 0, max = 1,
                                value = c(0.95,1), step = 0.01),
                ),
                mainPanel(
                    #width = 16,
                    p("This page only shows the correlation analysis results of target genes (HER2, MMP11, MYBL2 for TCGA-BRCA, NOTCH2, BRCA1, PDC for TCGA-COAD, HLA-A, ERBB3, HIF1A for TCGA-LUAD, and HLA-A, FNDC3B, ANXA1 for GSE102349-NPC) by candidate reference genes (HNRNPL, PCBP1, RER1, SF3A1, CIAO1, SRSF4, CNBP, MYL12B, UBC, TMBIM6, RPS27, EIF1, GAPDH, HEY1) reported in our manuscript."),
                    br(),
                    p("If you would like to analyze a specific reference gene panel using any target gene(s), please click 'Explore' tab."),
                    br(),
                    br(),
                    # Output: Verbatim text for data summary ----
                    #verbatimTextOutput("iRGvalid.description"),
                    
                    dataTableOutput("iRGvalid.result")
                )
            )
            
          
            
        ),
        tabPanel(
            "Explore",
            sidebarLayout(
                sidebarPanel(
                    
                    checkboxGroupInput("uc.project", h3("Select dataset(s)"),
                                choiceNames  = c("BLCA", "BRCA", "CESC", "COAD", "CHOL", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PCPG", "PRAD", "READ", "SARC", "STAD", "THCA", "THYM", "UCEC"),
                                choiceValues = c("TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-COAD", "TCGA-CHOL", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-STAD", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC"),
                                selected = "TCGA-COAD", inline = T),

                    textAreaInput("uc.panel",h3("Enter candidate reference gene(s)"), "CIAO1\nCNBP\nHEY1\nUBC", height = "150px"),
                    textInput("uc.target", h3("Enter target gene to normalize"), "HLA-A"),
                    actionButton("uc.showexpression", h4("Show expression")),
                    actionButton("uc.analysis", h4("Perform analysis"))
                ),

                
                mainPanel(
                    h3("Expression levels in log2(TPM). User can explore other panels."),
                    htmlOutput("uc.description"),
                    dataTableOutput("uc.expression.table"),
                    
                    # Output: Plot
                    column(12,plotOutput("uc.before.correlation")),
                    column(12,plotOutput("uc.after.correlation")),
                    
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
               "NPC" = all.projects.comb[,c(1,2,3,6,10)],
               "ALL" = all.projects.comb)
    })
    
    # output$iRGvalid.description = renderPrint({ input$IRGvalid.study.panel.size[1] })

    
    # Show 
    output$iRGvalid.result <- renderDataTable({
        datasetInput() %>% 
            filter(`Panel size` >= input$IRGvalid.study.panel.size[1] & `Panel size` <= input$IRGvalid.study.panel.size[2] ) %>%
            mutate(`Mean correlation` = rowMeans(dplyr::select(., 3:length(.)))  ) %>%
            filter(`Mean correlation` >= input$IRGvalid.study.mincor[1] & `Mean correlation` <= input$IRGvalid.study.mincor[2] )
        
    }, options = list(autoWidth = F, 
                      pageLength = 20, 
                      columnDefs = list(list( targets = c(1), width = '600px')) 
                      ) 
    )
    

    ############### User customed
    
    # For test only
    # output$uc.description = renderPrint({ c(input$uc.project, input$uc.target, input$uc.panel) })

    observeEvent(input$uc.showexpression, {
        referenece.genes = parseReferenceGene(input$uc.panel) 

        
        p.df.list <- lapply(input$uc.project, function(p){
            getGeneExpressionDf(p, input$uc.target, referenece.genes)
        })
        gene.df <- do.call(dplyr::bind_rows, p.df.list)
        

        output$uc.expression.table <- renderDataTable({
            #getGeneExpressionDf(input$uc.project, input$uc.target, referenece.genes) %>% 
            gene.df %>%
                tibble::add_column(Sample=rownames(.), .before = 1)
        }, options = list(autoWidth = F, 
                          pageLength = 10, 
                          columnDefs = list(list( targets = c(0), width = '300px')) 
                        ) 
        )
        
    })
                     
    
    observeEvent(input$uc.analysis, {
        target.gene = input$uc.target
        referenece.genes = parseReferenceGene(input$uc.panel) 

        if(length(referenece.genes)>10){
            showNotification("Please input no more than 10 genes.\nOne gene per line.")
        }
        
        # obtain gene expressino data
        # gene.exp.df <- getGeneExpressionDf(input$uc.project, target.gene, referenece.genes)
        
        p.df.list <- lapply(input$uc.project, function(p){
            getGeneExpressionDf(p, input$uc.target, referenece.genes)
        })
        gene.exp.df <- do.call(dplyr::bind_rows, p.df.list)
        gene.exp.df <- na.omit(gene.exp.df)
        
        # 有些基因被过滤掉了，getGeneExpressionDf函数没有考虑不在数据框的reference，所以这里更新reference
        referenece.genes <- intersect(referenece.genes, colnames(gene.exp.df))
        
        # normalizing
        norm.res <- normalize(gene.exp.df, target.gene, referenece.genes)
        
        
        # correlation analysis
        output$uc.before.correlation <- renderPlot( corrplot(norm.res$correlation.before, type = "upper", addCoef.col = "black", order = "hclust", title = paste0("Correlation between ", target.gene," (raw value) and reference genes"), mar=c(0,0,1,0))   )
        output$uc.expression.table <- NULL
        
        
        output$uc.beforeAndAfter.expression = renderDataTable(norm.res$target.df, options = list(pageLength = 10)  ) 
       
        output$uc.after.correlation <- renderPlot( norm.res$scatter + xlab(paste0("Raw value - ",target.gene) )  +  ylab(paste0("Normalized value - ",target.gene) ) )
        
        
        summary.word <- paste0(
            "The target genes is normalized against the combination of all the candidate gene(s) (", paste0(referenece.genes,sep="",collapse = ", "), "), ",
            "and the correlation (Rt value) is <font color=\"red\"><h3>", 
            norm.res$normalize.correlation.value, "</h3></font>\n--------",
            "<br>1. The correlations between target gene and reference genes are shown in the first figure.",
            "<br>2, The expression levels of pre- and post-normalization (normalized by ", paste0(referenece.genes,sep="",collapse = ", "),") are shown in the second figure. Detailed values are shown in the first table.",
            "<br>3, Ranked results of all reference gene combinations are listed in the second table. The best combination (<b>", unlist(norm.res$all.combination[1,3]), "</b>) has a Rt value of <b>", unlist(norm.res$all.combination[1,2]), "</b>.<br><br>"
        )
        
        output$uc.description = renderText({ summary.word  })
        
        output$uc.allcombination.res = renderDataTable( norm.res$all.combination, options = list(pageLength = 10) )
        
        
    })
    
    
}

# Run the application
shinyApp(ui = ui, server = server)


