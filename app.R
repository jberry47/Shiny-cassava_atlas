library("png")
library("grid")
library("ggplot2")
library(shiny)
library(shinyBS)
library(shinydashboard)
library(DT)
library(stringr)
library(stats)
library(mailR)
library(shinyjs)
library(reshape2)
library(plyr)

all.rna <- read.table("data/Rshiny_app_dataset_v2.txt",header=T,quote="",sep="\t")

rnaseq <- data.frame(all.rna$gene,all.rna$FEC,all.rna$Fibrous.Root,all.rna$Lateral.Bud,all.rna$Leaf,all.rna$Mid.Vein,all.rna$OES,all.rna$Petiole,all.rna$RAM,all.rna$SAM,all.rna$Stem,all.rna$Storage.Root)
colnames(rnaseq) <- c("gene","FEC",'Fibrous.Root','Lateral.Bud','Leaf','Mid.Vein','OES','Petiole','RAM','SAM','Stem','Storage.Root')

genefinder_data1 <- read.table("data/Rshiny_app_dataset_v3_withreps.txt",header=T,quote="",sep=",")
genefinder_data1 <- genefinder_data1[,1:39]
genefinder_data1$go <- all.rna$go
genefinder_data <- data.frame(genefinder_data1[,34],genefinder_data1[,1:32])
names(genefinder_data) <- c("gene",names(genefinder_data)[-1])


annot <- data.frame(seq(from=1,to=nrow(all.rna),by=1),all.rna$gene,all.rna$locus,all.rna$go,all.rna$annot)
colnames(annot) <- c("index","gene","locus","go","annot")

overlap <- read.table("data/overlapping.txt")
issues <- data.frame(seq(from=1,to=nrow(all.rna),by=1),all.rna$gene,all.rna$possible_issues)
issues$overlap<- all.rna$gene %in% overlap$V1
colnames(issues) <- c("index","gene","possible_issues","overlap")
issues$possible_issues <- as.logical(issues$possible_issues)
annot$status <- as.character(issues[,3] | issues[,4])

pcascores <- read.csv("data/som_pca.csv",header=T)
pcascores$node <- all.rna$node
pcascores$node <- as.factor(pcascores$node)

test1 <- readPNG("data/cassava_both_annot.png")
rtest <- rasterGrob(test1,interpolate = T)

starting_points <- data.frame(rbind(c(0,0),c(1,1)))
my_points <- data.frame(rbind(c(.665,0.095),c(.255,0.267),c(.255,0.92),c(.255,0.75),c(.255,0.641),c(.325,0.095),c(.255,0.532),c(.737,0.375),
  c(0.737,0.9),c(0.737,0.635),c(.255,0.425)))

rownames(my_points) <- c("FEC",'Fibrous.Root','Lateral.Bud','Leaf','Mid.Vein','OES','Petiole','RAM','SAM','Stem','Storage.Root')

pca_plot <- ggplot(pcascores[pcascores$node!=0,], aes(PC1, PC2,colour=node))+
  theme_bw() +
  geom_point(alpha=0.03)+
  theme(legend.position="none")

load("data/des.RData")
load("data/cd.RData")
load("data/bg.RData")
load("data/at.RData")

#********************************************************************************************************************
#********************************************************************************************************************
#********************************************************************************************************************

ui <- dashboardPage(skin="purple",
  dashboardHeader(
    title = "Bart Lab Cassava Atlas",
    titleWidth = 450
  ),
  dashboardSidebar(
    disable = T
  ),
  dashboardBody(
    tags$style(HTML("
      .tabbable > .nav > li > a[data-value='About'] {background-color: #676767;   color:white}
      .tabbable > .nav > li > a[data-value='Cas-Xam'] {background-color: #676767;   color:white}
      .tabbable > .nav > li > a[data-value='Plot'] {background-color: #676767;   color:white}
      .tabbable > .nav > li > a[data-value='Tissue-Specific Heatmap'] {background-color: #676767;  color:white}
      .tabbable > .nav > li > a[data-value='Gene Finder'] {background-color: #676767; color:white}
      .tabbable > .nav > li > a[data-value='Contact Us'] {background-color: #676767; color:white}
      .tabbable > .nav > li[class=active]    > a {background-color: #444444; color:white}
      .multicol{
      -webkit-column-count: 4; /* Chrome, Safari, Opera */
      -moz-column-count: 4; /* Firefox */
      column-count: 4;
      }
      .twocol{
      -webkit-column-count: 2; /* Chrome, Safari, Opera */
      -moz-column-count: 2; /* Firefox */
      column-count: 2;
      }
      .warning { 
      color: red;
      }"
    )
    ),
    tabsetPanel(
      tabPanel("About",
        br(),
        fluidRow(
          box(title="Welcome to Cassava Atlas",width=8,solidHeader = T,status = 'success',
            h3("Background"),
            p("This tool is brought to you by members of the",a(href="http://www.beckybartlab.org/",target='_blank',"Bart Lab"),"at the Donald Danforth Plant Science 
              Center in St. Louis MO. For more details on this tool and the RNAseq data used to generate 
              these analyses please read our paper",a("found here",target="_blank",href="http://onlinelibrary.wiley.com/doi/10.1111/nph.14443/full"),"."),
            p("Below you will find information on how to use two different gene expression tools that can
              be accessed through the tabs above. The Tissue-Specific Heatmap tool is designed to allow 
              users to visualize gene expression patterns for any gene of interest. The GeneFinder tool 
              allows users to specify a desired gene expression pattern and retrieve candidates. Questions 
              should be directed to jberry@danforthcenter.org and/or rbart@danforthcenter.org"),
            p("This webpage is created using open-source R Shiny developement tools and the souce code to this page is available",a("here",target='_blank',href="cassava_atlas_source.txt")),
            p("This work was funded by the Bill and Melinda Gates Foundation.")
            )
            )
        ,
        fluidRow(
          box(width=8,title="How to: ",solidHeader = T,status = "success",
            box(width=6,title = "Use the heatmap",collapsible = TRUE, collapsed = TRUE,
              p("This section is for when there is already a gene of interest and the question
                becomes where that gene is expressed. There are a couple different ways to find
                the gene to produce the heatmap (See searching). After the heatmap is generated, there a few customizable
                options available for the scale bar. The transform option will transform the FPKM's 
                for each tissue type and produce a new heatmap automatically with the new transformed
                values. The scale max option allows the user to set a maximum value to the color 
                gradient so the distribution of colors can be seen more easily. This is especially helpful when 
                the gene of interest is very highly expressed in one tissue type and the default
                color gradient becomes skewed to accommodate the large FPKM value."),
              hr(),
              h4("Searching"),
              p("Multiple terms separated by a comma are allowed"),
              h4("Manes"),
              p("* If the Manes name is already known, that is the 
                most efficient way to select the gene of interest. It should be noted that the name must begin with a 
                capital M and end with a period followed by the identification. For example Manes.02G126300"),
              hr(),
              h4("GO ID"),
              p("* You can also search for families of genes using the GO ID. It should be noted that the GO
                term must begin with capitals GO followed by a colon. For example GO:0008270"),
              hr(),
              h4("General"),
              p("* The least selective search is parsing keywords over the annotation. If the gene of interest
                is a zinc finger but the GO ID and Manes name are not known, entering zinc finger into the 
                search bar will return all genes that have zinc finger in the annotation.")
              )
            ,
            box(width=6,title="Use the Gene Finder",collapsible = TRUE,collapsed = TRUE,
              p("This tool allows the user to identify genes following a specified expression pattern. There 
                are three options for each tissue type that can be defined by the user with the criterion box. 
                If the cutoff for saying the gene is off is equivalent to FPKM < 1, enter 1 in that field and 
                similarly for equivalence for on, set the FPKM > 15. If the status for the tissue type is set 
                to either, a search is not done over that tissue type. In the end, the returned list of genes 
                will fit the criterion set by the user and the full dataset can be downloaded in CSV format."),
              hr(),
              h4("Secondary Search"),
              p("* After the expression pattern search is done and the table is generated, there is an option
                to further parse the data to extract only genes of interest. For example if the expression
                pattern of interest is on in the leaf blade and off in the petiole, a list if 1,711 genes
                is returned. A secondary parse of 'mitochondrial' will return 7 genes that have the
                expression pattern and the keyword in the annotation."),
              
              hr(),
              h4("Transferring to heatmap"),
              p("* There is developement in progress that will allow easy transfer for a gene found using
                the finder by a simple button click. However, the only way to do this currently is to copy 
                and paste the Manes name of the found gene into the search bar in the heatmap tab and go 
                from there.")
              ))
              )
              ),
      tabPanel("Tissue-Specific Heatmap", 
        fixedRow(
          br(),
          box(style = "overflow-y:scroll",width=10,title = "Generation Guide",solidHeader = T,status = 'success',collapsible = TRUE,
            h4("Step 1: Enter Search Terms"),
            column(5,
              p("Multiple terms separated by a comma is allowed"),
              textInput('gene_search_text', label="ex: zinc finger,b-box, GO:0005515", value = "",width=550),
              tags$div(class = "multicol", radioButtons("heat_radio", "Logical",c("And" = "AND", "Or"="OR"),selected = 'AND'))
            ),
            br(),
            br(),
            br(),
            actionButton("gene_search_button",label="Search",inline=TRUE),
            h6("To search, click the search button"),
            br(),
            hr(),
            h4("Step 2: Select a gene by clicking it"),
            dataTableOutput("gene_search_table"),
            uiOutput("download_found_genes"),
            uiOutput('is.problem'),
            hr(),
            h4("Step 3: Make Heatmap"),
            uiOutput('ready_for_heat')
          )
        ),
        uiOutput('heat_out')
      ),
      tabPanel("Gene Finder",
        br(),
        fluidRow(
          box(width=6,
            title = "Select expression status for tissue types",
            solidHeader = TRUE,
            status='success',
            wellPanel(
              tags$div(class = "multicol", radioButtons("leaf_sel", "Leaf",c("On" = "on", "Off"="off","Either"="either"),selected = 'either')),
              tags$div(class = "multicol",radioButtons("petiole_sel", "Petiole",c("On" = "on", "Off"="off","Either"="either"),selected = 'either')),
              tags$div(class = "multicol",radioButtons("mv_sel", "Mid Vein",c("On" = "on", "Off"="off","Either"="either"),selected = 'either')),
              tags$div(class = "multicol",radioButtons("lb_sel", "Lateral Bud",c("On" = "on", "Off"="off","Either"="either"),selected = 'either')),
              tags$div(class = "multicol",radioButtons("stem_sel", "Stem",c("On" = "on", "Off"="off","Either"="either"),selected = 'either')),
              tags$div(class = "multicol",radioButtons("sr_sel", "Storage Root",c("On" = "on", "Off"="off","Either"="either"),selected = 'either')),
              tags$div(class = "multicol",radioButtons("fr_sel", "Fibrous Root",c("On" = "on", "Off"="off","Either"="either"),selected = 'either')),
              tags$div(class = "multicol",radioButtons("fec_sel", "FEC",c("On" = "on", "Off"="off","Either"="either"),selected = 'either')),
              tags$div(class = "multicol",radioButtons("oes_sel", "OES",c("On" = "on", "Off"="off","Either"="either"),selected = 'either')),
              tags$div(class = "multicol",radioButtons("sam_sel", "SAM",c("On" = "on", "Off"="off","Either"="either"),selected = 'either')),
              tags$div(class = "multicol",radioButtons("ram_sel", "RAM",c("On" = "on", "Off"="off","Either"="either"),selected = 'either'))
            ),
            uiOutput("make_goButton")
          ),
          box(width = 2,
            title="Search Criterion",
            solidHeader = T,
            status = 'success',
            HTML('<table><tr><td>FPKM <   </td><td><textarea style="background-color: white;border-style:solid;" id="fpkmOff" rows="1" cols="1">2</textarea></td><td> = OFF</td></tr></table>'),
            HTML('<table><tr><td>FPKM >   </td><td><textarea style="background-color: white;border-style:solid;" id="fpkmOn" rows="1" cols="1">15</textarea></td><td> = ON</td></tr></table>')
          )),
        uiOutput("table_out")
      ),
      tabPanel("Cas-Xam",
        br(),
        fluidRow(
          box(width=10,title = "Select Subset",solidHeader = T,status = 'success',collapsible = TRUE,
              tabsetPanel(id = "tabs",
                          tabPanel("Subset by Treatment", value = "treatment_panel",
                                   br(),
                                   column(width = 6,
                                         wellPanel(
                                           column(4,
                                                  tags$div(radioButtons("mock8hr", "8hr: Mock",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')),
                                                  tags$div(radioButtons("xam6688hr", "8hr: Xam668",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')),
                                                  tags$div(radioButtons("xe8hr", "8hr: Xe",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')),
                                                  tags$div(radioButtons("xam668xe8hr", "8hr: Xam886 + Xe",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude'))
                                           ),
                                           column(4,
                                                  tags$div(radioButtons("mock24hr", "24hr: Mock",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')),
                                                  tags$div(radioButtons("xam66824hr", "24hr: Xam668",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')),
                                                  tags$div(radioButtons("xe24hr", "24hr: Xe",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')),
                                                  tags$div(radioButtons("xam668xe24hr", "24hr: Xam886 + Xe",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')) 
                                           ),
                                           column(4,
                                                  tags$div(radioButtons("mock50hr", "50hr: Mock",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')),
                                                  tags$div(radioButtons("xam66850hr", "50hr: Xam668",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')),
                                                  tags$div(radioButtons("xe50hr", "50hr: Xe",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')),
                                                  tags$div(radioButtons("xam668xe50hr", "50hr: Xam886 + Xe",c("Include" = "include", "Exclude"="exclude"),selected = 'exclude')) 
                                           ),
                                         numericInput("fc_cut", "Fold-change Cut-off", value = 2, width = 150),
                                         actionButton("casxam_sub_table_button", "Subset")
                                         )
                                     )
                          ),
                          tabPanel(value = "gene_panel", "Subset by Gene",
                                   br(),
                                   textInput('casxam_gene_search_text', label="ex. Manes.06G123400", value = "",width=350),
                                   actionButton("casxam_gene_search_button",label="Subset")
                          )
              )
          ),
          uiOutput("cdbg_sub"),
          uiOutput("casxam_gene_sub")        
        )
      ),
      tabPanel("Contact Us",
        br(),
        fluidRow(
          box(width = 6,style = "overflow-y:scroll",
            title = "Send us an email",
            solidHeader = T,
            status = 'success',
            textInput("name","Name",width = 400),
            textInput("from","Preferred Email",width = 400),
            textInput("subject","Subject",width = 400),
            selectInput("who", width = 400,
              label = "Topic",
              choices = c('General Inquiry', 'Problems using tool'),
              selected = "General Inquiry"),
            p(strong("Message")),
            HTML('<textarea style="background-color: white;border-style:solid;" id="body" rows="15" cols="80"></textarea>'),
            br(),br(),
            shinyjs::useShinyjs(),
            actionButton('sendEmailButton',"Send",icon("envelope-o")),
            uiOutput("sentSuccessful")
          )
        )
      )
              )
              )
        )

#********************************************************************************************************************
#********************************************************************************************************************
#********************************************************************************************************************

server <- function(input, output) {
  
  observeEvent(input$sendEmailButton, {
    is.error <- FALSE
    if(input$name=="" | input$from == "" | input$subject == "" | input$body == ""){
      is.error <- TRUE
    }
    output$sentSuccessful <- renderUI({
      if(is.error == TRUE){
        tags$div(class = "warning", p("Please include your name, email, a subject, and a message"))
      }else{h3("Thanks so much! We'll get back to you as soon as possible")}
    })
    if(is.error == FALSE){
      shinyjs::disable('sendEmailButton')
      sender <- "ddpsc_cassavaapp@yahoo.com"
      recipients <- switch(input$who,
        'General Inquiry'='RBart@danforthcenter.org',
        'Problems using tool'='RBart@danforthcenter.org')
      send.mail(from = sender,
        to = 'RBart@danforthcenter.org',
        subject = input$subject,
        cc= "jberry@danforthcenter.org",
        body = paste("from:",input$name,"\n","email:",input$from,"\n","topic:",input$who,"\n\n","message:","\n",input$body),
        smtp = list(host.name = "smtp.mail.yahoo.com", port = 465,
          user.name = "ddpsc_cassavaapp@yahoo.com",
          passwd = "STATS!jb47", ssl = TRUE),
        authenticate = TRUE,
        send = TRUE)
    }
  }
  )
  
  #********************************************************************************************************************
  # for heatmap
  #********************************************************************************************************************
  output$download_found_genes <- renderUI({
    s <- gene_search()
    downloadButton("get_all_data","Get all data for these genes (CSV)")
  })
  output$is.problem <- renderUI({
    s <- gene_sel()
    if(length(s)==0){
      return()
    }else{
      sub <- issues[as.numeric(s),]
      if(sub$possible_issues == 'True'){
        tags$div(class = "warning",p("WARNING: During assembly, multiple gene loci aligned to this gene name"))
      }
      if(sub$overlap == TRUE){
        tags$div(class = "warning",p("WARNING: This gene is called with low confidence"))
      }
      else{
        return()
      }
    }
    
  })
  
  output$ready_for_heat <- renderUI({
    s <- gene_search()
    if(s==""){
      return()
    }else if(is.null(input$gene_search_table_rows_selected)){
      return()
    }
    else{
      actionButton("make_heatmap_gene","Make heatmap with selected gene")
    }
  })
  
  gene_sel <-reactive({
    sub <- heat_sub_data()
    sub[input$gene_search_table_rows_selected,"index"]
  })
  
  make_heat<- eventReactive(input$make_heatmap_gene,{
    gene_sel()
  })
  
  output$heat_out <- renderUI({
    make_heat()
    fixedRow(
      column(4,
        box(width=12,title="Download Raw Data",solidHeader = T,status = "success",
          downloadButton("downloadHeatFPKM","Download All Data for Selected Gene (CSV)")
        ),
        box(width=12,title = "Heatmap Adjustments",status='success',solidHeader = T,
          tags$div(class = "multicol", radioButtons("color_sel", "Color for Heatmap Gradient",c("Blue" = "blue", "Red"="red","Purple"="purple","Orange"='orange',"Green"="green","Magenta"='magenta'),selected = 'green')),
          hr(),
          tags$div(class = "multicol", radioButtons("transform", "Transform",c('Identity'='identity',"Log10" = "log10", "Log2"="log2"),selected = 'identity')),
          hr(),
          textInput("scale_max","Scale Max",value="",width = 100)
        )
        ,
        box(width=12,title = "Tissue Expression PCA",status='info',solidHeader = T,
          tabsetPanel(
            tabPanel("Plot",
              plotOutput("plot2"
                ,height=350
                #width = 330)
              ),
              hr(),
              uiOutput("is.zero"),
              downloadButton("downloadDiversity","Download Expression PCA (PNG)")
            ),
            tabPanel("About",
              br(),
              p("This is a representation of the expression diversity over all the tissue types
                examined in this analysis. This plot is the result of principle component analysis (PCA)
                over all the tissues followed by clustering using self organizing maps (SOM). PC1 and PC2 are plotted against each other
                and the color coding is based on SOM classification. The large spot in the plot is 
                the selected gene and the color of the spot denotes which cluster the gene belongs to."),
              p("There is one spot for every gene in the dataset and if two genes are close together in 
                this space, then they show similar expression patterns. For instance, they could both be
                on in the leaf blade and petiole, and off in the mid vein and lateral bud."),
              p("If a warning message is displayed, that means the gene selected did not make the criteria to be 
                included in the PCA plot. The main filter is if the gene is not differentially expressed in at least
                one pairwise comparison, then it is removed. The last is if the gene is not close enough to any node 
                created in the SOM.")
              )
              )
              )
              ),
      box(style = "overflow-y:scroll",width=6,title = "Tissue Specific Expression", status = "info", solidHeader = TRUE,
        column(12,
          p("Mouse over the circles to see mean FPKM value"),
          hr(),
          plotOutput("plot1", 
            height = 425,width = 670,
            click = "plot1_click",
            hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce")
          ),
          uiOutput("hover_info"),
          hr(),
          downloadButton("downloadPlot","Download Heatmap (PNG)")
        )
      ),
      box(width = 6,height=650,title = "FPKM Boxplot", status = "info", solidHeader = TRUE,
        plotOutput("heat_boxplot"),
        br(),
        br(),
        br(),
        br(),
        br(),
        hr(),
        downloadButton("heat_box","Download Boxplot (PNG)")
      )
            )
  })
  
  Heat_boxplot <- reactive({
    s <- make_heat()
    sub <- genefinder_data[as.numeric(s),]
    my_scale <- input$transform
    my_transform <- function(x,type){
      switch(type,
        "identity" = x,
        "log10" = log10(x+1),
        'log2' = logb(x+1,2)
      )
    }
    sub1 <- data.frame(my_transform(as.numeric(sub[1,-1]),my_scale))
    sub1$group <- rep(row.names(my_points),each=3)[-33]
    names(sub1) <- c("FPKM","Group")
    sub1$Group <- ordered(sub1$Group,levels=c("OES","FEC","Fibrous.Root","Storage.Root","RAM","SAM","Leaf","Lateral.Bud","Mid.Vein","Petiole","Stem"))
    ggplot(sub1, aes(x=Group, y=FPKM,fill=Group))+
      ggtitle(paste(sub$gene,": FPKM Distribution"))+
      geom_boxplot(size=1,fatten=0) + 
      theme_bw()+ 
      theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")  + 
      ylab(paste(my_scale,"FPKM"))+
      xlab("Tissue Type")+
      theme(axis.text = element_text(size = 14),
        axis.title= element_text(size = 18))+
      scale_fill_manual(values=c("#33a02c","#b2df8a","#1f78b4","#a6cee3","#6a3d9a","#cab2d6","#ffff99","#ff7f00","#fdbf6f","#e31a1c","#fb9a99"))
  })
  output$heat_boxplot <- renderPlot(height=500,{
    Heat_boxplot()
  })
  
  output$heat_box <- downloadHandler(
    filename = function() {paste(heatName(),'-Boxplot', '.png', sep='')},
    content=function(file){
      png(file,height = 500,width = 670)
      print(Heat_boxplot())
      dev.off()
    },
    contentType='image/png'
  )
  
  Plot1 <- reactive({
    s <- make_heat()
    sub <- rnaseq[as.numeric(s),]
    my_scale <- input$transform
    my_transform <- function(x,type){
      switch(type,
        "identity" = x,
        "log10" = log10(x+1),
        'log2' = logb(x+1,2)
      )}
    hr()
    sub[1,2:12] <- my_transform(as.numeric(sub[1,2:12]),my_scale)
    top <- input$scale_max
    if(suppressWarnings(is.na(as.numeric(top)))){
      top <- max(as.numeric(sub[1,2:12]))+1
    }
    for(i in 2:12){
      if(sub[1,i] > as.numeric(top)){
        sub[1,i] <- as.numeric(top)
      }
    }
    my_col <- input$color_sel
    ggplot(data = starting_points, aes(X1, X2)) +
      ggtitle(paste('\n',sub$gene,": Tissue-Specific Gene Expression"))+
      scale_x_continuous(limits=c(0,1))+
      scale_y_continuous(limits=c(0,1))+
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())+
      annotation_custom(rtest)+
      scale_colour_gradient(low='black',high=my_col,paste(my_scale,"FPKM"),limits=c(0, as.numeric(input$scale_max)))+
      geom_point(data=my_points,aes(X1,X2,colour=as.numeric(sub[1,2:length(sub)])),size=6)
  })
  
  output$plot1 <- renderPlot({
    Plot1()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {paste(heatName(),'-heatmap', '.png', sep='')},
    content=function(file){
      png(file,height = 425,width = 670)
      print(Plot1())
      dev.off()
    },
    contentType='image/png'
  )
  
  output$hover_info <- renderUI({
    hover <- input$plot_hover
    point <- nearPoints(my_points, hover, threshold = 10, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    
    # create style property fot tooltipleft_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
      "left:", left_px + 20, "px; top:", top_px +20, "px;")
    
    s <- make_heat()
    sub1 <- rnaseq[as.numeric(s),]
    my_scale <- input$transform
    my_transform <- function(x,type){
      switch(type,
        "identity" = x,
        "log10" = log10(x+1),
        'log2' = logb(x+1,2)
      )
    }
    sub1[1,2:12] <- my_transform(as.numeric(sub1[1,2:12]),my_scale)
    
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(signif(sub1[,rownames(point)],4)))
    )
  })
  
  output$downloadHeatFPKM <- downloadHandler(
    filename = function() {paste(heatName(),'-FPKM','.csv', sep='')},
    content=function(file){
      s <- make_heat()
      sub1 <- genefinder_data1[as.numeric(s),]
      write.csv(sub1, file,row.names = F)
    }
  )
  
  output$get_all_data <- downloadHandler(
    filename = function() {paste(paste(as.character(sapply(strsplit(gene_search(),",")[[1]],function(i){gsub("^\\s+|\\s+$", "", i)})),collapse = "-"),'-Heatmap-All-Data.csv',sep = '')},
    content=function(file){
      s <- heat_sub_data()
      sub1 <- genefinder_data1[as.numeric(s$index),]
      write.csv(sub1, file,row.names = F)
    }
  )
  
  heatName <- reactive({
    s <- make_heat()
    sub <- pcascores[as.numeric(s),]
    return(sub$gene)
  })
  
  Plot2 <- reactive({
    s <- make_heat()
    sub <- pcascores[as.numeric(s),]
    if(sub$node == 0){
      pca_plot <- pca_plot+
        ggtitle(paste(sub$gene,"Expression Pattern"))+
        scale_colour_manual(values=c("#e7298a","#7570b3","#d95f02","#1b9e77"))
    }else{
      pca_plot <- pca_plot +
        geom_point(data=sub,aes(PC1,PC2,size=4))+
        ggtitle(paste(sub$gene,"Expression Pattern"))+
        scale_colour_manual(values=c("#e7298a","#7570b3","#d95f02","#1b9e77"))
      if(!is.null(input$arrow)){
        sub$PC1_ <- sub$PC1+0.6
        sub$PC2_ <- sub$PC2+0.6
        pca_plot <- pca_plot+
          geom_curve(data=sub,aes(xend=PC1+0.1,x=PC1_,yend=PC2+0.1,y=PC2_),color="black",arrow = arrow(length = unit(0.03, "npc")),curvature = 0)
      }
    }
    pca_plot
  })
  
  output$is.zero <- renderUI({
    s <- make_heat()
    sub <- pcascores[as.numeric(s),]
    
    if(sub$node == 0){
      tags$div(class = "warning",p("WARNING: This gene was not able to be classified and no point is created"))
    }else{
      checkboxGroupInput("arrow","Include Arrow",c("Yes"="yes"))
    }
  })
  
  output$plot2 <- renderPlot({
    Plot2()
  })
  
  output$downloadDiversity <- downloadHandler(
    filename = function() {paste(heatName(),'-PCA', '.png', sep='')},
    content=function(file){
      png(file,height = 350,width = 350)
      print(Plot2())
      dev.off()
    },
    contentType='image/png'
  )
  
  gene_search <- eventReactive(input$gene_search_button,{
    input$gene_search_text 
  })
  
  heat_sub_data <- eventReactive(input$gene_search_button,{
    s <- gene_search()
    s <- sapply(strsplit(s,",")[[1]],function(i){gsub("^\\s+|\\s+$", "", i)})
    my_logic <- data.frame(sapply(s,function(i){
      if(i == ""){
        return()
      }else if(str_detect(i,"GO:")==T){
        str_detect(annot$go,i)
      }else if(str_detect(i,"Manes.")==T){
        str_detect(annot$gene,i)
      }else{
        str_detect(tolower(annot$annot),tolower(i))
      }
    }))
    if(input$heat_radio == "AND"){
      sub <- annot[rowSums(my_logic)==length(s),]
    }else{
      sub <- annot[rowSums(my_logic) > 0,]
    }
    return(sub)
  })
  
  output$gene_search_table <- renderDataTable({
    sub <- heat_sub_data()
    test <- gsub(":","%3A",as.character(sub[,3]))
    test <- gsub("-","..",test)
    test <- sapply(test,function(i) paste("<a target='_blank' href='https://phytozome.jgi.doe.gov/jbrowse/index.html?data=genomes%2FMesculenta%2F&loc=",i,"&tracks=Transcripts%2CAlt_Transcripts%2CPASA_assembly%2CBlastx_protein%2CBlatx_Plant_protein&highlight='>","View in Phytozome","</a>",sep = ""))
    sub[,3] <- test
    
    datatable(sub,selection = 'single',rownames = F,escape = F,options = list(sDom  = '<"top">lrt<"bottom">ip')) %>% 
      formatStyle('status',  backgroundColor = styleEqual(c("FALSE","TRUE"),c("green","orange")),color = styleEqual(c("FALSE","TRUE"),c("green","orange")),fontWeight = 'bold')
    
  },server = T)
  
  
  
  
  #********************************************************************************************************************
  # for gene finder
  #********************************************************************************************************************
  output$make_goButton <- renderUI({
    actionButton('goButton',"Update Table", icon("refresh"))
  })
  
  output$table_out <- renderUI({
    tissue_choices()
    fluidRow(
      box(width=9,style = "overflow-y:scroll",
        solidHeader = TRUE,
        title='Annotations of Filtered Genes',
        status = 'info',
        dataTableOutput('table2'),
        br(),
        br(),
        hr(),
        downloadButton('x3', "Download All Data for Filtered Genes (CSV)"),
        br(),
        br(),
        p("To generate the heatmap for a particular gene found using the finder, copy and paste the 
          Manes name into the search bar in the heatmap page")
        # ,
        # actionButton('transfer','Transfer Selection to Heatmap')
        )
    )
  })
  
  tissue_choices <- eventReactive(input$goButton,{
    rbind(c("FEC",input$fec_sel),c("Fibrous.Root",input$fr_sel),c("Lateral.Bud",input$lb_sel),c("Leaf",input$leaf_sel),
      c("Mid.Vein",input$mv_sel),c("OES",input$oes_sel),c("Petiole",input$petiole_sel),c("RAM",input$ram_sel),
      c("SAM",input$sam_sel),c("Stem",input$stem_sel),c("Storage.Root",input$sr_sel))
  })
  
  gene_rows <- eventReactive(tissue_choices(),{
    sel <- tissue_choices()
    on_val <- input$fpkmOn
    
    if(suppressWarnings(is.na(as.numeric(on_val)))){
      on_val <- 10
    }else{
      on_val <- as.numeric(on_val)
    }
    
    off_val <- input$fpkmOff
    
    if(suppressWarnings(is.na(as.numeric(off_val)))){
      off_val <- 1
    }else{
      off_val <- as.numeric(off_val)
    }
    
    on <- sel[sel[,2]=='on',1]
    off <- sel[sel[,2]=='off',1]
    
    if(length(on) == 0){
      on_mat <- matrix(1,ncol=1,nrow=nrow(genefinder_data))
    }else{
      gsub_on <- genefinder_data[,rowSums(sapply(on,function(i) str_detect(names(genefinder_data),i)))>=1]
      on_mat <-  rowSums(gsub_on >on_val) == ncol(gsub_on)
    }
    on_final <- on_mat
    
    
    if(length(off)==0){
      off_mat <- matrix(1,ncol=1,nrow=nrow(genefinder_data))
    }else{
      gsub_off <- genefinder_data[,rowSums(sapply(off,function(i) str_detect(names(genefinder_data),i)))>=1]
      off_mat <- rowSums(gsub_off < off_val) == ncol(gsub_off)
    }
    off_final <- off_mat
    
    final <- (off_final+on_final)==2
    #annot.out <- annot[final,]
    return(final)
  })
  
  output$table2 = renderDataTable({
    indicies <- gene_rows()
    annot.out <- annot[indicies,][,-1]
    test <- gsub(":","%3A",as.character(annot.out[,2]))
    test <- gsub("-","..",test)
    test <- sapply(test,function(i) paste("<a target='_blank' href='https://phytozome.jgi.doe.gov/jbrowse/index.html?data=genomes%2FMesculenta%2F&loc=",i,"&tracks=Transcripts%2CAlt_Transcripts%2CPASA_assembly%2CBlastx_protein%2CBlatx_Plant_protein&highlight='>","View in Phytozome","</a>",sep = ""))
    annot.out[,2] <- test
    datatable(annot.out,selection = 'none',rownames = F,escape = F)%>% 
      formatStyle('status',  backgroundColor = styleEqual(c("FALSE","TRUE"),c("green","orange")),color = styleEqual(c("FALSE","TRUE"),c("green","orange")),fontWeight = 'bold')
  },server = T)
  
  
  genefinderOut <- eventReactive(tissue_choices(),{
    sel <- tissue_choices()
    if(sum(sel[,2]=='on')==0){
      on <- "None"
    }else{
      on <- sel[sel[,2]=='on',1]
    }
    if(sum(sel[,2]=='off')==0){
      off <- "None"
    }else{
      off <- sel[sel[,2]=='off',1]
    }
    return(paste('Genefinder--',"ON-",paste(on,collapse = ""),"_OFF-",paste(off,collapse = ""),".csv",sep=""))
  })
  
  output$x3 = downloadHandler(genefinderOut(), content = function(file) {
    sel <- tissue_choices()
    indicies <- gene_rows()
    gene_down <- genefinder_data1[indicies,]
    write.csv(gene_down, file,row.names = F)
  })
  
  #********************************************************************************************************************
  # for cas_xam
  #********************************************************************************************************************
  
  casxam_f <- function(cd, bg, des, v_treat, fc_cut){
    #sub <- cd[apply(cd, 1, function(i) all(c(i[5] %in% v_treat & i[6] %in% v_treat)) & (abs(as.numeric(i[10])) > fc_cut) & (as.numeric(i[13]) < 0.05)),]
    sub <- cd[with(cd, (sample_1 %in% v_treat & sample_2 %in% v_treat) & (abs(as.numeric(log2.fold_change.)) > fc_cut) & (as.numeric(q_value) < 0.05)),]
    sub$log.qvalue <- -log10(sub$q_value)
    sel_srr <- des$Run[des$treatment %in% v_treat]
    test <- bg[,sapply(colnames(bg),function(i){unlist(lapply(strsplit(i, "[.]"),function(j)j[2]))}) %in% sel_srr]
    colnames(test) <- unlist(lapply(strsplit(colnames(test), "[.]"),function(i) paste(i[2], des[des$Run == i[2],"treatment"],sep = ".")))
    test <- cbind(bg$gene_name, test)
    colnames(test)[1] <- "gene"
    cdbg1 <- join(sub, test, by = "gene", type = "inner")
    cdbg1 <- cdbg1[,-c(1)]
    return(cdbg1)
  }
  
  cdbg <- reactiveValues(data = NULL)
  casxam_annot <- reactiveValues(data = NULL)
  casxam_gene_search <- reactiveValues(data = NULL)
  annot_search <- reactiveValues(data = NULL)
  cdbg_annot <- reactiveValues(data = NULL)
  iterator <- reactiveValues(data = 0)
  
  observeEvent(input$casxam_sub_table_button, {
    print("test")
    id <- showNotification(h3("Subsetting data..."), duration = NULL)
    casxam_selections <- rbind(c("mock_8hr",input$mock8hr), c("mock_24hr",input$mock24hr), c("mock_50hr",input$mock50hr), 
                             c("Xam668_8hr",input$xam6688hr), c("Xam668_24hr",input$xam66824hr), c("Xam668_50hr",input$xam66850hr),
                             c("Xe_8hr",input$xe8hr),c("Xe_24hr",input$xe24hr), c("Xe_50",input$xe50hr), 
                             c("Xe(TAL20_Xam668)_8hr",input$xam668xe8hr), c("Xe(TAL20_Xam668)_24hr",input$xam668xe24hr),
                             c("Xe(TAL20_Xam668)_50hr",input$xam668xe50hr))
    v_treat <- casxam_selections[,1][which(casxam_selections[,2] == "include")]
    at1 <- at[,c("annot", "gene_name", "model", "gene")]
    cdbg$data <- casxam_f(cd, bg, des, v_treat, input$fc_cut)
    colnames(cdbg$data)[2] <- "gene_name"
    annot$data <- subset(at1, at1$gene_name %in% cdbg$data$gene_name)
    cdbg_annot$data <- join(cdbg$data, annot$data, by = "gene_name")
    cdbg_annot$data <- cdbg_annot$data[-c(6, 10, 11, 13, 14)]
    removeNotification(id)
  })
  
  observeEvent(input$casxam_gene_search_button, {
    casxam_gene_search$data <- subset(cd, cd$gene %in% input$casxam_gene_search_text)
    annot_search$data <- subset(at, at$gene_name %in% input$casxam_gene_search_text)
    iterator$data <- iterator$data + 1
  })
  
  observeEvent(input$tabs, {
    iterator$data <- iterator$data + 1
  })
  
  output$cdbg_sub <- renderUI({
    if(input$tabs == "treatment_panel"){
      if(!is.null(cdbg$data)){
        box(style = "overflow-y:scroll",width=12,title = "Subsetting by Treatment",solidHeader = T,status = 'success',collapsible = TRUE,
            br(),
            p("Select a gene from the data table to view a boxplot of the treatment comparison."),
            hr(),
            dataTableOutput("casxam_sub_table"),
            br(),
            downloadButton("casxam_data_download","Download Data (tsv)"),
            hr(),
            column(width = 6,
                   plotOutput("cdbg_boxplot"),
                   uiOutput("casxam_boxplot_download_ui")
            )
        )
      }
    }
    else if(input$tabs == "gene_panel"){  
      if(!is.null(casxam_gene_search$data)){
        box(style = "overflow-y:scroll",width=12,title = "Subsetting by Gene",solidHeader = T,status = 'success',collapsible = TRUE,
            br(),
            p("Select a gene from the data table to view a boxplot of the treatment comparison."),
            hr(),
            dataTableOutput("casxam_gene_sub_table"),
            br(),
            downloadButton("casxam_gene_data_download","Download Data (tsv)"),
            hr(),
            column(width = 6,
                   plotOutput("cdbg_boxplot"),
                   uiOutput("casxam_boxplot_download_ui")
            )
        )
      }
    }
  })
  
  output$casxam_gene_sub_table <- renderDataTable({
    datatable(annot_search$data,rownames = F,selection = "single", options = list(pageLength = 10, autoWidth = F))
  })
  
  output$casxam_sub_table <- renderDataTable({
    datatable(cdbg_annot$data,rownames = F,selection = "single", options = list(pageLength = 10, autoWidth = F))
  })
  
  casxam_gene_sel <- eventReactive(c((input$casxam_gene_sub_table_row_last_clicked + iterator$data), (input$casxam_sub_table_row_last_clicked)),{
    if(input$tabs == "treatment_panel"){  
      temp <- cdbg_annot$data[input$casxam_sub_table_row_last_clicked, "gene_name"]
    }else if(input$tabs == "gene_panel"){
      temp <- annot_search$data[input$casxam_gene_sub_table_row_last_clicked, "gene_name"]
    }
    box_dat <- bg[bg$gene_name == temp,]
    box_des <- des
    colnames(box_dat)[5:40] <- unlist(lapply(strsplit(colnames(box_dat)[5:40], "[.]"),function(i)i[2]))
    box_dat <- melt(box_dat, id.vars = colnames(box_dat[1:4]))
    colnames(box_dat)[5:6] <- c("sample","FPKM")
    box_des <- box_des[1:36,]
    sub <- box_dat[box_dat$gene_name == temp,]
    sub$meta <- sapply(sub$sample, function(i) box_des[box_des$Run == i,"treatment"])
    sub$treatment <- unlist(lapply(strsplit(sub$meta, "_"),function(i){paste(i[-length(i)],collapse = "")}))
    sub$time <- unlist(lapply(strsplit(sub$meta, "_"),function(i){i[length(i)]}))
    sub <- aggregate(data = sub, FPKM~sample + treatment + time + gene_name, FUN = "sum")
    sub$time_f <- factor(sub$time, levels = c("8hr","24hr","50hr"))
    sub
  })
  
  output$cdbg_boxplot <- renderPlot({
    if(!is.null(as.character(input$casxam_sub_table_row_last_clicked)) | !is.null(as.character(input$casxam_gene_sub_table_row_last_clicked))){
      if((input$tabs == "treatment_panel" & !is.null(input$casxam_sub_table_row_last_clicked) | (input$tabs == "gene_panel" & !is.null(input$casxam_gene_sub_table_row_last_clicked)))){ 
        ggplot(data = casxam_gene_sel(), aes(x = treatment, y = FPKM))+
          facet_grid(~time_f)+
          geom_boxplot()+
          theme_light()+
          xlab("")+
          ylab("FPKM")+
          ggtitle(casxam_gene_sel()$gene_name[1])+
          theme(axis.text = element_text(size = 12),
                axis.title= element_text(size = 18))+
          theme(plot.title = element_text(hjust = 0.5),
                strip.background=element_rect(fill="gray50"),
                strip.text.x=element_text(size=14,color="white"),
                strip.text.y=element_text(size=14,color="white"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
  })
  
  #sDom  = '<"top">lrt<"bottom">ip',
  output$casxam_data_download <- downloadHandler(
    filename = function() {"cas_xam_2014_table.tsv"},
    content = function(file){
      write.table(cdbg_annot$data,file,row.names = FALSE, quote = FALSE,sep = "\t")
    }
  )
  
  output$casxam_gene_data_download <- downloadHandler(
    filename = function() {"cas_xam_2014_table.tsv"},
    content = function(file){
      write.table(casxam_gene_search$data,file,row.names = FALSE, quote = FALSE,sep = "\t")
    }
  )
  
  casxam_boxplot <- reactive({
    ggplot(data = casxam_gene_sel(), aes(x = treatment, y = FPKM))+
      facet_grid(~time_f)+
      geom_boxplot()+
      theme_light()+
      xlab("")+
      ylab("FPKM")+
      ggtitle(casxam_gene_sel()$gene_name[1])+
      theme(axis.text = element_text(size = 12),
            axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
            strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$casxam_boxplot_download <- downloadHandler(
    filename = function() {"cas_xam_2014_boxplot.png"},
    content=function(file){
      ggsave(file,casxam_boxplot(),device = "png",width = 6,height = 4,dpi = 300)
    })
  
  output$casxam_boxplot_download_ui <- renderUI({
    if((input$tabs == "treatment_panel" & !is.null(input$casxam_sub_table_row_last_clicked)) | (input$tabs == "gene_panel" & !is.null(input$casxam_gene_sub_table_row_last_clicked))){
      downloadButton("casxam_boxplot_download", "Download Boxplot (png)")
    }
  })
}

shinyApp(ui, server)
