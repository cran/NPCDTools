# NPC & GNPC
library(gdata)
library(gmodels)
library(gplots)
library(gtools)
library(MASS)
library(mclust)
library(Matrix)  
library(shinythemes)
library(shiny)
library(NPCD)
library(GDINA)
library(CDM)

ui = navbarPage(
  tags$head(
    tags$style(HTML("
    .equal-height-box {
      background-color: #f5f5f5;
      border: 1px solid #ddd;
      border-radius: 4px;
      padding: 15px;
      height: 100%;
      min-height: 600px;
    }
    .box-title {
      font-weight: bold;
      font-size: 18px;
      margin-bottom: 10px;
      padding-bottom: 10px;
      border-bottom: 2px solid #ddd;
    }
  "))
  ),
  
  strong("Session 2: Nonparametric Classification Methods"), 
  theme = shinytheme("united"),
  tabPanel(
    strong("NPC & GNPC"),
    titlePanel(h2("Nonparametric Cognitive Diagnosis Workshop")),
    
    h3("Topic: Nonparametric Classification (NPC) and General NPC (GNPC) Methods"),
    br(),
    h4(p("1. This exercise is aimed to demonstrate how to use the NPC and the GNPC methods to classify data and then to
                                  compare the results with that obtained by fitting the data with a user-selected CDM.")),
    hr(),
    
    h4("1.1 Simulate data"),
    br(),
    
    h4(p("To demonstrate, let's use a Q-matrix containing 30 items requiring 5 attributes from the R package", code("GDINA"), "to generated data.")),
    hr(),
    
    h4(p("The function", code("simGDINA(N, Q, gs.parm = gs, model = 'GDINA',...)"), "is used to generate data.")),
    h4("Here, we need to insert the arguments 'N', 'gs.parm', and 'model'."),
    hr(),
    
fluidRow(
  # Block 1: Parameter & Model Selection
  column(3,
    div(class = "equal-height-box",
      div(class = "box-title", "PARAMETER & MODEL SELECTION"),
      
      h5(strong("Number of Examinees (N)")),
      fluidRow(
        column(7, sliderInput(inputId = "N",
                              label = NULL,
                              min = 1,
                              max = 500,
                              value = 100)),
        column(5, numericInput(inputId = "N_manual",
                               label = NULL,
                               value = 100,
                               min = 1,
                               max = 500))
      ),
      
      h5(strong("s, g ~ Unif(a,b)")),
      
      fluidRow(
        column(3, h6("a =")),
        column(5, sliderInput(inputId = "a",
                              label = NULL,
                              min = 0,
                              max = 0.5,
                              value = 0.1,
                              step = 0.01)),
        column(4, numericInput(inputId = "a_manual",
                               label = NULL,
                               value = 0.1,
                               min = 0,
                               max = 0.5,
                               step = 0.01))
      ),
      
      fluidRow(
        column(3, h6("b =")),
        column(5, sliderInput(inputId = "b",
                              label = NULL,
                              min = 0,
                              max = 0.5,
                              value = 0.2,
                              step = 0.01)),
        column(4, numericInput(inputId = "b_manual",
                               label = NULL,
                               value = 0.2,
                               min = 0,
                               max = 0.5,
                               step = 0.01))
      ),
      
      h5(strong("Select a Generating Model")),
      radioButtons("model", label = NULL,
                   choices = c("DINA" = "DINA", 
                              "DINO" = "DINO", 
                              "ACDM" = "ACDM", 
                              "RRUM" = "RRUM",
                              "G-DINA" = "GDINA"))
    )
  ),
  
  # Block 2: The Q-matrix
  column(3,
    div(class = "equal-height-box",
      div(class = "box-title", textOutput("caption.Q")),
      verbatimTextOutput("Q")
    )
  ),
  
  # Block 3: Data Generation
  column(6,
    div(class = "equal-height-box",
      div(class = "box-title", "DATA GENERATION"),
      
      actionButton("goButton1", "Generate data", class = "btn-primary btn-lg"),
      h5(textOutput("generateMessage"), style = "color:blue; margin-top: 10px;"), 
      br(),
      
      actionButton("goButton2", "Sample data", class = "btn-info btn-lg"), 
      br(),
      br(),
      verbatimTextOutput("sampleData")
    )
  )
),

br(),
hr(),
    
    #======= GDINA =====================
    h4("1.2 Fit the simulated data with a user-selected CDM and compute the PAR and AAR."),
    br(),
    
    h4(p("The function", code("GDINA(dat, Q, model = 'GDINA')"), "is used to fit the data with a selected model.")),
    h4(p(code("personparm()"), "can be used to extract the estimated attribute patterns.")),
    hr(),
    
fluidRow(
  # Block 1: Select a Fitted Model
  column(3,
         div(class = "equal-height-box",
             div(class = "box-title", "SELECT A FITTED MODEL"),
             
             radioButtons("fitted.model", 
                          label = NULL,
                          choices = c("DINA" = "DINA", 
                                      "G-DINA" = "GDINA"))
         )
  ),
  
  # Block 2: Model Fitting
  column(5,
         div(class = "equal-height-box",
             div(class = "box-title", "MODEL FITTING"),
             
             actionButton("goButton3", "Fit the data with the selected CDM", 
                          class = "btn-primary btn-lg"),
             h5(textOutput("generateMessage2"), style = "color:blue; margin-top: 10px;"), 
             br(),
             br(),
             
             actionButton("goButton4", "Sample estimated attribute profiles", 
                          class = "btn-info btn-lg"),
             br(),
             
             tags$style("#sample.att.gdina {color: blue;}"),
             verbatimTextOutput("sample.att.gdina")
         )
  ),
  
  # Block 3: Compute PAR and AAR
  column(3,
         div(class = "equal-height-box",
             div(class = "box-title", "COMPUTE PAR AND AAR"),
             
             actionButton("goButton5", "Compute PAR and AAR", class = "btn-success btn-lg"),
             br(),
             
             tags$style("#rate.gdina {color: blue;}"),
             verbatimTextOutput("rate.gdina")
         )
  )
),
br(),
hr(),
    
    #=== NPC =================
    h4("1.3 Classify the simulated data using the NPC method"),
    br(),
    
    h4(p("The", code("R"), "function in the package", code("NPCD"), "for running the NPC method is called", code("NPC"), ", formulated as")),
    
    h4(code("NPC = function (Y, Q, gate = c('AND', 'OR'), method = c('Hamming', 'Weighted', 'Penalized'), wg = 1, ws = 1)")),
    br(),
    
    h4(p(code("Y"),"= data")),
    h4(p(code("Q"),"= Q-matrix")),
    h4(p(code("gate"),"= conjunctive or disjunctoive structure, respectively")),
    h4(p(code("method"),"= distance measure")),
    h4(p(code("wg & ws"),"= weights assigned to the guessing and slipping parameters if the penalized Hamming distance is used")),
    hr(),
    br(),
    
fluidRow(
  # Block 1: Classify Using NPC
  column(3,
         div(class = "equal-height-box",
             div(class = "box-title", "CLASSIFY USING NPC"),
             
             h5(strong("Select a cognitive structure")),
             selectInput("gate", 
                         label = NULL,
                         choices = c("AND" = "AND", 
                                     "OR" = "OR")),
             
             h5(strong("Select a distance measure")),
             selectInput("method", 
                         label = NULL,
                         choices = c("Hamming distance" = "Hamming", 
                                     "Weighted Hamming distance" = "Weighted",
                                     "Penalized Hamming distance" = "Penalized"))
         )
  ),
  
  # Block 2: Model Fitting
  column(5,
         div(class = "equal-height-box",
             div(class = "box-title", "MODEL FITTING"),
             
             actionButton("goButton7", "Run the NPC method", 
                          class = "btn-primary btn-lg"),
             h5(textOutput("generateMessage3"), style = "color:blue; margin-top: 10px;"), 
             br(),
             br(),
             
             actionButton("goButton8", "Sample estimated attribute profiles", 
                          class = "btn-info btn-lg"),
             br(),
             br(),
             
             tags$style("#sample.npc.Att {color: blue;}"),
             verbatimTextOutput("sample.npc.Att")
         )
  ),
  
  # Block 3: Compute PAR and AAR
  column(3,
         div(class = "equal-height-box",
             div(class = "box-title", "COMPUTE PAR AND AAR"),
             
             actionButton("goButton9", "Compute the PAR", class = "btn-success btn-lg"),
             br(),
             h5(textOutput("PAR")),
             br(),
             
             actionButton("goButton10", "Compute the AAR", class = "btn-success btn-lg"),
             br(),
             h5(textOutput("AAR"))
         )
  )
),
br(),
hr(),

    #===== GNPC ================================
    h4("1.4 Classify the simulated data using the GNPC method"),
    br(),
    
    h4(p("The", code("R"), "function for executing the GNPC method is called", code("GNPC"), ", formulated as")),
    
    h4(p(code("GNPC = function(Y, Q, initial.dis= c('hamming', 'whamming'), initial.gate = c('AND', 'OR', 'Mix'))"))),
    br(),
    
    h4(p("Let's analyze the same simulated data set using the GNPC method.")),
    hr(),
    
fluidRow(
  # Block 1: Classify Using GNPC
  column(3,
         div(class = "equal-height-box",
             div(class = "box-title", "CLASSIFY USING GNPC"),
             
             h5(strong("Select a distance measure for the initial")),
             selectInput("distance", 
                         label = NULL,
                         choices = c("Hamming distance" = "Hamming", 
                                     "Weighted Hamming distance" = "Whamming")),
             
             h5(strong("Select a cognitive structure for the initial")),
             selectInput("start", 
                         label = NULL,
                         choices = c("Conjunctive" = "conjunctive", 
                                     "Disjunctive" = "disjunctive",
                                     "Mixed Structure" = "mix"))
         )
  ),
  
  # Block 2: Model Fitting
  column(5,
         div(class = "equal-height-box",
             div(class = "box-title", "MODEL FITTING"),
             
             actionButton("goButton11", "Run the GNPC method", 
                          class = "btn-primary btn-lg"),
             h5(textOutput("generateMessage4"), style = "color:blue; margin-top: 10px;"), 
             br(),
             br(),
             
             actionButton("goButton12", "Obtain the weights", 
                          class = "btn-info btn-lg"),
             br(),
             
             verbatimTextOutput("weight"),
             br(),
             
             actionButton("goButton13", "Sample estimated attribute profiles", 
                          class = "btn-info btn-lg"),
             br(),
             
             tags$style("#sampleEstAtt.gnpc {color: blue;}"),
             verbatimTextOutput("sampleEstAtt.gnpc")
         )
  ),
  
  # Block 3: Compute PAR and AAR
  column(3,
         div(class = "equal-height-box",
             div(class = "box-title", "COMPUTE PAR AND AAR"),
             
             actionButton("goButton14", "Compute the PAR", class = "btn-success btn-lg"),
             br(),
             
             tags$style("#PAR.gnpc {color: blue;}"),
             h5(textOutput("PAR.gnpc")),
             br(),
             
             actionButton("goButton15", "Compute the AAR", class = "btn-success btn-lg"),
             br(),
             
             tags$style("#AAR.gnpc {color: blue;}"),
             h5(textOutput("AAR.gnpc"))
         )
  )
),
br(),
hr()
  ),
  
  ##=== Simulation with replications ====
  tabPanel(
    strong("Simulation with Replications"), 
    
    titlePanel(h2("Nonparametric Cognitive Diagnosis Workshop")),
    h3("Topic: Nonparametric Classification (NPC) and General NPC (GNPC) Methods"),
    br(),
    hr(),
    
    h4("Simulations with Replications"),
    br(),
    hr(),
    
    fluidRow(
      column(3,
             #sidebarPanel(
               h4(strong("DATA GENERATION")),
               hr(),
               
               sliderInput(inputId = "N",
                           label = h4(strong("Number of Examinees (N)")),
                           min = 1,
                           max = 500,
                           value = 100),
               
               h5(strong("s, g ~ Unif(a,b)")),
               
               sliderInput(inputId = "a",
                           label = "a =",
                           min = 0,
                           max = 0.5,
                           value = 0.1),
               
               sliderInput(inputId = "b",
                           label = "b =",
                           min = 0,
                           max = 0.5,
                           value = 0.2),
               
             radioButtons("gen.model", h4(strong("Select a Generating Model")), choices=(c("DINA"= "DINA", 
                                                                                               "DINO"="DINO", 
                                                                                               "ACDM"="ACDM", 
                                                                                               "RRUM"="RRUM",
                                                                                               "G-DINA"="GDINA"))),
             
             radioButtons("structure", h4(strong("Select an attribute structure")), choices=(c("Uniform"= "Uniform",
                                                                                           "Higher-Order" = "Higher-Order"))),
        
             radioButtons("R", h4(strong("# Replications")), choices=(c("10"= 10,
                                                                                "20"= 20,
                                                                                "50"= 50))),

              # width = 3
             #),
      ),
      column(3,
             h4(strong("FIT WITH THE G-DINA MODEL")),
             hr(),
             
             radioButtons("fitted.structure", h4(strong("Select a fitted attribute structure")), choices=(c("Uniform"= "Uniform",
                                                                                                    "Higher-Order" = "Higher-Order"))),
             hr(),
             br(),
             
             h4(strong("CLASSIFY USING NPC")),
             hr(),
             
             selectInput("gate", h4(strong("Select a relation")), choices=c("AND"= "AND", 
                                                                                       "OR"="OR")),
             br(),
             
             selectInput("method", h4(strong("Select a distance measure")), choices=c("Hamming distance"= "Hamming", 
                                                                                      "Weighted Hamming distance"="Weighted",
                                                                                      "Penalized Hamming distance" = "Penalized")),
             br(),
             hr(),
             
             h4(strong("CLASSIFY USING GNPC")),
             hr(),
             
             selectInput("start", h4(strong("Select a cognitive structure for the initial")), choices=c("Conjunctive"= "conjunctive", 
                                                                                                        "Disjunctive"="disjunctive",
                                                                                                        "Mix" = "mix")),
             br(),
             
             selectInput("distance", h4(strong("Select a distance measure for the initial")), choices=c("Hamming distance"= "Hamming", 
                                                                                                        "Weighted Hamming distance"="Whamming")),
             br()
      ),
      column(4,
             actionButton("rep.goButton1", h4("Start the Simulation")),
             br(),
             
             h4(div(verbatimTextOutput("rep.sim"))),
             br()
      )
    ),
    br(),
    br()
    
    
  ),
  
  
  # real data analysis
  tabPanel(
    strong("Real Data Analysis"), 
    
    titlePanel(h2("Nonparametric Cognitive Diagnosis Workshop")),
    h3("Topic: Nonparametric Classification (NPC) and General NPC (GNPC) Methods"),
    br(),
    hr(),
    
    h4("1. Analyze the ECPE Data"),
    br(),
    
    h4(p("Data from a retired version of the Examination for the Certificate of Proficiency in English (ECPE) were previously analyzed 
          with the LCDM by Templin and Hoffman (2013). They used MMLE-EM relying on the implementation of the EM algorithm in Mplus. 
          To further probe these discrepant findings, the estimated attribute profiles of a subset of 5 examinees were inspected. 
          These five examinees had already been used by Templin and Hoffman as exemplary cases for a deeper analysis of their findings.")),
    br(),
    
    h4(p("The ECPE data set has 2922 examinees responding to 28 items. Each item requires at most 3 attributes. 
          In this session, we will fit the ECPE data with the G-DINA model and the GNPC method and take a close look at the estimates of the
          same 5 examinees as reported in Templin and Hoffman (2013).")),
    hr(),
    
    h4(p("In particular, the interest is in the inconsistent estimates of examinees' attribute patterns obtained by fitting the data
          with the parametric method.")),
    hr(),
    
    fluidRow(
      column(3,
             h5(textOutput("caption.ELI.Q")),
             h5(div(verbatimTextOutput("ELI.Q")))
      ),
      column(6,
             actionButton("goButton.realdata1", h4("Sample ECPE data - Examinees 1, 10, 14, 29, 33")),
             br(),
             
             h5(div(verbatimTextOutput("sample.ELI.data"))),
             br()
             )
    ),
    hr(),
    br(),
    #actionButton("goButton.realdata2", h4("Fit ECPE with the G-DINA model")),
    #br(),
    #br(),
    
    #actionButton("goButton.realdata3", h4("Analyze ECPE with the GNPC method")),
    #br(),
    #br()
    
    fluidRow(
      column(3,
             actionButton("goButton.realdata2", h4("Obtain their propotional correct")),
             br(),
             
             h5(div(verbatimTextOutput("prop.correct"),style = "color:blue")),
             br(),
             br()
      ),
      column(3,
             actionButton("goButton.realdata3", h4("Obtain their G-DINA estimates")),
             br(),
             
             h5(div(verbatimTextOutput("sample.gdina"),style = "color:blue")),
             br(),
             br()
      ),
      column(3,
             actionButton("goButton.realdata4", h4("Obtain their GNPC estimates")),
             br(),
             
             h5(div(verbatimTextOutput("sample.gnpc"),style = "color:blue")),
             br(),
             br()
      )
    ),

 
    # Randomly select 5 items
    h4("2. Randomly Select 5 Examinees from the ECPE Data Set"),
    h4(p("In addition to the analysis with the 5 examinees selected by Templin and Hoffman (2013), it is also of our interest in learning
                                  how the other estimates behave and whether the inconsistencies also occur. 
                                  This part of analysis allows us to randomly select 5 examinees and look into their attribute profile estimates.")),
    hr(),
    actionButton("goButton.realdata7", h4("Randomly select 5 examinees")),
    br(),
    
    h5(div(textOutput("random.select"),style = "color:blue")),
    hr(),
    
    actionButton("goButton.realdata8", h4("Their responses")),
    br(),
    
    h5(div(verbatimTextOutput("rsample.ELI.data"),style = "color:blue")),
    hr(),
    
    fluidRow(
      column(3,
             actionButton("goButton.realdata9", h4("Obtain their propotional correct")),
             br(),
             
             h5(div(verbatimTextOutput("rprop.correct"),style = "color:blue")),
             br()),
      
      column(3,
             actionButton("goButton.realdata10", h4("Obtain their G-DINA estimates")),
             br(),
             
             h5(div(verbatimTextOutput("rsample.gdina"),style = "color:blue")),
             br(),
             br()),
      
      column(3,
             actionButton("goButton.realdata11", h4("Obtain their GNPC estimates")),
             br(),
             
             h5(div(verbatimTextOutput("rsample.gnpc"),style = "color:blue")),
             br(),
             br())
    ),
    br(),
    br()
  )
)
  
server = function(input, output, session){
  # Sync N slider and numeric input
  observeEvent(input$N_manual, {
    updateSliderInput(session, "N", value = input$N_manual)
  })
  
  observeEvent(input$N, {
    updateNumericInput(session, "N_manual", value = input$N)
  })
  
  # Sync a slider and numeric input
  observeEvent(input$a_manual, {
    updateSliderInput(session, "a", value = input$a_manual)
  })
  
  observeEvent(input$a, {
    updateNumericInput(session, "a_manual", value = input$a)
  })
  
  # Sync b slider and numeric input
  observeEvent(input$b_manual, {
    updateSliderInput(session, "b", value = input$b_manual)
  })
  
  observeEvent(input$b, {
    updateNumericInput(session, "b_manual", value = input$b)
  })
  
  Q = reactive({GDINA::sim30GDINA$simQ})
  
  output$caption.Q <- renderText("The Q-matrix")
  output$Q <- renderPrint(Q())
  
  sim <- isolate(eventReactive(input$goButton1, {
    K=dim(Q())[2]
    J=dim(Q())[1]
    gs <- data.frame(guess=runif(J, input$a, input$b),slip=runif(J, input$a, input$b))
    sim <- GDINA::simGDINA(input$N, Q(), gs.parm = gs, model = input$model, 
                    gs.args = list(type = "random", mono.constraint = TRUE), item.names = TRUE)
    return(sim)
  }))
  
  Y = isolate(reactive({
    extract(sim(),what = "dat")
  }))
  
  att = isolate(reactive({
    extract(sim(),what = "attribute")
  }))
  
  sampleData <- isolate(eventReactive(input$goButton2, {
    Y=Y()
    Y[1:5,]
  }))
  
  output$sampleData = renderPrint({
    sampleData()
  })
  
  N_used <- eventReactive(input$goButton1, {
    input$N_manual
  })
  
  output$generateMessage <- renderText({
    sim()  # This triggers the data generation
    paste("Data generated successfully! N =", N_used())
  })
  
  #==== G-DINA =======================
  estAtt.gdina <- eventReactive(input$goButton3, {
    mod <- GDINA::GDINA(dat = Y(), Q = Q(), model = input$fitted.model, verbose = 0)
    estatt = as.matrix(personparm(mod))
  })
  
  sample.att.gdina <- eventReactive(input$goButton4, {
    (estAtt.gdina())[1:15,]
  })
  
  output$sample.att.gdina = renderPrint({
    sample.att.gdina()
  })
  
  rate.gdina = eventReactive(input$goButton5, {
    rate = c(NPCDTools::PAR(estAtt.gdina(), att()), NPCDTools::AAR(estAtt.gdina(), att()))
    names(rate) = c("PAR", "AAR")
    rate
  })
  
  output$rate.gdina = renderPrint({
    rate.gdina() 
  })
  
  fitted_model_used <- eventReactive(input$goButton3, {
    input$fitted.model
  })
  
  output$generateMessage2 <- renderText({
    req(estAtt.gdina())  # Only show message if fitting succeeded
    paste("Model fitting completed successfully! Model =", fitted_model_used())
  })

  
  #aar.gdina = isolate(eventReactive(input$goButton6, {
  #  AAR(estAtt.gdina(), att())
  #}))
  
  #output$AAR.gdina = renderPrint({
  #  aar.gdina()
  #})
  
  
  #==== NPC ====================
  estAtt <- isolate(eventReactive(input$goButton7, {
    est=NPCDTools::NPC(Y(), Q(), input$gate, input$method)
    as.matrix(est$alpha.est)
  }))
  
  output$generateMessage3 <- renderText({
    req(estAtt())
    paste("NPC completed successfully!")
  })
  
  sampleEstAtt <- isolate(eventReactive(input$goButton8, {
    estatt=estAtt()
    estatt[1:15,]
  }))
  
  output$sample.npc.Att = renderPrint({
    sampleEstAtt()
  })
  
  par = isolate(eventReactive(input$goButton9, {
    NPCDTools::PAR(estAtt(), att())
  }))
  
  output$PAR = renderText({
    par()
  })
  
  aar = isolate(eventReactive(input$goButton10, {
    NPCDTools::AAR(estAtt(), att())
  }))
  
  output$AAR = renderText({
    aar()
  })
  
  #========= GNPC ==========================
  out.gnpc <- isolate(eventReactive(input$goButton11, {
    NPCDTools::GNPC(Y(), Q(), input$distance, input$start)
  }))
  
  weight = isolate(eventReactive(input$goButton12, {
    out=out.gnpc()
    as.matrix((out$weight)[1:10,])
  }))
  
  output$weight = renderPrint({
    round(weight(),3)
  })
  
  sampleEstAtt.gnpc <- isolate(eventReactive(input$goButton13, {
    out=out.gnpc()
    as.matrix((out$att.est)[1:5,])
  }))
  
  output$sampleEstAtt.gnpc = renderPrint({
    sampleEstAtt.gnpc()
  })
  
  output$generateMessage4 <- renderText({
    req(out.gnpc())
    paste("GNPC completed successfully!")
  })
  
  par.gnpc = isolate(eventReactive(input$goButton14, {
    NPCDTools::PAR((out.gnpc())$att.est, att())
  }))
  
  output$PAR.gnpc = renderText({
    par.gnpc()
  })
  
  
  aar.gnpc = isolate(eventReactive(input$goButton15, {
    NPCDTools::AAR((out.gnpc())$att.est, att())
  }))
  
  output$AAR.gnpc = renderText({
    aar.gnpc()
  })

  ##==== Simulation with replications
  rep.sim <- isolate(eventReactive(input$rep.goButton1, {
    K=dim(Q())[2]
    J=dim(Q())[1]
    gs <- data.frame(guess=runif(J, input$a, input$b),slip=runif(J, input$a, input$b))
    lambda <- data.frame(a=rep(1,K),b=seq(-2,2,length.out=K))
    rate.gdina=rate.npc = rate.gnpc = NULL
    for (r in 1:input$R){
      if (input$structure=="Uniform"){
        sim <- GDINA::simGDINA(input$N, Q(), gs.parm = gs, model = input$gen.model, 
                        gs.args = list(type = "random", mono.constraint = T), item.names = TRUE)
      }else if (input$structure=="Higher-Order"){
        theta <- rnorm(input$N)
        sim <- GDINA::simGDINA(input$N, Q(), gs.parm = gs, model = input$gen.model, 
                        gs.args = list(type = "random", mono.constraint = T), att.dist="higher.order",
                        higher.order.parm = list(theta = theta,lambda = lambda), item.names = TRUE)
      }
      Y = extract(sim,what = "dat")
      true.att = extract(sim,what = "attribute")
      #==== G-DINA ======
      if (input$fitted.structure=="Uniform"){
        est1 <- GDINA::GDINA(dat = Y, Q = Q(), model = "GDINA", att.dist="higher.order", higher.order=list(model = "2PL"), verbose = 0)
      }
      else if (input$fitted.structure=="Higher-Order"){
        est1 <- GDINA::GDINA(dat = Y, Q = Q(), model = "GDINA", verbose = 0)
      }
      
      #est1 <- GDINA(dat = Y, Q = Q(), model = "GDINA", att.dist="higher.order", higher.order=list(model = "2PL"), verbose = 0)
      estatt.gdina = as.matrix(personparm(est1))
      rate.gdina = rbind(rate.gdina, c(NPCDTools::PAR(estatt.gdina, true.att), NPCDTools::AAR(estatt.gdina, true.att)))
      #=== NPC ==========
      est2=NPCDTools::NPC(Y, Q(), input$gate, input$method)
      estatt.npc = as.matrix(est2$alpha.est)
      rate.npc = rbind(rate.npc, c(NPCDTools::PAR(estatt.npc, true.att), NPCDTools::AAR(estatt.npc, true.att)))
      #=== GNPC ===========
      est3 = NPCDTools::GNPC(Y, Q(), input$distance, input$start)
      estatt.gnpc = as.matrix(est3$att.est)
      rate.gnpc = rbind(rate.gnpc, c(NPCDTools::PAR(estatt.gnpc, true.att), NPCDTools::AAR(estatt.gnpc, true.att)))
    }
    out = rbind(colMeans(rate.gdina), colMeans(rate.npc), colMeans(rate.gnpc))
    rownames(out) = c("G-DINA", "NPC", "GNPC")
    colnames(out) = c("PAR", "AAR")
    out
  }))
  
  output$rep.sim = renderPrint({rep.sim()})
  
  ##==== real data analysis =================================
  #ELI.data=reactive({
    #data=read.table("C:/Users/cychiu/Desktop/Shiny/ELI_data.txt", row.names=paste("sub.", 1:2922, sep=""))
  #  data=read.table("ELI_data.txt", row.names=paste("sub.", 1:2922, sep=""))
  #  data=as.matrix(data)
  #})
  #ELI.Q = reactive({
    #Q=read.table("C:/Users/cychiu/Desktop/Shiny/ELI_Q.txt", col.names=c("A1","A2","A3"))
  #  Q=read.table("ELI_Q.txt", col.names=c("A1","A2","A3"), row.names=paste("item.", 1:28, sep=""))
    #Q = cbind(Q[1:14,], Q[15:nrow(Q),])
  #  Q=as.matrix(Q)
  #})
  
  ELI.data = reactive({
    data = CDM::data.ecpe$data[,-1]
    data = as.matrix(data)
    rownames(data) = paste0("sub.",1:2922)
    data
  })
  
  ELI.Q = reactive({
    Q = CDM::data.ecpe$q.matrix
    Q=as.matrix(Q)
    Q
  })
  
  output$ELI.Q = renderPrint(ELI.Q())
  
  sample.ELI <- eventReactive(input$goButton.realdata1, {
    sample=(ELI.data())[c(1,10,14,29,33),]
    sample
  })
  
  output$sample.ELI.data = renderPrint(sample.ELI())
  
  # G-DINA
  prop.correct = isolate(eventReactive(input$goButton.realdata2, {
    round(((ELI.data())[c(1,10,14,29,33),]%*%ELI.Q())/matrix(rep(colSums(ELI.Q()),5),5,3,byrow=TRUE),2)
  }))
  
  output$prop.correct = renderPrint({
    prop.correct()
  })
  
  estAtt.ELI.gdina <- isolate(eventReactive(input$goButton.realdata3, {
    mod <- GDINA::GDINA(dat = ELI.data(), Q = ELI.Q(), model = "GDINA", verbose = 0)
    att.gdina = as.matrix(personparm(mod))
    rownames(att.gdina) = paste0("sub.",1:2922)
    colnames(att.gdina) = c("Skill1","Skill2","Skill3")
    att.gdina
  }))
  
  sample.gdina <- reactive({estAtt.ELI.gdina()[c(1,10,14,29,33),]})
  
  #sample.gdina = isolate(eventReactive(input$goButton.realdata5, {
  #  (estAtt.ELI.gdina())[c(1,10,14,29,33),]
  #}))
  
  output$sample.gdina = renderPrint({
    sample.gdina()})
  
  #  GNPC
  estAtt.ELI.gnpc = isolate(eventReactive(input$goButton.realdata4, {
    out=NPCDTools::GNPC(ELI.data(), ELI.Q(), "Hamming", "conjunctive")
    colnames(out$att.est)=c("A1", "A2", "A3")
    att.gnpc = as.matrix(out$att.est)
    rownames(att.gnpc) = paste0("sub.",1:2922)
    colnames(att.gnpc) = c("Skill1","Skill2","Skill3")
    att.gnpc
  }))
  
  sample.gnpc = reactive({(estAtt.ELI.gnpc())[c(1,10,14,29,33),]})
  
  output$sample.gnpc = renderPrint({
    sample.gnpc()
  })
  
  # Randomly select 5 items
  
  random.index <- eventReactive(input$goButton.realdata7, {
    sample(1:2922,5,prob=rep(1/2922,2922))
  })
  
  random.index.text = reactive({paste(random.index(),sep="   ", collapse=", ")})
  
  output$random.select = renderText({random.index.text()})
  
  rsample.ELI <- eventReactive(input$goButton.realdata8, {
    (ELI.data())[random.index(),]
  })
  
  output$rsample.ELI.data = renderPrint({
    rsample.ELI()
  })
  
  
  rprop.correct = isolate(eventReactive(input$goButton.realdata9, {
    round(((ELI.data())[random.index(),]%*%ELI.Q())/matrix(rep(colSums(ELI.Q()),5),5,3,byrow=TRUE),2)
  }))
  
  output$rprop.correct = renderPrint({
    rprop.correct()
  })
  
  rsample.gdina = isolate(eventReactive(input$goButton.realdata10, {
    (estAtt.ELI.gdina())[random.index(),]
  }))
  
  output$rsample.gdina = renderPrint({
    rsample.gdina()
  })
  
  rsample.gnpc = isolate(eventReactive(input$goButton.realdata11, {
    (estAtt.ELI.gnpc())[random.index(),]
  }))
  
  output$rsample.gnpc = renderPrint({
    rsample.gnpc()
  })
  

}
  
  
  shinyApp(ui, server)
  
