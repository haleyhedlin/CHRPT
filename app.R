## web  app to take as input a woman's information 
## and to output her predicted risk of several outcomes
## based on a competing risk model built using data from WHI

## online at https://hedlin.shinyapps.io/shiny/

## author: Haley Hedlin

# used the following site https://www.shinyapps.io/admin/#/dashboard
# update app online with the code rsconnect::deployApp('/Users/hhedlin/shiny')



library(shiny)
library(cmprsk)
options(shiny.deprecation.messages=FALSE)

load('data/WHIfit.Rdata')  
load('data/modfitRobbins5.Rdata')
load('data/modfitRobbins15.Rdata')
load('data/trpred10hist.Rdata')
load('data/trpred5hist.Rdata')
load('data/trpred15hist.Rdata')
## load data for the "event ever models"
load('data/modfitC10.Rdata')  
load('data/modfitC5.Rdata')
load('data/modfitC15.Rdata')
load('data/pred10Chist.Rdata')
load('data/pred5Chist.Rdata')
load('data/pred15Chist.Rdata')



pkyrfun = function(smkage, qsmkage, smking, smkyrs, smknw, a, smkct){
  if(smkage == "Less than 15") SMOKAGEn = 1
  if(smkage == "15-19") SMOKAGEn = 2
  if(smkage == "20-24") SMOKAGEn = 3
  if(smkage == "25-29") SMOKAGEn = 4
  if(smkage == "30-34") SMOKAGEn = 5
  if(smkage == "35-39") SMOKAGEn = 6
  if(smkage == "40-44") SMOKAGEn = 7
  if(smkage == "45-49") SMOKAGEn = 8
  if(smkage == "50 or older") SMOKAGEn = 9

  if(smknw == "0"){
    if(qsmkage == "Less than 15") QSMOKAGEn = 1
    if(qsmkage == "15-19") QSMOKAGEn = 2
    if(qsmkage == "20-24") QSMOKAGEn = 3
    if(qsmkage == "25-29") QSMOKAGEn = 4
    if(qsmkage == "30-34") QSMOKAGEn = 5
    if(qsmkage == "35-39") QSMOKAGEn = 6
    if(qsmkage == "40-44") QSMOKAGEn = 7
    if(qsmkage == "45-49") QSMOKAGEn = 8
    if(qsmkage == "50-54") QSMOKAGEn = 9
    if(qsmkage == "55-59") QSMOKAGEn = 10
    if(qsmkage == "60 or older") QSMOKAGEn = 11   
  }
  
  if(smking == "1" & smkyrs == "50 or more years" & 
     a < 61){
    if(smknw == "1"){
      if(smkage == "Less than 15"){
        YRS = a - 12.5
      }
      if(smkage!="15-19" & smkage!="50 or older"){
        YRS = a - ((SMOKAGEn + 1)*5 + 2)
      }
      if(smkage =="50 or older"){ 
        YRS = a - 50
      }
    }
    if(smknw == "0"){
      if( (QSMOKAGEn + 1)*5 <= a & 
          a <= (QSMOKAGEn + 2)*5 - 1 ){
        if(SMOKAGEn == "Less than 15"){
          YRS = ((a - (QSMOKAGEn + 1)*5)/2 + 
                   (QSMOKAGEn + 1)*5) - 12.5
        }
        if(2 <= SMOKAGEn & SMOKAGEn <= 8){
          YRS = ((a - (QSMOKAGEn + 1)*5)/2 + 
                   (QSMOKAGEn + 1)*5) - ((SMOKAGEn + 1)*5 + 2)
        }
        if(smkage  =="50 or older"){
          YRS = (a - (QSMOKAGEn + 1)*5)/2 + 
            (QSMOKAGEn + 1)*5 - 50
        }
        if(a > (QSMOKAGEn + 2)*5 - 1){
          if(smkage =="Less than 15"){
            YRS = (QSMOKAGEn + 1)*5 - 12.5
          }
          if(2 <= SMOKAGEn & SMOKAGEn <= 8){ 
            YRS = (QSMOKAGEn + 1)*5 - ((SMOKAGEn + 1)*5 + 2)
          }
          if(SMOKAGEn == 9){
            YRS = (QSMOKAGEn + 1)*5 - 50
          }
        }
      }
    }
  }
  if(smking == "1" & smkyrs == "50 or more years" & a >= 61){
      if(SMOKAGEn == 1) YRSX = a - 12.5
      if(2 <= SMOKAGEn & SMOKAGEn <= 8){
        YRSX = a - ((SMOKAGEn + 1)*5 + 2)
      }
      if(SMOKAGEn == 9) YRSX = a - 50
      YRS = max(50,YRSX)
      rm(YRSX)
  }
  
  if(smking == "1" & is.na(smkyrs) ){
      if(smknw == "1"){
        if(SMOKAGEn == 1) YRS = a - 12.5
        if(2 <= SMOKAGEn & SMOKAGEn <= 8){
          YRS = a - ((SMOKAGEn + 1)*5 + 2)
        }
        if(SMOKAGEn == 9) YRS = a - 50
      }
      if(smknw == "0" & SMOKAGEn <= QSMOKAGEn){
        if((QSMOKAGEn + 1)*5 <= a &
           a <= (QSMOKAGEn + 2)*5 - 1){
        if(SMOKAGEn == QSMOKAGEn & !is.na(smkage)){
            YRS = 2.5}
        if(SMOKAGEn == 1){
            YRS = ((a - (QSMOKAGEn + 1)*5)/2 +
                     (QSMOKAGEn + 1)*5) - 12.5
          }
        if(2 <= SMOKAGEn & SMOKAGEn<= 8){
            YRS = ((a - (QSMOKAGEn + 1)*5)/2 +
                   (QSMOKAGEn + 1)*5) - ((SMOKAGEn + 1)*5 + 2)
          }
        if(SMOKAGEn == 9){
            YRS = ((a - (QSMOKAGEn + 1)*5)/2 +
                     (QSMOKAGEn + 1)*5) - 50
          }
      }

      if(a > (QSMOKAGEn + 2)*5 - 1){
          if(SMOKAGEn == QSMOKAGEn & !is.na(smkage)) YRS = 2.5
          if(SMOKAGEn == 1) YRS = (QSMOKAGEn + 1)*5 - 12.5
          if(2 <= SMOKAGEn & SMOKAGEn <= 8){
            YRS = (QSMOKAGEn + 1)*5 - ((SMOKAGEn + 1)*5 + 2)
          }
          if(SMOKAGEn == 9) YRS = (QSMOKAGEn + 1)*5 - 50
        }
    }
  } 
  
  if(smking == "1" & smkyrs!="50 or more years"){
      if(smkyrs=="Less than 5 years") YRS=2.5
      if(smkyrs=="5-9 years") YRS=7
      if(smkyrs=="10-19 years") YRS=15
      if(smkyrs=="20-29 years") YRS=25
      if(smkyrs=="30-39 years") YRS=35
      if(smkyrs=="40-49 years") YRS=45
  }
  
  if(smkct == "Less than 1") PKYRS = 0.5/20*YRS
  if(smkct == "1-4") PKYRS = 2.5/20*YRS
  if(smkct == "5-14") PKYRS = 10/20*YRS
  if(smkct == "15-24") PKYRS = 20/20*YRS
  if(smkct == "25-34") PKYRS = 30/20*YRS
  if(smkct == "35-44") PKYRS = 40/20*YRS
  if(smkct == "45 or more") PKYRS = 50/20*YRS
  
  return(PKYRS)
}  

hordy5 <- 5*365.25
hordy10 <- 10*365.25
hordy15 <- 15*365.25

raceoptions = list(0,1,2)
names(raceoptions) = c("White","Black", "Other")
ethnicoptions = list(1,0)
names(ethnicoptions) = c("Hispanic","Non-Hispanic")
YNoptions = list(0,1)
names(YNoptions) = c("No","Yes")
age55options = list(0,1)
names(age55options) = c("Younger than 55","55 years or older")
oophnumoptions = list(1,2)
names(oophnumoptions) = c("One removed", "Both removed")
options12 = list(1,2)
names(options12) = c(1,"2 or more")
measoptions = list(1,2)
names(measoptions) = c("inches and pounds?","centimeters and kilograms?")

# Define UI for application 
# will have a button for the person to choose which outcome they want preds for
# then have tabs for the 5, 10, 15 

ui <- fluidPage(
  fluidRow(
    column(width = 5,
           
  numericInput(inputId="AGE", label="Age", value=65,
                            min=45, max=120),
                radioButtons(inputId="race", label="Race",
                     choices=raceoptions,selected="0", inline=TRUE),
                radioButtons(inputId="ethnic", label="Ethnicity",
                             choices=ethnicoptions,selected="0", inline=TRUE),
                radioButtons(inputId="DIAB",
                             label="Did a doctor ever say that you had sugar diabetes or high blood sugar when you were not pregnant?",
                             choices=YNoptions, inline=TRUE),
                radioButtons(inputId="HICHOLRP",
                             label="Did a doctor ever say that you had high cholesterol requiring pills?",
                             choices=YNoptions, inline=TRUE),       
                radioButtons(inputId="MIGRAINE",
                             label="Did a doctor ever say that you had migraine headaches?",
                             choices=YNoptions, inline=TRUE), 
                radioButtons(inputId="ATRIALFB",
                             label="Did a doctor ever say that you had atrial fibrillation (a type of irregular heart beat)?",
                             choices=YNoptions, inline=TRUE), 
                radioButtons(inputId="thy",
                             label="Did a doctor ever say that you had a thyroid problem?",
                             choices=YNoptions, inline=TRUE),
                conditionalPanel(condition = "input.thy==1",
                  radioButtons(inputId="undthy",
                             label="Did a doctor ever say that you had underactive thyroid?",
                             choices=YNoptions, inline=TRUE)),
                conditionalPanel(condition = "input.thy==1",
                  radioButtons(inputId="ovrthy",
                             label="Did a doctor ever say that you had overactive thyroid?",
                             choices=YNoptions, inline=TRUE)),                 
        ### this is not ideal - would like the age of the broken bone to appear only once someone says they did not break their hip
        ### alternatively, could set it so the age question gets asked whenever someone says they broke a bone even though I may not use it
                radioButtons("bkbone", "Have you ever broken a bone?",
                            YNoptions, inline=TRUE),
                
                    conditionalPanel(condition = "input.bkbone==1",
                                 radioButtons("bkhip", "Did you break your hip?", 
                                 YNoptions, inline=TRUE)),
                    conditionalPanel(condition = "input.bkbone==1 & input.bkhip==0",
                                 radioButtons("bkage", "At what age did you first break a bone?", 
                                              age55options, inline=TRUE)),
                radioButtons("htn", "Did a doctor ever say that you had hypertension or high blood pressure?",
                     YNoptions, inline=TRUE),
                  conditionalPanel(condition = "input.htn==1",
                         radioButtons("htncurr", "Do you currently take medication for high blood pressure?", 
                                      YNoptions, inline=TRUE)),
        selectInput("menarche","How old were you when you had your first menstrual period (menses)?",
                    choices=c("<11",11:14,">14"),selected=13),
        radioButtons("ooph01", "Did you ever have an operation to have one or both of your ovaries taken out?",
                     YNoptions, inline=TRUE),
          conditionalPanel(condition = "input.ooph01==1",
                         radioButtons("oophnum", "Was one or both ovaries removed?", 
                                      oophnumoptions, inline=TRUE)),
        radioButtons("BRSTBIOP", "Have you ever had a breast biopsy (where a doctor removes part or all of  breast lump to check for cancer)?",
                     YNoptions, inline=TRUE),
        radioButtons("preg","Have you ever been pregnant?",
                    YNoptions, inline=TRUE),
          conditionalPanel(condition = "input.preg==1",
                         selectInput("parity", "How many of the pregnancies lasted at least 6 months?", 
                                      choices=c(0:4,"5 or more"))),
          conditionalPanel(condition = "input.preg==1",
                         selectInput("agefbir", "How old were you at the end of the first of these pregnancies?", 
                                     choices=c("< 20 years", "20-29 years", "30 or older"))),
          conditionalPanel(condition = "input.preg==1",
                  radioButtons("BRSTFEED", "Did you breastfeed or nurse any children for at least one month?",
                      YNoptions, inline=TRUE)),     
  radioButtons("LACTDIET","Are you currently on a lactose-free (no milk or dairy foods) diet?",
               YNoptions, inline=TRUE),          radioButtons("ALCOHOL","During your entire life, have you had at least 12 drinks of any kind of alcoholic beverage? (One drink of alcohol is about equal to one can of beer, one glass of wine, or one shot of liquor.)?",
                                                              YNoptions, inline=TRUE), 
  conditionalPanel(condition = "input.ALCOHOL==1",
                   radioButtons("ALCNOW","Do you drink alcoholic beverages now?",
                                YNoptions, selected=1, inline=TRUE)),        
  conditionalPanel(condition = "input.ALCOHOL==1 & input.ALCNOW==1",
                   selectInput("alcoholepi", "How often do you drink alcoholic beverages?", 
                               choices=c("Never or less than once per month", 
                                         "< 1 drink per week", "1-6 drinks per week",
                                         "> 6 drinks per week"))),
  radioButtons("SMOKING","During your entire life, have you have you smoked at least 100 cigarettes?",
               YNoptions, inline=TRUE),
  conditionalPanel(condition = "input.SMOKING==1",
                   selectInput("SMOKAGE", "How old were you when you first started smoking cigarettes regularly?", 
                               choices=c("Less than 15", "15-19","20-24",
                                         "25-29", "30-34", "35-39", "40-44",
                                         "45-49", "50 or older"))),
  conditionalPanel(condition = "input.SMOKING==1",
                   radioButtons("SMOKNOW","Do you smoke cigarettes now?",
                                YNoptions, selected=1, inline=TRUE)),        
  conditionalPanel(condition = "input.SMOKING==1 & input.SMOKNOW==0",
                   selectInput("QSMOKAGE", "How old were you when you quit smoking regularly?", 
                               choices=c("Less than 15", "15-19","20-24",
                                         "25-29", "30-34", "35-39", "40-44",
                                         "45-49", "50-54", "55-59", "60 or older"))),
  conditionalPanel(condition = "input.SMOKING==1",
                   selectInput("smokct","On average, how many cigarettes do you (did you) usually smoke each day?",
                               choices=c("Less than 1", "1-4","5-14",
                                         "15-24", "25-34","35-44", "45 or more"))),        
  conditionalPanel(condition = "input.SMOKING==1",
                   selectInput("smokyrs","How many years have you been (were you) a regular smoker? Do not count the times you stayed off cigarettes",
                               choices=c("Less than 5 years", "5-9 years",
                                         "10-19 years", "20-29 years", 
                                         "30-39 years", "40-49 years", "50 or more years"))),
  radioButtons("aspirin","Are you currently taking aspirin, also known as acetylsalicylic acid (ASA)?",
               YNoptions, inline=TRUE),        
  selectInput("GENHEL","In general, how would you rate your health?",
              choices=c("Fair or poor", "Good", 
                        "Very good","Excellent"),
              selected="Very good")
    ),
  column(width = 5, offset = 1,
        radioButtons("momalv","Is your natural mother still alive?",
                     YNoptions, selected=1, inline=TRUE),    
          conditionalPanel(condition = "input.momalv==0",
                         numericInput("momdiedage2", "How old was she when she died?", 
                                      value=90,min=14, max=120)),
        radioButtons("dadalv","Is your natural father still alive?",
                     YNoptions, selected=1, inline=TRUE),    
          conditionalPanel(condition = "input.dadalv==0",
                         numericInput("daddiedage2", "How old was he when he died?", 
                                      value=90,min=14, max=120)),
        radioButtons("MIREL","Did your mother or father or full-blooded sisters, full-blooded brothers, daughters, or sons ever have a heart attack or myocardial infarction?",
                     YNoptions, inline=TRUE),      
          conditionalPanel(condition = "input.MIREL==1",
                         radioButtons("midad", "Did your father have a heart attack or myocardial infarction?", 
                                      YNoptions, inline=TRUE)),
            conditionalPanel(condition = "input.midad==1",
                             selectInput("midadage", "How old was he when he had his first heart attack?", 
                                         choices=c("< 55 years", "55-64 years", "65 or older", "Not sure"),
                                         selected="Not sure")),
          conditionalPanel(condition = "input.MIREL==1",
                               radioButtons("mimom", "Did your mother have a heart attack or myocardial infarction?", 
                                            YNoptions, inline=TRUE)),
              conditionalPanel(condition = "input.mimom==1",
                               selectInput("mimomage", "How old was she when she had her first heart attack?", 
                                           choices=c("< 55 years", "55-64 years", "65 or older", "Not sure"), 
                                           selected="Not sure")),
        radioButtons("STRKREL","Did your mother, father, full-blooded sisters, full-blooded brothers, daughters, or sons ever have a stroke?",
                     YNoptions, inline=TRUE),   
          conditionalPanel(condition = "input.STRKREL==1",
                         radioButtons("strkreln2", "How many of these relatives had a stroke?", 
                                      options12, inline=TRUE)),        
        radioButtons("BKBONREL", "Did your mother or father ever break or fracture a bone after 40 years of age?",
                     YNoptions, inline=TRUE),
          conditionalPanel(condition = "input.BKBONREL==1",
                         radioButtons("bkhipmom", "Did your mother break her hip after age 40?", 
                                      YNoptions, inline=TRUE)),
          conditionalPanel(condition = "input.BKBONREL==1",
                         radioButtons("bkhipdad", "Did your father break his hip after age 40?", 
                                      YNoptions, inline=TRUE)),
        radioButtons("CANCFREL","Did any of your female relatives ever have cancer? (For female relatives, please answer about your mother, full-blooded sisters, daughters, and grandmothers. Do not include aunts, cousins, or nieces.)",
                     YNoptions, inline=TRUE),    
          conditionalPanel(condition = "input.CANCFREL==1",
            radioButtons("brcafrel2","Did your mother, full-blooded sisters, daughters, or grandmothers ever have breast cancer?",
                     YNoptions, inline=TRUE)),  
          conditionalPanel(condition = "input.CANCFREL==1",
            radioButtons("colofrel2","Did your mother, full-blooded sisters, or daughters ever have cancer of the colon, rectum, intestine, or bowel?",
                     YNoptions, inline=TRUE)),  
        radioButtons("CANCMREL","Did any of your male relatives ever have cancer? (For male relatives, please answer about your father, full-blooded brothers, and sons. Do not include uncles, cousins, or nephews.)",
                     YNoptions, inline=TRUE),            
          conditionalPanel(condition = "input.CANCMREL==1",
            radioButtons("colomrel2","Did your father, full-blooded brothers, or sons ever have cancer of the colon, rectum, intestine, or bowel?",
                     YNoptions, inline=TRUE)),    

        #### This is not exactly the same as WHI
        ## too extensive to ask all the questions WHI does
        selectInput("tepi","How often each week (7 days) do you usually exercise?",
                    choices=c("0 days per week", "1 day per week", 
                              paste(2:4,"days per week"),"5 or more days per week"),
                    selected="0 days per week"),        
        selectInput("tmin","How long do you usually exercise at one time?",
                    choices=c("Less than 20 minutes", "20-39 minutes", 
                              "40-59 minutes","1 hour or more"),
                    selected="Less than 20 minutes"), 
        radioButtons("meas","Do you prefer entering",
                     measoptions),        
        conditionalPanel(condition = "input.meas==1",
                         numericInput(inputId="HEIGHTin", label="What is your height in inches?", 
                                      value=60,min=25, max=96),
                         numericInput(inputId="WEIGHTlb", label="What is your weight in pounds?", 
                                      value=150,min=50, max=500), 
                        numericInput(inputId="WAISTin", label="What is your waist circumference in inches?", 
                                       value=60,min=10, max=100),   
                        numericInput(inputId="HIPin", label="What is your hip circumference in inches?", 
                                       value=60,min=10, max=100) ),         
        conditionalPanel(condition = "input.meas==2",
                         numericInput(inputId="HEIGHTcm", label="What is your height in centimeters?", 
                                      value=180,min=70, max=300),
                         numericInput(inputId="WEIGHTkg", label="What is your weight in kilograms?", 
                                      value=75,min=25, max=250), 
      #### do we need to provide instructions on calcuting hip and waist??
                         numericInput(inputId="WAISTcm", label="What is your waist circumference in centimeter?", 
                                      value=180,min=30, max=300),   
                         numericInput(inputId="HIPcm", label="What is your hip circumference in centimeter?", 
                                      value=180,min=30, max=300) ),         
        numericInput(inputId="PULSE60", label="What is your 60-second resting pulse? That is, how many times does your heart beat in a minute?", 
                     value=160,min=30, max=400),         
      
        numericInput(inputId="SYST", label="What is your systolic blood pressure in mmHg?", 
                     value=120,min=30, max=300)
        
  )
  ),
  
  fluidRow(             
                actionButton("MI","Myocardial Infarction"),
                actionButton("strk","Stroke"),
                actionButton("lc","Lung Cancer"),
                actionButton("bc","Breast Cancer"),
                actionButton("crc","Colorectal Cancer"),
                actionButton("hip","Hip Fracture"),
                actionButton("death","Death")
  ),
  
  br(),
  
  fluidRow(                
                tabsetPanel(
                  tabPanel("5 year risk", textOutput("summary5"), 
                           plotOutput("hist5"), textOutput("summary5C"), 
                           plotOutput("hist5C")), 
                  tabPanel("10 year risk", textOutput("summary10"), 
                           plotOutput("hist10"), textOutput("summary10C"), 
                           plotOutput("hist10C")), 
                  tabPanel("15 year risk", textOutput("summary15"), 
                           plotOutput("hist15"), textOutput("summary15C"), 
                           plotOutput("hist15C")),
                   tabPanel("Interpretation of predictions",
                            fluidRow(
                              column(10,
                                     htmlOutput("interp"), offset = 1))
                  ),
                  selected="10 year risk"
                )
                
                
  )
                )


# Define server logic required to draw a histogram
server <- function(input, output) {

  rv <- reactiveValues(evt="Myocardial Infarction", outc10 = MI10h, 
                       vars10 = c(MIvars, 90, 92), fit10=MIfit,
                       ind10 = findInterval(hordy10,MIfit$uftime), 
                       outc5 = MI5h, vars5 = c(MIvars, 95:96), 
                       fit5=MIfit5, ind5 = findInterval(hordy5,MIfit5$uftime), 
                       outc15 = MI15h, vars15 = c(MIvars, 93:92), fit15=MIfit15, 
                       ind15 = findInterval(hordy15,MIfit15$uftime),
                       outc5C = MI5Ch, vars5C = c(MIvars, 90,99),
                       fit5C = MIfit5C, 
                       ind5C = findInterval(hordy5,MIfit5C$uftime),
                       outc10C = MI10Ch,vars10C = c(MIvars, 90,99),
                       fit10C = MIfit10C,
                       ind10C = findInterval(hordy10,MIfit10C$uftime),
                       outc15C = MI15Ch,vars15C = c(MIvars, 92), 
                       fit15C = MIfit15C,
                       ind15C = findInterval(hordy15,MIfit15C$uftime))
  observeEvent(input$MI, {rv$evt = "Myocardial Infarction"
                            rv$outc5 <- MI5h
                            rv$vars5 = c(MIvars, 95:96)   
                            rv$fit5 = MIfit5
                            rv$ind5 <- findInterval(hordy5,MIfit5$uftime)
                            rv$outc10 <- MI10h
                            rv$vars10 = c(MIvars, 90, 92)  
                            rv$fit10 = MIfit
                            rv$ind10 <- findInterval(hordy10,MIfit$uftime)
                            rv$outc15 <- MI15h
                            rv$vars15 = c(MIvars, 93:92)   
                            rv$fit15 = MIfit15
                            rv$ind15 <- findInterval(hordy15,MIfit15$uftime)
               rv$outc5C <- MI5Ch
               rv$vars5C = c(MIvars, 90,99) 
               rv$fit5C = MIfit5C
               rv$ind5C <- findInterval(hordy5,MIfit5C$uftime)
               rv$outc10C <- MI10Ch
               rv$vars10C = c(MIvars, 90,99)  
               rv$fit10C = MIfit10C
               rv$ind10C <- findInterval(hordy10,MIfit10C$uftime)
               rv$outc15C <- MI15Ch
               rv$vars15C = c(MIvars, 92)   
               rv$fit15C = MIfit15C
               rv$ind15C <- findInterval(hordy15,MIfit15C$uftime)})
  observeEvent(input$strk, {rv$evt = "Stroke" 
                             rv$outc5 <- strk5h
                            rv$vars5 = c(strkvars,97)  
                            rv$fit5 = strkfit5
                            rv$ind5 <- findInterval(hordy5,strkfit5$uftime)
                            rv$outc10 <- strk10h
                            rv$vars10 = rv$vars15 = strkvars  
                            rv$fit10 = strkfit
                            rv$ind10 <- findInterval(hordy10,strkfit$uftime)
                            rv$outc15 <- strk15h
                            rv$fit15 = strkfit15
                            rv$ind15 <- findInterval(hordy15,strkfit15$uftime)
                     rv$outc5C <- strk5Ch
                            rv$vars5C = c(strkvars, 100)  
                            rv$fit5C = strkfit5C
                            rv$ind5C <- findInterval(hordy5,strkfit5C$uftime)
                            rv$outc10C <- strk10Ch
                            rv$vars10C = c(strkvars, 100)   
                            rv$fit10C = strkfit10C
                            rv$ind10C <- findInterval(hordy10,strkfit10C$uftime)
                            rv$outc15C <- strk15Ch
                            rv$vars15C = c(strkvars, 101)  
                            rv$fit15C = strkfit15C
                            rv$ind15C <- findInterval(hordy15,strkfit15C$uftime) })
  observeEvent(input$lc, {rv$evt = "Lung Cancer"
                            rv$vars5 = rv$vars10 = rv$vars15 = lcvars 
                            rv$outc5 <- lc5h
                            rv$fit5 = lcfit5
                            rv$ind5 <- findInterval(hordy5,lcfit5$uftime)
                            rv$outc10 <- lc10h
                            rv$fit10 = lcfit
                            rv$ind10 <- findInterval(hordy10,lcfit$uftime)
                            rv$outc15 <- lc15h
                            rv$fit15 = lcfit15
                            rv$ind15 <- findInterval(hordy15,lcfit15$uftime)
                      rv$outc5C <- lc5Ch
                            rv$vars5C = rv$vars10C = rv$vars15C = lcvars  
                            rv$fit5C = lcfit5C
                            rv$ind5C <- findInterval(hordy5,lcfit5C$uftime)
                            rv$outc10C <- lc10Ch  
                            rv$fit10C = lcfit10C
                            rv$ind10C <- findInterval(hordy10,lcfit10C$uftime)
                            rv$outc15C <- lc15Ch
                            rv$fit15C = lcfit15C
                            rv$ind15C <- findInterval(hordy15,lcfit15C$uftime)})                            
  observeEvent(input$bc, {rv$evt = "Breast Cancer" 
                            rv$vars5 = rv$vars10 = rv$vars15 = bcvars  
                            rv$outc5 <- bc5h
                            rv$fit5 = bcfit5
                            rv$ind5 <- findInterval(hordy5,bcfit5$uftime)
                            rv$outc10 <- bc10h
                            rv$fit10 = bcfit
                            rv$ind10 <- findInterval(hordy10,bcfit$uftime)
                            rv$outc15 <- bc15h
                            rv$fit15 = bcfit15
                            rv$ind15 <- findInterval(hordy15,bcfit15$uftime)
                        rv$outc5C <- bc5Ch
                            rv$vars5C = rv$vars10C = rv$vars15C = bcvars  
                            rv$fit5C = bcfit5C
                            rv$ind5C <- findInterval(hordy5,bcfit5C$uftime)
                            rv$outc10C <- bc10Ch
                            rv$fit10C = bcfit10C
                            rv$ind10C <- findInterval(hordy10,bcfit10C$uftime)
                            rv$outc15C <- bc15Ch   
                            rv$fit15C = bcfit15C
                            rv$ind15C <- findInterval(hordy15,bcfit15C$uftime)})                           
  observeEvent(input$crc, {rv$evt = "Colorectal Cancer"
                           rv$vars5 = rv$vars10 = rv$vars15 = crcvars
                            rv$outc5 <- crc5h
                            rv$fit5 = crcfit5
                            rv$ind5 <- findInterval(hordy5,crcfit5$uftime)
                            rv$outc10 <- crc10h
                            rv$fit10 = crcfit
                            rv$ind10 <- findInterval(hordy10,crcfit$uftime)
                            rv$outc15 <- crc15h
                            rv$fit15 = crcfit15
                            rv$ind15 <- findInterval(hordy15,crcfit15$uftime)
                       rv$outc5C <- crc5Ch
                            rv$vars5C = rv$vars10C = rv$vars15C = crcvars   
                            rv$fit5C = crcfit5C
                            rv$ind5C <- findInterval(hordy5,crcfit5C$uftime)
                            rv$outc10C <- crc10Ch 
                            rv$fit10C = crcfit10C
                            rv$ind10C <- findInterval(hordy10,crcfit10C$uftime)
                            rv$outc15C <- crc15Ch
                            rv$fit15C = crcfit15C
                            rv$ind15C <- findInterval(hordy15,crcfit15C$uftime)}) 
  observeEvent(input$hip, {rv$evt = "Hip Fracture"
                           rv$outc5 <- hip5h
                            rv$vars5 = c(hipvars,98) 
                            rv$fit5 = hipfit5
                            rv$ind5 <- findInterval(hordy5,hipfit5$uftime)
                            rv$outc10 <- hip10h
                            rv$vars10 = c(hipvars,91)   
                            rv$fit10 = hipfit
                            rv$ind10 <- findInterval(hordy10,hipfit$uftime)
                            rv$outc15 <- hip15h
                            rv$vars15 = c(hipvars,94) 
                            rv$fit15 = hipfit15
                            rv$ind15 <- findInterval(hordy15,hipfit15$uftime)
                        rv$outc5C <- hip5Ch
                            rv$vars5C = rv$vars10C = rv$vars15C = hipvars   
                            rv$fit5C = hipfit5C
                            rv$ind5C <- findInterval(hordy5,hipfit5C$uftime)
                            rv$outc10C <- hip10Ch
                            rv$fit10C = hipfit10C
                            rv$ind10C <- findInterval(hordy10,hipfit10C$uftime)
                            rv$outc15C <- hip15Ch   
                            rv$fit15C = hipfit15C
                            rv$ind15C <- findInterval(hordy15,hipfit15C$uftime)})   
  observeEvent(input$death, {rv$evt = "Death"
                             rv$outc5 <- death5h
                            rv$vars5 = c(deathvars,95:98)  
                            rv$fit5 = deathfit5
                            rv$ind5 <- findInterval(hordy5,deathfit5$uftime)
                            rv$outc10 <- death10h
                            rv$vars10 = c(deathvars,90:92)  
                            rv$fit10 = deathfit
                            rv$ind10 <- findInterval(hordy10,deathfit$uftime)
                            rv$outc15 <- death15h
                            rv$vars15 = c(deathvars, 92:94)
                            rv$fit15 = deathfit15
                            rv$ind15 <- findInterval(hordy15,deathfit15$uftime) 
  ## removed death for "true competing risk" model because the predict command is different
  ## using crc here just as a placeholder to avoid errors when calculating pred                         
                            rv$outc5C <- crc5Ch
                            rv$vars5C = rv$vars10C = rv$vars15C = crcvars   
                            rv$fit5C = crcfit5C
                            rv$ind5C <- findInterval(hordy5,crcfit5C$uftime)
                            rv$outc10C <- crc10Ch 
                            rv$fit10C = crcfit10C
                            rv$ind10C <- findInterval(hordy10,crcfit10C$uftime)
                            rv$outc15C <- crc15Ch
                            rv$fit15C = crcfit15C
                            rv$ind15C <- findInterval(hordy15,crcfit15C$uftime)}) 
  
  pred = reactive({
    # initialize a vector of covariates from what the woman entered
    vec = data.frame(matrix(NA,nrow=1,ncol=length(deathfit$coef)))
    names(vec) = names(deathfit$coef) 
    vec$AGE = input$AGE
    vec$ethnic2_3 = ifelse(input$ethnic=="0" & input$race=="1", 
                           1, 0)
    vec$ethnic2_8 = ifelse(input$ethnic=="1" | input$race=="2", 1, 0)
    vec$DIAB = as.numeric(input$DIAB)
    vec$HICHOLRP = as.numeric(input$HICHOLRP)
    vec$MIGRAINE = as.numeric(input$MIGRAINE)
    vec$ATRIALFB = as.numeric(input$ATRIALFB)
    vec$undthy = as.numeric(input$undthy)
    vec$ovrthy = as.numeric(input$ovrthy)
    vec$fracthip55_1 = as.numeric(input$bkbone=="1" & input$bkage=="0" &
                                    input$bkhip=="0")
    vec$fracthip55_2 = as.numeric(input$bkbone=="1" & input$bkage=="1" &
                                    input$bkhip=="0")
    vec$fracthip55_3 = as.numeric(input$bkbone=="1" & input$bkhip=="1")
    vec$HTNTRT_1 = as.numeric(input$htn=="1" & input$htncurr=="0")
    vec$HTNTRT_2 = as.numeric(input$htn=="1" & input$htncurr=="1")
    vec$menarche_2 = as.numeric(input$menarche=="<11")
    vec$menarche_3 = as.numeric(input$menarche=="11")
    vec$menarche_4 = as.numeric(input$menarche=="12")
    vec$menarche_6 = as.numeric(input$menarche=="14")
    vec$menarche_7 = as.numeric(input$menarche==">14")
    vec$BRSTFEED = as.numeric(input$BRSTFEED)
    vec$ooph_1 = as.numeric(input$ooph01=="1" & input$oophnum=="1")
    vec$ooph_2 = as.numeric(input$ooph01=="1" & input$oophnum=="2")
    vec$BRSTBIOP = as.numeric(input$BRSTBIOP)
    vec$PARITY_0 = as.numeric(input$preg=="1" & input$parity=="0")
    vec$PARITY_1 = as.numeric(input$preg=="1" & input$parity=="1")
    vec$PARITY_2 = as.numeric(input$preg=="1" & input$parity=="2")
    vec$PARITY_3 = as.numeric(input$preg=="1" & input$parity=="3")
    vec$PARITY_4 = as.numeric(input$preg=="1" & input$parity=="4")
    vec$PARITY_5 = as.numeric(input$preg=="1" & input$parity=="5 or more")
    vec$agefbir3_1  = as.numeric(input$preg=="1" & input$agefbir=="< 20 years")
    vec$agefbir3_2  = as.numeric(input$preg=="1" & input$agefbir=="20-29 years")
    vec$agefbir3_3  = as.numeric(input$preg=="1" & input$agefbir=="30 or older")
    vec$momdiedage2_1 = as.numeric(input$momalv=="0" & input$momdiedage2 < 40)
    vec$momdiedage2_2 = as.numeric(input$momalv=="0" & input$momdiedage2 >= 40 & 
                                     input$momdiedage2 < 50)
    vec$momdiedage2_3 = as.numeric(input$momalv=="0" & input$momdiedage2 >= 50 & 
                                     input$momdiedage2 < 60)
    vec$momdiedage2_4 = as.numeric(input$momalv=="0" & input$momdiedage2 >= 60 & 
                                     input$momdiedage2 < 70)
    vec$momdiedage2_5 = as.numeric(input$momalv=="0" & input$momdiedage2 >= 70 & 
                                     input$momdiedage2 < 80)
    vec$momdiedage2_6 = as.numeric(input$momalv=="0" & input$momdiedage2 >= 80 & 
                                     input$momdiedage2 < 90)
    vec$momdiedage2_7 = as.numeric(input$momalv=="0" & input$momdiedage2 >= 90)
    vec$daddiedage2_1 = as.numeric(input$dadalv=="0" & input$daddiedage2 < 40)
    vec$daddiedage2_2 = as.numeric(input$dadalv=="0" & input$daddiedage2 >= 40 & 
                                     input$daddiedage2 < 50)
    vec$daddiedage2_3 = as.numeric(input$dadalv=="0" & input$daddiedage2 >= 50 & 
                                     input$daddiedage2 < 60)
    vec$daddiedage2_4 = as.numeric(input$dadalv=="0" & input$daddiedage2 >= 60 & 
                                     input$daddiedage2 < 70)
    vec$daddiedage2_5 = as.numeric(input$dadalv=="0" & input$daddiedage2 >= 70 & 
                                     input$daddiedage2 < 80)
    vec$daddiedage2_6 = as.numeric(input$dadalv=="0" & input$daddiedage2 >= 80 & 
                                     input$daddiedage2 < 90)
    vec$daddiedage2_7 = as.numeric(input$dadalv=="0" & input$daddiedage2 >= 90)
    vec$MIREL = as.numeric(input$MIREL=="1")
    vec$BKBONREL = as.numeric(input$BKBONREL=="1")
    vec$bkhipmom = as.numeric(input$BKBONREL=="1" & input$bkhipmom=="1")
    vec$bkhipdad = as.numeric(input$BKBONREL=="1" & input$bkhipdad=="1")
    vec$brcafrel2 = as.numeric(input$CANCFREL=="1" & input$brcafrel2=="1")
    vec$colofrel2 = as.numeric(input$CANCFREL=="1" & input$colofrel2=="1")
    vec$colomrel2 = as.numeric(input$CANCMREL=="1" & input$colomrel2=="1")
    vec$mimomage2_1 = as.numeric(input$MIREL=="1" & input$mimom=="1" &
                                   input$mimomage == "< 55 years")
    vec$mimomage2_2 = as.numeric(input$MIREL=="1" & input$mimom=="1" &
                                   input$mimomage == "55-64 years")
    vec$mimomage2_3 = as.numeric(input$MIREL=="1" & input$mimom=="1" &
                                   (input$mimomage == "65 or older" |
                                      input$mimomage == "Not sure"))
    vec$strkreln2_1 = as.numeric(input$STRKREL=="1" & input$strkreln2 == "1")
    vec$strkreln2_2 = as.numeric(input$STRKREL=="1" & 
                                   input$strkreln2 == "2")
    vec$CANCFREL = as.numeric(input$CANCFREL=="1" )
    vec$CANCMREL = as.numeric(input$CANCMREL=="1" )
    vec$midadage2_1 = as.numeric(input$MIREL=="1" & input$midad=="1" &
                                   input$midadage == "< 55 years")
    vec$midadage2_2 = as.numeric(input$MIREL=="1" & input$midad=="1" &
                                   input$midadage == "55-64 years")
    vec$midadage2_3 = as.numeric(input$MIREL=="1" & input$midad=="1" &
                                   (input$midadage == "65 or older" |
                                      input$midadage == "Not sure"))
    vec$LACTDIET = as.numeric(input$LACTDIET=="1")
    vec$TEPIWK = as.numeric(substr(input$tepi,1,2))
    vec$TEPIWK[vec$TEPIWK == 5] = 6
    tmin = 0
    tmin[input$tmin == "Less than 20 minutes"] = 10
    tmin[input$tmin == "20-39 minutes"] = 30
    tmin[input$tmin == "40-59 minutes"] = 50
    tmin[input$tmin == "1 hour or more"] = 70
    vec$TMINWK = vec$TEPIWK*tmin
    vec$ALCOHOL_2 = as.numeric(input$ALCOHOL=="1" & input$ALCNOW=="0")
    vec$ALCOHOL_3 = as.numeric(input$ALCOHOL=="1" & input$ALCNOW=="1" &
                                 input$alcoholepi=="Never or less than once per month")
    vec$ALCOHOL_4 = as.numeric(input$ALCOHOL=="1" & input$ALCNOW=="1" &
                                 input$alcoholepi=="< 1 drink per week")
    vec$ALCOHOL_5 = as.numeric(input$ALCOHOL=="1" & input$ALCNOW=="1" &
                                 input$alcoholepi=="1-6 drinks per week")
    vec$ALCOHOL_6 = as.numeric(input$ALCOHOL=="1" & input$ALCNOW=="1" &
                                 input$alcoholepi=="> 6 drinks per week")
    vec$PACKYRS =  0   #for never smoker
    if(input$SMOKING=="1"){
      vec$PACKYRS =
        pkyrfun(smkage=input$SMOKAGE, qsmkage=input$QSMOKAGE, a = input$AGE,
                smking=input$SMOKING, smknw=input$SMOKNOW,
                smkyrs=input$smokyrs, smkct=input$smokct)
    }
    vec$qsmokage2_4 = as.numeric(input$SMOKING=="1" & input$SMOKNOW=="0" &
                                   (input$QSMOKAGE=="Less than 15" |
                                      input$QSMOKAGE=="15-19" |
                                      input$QSMOKAGE=="20-24" |
                                      input$QSMOKAGE=="25-29") )
    vec$qsmokage2_6 = as.numeric(input$SMOKING=="1" & input$SMOKNOW=="0" &
                                   (input$QSMOKAGE=="30-34" |
                                      input$QSMOKAGE=="35-39") )
    vec$qsmokage2_8 = as.numeric(input$SMOKING=="1" & input$SMOKNOW=="0" &
                                   (input$QSMOKAGE=="40-44" |
                                      input$QSMOKAGE=="45-49") )
    vec$qsmokage2_10 = as.numeric(input$SMOKING=="1" & input$SMOKNOW=="0" &
                                    (input$QSMOKAGE=="50-54" |
                                       input$QSMOKAGE=="55-59"|
                                       input$QSMOKAGE=="60 or older") )
    vec$qsmokage2_12 = as.numeric(input$SMOKING=="1" & input$SMOKNOW=="1")
    vec$PULSE30 = input$PULSE60/2
    vec$HEIGHT = input$HEIGHTin*2.54
    vec$HEIGHT[input$meas==2] = input$HEIGHTcm
    vec$WEIGHT = input$WEIGHTlb*0.453592
    vec$WEIGHT[input$meas==2] = input$WEIGHTkg
    vec$WAIST = input$WAISTin*2.54
    vec$WAIST[input$meas==2] = input$WAISTcm
    vec$HIP = input$HIPin*2.54
    vec$HIP[input$meas==2] = input$HIPcm
    vec$BMI = vec$WEIGHT/(vec$HEIGHT/100)^2
    vec$WHR = vec$WAIST/vec$HIP
    vec$SYST = input$SYST
    vec$genhel_2 = as.numeric(input$GENHEL=="Very good")
    vec$genhel_3 = as.numeric(input$GENHEL=="Good")
    vec$genhel_4 = as.numeric(input$GENHEL=="Fair or poor")
    vec$aspirin = as.numeric(input$aspirin)
    # 10 yr intxn -- index 90-92
    vec$diabXsyst = vec$DIAB*vec$SYST  # in 5C, 10C too
    vec$fh551Xalc4 = vec$fracthip55_1*vec$ALCOHOL_4
    vec$diabXhtntrt2 = vec$DIAB*vec$HTNTRT_2   # in 15C too
    # 15 yr intxns -- index 93-94
    vec$diabXpkyrs = vec$DIAB*vec$PACKYRS
    vec$fh552Xalc6 = vec$fracthip55_2*vec$ALCOHOL_6
    # 5 yr intxns -- index 95-98
    vec$diabXqsmk6 = vec$DIAB*vec$qsmokage2_6
    vec$diabXqsmk8 = vec$DIAB*vec$qsmokage2_8
    vec$AFXqsmk4 = vec$ATRIALFB*vec$qsmokage2_4
    vec$bkrelXalc6 = vec$BKBONREL*vec$ALCOHOL_6
    # 10 yr intxns - true competing risk -- index 99-100 (and diabXsyst above)
    vec$htntrt1Xqsmk12 = vec$HTNTRT_1*vec$qsmokage2_12 # in 5C too
    vec$HTNTRT2Xqsmk4 = vec$HTNTRT_2*vec$qsmokage2_4
    # 15 yr intxns - true competing risk -- index 101 (and diabXhtntrt2 above)
    vec$hicholXhtntrt2 = vec$HICHOLRP*vec$HTNTRT_2
    # 5 yr intxns - true competing risk -- (htntrtXqsmk12, diabXsyst, htntrt2xqsmok4 above)
    
    ### make sure these match - otherwise I get an NA as output from predict
    
    list(pred5 = predict(rv$fit5,cov1=vec[rv$vars5])[rv$ind5,-1] ,
         pred10 = predict(rv$fit10,cov1=vec[rv$vars10])[rv$ind10,-1],
         pred15 = predict(rv$fit15,cov1=vec[rv$vars15])[rv$ind15,-1],
         pred5C = predict(rv$fit5C,cov1=vec[rv$vars5C])[rv$ind5C,-1] ,
         pred10C = predict(rv$fit10C,cov1=vec[rv$vars10C])[rv$ind10C,-1],
         pred15C = predict(rv$fit15C,cov1=vec[rv$vars15C])[rv$ind15C,-1]
         )
  })
  
  

  output$hist5 <- renderPlot({
    plot(rv$outc5, main="", xlab=paste("Probability of",rv$evt),
         ylab="Percent of women", xlim=c(0,0.2), col="gray", freq=FALSE)
    abline(v=pred()$pred5,lwd=2, col="blue")
  })  
  output$hist10 <- renderPlot({
    plot(rv$outc10, main="", xlab=paste("Probability of",rv$evt),
         ylab="Percent of women", xlim=c(0,0.2), col="gray", freq=FALSE)
   abline(v=pred()$pred10,lwd=2, col="blue")
  })
  output$hist15 <- renderPlot({
    plot(rv$outc15, main="", xlab=paste("Probability of",rv$evt),
         ylab="Percent of women", xlim=c(0,0.2), col="gray", freq=FALSE)
    abline(v=pred()$pred15,lwd=2, col="blue")
  })
  output$hist5C <- reactivePlot(function() {
    if(rv$evt!="Death"){
      plot(rv$outc5C, main="", xlab=paste("Probability of",rv$evt,"Ever"),
           ylab="Percent of women", xlim=c(0,0.2), col="gray", freq=FALSE)
      abline(v=pred()$pred5C,lwd=2, col="blue")  }
    if(rv$evt=="Death"){return(invisible())}
  }) 
  output$hist10C <- reactivePlot(function() {
    if(rv$evt!="Death"){
      plot(rv$outc10C, main="", xlab=paste("Probability of",rv$evt,"Ever"),
           ylab="Percent of women", xlim=c(0,0.2), col="gray", freq=FALSE)
      abline(v=pred()$pred10C,lwd=2, col="blue")  }
    if(rv$evt=="Death"){return(invisible())}
  })
  output$hist15C <- reactivePlot(function() {
      if(rv$evt!="Death"){
        plot(rv$outc15C, main="", xlab=paste("Probability of",rv$evt,"Ever"),
           ylab="Percent of women", xlim=c(0,0.2), col="gray", freq=FALSE)
        abline(v=pred()$pred15C,lwd=2, col="blue")  }
    if(rv$evt=="Death"){return(invisible())}
    })

  
  output$summary5 <- renderText({ 
    ifelse(pred()$pred5 < 0.005,
           paste0("A woman with these risk factors is predicted to have less than 0.5% chance of ", tolower(rv$evt)," occurring before another event within 5 years"),
           paste0("A woman with these risk factors is predicted to have a ", round(pred()$pred5*100,0),"% chance of ", tolower(rv$evt)," occurring before another event within 5 years") )
  })
  output$summary10 <- renderText({ 
    ifelse(pred()$pred10 < 0.005,
           paste0("A woman with these risk factors is predicted to have less than 0.5% chance of ", tolower(rv$evt)," occurring before another event within 10 years"),
           paste0("A woman with these risk factors is predicted to have a ", round(pred()$pred10*100,0),"% chance of ", tolower(rv$evt)," occurring before another event within 10 years") )
  })
  output$summary15 <- renderText({ 
    ifelse(pred()$pred15 < 0.005,
           paste0("A woman with these risk factors is predicted to have less than 0.5% chance of ", tolower(rv$evt)," occurring before another event within 15 years"),
           paste0("A woman with these risk factors is predicted to have a ", round(pred()$pred15*100,0),"% chance of ", tolower(rv$evt)," occurring before another event within 15 years") )
  })
  output$summary5C <- renderText({ 
    ifelse(rv$evt=="Death", "",
     ifelse(pred()$pred5C < 0.005,
             paste0("A woman with these risk factors is predicted to have less than 0.5% chance of ", tolower(rv$evt)," within 5 years"),
             paste0("A woman with these risk factors is predicted to have a ", round(pred()$pred5C*100,0),"% chance of ", tolower(rv$evt)," within 5 years") )
      )
 })
  output$summary10C <- renderText({ 
    ifelse(rv$evt=="Death", "",    
      ifelse(pred()$pred10C < 0.005,
             paste0("A woman with these risk factors is predicted to have less than 0.5% chance of ", tolower(rv$evt)," within 10 years"),
             paste0("A woman with these risk factors is predicted to have a ", round(pred()$pred10C*100,0),"% chance of ", tolower(rv$evt)," within 10 years") )
    )
  })
  output$summary15C <- renderText({ 
    ifelse(rv$evt=="Death", "",
       ifelse(pred()$pred15C < 0.005,
              paste0("A woman with these risk factors is predicted to have less than 0.5% chance of ", tolower(rv$evt)," within 15 years"),
              paste0("A woman with these risk factors is predicted to have a ", round(pred()$pred15C*100,0),"% chance of ", tolower(rv$evt)," within 15 years") )
    )
  })
  # removed decimal place based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3222170/

  output$interp = renderUI(HTML("<br/> <b>Understanding the Figures</b>  <br/> <br/> 
The first prediction and figure display a woman's 
probability of having the selected event (for example, stroke) before any of the 
other events considered (\"event first\") during the time period selected. The 
seven events considered are myocardial infarction, stroke, lung cancer, breast 
cancer, colorectal cancer, hip fracture, and death. <br/> <br/>
The second prediction and figure display a woman's probability of having the 
selected event at any time during the time period selected (\"event ever\"). <br/> <br/>
In both figures, the blue line shows the predicted chance of having the event 
within the time period selected (5, 10, or 15 years). The line is shown along 
with the distribution of predicted probabilities observed in the women who 
participated in the study used to build the models. <br/> <br/>
Note: Separate models are used to estimate the first and second predictions. As 
a result, it is possible that the predicted probability of experiencing an event 
ever could be lower than the predicted probability of experiencing the event 
first. <br/> <br/>
<b>About the model</b>  <br/> <br/>
Citation for publication to be added. <br/> <br/>
<b>Acknowledgment</b>  <br/> <br/>
The predictions are based on the data provided by women in the Women's Health 
Initiative  (www.whi.org). We gratefully acknowledge the dedicated efforts of 
the WHI participants and of 
key WHI Investigators and staff at the clinical centers and the Clinical 
Coordinating Center.  <br/> <br/>
<b>Sources of Funding</b>  <br/> <br/>
The researchers who developed this prediction app and the Women’s Health Initiative 
programs are funded by the National Heart, Lung, and Blood Institute, NIH, 
US Department of Health and Human Services through contracts, HHSN268201100046C, 
HHSN268201100001C, HHSN268201100002C, HHSN268201100003C, HHSN268201100004C.  <br/> <br/>
<b>Disclaimer</b>  <br/> <br/> 
This website provides risk predictions for postmenopausal women. These predictions are not a substitute for professional advice; 
health care providers who use the information provided here should exercise 
their own clinical judgment as to the advice they provide. Patients who use the 
predictions do so at their own risk. Individuals with any type of medical 
condition are specifically cautioned to seek professional medical advice before 
beginning any sort of health treatment. For medical concerns, including 
decisions about medications and other treatments, users should always consult 
their health care provider. The contents of the website, such as text, graphics 
and images are for informational purposes only. Your reliance upon information 
and content obtained by you at or through this site is solely at your own risk. 
We do not assume any liability or responsibility for damage or injury (including 
death) to you, other persons or property arising from any use of any product, 
information, idea or instruction contained in the content or services provided 
to you. We cannot and will not be held legally, financially, or medically 
responsible for decisions made using this risk calculator.                                      <br/>"))
 
  ## eventually add a citation to our methods
  ## The predictions were generated using the models described in (cite pub)
  
}

# Run the application 
shinyApp(ui = ui, server = server)
