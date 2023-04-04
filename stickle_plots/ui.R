

library(shiny)
library(shinythemes)

shinyUI(navbarPage("Imputation Plots",
                   theme = shinytheme("cerulean"),
                   tabPanel("Imputation Results",
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("variable", "Variable:",
                                             c("Standard Length" = "stl",
                                               "Interpelvic Spine Proportion" = "ips.prop",
                                               "Lower Pelvic Spine" = "lps",
                                               "Total Pectoral Girdle" = "tpg",
                                               "Maximum Ventral Angle" = "mav",
                                               "Maximum Caudal Angle" = "mcv",
                                               "Maximum Dorsal Fin" = "mdf",
                                               "Maximum Anal Fin" = "maf",
                                               "Minimum Dorsal Spine" = "mds",
                                               "Dorsal Spine 1" = "ds1",
                                               "Dorsal Spine 2" = "ds2",
                                               "Dorsal Spine 3" = "ds3",
                                               "Midpoint Pectoral Fin" = "mpt",
                                               "Lower Pectoral Fin" = "lpt",
                                               "Ectocoracoid Length" = "ect",
                                               "Cleithrum Length" = "cle",
                                               "Premaxillary Length" = "pmx"),
                                             inline = TRUE),

                              ),

                              mainPanel(
                                tabsetPanel(
                                  tabPanel("MICE Imputation", plotOutput("plot_mice")),
                                  tabPanel("mixgb Imputation", plotOutput("plot_mixgb"))
                                )
                              )
                            )
                   ),

                   tabPanel("Accuracy Plots",
                            mainPanel(
                              plotOutput("accuracy_plot_mice"),
                              verbatimTextOutput("auc_mice"),
                              plotOutput("accuracy_plot_mixgb"),
                              verbatimTextOutput("auc_mixgb"),
                              plotOutput("accuracy_plot"),
                              verbatimTextOutput("auc_missForest")
                            )
                   )
))
