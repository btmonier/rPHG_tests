#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   shiny_haplo_viewer.R
# Description:   HaploViewer using Shiny
# Author:        Brandon Monier
# Created:       2020-02-27 at 18:41:07
# Last Modified: 2020-02-27 at 18:48:39
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to create a Shiny
#    application that will project the number of consensus
#    haplotypes for a given reference range.
#--------------------------------------------------------------------

# === Preamble ======================================================

## Load packages ----
library(ggthemes)
library(plotly)
library(magrittr)
library(rPHG)
library(shiny)


## Make configuration file ----
rPHG:::configFileMaker(
    dbName = "phgSmallSeq.db",
    dbType = "sqlite",
    exportPath = paste0(getwd(), "/data")
)


## Get config file name
configFilePath <- paste0("data/configFilePHG.txt")


## Various graph building methods ----
method1 <- "GATK_PIPELINE"
method2 <- "CONSENSUS"



# === Get PHG data ==================================================

## Start logger ----
rPHG::startLogger(fullPath = paste0(getwd(), "/output"))


## Build the PHG (see rPHG graphBuilder code) ----
phgObject <- rPHG::graphBuilder(
    configFile = configFilePath,
    methods = method2 # CONSENSUS
)

## Get number of haplotypes in range ----
numHaplo <- rPHG::numHaploPerRange(phgObject = phgObject)



# === Shiny application =============================================

## User interface ----
ui <- fluidPage(
    ## Title
    shiny::h4("HaploViewer v0.0.1"),

    shiny::mainPanel(
        ## Drop down for chromosomes
        selectInput(
            inputId = "chrom_select",
            label = "Select Chromosome",
            choices = NULL
        ),

        ## Interactive plot declaration
        plotly::plotlyOutput("hapPlot")
    )
)


## Server generation ----
server <- function(input, output, session) {

    ## Get levels of chromosomes
    choices_chrom <- reactive({
        as.vector(levels(numHaplo$seqnames))
    })

    ## Update chromosome selection
    observe({
        updateSelectInput(
            session = session,
            inputId = "chrom_select",
            choices = choices_chrom()
        )
    })

    ## Haploplot
    output$hapPlot <- plotly::renderPlotly({
        tmp <- as.data.frame(numHaplo)

        # Shape proportions
        yfrac <- 0.1
        xfrac <- 0.001

        # Add shape data
        tmp$med <- apply(tmp[, 3:4], 1, stats::median)
        tmp$color <- "#91baff"
        tmp[seq(1, nrow(tmp), by = 2),]$color <- "#3e619b"

        tmp$numHaplotypes <- sample(x = 1:5, size = nrow(tmp), replace = TRUE)

        # Get limit data
        xbeg <- min(tmp$start)
        xend <- max(tmp$end)
        yend <- max(tmp$numHaplotypes)

        p <- ggplot(tmp) +
            geom_rect(
                mapping = aes(
                    xmin = start,
                    xmax = end,
                    ymin = 0,
                    ymax = max(tmp$numHaplotypes)
                ),
                fill = tmp$color,
                alpha = 0.4
            ) +
            geom_point(aes(med, numHaplotypes, text = paste0("Ref Range ID: ", refRange_id))) +
            geom_path(aes(med, numHaplotypes)) +
            xlab("Physical Position (bp)") +
            ylab("Number of Haplotypes")
        ggplotly(p, tooltip = "text")

        # print(
        #     ggplotly(
        #         ggplot(data = tmp) +
        #             ylim(-(yend * yfrac), yend) +
        #             scale_x_continuous(limits = c(xbeg, xend)) +
        #             geom_rect(
        #                 mapping = aes(
        #                     xmin = .data$start,
        #                     xmax = .data$end,
        #                     ymin = 0,
        #                     ymax = -(yend * yfrac)
        #                 ),
        #                 fill = tmp$color
        #             ) +
        #             geom_path(aes(x = .data$med, y = .data$numHaplotypes)) +
        #             geom_point(
        #                 aes(
        #                     x = .data$med,
        #                     y = .data$numHaplotypes
        #                 ),
        #                 size = 1
        #                 ) +
        #             facet_grid(seqnames ~ .) +
        #             xlab("Physical Position (bp)") +
        #             ylab("Number of Haplotypes")
        #     )
        # )
    })

}


## Launch the application ----
shinyApp(ui = ui, server = server)


