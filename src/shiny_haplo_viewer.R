#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   shiny_haplo_viewer.R
# Description:   HaploViewer using Shiny
# Author:        Brandon Monier
# Created:       2020-02-27 at 18:41:07
# Last Modified: 2020-02-27 at 19:24:21
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to create a Shiny
#    application that will project the number of consensus
#    haplotypes for a given reference range.
#
#    NOTE: run `btm_rphg_test.R` first for actual PHG data
#--------------------------------------------------------------------

# === Preamble ======================================================

## Load packages ----
library(ggthemes)
library(plotly)
library(magrittr)
library(shiny)



# === Toy data generation ===========================================

## Parameters ----
rrSamples <- 300


## Random generation of samples (reference ranges) ----
seqnames <- rep(paste0("1", LETTERS[1:3]), rrSamples / 3)
seqnames <- seqnames[order(seqnames)]

w <- sample(250:1500, rrSamples, TRUE)
start <- seq(1, 5e6, by = 3000)[1:300]
end <- start + w


## Data frame generation (S4 Vectors::DataFrame) ----
testHaplo <- S4Vectors::DataFrame(
    refRange_id   = factor(paste0("R", seq_len(rrSamples))),
    seqnames      = factor(seqnames),
    start         = as.integer(start),
    end           = as.integer(end),
    width         = as.integer((end - start) + 1),
    numHaplotypes = as.integer(sample(1:10, size = rrSamples, replace = TRUE))
)
rm(w, start, end, seqnames, rrSamples)



# === Shiny application =============================================

## User interface ----
ui <- fluidPage(
    ## Title
    shiny::h4("HaploViewer v0.0.1"),

    shiny::mainPanel(
        ## Drop down for chromosomes
        # selectInput(
        #     inputId = "chrom_select",
        #     label = "Select Chromosome",
        #     choices = NULL
        # ),
        selectInput(
            "chrom_select",
            "Select Chromosome",
            choices = levels(testHaplo$seqnames),
            selected = levels(testHaplo$seqnames)[1]
        ),

        ## Interactive plot declaration
        div(plotly::plotlyOutput("hapPlot", height = "50%"), align = "center")
    )
)


## Server generation ----
server <- function(input, output, session) {

    ## Get levels of chromosomes
    choices_chrom <- reactive({
        as.vector(levels(testHaplo$seqnames))
    })

    ## Haploplot
    output$hapPlot <- plotly::renderPlotly({
        tmp <- as.data.frame(testHaplo)
        tmp <- tmp[which(tmp$seqnames == as.character(input$chrom_select)), ]

        # Shape proportions
        yfrac <- 0.1
        xfrac <- 0.001

        # Add shape data
        tmp$med <- apply(tmp[, 3:4], 1, stats::median)
        tmp$color <- "#91baff"
        tmp[seq(1, nrow(tmp), by = 2),]$color <- "#3e619b"

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
                    ymax = max(tmp$numHaplotypes),
                    text = paste0(
                        "<b>Ref Range ID: </b>", refRange_id, "\n",
                        "<b>Start: </b>", start, "\n",
                        "<b>Stop: </b>", end
                    )
                ),
                fill = tmp$color,
                alpha = 0.5
            ) +
            geom_point(
                mapping = aes(
                    med,
                    numHaplotypes,
                    text = paste0("<b>No. Haplotypes: </b>", numHaplotypes)
                )
            ) +
            geom_path(aes(med, numHaplotypes)) +
            xlab("Physical Position (bp)") +
            ylab("Number of Haplotypes")

        ## GG Plotly
        m <- list(
            l = 100,
            r = 10,
            b = 50,
            t = 50,
            pad = 2
        )
        ggplotly(p, tooltip = "text") %>%
            config(displayModeBar = FALSE) %>%
            layout(
                autosize = F,
                xaxis = list(
                    # rangeselector = list(
                    #     buttons = list(
                    #         list(
                    #             count = 3,
                    #             label = "3 mo",
                    #             step = "month",
                    #             stepmode = "backward"
                    #         ),
                    #         list(
                    #             count = 6,
                    #             label = "6 mo",
                    #             step = "month",
                    #             stepmode = "backward"
                    #         ),
                    #         list(
                    #             count = 1,
                    #             label = "1 yr",
                    #             step = "year",
                    #             stepmode = "backward"
                    #         ),
                    #         list(
                    #             count = 1,
                    #             label = "YTD",
                    #             step = "year",
                    #             stepmode = "todate"),
                    #         list(step = "all")
                    #     )
                    # ),
                    width = 10000,
                    rangeslider = list(),
                    margin = m
                )
            )
    })

}


## Launch the application ----
shinyApp(ui = ui, server = server)


