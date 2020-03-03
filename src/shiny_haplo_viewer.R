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
rrSamples <- 10e4


## Random generation of samples (reference ranges) ----
seqnames <- rep(paste0("1", LETTERS[1:3]), rrSamples / 3)
seqnames <- seqnames[order(seqnames)]

w <- sample(2500:15000, rrSamples, TRUE)
start <- seq(1, 500e7, by = 3000)
start <- sample(start, size = rrSamples, replace = FALSE)
start <- start[order(start)]
end <- start + w


## Data frame generation (S4 Vectors::DataFrame) ----
testHaplo <- S4Vectors::DataFrame(
    refRange_id   = factor(paste0("R", seq_len(rrSamples))),
    seqnames      = factor(seqnames),
    start         = start,
    end           = end,
    width         = (end - start) + 1,
    numHaplotypes = sample(1:10, size = rrSamples, replace = TRUE)
)
rm(w, start, end, seqnames, rrSamples)



# === Shiny application =============================================

## User interface ----
ui <- fluidPage(
    ## Title
    shiny::h4("HaploViewer v0.0.3"),

    ## Interactive plot declaration
    div(plotly::plotlyOutput("hapPlot"), align = "center"),

    ## Drop down for chromosomes
    shiny::column(
        width = 3,
        shiny::selectInput(
            inputId  = "chrom_select",
            label    = "Select Chromosome",
            choices  = levels(testHaplo$seqnames),
            selected = levels(testHaplo$seqnames)[1]
        )
    ),
    shiny::column(
        width = 1,
        style = "margin-top: 25px; margin-right: 0px",
        shiny::actionButton(
            inputId = "zoom_in",
            label = "",
            icon = icon("search-plus")
        )
    ),
    shiny::column(
        width = 1,
        style = "margin-top: 25px; margin-left: 0px",
        shiny::actionButton(
            inputId = "zoom_out",
            label = "",
            icon = icon("search-minus")
        )
    ),
    shiny::column(
        width = 1,
        style = "margin-top: 25px; margin-right: 0px",
        shiny::actionButton(
            inputId = "mov_left",
            label = "",
            icon = icon("arrow-left")
        )
    ),
    shiny::column(
        width = 1,
        style = "margin-top: 25px; margin-left: 0px",
        shiny::actionButton(
            inputId = "mov_right",
            label = "",
            icon = icon("arrow-right")
        )
    )
)


## Server generation ----
server <- function(input, output, session) {

    ## Get levels of chromosomes
    choices_chrom <- reactive({
        as.vector(levels(testHaplo$seqnames))
    })


    ## Zoom in
    counter <- reactiveValues(countervalue = 1)
    observeEvent(input$zoom_in, {
        counter$countervalue <- counter$countervalue - 0.1     # if the add button is clicked, increment the value by 1 and update it
    })
    observeEvent(input$zoom_out, {
        counter$countervalue <- counter$countervalue + 0.1  # if the sub button is clicked, decrement the value by 1 and update it
    })


    ## Nav left and right
    counter_mov <- reactiveValues(countervalue = 0)
    observeEvent(input$mov_left, {
        counter$countervalue <- counter$countervalue - 3e3     # if the add button is clicked, increment the value by 1 and update it
    })
    observeEvent(input$mov_right, {
        counter$countervalue <- counter$countervalue + 3e3  # if the sub button is clicked, decrement the value by 1 and update it
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


        # plotly test ----
        fig <- plot_ly(tmp, x = ~med)
        fig <- fig %>%
            add_trace(
                y = ~numHaplotypes,
                name = "Number of haplotypes",
                type = "scatter",
                mode = "markers",
                text = paste0(
                    "<b>No. Haplotypes: </b>", tmp$numHaplotypes, "\n",
                    "<b>Ref Range ID: </b>", tmp$refRange_id, "\n",
                    "<b>Start: </b>", tmp$start, "\n",
                    "<b>Stop: </b>", tmp$end
                ),
                hoverinfo = "text"
            )

        # plotly doesn't know vectorization?
        shape_ls <- list()
        for (i in seq_len(nrow(tmp))) {
            shape_ls[[i]] <- list(
                type = "rect",
                line = list(color = tmp$color[i]),
                fillcolor = tmp$color[i],
                x0 = tmp$start[i],
                x1 = tmp$end[i],
                y0 = 0,
                y1 = max(tmp$numHaplotypes),
                opacity = 0.4
            )
        }

        fig <- fig %>%
            layout(
                shapes = shape_ls,
                xaxis = list(
                    title = "Physical position (bp)",
                    rangeselector = list(
                        buttons = list(
                            count = 1000,
                            label = "TEST"
                        )
                    ),
                    rangeselector = list(
                        buttons = list(
                            list(
                                count = 3,
                                label = "3 mo",
                                step = "month",
                                stepmode = "backward"
                            )
                        )
                    )
                ),
                yaxis = list(
                    title = "Number of haplotypes"
                )
            ) %>%
            rangeslider(
                start = (0 + counter_mov$countervalue) * counter$countervalue,
                end = (1e5 + counter_mov$countervalue) * counter$countervalue
            ) %>%
            config(displayModeBar = FALSE) %>%
            toWebGL()
        fig

    })

}


## Launch the application ----
shinyApp(ui = ui, server = server)





# === DEBUG =========================================================
        # p <- ggplot(tmp) +
        #     geom_rect(
        #         mapping = aes(
        #             xmin = start,
        #             xmax = end,
        #             ymin = 0,
        #             ymax = max(tmp$numHaplotypes)
        #         ),
        #         fill = tmp$color,
        #         alpha = 0.5
        #     ) +
        #     geom_point(
        #         mapping = aes(
        #             med,
        #             numHaplotypes,
        #             text = paste0(
        #                 "<b>No. Haplotypes: </b>", numHaplotypes, "\n",
        #                 "<b>Ref Range ID: </b>", refRange_id, "\n",
        #                 "<b>Start: </b>", start, "\n",
        #                 "<b>Stop: </b>", end
        #             )
        #         )
        #     ) +
        #     geom_path(aes(med, numHaplotypes)) +
        #     xlab("Physical Position (bp)") +
        #     ylab("Number of Haplotypes")
        #
        # ## GG Plotly
        # m <- list(
        #     l = 100,
        #     r = 10,
        #     b = 50,
        #     t = 50,
        #     pad = 2
        # )
        # ggplotly(p, tooltip = "text") %>%
        #     config(displayModeBar = FALSE) %>%
        #     layout(
        #         autosize = FALSE,
        #         xaxis = list(
        #             rangeselector = list(
        #                 buttons = list(
        #                     # list(
        #                     #     count = 3,
        #                     #     label = "3 mo",
        #                     #     step = "month",
        #                     #     stepmode = "backward"
        #                     # ),
        #                     list(
        #                         count = 1000,
        #                         label = "YTD",
        #                         # step = "linear",
        #                         stepmode = "backward"
        #                     )
        #                 )
        #             ),
        #             # width = 10000,
        #             # rangeslider = list(start = 0, end = 100),
        #             margin = m
        #         )
        #     ) %>%
        #     rangeslider(start = 0, end = 1e5)
