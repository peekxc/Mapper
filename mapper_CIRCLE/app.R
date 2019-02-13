
library("shiny")
library("Mapper")
library("grapher")

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel(""),
    fluidPage(
        fluidRow( grapher::grapherOutput("grapher"), style = "height: 60vh !important;"), 
        fluidRow(
            column(6, plotOutput(outputId = "data_out", brush = "data_brush")),
            column(6, plotOutput(outputId = "filter_out", brush = brushOpts(id = "filter_brush", direction = "x"))), 
            style = "height: 40vh !important;"
        )
    )
)

## Generate noisy points around the perimeter of a circle 
data("noisy_circle", package = "Mapper")
X <- as.data.frame(noisy_circle)
colnames(X) <- c("x", "y")

## Filter points
left_pt <- noisy_circle[which.min(noisy_circle[, 1]),]
f_x <- matrix(apply(noisy_circle, 1, function(pt) (pt - left_pt)[1]))
f_x_jitter <- structure(data.frame(cbind(f_x, runif(nrow(noisy_circle)))), names = c("f(X)", "Jitter"))
 

## Make the mapper
m <- mapper(X = noisy_circle, filter_values = f_x, 
            cover_params = list(typename="restrained rectangular", 
                                number_intervals=8L, percent_overlap=50),
            measure = "euclidean", 
            cluster_params = list(cl="single", num_bins=10L), return_reference = TRUE)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$grapher <- grapher::renderGrapher({ m$as_grapher() %>% grapher::enableForce() })
    
    ## Filter plot 
    output$filter_out <- renderPlot({
        brush_data <- brushedPoints(X, input$data_brush, allRows = TRUE)
        pt_color <- bin_color(f_x_jitter[["f(X)"]]) 
        if (all(brush_data[["selected_"]] == FALSE) || is.null(brush_data)){
            # plot(f_x_jitter, pch = 20, cex = 1.25, col = pt_color)
            ggplot() + 
                geom_point(data = f_x_jitter, aes(x=`f(X)`,y=Jitter), color = pt_color) + 
                ggplot2::xlab("f(X)") + ggplot2::ggtitle("Result of map") + theme_bw() + 
                theme(axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)) + 
                theme(plot.title = element_text(size=35)) + 
                theme(plot.title = element_text(hjust = 0.5))
        } else {
            # selected_idx <- which(brush_data[["selected_"]])
            pt_color[!brush_data[["selected_"]]] <- adjustcolor("black", alpha.f = 0.25)
            # plot(f_x_jitter, pch = 20, cex = 1.25, col = pt_color)
            ggplot() + 
                geom_point(data = f_x_jitter, aes(x=`f(X)`,y=Jitter), color = pt_color) + 
                ggplot2::xlab("f(X)") + ggplot2::ggtitle("Result of map") + theme_bw() + 
                theme(axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)) + 
                theme(plot.title = element_text(size=35)) + 
                theme(plot.title = element_text(hjust = 0.5))
        }
    })
    
    output$data_out <- renderPlot({
        brush_data <- brushedPoints(f_x_jitter, input$filter_brush, allRows = TRUE)
        pt_color <- bin_color(f_x_jitter[["f(X)"]]) 
        
        if (all(brush_data[["selected_"]] == FALSE) || is.null(brush_data)){
            ggplot(X, aes(x=x, y=y)) + geom_point(color="black") + theme_bw() + 
                ggplot2::ggtitle("The data") + 
                theme(axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)) + 
                theme(plot.title = element_text(size=35)) + 
                theme(plot.title = element_text(hjust = 0.5))
        } else {
            selected_idx <- which(brush_data[["selected_"]])
            pt_color[!brush_data[["selected_"]]] <- adjustcolor("black", alpha.f = 0.25)
            ggplot(X, aes(x=x, y=y)) + geom_point(color=pt_color) + theme_bw() + 
                ggplot2::ggtitle("The data") + 
                theme(axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)) + 
                theme(plot.title = element_text(size=35)) + 
                theme(plot.title = element_text(hjust = 0.5))
        }

        # if (all(brush_data[["selected_"]] == FALSE)){
        #     ggplot(X, aes(x=x, y=y)) + geom_point(color=pt_color) + theme_bw()
        # } else {
        #     selected_idx <- which(brush_data[["selected_"]])
        #     pt_color[!selected_idx] <- adjustcolor("black", alpha.f = 0.25)
        #     ggplot(X, aes(x=x, y=y)) + geom_point(color=pt_color) + theme_bw()
        # }
    })
    
    vertex_filter_val <- sapply(sapply(m$vertices, function(v_idx){ 
        apply(as.matrix(m$cover$filter_values[v_idx,]), 1, mean)
    }), mean)
    default_colors <- bin_color(vertex_filter_val)
    observeEvent(input$filter_brush, {
        req(input$filter_brush)
        brush_data <- brushedPoints(f_x_jitter, input$filter_brush, allRows = TRUE)
        v_col <- switch(as.character(all(brush_data[["selected_"]] == FALSE)), 
                        "TRUE" = default_colors, 
                        "FALSE" = {
                            selected_idx <- which(brush_data[["selected_"]])
                            v_weight <- sapply(m$vertices, function(v) sum(v %in% selected_idx))
                            bin_color(v_weight, col_pal = viridis::inferno(100, begin = 0.15, end = 0.8, alpha = seq(0.65, 1, length.out = 100)))
                        })
        grapher::setNodeColor("grapher", color = v_col)
    })
    
    
    observeEvent(input$data_brush, {
        brush_data <- brushedPoints(X, input$data_brush, allRows = TRUE)
        v_col <- switch(as.character(all(brush_data[["selected_"]] == FALSE)), 
                        "TRUE" = default_colors, 
                        "FALSE" = {
                            selected_idx <- which(brush_data[["selected_"]])
                            v_weight <- sapply(m$vertices, function(v) sum(v %in% selected_idx))
                            bin_color(v_weight, col_pal = viridis::inferno(100, begin = 0.15, end = 0.8, alpha = seq(0.65, 1, length.out = 100)))
                        })
        grapher::setNodeColor("grapher", color = v_col)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
