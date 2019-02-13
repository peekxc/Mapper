library("shiny")
library("Mapper")
library("promises")

## Load chemical diabetes data set
library(locfit)
data("chemdiab", package = "locfit")

## L-Infinity eccentricity 
X <- as.matrix(chemdiab[, -6])
f_ecc <- matrix(apply(as.matrix(dist(X)), 1, max))

## Make the mapper 
m <- Mapper::MapperRef$new(X = X)
m$use_cover(f_ecc, typename = "fixed rectangular", number_intervals = 8L, percent_overlap = 35)
m$use_clustering_algorithm(cl = "single", num_bins = 10L)
m$use_distance_measure(measure = "euclidean")
m$compute_k_skeleton(k = 1L)

## Enable multiscale
m$enable_multiscale()
m$update_mapper(percent_overlap = 35)

ui <- fluidPage(
    titlePanel("Multiscale Mapper example"),
    sidebarLayout(
        sidebarPanel(
            sliderInput("overlap", label="Percent overlap", min = 0, max = 50, step = 0.1, value = m$cover$percent_overlap)
        ),
        mainPanel(
           grapher::grapherOutput("grapher")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$grapher <- grapher::renderGrapher({ m$as_grapher() })
    
    ## Observer to detect node changes
    idx_set <- m$cover$index_set
    ls_idx <- m$cover$level_sets_to_compare()
    # c_network <- reactiveVal()
    # observeEvent(input$network, { c_network <- input$network })
    observeEvent(input$overlap, {
        req(input$overlap)
        update_stats <- m$update_mapper(percent_overlap = input$overlap, stats = TRUE)
        grapher::getNetwork("grapher") ## updates the input$network
        # future::future({ c_network() }) %...>% (function(c_net) {
        #     print(c_net)
        # })
        ## Need to update the update function to allow specific level sets to be updated
        observeEvent(input$network, {
            req(input$network)
            nodes_xy <- as.matrix(data.table::rbindlist(input$network$nodes)[, c("x", "y")])
            links <- data.table::rbindlist(lapply(c_net$links, function(link){ link[c("from", "to", "color")]}))
            for (idx in idx_set){
                c_old_vertices <- update_stats$old_ls_map[[idx]]
                centroid <- colMeans(nodes_xy[with(nodes_df, { id %in% c_old_vertices}),])
                grapher::removeNodes("grapher", node_ids = c_old_vertices)
                c_new_vertices <- update_stats$new_ls_map[[idx]]
                m$simplicial_complex$adjacent_vertices(0)
                grapher::insertNodes("grapher", node_config = )
                apply(ls_idx, 1, function(li) idx %in% li)
            }
        })


# 
#         grapher::removeNodes()
#         grapher::setNodes()
#         wut$new_ls_map$`(1)`
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
