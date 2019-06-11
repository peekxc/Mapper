## shiny_modules.R
## Author: Matt Piekenbrock

#' MapperVis
#' @description 'MapperVis' is a \code{\link[modules]{module}} containing \code{\link[shiny::callModule]{shiny modules}}
#' that simplify making shiny applications with Mapper.
#' @import modules
MapperVis <- modules::module({
  modules::import("utils", attach = FALSE)
  modules::import("shiny", attach = FALSE)
  
  MapperDataInput <- function(id, datasets=NULL){
    ns <- NS(id)
    data_supplied <- !(missing(datasets) || is.null(datasets))
    if (!data_supplied){
      extdata <- file.path(system.file(package="Mapper"), "extdata")
      datasets <- list(
        "circle" = function(){ readRDS(file.path(extdata, "noisy_circle.rds")) }, 
        "WVS" = function() { readRDS(file.path(extdata, "wvs_us_wave6_ex.rds")) }
      )
      shiny::selectInput(ns("data_name"), label = "Data", choices = datasets, selected = utils::head(names(datasets), 1))
    } else {
      functor <- function(symbol){ function(){ eval(parse(text=symbol), envir = .GlobalEnv) } }
      datasets <- structure(lapply(ls(.GlobalEnv), functor), names = ls(.GlobalEnv))
      shiny::selectizeInput(shiny::ns("data_name"), label = "Data", choices = datasets, multiple = FALSE, selectize = TRUE)
    }
  }

  MapperData <- function(input, output, session){
    return(shiny::reactive({ input$data_name }))
  }

  MapperFilterInput <- function(){
    shiny::selectInput(shiny::ns("filter_name"), label = "Filter", choices = c("PCA", "ICA", "Eccentricity", "MDS", "LE", "tSNE", "UMAP", "Custom"))
  }

  ShinyMapperInput <- function(id){
    # stopifnot(is.list(X), !is.null(names(X)))
    ns <- shiny::NS(id)
    shiny::tagList(
      # selectInput(ns("data_name"), label = "Input data", choices = names(X), selected = head(names(X), 1)),
      shiny::actionButton(ns("construct_mapper"), label = "Construct Mapper", width = "100%"),
      shiny::selectInput(ns("filter_name"), label = "Filter", choices = c("PCA", "ICA", "Eccentricity", "MDS", "LE", "tSNE", "UMAP", "Custom")),
      shiny::uiOutput(ns("filter_options")),
      shiny::selectInput(ns("cover_type"), label="Cover", 
                  choices = c("fixed interval", "restrained interval", "adaptive", "ball"), 
                  selected = "fixed interval"),
      shiny::uiOutput(ns("cover_options")),
      shiny::selectInput(ns("metric"), label = "Metric", choices=proxy::pr_DB$get_entry_names(), 
                  selected = "Euclidean"),
      shiny::selectInput(ns("clustering"), label="Clustering algorithm (or linkage criterion)", 
                  choices=c("ward.D", "ward.D2", "single", "complete", 
                            "average", "mcquitty", "median", "centroid", "custom"), 
                  selected = "single"),
      shiny::selectInput(ns("cutoff_heuristic"), label = "Cutoff Heuristic", 
                  choices = c("histogram", "continuous"), 
                  selected="continuous")
    )
  }

  ShinyCoverInput <- function(id){
    ns <- shiny::NS(id)
    shiny::tagList(
      shiny::conditionalPanel(condition = "(input.cover_type == 'fixed interval') || (input.cover_type == 'restrained interval')", 
        shiny::tagList(
          shiny::sliderInput(ns("number_intervals"), label = "Number intervals", min = 1, max = 50, value = 10, step = 1),
          shiny::sliderInput(ns("percent_overlap"), label = "Percent overlap", min = 0, max = 100, value = 50, step = 1) 
        )
      ), 
      shiny::conditionalPanel(condition = "input.cover_type == 'adaptive interval'", 
        shiny::tagList(
          shiny::sliderInput(ns("number_intervals"), label = "Number intervals", min = 1, max = 50, value = 10, step = 1),
          shiny::sliderInput(ns("percent_overlap"), label = "Percent overlap", min = 0, max = 100, value = 50, step = 1), 
          shiny::numericInput(ns("quantile_method"), label = "Quantile method", min = 1, max = 9, step = 1L, value = 7) 
        )
      ), 
      shiny::conditionalPanel(condition = "input.cover_type == 'ball'", 
        shiny::tagList(
          shiny::sliderInput(ns("epsilon"), label = "Radius (normalized)", min = 0, max = 1, step = 0.001)
        )
      )
    )
  }

  ShinyMapper <- function(input, output, session, dataset, cover, M = NULL){
    mapper_supplied <- !missing(M)
    if (mapper_supplied){ stopifnot(is(M, "MapperRef")) }
    M <- if (mapper_supplied) { M } else { methods::new(MapperRef) }
    MapperObj <- shiny::reactive({
      M <- MapperRef$new(dataset)
      M$cover <- cover()
      M$use_distance_measure(input$metric)
      M$use_clustering_algorithm(cl = input$clustering, cutoff_method = input$cutoff_heuristic)
      return(M$construct_k_skeleton(k=1L))
    })
    return(MapperObj)
  }

  ShinyMultiscaleMapperUI <- function(id, smaps){
    ns <- shiny::NS(id)
    bd_idx <- grep(smaps, pattern = "[#] .*")
    bd <- as.numeric(sapply(strsplit(smaps[bd_idx], split = " "), function(x) utils::tail(x, 1)))
    shiny::sliderInput(ns("eps"), label = "eps", min = min(bd), max = max(bd[bd != Inf]), value = min(bd))
  }

  ShinyMultiscaleMapper <- function(input, output, session, initial_complex, smaps, pp_shiny){
    base_complex <- simplextree::simplex_tree()
    base_complex$deserialize(initial_complex$serialize()) ## copy 
    bd_idx <- grep(smaps, pattern = "[#] .*")
    simplicial_maps <- smaps[-bd_idx]
    pp_mm <- vector("list", length(simplicial_maps))
    pp_mm_rev <- vector("list", length(simplicial_maps))
    cc <- 1L
    for (si_map in simplicial_maps){
      op <- strsplit(si_map, split = " ")[[1]]
      if (op[1] == "i"){
        simplex <- as.integer(op[-1])
        base_complex$insert(simplex)
        if (length(simplex) == 1){
          new_node <- list(nodes=data.frame(id=simplex, x=0.50, y=0.50))
          pp_mm[[cc]] <- list(call=pixiplex::insertNodes, params = new_node)
          pp_mm_rev[[cc]] <- list(call=pixiplex::removeNodes, params = list(node_ids=simplex))
        } 
        else if (length(simplex) == 2){
          new_link <- list(links=data.frame(source=simplex[1], target=simplex[2]))
          pp_mm[[cc]] <- list(call=pixiplex::insertLinks, params = new_link)
          pp_mm_rev[[cc]] <- list(call=pixiplex::removeLinks, params = new_link)
        }
        cc <- cc + 1L
      }
      else if (op[1] == "c"){
        from <- as.integer(op[2:(which(op == "t")-1)])
        to <- as.integer(utils::tail(op, 1))
        params <- list()
        params_rev <- list()
        ii <- 1L
        for (tau in from){
          tau_cofaces <- base_complex$ltraverse(tau, identity, "cofaces")
          mapped <- lapply(tau_cofaces, function(coface){
            coface[coface == tau]  <- to
            return(sort(unique(coface)))
          })
          ## Insert mapped simplices 
          for (new_simplex in mapped){
            if (length(new_simplex) == 1){
              new_node <- list(nodes = data.frame(id=new_simplex, x=0.5, y=0.5))
              params[[ii]] <- list(call=pixiplex::insertNodes, params = new_node)
              params_rev[[ii]] <- list(call=pixiplex::removeNodes, params = list(node_ids=new_simplex))
            } else if (length(new_simplex) == 2){
              new_link <- list(links = data.frame(source = new_simplex[1], target = new_simplex[2]))
              params[[ii]] <- list(call=pixiplex::insertLinks, params = new_link)
              params_rev[[ii]] <- list(call=pixiplex::removeLinks, params = list(links=new_link))
            }
            ii <- ii + 1L
          }
          ## Remove old simplices
          for (cof in tau_cofaces){
            if (length(cof) == 1){
              params[[ii]] <- list(call=pixiplex::removeNodes, params = list(node_ids = cof))
              params_rev[[ii]] <- list(call=pixiplex::insertNodes, params = list(nodes = data.frame(id=cof, x=0.5, y=0.5)))
            } else if (length(cof) == 2){
              new_link <- data.frame(source = cof[1], target = cof[2])
              params[[ii]] <- list(call=pixiplex::removeLinks, params = list(links = force(new_link)))
              params_rev[[ii]] <- list(call=pixiplex::insertLinks, params = list(links = new_link))
            }
            ii <- ii + 1L
          }
          base_complex$collapse(tau, to, to)
        } # for (tau in from)
        pp_mm[[cc]] <- params
        pp_mm[[cc]] <- params_rev
        cc <- cc + 1L
      }
    }
    
    ## Attach observer to eps slider that sends changing signals to the visualization when moved
    bd <- as.numeric(sapply(strsplit(smaps[bd_idx], split = " "), function(x) utils::tail(x, 1)))
    idx <- 1L
    observeEvent(input$eps, {
      req(pp_shiny$pp_vis())
      #browser()
      ub_idx <- ifelse(input$eps < min(bd), 0, max(which(input$eps >= bd)))
      if (idx == ub_idx){ return(); }
      else if (idx < ub_idx){
        while(idx <= ub_idx){
          if (!is.null(pp_mm[[idx]])){
            do.call(pp_mm[[idx]]$call, modifyList(list(id = pp_shiny$pp_vis()), pp_mm[[idx]]$params))
          }
          idx <<- idx + 1L
        }
      } else {
        while(idx > ub_idx){
          if (!is.null(pp_mm_rev[[idx]])){
            do.call(pp_mm_rev[[idx]]$call, modifyList(list(id = pp_shiny$pp_vis()), pp_mm_rev[[idx]]$params))
          }
          idx <<- idx - 1L
        }
      }
    })
    current_complex <- simplextree::simplex_tree()
    current_complex$deserialize(initial_complex$serialize()) ## copy 
    filtration <- Simpers::filtration(x = current_complex, sm = smaps)
    
    ## Final result
    result <- list(
      current_eps = reactive({ input$eps }), 
      filtration = reactive({ filtration }),
      eps = reactive({ bd })
    )
    return(result)
  } # ShinyMultiscaleMapper
}) # MapperVis module



