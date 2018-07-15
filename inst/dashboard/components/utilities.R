# Call this function with an input (such as `textInput("text", NULL, "Search")`) if you
# want to add an input to the navbar
dropDownItems <- function(inputs) {
  for(i in 1:length(inputs)){
    input_el <- inputs[[i]]
    input_el$attribs$class <- paste(input_el$attribs$class, "dropdown_item")
    inputs[[i]] <- input_el
  }
  inputs
}

## Utility function to make a fancy datatables object with all the cool extensions
fancytable <- function(dt){
  DT::datatable(dt, fillContainer = TRUE,
                extensions = c("Scroller", "Buttons", "ColReorder", "FixedColumns", "Scroller"),
                options = list(dom = 'Bfrtip', buttons = I('colvis'), scrollX = TRUE, colReorder = TRUE, fixedColumns = TRUE, deferRender = TRUE, scrollY = 200, scroller = TRUE)
  )
}

## Function to create a table incase a more useful one is available, default to one in mapper object otherwise
# createTable <- function(M, column_names=NULL){
#   X_df <- as.data.frame(M$X)
#   colnames(X_df) <- column_names
#   return(X_df)
# }
