
#' Create DiNeR object
#'
#' This function creates a DiNeR object
#'
#' @param edge_list A list of data frames; each data frame is exactly two columns, representing the two nodes connected
#' @return A DiNeR object consisting of a node adjacency list for each sample and XY coordinates for nodes.
#' @export
create_diff_network_object <- function(edge_list){

  # error checking to assure we have exactly 2 columns to identify node-pair relationships from
  adequate_columns = unlist(lapply(edge_list, function(x){
    if (ncol(x) == 2){
      return (0)
    }  else if (ncol(x) < 2){
      return (-1)
    } else{
      return (1)
    }
  }))

  # cannot continue if any of the lists don't have at least two columns
  if (-1 %in% adequate_columns){
    stop( paste("The following edge lists do not have enough columns: ", paste(which (adequate_columns == -1), collapse=",") ))
  }

  # can continue if the lists have too many columns, but output may not be accurate
  if (1 %in% adequate_columns){
    warning( paste("The following edge lists have too many columns: ", paste(which (adequate_columns == 1), collapse=",") , ". Using the first 2 columns to determine node relationships."))
  }


  # identify all distinct node types
  node_labels <- unique(c(unlist(lapply(edge_list, "[[", 1)),
                          unlist(lapply(edge_list, "[[", 2))))

  # for now: square layout
  n_row_col <- ceiling(sqrt(length(node_labels)))

  spacing = 10

  # the dataframe with x y coords for nodes
  node_xy_df <- data.frame(matrix(nrow=0,ncol=3), stringsAsFactors=FALSE)

  for (row in 1:n_row_col){
    for (col in 1:n_row_col){

      if (col + (row-1)*n_row_col > length(node_labels)){
        break
      }

      node_label <- node_labels[col + (row-1)*n_row_col]

      x=0
      y=row*spacing

      if (row %% 2 == 1){
        x = col*spacing
      } else{
        x = col*spacing + spacing/2
      }

      node_xy_df <- rbind(node_xy_df, c(node_label,x,y), stringsAsFactors=FALSE)
    }
  }

  colnames(node_xy_df) <- c("Node", "X", "Y")
  node_xy_df$X <- as.numeric(node_xy_df$X)
  node_xy_df$Y <- as.numeric(node_xy_df$Y)
  node_xy_df$radius = 1

  edge_list = lapply(edge_list, function(x){
    tmp = x
    colnames(tmp) = c("N1", "N2")
    return(tmp)
  })

  return(list(
    "EdgeList" = edge_list,
    "NodeXY" = node_xy_df
  ))
}

#' Add features for individual networks
#'
#' This function adds features for individuals (i.e. for each network) to the input DiNeR object
#'
#' @param diff_network_object The DiNeR object to add individual features to
#' @param df A data.frame containing the features to add. Each column represents a feature, and number of rows should be equal to the number of networks in the object.
#' @return A DiNeR object with the individual features added
#' @export
add_individual_features <- function(diff_network_object, df){
  # sanity check
  if (nrow(df) != length(diff_network_object$EdgeList)){
    stop("Feature DF must contain info for all networks")
  }

  # if we already have features, combine
  if ("IndividualFeatures" %in% names(diff_network_object)){
    diff_network_object$IndividualFeatures = cbind(diff_network_object$IndividualFeatures, df)
  } else { # otherwise....
    diff_network_object$IndividualFeatures = df
    row.names(df) = 1:nrow(df)
  }

  return (diff_network_object)
}

#' Add edge features
#'
#' This function adds features for the edges in each network
#'
#' @param diff_network_object The DiNeR object to add edge features to
#' @param df_list A list of data.frame objects. Within each data.frame, each column represents a feature, and the rows signify the edge the feature belongs to.
#' @return A DiNeR object with the edge features added
#' @export
add_edge_features <- function(diff_network_object, df_list){
  # if we already have features, combine
  if ("EdgeFeatures" %in% names(diff_network_object)){
    if (length(df_list) != length(diff_network_object$EdgeList)){
      stop("DF lists not the same size")
    }
    diff_network_object$EdgeFeatures = lapply(1:length(df_list), function(i){
      cbind(diff_network_object$EdgeFeatures[[i]], df_list[[i]])
    })
  } else { # otherwise....
    diff_network_object$EdgeFeatures = df_list
  }

  return (diff_network_object)


}

#' Add node features
#'
#' This function adds features for the node in each network
#'
#' @param diff_network_object The DiNeR object to add node features to
#' @param df_list A list of data.frame objects. Within each data.frame, each column represents a feature, and the rows signify the node the feature belongs to.
#' @return A DiNeR object with the edge features added
#' @export
add_node_features <- function(diff_network_object, df){
  # if we already have features, combine
  if ("NodeFeatures" %in% names(diff_network_object)){
    if (nrow(df) != nrow(diff_network_object$NodeXY)){
      stop(paste("Number of rows in DF (", nrow(df)  ,") is not equal to the number of nodes, (",nrow(diff_network_object$NodeXY), ") in object." ))
    }
    # IMPORTANT: must be in the same order!!!
    diff_network_object$NodeFeatures = cbind(diff_network_object$NodeFeatures, df)
  } else { # otherwise....
    diff_network_object$NodeFeatures = df
  }
}
