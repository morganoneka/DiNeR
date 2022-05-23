#' Split loops
#'
#' This function, used internally by DiNeR, takes a list of edges and returns a list with edges that are loops (go to and from the same node) and edges that are not loops.
#'
#' @param edge_info The edge info (i.e. the node adjacency list) for an individual network
#' @return A DiNeR object with the edge features added
#' @export
split_loops <- function(edge_info){
  edges_loops <- edge_info[which(edge_info$X.x == edge_info$X.y & edge_info$Y.x == edge_info$Y.y ), ]
  edges_loops$X.y <- edges_loops$X.y + .2
  edges_loops$Y.y <- edges_loops$Y.y + .2

  edges_nonloops <- edge_info[setdiff(1:nrow(edge_info),which(edge_info$X.x == edge_info$X.y & edge_info$Y.x == edge_info$Y.y )), ]

  return(list(edges_loops, edges_nonloops))
}


#' Get node IDs
#'
#' This function returns the names of all the node IDs present in this list of networks
#'
#' @param diff_network_object The DiNeR object to get node IDs for
#' @return A list of all node IDs present in the DiNeR object
#' @export
get_node_ids <- function(diff_network_object){
  return(diff_network_object[[2]]$Node)
}

#' Fill in missing values
#'
#' This function is used primarily by `make_piecharts`. A DiNeR object contains one row per edge (i.e., the edge between A and B is saved as A,B) but symmetric visualizations such as `make_piecharts` need an entry for both A,B and B,A. This function creates a new edge list with these symmetric entries.
#'
#' @param edge_list_df
#' @param node_ids
#' @param edge_vars
#' @return An edge list expanded to include entries for both symmetries of each interaction
#' @export
fill_missing_vals <- function(edge_list_df, node_ids, edge_vars){
  for (i in 1:length(node_ids)){
    for (j in i:length(node_ids)){
      forward = which(edge_list_df[,edge_vars$node1] == node_ids[i] & edge_list_df[,edge_vars$node2] == node_ids[j])
      backward = which(edge_list_df[,edge_vars$edge3] == node_ids[i] & edge_list_df[,edge_vars$node1] == node_ids[j])

      if (length(forward) == 0 & length(backward) == 0){
        edge_list_df <- rbind(edge_list_df[,unlist(edge_vars)], c(node_ids[i], node_ids[j], 0, ""))
      }
    }
  }

  return(edge_list_df)
}

#' Sort edge list
#'
#' This function sorts an edge list so that the node pairs are in alphabetical order, i.e. the row "B,A" becomes "A,B", then sorts the entire list
#'
#' @param edge_list_df The edge list to be sorted
#' @param edge_vars
#' @return The edge list with nodes in alphabetical order and entries sorted
#' @export
sort_list <- function(edge_list_df, edge_vars){
  edge_list_in_order <- edge_list_df[which(edge_list_df[,edge_vars$node1] <= edge_list[,edge_vars$node2]),]
  edge_list_out_of_order <- edge_list_df[which(edge_list_df[,edge_vars$node1] > edge_list[,edge_vars$node2]),]
  cnames <- colnames(edge_list_out_of_order)
  cnames[c(which(cnames == edge_vars$node1), which(cnames == edge_vars$node2))] = c(edge_vars$node2, edge_vars$node1)
  colnames(edge_list_out_of_order) <- cnames
  fixed <- rbind(edge_list_in_order[,unlist(edge_vars)], edge_list_out_of_order[,unlist(edge_vars)])
  fixed[order(fixed[,edge_vars$node1], fixed[,edge_vars$node2]),]
  return(fixed)
}

#' Split networks by network group
#'
#' Split networks into groups based on values for an individual feature
#'
#' @param diff_network_object The DiNeR object to split
#' @param group_label The name of the feature used to split into groups
#' @return The edge list with nodes in alphabetical order and entries sorted
#' @export
split_by_group <- function(diff_network_object, group_label){

  # TODO: check if group_label is in individual features

  groups = unique(diff_network_object$IndividualFeatures[,group_label])


  # TODO: return list where index is group name, value is diff network object just for that group
  group_dnos = lapply(groups, function(g){

    # indices of individuals in this group
    indices = which(diff_network_object$IndividualFeatures[,group_label] == g)

    sub_dno = list(
      EdgeList = diff_network_object$EdgeList[indices],
      NodeXY = diff_network_object$NodeXY,
      IndividualFeatures = diff_network_object$IndividualFeatures[indices,],
      EdgeFeatures = diff_network_object$EdgeFeatures[indices]
    )

    # sub_dno$NodeXY[,"Index"] = 1:length(indices)

    return(sub_dno)

  })

  return(group_dnos)
}

#' Get Possible Values for Edge Feature
#'
#' Identify all values an edge feature can have
#'
#' @param diff_network_object The DiNeR object to get feature names for
#' @param edge_feature The name of the feature we want to get values for
#' @return A list of all possible values for the specified edge feature
#' @export
get_possible_values_edge_feature <- function(diff_network_object, edge_feature=""){
  if (edge_feature == ""){
    stop("No value was entered for edge_feature. Please identify the column name of the feature you would like to obtain values for. ")
  }
  do.call(rbind, diff_network_object$EdgeFeatures)[,edge_feature] %>% unique()
}

#' Get Possible Values for Individual Feature
#'
#' Identify all values an individual featur
#'
#' @param diff_network_object The DiNeR object to get feature names for
#' @param edge_feature The name of the feature we want to get values for
#' @return A list of all possible values for the specified individual feature
#' @export
get_possible_values_individual <- function(diff_network_object, individual_feature=""){
  if (individual_feature == ""){
    stop("No value was entered for individual_feature Please identify the column name of the feature you would like to obtain values for. ")
  }
  diff_network_object$IndividualFeatures[,individual_feature] %>% unique()
}
