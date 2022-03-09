# helper function to split loops and non-loops
split_loops <- function(edge_info){
  edges_loops <- edge_info[which(edge_info$X.x == edge_info$X.y & edge_info$Y.x == edge_info$Y.y ), ]
  edges_loops$X.y <- edges_loops$X.y + .2
  edges_loops$Y.y <- edges_loops$Y.y + .2

  edges_nonloops <- edge_info[setdiff(1:nrow(edge_info),which(edge_info$X.x == edge_info$X.y & edge_info$Y.x == edge_info$Y.y )), ]

  return(list(edges_loops, edges_nonloops))
}

# get the name of all the possible nodes
get_node_ids <- function(diff_network_object){
  return(diff_network_object[[2]]$Node)
}

# fill in blanks
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

# a function which splits edge lists by a group
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
}

# get all possible values for an edge feature
get_possible_values_edge_feature <- function(diff_network_object, edge_feature=""){
  do.call(rbind, diff_network_object$EdgeFeatures)[,edge_feature] %>% unique()
}

# get all possible values for an individual grouping
get_possible_values_individual <- function(diff_network_object, individual_feature=""){
  diff_network_object$IndividualFeatures[,individual_feature] %>% unique()
}
