library(igraph)
library(ggplot2)
library(ggforce)
library(reshape2)
library(ggrepel)
library(ggnewscale)
library(scatterpie)
library(ggpubr)

############################ DIFF NETWORK OBJECT ######################
create_diff_network_object <- function(edge_list){

  
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
  
  # ggplot() + geom_circle(aes(x0 = X, y0 = Y, r = radius, fill = Node), data=node_xy_df)+
    # coord_fixed()
  
  return(list(
    "EdgeList" = edge_list,
    "NodeXY" = node_xy_df
    ))
}

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

add_edge_features <- function(diff_network_object, df_list){
  if (length(df_list) != length(diff_network_object$EdgeList)){
    stop("DF lists not the same size")
  }
  
  # if we already have features, combine
  if ("EdgeFeatures" %in% names(diff_network_object)){
    diff_network_object$EdgeFeatures = lapply(1:length(df_list), function(i){
      cbind(diff_network_object$EdgeFeatures[[i]], df_list[[i]])
    })
  } else { # otherwise....
    diff_network_object$EdgeFeatures = df_list
  }
  
  return (diff_network_object)
  
  
}

add_node_features <- function(diff_network_object){
  
}

############################ HELPER FUNCTIONS ######################
 
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

############################ NETWORK VISUALIZATION ######################
plot_individual <-function(diff_network_object,idx=1,edge_weight="Weight", edge_color="Type", colors=NULL){
  
  node1="N1"
  node2="N2"
  
  edge_info <- cbind(diff_network_object$EdgeList[[idx]], diff_network_object$EdgeFeatures[[idx]])
  node_xy_df <- diff_network_object$NodeXY
  
  # get node1 coords
  edge_info <- merge(edge_info, node_xy_df, by.x=node1, by.y="Node", all.x=TRUE, all.y=FALSE)
  edge_info <- merge(edge_info, node_xy_df, by.x=node2, by.y="Node", all.x=TRUE, all.y=FALSE)
  
  names(edge_info)[names(edge_info) == edge_weight] <- 'Weight'
  names(edge_info)[names(edge_info) == edge_color] <- 'Type'
  
  
  edge_info$size <- 1.5
  
  
  # split based on edge_type
  split_data <- split(edge_info,diff_network_object$EdgeFeatures[[idx]][,edge_color])
  g <- ggplot()
  
  if(is.null(colors)){
    colors = as.data.frame(matrix(sample((rainbow(length(split_data)*2))), ncol=2, nrow=length(split_data)))
    colnames(colors) = c("high", "low")
  }
  
  for (i in 1:length(split_data)){
    # logic for self loops
    splits <- split_loops(split_data[[i]])
    
    edges_loops <- splits[[1]]
    edges_nonloops <- splits[[2]]
    
    weight_name = paste(names(split_data)[i], "Weight", sep="_")
    
    edges_loops[,weight_name] <- edges_loops$Weight
    edges_nonloops[,weight_name] <- edges_nonloops$Weight
    
    g <- g + new_scale_color() + geom_curve(
      aes_string(x = "X.x", y = "Y.x", xend = "X.y", yend = "Y.y", color=weight_name, size="size"), alpha=0.8,
      data = edges_nonloops) +
      geom_curve(
        aes_string(x = "X.x", y = "Y.x", xend = "X.y", yend = "Y.y", color=weight_name, size="size"), alpha=0.8,  curvature = 50,angle = 270,
        data = edges_loops) + scale_color_gradient(low=colors[i,"low"],high=colors[i,"high"]) 
  }
  
  
  g <- g + geom_circle(aes(x0 = X, y0 = Y, r = radius, fill = Node), data=node_xy_df) + 
    geom_label_repel(aes(x=X,y=Y,label=Node),hjust=0, vjust=0, data=node_xy_df) +
    guides(fill=guide_legend(title="Cell Type")) + guides(size = FALSE) +
    ylim(c(0, max(node_xy_df$Y+10))) + xlim(c(0, max(node_xy_df$X+10))) +
    coord_fixed() + theme_void()
  
  g
}

#TODO: consensus edge

#TODO: unique edges