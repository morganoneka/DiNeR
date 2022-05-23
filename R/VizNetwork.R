#' Plot Individual
#'
#' Create a network plot for one network
#'
#' @param diff_network_object The DiNeR object to create the plot for
#' @param idx Which index to plot
#' @param edge_weight The name of the edge feature to be used for plotting edge weight
#' @param edge_color The name of the edge feature to be used for plotting edge color
#' @param colors A dataframe of colors to be used for edge colors
#' @param nodetype A label for the node type legend
#' @return The network visualization
#' @export
plot_individual <-function(diff_network_object,idx=1,edge_weight=NULL, edge_color=NULL, colors=NULL, nodetype="Node Type"){

  node1="N1"
  node2="N2"

  edge_info <- cbind(diff_network_object$EdgeList[[idx]], diff_network_object$EdgeFeatures[[idx]])
  node_xy_df <- diff_network_object$NodeXY

  # get node1 coords
  edge_info <- merge(edge_info, node_xy_df, by.x=node1, by.y="Node", all.x=TRUE, all.y=FALSE)
  edge_info <- merge(edge_info, node_xy_df, by.x=node2, by.y="Node", all.x=TRUE, all.y=FALSE)


  if (!is.null(edge_weight)){
    names(edge_info)[names(edge_info) == edge_weight] <- 'Weight'
  } else{
    edge_info$Weight = 1
  }
  edge_info$size <- 1.5
  split_data = NULL


  if (!is.null(edge_color)){
    names(edge_info)[names(edge_info) == edge_color] <- 'Type'
    split_data <- split(edge_info,diff_network_object$EdgeFeatures[[idx]][,edge_color])

  } else{
    edge_info$Type = "Edge"
    split_data <- list(edge_info)
  }





  # split based on edge_type
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

    weight_name = paste("Weight", names(split_data)[i], sep="_")

    edges_loops[,weight_name] <- edges_loops$Weight
    edges_nonloops[,weight_name] <- edges_nonloops$Weight

    g <- g + new_scale_color() + geom_curve(
      aes_string(x = "X.x", y = "Y.x", xend = "X.y", yend = "Y.y", color=weight_name, size="size"), alpha=0.8,
      data = edges_nonloops) +
      geom_curve(
        aes_string(x = "X.x", y = "Y.x", xend = "X.y", yend = "Y.y", color=weight_name, size="size"), alpha=0.8,  curvature = 50,angle = 270,
        data = edges_loops) + scale_color_gradient(low=colors[i,"low"],high=colors[i,"high"])

    if (is.null(edge_color)){
      # g = g + guides(color = FALSE)
    }

    if (is.null(edge_weight)){
      g = g+ guides(size = FALSE)
    }
  }


  g <- g + geom_circle(aes(x0 = X, y0 = Y, r = radius, fill = Node), data=node_xy_df) +
    geom_label_repel(aes(x=X,y=Y,label=Node),hjust=0, vjust=0, data=node_xy_df) +
    guides(fill=guide_legend(title=nodetype)) + guides(size = FALSE) +
    ylim(c(0, max(node_xy_df$Y+10))) + xlim(c(0, max(node_xy_df$X+10))) +
    coord_fixed() + theme_void()



  g
}

#' Consensus Edge
#'
#' Create a network plot where edge weight is equal to the number of networks in the object that have an edge between the two nodes
#'
#' @param diff_network_object The DiNeR object to create the plot for
#' @param colors A list of colors to be used for edge colors; there should be one color for "low" and one for "high"
#' @param nodetype A label for the node type legend
#' @return The network visualization
#' @export
consensus_edge <- function(diff_network_object, colors=NULL, nodetype="Node Type"){

  node1="N1"
  node2="N2"


  node_xy_df <- diff_network_object$NodeXY

  edge_list_labels = c(node1,node2)

  # combine the edge lists from all individuals
  edge_list <- do.call(rbind,lapply(diff_network_object$EdgeList, FUN = function(x){x[,edge_list_labels]} ))

  edge_list_in_order <- edge_list[which(edge_list[,node1] <= edge_list[,node2]),]
  edge_list_out_of_order <- edge_list[which(edge_list[,node1] > edge_list[,node2]),]
  colnames(edge_list_out_of_order) <- c(node2,node1)


  occurrences <- melt(table(rbind(edge_list_in_order, edge_list_out_of_order)))
  occurrences <- occurrences[which(occurrences$value >0),]

  occurrences <- merge(occurrences, node_xy_df, by.x=node1, by.y="Node", all.x=TRUE, all.y=FALSE)
  occurrences <- merge(occurrences, node_xy_df, by.x=node2, by.y="Node", all.x=TRUE, all.y=FALSE)
  occurrences$size = 1.5

  splits_occurrences <- split_loops(occurrences)

  edges_loops_occurrences <- splits_occurrences[[1]]
  edges_nonloops_occurrences <- splits_occurrences[[2]]


  if(is.null(colors)){
    colors = sample(rainbow(2))
  }

  ggplot() +
    geom_curve(
      aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, size=size, color=abs(value), alpha=0.8),
      data = edges_nonloops_occurrences
    ) +
    geom_curve(
      aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, size=size, color=abs(value)),alpha=0.8,  curvature = 50,angle = 270,
      data = edges_loops_occurrences
    ) + scale_color_gradient(low=colors[1],high=colors[2]) +
    guides(size=FALSE) + guides(alpha=FALSE) +
    guides(color=guide_colorbar(title="Number of Individuals with Edge")) +
    geom_circle(aes(x0 = X, y0 = Y, r = radius, fill = Node), data=node_xy_df) +geom_label_repel(aes(x=X,y=Y,label=Node),hjust=0, vjust=0, data=node_xy_df) +
    guides(fill=guide_legend(title=nodetype)) +
    ylim(c(0, max(node_xy_df$Y+10))) + xlim(c(0, max(node_xy_df$X+10))) +
    coord_fixed() + theme_void()
}

#TODO: unique edges

subset_network_plot <- function(network_plot, nodename){
  q <- ggplot_build(network_plot)

  # get index for df within ggplot_build that works with nodes (last one)
  node_idx = length(q$data)

  # identify location of this cell
  celltype_loc = q$data[[node_idx]][which(q$data[[node_idx]]$label == nodename),c("y","x")]

  for (i in 1:(length(q$data)-1)){
    # "yend" is only a column when we are dealing with lines
    if ("yend" %in% colnames(q$data[[i]])){
      relevant_lines = c(which(q$data[[i]]$y == celltype_loc$y & q$data[[i]]$x == celltype_loc$x),
                         which(q$data[[i]]$yend == celltype_loc$y & q$data[[i]]$xend == celltype_loc$x))
      color_title = colnames(q$data[[i]])[grep("colour", colnames(q$data[[i]]))[1]]
      q$data[[i]][setdiff(1:nrow(q$data[[i]]),relevant_lines),color_title] = "#D3D3D3"
      q$data[[i]]$alpha = 1
      q$data[[i]][setdiff(1:nrow(q$data[[i]]),relevant_lines),"alpha"] = 0.2
    }

  }


  q <- ggplot_gtable(q)
  grid::grid.draw(q)
}
