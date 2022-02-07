library(igraph)
library(ggplot2)
library(ggforce)
library(reshape2)
library(ggrepel)
library(ggnewscale)
library(scatterpie)
library(ggpubr)
library(ggthemes)

############################ DIFF NETWORK OBJECT ######################
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
    
    sub_dno$NodeXY[,"Index"] = 1:length(indices)
    
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

#TODO: allow calculation for consensus edge for ALL vs. just a subset 
# (i.e. consensus for certain indices, or patients with a certain label)
consensus_edge <- function(diff_network_object, edge_weight="Weight", edge_color="Type"){
 
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
  
  splits_occurrences <- split_loops(occurrences)
  
  edges_loops_occurrences <- splits_occurrences[[1]]
  edges_nonloops_occurrences <- splits_occurrences[[2]]
  
  ggplot() + 
    geom_curve(
      aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, size=abs(value), color=abs(value)),
      data = edges_nonloops_occurrences
    ) +
    geom_curve(
      aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, size=abs(value), color=abs(value)),  curvature = 50,angle = 270,
      data = edges_loops_occurrences
    ) +
    guides(size=FALSE) +
    guides(color=guide_colorbar(title="Number of Individuals with Edge")) +
    geom_circle(aes(x0 = X, y0 = Y, r = radius, fill = Node), data=node_xy_df) +geom_label_repel(aes(x=X,y=Y,label=Node),hjust=0, vjust=0, data=node_xy_df) +
    guides(fill=guide_legend(title="Cell Type")) + 
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
    relevant_lines = c(which(q$data[[i]]$y == celltype_loc$y & q$data[[i]]$x == celltype_loc$x),
                       which(q$data[[i]]$yend == celltype_loc$y & q$data[[i]]$xend == celltype_loc$x))
    color_title = colnames(q$data[[i]])[grep("colour", colnames(q$data[[i]]))[1]]
    q$data[[i]][setdiff(1:nrow(q$data[[i]]),relevant_lines),color_title] = "#D3D3D3"
    q$data[[i]]$alpha = 1
    q$data[[i]][setdiff(1:nrow(q$data[[i]]),relevant_lines),"alpha"] = 0.2
  }
  

  q <- ggplot_gtable(q)
  grid::grid.draw(q)
}


############################ PIE CHART VISUALIZATION ######################
make_piecharts <- function(diff_network_object, colors=c(), edge_feature=""){
  
  #TODO test to see if edge_Feature is in columns for diff_network_object$EdgeFeature
  
  node1="N1"
  node2="N2"
  
  node_xy_df <- diff_network_object$NodeXY
  
  edge_list_labels = c(node1,node2)
  
  # for each individual, combine edge list with edge features
  edge_list <- do.call(rbind,lapply(1:length(diff_network_object$EdgeList), FUN = function(x){
    cbind(diff_network_object$EdgeList[[x]][,c(node1,node2)], diff_network_object$EdgeFeatures[[x]][,edge_feature])
  } ))
  
  # combine in-order and reverse-order edge lists
  edge_list_in_order <- edge_list[which(edge_list[,node1] <= edge_list[,node2]),]
  edge_list_out_of_order <- edge_list[which(edge_list[,node1] > edge_list[,node2]),]
  colnames(edge_list_out_of_order) <- c(node2,node1,  edge_feature)
  colnames(edge_list_in_order) <- c(node1,node2,  edge_feature)
  
  # get counts for each node-node pair and possible edge class
  occurrences <- melt(table(rbind(edge_list_in_order, edge_list_out_of_order)))
  
  # create new coordinate system for pie chart
  celltypes <- unlist(diff_network_object$NodeXY$Node)
  coords <- as.data.frame(cbind(celltypes, 1:length(celltypes)), stringsAsFactors = FALSE)
  coords$V2 = as.numeric(coords$V2)
  
  # combine edge counts with coordinate system
  combo = merge(occurrences, coords, by.x=node1, "celltypes", all=FALSE)
  combo = merge(combo, coords, by.x=node2, "celltypes", all=FALSE)
  colnames(combo) <- c(node1,node2,  edge_feature, "value", "X", "Y")
  combo_flipped  = merge(occurrences, coords, by.x=node2, "celltypes", all=FALSE)
  combo_flipped = merge(combo_flipped, coords, by.x=node1, "celltypes", all=FALSE)
  colnames(combo_flipped) <- c(node1,node2, edge_feature, "value", "X", "Y")
  combo = rbind(combo, combo_flipped)
  combo = combo[which(combo$value >0),]
  
  # how many individuals do we have
  num_patients = length(diff_network_object[[1]])
  
  # convert from factor to character 
  combo[,node1] = as.character(combo[,node1])
  combo[,node2] = as.character(combo[,node2])
  
  # identify all possible interaction types
  interaction_types = unique(unlist(lapply(diff_network_object$EdgeFeatures, "[[", edge_feature)))
  
  # adds rows for missing pairs (i.e. if there isn't a depleted edge between A and C, adds a row that's like A C depleted 0)
  for (cell1 in celltypes){
    for (cell2 in celltypes){
      for (inxtype in c(interaction_types)){

        
        w = which(combo[,node1] == cell1 & combo[,node2] == cell2 & combo[,edge_feature] == inxtype)
        # print(cell1,cell2,inxtype)
        # print(w)
        if (length(w) == 0){
          # print(cell1,cell2,inxtype)
          # print("UGH")
          combo = rbind(combo, c(cell1,cell2,inxtype,0,coords[which(coords$celltypes == cell1), "V2"],coords[which(coords$celltypes == cell2), "V2"]), stringsAsFactors=FALSE)
        }
      }
    }
  }
  
  #TODO: fix to not need this later
  combo = combo[!duplicated(combo[,c(node1,node2,edge_feature)]),]
  
  # convert to the correct data types
  combo$X = as.numeric(combo$X)
  combo$Y = as.numeric(combo$Y)
  combo[,edge_feature] = as.character(combo[,edge_feature])
  
  # convert value to percentage
  combo$value = as.numeric(combo$value)
  combo$value = combo$value / num_patients
  
  # combo_clean <- dcast(unique(combo), node1 + node2 ~ edge_feature, value.var="value")
  # TODO: dcast is giving an issue so see what's up
  combo_clean <- dcast(unique(combo), as.formula(paste(paste(node1, node2, sep="+"), "~", edge_feature)), value.var="value")
  combo_clean = merge(combo_clean, coords, by.x=node1, "celltypes", all=FALSE)
  combo_clean = merge(combo_clean, coords, by.x=node2, "celltypes", all=FALSE)
  colnames(combo_clean) <- c(node1,node2, interaction_types, "X", "Y")
  combo_clean$Neither = 1 - rowSums(combo_clean[,interaction_types])/num_patients
  combo_clean$radius = 0.25
  
  if (length(colors) != (length(interaction_types)+1)){
    colors = rainbow(length(interaction_types)+1)
  }
  
  
  ggplot() + geom_scatterpie(aes_string(x="X", y="Y", r="radius"), data=combo_clean, cols=c(interaction_types, "Neither")) + coord_equal() +
    scale_x_discrete(limits = as.character(0:(length(celltypes)-1) + 0.5), breaks=as.character(0:(length(celltypes)-1) + 0.5), labels=celltypes) +
    scale_y_discrete(limits = as.character(0:(length(celltypes)-1) + 0.5), breaks=as.character(0:(length(celltypes)-1) + 0.5), labels=celltypes) + 
    theme_pubr() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7), axis.text.y = element_text(size=7), axis.title.x = element_blank(), axis.title.y = element_blank())  + 
    scale_fill_manual(values=colors) + labs(fill = "Interaction Type")
  
}

make_piecharts_by_group <- function(diff_network_object, colors=c(), edge_feature="", group_label=""){
  
  node1 = "N1"
  node2 = "N2"
  
  # split into groups 
  split_dnos = split_by_group(diff_network_object, group_label)
  
  # get possible values for groups 
  #TODO make sure this is always in the same order?
  groups = get_possible_values_individual(diff_network_object, group_label)
  
  # get possible edge feature values
  edge_features = get_possible_values_edge_feature(diff_network_object, edge_feature)
  
  tables = lapply(split_dnos, function(x){
    # get one dataframe with all edges and their features
    g1 = cbind(do.call(rbind,x$EdgeList), do.call(rbind,x$EdgeFeatures))
    
    # create a column for each interaction
    g1$Interaction = unlist(lapply(1:nrow(g1), function(i) paste(g1[i,c(node1,node2)],collapse=" / ")))
    
    # convert to distribution table
    dist_table = as.data.frame(table(g1[,c( "Interaction",edge_feature)])) %>% spread(edge_feature, Freq)
    
    # number of "other" individuals: 
    dist_table$Other = length(x$EdgeList) - rowSums(dist_table[,setdiff(colnames(dist_table), "Interaction")])
    return(dist_table)
  })
  
  get_row <- function(row){
    return(lapply(1:length(tables), function(x){
      tmp = tables[[x]]
      idx = which(tmp$Interaction == row)
      if (length(idx >0)){
        return(tmp[idx,])
      } else{
        return(c(row, 0,0,length(split_dnos[[x]]$EdgeList)))
      }

    }))
  }

  
  all_inx <- as.character(unique(do.call(rbind,tables)$Interaction))
  
  # create matrix to store distances
  #TODO: this will only happen if flag says we sort by distance vs. 
  inx_dist <- as.data.frame(matrix(nrow=0,ncol=2), stringsAsFactors = FALSE)
  
  # scatterpie coords
  scatter_coords <- as.data.frame(matrix(nrow=0,ncol=6), stringsAsFactors = FALSE)
  
  for (idx in 1:length(all_inx)){
    # get interaction counts for all groups
    inx_numbers = do.call(rbind,get_row(all_inx[idx])) %>% as.data.frame()
    colnames(inx_numbers) = c("Interaction", edge_features, "Other")
    inx_numbers$Group = groups
    inx_numbers$X = 1:nrow(inx_numbers)
    scatter_coords = rbind(scatter_coords, inx_numbers)
    
    # calculate all of the probability vectors
    prob_vec = lapply(1:nrow(inx_numbers), function(x){
      vec = inx_numbers[x, c(edge_features, "Other")]
      return(as.numeric(vec) / sum(as.numeric(vec)))
    })
    
    # calculate all the pairwise distances
    dists = lapply(1:(nrow(inx_numbers)-1), function(x){
      unlist(lapply((x+1):nrow(inx_numbers), function(y){
        return( sqrt(sum((prob_vec[[x]] - prob_vec[[y]])^2)) )
      }))
    })
    
    dist_sum = sum(unlist(dists))
    inx_dist = rbind(inx_dist, c(all_inx[[idx]], dist_sum), stringsAsFactors=FALSE)
    
  }
  
  inx_order = cbind(inx_dist[order(inx_dist[,2]),], 1:nrow(inx_dist))
  colnames(inx_order) <- c("Interaction", "Distance", "Y")
  
  merged_data = merge(scatter_coords,inx_order, by="Interaction")
  merged_data$Radius = 0.4
  
  merged_data[,c(edge_features, "Other")] = sapply(merged_data[,c(edge_features, "Other")], as.numeric)
  
  if (length(colors) != (length(interaction_types)+1)){
    colors = rainbow(length(interaction_types)+1)
  }
  
  ggplot() + geom_scatterpie(aes(x=X, y=Y, r=Radius), data=merged_data, cols=c(edge_features, "Other"))  + coord_equal() +
    scale_x_discrete(limits = as.character(0:(length(groups)-1) + 0.5), breaks=as.character(0:(length(groups)-1) + 0.5), labels=groups, position="top") +
    scale_y_discrete(limits = as.character(0:(length(all_inx)-1) + 0.5), breaks=as.character(0:(length(all_inx)-1) + 0.5), labels=unique(merged_data[order(merged_data$Y, decreasing = FALSE), "Interaction"]))+ 
    theme_pubr() + theme(axis.text.x = element_text(angle = 45, hjust = 0, size=9), axis.text.y = element_text(size=9), axis.title.x = element_blank(), axis.title.y = element_blank())  + 
    scale_fill_manual(values=colors) + labs(fill = "Interaction Type")
}

############################ IGRAPH STUFF
add_graphs <- function(diff_network_object, edge_weight=""){
  
  graphs = lapply(1:length(diff_network_object$EdgeList), function(i){
    
    g = graph_from_edgelist(diff_network_object$EdgeList[[i]] %>% as.matrix(), directed=FALSE)
    
    if (edge_weight != ""){
      E(g)$weight = diff_network_object$EdgeFeatures[[i]][,edge_weight]
    } 
    
    return(g)
  })
  
  diff_network_object$Graphs = graphs
  
  return(diff_network_object)
}
