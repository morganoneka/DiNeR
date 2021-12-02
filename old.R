library(igraph)
library(ggplot2)
library(ggforce)
library(reshape2)
library(ggrepel)
library(ggnewscale)
library(scatterpie)
library(ggpubr)

### DIFF NETWORK OBJECT
create_diff_network_object <- function(df_list, node1="N1", node2="N2", edge_weight="Weight", edge_color="Type"){
  edge_weight = edge_weight
  edge_color = edge_color
  edge_list_labels = c(node1,node2)
  
  edge_vars = list(
    node1 = node1,
    node2 = node2,
    edge_weight = edge_weight,
    edge_color = edge_color
  )
  
  # identify all distinct node types
  node_labels <- unique(c(unlist(lapply(df_list, "[[", edge_list_labels[1])), 
                          unlist(lapply(df_list, "[[", edge_list_labels[2]))))
  
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
  
  return(list(df_list,node_xy_df, edge_vars))
}

#TODO
add_node_features <- function(){
  
}

#TODO
add_sample_features <- function(){
  
}

### HELPER FUNCTIONS
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


### MAIN FUNCTIONS
plot_individual <- function(diff_network_object,idx=1, node1="N1", node2="N2", edge_weight="Weight", edge_color="Type"){
  
  edge_info <- diff_network_object[[1]][[idx]]
  node_xy_df <- diff_network_object[[2]]
  
  # get node1 coords
  edge_info <- merge(edge_info, node_xy_df, by.x=node1, by.y="Node", all.x=TRUE, all.y=FALSE)
  edge_info <- merge(edge_info, node_xy_df, by.x=node2, by.y="Node", all.x=TRUE, all.y=FALSE)
  
  names(edge_info)[names(edge_info) == edge_weight] <- 'Weight'
  names(edge_info)[names(edge_info) == edge_color] <- 'Type'
  
  # logic for self loops
  splits <- split_loops(edge_info)
  
  edges_loops <- splits[[1]]
  edges_nonloops <- splits[[2]]
  
  
  ggplot() + 
    geom_curve(
      aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, color=Type, size=abs(Weight)),
      data = edges_nonloops, alpha=0.8
    ) +
    geom_curve(
      aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, color=Type, size=abs(Weight)),  curvature = 50,angle = 270,
      data = edges_loops, alpha=0.8
    ) +
    geom_circle(aes(x0 = X, y0 = Y, r = radius, fill = Node), data=node_xy_df) +geom_label_repel(aes(x=X,y=Y,label=Node),hjust=0, vjust=0, data=node_xy_df) +
    guides(fill=guide_legend(title="Cell Type")) +
    guides(size=guide_legend(title="Interaction Strength")) + 
    ylim(c(0, max(node_xy_df$Y+10))) + xlim(c(0, max(node_xy_df$X+10))) +
    coord_fixed() + theme_void()
}

### TODO: get working
plot_individual_alt <-function(diff_network_object,idx=1, node1="N1", node2="N2", edge_weight="Weight", edge_color="Type"){
  
  edge_info <- diff_network_object[[1]][[idx]]
  node_xy_df <- diff_network_object[[2]]
  
  # get node1 coords
  edge_info <- merge(edge_info, node_xy_df, by.x=node1, by.y="Node", all.x=TRUE, all.y=FALSE)
  edge_info <- merge(edge_info, node_xy_df, by.x=node2, by.y="Node", all.x=TRUE, all.y=FALSE)
  
  names(edge_info)[names(edge_info) == edge_weight] <- 'Weight'
  names(edge_info)[names(edge_info) == edge_color] <- 'Type'
  
  
  edge_info$size <- 1.5
  
  
  # split based on edge_type
  split_data <- split(edge_info,edge_info$Type)
  g <- ggplot()
  
  # TODO: generate colors based on different groups - need high and low value
  colors <- data.frame(high=c("red", "gray50"), low=c("blue", "gray90"))
  for (i in 1:length(split_data)){
    # logic for self loops
    splits <- split_loops(split_data[[i]])
    
    edges_loops <- splits[[1]]
    edges_nonloops <- splits[[2]]
    
    # g = g 
    # + scale_color_gradient(high="#ff0000", low="#0000ff")
    # scale_color_gradient(high=colors[i,"high"],low=colors[i,"low"])
    
    # TODO: need to generalize 
    # thanks to: https://github.com/eliocamp/ggnewscale 
    if (i == 1){
      print("i = 1")
      edges_loops[,"Depleted Weight"] <- edges_loops$Weight
      edges_nonloops[,"Depleted Weight"] <- edges_nonloops$Weight
      g = g + geom_curve(
        aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, color=`Depleted Weight`, size=size), alpha=0.8,
        data = edges_nonloops) +
        geom_curve(
          aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, color=`Depleted Weight`, size=size), alpha=0.8,  curvature = 50,angle = 270,
          data = edges_loops) + scale_color_gradient(low="#191970",high="#ADD8E6") 
    } else{
      print("i = 2")
      edges_loops[,"Enriched Weight"] <- edges_loops$Weight
      edges_nonloops[,"Enriched Weight"] <- edges_nonloops$Weight
      g = g + new_scale_color() + geom_curve(
        aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, color=`Enriched Weight`, size=size), alpha=0.8,
        data = edges_nonloops) +
        geom_curve(
          aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, color=`Enriched Weight`, size=size), alpha=0.8,  curvature = 50,angle = 270,
          data = edges_loops)   + scale_color_gradient(high="#9A2A2A",low="#FAA0A0") 
    }
  }
  
  
  g = g + geom_circle(aes(x0 = X, y0 = Y, r = radius, fill = Node), data=node_xy_df) + geom_label_repel(aes(x=X,y=Y,label=Node),hjust=0, vjust=0, data=node_xy_df) +
    guides(fill=guide_legend(title="Cell Type")) + guides(size = FALSE) +
    ylim(c(0, max(node_xy_df$Y+10))) + xlim(c(0, max(node_xy_df$X+10))) +
    coord_fixed() + theme_void()
  
  g
}

consensus_edge <- function(diff_network_object, node1="N1", node2="N2", edge_weight="Weight", edge_color="Type"){
  
  # adj_mx <- table(do.call(rbind,lapply(df_list, FUN = function(x){x[,edge_list_labels]} )))
  # adj_mx_2 <- table(do.call(rbind,lapply(df_list, FUN = function(x){x[,c(edge_list_labels[2], edge_list_labels)]} )))
  # melt(adj_mx)
  
  node_xy_df <- diff_network_object[[2]]
  
  edge_list_labels = c(node1,node2)
  
  edge_list <- do.call(rbind,lapply(diff_network_object[[1]], FUN = function(x){x[,edge_list_labels]} ))
  
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
    guides(size=guide_legend(title="Number of Patients w/ Interaction")) +
    guides(color=guide_legend(title="Number of Patients w/ Interaction")) +
    geom_circle(aes(x0 = X, y0 = Y, r = radius, fill = Node), data=node_xy_df) +geom_label_repel(aes(x=X,y=Y,label=Node),hjust=0, vjust=0, data=node_xy_df) +
    guides(fill=guide_legend(title="Cell Type")) + 
    ylim(c(0, max(node_xy_df$Y+10))) + xlim(c(0, max(node_xy_df$X+10))) +
    coord_fixed() + theme_void()
} 

### UNFINISHED
#TODO: i think this is about done but double check
diff_edge <- function(diff_network_object){
  
  node_ids <- get_node_ids(diff_network_object)
  edge_vars = diff_network_object[[3]]
  
  
  mx_diff <- data.frame(matrix(nrow=0, ncol=3))
  for (i in 1:n_row_col){
    for (j in 1:n_row_col){
      # network similarity
      edge_list_1 <- sort_list(fill_missing_vals(diff_network_object[[1]][[i]], node_ids, edge_vars),edge_vars)
      edge_list_2 <- sort_list(fill_missing_vals(diff_network_object[[1]][[j]], node_ids, edge_vars),edge_vars)
      
      mx_diff <- rbind(mx_diff,c(i,j,sum(sum(abs(as.numeric(edge_list_1[,edge_vars$edge_weight]) - as.numeric(edge_list_2[,edge_vars$edge_weight]))))))
    }
  }
  
  colnames(mx_diff) <- c("RegionA", "RegionB", "Difference")
  ggplot(mx_diff, aes(RegionA,RegionB, fill=Difference)) + geom_raster() 
}

#TODO: finish me
scale_edge_weight <- function(diff_network_object){
  # identify the maximum edge weight (maximum absolute value)
  weights <- lapply(diff_network_object[[1]], "[[", diff_network_object[[3]]$edge_weight)
  heaviest <- max(abs(unlist(weights)))
  
  weighted_mx <- lapply(diff_network_object[[1]], function(x, mx, label){
    x[,label] = x[,label] / heaviest
    return(x)
  }, heaviest, diff_network_object[[3]]$edge_weight)
  
  diff_network_object[[1]] = weighted_mx
  return(diff_network_object)
}

#TODO: finish me
calc_hub_score <- function(diff_network_object){
  x <- lapply(diff_network_object[[1]], fill_missing_vals, diff_network_object[[3]], get_node_ids(diff_network_object))
}

plot_hub_score_heatmap <- function(diff_network_object){
  # lapply(diff_network_object[[1]], function(x){
  #  return(x[,c(diff_network_object[[3]]$node1, diff_network_object[[3]]$node2, diff_network_object[[3]]$edge_weight)]) 
  # })
  # adj_mxs <- as.matrix(get.adjacency(graph.data.frame(edges)))
  # graphs <- lapply(, graph_from_adjacency_matrix, weighted=T)
  
}

# TODO: fix default colors
make_piecharts <- function(diff_network_object, colors=c()){
  
  node_xy_df <- diff_network_object[[2]]
  edge_vars = diff_network_object[[3]]
  
  edge_list_labels = c(edge_vars[["node1"]],edge_vars[["node2"]])
  
  edge_list <- do.call(rbind,lapply(diff_network_object[[1]], FUN = function(x){x[,c(edge_list_labels, edge_vars[["edge_color"]])]} ))
  
  edge_list_in_order <- edge_list[which(edge_list[,edge_vars[["node1"]]] <= edge_list[,edge_vars[["node2"]]]),]
  edge_list_out_of_order <- edge_list[which(edge_list[,edge_vars[["node1"]]] > edge_list[,edge_vars[["node2"]]]),]
  colnames(edge_list_out_of_order) <- c(edge_vars[["node2"]],edge_vars[["node1"]],  edge_vars[["edge_color"]])
  
  occurrences <- melt(table(rbind(edge_list_in_order, edge_list_out_of_order)))
  
  celltypes <- unlist(diff_network_object[[2]]$Node)
  coords <- as.data.frame(cbind(celltypes, 1:length(celltypes)), stringsAsFactors = FALSE)
  coords$V2 = as.numeric(coords$V2)
  
  combo = merge(occurrences, coords, by.x=edge_vars[["node1"]], "celltypes", all=FALSE)
  combo = merge(combo, coords, by.x=edge_vars[["node2"]], "celltypes", all=FALSE)
  
  colnames(combo) <- c(edge_vars[["node1"]],edge_vars[["node2"]],  edge_vars[["edge_color"]], "value", "X", "Y")
  # combo[is.na(combo$value),"value"] = 0
  
  combo_flipped  = merge(occurrences, coords, by.x=edge_vars[["node2"]], "celltypes", all=FALSE)
  combo_flipped = merge(combo_flipped, coords, by.x=edge_vars[["node1"]], "celltypes", all=FALSE)
  
  colnames(combo_flipped) <- c(edge_vars[["node1"]],edge_vars[["node2"]], edge_vars[["edge_color"]], "value", "X", "Y")
  
  combo = rbind(combo, combo_flipped)
  combo = combo[which(combo$value >0),]
  
  num_patients = length(diff_network_object[[1]])
  
  combo[,edge_vars[["node1"]]] = as.character(combo[,edge_vars[["node1"]]])
  combo[,edge_vars[["node2"]]] = as.character(combo[,edge_vars[["node2"]]])
  
  interaction_types = unique(unlist(lapply(diff_network_object[[1]], "[[", edge_vars[["edge_color"]])))
  
  for (cell1 in celltypes){
    for (cell2 in celltypes){
      for (inxtype in c(interaction_types)){
        
        w = which(combo$cell_1 == cell1 & combo$cell_2 == cell2 & combo$InxType == inxtype)
        # print(w)
        if (length(w) == 0){
          combo = rbind(combo, c(cell1,cell2,inxtype,0,coords[which(coords$celltypes == cell1), "V2"],coords[which(coords$celltypes == cell2), "V2"]), stringsAsFactors=FALSE)
        }
      }
    }
  }
  
  combo$X = as.numeric(combo$X)
  combo$Y = as.numeric(combo$Y)
  combo[,edge_vars[["edge_color"]]] = as.character(combo[,edge_vars[["edge_color"]]])
  combo$value = as.numeric(combo$value)
  combo$value = combo$value / num_patients
  
  # combo_clean <- dcast(unique(combo), edge_vars[["node1"]] + edge_vars[["node2"]] ~ edge_vars[["edge_color"]], value.var="value")
  # TODO: dcast is giving an issue so see what's up
  combo_clean <- dcast(unique(combo), as.formula(paste(paste(edge_vars[["node1"]], edge_vars[["node2"]], sep="+"), "~", edge_vars[["edge_color"]])), value.var="value")
  combo_clean = merge(combo_clean, coords, by.x=edge_vars[["node1"]], "celltypes", all=FALSE)
  combo_clean = merge(combo_clean, coords, by.x=edge_vars[["node2"]], "celltypes", all=FALSE)
  colnames(combo_clean) <- c("cell_1", "cell_2", interaction_types, "X", "Y")
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

# TODO: make generalized
make_heatmap <- function(x){
  # heatmap stuff
  
  # combo2 <- combo
  # combo2$X = combo$X - 1
  # 
  # combo3 <- combo 
  # combo3$Y = combo$Y + 1
  # 
  # combo4 <- combo
  # combo4$X = combo$X - 1
  # combo4$Y = combo$Y + 1
  # 
  # 
  # lower_tri = rbind(rbind(combo,combo2),combo3)
  # lower_tri$ID = paste(lower_tri$cell_1, lower_tri$cell_2, sep= " ")
  # lower_tri = lower_tri[which(lower_tri$InxType == "Depletion"),]
  # 
  # 
  # upper_tri = rbind(rbind(combo4,combo3), combo2)
  # upper_tri$ID = paste(upper_tri$cell_1, upper_tri$cell_2, sep= " ")
  # upper_tri = upper_tri[which(upper_tri$InxType == "Enrichment"),]
  # 
  # 
  # ggplot() + geom_polygon(data=lower_tri, aes(x=X, y=Y, fill=value, group=ID)) + scale_fill_gradient(high="#191970",low="#ADD8E6") + 
  #   new_scale_fill() + geom_polygon(data=upper_tri, aes(x=X, y=Y, fill=value, group=ID))  + scale_fill_gradient(high="#9A2A2A",low="#FAA0A0") +
  #   scale_x_discrete(limits = as.character(0:(length(celltypes)-1) + 0.5), breaks=as.character(0:(length(celltypes)-1) + 0.5), labels=celltypes) +
  #   scale_y_discrete(limits = as.character(0:(length(celltypes)-1) + 0.5), breaks=as.character(0:(length(celltypes)-1) + 0.5), labels=celltypes) + 
  #   coord_fixed() + theme_pubr() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  # 
  # 
}


