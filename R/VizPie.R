
#' Make pie charts
#'
#' Create symmetric pie chart plot for edge class
#'
#' @param diff_network_object The DiNeR object to create the plot for
#' @param edge_feature The name of the feature we want to get edge class composition for
#' @return The pie chart visualization
#' @export
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

#' Make pie charts
#'
#' Create piecharts of edge class for each edge, split by a group
#'
#' @param diff_network_object The DiNeR object to create the plot for
#' @param colors The colors to use. This list must have as many colors as there are possible values for the edge feature, plus an additional color for "Other"
#' @param edge_feature The name of the feature we want to get edge class composition for
#' @param rows OPTIONAL: The names of the interactions we want to visualize. Each entry is styled as: "Node1 / Node2".
#' @return The pie chart visualization
#' @export
make_piecharts_by_group <- function(diff_network_object, colors=c(), edge_feature="", group_label="", rows=0){

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
  if (rows != 0){
    all_inx = all_inx[all_inx %in% rows]
  }

  # create matrix to store distances
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


  if (length(colors) != (length(edge_features)+1)){
    colors = rainbow(length(edge_features)+1)
  }


  ggplot() + geom_scatterpie(aes(x=X, y=Y, r=Radius), data=merged_data, cols=c(edge_features, "Other"))  + coord_equal() +
    scale_x_discrete(limits = as.character(0:(length(groups)-1) + 0.5), breaks=as.character(0:(length(groups)-1) + 0.5), labels=groups, position="top") +
    scale_y_discrete(limits = as.character(0:(length(all_inx)-1) + 0.5), breaks=as.character(0:(length(all_inx)-1) + 0.5), labels=unique(merged_data[order(merged_data$Y, decreasing = FALSE), "Interaction"]))+
    theme_pubr() + theme(axis.text.x = element_text(angle = 45, hjust = 0, size=9), axis.text.y = element_text(size=9), axis.title.x = element_blank(), axis.title.y = element_blank())  +
    scale_fill_manual(values=colors) + labs(fill = "Interaction Type")



}
