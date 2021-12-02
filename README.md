# DiNeR (Differential Network Visualization in R)

This is code to quantify and visualize differences in networks that have a shared node set. This was designed with biological networks, such as genetic coexpression networks, in mind.

## Data pre-processing
Requires a list of `data.frame` objects representing the edge list for a graph. Each edge can also be given weights or types.

```r
test_df_1 <- data.frame(E1 = c("A","B","C"), E2 = c("A","A","F"), Weight=c(2.78, -5.12, -1.28), Type=c("Enriched", "Depleted", "Depleted"), stringsAsFactors=FALSE)
test_df_2 <- data.frame(E1 = c("E","A","A","H","C","F"), E2 = c("D","B","F","G","A","F"), Weight=c(2.78, -5.12, 1.28, 6.7, 1.9, -8.3), Type=c("Enriched", "Depleted", "Enriched", "Enriched", "Enriched", "Depleted"), stringsAsFactors=FALSE)
test_df_3 <- data.frame(E1 = c("A","B","F","C"), E2 = c("G","A","B","D"), Weight=c(1,2,3,4), Type=c("Enriched", "Enriched", "Enriched", "Enriched"), stringsAsFactors=FALSE)

df_list <- list(test_df_1, test_df_2, test_df_3)
```

## How to run
First, run `create_diff_network_object`. Then, this object is passed to any subsequent function calls.

```r
diff_network_object <- create_diff_network_object(df_list)
```

There are several parameters to pass to this function
- **node1** and **node2** specify the name of the two columns containing the node IDs.
- **weight** specifies the weight or strength of the edge.
- **type** specifies the type of edge.

## Individual network plot
The `plot_individual` function plots an individual network.

```r
plot_patient(diff_network_object)
```

This can be run on all networks using `lapply`.

```r
lapply(1:length(df_list), function(x){plot_patient(diff_network_object, x)})
```

For each individual network, the nodes will be in the same location, allowing for easier comparison of individual networks.

## Consensus edge network plot
The `consensus_edge` function creates one plot identifying the number of networks that share any given edge.


## Differential edge heatmap
The `diff_edge` function creates one heatmap showing how many networks share a given edge. This is an alternative to `consensus_edge` that is perhaps a bit more practical for large networks.




## Future functionality
- Integrating node information
- Additional layouts for network plots
- Color scheme customization 
