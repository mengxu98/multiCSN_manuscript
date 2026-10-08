# Project-local rendering of the saved Fig. 5 network. No inference is performed.
# Keep this helper independent of optional multiCSN plotting APIs: the pinned
# inference package does not implement the state/shape arguments used here.
plot_focal_network <- function(network_table, tf_nodes, node_state, state_colors,
                               node_size, label_nodes, seed = 42,
                               label.size = 3.5, edge_alpha = 0.55,
                               legend_position = "none") {
  edges <- as.data.frame(network_table)
  graph <- igraph::graph_from_data_frame(edges[, c("regulator", "target")],
                                        directed = FALSE)
  set.seed(seed)
  xy <- igraph::layout_with_fr(graph)
  nodes <- data.frame(name = igraph::V(graph)$name, x = xy[, 1], y = xy[, 2])
  nodes$state <- factor(unname(node_state[nodes$name]), levels = names(state_colors))
  nodes$kind <- ifelse(nodes$name %in% tf_nodes, "Focal TF", "Target")
  nodes$size <- unname(node_size[nodes$name])
  nodes$size[is.na(nodes$size)] <- 3
  edges$x <- nodes$x[match(edges$regulator, nodes$name)]
  edges$y <- nodes$y[match(edges$regulator, nodes$name)]
  edges$xend <- nodes$x[match(edges$target, nodes$name)]
  edges$yend <- nodes$y[match(edges$target, nodes$name)]
  edges$interaction <- ifelse(edges$weight > 0, "Activation", "Repression")
  ggplot2::ggplot() +
    ggplot2::geom_curve(data = edges, ggplot2::aes(x, y, xend = xend, yend = yend,
      colour = interaction, linewidth = abs(weight)), curvature = 0.08, alpha = edge_alpha) +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 0.85)) +
    ggplot2::scale_colour_manual(values = c(Activation = "#1F78B4", Repression = "#E31A1C")) +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x, y, fill = state, shape = kind,
      size = size), colour = "black", stroke = 0.35) +
    ggplot2::scale_fill_manual(values = state_colors) +
    ggplot2::scale_shape_manual(values = c("Focal TF" = 23, Target = 21)) +
    ggplot2::scale_size_identity() +
    ggrepel::geom_text_repel(data = nodes[nodes$name %in% label_nodes, ],
      ggplot2::aes(x, y, label = name), size = label.size,
      segment.colour = NA, max.overlaps = Inf, seed = seed, bg.color = "white") +
    ggplot2::theme_void() + ggplot2::theme(legend.position = legend_position)
}
