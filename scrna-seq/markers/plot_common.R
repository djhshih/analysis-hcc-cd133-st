# this script contains common functions for plot.R

summarize_gset <- function(sce, gset, groups) {
  
  idx <- match(tolower(gset), tolower(rowData(sce)$Symbol));
  valid <- !is.na(idx);
  message("matched fraction: ", mean(valid))
  idx.valid <- idx[valid];
  
  expr <- logcounts(sce[idx.valid, ]);
  
  ms <- scoreMarkers(expr, groups);
  
  d <- do.call(rbind, mapply(
    function(m, k) {
      data.frame(
        group = k,
        gene = rownames(m),
        logfc = m[, "self.average"] - m[, "other.average"],
        p_expressed = m[, "self.detected"]
      )
    },
    ms, names(ms),
    SIMPLIFY = FALSE
  ));
  rownames(d) <- NULL;
  
  d
}

plot_graph <- function(d, title) {
    ggplot(d,
           aes(
             x = group, y = gene,
             colour = bound(logfc, c(-1.5, 1.5)),
             size = p_expressed * 100
           )
    ) +
      theme_classic() +
      labs(title = title) +
      geom_point() +
      facet_grid(gset ~ ., scales="free_y", space="free", switch="y") +
      scale_colour_gradient2(breaks = (-1):1,
                             low="royalblue4", 
                             mid="grey90",
                             high="orangered2"
      ) +
      guides(
        colour = guide_colourbar("log FC"),
        size = guide_legend("% expressed")
      ) +
      theme(
        strip.text.y.left = element_text(angle = 0, hjust=0),
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1)
      ) +
      scale_size_continuous(range = c(-0.5, 3.5)) +
      scale_y_discrete(limits = rev) +
      xlab("") + ylab("")
}
