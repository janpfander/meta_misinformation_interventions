# set general theme for plots
plot_theme <- theme_minimal(base_size = 12) +
  theme(# Bold title
    plot.title = element_text(face = "bold", size = rel(0.7)),
    # Plain subtitle that is grey
    plot.subtitle = element_text(face = "plain", size = rel(1), color = "grey70"),
    # Bold legend titles
    legend.title = element_text(face = "bold"),
    # Bold axis titles
    axis.title = element_text(face = "bold"), 
    # Bold legend labels
    legend.text = element_text(face = "bold")
    )