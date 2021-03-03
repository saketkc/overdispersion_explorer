MeanVarFit <- function(counts) {
  means <- sparseMatrixStats::rowMeans2(counts)
  variance <- sparseMatrixStats::rowVars(counts)
  df <- data.frame(mean = means, variance = variance)
  rownames(df) <- rownames(counts)
  df["var_alpha_0"] <- df$mean
  for (alpha in c(0.01, 0.1, 1)) {
    df[paste0("var_alpha_", alpha)] <- df$mean + alpha * df$mean^2
  }
  return(df)
}

ZerosFractionFit <- function(counts) {
  n_cells <- dim(counts)[2]
  means <- sparseMatrixStats::rowMeans2(counts)
  # variance <- sparseMatrixStats::rowVars(counts)
  zeros.df <- data.frame(
    mean = means,
    observed_zeros = apply(counts == 0, 1, sum) / n_cells,
    gene = rownames(counts)
  )

  zeros.df["expected_zeros_nb0"] <- exp(-zeros.df$mean) #
  for (alpha in c(0.01, 0.1, 1)) {
    zeros.df[paste0("expected_zeros_nb", alpha)] <- (1 / (alpha * means + 1))^(1 / alpha)
  }

  return(zeros.df)
}

PlotMeanVar <- function(cm) {
  meanvarfit <- MeanVarFit(cm)
  meanvarfit$gene <- rownames(meanvarfit)
  meanvarfit_melt <- melt(meanvarfit, id.vars = c("gene", "mean"))
  meanvarfit_melt$variable <- factor(
    x = as.character(meanvarfit_melt$variable),
    levels = c(
      "variance", "var_alpha_0",
      "var_alpha_0.01", "var_alpha_0.1", "var_alpha_1"
    )
  )
  p <- ggplot(meanvarfit_melt[meanvarfit_melt$variable == "variance", ], aes(mean, value)) +
    geom_scattermore(size = 0.71) +
    geom_line(
      data = meanvarfit_melt[meanvarfit_melt$variable != "variance", ],
      mapping = aes(mean, value,
        group = variable,
        color = variable, linetype = variable
      ), size = 0.9
    ) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    xlab("Mean") +
    ylab("Variance")
  labels <- c(
    expression(mu),
    expression(paste(mu, "+", 0.01, mu^2)),
    expression(paste(mu, "+", 0.1, mu^2)),
    expression(paste(mu, "+", mu^2))
  )
  # p <- p +   annotation_logticks()
  p <- p +
    scale_color_manual(
      values = c(
        "blue", "#fdcc8a",
        "#fc8d59", "#d7301f"
      ),
      name = "",
      labels = labels
    ) +
    scale_linetype_manual(
      name = "",
      values = c(1, 2, 4, 5),
      labels = labels
    ) +
    # theme(legend.position = c(0, 1),
    #      legend.justification = c(-0.1, 1), legend.background=element_blank())
    theme(
      legend.position = "bottom", legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.text = element_text(size = 12)
    )

  return(p)
}


PlotZerosFit <- function(cm) {
  zerosfit <- ZerosFractionFit(cm)
  zerosfit$gene <- rownames(zerosfit)
  zerosfit_melt <- melt(zerosfit, id.vars = c("gene", "mean"))
  zerosfit_melt$variable <- factor(
    x = as.character(zerosfit_melt$variable),
    levels = c(
      "observed_zeros", "expected_zeros_nb0",
      "expected_zeros_nb0.01", "expected_zeros_nb0.1", "expected_zeros_nb1"
    )
  )
  p <- ggplot(zerosfit_melt[zerosfit_melt$variable == "observed_zeros", ], aes(mean, value)) +
    # geom_point(size = 0.71) +
    geom_scattermore(size = 0.71) +
    geom_line(
      data = zerosfit_melt[zerosfit_melt$variable != "observed_zeros", ],
      mapping = aes(mean, value,
        group = variable,
        color = variable, linetype = variable
      ), size = 0.9
    ) +
    xlab("Mean") +
    ylab("Fraction of Zeros") +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  labels <- c(
    expression(mu),
    expression(paste(mu, "+", 0.01, mu^2)),
    expression(paste(mu, "+", 0.1, mu^2)),
    expression(paste(mu, "+", mu^2))
  )
  p <- p +
    scale_color_manual(
      values = c(
        "blue", "#fdcc8a",
        "#fc8d59", "#d7301f"
      ),
      name = "",
      labels = labels
    ) +
    scale_linetype_manual(
      name = "",
      values = c(1, 2, 4, 5),
      labels = labels
    ) +
    # theme(legend.position = c(0, 1),
    #      legend.justification = c(-0.1, 1), legend.background=element_blank())
    theme(
      legend.position = "bottom", legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.text = element_text(size = 12)
    )

  return(p)
}

residualVarPlot <- function(gene_var, col = "residual_variance") {
  col <- gsub("-", "_", col)
  top20 <- subset(gene_var, rank(-gene_var[, col]) <= 30)$gene

  p <- ggplot(gene_var, aes_string("gmean", col)) +
    # geom_point(size = 0.6, shape = 16, alpha = 0.5) +
    # geom_point(data = subset(gene_var, gene %in% top20), size = 0.6, shape = 16, alpha = 1.0, color = "deeppink") +
    geom_scattermore(size = 0.6, shape = 16, alpha = 0.5) +
    geom_scattermore(data = subset(gene_var, gene %in% top20), size = 0.6, shape = 16, alpha = 1.0, color = "deeppink") +
    geom_hline(yintercept = 1, color = "red") +
    geom_smooth(method = "loess", span = 0.1, size = 1) +
    # geom_smooth(method = "gam", size = 1) +
    scale_y_continuous(trans = "sqrt", breaks = c(0, 1, 50, 100, 150), limits = c(0, 230)) +
    scale_x_continuous(trans = "log10", breaks = c(0.001, 0.01, 0.1, 1, 10), labels = MASS::rational) +
    # scale_y_continuous(trans='log1p') +
    # facet_wrap(~ model, ncol=3, scales = 'free_y') +
    xlab("Gene mean") +
    ylab("Residual variance") +
    geom_text_repel(
      data = subset(gene_var, gene %in% top20), aes(label = gene), color = "gray25",
      size = 1.8,
      nudge_y = 230 - subset(gene_var, gene %in% top20)[, col],
      direction = "x",
      angle = 90,
      vjust = 0.5,
      hjust = 0.5,
      segment.size = 0.2,
      segment.alpha = 0.2
    )
  return(p)
}
