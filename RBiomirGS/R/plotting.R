#' @title rbiomirgs_volcano
#'
#' @description volcano plot function
#' @param gsadfm Input dataframe. The dataframe is the output from \code{\link{rbiomirgs_logistic}}.
#' @param topgsLabel If to display the top GS identification on the plot. Default is \code{FALSE}.
#' @param n Number of top GS name to display. Only set this argument when \code{topgsLabel = TRUE}. Default is \code{5}.
#' @param gsVar The name of the variable for GS name. Only set this argument when \code{geneName = TRUE}. Default is \code{GS}.
#' @param padding When \code{topgsLabel = TRUE}, to set the distance between the dot and the gene symbol. Default is \code{0.5}.
#' @param gsLabelSize When \code{topgsLabel = TRUE}, to set the label size. Default is \code{5}.
#' @param logoddsratio Threshold for logoddsratio for volcano plot. Default is \code{0}.
#' @param fdr If to use fdr corrected p value for plotting. Default is \code{TRUE}.
#' @param q_value Threshold for the p value. Default is \code{0.05}.
#' @param Title Figure title. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param yLabel Y-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param symbolSize Size of the symbol. Default is \code{2}.
#' @param sigColour Colour of the significant genes or probes. Default is \code{"red"}.
#' @param nonsigColour Colour of the non-significant genes or probes. Default is \code{"gray"}.
#' @param xTxtSize Font size for the x-axis text. Default is \code{10}.
#' @param yTxtSize Font size for the y-axis text. Default is \code{10}.
#' @param plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @details The function will use uncorrected p value for thresholding if no significant results were found under fdr. A message will display.
#' @return A \code{pdf} file of the volcano plot.
#' @import ggplot2
#' @importFrom grid grid.newpage grid.draw
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom ggrepel geom_text_repel
#' @examples
#' \dontrun{
#'
#' rbiomirgs_volcano(gsadfm = mirna_mrna_gs, n = 5)
#'
#' }
#' @export
rbiomirgs_volcano <- function(gsadfm,
                             topgsLabel = FALSE, n = 5, gsVar = "GS", padding = 0.5, gsLabelSize = 5,
                             logoddsratio = 0, fdr = TRUE, q_value = 0.05,
                             Title = NULL, xLabel = "Log odds ratio", yLabel = "-log10(p value)",
                             symbolSize = 2, sigColour = "red", nonsigColour = "gray",
                             xTxtSize = 10, yTxtSize =10,
                             plotWidth = 170, plotHeight = 150,
                             boxplot = TRUE){

  # set up the data frame
  tmpdfm <- gsadfm

  # set the cutoff
  if (fdr){
    if (length(which(tmpdfm$adj.p.val < q_value)) == 0){
      message(cat("No significant result was found under fdr correction. Proceed thresholding is conducted on raw p value."))
      pcutoff <- q_value
    } else {
      pcutoff <- max(tmpdfm[tmpdfm$adj.p.val < q_value, ]$p.value)
    }
  } else {
    pcutoff <- q_value
  }
  cutoff <- as.factor(abs(tmpdfm$coef) >= logoddsratio & tmpdfm$p.value < pcutoff)

  # plot
  loclEnv <- environment()
  plt <- ggplot(tmpdfm, aes(x = coef, y = -log10(p.value)), environment = loclEnv) +
    geom_point(alpha = 0.4, size = symbolSize, aes(colour = cutoff)) +
    scale_color_manual(values = c(nonsigColour, sigColour)) +
    ggtitle(Title) +
    scale_y_continuous(expand = c(0.02, 0)) +
    xlab(xLabel) +
    ylab(yLabel) +
    geom_vline(xintercept = logoddsratio, linetype = "dashed") +
    geom_vline(xintercept = - logoddsratio, linetype = "dashed") +
    geom_hline(yintercept = - log10(pcutoff), linetype = "dashed") +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          legend.title = element_blank(),
          axis.text.x = element_text(size = xTxtSize),
          axis.text.y = element_text(size = yTxtSize, hjust = 0.5))

  if (topgsLabel){
    tmpfltdfm <- tmpdfm[abs(tmpdfm$coef) >= logoddsratio & tmpdfm$p.value < pcutoff, ]
    tmpfltdfm <- tmpfltdfm[order(tmpfltdfm$p.value), ]
    plt <- plt + geom_text_repel(data = head(tmpfltdfm, n = n),
                                 aes(x = coef, y = -log10(p.value), label = head(tmpfltdfm, n = n)[, gsVar]),
                                 point.padding = unit(padding, "lines"), size = gsLabelSize)
  }

  grid.newpage()

  # extract gtable
  pltgtb <- ggplot_gtable(ggplot_build(plt))

  # add the right side y axis
  Aa <- which(pltgtb$layout$name == "axis-l")
  pltgtb_a <- pltgtb$grobs[[Aa]]
  axs <- pltgtb_a$children[[2]]
  axs$widths <- rev(axs$widths)
  axs$grobs <- rev(axs$grobs)
  axs$grobs[[1]]$x <- axs$grobs[[1]]$x - unit(1, "npc") + unit(0.08, "cm")
  Ap <- c(subset(pltgtb$layout, name == "panel", select = t:r))
  pltgtb <- gtable_add_cols(pltgtb, pltgtb$widths[pltgtb$layout[Aa, ]$l], length(pltgtb$widths) - 1)
  pltgtb <- gtable_add_grob(pltgtb, axs, Ap$t, length(pltgtb$widths) - 1, Ap$b)

  # export the file and draw a preview
  ggsave(filename = paste(deparse(substitute(gsadfm)),".volcano.pdf", sep = ""), plot = pltgtb,
         width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
  grid.draw(pltgtb) # preview
}


#' @title rbiomirgs_bar
#'
#' @description Bar graph function
#' @param gsadfm Input dataframe. The dataframe is the output from \code{\link{rbiomirgs_logistic}}.
#' @param gs.name Whether or not to display gene set names on the y-axis. Default is \code{FALSE}.
#' @param n Number of gene sets to plot. Options are \code{"all"} and an integer number. Default is \code{"all"}.
#' @param signif_only Whether to plot only the significantly enriched GS. Default is \code{FALSE}.
#' @param p_adj Set when \code{signif_only = TRUE}, whether to use adjusted p value for thresholding. Default is \code{TRUE}.
#' @param q_value Set when \code{signif_only = TRUE}, the significance threshold. Default is \code{0.05}.
#' @param y.rightside Whether or not to display right side y-axis. Default is \code{FALSE}.
#' @param Title Title of the plot. Default is \code{NULL}.
#' @param xLabel X-axis label. Default is \code{"Log  odds ratio"}.
#' @param yLabel Y-axis label. Default is \code{NULL}.
#' @param plotWidth Set the width of the plot. Default is \code{250}.
#' @param plotHeight Set the height of the plot. Default is \code{230}.
#' @details The function produces either a bar graph or a vertical version. When \code{gs.name = FALSE}, the function produces a plot with x-axis as gene set. When \code{gs.name = TRUE}, the plot will have gene set names displayed on the y-axis, and x-axis will be used for log oods ratio. The error bar is standard error of the coefficient (log odds ratio).
#' @return Outputs a \code{pdf} boxplot figure.
#' @import ggplot2
#' @importFrom grid grid.newpage grid.draw
#' @examples
#' \dontrun{
#'
#' # horizontal version
#' rbiomirgs_bar(gsadfm = mirna_mrna_GS, n = "all", y.rideside = TRUE, yTxtSize = 8, plotWidth = 300, plotHeight = 200, xLabel = "Gene set", yLabel = "Log odds ratio")
#'
#' # vertical version
#'  rbiomirgs_bar(gsadfm = mirna_mrna_GS, n = 60, yAxis = TRUE,
#'  yTxtSize = 8, plotWidth = 300, plotHeight = 200, xLabel = "Log odds ratio")
#'
#' }
#' @export
rbiomirgs_bar <- function(gsadfm,
                                gs.name = FALSE,  n = "all", signif_only = FALSE, p_adj = TRUE, q_value = 0.05,
                                y.rightside = FALSE,
                                Title = NULL, xLabel = "Log odds ratio", yLabel = NULL,
                                errorbarWidth = 0.2,
                                xTxtSize = 10, yTxtSize = 10,
                                plotWidth = 250, plotHeight = 230){

  # prepare the dataframe for ggplot2
  if (signif_only){ # if significant GS only
    if (p_adj){
      dfm <- gsadfm[gsadfm$adj.p.val < q_value, ]
    } else {
      dfm <- gsadfm[gsadfm$p.value < q_value, ]
    }
  } else {
    dfm <- gsadfm
  }

  dfm <- dfm[order(abs(dfm[, "coef"]), decreasing = TRUE), ] # sort according to the absolute value of the coef

  # plot
  # prepare plotting dataframe
  if (n != "all"){
    pltdfm <- head(dfm, n)
  } else {
    pltdfm <- dfm
  }

  pltdfm <- pltdfm[order(pltdfm$coef, decreasing = TRUE), ] # sort accroding to the raw value of the coef

  pltdfm$GS <- factor(pltdfm$GS, levels = pltdfm$GS)

  # plotting
  loclEnv <- environment()
  baseplt <- ggplot(pltdfm, aes(x = GS, y = coef), environment = loclEnv) +
    geom_bar(position = "dodge", stat = "identity", color = "black", fill = "gray66")+
    scale_x_discrete(expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    ggtitle(Title) +
    geom_hline(yintercept = 0) +
    geom_errorbar(aes(ymin = coef - std.err, ymax = coef + std.err), width = errorbarWidth,
                  position = position_dodge(0.9))

  if (gs.name){
    plt <- baseplt +
      xlab(yLabel) + # the arguments for x and y labls are switched as the figure will be rotated
      ylab(xLabel) + # the arguments for x and y labls are switched as the figure will be rotated
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(hjust = 0.5),
            legend.position = "bottom",
            legend.title = element_blank(),
            axis.text.x = element_text(size = xTxtSize, angle = 0, hjust = 0.5), # x and y not reversed as they are not associated with the roation of the axes.
            axis.text.y = element_text(size = yTxtSize, hjust = 1)) +
      coord_flip()
  } else { # the axes will not be flipped if not to display GS names
    plt <- baseplt +
      xlab(xLabel) +
      ylab(yLabel) +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(hjust = 0.5),
            legend.position = "bottom",
            legend.title = element_blank(),
            axis.text.y = element_text(size = xTxtSize, angle = 0, hjust = 0.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }

  # prepare for preview
  grid.newpage()

  if (y.rightside){   ## add the right-side y axis
    # extract gtable
    pltgtb <- ggplot_gtable(ggplot_build(plt))

    # add the right side y axis
    Aa <- which(pltgtb$layout$name == "axis-l")
    pltgtb_a <- pltgtb$grobs[[Aa]]
    axs <- pltgtb_a$children[[2]]
    axs$widths <- rev(axs$widths)
    axs$grobs <- rev(axs$grobs)
    axs$grobs[[1]]$x <- axs$grobs[[1]]$x - unit(1, "npc") + unit(0.08, "cm")
    Ap <- c(subset(pltgtb$layout, name == "panel", select = t:r))
    pltgtb <- gtable_add_cols(pltgtb, pltgtb$widths[pltgtb$layout[Aa, ]$l], length(pltgtb$widths) - 1)
    pltgtb <- gtable_add_grob(pltgtb, axs, Ap$t, length(pltgtb$widths) - 1, Ap$b)
  } else {
    pltgtb <- ggplot_gtable(ggplot_build(plt))
  }

  # export the file and draw a preview
  ggsave(filename = paste(deparse(substitute(gsadfm)),".bar.plot.pdf", sep = ""), plot = pltgtb,
         width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
  grid.draw(pltgtb) # preview
}
