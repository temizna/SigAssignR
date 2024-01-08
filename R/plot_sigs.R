#' Plot signature contribution barplot
#'
#' Plot contribution of known signatures. Adopted from plot_contribution 
#' from MutationalPatterns
#'
#' @param fit Signature contribution matrix (Fit$coeffecient from assign_signatures)
#' @param signatures Signature matrix. Default is NA
#' @param mode "relative" or "absolute"; to plot the relative contribution or
#' absolute number of mutations, default = "relative"
#' @return Stacked barplot with contribution of each signature for each sample
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#'
#' @examples
#'
#' 
#' ## Plot the relative contribution
#' plot_sigs(fit$coefficient, mode="relative")
#'
#' ## Plot the absolute contribution.
#' ## When plotting absolute NMF results, the signatures need to be included.
#' plot_sigs(fit$coefficient,
#'   mode = "absolute"
#' )
#'
#'
#' @export
plot_sigs<-function(fit, mode = c("relative", "absolute"), signatures=NA){
  mode <- match.arg(mode)
  if(!(is.na(signatures))){
    fit=fit[,signatures]
  }
  #if there are less than 10 mutations assgined to a signature, remove those 
  # signatures from the plot
  zero_sigs=which(colSums(fit)<10)
  fit=fit[,-zero_sigs]
  Sample <- Contribution <- Signature <- NULL
  tb <-t(fit) %>% as.data.frame() %>% tibble::rownames_to_column("Signature") %>% 
    tidyr::pivot_longer(-Signature, names_to = "Sample", 
                        values_to = "Contribution") %>% dplyr::mutate(Sample = factor(Sample, 
                        levels = unique(Sample)), Signature = factor(Signature, 
                        levels = unique(Signature)))
  if (mode == "absolute") {
    bar_geom <- geom_bar(position=position_stack(reverse = TRUE),stat = "identity", colour = "black")
    y_lab <- "Absolute contribution \n (no. mutations)"
  } else if (mode == "relative") {
    bar_geom <- geom_bar(position = "fill", stat = "identity", 
                         colour = "black")
    y_lab <- "Relative contribution"
  }
  
  present_sigs <- tb %>% dplyr::filter(Contribution != 0) %>% 
    dplyr::pull(Signature) %>% unique()
  plot <- ggplot(tb, aes(x = Sample, y = Contribution, fill = Signature)) + 
    bar_geom + labs(x = "", y = y_lab) + scale_fill_discrete(breaks = present_sigs) + 
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), 
                       panel.grid.major.y = element_blank())+
    theme(axis.title.x=element_blank(),
                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
