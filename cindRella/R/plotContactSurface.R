# a script which gets the contact matrix of a shoe allows you to plot it

#' A shoe plotting function
#'
#' This function is for making a ggplot2 plot of a shoe. It gives options to include accidentals and their labels.
#' @param shoenumber The number of the shoe to be plotted.
#' @param show_accidentals A boolean for whether or not to include the accidentals in the plot. Defaults to FALSE.
#' @param shoe_colour A string indicating which colour the contact surface should be. Defaults to chocolate1
#' @param accidental_colour A string indicating which colour the accidentals should be. Defaults to blue.
#' @param include_axes A boolean for whether or not to include axes ticks and scales in the plot. Defaults to FALSE.
#' @keywords plot shoe
#' @export
#' @examples
#' plotContactSurface()


plotContactSurface <- function(shoenumber, show_accidentals = FALSE, shoe_colour = "chocolate1", accidental_colour = "blue", include_axes = F){
  surface <- subset(contacts_example, shoenum == shoenumber)
  w <- ggplot2::ggplot(surface)
  p <- w + geom_raster(aes(x, y), fill = shoe_colour) + guides(fill = F, alpha = F)
  if(show_accidentals){
      p <- p + geom_point(data = subset(accidentals_example, (shoenum == shoenumber)), aes(x = x, y = y), color = accidental_colour)
  }
  return(default_shoe_axes(p, include_axes) + theme(panel.background = element_rect(fill='white', colour='white')))
}

#' Defaulting the shoe axes for good scaling.
#'
#' A function for defaulting the axes when plotting so the shoe is scaled well. Also gives option for removing axes.
#' @param theplot A ggplot object representing a shoe to be scaled
#' @param include_axes Whether or not to include the scale, axis labels, and axis ticks. A boolean which defaults to FALSE.
#' @keywords default
#' @export
#' @examples
#' default_shoe_axes()
default_shoe_axes <- function(theplot, include_axes = FALSE){
  if(include_axes == FALSE){
    theplot + scale_y_continuous(labels=NULL, breaks = NULL, limits = c(0, 200)) + scale_x_continuous(labels = NULL, breaks = NULL, limits = c(0, 100)) + xlab("") + ylab("")
  }
  else{
   theplot + scale_y_continuous(limits = c(0, 200)) + scale_x_continuous(limits = c(0, 100))
  }

}

dfposterior <- function(fit, contact, w_indices, zoptions, nzoptions, optiondiffsx, optiondiffsy, grid_memberships, chain_indices){
  nonzeros <- which(is.na(grid_memberships) == F)
  q <- fit$q[chain_indices]
  w <- matrix(ncol = length(q), nrow = length(fit$w[[1]]))
  Phi <- matrix(ncol = length(q),nrow = length(fit$phi[[1]]))
  kernelx <- matrix(ncol = length(q), nrow = length(fit$kernelx[[1]]))
  kernely <- matrix(ncol = length(q), nrow = length(fit$kernely[[1]]))
  for(i in 1:length(q)){
    w[,i] <- fit$w[[chain_indices[i]]]
    Phi[,i] <- fit$phi[[chain_indices[i]]]
    kernelx[,i] <- fit$kernelx[[chain_indices[i]]]
    kernely[,i] <- fit$kernely[[chain_indices[i]]]
  }
  xx <- plot_posterior(contact, w, Phi, kernelx, kernely,q, w_indices, zoptions, nzoptions, optiondiffsx, optiondiffsy)
  x = rep(seq(from = 0.005,to = 0.995, length.out = 100), 200)
  y = rep(seq(from = 0.0025,to = 0.9975, length.out = 200), each = 100)
  dens <- rep(0, 100 * 200)
  dens[nonzeros] <- exp(xx[[3]])
  x <- x * 100
  y <- y * 200
  dat <- data.frame(x, y, density = dens)
  return(dat)
}


