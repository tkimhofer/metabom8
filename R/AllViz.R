#' Simple plotting of a single NMR spectrum
#' @export
#' @param ppm ppm vector.
#' @param x NMR spectrum.
#' @param shift chemical shift region to be plotted.
#' @param add Logical indicating if spectrum should be added to a current plot generated with \code{spec()} or \code{matspec()}.
#' @param interactive logical, interactive version (plotly package)
#' @param name string, name of trace (only used in interactive mode)
#' @param mode string, plot mode for interactive version: 'lines', 'lines+markers' or 'markers' (see Details)
#' @param ... Additional parameters to be passed on to the graphics generic plot function.
#' @seealso  \code{\link{matspec}} \code{\link{plot}}
#' @aliases spec
#' @details Low-level plotting function for a single NMR spectrum (base graphics). Interactive visualisation in mode of markers or line+markers will take longer (mode='lines' is recommended).
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom graphics points
#' @importFrom plotly plot_ly add_lines layout %>% add_trace
#' @family visNMR
#' @section
spec <- function(x, ppm, shift = c(0, 11), add = FALSE, interactive=TRUE, name='A', mode='lines', ...) {

  if(!is.null(ncol(x))){ stop('More than one spectrum provided.')}
  if(length(x)!=length(ppm)){stop('X and ppm don\'t match.')}
  idx <- get.idx(shift, ppm)

  if(interactive){
    sp=data.frame(ppm=ppm[idx], spec=x[idx])
    if(add){
      if(!exists('.ind_interactive', envir = parent.frame())){stop('First create interactive plot (set add=FALSE), then add desired spectrum')}
      p= get(".p", envir = parent.frame())
      p=p %>% add_trace(data=sp, y = ~spec, name = 'trace 1', mode = mode)
      assign(".p", p, envir =   parent.frame())
      return(p)
    }else{
     # if(!exists('mm8_plot', mode='environment')) mm8_plot <- new.env(parent = globalenv())
      assign(".ind_interactive", TRUE, envir =  globalenv())
      x <- list(title = "&delta;<sup>1</sup>H (ppm)", autorange="reversed")
      y <- list( title = "Intensity")
      p=plot_ly(data=sp, x=~ppm, y=~spec,   name='A', type = 'scatter', mode = mode,  hovertemplate = '%{x} ppm<extra></extra>') %>% layout(xaxis = x, yaxis=y)
      assign(".p", p, envir  = parent.frame())
    }
   return(p)
  }


  if(exists('.ind_interactive', envir = parent.frame())){ assign(".ind_interactive", FALSE, envir =  globalenv())}
  if (add) {
    points(ppm[idx], x[idx], type = "l", ...)
  } else {
      plot(ppm[idx], as.numeric(x[idx]), type = "l", xlim = rev(range(ppm[idx])), xlab = "ppm", ylab = "Intensity", ...)
    }

}
# spec(Xb[1,], ppm, shift=c(-1,11), name='A', interactive=T, mode='lines')
# spec(X[4,], ppm, shift=c(-1,11), name='C', add=T, interactive=T, mode='lines')





#' Simple plotting of multiple NMR spectra overlayed
#' @export
#' @param ppm ppm vector.
#' @param X NMR matrix with spectra represented in rows.
#' @param shift Chemical shift region to be plotted (in ppm).
#' @param interactive logical, interactive version (plotly package)
#' @param ... Additional parameters to be passed on to the graphics generic plot function.
#' @seealso \code{\link{spec}} \code{\link{plot}}
#' @aliases matspec
#' @details Low-level plotting function for NMR spectra, interactive plotting with ggplotly
#' @importFrom graphics matplot matpoints
# #' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @family visNMR
#' @section
matspec <- function(X, ppm, shift = c(0, 9.5), interactive=TRUE, ...) {
  if(is.null(ppm)){ppm=as.numeric(colnames(X)); } else{
    if(!.check_X_ppm(X, ppm)) stop('Non-matching dimensions X matrix and ppm vector or missing values in ppm.')
  }

  idx <- get.idx(shift, ppm)

  if(interactive){
    df=melt(X[,idx])
    x <- list(title = "&delta;<sup>1</sup>H (ppm)", autorange="reversed")
    y <- list( title = "Intensity")

    cols <- colorRampPalette(brewer.pal(8, "Set2"))(nrow(X))
    df$col=rep(cols, length(idx))

    p=plot_ly(data=df, x = ~Var2, y = ~value, color = ~col, name=~Var1, hovertemplate = '%{x} ppm<extra></extra>') %>% layout(xaxis = x, yaxis=y) %>% add_lines()
      return(p)
    }

    matplot(ppm[idx], t(X[, idx]), type = "l", xlim = rev(range(ppm[idx])), xlab = "ppm", ylab = "Intensity", ...)

}

# # Xb=bline(X)
# # matspec(Xb, ppm, shift=c(-0.5,4), interactive=T)
#
# idx=get.idx(c(6,8), ppm)
# plot(apply(Xb[,idx], 1, sum))



#' @title  Higher level plotting function to overlay NMR spectra (ggplot2 based)
#' @aliases specOverlay
#' @export
#' @param X Input NMR data matrix with row representing spectra.
#' @param ppm ppm vector with its length equals to \code{nrow(X)}.
#' @param shift Chemical shift area to be plotted. This should be kept as small as possible (see Details).
#' @param an List with one to three elements specifying facetting, colour and linetype (see Details).
#' @param alp Alpha value for lines (number between 0 and 1 whereas 0 is fully transparent and 1 is fully opaque).
#' @param title Plot title.
#' @param size Line width (0.5 is a good start).
#' @param ... Additional paramters passed on to ggplot's facet function.
#' @description  Plotting overlayed NMR specra. This function is based on ggplot2, a high-level plotting R package. For large ppm ranges the computation time is relatively long, so the chemical \code{shift} range should be as small as possible. For list argument \code{an}, the first element describes the colour and must be defined (even if it is only a single value). If colour and line width are specified, then at least one list elements of \code{an} must have the same length as \code{X}.
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes_string scale_x_reverse ggtitle xlab facet_grid theme_bw theme element_text geom_line scale_colour_gradientn
#' @importFrom colorRamps matlab.like2
#' @importFrom scales breaks_pretty
#' @importFrom stats as.formula
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @section
specOverlay <- function(X, ppm=NULL, shift = c(-0.01, 0.01), an = list("facet", "col", "ltype"), alp = 0.7, size = 0.5, title = "", ...) {

  if(is.null(ppm)){ppm=as.numeric(colnames(X)); } else{
    if(!.check_X_ppm(X, ppm)) stop('Non-matching dimensions X matrix and ppm vector or missing values in ppm.')
  }

  if (is.null(names(an))) {
    cat("No facet, colour and linetype names given. See an argument in ?specOverlay\n")
    names(an) <- paste("an", seq_len(length(an)), sep = "")
  }
  le.arg <- paste(length(an))
  if ("" %in% names(an)) {
    idx <- which(names(an) == "")
    names(an)[idx] <- paste("an", idx, sep = "")
  }
  names(an) <- gsub(" ", ".", names(an))
  idx <- get.idx(shift, ppm)
  specs <- X[, idx]
  colnames(specs) <- paste("Idx", idx, sep = "_")
  # create dataframe for ggplot function
  df <- data.frame(do.call(cbind.data.frame, an), ID = seq_len(nrow(specs)), alp, specs)
  colnames(df)[seq_len(le.arg)] <- names(an)
  df <- melt(df, id.vars = c("alp", "ID", names(an)))
  df$variable <- ppm[as.numeric(gsub("Idx_", "", df$variable))]
  # initiate generic ggplot object
  g <- ggplot() + scale_x_reverse(breaks = seq(shift[1], shift[2], by = abs(diff(shift))/20), name = expression(delta ~ {
  }^1 * H ~ "(ppm)")) + scale_y_continuous(breaks = breaks_pretty(), name = "Intensity") + ggtitle(title) + facet_grid(as.formula(paste(names(an)[1],
                                                                                                                                        "~ ."))) + theme_bw() + theme(axis.text = element_text(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1))
  # add colour and line type
  switch(le.arg, `1` = {
    g <- g + geom_line(data = df, aes_string(x = "variable", y = "value", group = "ID"), colour = "black", alpha = alp, size = size)
  }, `2` = {
    g <- g + geom_line(data = df, aes_string(x = "variable", y = "value", group = "ID", colour = names(an)[2]), alpha = alp, size = size)
    # add multi-colour gradient if colour vector is not factor/char
    col.cat <- is.factor(an[[2]]) | is.character(an[[2]]) | is.logical(an[[2]])
    if (!col.cat) {
      g <- g + scale_colour_gradientn(colors = matlab.like2(length(an[[2]])))
    }
  }, `3` = {
    an[[3]] <- factor(an[[3]])
    g <- g + geom_line(data = df, aes_string(x = "variable", y = "value", group = "ID", colour = names(an)[2], linetype = names(an)[3]), alpha = alp,
                       size = size)
    # add multi-colour gradient if colour vector is not factor/char
    col.cat <- is.factor(an[[2]]) | is.character(an[[2]]) | is.logical(an[[2]])
    if (!col.cat) {
      g <- g + scale_colour_gradientn(colors = matlab.like2(length(an[[2]])))
    }
  })
  return(g)
}


#' Overlay PCA or OPLS loadings with spectra
#' @export
#' @param mod PCA or OPLS model generated via \emph{MetaboMate} package functions.
#' @param shift ppm region to visualise.
#' @param pc index of principal component to visualise, set to 1 if input model is OPLS
#' @param type Type of loadings visualisation, either \code{'Statistical reconstruction'} or \code{'Backscaled'} (see Details).
#' @param an List with one to three elements specifying facetting, colour and linetype (see Details).
#' @param alp Alpha value for spectral lines.
#' @param title Plot title.
#' @param size plot line width.
#' @param r_scale logical, adjust limits of color gradient to 0 and 1 (only applies for type stat reconstruction)
#' @description  Plotting overlayed NMR spectra. This function is based on ggplot2, a high-level plotting R package. For high ppm ranges computation time is relatively, so the range of input argument \code{shift} should be as small as possible. List argument \code{an} must have the first element define, even if it is only a single value. If colour and line width is specified, then at least one list elements of \code{an} must have the same length as \code{X}.
#' @details OPLS: If \code{type='Statistical reconstruction'} the function calculates the covariance (y axis) and Pearson's correlation (colouring) of the predictive OPLS scores with each X variable (x axis is ppm variable). If \code{type='Backscaled'} the OPLS loadings are backscaled with X feature standard deviations. Results are plotted over ppm, coloured according to OPLS model weights. Often, the latter method visualises model importance more robust due to the presence of false positive correlations. PCA: Function always calculates the statistical reconstruction.
# @seealso \code{\link{plotload}} \code{\link{specOverlay}} \code{\link[=OPLS_MetaboMate-class]{OPLS_MetaboMate}} \code{\link{opls}} \code{\link[=PCA_MetaboMate-class]{PCA_MetaboMate}} \code{\link{pca}}
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes_string scale_x_reverse ggtitle xlab ylab facet_grid theme_bw theme element_text geom_line scale_colour_gradientn
#' @importFrom colorRamps matlab.like2
#' @importFrom stats as.formula
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @family dataviz
#' @section
specload <- function(mod, shift = c(0, 10), an, alp = 0.3, size = 0.5, pc = 1, type = "Backscaled", title = "", r_scale=FALSE) {

  if(!class(mod)[1] %in% c('OPLS_metabom8', 'PCA_metabom8')) {stop('Need metabom8 PCA or OPLS object.')}

  # calc loadings
  if (grepl("st|recon", type, ignore.case = TRUE)) {
    type <- "Statistical reconstruction"
  } else {
    type <- "Backscaled"
  }

  X=mod@X
  ppm=as.numeric(colnames(mod@X))


  idx <- get.idx(shift, ppm)
  if (length(idx) < 3) {
    stop("Shift area not in ppm variable.")
  }



  idx <- get.idx(shift, ppm)


  if (type == "Statistical reconstruction") {

    # extrat data
    #t_mod <- .viz_df_helper(mod, pc, an=NA, type='t')


    df_l <- .load_stat_reconstr_nmr(mod, pc, X, idx, ppm)
    if(r_scale) {raCol =  c(0, 1)}else{  raCol <- c(0, max(abs(df$cor))) }

    y <- df_l$cov
    cols <- df_l$cor

  }
  if (type == "Backscaled") {
    # extract data
    #p_mod=.viz_df_helper(mod, pc, an=NA, type='p')

    # backscaling
    df_l=.load_backscaled_nmr(mod, pc, idx, ppm)

    y <- df_l$p_bs
    cols <- abs(df_l[,1])

  }
#
#
#   #####################
#   # plot specs

  specs <- X[, idx]
  limY <- range(specs)

  colnames(specs) <- paste("ppm", ppm[idx], sep = "_")

  an=.check_an_viz(an, mod)

  df <- data.frame(ID = seq(nrow(specs)), do.call(cbind.data.frame, an), specs, row.names = NULL)
  df <- melt(df, id.vars = c("ID", names(an)))
  df$variable <- as.numeric(gsub("^\\.", "-", gsub("ppm_", "", df$variable)))
  colnames(df)[match(names(an)[1], colnames(df))]='facet'


  # scale loadings to appear in spectral range

  # scale colour line to dimensions of spectrum intensity
  cv1 <- (minmax(y) * (limY[2]/3)) + limY[2] * 0.67
  if (max(cv1) > limY[2]) {
    cv1 <- cv1 - abs(max(cv1 - limY[2]))
  }

  fac_lev <- unique(an[[1]])
  # define loadings
  df1 <- data.frame(alp, ID = nrow(X) + 1, facet = fac_lev[1], Group = "load", ppm = ppm[idx], Intensity = cv1, load = cols)

  # add loadings for each facet level
  if(length(fac_lev)>1){
    for (i in 2:length(fac_lev)) {
      df1 <- rbind(df1, data.frame(alp, ID = nrow(X) + i, facet = fac_lev[i], Group = "load", ppm = ppm[idx], Intensity = cv1, load = cols))
    }
  }

  g <- ggplot() +
    geom_line(data = df1, aes_string("ppm", "Intensity", color = "load", group = "ID"), size = 0.8) +
    geom_line(data = df, aes_string("variable", "value", group = "ID"), alpha = alp, size = 0.1) +
    facet_grid(facet ~ .)+
    scale_x_reverse(breaks = round(seq(shift[1], shift[2], by = abs(diff(shift))/20), 3)) +
    scale_y_continuous(limits = limY) +
    ggtitle(title) +
    xlab(expression(delta ~ {}^1 * H ~ "(ppm)")) +
    ylab("Intensity (AU)") + facet_grid(facet ~ .) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1))



  if (type == "Statistical reconstruction") {
    g <- g + scale_colour_gradientn(colors = matlab.like2(10), name = "cor(t,x)", limits = raCol)
  } else {
    g <- g + scale_colour_gradientn(colors = matlab.like2(10), name = expression(abs ~ w[pred * "," ~ sc]))
  }

  return(g)
}



#' Plotting PCA or OPLS loadings
#' @export
#' @param mod PCA or OPLS model generated via \emph{metabom8} package functions.
#' @param shift ppm region to visualise.
#' @param pc index of principal component to visualise, set to 1 if input model is OPLS
#' @param type Type of loadings visualisation, either \code{'Statistical reconstruction'} or \code{'Backscaled'} (see Details).
#' @param title Plot title.
#' @param r_scale logical, adjust limits of color gradient to 0 and 1 (only applies for type stat reconstruction)
#' @details OPLS: If \code{type='Statistical reconstruction'} the function calculates the covariance (y axis) and Pearson's correlation (colouring) of the predictive OPLS scores with each X variable (x axis is ppm variable). If \code{type='Backscaled'} the OPLS loadings are backscaled with X feature standard deviations. Results are plotted over ppm, coloured according to OPLS model weights. Often, the latter method visualises model importance more robust due to the presence of false positive correlations. PCA: Function always calculates the statistical recostruction.
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @references Cloarec, O., \emph{et al.} (2005). Evaluation of the Orthogonal Projection on Latent Structure Model Limitations Caused by Chemical Shift Variability and Improved Visualization of Biomarker Changes in 1H NMR Spectroscopic Metabonomic Studies. \emph{Analytical Chemistry} 77.2, 517-26.
#' @importFrom stats cor cov
#' @importFrom ggplot2 ggplot geom_line scale_x_reverse ggtitle xlab ylab theme_bw ggtitle aes_string scale_colour_gradientn geom_point
#' @importFrom colorRamps matlab.like2
#' @importFrom scales breaks_pretty
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @family dataviz
#' @section

plotload <- function(mod, shift = c(0, 10), pc = 1, type = "Backscaled", title = NULL, r_scale=FALSE) {

  if(!class(mod)[1] %in% c('OPLS_metabom8', 'PCA_metabom8')) {stop('Need metabom8 PCA or OPLS object.')}

  if (grepl("st|recon", type, ignore.case = TRUE)) {
    type <- "Statistical reconstruction"
  } else {
    type <- "Backscaled"
  }

  X=mod@X
  ppm=as.numeric(colnames(mod@X))

  idx <- get.idx(shift, ppm)


  if (type == "Statistical reconstruction") {
    # extrat data
    #p_mod <- .viz_df_helper(mod, pc, an=NA, type='p')

    df <- .load_stat_reconstr_nmr(mod, pc, X, idx, ppm)

    if(r_scale) {raCol =  c(0, 1)}else{  raCol <- c(0, max(abs(df$cor))) }

    g <- ggplot(df, aes_string("ppm", "cov", colour = "cor")) +
      geom_line() +
      scale_x_reverse(breaks = breaks_pretty(n = 15)) +
      scale_colour_gradientn(colors = matlab.like2(10),
                            limits = raCol,
                            name = "r") +
      labs(title=title, x=expression(delta ~ {}^1 * H ~ (ppm)), y="cov(t,x)",
           caption = paste(gsub('_metabom8','', class(mod)[1]), '-', mod@type, 'component', pc)) +
      theme_bw()
  }
  if (type == "Backscaled") {
    # extract data
    #p_mod=.viz_df_helper(mod, pc, an=NA, type='p')

    # backscaling
    df=.load_backscaled_nmr(mod, pc, idx, ppm)

    g <- ggplot(df, aes_string("ppm", "p_bs", colour = "p_abs")) +
      geom_line() +
      scale_x_reverse(breaks = breaks_pretty(n = 15)) +
      scale_colour_gradientn(colors = matlab.like2(10),
                            name = expression('|p'['sc']*'|')) +
      labs(title=title, x=expression(delta ~ {}^1* H ~ (ppm)), y=expression(p*'*'*sigma[x]),
           caption = paste(gsub('_metabom8','', class(mod)[1]), '-', mod@type, 'component', pc)) +
      theme_bw()
  }
  return(g)
}



#' Calculating distance to the model in X space
#' @export
#' @param mod OPLS model of type \code{OPLS_metabom8}.
#' @param plot Logical indicating if results should be visualised.
#' @return The projection distance of each observation in the model (\code{DModX}).
#' @references Bylesjö M., \emph{et al.} (2002) OPLS discriminant analysis: combining the strengths of PLS-DA and SIMCA classification. \emph{Journal of Chemometrics}, 20, 341-51.
#' @references Wold S. (1976) Pattern recognition by means of disjoint principal components models.  \emph{Pattern Recognition}, 8, 127-39.
#' @seealso \code{\link{opls}}
#' @importFrom ggplot2 ggplot aes_string geom_point scale_colour_gradientn geom_hline xlab scale_y_continuous theme_bw theme element_blank element_text geom_segment
#' @importFrom colorRamps matlab.like
#' @importFrom stats t.test sd
#' @author Torben Kimhofer \email{torben.kimhofer@@murdoch.edu.au}
#' @family OPLS model validation functions
#' @section
# E=residual Matrix N=number of samples K=number of variables A=number of model components A0= (1 if mean centred, 0 otherwise)

dmodx <- function(mod, plot = TRUE) {
  if (class(mod)[1] != "OPLS_metabom8") {
    stop("Please provide a OPLS_metabom8 object.")
  }
  E <- mod@X_res
  N <- nrow(E)
  K <- ncol(E)
  A <- ncol(mod@t_pred)  # in case of OPLS-DA (alwasy one predictive component)
  if (mod@Parameters$center) {  A0 <- 1 } else { A0 <- 0 }
  # loop over all observations in residual matrix, calc SS residuals / observations and normalise by TSS
  ss_res <- apply(E, 1, function(x) sum(x^2))
  dmodX <- sqrt(ss_res/(K - A))/sqrt(sum(ss_res)/((N - A - A0) * (K - A)))
  tt <- t.test(dmodX, alternative = "less")
  ci95 <- tt$conf.int[2] + 2 * sd(dmodX)
  df <- data.frame(col = mod@t_pred_cv[,1], ID = seq_len(length(dmodX)), DmodX = dmodX, passedT.test = dmodX < tt$conf.int[2] + 2 * sd(dmodX))
  if(plot) {
    df$Y=mod@Y$ori
    g <- ggplot(data = df) +
      geom_segment(aes_string(x = "ID", xend = "ID", y = "min(dmodX)-0.1", yend = "DmodX"), colour = "gray60", size = 0.1) +
      geom_point(aes_string(x = "ID", y = "DmodX", colour = "Y")) +
      #scale_colour_gradientn(colours = matlab.like2(10), name = expression(t[pred])) +
      scale_y_continuous(limits = c(min(dmodX) - 0.1, max(c(dmodX, ci95)) + 0.2), name = "DModX", expand = c(0, 0)) +
      geom_hline(yintercept = ci95,  linetype = 2, colour = "black") + xlab("Sample index") + theme_bw() +
      theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(), axis.text = element_text(colour = "black"))+
      labs(caption = 'Dashed line indicates uppler limit of 95% CI')
    if(mod@type == 'R'){
     g<- g +
       scale_colour_gradientn(colours = matlab.like2(10), name = expression(t[pred]))
    }
    plot(g)
  }
  return(df[, -1])
}

