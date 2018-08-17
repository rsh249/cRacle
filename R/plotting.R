#' @import wesanderson
NULL


#' Plot PDF curve of given type and variable for single density object.
#' 
#' Using an object from the densform() function. Plot a single PDF curve in a new plot window.
#' @param dens.ob An object derived from the cRacle::dens_ob() function.
#' @param var A character string that matches one of the layer names in the source raster object.
#' @param col A color declaration. Default is a random color.
#' @param type A character string of value either ".kde" for a Kernel Density Estimator curve, or ".gauss" for a Gaussian (normal) curve. All other values will result in errors.
#' @param w TRUE or FALSE to show weighted probability functions
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr); data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' extr.sub = subset(extr.raw, extr.raw$tax == extr.raw[5,'tax']);
#' dens.sub = densform(extr.sub, clim = climondbioclim, bw = 'nrd0', n = 512);
#' densplot(dens.sub, names(climondbioclim[[1]]));
#' }

densplot <- function(dens.ob, var, col = sample(grDevices::colours()), type = ".kde", w=FALSE) {
  
  varx <- paste(var, "x", sep = ".");
  varw <- paste(var, "w", sep = ".");
  
  graphics::par(mar= c(5,4,4,4) + 0.3);
  tempvarlist <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "MAT", "MaximumT", "MinimumT");
  if(var %in% tempvarlist){by = 10}else{by = 1};
  var <- paste(var, type, sep = "");
  if(w ==TRUE){
    to <- max(stats::na.omit(dens.ob[[varx]]));
    from <- min(stats::na.omit(dens.ob[[varx]]));
    num = length(dens.ob[[varx]]);
    lby = (to - from)/num;
    dens.ob[[var]] = dens.ob[[var]]^dens.ob[[varw]];
    den.area <- sum(stats::na.omit(dens.ob[[var]]))*lby;
    dens.ob[[var]] = dens.ob[[var]]/den.area
  }
  graphics::plot(dens.ob[[varx]]/by, dens.ob[[var]], xlab = "", ylab = "", ylim = c(0, 3.5*max(stats::na.omit(dens.ob[[var]]))), type = "l", lwd = 3, col = col, frame.plot=F, axes = F);
  graphics::axis(side = 2, at = pretty(c(0, 2.5*max(stats::na.omit(dens.ob[[var]])))));
  graphics::axis(side = 1, at = pretty(range(stats::na.omit(dens.ob[[varx]])/by)));
  graphics::mtext(var, side = 1, line =3);
  graphics::mtext("Probability Density Estimation", side = 2, line = 3);
};



#' Plot PDF curves of given type and variable for density list object
#' 
#' Using an object from the cRacle::dens_obj() function. Plot a series PDF curves in a new plot window.
#' @param dens.oblist An object derived from the cRacle::dens_obj() function.
#' @param var A character string that matches one of the layer names in the source raster object.
#' @param col A color vector of the same length as the dens.oblist object. Default is 'wes_palette("FantasticFox1")'.
#' @param type A character string of value either ".kde" for a Kernel Density Estimator curve, or ".gauss" for a Gaussian (normal) curve. All other values will result in errors.
#' @param l.pos  Legend position. Recommend 'topleft' or 'topright'. Default is 'topleft'.
#' @param l.cex  cex setting for legend. Default is 0.8.
#' @param w TRUE or FALSE to show weighted probability functions
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr); data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' }

multiplot <- function(dens.oblist, var, col = wes_palette("FantasticFox1"), type = ".kde", l.pos = 'topleft', l.cex = 0.8, w = FALSE){ 
  arr.dens.ob = dens.oblist;
  varx <- paste(var, "x", sep = ".");
  vart = paste(var, type, sep = '');
  
  current <- arr.dens.ob[[1]];
  densplot(current, var, col[1], type = type, w=w);
  max.x.hold = list(max(current[[varx]]));
  max.y.hold = list(max(current[[vart]]));
  names.hold = as.character(current[["name"]]);
  for(i in 2:length(arr.dens.ob)){
    current <- arr.dens.ob[[i]];
    addplot(current, var, col[i], type = type, w=w);
    #		max.x.hold = c(max.x.hold, max(current[[varx]]));
    #		max.y.hold = c(max.y.hold, max(current[[vart]]));
    names.hold = c(names.hold, as.character(current[["name"]]));
  };
  #	max.x <- mean(as.numeric(as.character(max.x.hold)));
  #	max.y <- mean(as.numeric(as.character(max.y.hold)));
  graphics::legend(l.pos, legend = as.character(names.hold), lty=1, lwd=2, cex=l.cex, col = col, box.col=NA);
};

#' Adds a single PDF plot to already open plot
#' 
#' Using an object from the cRacle::densform() or cRacle::and_fun() or cRacle::or_fun() functions, add a PDF over the existing plot. Useful for visualizing joint likelihood curves.
#' @param dens.ob An object derived from the cRacle::densform() or cRacle::and_fun() or cRacle::or_fun() functions.
#' @param var A character string that matches one of the layer names in the source raster object.
#' @param col A color declaration. Default is a random color.
#' @param type A character string of value either ".kde" for a Kernel Density Estimator curve, or ".gauss" for a Gaussian (normal) curve. All other values will result in errors.
#' @param w TRUE or FALSE to show weighted probability functions
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr); data(climondbioclim);
#' #extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' #dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' #multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' #or <- or_fun(dens.list.raw);
#' #addplot(or, names(climondbioclim[[1]]), col ='black');
#' #and <- and_fun(dens.list.raw);
#' #addplot(and, names(climondbioclim[[1]]), col ='black');
#' }

addplot <- function(dens.ob, var, col = sample(grDevices::colours()), type = ".kde", w=FALSE) {
  varx <- paste(var, "x", sep = ".");
  varw <- paste(var, "w", sep = ".");
  
  tempvarlist <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "MAT", "MaximumT", "MinimumT");
  if(var %in% tempvarlist){by = 10}else{by = 1};
  var <- paste(var, type, sep = "");
  if(w ==TRUE){
    to <- max(dens.ob[[varx]]);
    from <- min(dens.ob[[varx]]);
    num = length(dens.ob[[varx]]);
    lby = (to - from)/num;
    dens.ob[[var]] = dens.ob[[var]]^dens.ob[[varw]];
    den.area <- sum(dens.ob[[var]])*lby;
    dens.ob[[var]] = dens.ob[[var]]/den.area
  }
  graphics::points(dens.ob[[varx]]/by, dens.ob[[var]], type = "l", lwd = 3, col = col);
};
