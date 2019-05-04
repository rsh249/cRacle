
#' P(A | B) = P(A) + P(B)
#' 
#' Using an object from the cRacle::dens_obj() function. Create a single density object (i.e., like that produced by cRacle::densform()) where the probability curves correspond to the probability density function of any one taxon/species from the original set occurring. This is not actually used in the implementation of finding the maximum joint likelihood in a CRACLE analysis, but is a good companion to the cRacle::and_fun() function.
#' @param dens.oblist An object derived from the cRacle::dens_ob() function.
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");

#' #OR:
#' data(distr);
#' data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' or <- or_fun(dens.list.raw);
#' addplot(or, names(climondbioclim[[1]]), col ='black');
#' }

or_fun <- function(dens.oblist){
  varlist <- names(dens.oblist[[1]]);
  varlist <- (varlist[1:((length(varlist)-1)/6)]);
  varlist <- sub(".kde", "", varlist);
  
  field <- list();
  gfield <- list();
  xfield <- list();
  meanadjust <- list();
  variances <- list();
  name = "ADDITION";
  for (n in 1:length(varlist)){
    var = varlist[n];
    varx <- paste(var, "x", sep = ".");
    vargauss <- paste(var, "gauss", sep = ".");
    varkde <- paste(var, "kde", sep = ".");
    
    varmean <- paste(var, "mean", sep = ".");
    varsd <- paste(var, "sd", sep = ".");
    meanlist <- list();
    sdlist <- list();
    dens.obcurr <- dens.oblist[[1]];
    to <- max(dens.obcurr[[varx]]);
    from <- min(dens.obcurr[[varx]]);
    num = length(dens.obcurr[[varx]]);
    by = (to - from)/num;
    meanlist[[1]] <- as.numeric(dens.obcurr[[varmean]]);
    sdlist[[1]] <- as.numeric(dens.obcurr[[varsd]])^2;
    prod <- as.numeric(dens.obcurr[[varkde]]);
    prod.gauss <- as.numeric(dens.obcurr[[vargauss]]);
    
    for(i in 2:length(dens.oblist)){
      dens.obnow <- dens.oblist[[i]];
      prod <- prod + (as.numeric(dens.obnow[[varkde]]));
      prod.gauss <- prod.gauss + (as.numeric(dens.obnow[[vargauss]]));
      prod.area <- sum(prod)*by;
      prod <- prod / prod.area;
      prod.gauss.area <- sum(prod.gauss)*by;
      prod.gauss <- prod.gauss / prod.gauss.area;
      meanlist[[i]] <- as.numeric(dens.obnow[[varmean]]);
      sdlist[[i]] <- as.numeric(dens.obnow[[varsd]])^2;
    };
    prod.area <- sum(prod)*by;
    prod <- prod / prod.area;
    prod.gauss.area <- sum(prod.gauss)*by;
    prod.gauss <- prod.gauss / prod.gauss.area;
    field[[n]] <- prod;
    gfield[[n]] <- prod.gauss;
    xfield[[n]] <- dens.obcurr[[varx]];
    meanadjust[[n]] <- as.numeric(meanlist)/as.numeric(sdlist);
    variances[[n]] <- 1/as.numeric(sdlist);
  };
  meansum <- lapply(meanadjust, sum);
  varisum <- lapply(variances, sum);
  wmeans <- mapply("/", meansum, varisum);
  wsd <- mapply("/", 1, varisum);
  wsd <- lapply(wsd, sqrt);
  field <- data.frame(field);
  gfield <- data.frame(gfield);
  xfield <- data.frame(xfield);
  colnames(field) <- (paste(varlist, "kde", sep = "."));
  colnames(gfield) <- paste(varlist, "gauss", sep = ".");
  names(wmeans) <- paste(varlist, "mean", sep = ".");
  names(wsd) <- paste(varlist, "sd", sep = ".");
  colnames(xfield) <- (paste(varlist, "x", sep = "."));
  name = data.frame(name);
  colnames(name) <- "name";
  fin <- c(field, xfield, gfield, wmeans, wsd, name);
  fin <- .makeaucone(fin);
  return(fin);
};

#and_fun = compiler::cmpfun(and_fun);

#' P(A | B) = P(A) * P(B)
#' 
#' Using an object from the cRacle::dens_obj() function. Create a single density object (i.e., like that produced by cRacle::densform()) where the probability curves correspond to the probability density function of ALL taxa/species from the original set occurring. 
#' @param dens.oblist An object derived from the cRacle::dens_ob() function.
#' @param w Weight importance of probability functions
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr);
#' data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' or <- or_fun(dens.list.raw);
#' addplot(or, names(climondbioclim[[1]]), col ='black');
#' and <- and_fun(dens.list.raw);
#' addplot(and, names(climondbioclim[[1]]), col ='black');
#' }

and_fun <- function(dens.oblist, w = FALSE){
  dens.oblist <- .scramble(dens.oblist);
  varlist <- names(dens.oblist[[1]]); #print(varlist)
  varlist <- (varlist[1:((length(varlist)-1)/6)]) ;
  varlist <- sub(".kde", "", varlist);
  
  field <- list();
  gfield <- list();
  xfield <- list();
  meanadjust <- list();
  variances <- list();
  name = "PRODUCT";
  for (n in 1:length(varlist)){ #print(varlist[n])
    var = varlist[n]; 
    varx <- paste(var, "x", sep = ".");
    varkde <- paste(var, "kde", sep = ".");
    
    vargauss <- paste(var, "gauss", sep = ".");
    varmean <- paste(var, "mean", sep = ".");
    varsd <- paste(var, "sd", sep = ".");
    varw <- paste(var, "w", sep = '.');
    meanlist <- list();
    sdlist <- list();
    
    dens.obcurr <- dens.oblist[[1]];
    
    if(w == TRUE){ 
      we = dens.obcurr[[varw]]; #print(we);
    } else {
      we <- 1;
    }
    
    to <- max(stats::na.omit(dens.obcurr[[varx]]));
    from <- min(stats::na.omit(dens.obcurr[[varx]]));
    num = length((dens.obcurr[[varx]]));
    by = (to - from)/num;
    meanlist[[1]] <- as.numeric(dens.obcurr[[varmean]]);
    sdlist[[1]] <- as.numeric(dens.obcurr[[varsd]])^2;
    prod <- as.numeric(dens.obcurr[[varkde]]) ^ we;
    prod <- prod*by;
    prod.gauss <- as.numeric(dens.obcurr[[vargauss]])*by;
    for(i in 2:length(dens.oblist)){
      dens.obnow <- dens.oblist[[i]];
      if(w == TRUE){
        we = dens.obnow[[varw]]; #print(we);
      } else {
        we = 1;
      }
      if(sum(stats::na.omit(dens.obnow[[varkde]]*by)) == 0) {next;}
      prod <- prod * ((as.numeric(dens.obnow[[varkde]])*by) ^ we);
      prod.area <- sum(stats::na.omit(prod))*by;
      prod <- prod / prod.area;
      prod.gauss <- prod.gauss * ((as.numeric(dens.obnow[[vargauss]])*by)^we);
      prod.gauss.area <- sum(stats::na.omit(prod.gauss))*by;
      prod.gauss <- prod.gauss / prod.gauss.area;
      meanlist[[i]] <- as.numeric(dens.obnow[[varmean]]);
      sdlist[[i]] <- as.numeric(dens.obnow[[varsd]])^2;
    };
    prod.area <- sum(stats::na.omit(prod))*by;
    prod <- prod / prod.area;
    prod.gauss.area <- sum(stats::na.omit(prod.gauss))*by;
    prod.gauss <- prod.gauss / prod.gauss.area;
    field[[n]] <- prod;
    gfield[[n]] <- prod.gauss;
    xfield[[n]] <- dens.obcurr[[varx]];
    meanadjust[[n]] <- as.numeric(unlist(meanlist))/as.numeric(unlist(sdlist));
    variances[[n]] <- 1/as.numeric(unlist(sdlist));
  };
  meansum <- lapply(meanadjust, sum);
  varisum <- lapply(variances, sum);
  wmeans <- mapply("/", meansum, varisum); ## Does this make sense for the weighted standard deviation calculation?
  wsd <- mapply("/", 1, varisum); ## wsd = sqrt(1(sum(1/sd)))???
  wsd <- lapply(wsd, sqrt);
  field <- data.frame(field);
  gfield <- data.frame(gfield);
  xfield <- data.frame(xfield);
  colnames(field) <- paste(varlist, "kde", sep = ".");
  colnames(gfield) <- paste(varlist, "gauss", sep = ".");
  names(wmeans) <- paste(varlist, "mean", sep = ".");
  names(wsd) <- paste(varlist, "sd", sep = ".");
  colnames(xfield) <- (paste(varlist, "x", sep = "."));
  name = data.frame(name);
  colnames(name) <- "name";
  fin <- c(field, xfield, gfield, wmeans, wsd, name);
  fin <- .makeaucone(fin);
  return(fin);
};



#get_optim() takes an object output from the densform function or and_fun or or_fun and finds optimal values for each PDF
#' Find PDF optim(a)um
#' 
#' Using an object from the cRacle::dens_obj() function. Create a single density object (i.e., like that produced by cRacle::densform()) where the probability curves correspond to the probability density function of ALL taxa/species from the original set occurring. 
#' @param dens.ob An object derived from the cRacle::dens_ob() function.
#' @export
#' @examples \dontrun{
#' #distr <- read.table('test_mat.txt', head=T, sep ="\t");
#' #OR:
#' data(distr);
#' data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' dens.list.raw <- dens_obj(extr.raw, clim = climondbioclim, bw = 'nrd0', n = 1024);
#' multiplot(dens.list.raw, names(climondbioclim[[1]]));
#' and <- and_fun(dens.list.raw);
#' addplot(and, names(climondbioclim[[1]]), col ='black');
#' optim.and <- get_optim(and);
#' abline(v=optim.and$means[paste(names(climondbioclim[[1]]), 'mean', sep = ".") ])
#' }

get_optim <- function(dens.ob){
  dens.ob1 <- dens.ob;
  varlist <- names(dens.ob1);
  varlist <- (varlist[1:((length(varlist)-1)/5)]);
  varlist <- sub(".kde", "", varlist);
  conintkde <- list();
  conintgauss <- list();
  dirconint <- list();
  origk <- list();
  origg <- list();
  means <- list();
  sds <- list();
  for (j in 1:length(varlist)){
    
    var = varlist[[j]];
    varx <- paste(var, "x", sep = ".");
    #print(varx)
    vargauss <- paste(var, "gauss", sep = ".");
    varkde <- paste(var, "kde", sep = ".");
    
    varmean <- paste(var, "mean", sep = ".");
    varsd <- paste(var, "sd", sep = ".");
    cumulkde <- list();
    cumulgauss <- list();
    cikde <- list(0,0);
    cigauss <- list(0,0);
    runkde <- 0;
    rungauss <- 0;
    to <- max(stats::na.omit(dens.ob1[[varx]]));
    from <- min(stats::na.omit(dens.ob1[[varx]]));
    num = length(stats::na.omit(dens.ob1[[varx]]));
    by = (to - from)/num;
    for (i in 1:length(dens.ob1[[varkde]])){
      runkde = runkde + (dens.ob1[[varkde]][i]*by);
      cumulkde[[i]] <- runkde;
      if(is.na(cumulkde[[i]])){cumulkde[[i]] = 0;}
      if(i==1){ 
        if(cumulkde[[i]] >= 0.025){
          cikde[[1]] <- dens.ob1[[varx]][i];
        };
        if(cumulkde[[i]] >= 0.975){
          cikde[[1]] <- dens.ob1[[varx]][i];
          cikde[[2]] <- dens.ob1[[varx]][i];
        };
      } else {
        if(cumulkde[[i-1]] < 0.025 && cumulkde[[i]] >= 0.025){
          cikde[[1]] <- dens.ob1[[varx]][i];
        };
        if(cumulkde[[i-1]] < 0.975 && cumulkde[[i]] >= 0.975){
          cikde[[2]] <- dens.ob1[[varx]][i];
        };
      };
      rungauss = rungauss + (dens.ob1[[vargauss]][i]*by);
      cumulgauss[[i]] <- rungauss;
      if(is.na(cumulgauss[[i]])){cumulgauss[[i]] = 0;}
      
      if(i==1){ 
        if(cumulgauss[[i]] >= 0.025){
          cigauss[[1]] <- dens.ob1[[varx]][i];
        };
        if(cumulgauss[[i]] >= 0.975){
          cigauss[[1]] <- dens.ob1[[varx]][i];
          cigauss[[2]] <- dens.ob1[[varx]][i];
        };
      } else {
        if(cumulgauss[[i-1]] < 0.025 && cumulgauss[[i]] >= 0.025){
          cigauss[[1]] <- dens.ob1[[varx]][i];
        };
        if(cumulgauss[[i-1]] < 0.975 && cumulgauss[[i]] >= 0.975){
          cigauss[[2]] <- dens.ob1[[varx]][i];
        };
      };
    };
    logkde <- ifelse(dens.ob1[[varkde]]>0, log(dens.ob1[[varkde]]*by), -Inf);
    #	print(varkde);
    origkde <- subset(dens.ob1[[varx]], logkde >= max(stats::na.omit(logkde))*1.01);
    origk[[j]] <- c(min(stats::na.omit(origkde)), max(stats::na.omit(origkde)));
    loggauss <- ifelse(dens.ob1[[vargauss]]>0, log(dens.ob1[[vargauss]]*by), -Inf);
    #	loggauss <- log(dens.ob1[[vargauss]]*by)
    origgauss <- subset(dens.ob1[[varx]], loggauss >= max(stats::na.omit(loggauss))*1.01); 
    origg[[j]] <- c(min(stats::na.omit(origgauss)), max(stats::na.omit(origgauss)));
    conintkde[[j]] <- c(cikde[[1]], cikde[[2]]);
    conintgauss[[j]] <- c(cigauss[[1]], cigauss[[2]]);
    dirconint[[j]] <- c((dens.ob1[[varmean]] - 1.96*dens.ob1[[varsd]]), (dens.ob1[[varmean]]+1.96*dens.ob1[[varsd]]));
    means[[j]] <- dens.ob1[[varmean]];
    sds[[j]] <- dens.ob1[[varsd]];
  };
  conintkde <- data.frame(conintkde);
  conintgauss <- data.frame(conintgauss);
  origk <- data.frame(origk);
  origg <- data.frame(origg);
  dirconint <- data.frame(dirconint);
  means <- data.frame(means);
  sds <- data.frame(sds);
  colnames(conintkde) <- paste(varlist, "cikde", sep = ".");
  colnames(conintgauss) <- paste(varlist, "cigauss", sep = ".");
  colnames(origk) <- paste(varlist, "origkde", sep = ".");
  colnames(origg) <- paste(varlist, "origgauss", sep = ".");
  colnames(dirconint) <- paste(varlist, "cidir", sep = ".");
  colnames(means) <- paste(varlist, "mean", sep = ".");
  colnames(sds) <- paste(varlist, "sd", sep = ".");
  ret <- list(conintkde, conintgauss, origk, origg, dirconint, means, sds);
  names(ret) <- c("conintkde", "conintgauss", "origk", "origg", "dirconint", "means", "sds");
  return(ret);
};


#makes area under any PDF curve equal 1 (good for standardizing curves to be compared). HIDDEN!
.makeaucone <- function(dens.ob1, var){ var <- names(dens.ob1);
var <- (var[1:((length(var)-1)/6 )]);
var <- sub('.kde', '', var)
for(i in 1:length(var)){
  varnow <- var[[i]];
  varx <- paste(var[[i]], "x", sep = ".");
  gauss <- paste(var[[i]], "gauss", sep = ".");
  kde <- paste(var[[i]], "kde", sep = ".");
  
  to <- max(subset(dens.ob1[[varx]], !is.na(dens.ob1[[kde]])));
  from <- min(subset(dens.ob1[[varx]], !is.na(dens.ob1[[kde]])));
  num <- length(subset(dens.ob1[[varx]], !is.na(dens.ob1[[kde]])));
  by = (to - from)/num;
  do <- sum(stats::na.omit(dens.ob1[[kde]]))*by;
  if(do == '0'){
    return(0)
  }
  do.gauss <- sum(stats::na.omit(dens.ob1[[gauss]]))*by;
  dens.ob1[[kde]] <- dens.ob1[[kde]]/do;
  dens.ob1[[gauss]] <- dens.ob1[[gauss]]/do.gauss;
};
return(dens.ob1);
};

#scramble() reorders pdfs. No real reason to do this as order does not matter for the operations being done here.
.scramble <- function(x, k=3) {
  x.s <- seq_along(x);
  y.s <- sample(x.s);
  idx <- unlist(split(y.s, (match(y.s, x.s)-1) %/% k), use.names = FALSE);
  x[idx];
};
