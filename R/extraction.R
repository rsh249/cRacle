#' @import grDevices
#' @import doParallel
#' @import foreach
#' @import iterators
#' @import doSNOW
#' @useDynLib cRacle
NULL



#' Extract environmental data
#'
#' This function is a feature-added wrapper for raster::extract();
#' @param data Distribution data. A data.frame that should include (at least) a column for species/taxon name named 'tax', latitude named 'lat', longitude named 'lon', and an optional column 'sub' if it is necessary to define subgroups (i.e., if 'tax' corresponds to genera but sampling needs to know which records belong to which species. See param 'schema' below).
#' @param clim A raster object (see raster::raster() and raster::stack() documentation for reading raster files into R).
#' @param schema A string of value "raw", "flat", or "species" to define the sampling protocol. In "raw", all records are counted (including duplicate exact localities). In "flat", all unique localities will be counted where a unique locality is defined as a raster grid cell. Under the "flat" sampling strategy two records in the same raster grid cell will be counted as one. The option "species", only applies when taxa are identified as genera and species identities are represented in the "sub" column of the data object. In "species", each unique locality is counted for each species within the group (taxon). This weighs more diverse localities higher. Default is "raw".
#' @param factor An integer value for the methods "flat" and "spec" to increase the systematic sampling grid size to courser resolutions than the given climate grid. The value of factor corresponds to the number of rows and columns to aggregate into each courser grid cell. Default is 0 which will not be processed.
#' @param rm.outlier TRUE or FALSE. Indicate whether to remove points that are climatic outliers for at least one variable given a normal 95 percent confidence interval.
#' @param alpha Confidence level (i.e., 0.05) for clipping out outlier records.
#' @param nmin Minimum number of records allowed. Taxa or groups with fewer records will not be returned.
#' @export
#' @examples
#' data(distr);
#' data(climondbioclim);
#' extr.raw = extraction(data=distr, clim= climondbioclim, schema='raw');
#' extr.flat = extraction(data=distr, clim= climondbioclim, schema='flat');
#' extr.spec = extraction(data=distr, clim= climondbioclim, schema='species');
extraction <- function(data, clim, schema = "raw", factor = 0, rm.outlier = FALSE,  alpha = 0.01, nmin = 5){

  if(nrow(data) < 5){cat('ERR: Too few records\n'); return(NULL);}

  mat.larr <- data;
  phytoclim <- clim;

  extr.larr <- terra::extract(phytoclim, cbind(mat.larr$lon, mat.larr$lat), cells=T);
  count = nrow(stats::na.omit(extr.larr))
  if(count<5){
    cat("ERR: Records out of study area\n")
    return(NULL)
  }


  extr.larr <- cbind(mat.larr, extr.larr);

  if(schema != 'raw'){
    if(factor == 0){} else {
      r2 <- terra::aggregate(phytoclim, fact = factor, fun="mean");
      tmp.ext <- terra::extract(r2, cbind(mat.larr$lon, mat.larr$lat), cells=T);
      extr.larr[,"rawcell"] = extr.larr[,'cell'];
      extr.larr[,'cell'] = tmp.ext[,'cell'];

    }


  }
  extr.larr <- stats::na.omit(extr.larr)
  if(schema == "raw"){
    holder <- data.frame();
    tlist <- unique(extr.larr$tax);
    for(i in 1:length(tlist)){
      set <- subset(extr.larr, extr.larr$tax == tlist[i]);

      #if(length(set[,1])>=5){
      holder <- rbind(holder, set);
      #print(length(holder[,1]))

      #}
    }
    extr.larr = holder;
  } else {
    holder <- data.frame();
    tlist <- unique(extr.larr$tax);
    for(i in 1:length(tlist)){
      set <- subset(extr.larr, extr.larr$tax == tlist[i]);
      if(schema == "flat"){
        sub = set

        sub <- sub[!duplicated(sub[,"cell"]),];
        #if(length(sub[,1])>=5){
        holder <- rbind(holder, sub);
        #	}
      }
      if(schema == "species"){
        glist <- unique(set$sub);
        for(n in 1:length(glist)){
          sub <- subset(set, set$sub == glist[n]);
          sub <- sub[!duplicated(sub[,"cell"]),];
          if(length(sub[,1])>=5){
            holder <- rbind(holder, sub);
          }
        }
      }
    }
    extr.larr <- holder;
  }
  #if(schema != 'raw'){

  #  extr.larr[,'cell'] = extr.larr[,ncol(extr.larr)];

   # extr.larr = extr.larr[,-ncol(extr.larr)];

  #}
  #print("EXTRACTION MONITOR:")
  #  print(length(holder[,1]));

  #print(length(extr.larr[,1]));
  head = ncol(mat.larr)+1
  #print(head)

  #extr.larr[,1] = as.numeric(as.character(extr.larr[,1]))
  if(rm.outlier== TRUE){
    df = as.data.frame(extr.larr)
    for(nn in 1:terra::nlyr(phytoclim)){
      print(head+nn)

      n.mean <- mean(as.numeric(df[,head+nn]), na.rm=T);
      n.sd <- stats::sd(as.numeric(df[,(head+nn)]));
      rn <- nrow(df);
      t = stats::qt((1-(alpha/2)), rn-1);
      minci = n.mean-(t*n.sd);
      maxci = n.mean+(t*n.sd);
      extr.larr <- subset(df, df[,(head+nn)] >= minci);
      extr.larr <- subset(df, df[,(head+nn)] <= maxci);

    }

  }
  extr.larr = df
  t.list = unique(extr.larr$tax);

  if(length(t.list)>1){
    hold = list()

    for(zz in 1:length(t.list)){
      sub <- subset(extr.larr, extr.larr$tax == t.list[[zz]]);
      sub <- stats::na.omit(sub)
      if(nrow(sub) < nmin){

      } else {
        hold[[zz]] = sub
      }
    }
    #colnames(hold) = colnames(extr.larr);
    hold = dplyr::bind_rows(hold)
  } else {
    hold = extr.larr;
  }


  if(nrow(hold) == 0){return(NULL)}
  return(hold);
};


