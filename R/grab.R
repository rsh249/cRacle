#' @import jsonlite
#' @import urltools
#' @import rinat
NULL

#' Download distribution data directly from GBIF API
#'
#' This function requests data from the GBIF database for a single taxon using the GBIF callback API.
#'
#' @param taxon A string of the form 'genus species' or 'genus'.
#' @param maxrec Maximum number of records to download. 
#' Under 200 really doesn't mean anything because a single page (200) 
#' of results is returned and all records are kept.
#' @export
#' @examples \dontrun{
#' abies <- gbif_get('Abies');
#' }

gbif_get <- function(taxon, maxrec = 200000) {
  # require('jsonlite');
  #require('urltools')
  n = 0
  
  round = 0
  
  hold = list()
  
  offset = 0
  tori = taxon;
  taxon = urltools::url_encode(taxon)
  
  while (n < 1) {
    html_str = paste(
      "https://api.gbif.org/v1/occurrence/search?scientificName=",
      taxon,
      "&limit=300&offset=",
      offset,
      sep = ''
    )
    
    jsonget = jsonlite::fromJSON(html_str)
    
    round = round + 1
    
    if (is.null(nrow(jsonget$results))) {
      return(NULL)
    } else {
      hold[[round]] = jsonget$results
      
    }
    
    if (jsonget$endOfRecords == TRUE) {
      n = 1
      
    } else {
      offset = offset + 300
      
    }
    if (offset > maxrec) {
      break
      
    }
    
    ##DO filtering step:
    #exclude fossils
    #return(hold[[1]])
    foss = grep('FOSSIL_', hold[[round]]$basisOfRecord, ignore.case=TRUE);
    cult = grep('cultivat', hold[[round]]$locality, ignore.case=TRUE);
    gard = grep('gard', hold[[round]]$locality, ignore.case=TRUE);
    botan = grep('botan', hold[[round]]$locality, ignore.case=TRUE)
    if(length(foss) != 0){
  #    print(1)
      hold[[round]] = hold[[round]][-foss,]
    }
    if(length(cult) != 0){
    #  print(2)
      
      hold[[round]] = hold[[round]][-cult,]
    }
    if(length(gard) != 0){
     # print(3)
      
      hold[[round]] = hold[[round]][-gard,]
    }
    if(length(botan) != 0){
      #print(4)
      
      hold[[round]] = hold[[round]][-botan,]
      
    }
    
    
    
  }
  #return(hold)
  
  
  # cols = c('key',
  # 'genus',
  # 'specificEpithet',
  # 'decimalLatitude',
  # 'decimalLongitude');
  # 
    if(length(df)>0){
    df = cbind(hold[[1]]$key, rep(tori, nrow(hold[[1]])), hold[[1]]$decimalLatitude, hold[[1]]$decimalLongitude)
    } else {return(NULL);}  
    if(length(hold)>1){
      for(i in 2:length(hold)){
        nex = cbind(hold[[i]]$key, rep(tori, nrow(hold[[i]])), hold[[i]]$decimalLatitude, hold[[i]]$decimalLongitude)
        df = rbind(df, nex)
      }
    } 
  colnames(df) = c('ind_id', 'tax', 'lat', 'lon')
  df = as.data.frame(df)
  return(df)
  
#   if (sum(cols %in% names(hold[[1]])) == length(cols)) {
#     
#     df = hold[[1]][, c(cols)]
#     if (length(hold) > 1) {
#       for (n in 2:length(hold)) {
#         # print(n);
#         if (sum(cols %in% colnames(hold[[n]])) == length(cols)) {
#           nex = hold[[n]][, c(cols)]
#           df = rbind(df, nex)
#         } else {
#           
#           next
#         }
#         
#       } 
#       
#     }
#     
#     df[, 2] = paste(df[, 2], df[, 3])
#     
#     df = df[, -3]
#     
#     
#     colnames(df) = c('ind_id', 'tax', 'lat', 'lon')
#     #df$tax = rep(tori, nrow(df));
#     return(df)
#   } else {
#     
#     print("this")
# #    return(NULL)
#     
#   }
}

#' Download distribution data from repositories
#'
#' This function requests data from the GBIF database for a single taxon using the GBIF callback API.
#'
#' @param taxon A string of the form 'genus species' or 'genus'.
#' @param maxrec Maximum number of records to download.
#' @param repo Which data repositor(ies) should be searched. Accepts 'gbif', 'inat', and/or 'bison'. Defaults to 'gbif'.
#' @export
#' @examples \dontrun{
#' abies <- get_dist_all('Abies', maxrec = 1000);
#' }
get_dist_all <- function(taxon, maxrec = 19999, repo=c('gbif')) {
  ###GET DATA
  #GET GBIF DATA direct
  gbif = cbind(1,1,1,1);
  if('gbif' %in% repo | "GBIF" %in% repo){
  tryCatch({
      gbif <- gbif_get(taxon, maxrec = maxrec)
  },
  error = function(cond) {
    message(paste("GBIF", cond))
    return(NULL)
  })
  }
  
  
  #GET BIEN DATA
  bien = cbind(1,1,1,1);
  # tryCatch({
  # bien <-
  #   BIEN::BIEN_occurrence_species(
  #     species = taxon,
  #     native.status = TRUE,
  #     only.new.world = TRUE
  #   )
  # },
  # error = function(cond) {
  #   message(paste("BIEN", cond))
  #   return(NULL)
  # })
  
  #get bison data
  bison = cbind(1,1,1,1)
  if('bisont' %in% repo){
  tryCatch({
    bison <- cRacle::get_bison(taxon, maxrec = maxrec)
    
  },
  error = function(cond) {
    message(paste("BISON", cond))
    return(NULL)
  })
  }
  
  #get inaturalist data **research grade only**
  inatr = cbind(1,1,1,1);
  if('inat' %in% repo){
  tryCatch({
    inatr = cRacle::inat(taxon, maxrec = maxrec)
  }, 
  error = function(cond) {
    message(paste("inat", cond))
    return(NULL)
  })
  }
  
  
  cnames <- c('ind_id', 'tax', 'lat', 'lon')
  
  if (nrow(inatr) > 5) {
    inatr = inatr[, c('id', 'scientific_name', 'latitude', 'longitude')]
    colnames(inatr) = cnames
  } else {
    inatr = NA
  }
  if (nrow(bison) > 5) {
    bison = bison[, c('occurrenceID',
                      'name',
                      'decimalLatitude',
                      'decimalLongitude')]
    colnames(bison) = cnames
  } else {
    bison = NA
  }
  #if (nrow(gbif) > 5) {
    #gbif = gbif[, c('key',
    #                   'genus',
    #                   'specificEpithet',
    #                   'decimalLatitude',
    #                   'decimalLongitude')]
    #gbif[, 2] = paste(gbif[, 2], gbif[, 3], sep = ' ')
    #gbif = gbif[,-3]
    #colnames(gbif) = cnames
  #} else {
  #  gbif = NA
  #}
  
  if (nrow(bien) > 5) {
    bien = bien[, c('datasource_id',
                    'scrubbed_species_binomial',
                    'latitude',
                    'longitude')]
    colnames(bien) = cnames
  } else {
    bien = NA
  }
  data <- rbind(inatr, bison, gbif, bien) ## Consider using plyr::rbind.fill here
  if(nrow(data)<5){return(NULL)} #effectively nothing returned and catches a NULL error from the rbind above.
  data$lat <- as.numeric(as.character(data$lat))
  data$lon <- as.numeric(as.character(data$lon))
  #data = subset(data, data$tax == taxon)
  data$tax = rep(taxon, nrow(data))
  data = stats::na.omit(data)
  return(data)
}

#' Download distribution data, filter, and merge with climate or environmental
#'
#' getextr is a function that gets GBIF data and extracts climate or environmental 
#' data for each occurrence. This is a whole workflow for distribution 
#' data acquisition and value addition that draws on several other functions in cRacle
#' including gbif_get and extraction. Parallel option is useful for speeding up data collection for
#' many species when computational resources are available.
#' 
#' @param x A taxon name or list of taxon names. It is sometimes good to 
#' test these on the cRacle::get_gbif() function first.
#' @param maxrec Maximum number of records to download.
#' @param clim A raster object of climate or other environmental data to extract from.
#' @param repo Pass to get_dist_all
#' @param schema To be passed to cRacle::extraction
#' @param rm.outlier To be passed to cRacle::extraction
#' @param factor To be passed to cRacle::extraction
#' @param alpha To be passed to cRacle::extraction
#' @param nmin To be passed to cRacle::extraction
#' @param parallel TRUE or FALSE. Should this be executed in parallel.
#' @param nclus If parallel == TRUE then how many cores should be used? Default is 4.
#' 
#' @export
#' @examples \dontrun{
#' abies <- getextr(c('Abies fraseri', 'Abies lasiocarpa', 'Pinus strobus'), 
#' clim = clim, maxrec=500, 
#' schema= 'flat', rm.outlier = TRUE, 
#' alpha=0.01, factor = 2, nmin = 5, parallel=FALSE, nclus = 4);
#' }
#' 
getextr = function(x, clim = clim, maxrec=500, schema= 'flat', repo=c('gbif'),
                   rm.outlier = FALSE, alpha=0.01,  
                   factor = 4, nmin = 5, parallel=FALSE, nclus = 4){
  
  clim = clim;
  maxrec = maxrec;
  schema = schema;
  rm.outlier = rm.outlier;
  alpha = alpha;
  factor = factor;
  nmin = nmin;
  parallel = parallel;
  nclus = nclus;
  repo=repo;
  
  
  subfun = function(x){
    ex = list();
    for(i in 1:length(x)){
      print(x[i]);
      ex[[i]] = NULL;
      dat2 = cracle::get_dist_all(x[i], maxrec = maxrec, repo=repo)
      print(nrow(dat2));
      if(is.null(dat2)){ ex[[i]]=NULL; next; }
      dat2 = stats::na.omit(dat2);
      if(any(is.na(dat2))){ ex[[i]]=NULL; next;}
      if(nrow(dat2)<nmin){ ex[[i]]=NULL; next; }
      ex.hold = cRacle::extraction(dat2, clim, 
                                       schema = schema, 
                                       rm.outlier = rm.outlier, 
                                       alpha = alpha, 
                                       factor = factor, 
                                       nmin = nmin);
      if(length(ex.hold) == 0){ ex[[i]] = NULL;next;} else {
        ex.hold$tax = rep(x[i], nrow(ex.hold))
        ex[[i]] = ex.hold;
      }
    }
    
    ex = stats::na.omit(ex);
    #	if(any(is.null(ex))){ return(NULL); }
    if(length(ex) == 0) { return(NULL); }
    
    ex2 = plyr::rbind.fill(ex[[1]]);
    if(length(ex)>1){
      for(k in 2:length(ex)){
        ex2 = plyr::rbind.fill(ex2, ex[[k]]);
      }
    } else { return(ex[[1]]); }
    print(nrow(ex2))
    return(ex2);
  }
  
  
  if(parallel==FALSE){
    return(subfun(x));
  } else {
    clim = clim;
    maxrec = maxrec;
    schema = schema;
    rm.outlier = rm.outlier;
    alpha = alpha;
    factor = factor;
    nmin = nmin;
    parallel = parallel;
    nclus = nclus;
    
    cl = parallel::makeCluster(nclus, type = "SOCK", outfile = '')
    parallel::clusterExport(cl, varlist = c('clim',  'maxrec', 'nmin', 'schema', 'rm.outlier', 'alpha', 'factor' ), envir = environment())
    splits = parallel::clusterSplit(cl, x);
    extr = parallel::parLapply(cl, splits, subfun);
    parallel::stopCluster(cl);
    #return(extr);
    
    ##code below here not executed and problematic::
    extall = plyr::rbind.fill(extr[[1]]); ##Need to check that this object is OK as below.
    
    for(k in 2:length(extr)){
      
      if(is.null(extr[[k]])){} else {
        if(length(ncol(extr[[k]]))==0){} else {
          extall=plyr::rbind.fill(extall, extr[[k]]);
        }
      }
    }
    return(extall);
    
  }
}


##Add summary stats function calling stats::aggregate
#sum.stats = stats::aggregate(g$annualPET, by = list(g$tax), median); #for example


#' Download distribution data from iNaturalist
#' 
#' This function requests data from the iNaturalist database for a single species that
#' wraps the data access function(s) from the rinat library.
#' 
#' @param taxon A string of the form 'genus species'.
#' @param maxrec Limit on number of records to download.
#' @export
#' @examples \dontrun{
#' abies <- inat('Abies fraseri', 10000);
#' }

inat <- function(taxon, maxrec = 10000){
  #require(rinat);
  #8/17/2018 ERROR NOTE: Getting errors "your search returned too many results" even when setting maxresults to <100.
  di = rinat::get_inat_obs(taxon_name=taxon, maxresults=maxrec, quality='research')
  return(di);
}

#' Download distribution data from BISON (https://bison.usgs.gov)
#' 
#' This function requests data from the BISON database for a single species. 
#' Note that BISON requires exact name matching to binomial. Searching on a
#' genus, for example, will match only records that are enterred with that name
#' only, NOT all records in that genus identified by a binomial. Also note that
#' BISON has no name correction so mispelling and errors are likely.
#' 
#' @param taxon A string of the form 'genus species'.
#' @param maxrec Limit on number of records to download.
#' @export
#' @examples \dontrun{
#' abies <- get_bison('Abies fraseri', 10000);
#' }



get_bison <- function(taxon, maxrec = 10000) {
  #require('jsonlite');
  #require('urltools')
  taxon = urltools::url_encode(taxon)
  
  ###Run test query here to get total number of records so can adjust maxrecs
  html_test = paste("https://bison.usgs.gov/api/search.json?species=", 
                    taxon, "&type=scientific_name&start=0&count=1", sep = '');
  #print(html_test); 
  
  jsontest = jsonlite::fromJSON(html_test); 
  #return(jsontest)
  ntotal = sum(unlist(jsontest$occurrences$legend))
  if(maxrec>ntotal){maxrec=ntotal}
  dat = matrix(ncol=8, nrow=maxrec+100);
  dat = as.data.frame(dat);
  tdat = jsontest$data;
  colnames(dat)=colnames(tdat);
  if(maxrec<500){recs = maxrec}else{recs=500}
  if(maxrec>500){
    i = maxrec/500; print(i);
    for(z in 0:i){
      n=z*500;
      if(n>maxrec){break}
      if(n+500>maxrec){recs = maxrec-n}
      html_str = paste("https://bison.usgs.gov/api/search.json?species=", 
                       taxon, "&type=scientific_name&start=", n, "&count=", 
                       recs, sep = '');
      # print(html_str); 
      
      jsonget = jsonlite::fromJSON(html_str); 
      #return(jsonget)
      d = jsonget$data;
      #d2 = subset(d, dat$geo != "No");
      #print(d)
      # print(n)
      # print(n+nrow(d))
      # print(nrow(dat))
      dat[(n+1):(n+nrow(d)),]=d;
      if(class(dat)=='list' & n == 0){
        stop('no data from bison\n')
        return(NA)
      }
    }
    
  } else {
    n=0
    html_str = paste("https://bison.usgs.gov/api/search.json?species=", 
                     taxon, "&type=scientific_name&start=", n, "&count=", 
                     recs, sep = '');
    # print(html_str); 
    
    jsonget = jsonlite::fromJSON(html_str); 
    #return(jsonget)
    d = jsonget$data;
    #d2 = subset(d, dat$geo != "No");
    # print(d)
    dat[n+1:(n+nrow(d)),]=d;
    if(class(dat)=='list'){
      stop('no data from bison\n')
      return(NA)
    }
    
  }  
  
  ##Opportunity to provide filtering options here
  ##geo==TRUE
  
  ##fossil==FALSE
  
  
  return(dat)
}  







