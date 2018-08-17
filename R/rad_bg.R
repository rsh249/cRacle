#' Generate backround data within a radius of occurrence records
#' 
#' This function generates a sample of background coordinates that are within n kilometers of one occurrene record up to x background points per record.
#' @param coords A two-column matrix or data.frame of coordinates with the first column being longitude and the second being latitude.
#' @param clim A raster object to extract background data from.
#' @param radius A distance in km to define a radius around each occurrence to sample from.
#' @param n Number of background points to generate per occurrence point. Note that background points will be thinned so that only one occurs per grid cell, so the total number returned will be less than n * length(coords[,1]).
#' @export
#' @examples \dontrun{
#' data(distr);
#' bg.ext <- rad_bg(distr[,4:3], climondbioclim, radius=100, n = 50)
#' }

rad_bg <- function(coords, clim, radius, n){
  extr = coords;
  ##Premise: IF a species PDF - P(sp) - is approximately equal to the background PDF
  #P(bg) then P(condi) = P(sp)/P(bg) should be a uniform distribution and
  #not contribute to the overall estimation of climate.
  
  bg.mat <- matrix(ncol = 2, nrow = n * nrow(extr));
  for(i in 1:length(extr[,1])){
    # print(i)
    for(zz in 1:n){
      # print(paste('zz', zz))
      dir = sample(1:360, 1)
      dist = sample(1:radius, 1);
      bg.mat[(i*zz),1:2] = findcoord(extr[i,1], extr[i,2], dist, dir)
      
    }
  }
  bg.mat = cbind(rep(1111, length(bg.mat[,1])), rep('bg', length(bg.mat[,1])), bg.mat[,2], bg.mat[,1])
  #  bg.mat = stats::na.omit(bg.mat)
  colnames(bg.mat) = c('ind_id', 'tax', 'lat', 'lon')
  bg.mat <- data.frame(bg.mat)
  bg.mat$lon = as.numeric(as.character(bg.mat$lon))
  bg.mat$lat = as.numeric(as.character(bg.mat$lat))
  # bg.mat = unique(bg.mat);
  bg.ext <- raster::extract(clim, bg.mat[,4:3], cellnumbers=T)
  bg.ext <- cbind(bg.mat, bg.ext);
  bg.ext <- stats::na.omit(bg.ext);
  
  bg.ext = bg.ext[!duplicated(bg.ext[,"cells"]),]
  
  return(bg.ext);
  
}
