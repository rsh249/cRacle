#!/usr/bin/R

#####
#Demo CRACLE implementation
####

#for a single locality with the species:
data(climondbioclim) #packaged data from CliMond for the Eastern US. Should replace with climate grids of your choice
sp_list = c(
  "Quercus rubra",
  "Pinus strobus",
  "Vaccinium angustifolium",
  "Betula papyrifera",
  "Rhododendron maximum",
  "Toxicodendron radicans"
)

extr = getextr(
  sp_list,
  climondbioclim,
  maxrec = 200,
  #Should set higher in practice!
  schema = 'flat',
  factor = 2,
  #adjust to spatially thin. Use a higher number (4 to 8) for spatial grids of 2.5 arcminutes or less.
  nmin = 10,
  parallel = FALSE
)

densall = dens_obj(extr,
                   climondbioclim,
                   manip = 'condi',
                   kern = 'gaussian',
                   bg.n = 40) #consider using parallel=TRUE if bg.n>500

and = and_fun(densall)
optim = get_optim(and)

print(optim$origk) #As implemented in Harbert and Nixon, 2015
print(optim$conintkde) #with 95% confidence intervals on the likelihood distribution.

n = 1
multiplot(densall, names(climondbioclim)[[n]])
addplot(and, names(climondbioclim)[[n]])

abline(v = median(optim$origk[, n]), col = 'green') #Midpoint will be very close to the optimum. Original method of Harbert and Nixon, 2015 returned a top 1% range to introduce some slop
