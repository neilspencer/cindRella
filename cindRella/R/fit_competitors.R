fit_contact <- function(traincontacts, naccidentalstrain, traindat, zhometrain, ncategories){
  ntype <- matrix(nrow = ncategories, ncol = ntrains)
  counts_per_shoe <- naccidentalstrain
  count_on_js <- rep(0, ncategories)
  for(i in 1:ntrains){
    zs <- zhometrain[traindat$shoenum == train[i]]
    for(j in 1:ncategories){
      ntype[j,i] = sum(traincontacts[i,] == (j-1))
      count_on_js[j] <- count_on_js[j] + sum(traincontacts[i, zs] == (j-1))
    }
  }


  tominimize <- function(logphi, phis, j, count_on_js, naccidentals, ntype, nshoes){
    phis[j] <- exp(logphi)
    res <- count_on_js[j]*logphi
    for(i in 1:nshoes){
      res <- res - naccidentals[i] * log(sum(phis * ntype[,i]))
    }
    return(-res)
  }


  phis <- runif(32)
  for(qq in 1:30){
    for(j in 1:32){
      xx <- optim(log(phis[j]), tominimize, gr = NULL, phis, j=j, count_on_js, naccidentalstrain, ntype, ntrains)
      phis[j] <- exp(xx$par)
    }
  }
  return(phis)
}

fit_kde <- function(naccidentalstrain, traindat, grid_memberships){
  thekde <- kde2d(100 * traindat$x, 200 * traindat$y, n = c(100,200), lims = c(0,100,0,200))
  gridpredicts <- 20000 * (thekde$z*(1 - is.na(grid_memberships)))/sum(thekde$z * (1 - is.na(grid_memberships)))
  return(gridpredicts)
}

