
# demonstration of constructing the necessary objects
# for running the algorithm

data("contacts_example") # load the example contact surfaces from the package

nshoes <- length(unique(contacts_example$shoenum)) # count the number of examples

ngridx <- 100 # the number of horizontal gridpoints
ngridy <- 200 # the number of vertical gridpoints
ngrid <- ngridx * ngridy # the total number of gridpoints

contactmat_example <- matrix(0, nrow = nshoes, ncol = ngrid) # to be a matrix to
# that conveniently stores all contact surfaces as {0,1} row vectors

for(i in 1:nshoes){
  shoenow <- subset(contacts_example, shoenum == i)
  gridnums <- shoenow$x + (shoenow$y - 1) * ngridx # the gridpoints with contact
  contactmat_example[i, gridnums] <- 1
}

contactcats_example <- matrix(0, nrow = nshoes, ncol = ngrid) # to be a matrix to
# that conveniently stores the phi value for all contact {0,1,...32} row vector
# entries


# fill in the contact categories being careful of not overrunning edges
for(i in 1:nrow(contactmat_example)){
  for(j in 1:ncol(contactmat_example)){
    if(j > ngridx){
      contactcats_example[i,j] <- contactcats_example[i,j] + contactmat_example[i, j - ngridx]
    }
    if(j %% ngridx != 1){
      contactcats_example[i,j] <- contactcats_example[i,j] + 2 * contactmat_example[i, j - 1]
    }

    contactcats_example[i,j] <- contactcats_example[i,j] + 4 * contactmat_example[i,j]

    if(j %% ngridx != 0){
      contactcats_example[i,j] <- contactcats_example[i,j] + 8 * contactmat_example[i, j + 1]
    }

    if(j <= (ngrid - ngridx)){
      contactcats_example[i,j] <- contactcats_example[i,j] + 16 * contactmat_example[i, j + ngridx]
    }

  }
}


data("grid_memberships") # the w checkerboard as defined in the paper


# Iterates over which gridpoints are actually considered by our prior,
# indexes them in w_memberships
kk <- 1
for(i in 1:(ngrid)){
  if(grid_memberships[i] != 0){
    grid_memberships[i] <- kk
    kk <- kk + 1
  }
  else{
    grid_memberships[i] <- NA
  }
}
ngrid2 <- length(unique(grid_memberships)) - 1 # number of actually considered
# gridpoints

# for each gridpoints, determine the possible auxiliary variable Z assignments
# corresponding to an accidental occuring there
Kx <- 3 # kernel width in x direction
Ky <- 3 # kernel width in y direction
zoptions <- matrix(NA, nrow = ngrid2, ncol = (2*Kx + 1) * (2*Ky + 1))
# the x offset for each z
optiondiffsx <- matrix(NA, nrow = ngrid2, ncol = (2*Kx + 1) * (2*Ky + 1))
# the y offset for each z
optiondiffsy <- matrix(NA, nrow = ngrid2, ncol = (2*Kx + 1) * (2*Ky + 1))
nzoptions <- rep(NA, ngrid2) #total number of options

kk <- 1
for(a in 1:(ngrid)){
  indx <- 1
  if(is.na(grid_memberships[a]) == F){
    for(x in (-Kx:Kx)){
      if((((a - 1)%%ngridx + 1) + x > 0) & (((a - 1)%%ngridx + 1) + x <= ngridx)){
        for(y in (-Ky:Ky)){
          if((ceiling(a/ngridx) + y > 0) & (ceiling(a/ngridx) + y <= ngridy)){
            if(is.na(grid_memberships[a + x + ngridx * y]) == F){
              zoptions[kk, indx] <- grid_memberships[a + x + ngridx * y] - 1
              optiondiffsy[kk, indx] <- abs(y)
              optiondiffsx[kk, indx] <- abs(x)
              indx <- indx + 1
            }
          }
        }
      }
    }
    nzoptions[kk] <- sum(!is.na(zoptions[kk,]))
    kk <- kk + 1
  }
}


# create the coarsened grid of our prior for w

nxsquares <- 10
nysquares <- 20
nsquares <- nxsquares * nysquares

# create a checkerboard of J's which each J contains a 10 by 10 grid
w_indices <- rep(rep(1:nxsquares, each = ngridx/nxsquares), ngridy) + nxsquares * rep(0:(nysquares - 1), each = ngrid/nysquares)
w_indices <- w_indices[is.na(grid_memberships) == F]


# create the correlations for the logNormal prior on w
adja1 <- 0 * diag(nsquares)
adja2 <- 0 * diag(nsquares)
for(i in 2:nsquares){
  if((i %% nxsquares) != 0){
    adja1[i,i-1] <- 1
    adja1[i-1,i] <- 1
  }
  if((i > nxsquares)){
    adja2[i, i - nxsquares] <- 1
    adja2[i - nxsquares, i] <- 1
  }
}

Ax <- adja1[which((1:nsquares) %in% w_indices),]
Ax <- Ax[,which((1:nsquares) %in% w_indices)]
Ay <- adja2[which((1:nsquares) %in% w_indices),]
Ay <- Ay[,which((1:nsquares) %in% w_indices)]

w_indices <- as.numeric(as.factor(w_indices)) - 1 # convert to 1= 0 for C++





