set.seed(1)
ntrains <- 20

data("accidentals_example") # load the example accidentals

# Set up training and testing
ords <- sample(1:nshoes)
train <- sort(ords[1:ntrains])
test <- sort(ords[(ntrains + 1):nshoes])
ntests <- length(test)

traindat <- subset(accidentals_example, shoenum %in% train)
testdat <- subset(accidentals_example, shoenum %in% test)
traincontacts <- contactcats_example[train, is.na(grid_memberships) == F]
testcontacts <- contactcats_example[test, is.na(grid_memberships) == F]

xtrain <- traindat$x
ytrain <- traindat$y
xtest <- testdat$x
ytest <- testdat$y

zhometrain <- grid_memberships[xtrain + ngridx * (ytrain - 1)] - 1
zhometest <- grid_memberships[xtest + ngridx * (ytest - 1)] - 1


naccidentalstrain <- ddply(traindat, ~shoenum, nrow)$V1 # number of accidentals per shoe
naccidentalstest <- ddply(testdat, ~shoenum, nrow)$V1 # number of accidentals per shoe

naccidentalsalltrain <- nrow(traindat)
naccidentalsalltest <- nrow(testdat)

ncategories <- 32
phi_init <- (1:ncategories)/ncategories

# parameters of distribution of q
qa <- 2
qb <- 2

# covariances in x and y direction
rhox_init <- 0.2
rhoy_init <- 0.2

shoenums <- as.numeric(as.factor(traindat$shoenum)) - 1
phi_init <- (1:ncategories)/ncategories # initialize phi
q_init <- 1.0 # initialize q

# initialize px and py
sigma <- 2
py_init <- rnorm(Ky + 1, 0, sigma)
px_init <- rnorm(Kx + 1, 0, sigma)

#initialize w
nws <- length(unique(w_indices))
w_init <- exp(log(table(c(w_indices[zhometrain + 1] + 1, 1:nws)))-4)

naccidentals <- rep(NA, nshoes)
second_params <- rep(NA, nshoes)
for(i in 1:ntrains){
  second_params[i] <- sum(q_init * w_init[w_indices + 1] * phi_init[1 + traincontacts[i,]])
}
u_init <- rgamma(ntrains, naccidentalstrain, second_params)
z_init <- rep(0, length(xtrain))
for(i in 1:length(xtrain)){
  z_init[i] <- sample(0:(nzoptions[zhometrain[i] + 1] - 1),1)
}



nsim <- 100 # chains have just 100 iterations for demonstration.

# fit our model
fit1 <- mcmc_inference(shoenums, traincontacts, ncategories, ntrains, z_init, w_init, w_indices, zhometrain, zoptions, nzoptions, optiondiffsx, optiondiffsy, phi_init, Ax, Ay,
                          Kx, Ky, px_init, py_init, sigma,  u_init,  q_init, qa, qb, rhox_init, rhoy_init, nsim, print_output = 0)
# fit our model with no kernel
fit2 <- mcmc_inference_noZ(shoenums, traincontacts, ncategories, ntrains, w_init, w_indices, zhometrain, phi_init, Ax, Ay, u_init,
                   q_init, qa,  qb, rhox_init, rhoy_init, nsim, print_output = 0)

# fit our model with uniform w
fit3 <- mcmc_inference_now(shoenums, traincontacts,ncategories, ntrains, z_init, zhometrain, zoptions, nzoptions, optiondiffsx, optiondiffsy,
                           phi_init, Ax, Ay, Kx, Ky, px_init, py_init, sigma,
                           u_init, q_init, qa, qb, nsim, print_output = 0)

# fit our model with no scores
fit4 <- mcmc_inference_noepsilon(shoenums,traincontacts, ncategories, ntrains, z_init, w_init,
                         w_indices, zhometrain, zoptions, nzoptions,  optiondiffsx, optiondiffsy,
                         phi_init, Ax, Ay, Kx, Ky, px_init, py_init, sigma, rhox_init,
                         rhoy_init, nsim, print_output = 0)

# fit our model without scores or a kernel
fit5 <- mcmc_inference_noepsilonnoZ(shoenums, traincontacts, ncategories, ntrains, w_init, w_indices, zhometrain, phi_init, Ax, Ay, rhox_init,
                                    rhoy_init, nsim, print_output = 0)

# fit our model with uniform Phi
fit6 <- mcmc_inference_noPhi(shoenums, ntrains, z_init, w_init, w_indices, zhometrain,
                             zoptions, nzoptions, optiondiffsx, optiondiffsy,
                             Ax, Ay, Kx, Ky, px_init, py_init, sigma, u_init, q_init, qa, qb,
                             rhox_init, rhoy_init, nsim ,  print_output = 0)

#fit the simple contact model via optimization (ignore warnings)
suppressWarnings(exprphis <- fit_contact(traincontacts, naccidentalstrain, traindat, zhometrain, ncategories))

# fit the kde model
gridpredicts <- fit_kde(naccidentalstrain, traindat, grid_memberships)
