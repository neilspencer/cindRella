# extact all fits and evaluate the held out likelihood on the test set
nsamps <- nsim
Phi1 <- matrix(ncol = nsamps, nrow = ncategories)
Phi2 <- matrix(ncol = nsamps, nrow = ncategories)
Phi3 <- matrix(ncol = nsamps, nrow = ncategories)
Phi4 <- matrix(ncol = nsamps, nrow = ncategories)
Phi5 <- matrix(ncol = nsamps, nrow = ncategories)
kernelx1 <- matrix(ncol = nsamps, nrow = Kx + 1)
kernely1 <- matrix(ncol = nsamps, nrow = Ky + 1)
kernelx3 <- matrix(ncol = nsamps, nrow = Kx + 1)
kernely3 <- matrix(ncol = nsamps, nrow = Ky + 1)
kernelx4 <- matrix(ncol = nsamps, nrow = Kx + 1)
kernely4 <- matrix(ncol = nsamps, nrow = Ky + 1)
kernelx6 <- matrix(ncol = nsamps, nrow = Kx + 1)
kernely6 <- matrix(ncol = nsamps, nrow = Ky + 1)
q1 <- fit1$q
q2 <- fit2$q
q3 <- fit3$q
q6 <- fit6$q

W1 <- matrix(ncol = nsamps, nrow = length(fit1$w[[1]]))
W2 <- matrix(ncol = nsamps, nrow = length(fit2$w[[1]]))
W4 <- matrix(ncol = nsamps, nrow = length(fit4$w[[1]]))
W5 <- matrix(ncol = nsamps, nrow = length(fit5$w[[1]]))
W6 <- matrix(ncol = nsamps, nrow = length(fit6$w[[1]]))


for(i in 1:nsamps){
  W1[,i] <- fit1$w[[i]]
  W2[,i] <- fit2$w[[i]]
  W4[,i] <- fit4$w[[i]]
  W5[,i] <- fit5$w[[i]]
  W6[,i] <- fit6$w[[i]]
  Phi1[,i] <- fit1$phi[[i]]
  Phi2[,i] <- fit2$phi[[i]]
  Phi3[,i] <- fit3$phi[[i]]
  Phi4[,i] <- fit4$phi[[i]]
  Phi5[,i] <- fit5$phi[[i]]
  kernelx1[,i] <- fit1$kernelx[[i]]
  kernely1[,i] <- fit1$kernely[[i]]
  kernelx3[,i] <- fit3$kernelx[[i]]
  kernely3[,i] <- fit3$kernely[[i]]
  kernelx4[,i] <- fit4$kernelx[[i]]
  kernely4[,i] <- fit4$kernely[[i]]
  kernelx6[,i] <- fit6$kernelx[[i]]
  kernely6[,i] <- fit6$kernely[[i]]
}

heldoutprobs1 <- rep(NA, ntests)
heldoutprobs2 <- rep(NA, ntests)
heldoutprobs3 <- rep(NA, ntests)
heldoutprobs4 <- rep(NA, ntests)
heldoutprobs5 <- rep(NA, ntests)
heldoutprobs6 <- rep(NA, ntests)

for(i in 1:ntests){
  shoen <- i
  zhomen <- zhometest[testdat$shoenum == test[shoen]]
  contactn <- testcontacts[shoen,]

  xx <- heldoutprob_importance(contactn, ncategories, W1, Phi1, Kx, Ky, kernelx1,  kernely1, q1, w_indices, zhomen, zoptions, nzoptions, optiondiffsx,  optiondiffsy, nsim_per = 1)
  heldoutprobs1[i] <- exp((logSumExp(c(xx)) - log(length(c(xx))))/length(zhomen)) * 20000/length(contactn)

  xx <- heldoutprob_importance_noZ(contactn, ncategories, W2, Phi2, q2, w_indices, zhomen, nsim_per = 1)
  heldoutprobs2[i] <- exp((logSumExp(c(xx)) - log(length(c(xx))))/length(zhomen)) * 20000/length(contactn)

  xx <- heldoutprob_importance_now(contactn, ncategories,  Phi3, Kx, Ky, kernelx3,  kernely3, q3, zhomen, zoptions, nzoptions, optiondiffsx,  optiondiffsy, nsim_per = 1)
  heldoutprobs3[i] <- exp((logSumExp(c(xx)) - log(length(c(xx))))/length(zhomen)) * 20000/length(contactn)

  xx <- heldoutprob_importance_noepsilon(contactn, ncategories, W4, Phi4, Kx, Ky, kernelx4,  kernely4, w_indices, zhomen, zoptions, nzoptions, optiondiffsx,  optiondiffsy)
  heldoutprobs4[i] <- exp((logSumExp(c(xx)) - log(length(c(xx))))/length(zhomen)) * 20000/length(contactn)

  xx <- heldoutprob_importance_noepsilonnoZ(contactn, ncategories, W5, Phi5, w_indices, zhomen)
  heldoutprobs5[i] <- exp((logSumExp(c(xx)) - log(length(c(xx))))/length(zhomen)) * 20000/length(contactn)

  xx <- heldoutprob_importance_noPhi(contactn, W6, Kx, Ky, kernelx6,  kernely6, q6, w_indices, zhomen, zoptions, nzoptions, optiondiffsx,  optiondiffsy, nsim_per = 1)
  heldoutprobs6[i] <- exp((logSumExp(c(xx)) - log(length(c(xx))))/length(zhomen)) * 20000/length(contactn)
}


heldoutprobscontact<- rep(NA, ntests)
for(i in 1:ntests){
  zs <- zhometest[testdat$shoenum == test[i]]

  probs <- phis[testcontacts[i,] + 1]/sum(phis[testcontacts[i,] + 1])
  heldoutprobscontact[i] <- exp(sum(log(probs * 20000)[zs])/length(zs))
}


heldoutprobskde <- rep(NA, ntests)
for(i in 1:ntests){
  xnews <- xtest[testdat$shoenum == test[i]]
  ynews <- ytest[testdat$shoenum == test[i]]
  kderess <- c()
  for(j in 1:length(xnews)){
    kderess[j] <- log(gridpredicts[floor(xnews[j]), floor(ynews[j])])
  }
  heldoutprobskde[i] <- exp(sum(kderess)/naccidentalstest[i])
}

