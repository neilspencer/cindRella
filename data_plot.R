heldoutprobsunif <- rep(20000/sum(1- is.na(grid_memberships)), ntests)
# values for uniform contact surface

# create example plots comparing performance of different methods
df <- data.frame(index = rep(NA, 9 * ntests), performance = rep(NA, 9 * ntests), model = rep(NA, 9 * ntests))
order1 <- sort(heldoutprobs1, index.return = T)$ix

df$performance <- c(heldoutprobs1[order1], heldoutprobs2[order1], heldoutprobs3[order1], heldoutprobs4[order1], heldoutprobs5[order1], heldoutprobs6[order1], heldoutprobscontact[order1], heldoutprobskde[order1], heldoutprobsunif[order1])
df$index <- rep(c(1:ntests), 9)
df$model <- rep(c("Ours", "-kernel", "-w", "-scores", "-kernel/score","-contact", "Contact", "KDE", "Uniform"), each = ntests)


competitorsdf <- subset(df, model %in% c("Ours", "Contact", "KDE", "Uniform"))
variantsdf <- subset(df, !(model %in%  c("Contact", "KDE", "Uniform")))

competitorsplot <- ggplot(competitorsdf, aes(x = index, y = performance, color = model)) + geom_line()  +ylab("Heldout Probability per Accidental on Shoe") +xlab("Shoe Index (Sorted by Our Model's Performance)")
competitorsplot <- competitorsplot + theme_bw() # comparison of our method to competitor models

variantsplot <- ggplot(variantsdf, aes(x = index, y = performance, color = model)) + geom_line()  +ylab("Heldout Probability per Accidental on Shoe") +xlab("Shoe Index (Sorted by Our Model's Performance)")
variantsplot <- variantsplot + theme_bw() # comparison of all variants of our model

# plot the posterior predictive for the last shoe in the tests
predictdat <- dfposterior(fit1, testcontacts[1,], w_indices, zoptions, nzoptions, optiondiffsx, optiondiffsy, grid_memberships, chain_indices = 1:nsim)

predictplot <- ggplot(data = predictdat, aes(x, y)) + geom_raster(aes(fill = density)) + scale_fill_gradientn(name = "Posterior Density", colours=rev(brewer.pal(10,"Spectral"))) + theme_bw()
predictplot <- predictplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + xlab("h") + ylab("v")


