# --------------- STEP 4: Cross validation results ----------------- # (may 24th start)
repetitions <- 2
# Getting the results together (held out indicies and predictions)
all_indices <- array(NA, dim = c(repetitions, 50, 2))
our_preds <- array(NA, dim = c(repetitions, nB, nP))

for (rr in 1 : repetitions) {
  load(paste0(result_path, 'cv_indices.may.28.ob1_', rr, '.dat'))
  load(paste0(result_path, 'pred.may.28.ob1_', rr, '.dat'))
  all_indices[rr, , ] <- cv_indices
  our_preds[rr, , ] <- pred
}

# Predictions of the held out data from both models:
pred <- array(NA, dim = c(repetitions, 50, 2))
for (rr in 1 : repetitions) {
  for (ii in 1 : 50) {
    pred[rr, ii, 1] <- our_preds[rr, all_indices[rr, ii, 1], all_indices[rr, ii, 2]]
    # pred[rr, ii, 2] <- alt_preds[rr, all_indices[rr, ii, 1], all_indices[rr, ii, 2]]
  }
}
# may 24th fin.


### BELOW IS GEORGIA'S CODE (for plotting)
# Average and median probability of interaction based on the two models in the
# overall data:
overall_mean <- cbind(apply(our_preds, 1, mean), apply(our_preds, 1, mean))
overall_median <- cbind(apply(our_preds, 1, median), apply(our_preds, 1, median))

# Average and median in the held out data.
pred_mean <- apply(pred, c(1, 3), mean)
pred_median <- apply(pred, c(1, 3), median)

# Creating the data frame we will plot:
plot_dta <- abind::abind(pred_mean / overall_mean, pred_median / overall_median, along = 3)
dimnames(plot_dta)[2 : 3] <- list(method = c('Latent Factors', 'Covariates'),
                                  stat = c('Prediction mean / Overall mean',
                                           'Prediction median / Overall median'))
names(dimnames(plot_dta)) <- c('Iteration', 'Method', 'Statistic')
plot_dta <- reshape2::melt(plot_dta)

plot_dta2 <- plot_dta %>% filter(Method == "Latent Factors")

# Plotting cross validation results:
ggplot(data = plot_dta2) +
  geom_boxplot(aes(x = Method, y = value)) +
  facet_wrap(~ Statistic, scales = 'free_y') +
  theme_bw() +
  ylab('') +
  xlab('') +
  ggtitle('Out of sample performance', subtitle = 'Predictions for held-out recorded interactions') +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = function(x) c(1, x[2]), n.breaks = 6)
