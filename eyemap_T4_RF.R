# UPDATE 2021-09-24,  use T4_v2.R for PDs, load "data/T4_gallery.RData"
# this for interpolation of T4_RF at all retina locations

# UPDATE 2021-01-12
# test lm vs npreg vs npreg_1D --> use npreg
# med_xyz nbhd 9um --> 20um

# UPDATE 2021-01-05
# no aux for med->lens

# UPDATE 2022-07-25
# use PD calculated from eyemap_T4_v2

# When SAVE, switch from T4 orientation to PD

# init --------------------------------------------------------------------

library(natverse)
library(tidyverse)
library(RColorBrewer)
library(np) # kernel regression to find the base pt on the eye

# clean everythign up.
rm(list=ls())
#close any open rgl windows
while (rgl.cur() > 0) { rgl.close() }
while (dev.cur() > 1) { dev.off() }

source("eyemap_func.R")

load('data/eyemap.RData')
load('data/hexnb_ind_dist.RData')

load("data/T4_gallery.RData")

# npreg vector field on lens, T4b, use global np with fixed bw  ------------------------

LL <- 2

RF_lens <- lens_type[[LL]][, c('x0','y0','z0','xd','yd','zd')] %>% as.matrix()
colnames(RF_lens) <- NULL
RF_com <- lens_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
colnames(RF_com) <- NULL

## ## 
# global -- bwtype= 'fixed', regtype= 'll'
np_eval <- ucl_rot_sm
vf_pred <- matrix(nrow = dim(np_eval)[1], ncol = 3)
for (k in 1:3) {
  # npdata <- data.frame(mc = RF_lens[,1:3], ec = RF_lens[,3+k]) #TODO, use com instead
  npdata <- data.frame(mc = RF_com, ec = RF_lens[,3+k]) #use com
  bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'fixed', regtype= 'll')
  # bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
  model_np <- npreg(bw)
  # model_np <- npreg(ec ~ mc.1 + mc.2 + mc.3, data = npdata)
  for (j in 1:nrow(np_eval)) {
    np_eval_one <- data.frame(mc.1 = np_eval[j,1], mc.2 = np_eval[j,2], mc.3 = np_eval[j,3])
    vf_pred[j,k] <- predict(model_np, newdata = np_eval_one)
  }
}
RF_lens_T4b_pred <- 2*(vf_pred - ucl_rot_sm) + ucl_rot_sm # restore length


## ##
# local -- bwtype= 'adaptive_nn', regtype= 'll'
# reason not to do this is that T4s are quite non-uniform distributed. 
# unclear how many to include; region systematic bias ?? 

# smoothing

# ## use delsgs, neighbor >=5 & com-dist < mean(el)/4, cf "eyemap_lens_med_lop.R"
# pts <- ucl_rot_sm
# vv <- RF_lens_T4b_pred - ucl_rot_sm
# RF_lens_T4b_pred_sm <- vv
# Nnb <- 5
# 
# el <- pts[delsgs_mod$ind1,] - pts[delsgs_mod$ind2,]
# el <- sqrt(rowSums(el^2))
# dcom <- mean(el)/4
# 
# for (j in 1:nrow(pts)) {
#   ind_nb <- c(delsgs_mod[delsgs_mod$ind1 %in% j, "ind2"],
#               delsgs_mod[delsgs_mod$ind2 %in% j, "ind1"])
#   if (length(ind_nb) >= 5) {
#     dd <- (pts[j,] - colMeans(pts[ind_nb,]))^2 %>% sum() %>% sqrt()
#     if (dd < dcom) {
#       RF_lens_T4b_pred_sm[j,] <- vv[j,]*0.5 + colMeans(vv[ind_nb,])*0.5
#     }
#   }
# }
# RF_lens_T4b_pred_sm <- RF_lens_T4b_pred_sm + ucl_rot_sm
# RF_lens_T4b_pred_sm <- sweep(RF_lens_T4b_pred_sm,1,sqrt(rowSums(RF_lens_T4b_pred_sm^2)),'/') #normalize


## ## use nb_ind to smooth
vv <- RF_lens_T4b_pred - ucl_rot_sm
vvnew <- vv
for (m in 1:nrow(nb_ind)) {
  if (!any(is.na(nb_ind[m,]))) {
    vvnew[nb_ind[m,1],] <- vv[nb_ind[m,1],]/2 + colMeans(vv[nb_ind[m,-1],])/2
  }
}
RF_lens_T4b_pred_sm <- ucl_rot_sm + vvnew
RF_lens_T4b_pred_sm <- sweep(RF_lens_T4b_pred_sm,1,sqrt(rowSums(RF_lens_T4b_pred_sm^2)),'/') #normalize


## ## no smoothing
# RF_lens_T4b_pred_sm <- RF_lens_T4b_pred

# PLOT
vf_pred <- RF_lens_T4b_pred_sm
# vf_pred <- RF_lens_T4b_pred
nopen3d()
points3d(ucl_rot_aux[1:Npt,], col = "grey", size = 5)
for (j in 1:nrow(RF_lens)) {
  arrow3d(RF_lens[j,1:3], RF_lens[j,4:6], theta = pi / 12, n = 4, col = "red", type = "rotation")
  # text3d(RF_lens[j,1:3], texts = paste(j), adj = -1)
}
for (j in 1:Npt) {
  arrow3d(ucl_rot_aux[j,], vf_pred[j,1:3], theta = pi / 12, n = 4, col = "blue", type = "rotation")
}


# npreg vector field on lens, T4d ------------------------

LL <- 4

RF_lens <- lens_type[[LL]][, c('x0','y0','z0','xd','yd','zd')] %>% as.matrix()
colnames(RF_lens) <- NULL
RF_com <- lens_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
colnames(RF_com) <- NULL


# - np
np_eval <- ucl_rot_sm
vf_pred <- matrix(nrow = dim(np_eval)[1], ncol = 3)
# for (j in 1:dim(np_eval)[1]) {
#   np_eval_one <- data.frame(mc.1 = np_eval[j,1], mc.2 = np_eval[j,2], mc.3 = np_eval[j,3])
#   for (k in 1:3) {
#     npdata <- data.frame(mc = RF_lens[,1:3], ec = RF_lens[,3+k])
#     model_np <- npreg(ec ~ mc.1 + mc.2 + mc.3, data = npdata)
#     vf_pred[j,k] <- predict(model_np, newdata = np_eval_one)
#   }
# }
for (k in 1:3) {
  npdata <- data.frame(mc = RF_com, ec = RF_lens[,3+k] ) #use com
  bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'fixed', regtype= 'll')
  # bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
  model_np <- npreg(bw)
  for (j in 1:nrow(np_eval)) {
    np_eval_one <- data.frame(mc.1 = np_eval[j,1], mc.2 = np_eval[j,2], mc.3 = np_eval[j,3])
    vf_pred[j,k] <- predict(model_np, newdata = np_eval_one)
  }
}
RF_lens_T4d_pred <- 2*(vf_pred - ucl_rot_sm) + ucl_rot_sm # restore length (approx)


# - smooth, use nb_ind
vv <- RF_lens_T4d_pred - ucl_rot_sm
vvnew <- vv
for (m in 1:nrow(nb_ind)) {
  if (!any(is.na(nb_ind[m,]))) {
    vvnew[nb_ind[m,1],] <- vv[nb_ind[m,1],]/2 + colMeans(vv[nb_ind[m,-1],])/2
  }
}
RF_lens_T4d_pred_sm <- ucl_rot_sm + vvnew
RF_lens_T4d_pred_sm <- sweep(RF_lens_T4d_pred_sm,1,sqrt(rowSums(RF_lens_T4d_pred_sm^2)),'/') #normalize


# PLOT
vf_pred <- RF_lens_T4d_pred
# vf_pred <- RF_lens_T4d_pred_sm
nopen3d()
points3d(ucl_rot_aux[1:Npt,], col = "grey", size = 5)
for (j in 1:nrow(RF_lens)) {
  arrow3d(RF_lens[j,1:3], RF_lens[j,4:6], theta = pi / 12, n = 4, col = "red", type = "rotation")
}
for (j in 1:Npt) {
  # for (j in Npt:nrow(vf_pred)) {
  arrow3d(ucl_rot_aux[j,], vf_pred[j,1:3], theta = pi / 12, n = 4, col = "blue", type = "rotation")
}




# all T4 direction assuming antiparallelism  ------------------------------------------------------

# RF_lens_T4_pred <- cbind(ucl_rot_aux*2 - RF_lens_T4b_pred,
#                          RF_lens_T4b_pred,
#                          ucl_rot_aux*2 - RF_lens_T4d_pred,
#                          RF_lens_T4d_pred)
# # reverse "re-define"
# RF_lens_T4_pred <- RF_lens_T4_pred[1:Npt,]
# # ucl_rot <- ucl_rot_aux[1:Npt,]
# # med_xyz <- med_xyz[1:Npt,]

# This is RF or preferred direction (PD), opposite of predicated orientation of T4 dendrites
RF_lens_T4_pred <- cbind(RF_lens_T4b_pred,
                         ucl_rot_sm*2 - RF_lens_T4b_pred,
                         RF_lens_T4d_pred,
                         ucl_rot_sm*2 - RF_lens_T4d_pred )

RF_lens_T4_pred_sm <- cbind(RF_lens_T4b_pred_sm,
                         ucl_rot_sm*2 - RF_lens_T4b_pred_sm,
                         RF_lens_T4d_pred_sm,
                         ucl_rot_sm*2 - RF_lens_T4d_pred_sm )


# m <- 2
# nopen3d()
# points3d(ucl_rot_sm, col = "grey", size = 5)
# for (j in 1:dim(RF_lens_T4_pred)[1]) {
#   arrow3d(ucl_rot_sm[j,], RF_lens_T4_pred[j,(3*(m-1)+1):(3*(m-1)+3)], theta = pi / 9, n = 4, col = "gray", type = "rotation")
#   arrow3d(ucl_rot_sm[j,], RF_lens_T4_pred_sm[j,(3*(m-1)+1):(3*(m-1)+3)], theta = pi / 9, n = 8, col = "blue", type = "rotation")
# }

# SAVE, RF_lens_, RF_med, ucl_rot
# save(ucl_rot,ucl_rot_sm,med_xyz, RF_lens_T4_pred, RF_lens_T4b, RF_med_T4b, RF_lens_T4d, RF_med_T4d,
#      file = "data/RF_pred.RData")
# save(ucl_rot,ucl_rot_sm,med_xyz, RF_lens_T4_pred, RF_lens_T4_pred_sm, RF_lens_T4b, RF_med_T4b, RF_lens_T4d, RF_med_T4d,
#      file = "data/T4_RF_pred.RData")
# save(ucl_rot,ucl_rot_sm,med_xyz, RF_lens_T4_pred, RF_lens_T4_pred_sm,
#      file = "data/T4_RF_pred.RData") #2021

save(RF_lens_T4_pred, RF_lens_T4_pred_sm, file = "data/T4_RF_pred.RData") #2023


