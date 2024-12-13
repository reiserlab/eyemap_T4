# cf. /ReiserGroup/p_eyemap/microCT

library(natverse)
library(tidyverse)


library(clue) # hungarian
library(alphashape3d) # ashape3d
library(alphahull)
library(RColorBrewer) #palette
library(reshape2) #melt
library(deldir)


# clean everythign up.
rm(list=ls())
 
# # cross product 
# cross3D <- function(a, b){
#   prod <- matrix(c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1]), ncol = 1)
#   return(prod)
# }

source("eyemap_func.R")



# - S2 grid for kde
S2grid <- matrix(c(0,0), ncol = 2) # [theta phi] in deg
darc <- pi/180*4
thetaN <- 180/4 - 1
for (ti in 1:thetaN) {
  theta <- ti * darc / 1
  ddeg <- darc / sin(theta) / pi *180
  phiN <- floor(180 / ddeg)
  tpmat <- cbind(rep(theta/pi*180, 2*phiN+1), seq(-phiN, phiN, by = 1)*ddeg )
  S2grid <- rbind(S2grid, tpmat)
}
S2grid <- rbind(S2grid, c(180,0)) #south pole
darea <- 4*pi/nrow(S2grid)
S2grid <- cbind(S2grid, Mollweide(S2grid))
colnames(S2grid) <- c('t','p','xM','yM')

# Mollweide guidelines ---------------------------------------------------------------

# load("data/eyemap.RData")
# 
# ucl_rot_Mo <- ucl_rot_sm
# colnames(ucl_rot_Mo) <- c('x','y','z')
# ucl_rot_Mo %<>% as_tibble() %>%  
#   mutate(y = -y) %>%
#   mutate(theta = acos(z)) %>%
#   mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
#   mutate(t = theta / pi * 180, p = phi/pi*180) %>%  
#   mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
#   as.data.frame()
# ucl_rot_Mo <- Mollweide(ucl_rot_Mo[,c('t', 'p')])
# colnames(ucl_rot_Mo) <- c('xM','yM')
# rownames(ucl_rot_Mo) <- rownames(ucl_rot_sm)

# Mollweide guidelines
Mollweide_ori <- c(0,0)
Mollweide_mul <- 1

bkgd_eq <- Mollweide(cbind(90, seq(-180, 180, by = 2)))
bkgd_eq <- sweep(bkgd_eq*Mollweide_mul, 2, Mollweide_ori, '+')
colnames(bkgd_eq) <- c('xM','yM')
bkgd_eq_p45 <- Mollweide(cbind(45, seq(-180, 180, by = 2)))
bkgd_eq_p45 <- sweep(bkgd_eq_p45*Mollweide_mul, 2, Mollweide_ori, '+')
colnames(bkgd_eq_p45) <- c('xM','yM')
bkgd_eq_m45 <- Mollweide(cbind(135, seq(-180, 180, by = 2)))
bkgd_eq_m45 <- sweep(bkgd_eq_m45*Mollweide_mul, 2, Mollweide_ori, '+')
colnames(bkgd_eq_m45) <- c('xM','yM')
bkgd_mer_ww <- Mollweide(cbind(seq(0, 180, by = 10), rep(-180,19)))
bkgd_mer_w <- Mollweide(cbind(seq(180, 0, by = -10), rep(-90,19)))
bkgd_mer <- Mollweide(cbind(seq(0, 180, by = 10), rep(0,19)))
bkgd_mer_e <- Mollweide(cbind(seq(180, 0, by = -10), rep(90,19)))
bkgd_mer_ee <- Mollweide(cbind(seq(0, 180, by = 10), rep(180,19)))
bkgd_mer <- rbind(bkgd_mer_ww, bkgd_mer_w, bkgd_mer,bkgd_mer_e,bkgd_mer_ee)
bkgd_mer <- sweep(bkgd_mer*Mollweide_mul, 2, Mollweide_ori, '+')
colnames(bkgd_mer) <- c('xM','yM')

bkgd_chull <- rbind(bkgd_mer_ww, bkgd_mer_ee[seq(nrow(bkgd_mer_ee),1,by=-1),])
colnames(bkgd_chull) <- c('xM','yM')

# proc 2023 test eyemap -----------------------------------------------------------

# fn <- c('20230705_70EtOH_')
# 
# for (f in 1:length(fn)) {
#   lens_left <- read.csv(paste("../microCT/2023_test/", fn[f], "Left_Lens.csv", sep=''), header=T)
#   lens_right <- read.csv(paste("../microCT/2023_test/", fn[f], "Right_Lens.csv", sep=''), header=T)
#   cone_left <- read.csv(paste("../microCT/2023_test/", fn[f], "Left_PR.csv", sep=''), header=T)
#   cone_right <- read.csv(paste("../microCT/2023_test/", fn[f], "Right_PR.csv", sep=''), header=T)
#   
#   lens <- rbind(lens_left, lens_right)
#   cone <- rbind(cone_left, cone_right)
# 
#   # extract
#   lens <- lens[, c('Position.X','Position.Y','Position.Z')]
#   colnames(lens) <- c("x", "y", "z")
#   # lens <- lens[order(lens$x, lens$y,lens$z),]
#   # rownames(lens) <- seq(1,length(lens$x), by = 1)
#   
#   cone <- cone[, c('Position.X','Position.Y','Position.Z','Default.Labels')]
#   colnames(cone) <- c("x", "y", "z", "label")
#   # cone <- cone[order(cone$x, cone$y,cone$z),]
#   # rownames(cone) <- seq(1,length(cone$x), by = 1)
#   ind_Down <- which(cone$label == "DOWN")
#   ind_Up <- which(cone$label == "UP")
#   cone <- cone[,c("x", "y", "z")]
#   
#   # - alignment
#   # equator
#   v_up <- colMeans(cone[ind_Up,]) - colMeans(cone[ind_Down,])
#   cone_equator <- cone[c(ind_Up, ind_Down),]
#   pt_theta90 <- cone_equator
#   # lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z)
#   # vn <- c(-1, lm_plane$coefficients[2], lm_plane$coefficients[3]) #norm
#   vn <- prcomp(pt_theta90)$rotation[,3]
#   vn <- vn/sqrt(sum(vn^2))
#   # which side is up
#   if (sum(vn * v_up) < 0) {
#     vn <- -vn
#   }
#   # pt_center <- colMeans(cone)
#   
#   # front via pca, use both eyes
#   eye_pca <- prcomp(rbind(lens, cone))
#   vf <- eye_pca$rotation[,3]
#   vf <- vf - c(vf %*% vn) * vn
#   vf <- vf / sqrt(sum(vf^2))
#   if (vf %*% (colMeans(lens) - eye_pca$center) < 0) {
#     vf <- - vf
#   }
#   
#   vleft <- cross3D(vn, vf)
#   
#   # # DEBUG
#   # nopen3d()
#   # points3d(cone)
#   # points3d(cone[ind_Up,], col='red', size =10)
#   # points3d(cone[ind_Down,], col='blue', size =10)
#   # # arrow3d(eye_pca$center, eye_pca$center + eye_pca$rotation[,2] * 100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   # arrow3d(eye_pca$center, eye_pca$center + vn * 100, theta = pi/6,n = 4, col="red", type = "rotation")
#   # arrow3d(eye_pca$center, eye_pca$center + vf * 100, theta = pi/6,n = 4, col="magenta", type = "rotation")
# 
#   
#   # - remove shifted equator pts
#   # centering
#   mrot <- cbind(vf,vleft,vn)
#   lens_tmp <- sweep(as.matrix(lens),2,eye_pca$center,'-') %*% mrot
#   cone_tmp <- sweep(as.matrix(cone),2,eye_pca$center,'-') %*% mrot
#   ind_Up2 <- ind_Up[cone_tmp[ind_Up,1] > 0]
#   ind_Down2 <- ind_Down[cone_tmp[ind_Down,1] > 0]
#   
#   # alignment again
#   v_up <- colMeans(cone[ind_Up2,]) - colMeans(cone[ind_Down2,])
#   cone_equator <- cone[c(ind_Up2, ind_Down2),]
#   pt_theta90 <- cone_equator
#   # lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z)
#   # vn <- c(-1, lm_plane$coefficients[2], lm_plane$coefficients[3]) #norm
#   vn <- prcomp(pt_theta90)$rotation[,3]
#   vn <- vn/sqrt(sum(vn^2))
#   # which side is up
#   if (sum(vn * v_up) < 0) {
#     vn <- -vn
#   }
#   # front via pca, use both eyes
#   eye_pca <- prcomp(rbind(lens, cone))
#   vf <- eye_pca$rotation[,3]
#   vf <- vf - c(vf %*% vn) * vn
#   vf <- vf / sqrt(sum(vf^2))
#   if (vf %*% (colMeans(lens) - eye_pca$center) < 0) {
#     vf <- - vf
#   }
#   # separate left and right
#   vleft <- cross3D(vn, vf)
#   ind_left_lens <- c(t(vleft) %*% t(sweep(lens,2,eye_pca$center,'-')) > 0)
#   ind_left_cone <- c(t(vleft) %*% t(sweep(cone,2,eye_pca$center,'-')) > 0)
#   
#   
#   # - matching, number of points SHOULD be equal
#   L <- dim(cone)[1]
#   mDist <- matrix(ncol = L, nrow = L)
#   dr <- cone
#   dc <- lens
#   for (j in 1:L) {
#     mDist[j,] <- as.matrix(dist(rbind(dr[j,], dc)))[-1,1]
#   }
#   i_match <- solve_LSAP(mDist)
#   
#   # PLOT
#   nopen3d()
#   points3d(lens[ind_left_lens,], size = 5, color = 'magenta')
#   points3d(lens[!ind_left_lens,], size = 5, color = 'green')
#   points3d(cone, size = 5, color = 'blue')
#   points3d(cone[c(ind_Down2, ind_Up2),], size=10, col='gold2')
#   points3d(lens[i_match[c(ind_Down2, ind_Up2)],], size=10, col='cyan')
#   arrow3d(eye_pca$center, eye_pca$center + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   # arrow3d(c(0,0,0), c(0,0,0) + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   planes3d(vn[1],vn[2],vn[3], -sum(eye_pca$center*vn), alpha = 0.5)
#   # planes3d(0,0,1, 0, alpha = 0.3)
#   arrow3d(eye_pca$center, eye_pca$center + vf*100, theta = pi/6,n = 4, col="brown", type = "rotation")
#   # arrow3d(c(0,0,0), c(0,0,0) + vf*100, theta = pi/6,n = 4, col="brown", type = "rotation")
#   # for (j in 1:L) {
#   #   lines3d(rbind(dr[j,], dc[i_match[j],]), col = "grey")
#   # }
#   linemat <- matrix(t(cbind(cone, lens[i_match,])), ncol = 3, byrow = T)
#   segments3d(linemat, color = "grey")
#   # title3d(fn[f])
#   
#   # index of lens
#   ind_Up_lens <- i_match[ind_Up2]
#   ind_Down_lens <- i_match[ind_Down2]
#   
#   # # DEBUG, equator viewing dir on a plane ?
#   # nopen3d()
#   # planes3d(vn[1],vn[2],vn[3], -sum(eye_pca$center*vn), alpha = 0.5)
#   # points3d(cone[c(ind_Up2, ind_Down2),], size=6)
#   # # points3d(lens[c(ind_Up_lens, ind_Down_lens),], size=6, col='blue')
#   # points3d(cone[c(ind_Up2, ind_Down2),] + 10*(lens[c(ind_Up_lens, ind_Down_lens),] - cone[c(ind_Up2, ind_Down2),]), size=8, col='red')
#   # arrow3d(eye_pca$center, eye_pca$center + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   # title3d(fn[f])
#   
#   # - centering
#   mrot <- cbind(vf,vleft,vn)
#   lens <- sweep(as.matrix(lens),2,eye_pca$center,'-') %*% mrot
#   cone <- sweep(as.matrix(cone),2,eye_pca$center,'-') %*% mrot
#   colnames(lens) <- c('x','y','z')
#   colnames(cone) <- c('x','y','z')
#   # vf <- vf %*% mrot
#   # vleft <- c(vleft) %*% mrot
#   # vn <- vn %*% mrot
#   vf <- c(1,0,0)
#   vleft <- c(0,1,0)
#   vn <- c(0,0,1)
#   
#   # - direction unit vec
#   ucl <- lens[i_match,] - cone #use cone index
#   ucl_norm <- sqrt(rowSums(ucl^2))
#   ucl <- ucl / ucl_norm
#   
#   # rotation to canonical coordinate, this dose NOT make sense
#   
#   ### ### vn -> +z, vf -> +x, this way, equatorial ucl not on equatorial plane given lens positions
#   mrot <- cbind(vf,vleft,vn)
#   ucl_rot <- as.matrix(ucl) %*% mrot
#   
#   
#   ### ### based on equator and ucl
#   # ucl_pca <- prcomp(ucl[c(ind_Up2, ind_Down2), ]) # can NOT do shift
#   # vn2 <- ucl_pca$rotation[,3]
#   
#   pt_theta90 <- as.data.frame(ucl[c(ind_Up2, ind_Down2), ])
#   lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z + 0)
#   vn2 <- c(-1, lm_plane$coefficients[1], lm_plane$coefficients[2]) #norm
#   vn2 <- vn2/sqrt(sum(vn2^2))
#   # which side is up
#   if (vn2 %*% c(0,0,1) < 0) {
#     vn2 <- -vn2
#   }
#   
#   vleft2 <- c(0,1,0) - c(vn2 %*% c(0,1,0)) * vn2
#   vleft2 <- vleft2/sqrt(sum(vleft2^2))
#   vf2 <- cross3D(vleft2, vn2)
#   
#   mrot <- cbind(vf2,vleft2,vn2)
#   # ucl_rot2 <- sweep(as.matrix(ucl),2, ucl_pca$center, '-') %*% mrot 
#   ucl_rot2 <- as.matrix(ucl) %*% mrot
#   colnames(ucl_rot2) <- c('x','y','z')
#   
#   # ucl_rot2[,2] <- -ucl_rot2[,2] # right eye, phi angle viewed from outside,
#   
#   # # PLOT, unit sphere, rotated
#   # nopen3d()
#   # points3d(ucl_rot2[ind_left_cone, ], size = 5, color = 'magenta')
#   # points3d(ucl_rot2[!ind_left_cone,], size = 5, color = 'green')
#   # points3d(ucl_rot2[c(ind_Up2, ind_Down2),], size = 10, color = 'gray20')
#   # spheres3d(0,0,0,0.99, col='grey', alpha=0.3, lit=F)
#   # arrow3d(c(0,0,0), vn2*1, theta = pi/6,n = 4, col="blue", type = "rotation")
#   # planes3d(0,0,1, 0, alpha = 0.3)
#   # planes3d(0,1,0, 0, alpha = 0.3)
#   # axes3d(c('x','y','z')); title3d('','','x','y','z')
#   # title3d(fn[f])
#   
#   # save(lens, cone, ind_Up, ind_Down, ind_left_lens, ind_left_cone,
#   #      eye_pca, vn, vf, vleft, L, i_match, ucl,
#   #      vn2, vf2, vleft2, ucl_rot2,
#   #      file = paste("../microCT/data/", fn[f], ".RData", sep='') )
#   
#   save(lens, cone, ind_Up, ind_Down, ind_left_lens, ind_left_cone,
#        eye_pca, vn, vf, vleft, ucl_rot, L, i_match, ucl,
#        ind_Up2, ind_Down2, vn2, vf2, vleft2, ucl_rot2,
#        file = paste("../microCT/data/eq_corr/", fn[f], ".RData", sep='') )
# }


# proc 2023 test (incl. LM) eyemap, one eye -----------------------------------------------------------

# fn <- c('20230727_EtOH_1PTA_')
# # fn <- c('LM_1_')
# 
# for (f in 1:length(fn)) {
#   lens_left <- read.csv(paste("../microCT/2023_test/", fn[f], "Left_Lens.csv", sep=''), header=T)
#   lens_right <- read.csv(paste("../microCT/2023_test/", fn[f], "Right_Lens.csv", sep=''), header=T)
#   cone_right <- read.csv(paste("../microCT/2023_test/", fn[f], "Right_PR.csv", sep=''), header=T)
#   
#   # one eye
#   lens <- lens_right
#   cone <- cone_right
# 
#   # extract
#   lens <- lens[, c('Position.X','Position.Y','Position.Z')]
#   colnames(lens) <- c("x", "y", "z")
#   
#   # use left lens with UP/DOWN with right cone for coord system
#   lens_left <- lens_left[, c('Position.X','Position.Y','Position.Z','Default.Labels')]
#   colnames(lens_left) <- c("x", "y", "z", "label")
#   ind_Down_ll <- which(lens_left$label == "DOWN")
#   ind_Up_ll <- which(lens_left$label == "UP")
#   lens_left <- lens_left[,1:3]
#   
#   # lens <- lens[order(lens$x, lens$y,lens$z),]
#   # rownames(lens) <- seq(1,length(lens$x), by = 1)
#   
#   cone <- cone[, c('Position.X','Position.Y','Position.Z','Default.Labels')]
#   colnames(cone) <- c("x", "y", "z", "label")
#   # cone <- cone[order(cone$x, cone$y,cone$z),]
#   # rownames(cone) <- seq(1,length(cone$x), by = 1)
#   ind_Down <- which(cone$label == "DOWN")
#   ind_Up <- which(cone$label == "UP")
#   cone <- cone[,c("x", "y", "z")]
#   
#   # - alignment
#   # equator
#   v_up <- colMeans(rbind(cone[ind_Up,],lens_left[ind_Up_ll,])) - 
#     colMeans(rbind(cone[ind_Down,], lens_left[ind_Down_ll,]))
#   cone_equator <- rbind(cone[c(ind_Up, ind_Down),], lens_left[c(ind_Up_ll,ind_Down_ll),])
#   pt_theta90 <- cone_equator
#   # lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z)
#   # vn <- c(-1, lm_plane$coefficients[2], lm_plane$coefficients[3]) #norm
#   vn <- prcomp(pt_theta90)$rotation[,3]
#   vn <- vn/sqrt(sum(vn^2))
#   # which side is up
#   if (sum(vn * v_up) < 0) {
#     vn <- -vn
#   }
#   # pt_center <- colMeans(cone)
#   
#   # front via pca, use both eyes
#   eye_pca <- prcomp(rbind(lens, lens_left))
#   vf <- eye_pca$rotation[,2]
#   vf <- vf - c(vf %*% vn) * vn
#   vf <- vf / sqrt(sum(vf^2))
#   # if (vf %*% (colMeans(lens) - eye_pca$center) < 0) {
#   if (vf %*% (colMeans(lens) - colMeans(cone)) < 0) {
#     vf <- - vf
#   }
#   
#   vleft <- cross3D(vn, vf)
#   
#   # DEBUG
#   nopen3d()
#   points3d(cone)
#   points3d(lens_left)
#   points3d(lens_left[ind_Up_ll,], col='red', size =10)
#   points3d(lens_left[ind_Down_ll,], col='blue', size =10)
#   points3d(cone[ind_Up,], col='red', size =10)
#   points3d(cone[ind_Down,], col='blue', size =10)
#   # arrow3d(eye_pca$center, eye_pca$center + eye_pca$rotation[,2] * 100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   arrow3d(eye_pca$center, eye_pca$center + vn * 100, theta = pi/6,n = 4, col="red", type = "rotation")
#   arrow3d(eye_pca$center, eye_pca$center + vf * 100, theta = pi/6,n = 4, col="magenta", type = "rotation")
# 
#   
#   # # separate left and right
#   # vleft <- cross3D(vn, vf)
#   # ind_left_lens <- c(t(vleft) %*% t(sweep(lens,2,eye_pca$center,'-')) > 0)
#   # ind_left_cone <- c(t(vleft) %*% t(sweep(cone,2,eye_pca$center,'-')) > 0)
#   
#   
#   # - matching, number of points SHOULD be equal
#   L <- dim(cone)[1]
#   mDist <- matrix(ncol = L, nrow = L)
#   dr <- cone
#   dc <- lens
#   for (j in 1:L) {
#     mDist[j,] <- as.matrix(dist(rbind(dr[j,], dc)))[-1,1]
#   }
#   i_match <- solve_LSAP(mDist)
#   
#   # PLOT
#   nopen3d()
#   points3d(lens[i_match[c(ind_Down)],], size=5, col='magenta')
#   points3d(lens[i_match[c(ind_Up)],], size=5, col='cyan')
#   points3d(lens, size = 5, color = 'gray')
#   points3d(cone, size = 5, color = 'brown')
#   arrow3d(eye_pca$center, eye_pca$center + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   # arrow3d(c(0,0,0), c(0,0,0) + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   planes3d(vn[1],vn[2],vn[3], -sum(eye_pca$center*vn), alpha = 0.5)
#   # planes3d(0,0,1, 0, alpha = 0.3)
#   arrow3d(eye_pca$center, eye_pca$center + vf*100, theta = pi/6,n = 4, col="brown", type = "rotation")
#   # arrow3d(c(0,0,0), c(0,0,0) + vf*100, theta = pi/6,n = 4, col="brown", type = "rotation")
#   # for (j in 1:L) {
#   #   lines3d(rbind(dr[j,], dc[i_match[j],]), col = "grey")
#   # }
#   linemat <- matrix(t(cbind(cone, lens[i_match,])), ncol = 3, byrow = T)
#   segments3d(linemat, color = "grey")
#   # title3d(fn[f])
#   
#   # # index of lens
#   # ind_Up_lens <- i_match[ind_Up2]
#   # ind_Down_lens <- i_match[ind_Down2]
#   
#   # # DEBUG, equator viewing dir on a plane ?
#   # nopen3d()
#   # planes3d(vn[1],vn[2],vn[3], -sum(eye_pca$center*vn), alpha = 0.5)
#   # points3d(cone[c(ind_Up2, ind_Down2),], size=6)
#   # # points3d(lens[c(ind_Up_lens, ind_Down_lens),], size=6, col='blue')
#   # points3d(cone[c(ind_Up2, ind_Down2),] + 10*(lens[c(ind_Up_lens, ind_Down_lens),] - cone[c(ind_Up2, ind_Down2),]), size=8, col='red')
#   # arrow3d(eye_pca$center, eye_pca$center + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   # title3d(fn[f])
#   
#   # - centering
#   mrot <- cbind(vf,vleft,vn)
#   lens <- sweep(as.matrix(lens),2,eye_pca$center,'-') %*% mrot
#   cone <- sweep(as.matrix(cone),2,eye_pca$center,'-') %*% mrot
#   colnames(lens) <- c('x','y','z')
#   colnames(cone) <- c('x','y','z')
#   # vf <- vf %*% mrot
#   # vleft <- c(vleft) %*% mrot
#   # vn <- vn %*% mrot
#   vf <- c(1,0,0)
#   vleft <- c(0,1,0)
#   vn <- c(0,0,1)
#   
#   # # smoothing 1
#   # lens_orig <- lens
#   # cone_orig <- cone
#   # 
#   # ls <- reghex(cone)
#   # nb_ind <- ls[[1]]
#   # 
#   # cone_new <- matrix(nrow = nrow(cone), ncol = ncol(cone))
#   # lens_new <- matrix(nrow = nrow(cone), ncol = ncol(cone))
#   # for (m in 1:nrow(nb_ind)) {
#   #   cone_new[nb_ind[m,1],] <- cone[nb_ind[m,1],]/2 + colMeans(cone[nb_ind[m,-1],])/2
#   #   lens_new[i_match[nb_ind[m,1]],] <- lens[i_match[nb_ind[m,1]],]/2 + colMeans(lens[i_match[nb_ind[m,-1]],])/2
#   # }
#   # cone <- cone_new
#   # lens <- lens_new
#   
#   
#   # - direction unit vec
#   ucl <- lens[i_match,] - cone #use cone index
#   ucl_norm <- sqrt(rowSums(ucl^2))
#   ucl <- ucl / ucl_norm
#   
#   # rotation to canonical coordinate, 
#   # vn -> +z, vf -> +x, this way, equatorial ucl not on equatorial plane given lens positions
#   mrot <- cbind(vf,vleft,vn)
#   ucl_rot <- as.matrix(ucl) %*% mrot
#   
#   
#   # smoothing 2
#   ucl_rot_orig <- ucl_rot
# 
#   ls <- reghex(cone)
#   nb_ind <- ls[[1]]
#   
#   ucl_rot_new <- matrix(nrow = nrow(cone), ncol = ncol(cone))
#   for (m in 1:nrow(nb_ind)) {
#     ucl_rot_new[nb_ind[m,1],] <- ucl_rot[nb_ind[m,1],]/2 + colMeans(ucl_rot[nb_ind[m,-1],])/2
#   }
#   ucl_rot <- ucl_rot_new
#   
# 
#   # sphere
#   nopen3d()
#   points3d(ucl_rot)
#   points3d(ucl_rot[ind_Down,], col='magenta', size=10)
#   points3d(ucl_rot[ind_Up,], col='cyan', size=10)
#   spheres3d(0,0,0, 0.99, col='grey90', alpha=1, lit=F)
#   planes3d(0,1,0, 0, alpha = 0.2)
#   
#   # plot Mollweide
#   ucl_rot_Mo <- ucl_rot
#   colnames(ucl_rot_Mo) <- c('x','y','z')
#   ucl_rot_Mo %<>% as_tibble() %>%  
#     mutate(y = -y) %>%
#     mutate(theta = acos(z)) %>%
#     mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
#     mutate(t = theta / pi * 180, p = phi/pi*180) %>%  
#     mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
#     as.data.frame()
#   ucl_rot_Mo <- Mollweide(ucl_rot_Mo[,c('t', 'p')])
#   colnames(ucl_rot_Mo) <- c('xM','yM')
#   rownames(ucl_rot_Mo) <- rownames(ucl_rot)
#   
#   # need run Mollweide setup
#   
#   df <- as.data.frame(ucl_rot_Mo)
#   windows(width = 16, height = 8)
#   # pdf("viewing_angles_bi.pdf")
#   plt <- plt_Mo + 
#     geom_point(data=df, aes(x=xM, y=yM), size=1)
#   
#   plt
#   # save(lens, cone, ind_Up, ind_Down, ind_left_lens, ind_left_cone,
#   #      eye_pca, vn, vf, vleft, ucl_rot, L, i_match, ucl,
#   #      file = paste("../microCT/data/eq_corr/", fn[f], ".RData", sep='') )
# }


# proc 2023 eyemap -----------------------------------------------------------

# load position data, match lens to cone/tip to compute directions
# smooth positions and directions (almost redundant)

# fn <- c(
#   # '20230926',
#   "20231107",
#   # # '20230925_2',
#   # # "20240125",
#   # "20240206",
#   # "20240510",
#   # "20240513",
#   # "20240520",
#   # "20240522",
#   # "20240524",
#   "20240530",
#   # "20240530_dried",
#   # "20240612",
#   "20240701"
#   )

fn <- c(
  '20240701',
  '20231107',
  '20240530'
)

for (f in 1:length(fn)) {
  # lens <- read.csv(paste("../microCT/2023_eyemap/", fn[f], "_position/lens.csv", sep=''), header=T) 
  # cone <- read.csv(paste("../microCT/2023_eyemap/", fn[f], "_position/cone.csv", sep=''), header=T) 
  lens <- read.csv(paste("data/microCT/", fn[f], "_position/lens.csv", sep=''), header=T) 
  cone <- read.csv(paste("data/microCT/", fn[f], "_position/cone.csv", sep=''), header=T) 
  
  # extract
  lens <- lens[, c('Position.X','Position.Y','Position.Z')]
  colnames(lens) <- c("x", "y", "z")
  lens <- lens[order(lens$x, lens$y,lens$z),]
  rownames(lens) <- seq(1,length(lens$x), by = 1)
  
  cone <- cone[, c('Position.X','Position.Y','Position.Z','Default.Labels')]
  colnames(cone) <- c("x", "y", "z", "label")
  cone <- cone[order(cone$x, cone$y,cone$z),]
  rownames(cone) <- seq(1,length(cone$x), by = 1)
  ind_Down <- which(cone$label == "DOWN")
  ind_Up <- which(cone$label == "UP")
  cone <- cone[,c("x", "y", "z")]
  
  # - alignment
  # equator
  v_up <- colMeans(cone[ind_Up,]) - colMeans(cone[ind_Down,])
  cone_equator <- cone[c(ind_Up, ind_Down),]
  pt_theta90 <- cone_equator
  # lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z)
  # vn <- c(-1, lm_plane$coefficients[2], lm_plane$coefficients[3]) #norm
  vn <- prcomp(pt_theta90)$rotation[,3]
  vn <- vn/sqrt(sum(vn^2))
  # which side is up
  if (sum(vn * v_up) < 0) {
    vn <- -vn
  }
  # pt_center <- colMeans(cone)
  
  # front via pca, use both eyes
  eye_pca <- prcomp(rbind(lens, cone))
  vf <- eye_pca$rotation[,3]
  vf <- vf - c(vf %*% vn) * vn
  vf <- vf / sqrt(sum(vf^2))
  if (vf %*% (colMeans(lens) - eye_pca$center) < 0) {
    vf <- - vf
  }
  vleft <- cross3D(vn, vf)

  # - remove shifted equator pts
  # if (FALSE) {
  #   # centering
  #   mrot <- cbind(vf,vleft,vn)
  #   lens_tmp <- sweep(as.matrix(lens),2,eye_pca$center,'-') %*% mrot
  #   cone_tmp <- sweep(as.matrix(cone),2,eye_pca$center,'-') %*% mrot
  #   ind_Up2 <- ind_Up[cone_tmp[ind_Up,1] > 0]
  #   ind_Down2 <- ind_Down[cone_tmp[ind_Down,1] > 0]
  #   
  #   # alignment again
  #   v_up <- colMeans(cone[ind_Up2,]) - colMeans(cone[ind_Down2,])
  #   cone_equator <- cone[c(ind_Up2, ind_Down2),]
  #   pt_theta90 <- cone_equator
  #   # lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z)
  #   # vn <- c(-1, lm_plane$coefficients[2], lm_plane$coefficients[3]) #norm
  #   vn <- prcomp(pt_theta90)$rotation[,3]
  #   vn <- vn/sqrt(sum(vn^2))
  #   # which side is up
  #   if (sum(vn * v_up) < 0) {
  #     vn <- -vn
  #   }
  #   # front via pca, use both eyes
  #   eye_pca <- prcomp(rbind(lens, cone))
  #   vf <- eye_pca$rotation[,3]
  #   vf <- vf - c(vf %*% vn) * vn
  #   vf <- vf / sqrt(sum(vf^2))
  #   if (vf %*% (colMeans(lens) - eye_pca$center) < 0) {
  #     vf <- - vf
  #   }
  # } else {
  #   ind_Up2 <- ind_Up
  #   ind_Down2 <- ind_Down
  # }
  
  # - separate left and right
  vleft <- cross3D(vn, vf)
  ind_left_lens <- c(t(vleft) %*% t(sweep(lens,2,eye_pca$center,'-')) > 0)
  ind_left_cone <- c(t(vleft) %*% t(sweep(cone,2,eye_pca$center,'-')) > 0)
  
  # - matching, number of points SHOULD be equal
  L <- dim(cone)[1]
  Lr <- L
  mDist <- matrix(ncol = L, nrow = L)
  # Lc <- dim(cone)[1]
  # Lr <- dim(lens)[1] # larger
  # mDist <- matrix(ncol = Lc, nrow = Lr)
  dr <- cone
  dc <- lens
  for (j in 1:Lr) {
    mDist[j,] <- as.matrix(dist(rbind(dr[j,], dc)))[-1,1]
  }
  i_match <- solve_LSAP(mDist)
  
  # PLOT
  nopen3d()
  points3d(lens[ind_left_lens,], size = 5, color = 'magenta')
  points3d(lens[!ind_left_lens,], size = 5, color = 'green')
  points3d(cone, size = 5, color = 'blue')
  points3d(cone[c(ind_Down2, ind_Up2),], size=10, col='gold2')
  points3d(lens[i_match[c(ind_Down2, ind_Up2)],], size=10, col='cyan')
  # points3d(lens[match(c(ind_Down2, ind_Up2), i_match),], size=10, col='cyan')
  arrow3d(eye_pca$center, eye_pca$center + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
  # arrow3d(c(0,0,0), c(0,0,0) + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
  planes3d(vn[1],vn[2],vn[3], -sum(eye_pca$center*vn), alpha = 0.5)
  # planes3d(0,0,1, 0, alpha = 0.3)
  arrow3d(eye_pca$center, eye_pca$center + vf*100, theta = pi/6,n = 4, col="brown", type = "rotation")
  # arrow3d(c(0,0,0), c(0,0,0) + vf*100, theta = pi/6,n = 4, col="brown", type = "rotation")
  # for (j in 1:L) {
  #   lines3d(rbind(dr[j,], dc[i_match[j],]), col = "grey")
  # }
  linemat <- matrix(t(cbind(dr, dc[i_match,])), ncol = 3, byrow = T)
  segments3d(linemat, color = "grey")
  # title3d(fn[f])
  
  # index of lens
  ind_Up_lens <- i_match[ind_Up2]
  ind_Down_lens <- i_match[ind_Down2]
  
  # # DEBUG, equator viewing dir on a plane ?
  # nopen3d()
  # planes3d(vn[1],vn[2],vn[3], -sum(eye_pca$center*vn), alpha = 0.5)
  # points3d(cone[c(ind_Up2, ind_Down2),], size=6)
  # # points3d(lens[c(ind_Up_lens, ind_Down_lens),], size=6, col='blue')
  # points3d(cone[c(ind_Up2, ind_Down2),] + 10*(lens[c(ind_Up_lens, ind_Down_lens),] - cone[c(ind_Up2, ind_Down2),]), size=8, col='red')
  # arrow3d(eye_pca$center, eye_pca$center + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
  # title3d(fn[f])
  
  # - centering
  mrot <- cbind(vf,vleft,vn)
  lens <- sweep(as.matrix(lens),2,eye_pca$center,'-') %*% mrot
  cone <- sweep(as.matrix(cone),2,eye_pca$center,'-') %*% mrot
  colnames(lens) <- c('x','y','z')
  colnames(cone) <- c('x','y','z')
  # vf <- vf %*% mrot
  # vleft <- c(vleft) %*% mrot
  # vn <- vn %*% mrot
  vf <- c(1,0,0)
  vleft <- c(0,1,0)
  vn <- c(0,0,1)
  
  # - smoothing  (separately --> not much difference)
  lens_sm <- lens
  cone_sm <- cone

  # left lens
  pts <- lens[ind_left_lens,]
  ls <- reghex(pts)
  nb_ind <- ls[[1]]
  ptsnew <- pts
  for (m in 1:nrow(nb_ind)) {
    if (!any(is.na(nb_ind[m,]))) {
      ptsnew[nb_ind[m,1],] <- pts[nb_ind[m,1],]/2 + colMeans(pts[nb_ind[m,-1],])/2
    }
  }
  lens_sm[ind_left_lens,] <- ptsnew

  # left cone
  pts <- cone[ind_left_cone,]
  ls <- reghex(pts)
  nb_ind <- ls[[1]]
  ptsnew <- pts
  for (m in 1:nrow(nb_ind)) {
    if (!any(is.na(nb_ind[m,]))) {
      ptsnew[nb_ind[m,1],] <- pts[nb_ind[m,1],]/2 + colMeans(pts[nb_ind[m,-1],])/2
    }
  }
  cone_sm[ind_left_cone,] <- ptsnew

  # right lens
  pts <- lens[!ind_left_lens,]
  ls <- reghex(pts)
  nb_ind <- ls[[1]]
  ptsnew <- pts
  for (m in 1:nrow(nb_ind)) {
    if (!any(is.na(nb_ind[m,]))) {
      ptsnew[nb_ind[m,1],] <- pts[nb_ind[m,1],]/2 + colMeans(pts[nb_ind[m,-1],])/2
    }
  }
  lens_sm[!ind_left_lens,] <- ptsnew

  # right cone
  pts <- cone[!ind_left_cone,]
  ls <- reghex(pts)
  nb_ind <- ls[[1]]
  ptsnew <- pts
  for (m in 1:nrow(nb_ind)) {
    if (!any(is.na(nb_ind[m,]))) {
      ptsnew[nb_ind[m,1],] <- pts[nb_ind[m,1],]/2 + colMeans(pts[nb_ind[m,-1],])/2
    }
  }
  cone_sm[!ind_left_cone,] <- ptsnew
  
  lens0 <- lens
  cone0 <- cone
  lens <- lens_sm
  cone <- cone_sm
  
  # - direction unit vec
  ucl <- lens[i_match,] - cone #use cone index
  # ucl <- lens_sm[i_match,] - cone_sm #use cone index
  ucl_norm <- sqrt(rowSums(ucl^2))
  ucl <- ucl / ucl_norm
  
  # - rotation to canonical coordinate
  
  # vn -> +z, vf -> +x, this way, equatorial ucl not on equatorial plane given lens positions
  mrot <- cbind(vf,vleft,vn)
  ucl_rot <- as.matrix(ucl) %*% mrot
  
  # # alt,  based on equator and ucl, note: can NOT add shift
  # ucl_pca <- prcomp(ucl[c(ind_Up2, ind_Down2), ]) #
  # vn2 <- ucl_pca$rotation[,3]
  # 
  # pt_theta90 <- as.data.frame(ucl[c(ind_Up2, ind_Down2), ])
  # lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z + 0)
  # vn2 <- c(-1, lm_plane$coefficients[1], lm_plane$coefficients[2]) #norm
  # vn2 <- vn2/sqrt(sum(vn2^2))
  # # which side is up
  # if (vn2 %*% c(0,0,1) < 0) {
  #   vn2 <- -vn2
  # }
  # 
  # vleft2 <- c(0,1,0) - c(vn2 %*% c(0,1,0)) * vn2
  # vleft2 <- vleft2/sqrt(sum(vleft2^2))
  # vf2 <- cross3D(vleft2, vn2)
  # 
  # mrot <- cbind(vf2,vleft2,vn2)
  # # ucl_rot2 <- sweep(as.matrix(ucl),2, ucl_pca$center, '-') %*% mrot 
  # ucl_rot2 <- as.matrix(ucl) %*% mrot
  # colnames(ucl_rot2) <- c('x','y','z')
  # 
  # # PLOT, unit sphere, rotated
  # nopen3d()
  # points3d(ucl_rot2[ind_left_cone, ], size = 5, color = 'magenta')
  # points3d(ucl_rot2[!ind_left_cone,], size = 5, color = 'green')
  # points3d(ucl_rot2[c(ind_Up2, ind_Down2),], size = 10, color = 'gray20')
  # spheres3d(0,0,0,0.99, col='grey', alpha=0.3, lit=F)
  # arrow3d(c(0,0,0), vn2*1, theta = pi/6,n = 4, col="blue", type = "rotation")
  # planes3d(0,0,1, 0, alpha = 0.3)
  # planes3d(0,1,0, 0, alpha = 0.3)
  # axes3d(c('x','y','z')); title3d('','','x','y','z')
  # title3d(fn[f])
  
  
  # -- smoothing (already done earlier)
  ucl_rot_sm <- ucl_rot
  
  # separate left and right ommatidia directions
  ii <- cone[,2] > 0 
  # left
  cc <- cone[ii,]
  ucl_rot_1 <- ucl_rot[ii, ]
  ls <- reghex(cc)
  nb_ind <- ls[[1]]
  ucl_rot_new <- ucl_rot_1
  for (m in 1:nrow(nb_ind)) {
    if (!any(is.na(nb_ind[m,]))) {
      ucl_rot_new[nb_ind[m,1],] <- ucl_rot_1[nb_ind[m,1],]/2 + colMeans(ucl_rot_1[nb_ind[m,-1],])/2
    }
  }
  ucl_rot_sm[ii, ] <- ucl_rot_new
  
  # right
  cc <- cone[!ii,]
  ucl_rot_1 <- ucl_rot[!ii, ]
  ls <- reghex(cc)
  nb_ind <- ls[[1]]
  ucl_rot_new <- ucl_rot_1
  for (m in 1:nrow(nb_ind)) {
    if (!any(is.na(nb_ind[m,]))) {
      ucl_rot_new[nb_ind[m,1],] <- ucl_rot_1[nb_ind[m,1],]/2 + colMeans(ucl_rot_1[nb_ind[m,-1],])/2
    }
  }
  ucl_rot_sm[!ii, ] <- ucl_rot_new
  
  # normalize
  ucl_rot_sm <- sweep(ucl_rot_sm,1,sqrt(rowSums(ucl_rot_sm^2)),'/') #normalize
  
  # save(lens, cone, ind_Up, ind_Down, ind_left_lens, ind_left_cone,
  #      eye_pca, vn, vf, vleft, L, i_match, ucl,
  #      vn2, vf2, vleft2, ucl_rot2,
  #      file = paste("../microCT/data/", fn[f], ".RData", sep='') )
  
  # note, ind_Up2/etc are corrected eq labels, not necessary if manually corr
  # but vn2/vf2/etc are different, they're values based on ucl
  save(lens, cone, ind_Up, ind_Down, ind_left_lens, ind_left_cone,
       eye_pca, vn, vf, vleft, ucl_rot, L, i_match, ucl,
       ind_Up2, ind_Down2, ucl_rot_sm,
       file = paste("data/microCT/", fn[f], ".RData", sep='') )
  
}

# proc 2021 eyemap -----------------------------------------------------------

# fn <- c('12102019_female', # match EM
#         '12062019_female',
#         '12052019_female',
#         '033021_male',  # no overlap
#         '033021_male_2', # ok
#         '033021_female', # bad edge
#         '033021_female_2', # large overlap
#         '033121_female') # bad edge
# 
# for (f in 1:length(fn)) {
#   lens <- read.csv(paste("../microCT/data/", fn[f], "_lens_Detailed.csv", sep=''), header=T) 
#   cone <- read.csv(paste("../microCT/data/", fn[f], "_tip_Detailed.csv", sep=''), header=T) 
#   
#   # extract
#   lens <- lens[, c('Position.X','Position.Y','Position.Z')]
#   colnames(lens) <- c("x", "y", "z")
#   lens <- lens[order(lens$x, lens$y,lens$z),]
#   rownames(lens) <- seq(1,length(lens$x), by = 1)
#   
#   cone <- cone[, c('Position.X','Position.Y','Position.Z','PR_tips_equator')]
#   colnames(cone) <- c("x", "y", "z", "label")
#   cone <- cone[order(cone$x, cone$y,cone$z),]
#   rownames(cone) <- seq(1,length(cone$x), by = 1)
#   ind_Down <- which(cone$label == "down")
#   ind_Up <- which(cone$label == "up")
#   cone <- cone[,c("x", "y", "z")]
#   
#   # - alignment
#   # equator
#   v_up <- colMeans(cone[ind_Up,]) - colMeans(cone[ind_Down,])
#   cone_equator <- cone[c(ind_Up, ind_Down),]
#   pt_theta90 <- cone_equator
#   # lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z)
#   # vn <- c(-1, lm_plane$coefficients[2], lm_plane$coefficients[3]) #norm
#   vn <- prcomp(pt_theta90)$rotation[,3]
#   vn <- vn/sqrt(sum(vn^2))
#   # which side is up
#   if (sum(vn * v_up) < 0) {
#     vn <- -vn
#   }
#   # pt_center <- colMeans(cone)
#   
#   # front via pca, use both eyes
#   eye_pca <- prcomp(rbind(lens, cone))
#   vf <- eye_pca$rotation[,3]
#   vf <- vf - c(vf %*% vn) * vn
#   vf <- vf / sqrt(sum(vf^2))
#   if (vf %*% (colMeans(lens) - eye_pca$center) < 0) {
#     vf <- - vf
#   }
#   
#   vleft <- cross3D(vn, vf)
#   
#   # # DEBUG
#   # nopen3d()
#   # points3d(cone)
#   # points3d(cone[ind_Up,], col='red', size =10)
#   # points3d(cone[ind_Down,], col='blue', size =10)
#   # # arrow3d(eye_pca$center, eye_pca$center + eye_pca$rotation[,2] * 100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   # arrow3d(eye_pca$center, eye_pca$center + vn * 100, theta = pi/6,n = 4, col="red", type = "rotation")
#   # arrow3d(eye_pca$center, eye_pca$center + vf * 100, theta = pi/6,n = 4, col="magenta", type = "rotation")
# 
#   
#   # - remove shifted equator pts
#   # centering
#   mrot <- cbind(vf,vleft,vn)
#   lens_tmp <- sweep(as.matrix(lens),2,eye_pca$center,'-') %*% mrot
#   cone_tmp <- sweep(as.matrix(cone),2,eye_pca$center,'-') %*% mrot
#   ind_Up2 <- ind_Up[cone_tmp[ind_Up,1] > 0]
#   ind_Down2 <- ind_Down[cone_tmp[ind_Down,1] > 0]
#   
#   # alignment again
#   v_up <- colMeans(cone[ind_Up2,]) - colMeans(cone[ind_Down2,])
#   cone_equator <- cone[c(ind_Up2, ind_Down2),]
#   pt_theta90 <- cone_equator
#   # lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z)
#   # vn <- c(-1, lm_plane$coefficients[2], lm_plane$coefficients[3]) #norm
#   vn <- prcomp(pt_theta90)$rotation[,3]
#   vn <- vn/sqrt(sum(vn^2))
#   # which side is up
#   if (sum(vn * v_up) < 0) {
#     vn <- -vn
#   }
#   # front via pca, use both eyes
#   eye_pca <- prcomp(rbind(lens, cone))
#   vf <- eye_pca$rotation[,3]
#   vf <- vf - c(vf %*% vn) * vn
#   vf <- vf / sqrt(sum(vf^2))
#   if (vf %*% (colMeans(lens) - eye_pca$center) < 0) {
#     vf <- - vf
#   }
#   # separate left and right
#   vleft <- cross3D(vn, vf)
#   ind_left_lens <- c(t(vleft) %*% t(sweep(lens,2,eye_pca$center,'-')) > 0)
#   ind_left_cone <- c(t(vleft) %*% t(sweep(cone,2,eye_pca$center,'-')) > 0)
#   
#   
#   # - matching, number of points SHOULD be equal
#   L <- dim(cone)[1]
#   mDist <- matrix(ncol = L, nrow = L)
#   dr <- cone
#   dc <- lens
#   for (j in 1:L) {
#     mDist[j,] <- as.matrix(dist(rbind(dr[j,], dc)))[-1,1]
#   }
#   i_match <- solve_LSAP(mDist)
#   
#   # PLOT
#   nopen3d()
#   points3d(lens[ind_left_lens,], size = 5, color = 'magenta')
#   points3d(lens[!ind_left_lens,], size = 5, color = 'green')
#   points3d(cone, size = 5, color = 'blue')
#   points3d(cone[c(ind_Down2, ind_Up2),], size=10, col='gold2')
#   points3d(lens[i_match[c(ind_Down2, ind_Up2)],], size=10, col='cyan')
#   arrow3d(eye_pca$center, eye_pca$center + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   # arrow3d(c(0,0,0), c(0,0,0) + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   planes3d(vn[1],vn[2],vn[3], -sum(eye_pca$center*vn), alpha = 0.5)
#   # planes3d(0,0,1, 0, alpha = 0.3)
#   arrow3d(eye_pca$center, eye_pca$center + vf*100, theta = pi/6,n = 4, col="brown", type = "rotation")
#   # arrow3d(c(0,0,0), c(0,0,0) + vf*100, theta = pi/6,n = 4, col="brown", type = "rotation")
#   # for (j in 1:L) {
#   #   lines3d(rbind(dr[j,], dc[i_match[j],]), col = "grey")
#   # }
#   linemat <- matrix(t(cbind(cone, lens[i_match,])), ncol = 3, byrow = T)
#   segments3d(linemat, color = "grey")
#   # title3d(fn[f])
# 
#   # index of lens
#   ind_Up_lens <- i_match[ind_Up2]
#   ind_Down_lens <- i_match[ind_Down2]
#   
#   # # DEBUG, equator viewing dir on a plane ?
#   # nopen3d()
#   # planes3d(vn[1],vn[2],vn[3], -sum(eye_pca$center*vn), alpha = 0.5)
#   # points3d(cone[c(ind_Up2, ind_Down2),], size=6)
#   # # points3d(lens[c(ind_Up_lens, ind_Down_lens),], size=6, col='blue')
#   # points3d(cone[c(ind_Up2, ind_Down2),] + 10*(lens[c(ind_Up_lens, ind_Down_lens),] - cone[c(ind_Up2, ind_Down2),]), size=8, col='red')
#   # arrow3d(eye_pca$center, eye_pca$center + vn*100, theta = pi/6,n = 4, col="blue", type = "rotation")
#   # title3d(fn[f])
#   
#   # - centering
#   mrot <- cbind(vf,vleft,vn)
#   lens <- sweep(as.matrix(lens),2,eye_pca$center,'-') %*% mrot
#   cone <- sweep(as.matrix(cone),2,eye_pca$center,'-') %*% mrot
#   colnames(lens) <- c('x','y','z')
#   colnames(cone) <- c('x','y','z')
#   # vf <- vf %*% mrot
#   # vleft <- c(vleft) %*% mrot
#   # vn <- vn %*% mrot
#   vf <- c(1,0,0)
#   vleft <- c(0,1,0)
#   vn <- c(0,0,1)
#   
#   # - direction unit vec
#   ucl <- lens[i_match,] - cone #use cone index
#   ucl_norm <- sqrt(rowSums(ucl^2))
#   ucl <- ucl / ucl_norm
#   
#   # rotation to canonical coordinate, this dose NOT make sense
#   
#   ### ### vn -> +z, vf -> +x, this way, equatorial ucl not on equatorial plane given lens positions
#   mrot <- cbind(vf,vleft,vn)
#   ucl_rot <- as.matrix(ucl) %*% mrot
#   
#   
#   ### ### based on equator and ucl
#   # ucl_pca <- prcomp(ucl[c(ind_Up2, ind_Down2), ]) # can NOT do shift
#   # vn2 <- ucl_pca$rotation[,3]
#   
#   pt_theta90 <- as.data.frame(ucl[c(ind_Up2, ind_Down2), ])
#   lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z + 0)
#   vn2 <- c(-1, lm_plane$coefficients[1], lm_plane$coefficients[2]) #norm
#   vn2 <- vn2/sqrt(sum(vn2^2))
#   # which side is up
#   if (vn2 %*% c(0,0,1) < 0) {
#     vn2 <- -vn2
#   }
#   
#   vleft2 <- c(0,1,0) - c(vn2 %*% c(0,1,0)) * vn2
#   vleft2 <- vleft2/sqrt(sum(vleft2^2))
#   vf2 <- cross3D(vleft2, vn2)
#   
#   mrot <- cbind(vf2,vleft2,vn2)
#   # ucl_rot2 <- sweep(as.matrix(ucl),2, ucl_pca$center, '-') %*% mrot 
#   ucl_rot2 <- as.matrix(ucl) %*% mrot
#   colnames(ucl_rot2) <- c('x','y','z')
#   
#   # ucl_rot2[,2] <- -ucl_rot2[,2] # right eye, phi angle viewed from outside,
#   
#   # # PLOT, unit sphere, rotated
#   # nopen3d()
#   # points3d(ucl_rot2[ind_left_cone, ], size = 5, color = 'magenta')
#   # points3d(ucl_rot2[!ind_left_cone,], size = 5, color = 'green')
#   # points3d(ucl_rot2[c(ind_Up2, ind_Down2),], size = 10, color = 'gray20')
#   # spheres3d(0,0,0,0.99, col='grey', alpha=0.3, lit=F)
#   # arrow3d(c(0,0,0), vn2*1, theta = pi/6,n = 4, col="blue", type = "rotation")
#   # planes3d(0,0,1, 0, alpha = 0.3)
#   # planes3d(0,1,0, 0, alpha = 0.3)
#   # axes3d(c('x','y','z')); title3d('','','x','y','z')
#   # title3d(fn[f])
#   
#   # save(lens, cone, ind_Up, ind_Down, ind_left_lens, ind_left_cone,
#   #      eye_pca, vn, vf, vleft, L, i_match, ucl,
#   #      vn2, vf2, vleft2, ucl_rot2,
#   #      file = paste("../microCT/data/", fn[f], ".RData", sep='') )
#   
#   save(lens, cone, ind_Up, ind_Down, ind_left_lens, ind_left_cone,
#        eye_pca, vn, vf, vleft, ucl_rot, L, i_match, ucl,
#        ind_Up2, ind_Down2, vn2, vf2, vleft2, ucl_rot2,
#        file = paste("../microCT/data/eq_corr/", fn[f], ".RData", sep='') )
# }


# nb_ind and ioa, 20240701 ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240701.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
# ls <- reghex(pt, j7 = c(427, 448, 462, 440, 406, 393, 413) ) #[home, +p,+v,+q,-p,-v,-q] 
ls <- reghex(pt, lefteye = T )
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
# ls <- reghex(pt, j7 = c(406, 430, 441, 415, 382, 372, 395))
ls <- reghex(pt, lefteye = F)
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
nb_ind <- nb_ind_left
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1],acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}
nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1], acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) 
  )
}
nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240701_nb.RData") )


# nb_ind and ioa, 20240612 ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240612.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(391, 410, 425, 406, 371, 356, 377) ) #[home, +p,+v,+q,-p,-v,-q] 
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(375, 398, 409, 384, 352, 342, 365))
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
nb_ind <- nb_ind_left
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1],acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}
nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1], acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) 
  )
}
nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240612_nb.RData") )


# nb_ind and ioa, 20240530_dried ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240530_dried.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(464, 487, 501, 474, 436, 427, 453) ) #[home, +p,+v,+q,-p,-v,-q] 
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(328, 332, 361, 356, 323, 295, 298))
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
nb_ind <- nb_ind_left
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1],acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}
nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1], acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) 
  )
}
nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240530_dried_nb.RData") )


# nb_ind and ioa, 20240530 ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240530.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(369, 389, 401, 381, 349, 335, 355) ) #[home, +p,+v,+q,-p,-v,-q] 
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(373, 394, 407, 385, 353, 339, 364))
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
nb_ind <- nb_ind_left
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1],acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}
nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1], acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) 
  )
}
nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240530_nb.RData") )

# nb_ind and ioa, 20240524 ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240524.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(428, 447, 464, 446, 413, 392, 412) ) #[home, +p,+v,+q,-p,-v,-q] 
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(411, 425, 447, 435, 400, 377, 389))
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
# dd <- sweep(dd, 1, sqrt(rowSums(dd^2)), '/')
nb_ind <- nb_ind_left

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1],acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}

nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  # nb_dist_ucl[j,] <- c(nb_ind[j,1],
  #                      sqrt(rowSums(sweep(dd[nb_ind[j,-1], ], 2, dd[nb_ind[j,1],,drop=F])^2)) )
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1], acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) 
  )
}

nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240524_nb.RData") )


# nb_ind and ioa, 20240522 ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240522.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(331, 347, 365, 350, 317, 300, 314) ) #[home, +p,+v,+q,-p,-v,-q]
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

# only for this scan, remove that extra lens-cone
# ii <- which(!(seq(1,780) %in% nb_ind_left[,1]))
# nb_ind_left <- nb_ind_left[!is.na(nb_ind_left[,1]),]
# ind_xy_left <- ind_xy_left[!is.na(ind_xy_left[,1]),]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(364, 377, 396, 378, 346, 329, 345) )
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
# dd <- sweep(dd, 1, sqrt(rowSums(dd^2)), '/')
nb_ind <- nb_ind_left
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  if (!is.na(nb_ind[j,1])) { #for na rows
    nb_dist_ucl[nb_ind[j,1],] <- c(
      nb_ind[j,1],acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
    )
  }
}
nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  # nb_dist_ucl[j,] <- c(nb_ind[j,1],
  #                      sqrt(rowSums(sweep(dd[nb_ind[j,-1], ], 2, dd[nb_ind[j,1],,drop=F])^2)) )
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1], acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) 
  )
}
nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240522_nb.RData") )


# nb_ind and ioa, 20240513 ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240520.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(370, 391, 406, 387, 350, 335, 355) ) #[home, +p,+v,+q,-p,-v,-q] 20240125
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(353, 368, 388, 375, 340, 317, 333) ) # 20240125
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
# dd <- sweep(dd, 1, sqrt(rowSums(dd^2)), '/')
nb_ind <- nb_ind_left

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1],acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}

nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  # nb_dist_ucl[j,] <- c(nb_ind[j,1],
  #                      sqrt(rowSums(sweep(dd[nb_ind[j,-1], ], 2, dd[nb_ind[j,1],,drop=F])^2)) )
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1], acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) 
  )
}

nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240520_nb.RData") )

# nb_ind and ioa, 20240513 ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240513.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(363, 370, 397, 390, 355, 329, 338) ) #[home, +p,+v,+q,-p,-v,-q] 20240125
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(342, 362, 377, 358, 321, 307, 325) ) # 20240125
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
# dd <- sweep(dd, 1, sqrt(rowSums(dd^2)), '/')
nb_ind <- nb_ind_left

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1],acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}

nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  # nb_dist_ucl[j,] <- c(nb_ind[j,1],
  #                      sqrt(rowSums(sweep(dd[nb_ind[j,-1], ], 2, dd[nb_ind[j,1],,drop=F])^2)) )
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1], acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) 
    )
}

nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240513_nb.RData") )

# nb_ind and ioa, 20240510 ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240510.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(469, 485, 501, 483, 447, 433, 449) ) #[home, +p,+v,+q,-p,-v,-q] 
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(378, 391, 414, 400, 366, 344, 358) ) 
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
# dd <- sweep(dd, 1, sqrt(rowSums(dd^2)), '/')
nb_ind <- nb_ind_left

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1],acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}

nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  # nb_dist_ucl[j,] <- c(nb_ind[j,1],
  #                      sqrt(rowSums(sweep(dd[nb_ind[j,-1], ], 2, dd[nb_ind[j,1],,drop=F])^2)) )
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1], acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) 
  )
}

nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240510_nb.RData") )

# nb_ind and ioa, 20240206 ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240206.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(388, 400, 423, 409, 375, 352, 364) ) #[home, +p,+v,+q,-p,-v,-q] 20240125
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(434, 457, 468, 446, 410, 401, 421) ) # 20240125
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
# dd <- sweep(dd, 1, sqrt(rowSums(dd^2)), '/')
nb_ind <- nb_ind_left

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(nb_ind[j,1],
                                 acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}

nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  # nb_dist_ucl[j,] <- c(nb_ind[j,1],
  #                      sqrt(rowSums(sweep(dd[nb_ind[j,-1], ], 2, dd[nb_ind[j,1],,drop=F])^2)) )
  nb_dist_ucl[nb_ind[j,1],] <- c(nb_ind[j,1],
                                 acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) )
}

nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240206_nb.RData") )

# nb_ind and ioa, 20240125 ------------------------------------------------

load(paste0("../microCT/2023_eyemap/20240125.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(234, 248, 211, 200, 228, 261, 273) ) #[home, +p,+v,+q,-p,-v,-q] 20240125
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(310, 325, 284, 272, 296, 337, 349) ) # 20240125
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
# dd <- sweep(dd, 1, sqrt(rowSums(dd^2)), '/')
nb_ind <- nb_ind_left

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(nb_ind[j,1],
                                 acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}

nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  # nb_dist_ucl[j,] <- c(nb_ind[j,1],
  #                      sqrt(rowSums(sweep(dd[nb_ind[j,-1], ], 2, dd[nb_ind[j,1],,drop=F])^2)) )
  nb_dist_ucl[nb_ind[j,1],] <- c(nb_ind[j,1],
                                 acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) )
}

nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20240125_nb.RData") )

# nb_ind and ioa, 20231107 ----------------------------------------------------------

load(paste0("../microCT/2023_eyemap/20231107.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(497, 504, 536, 526, 490, 461, 468) ) #[home, +p,+v,+q,-p,-v,-q] 20231107
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(395, 388, 424, 430, 401, 366, 360) ) # 2023107
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]


# # DEBUG
# nopen3d()
# points3d(pt,size=7)
# text3d(pt[ind_xy_right[,1],], texts = paste(ind_xy_right[,2],ind_xy_right[,3],sep='|'), adj=1.3)
# points3d(pt[ind_xy_right[ind_xy_right[,3] %% 3 == 0,1],], size=10, col='magenta')

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
# dd <- sweep(dd, 1, sqrt(rowSums(dd^2)), '/')
nb_ind <- nb_ind_left

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(nb_ind[j,1],
                                 acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}

nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  # nb_dist_ucl[j,] <- c(nb_ind[j,1],
  #                      sqrt(rowSums(sweep(dd[nb_ind[j,-1], ], 2, dd[nb_ind[j,1],,drop=F])^2)) )
  nb_dist_ucl[nb_ind[j,1],] <- c(nb_ind[j,1],
                                 acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) )
}

nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20231107_nb.RData") )

# nb_ind and ioa, 20230926 ----------------------------------------------------------

load(paste0("../microCT/2023_eyemap/20230926.RData"))

# - hex coord and ioa
pt <- lens[ind_left_lens, ]
ls <- reghex(pt, j7 = c(415, 403, 441, 453, 426, 388, 380) ) #[home, +p,+v,+q,-p,-v,-q]
nb_ind_left <- ls[[1]]
ind_xy_left <- ls[[2]]

pt <- lens[!ind_left_lens, ]
ls <- reghex(pt, j7 = c(647, 693, 697, 648, 603, 600, 646) )
nb_ind_right <- ls[[1]]
ind_xy_right <- ls[[2]]


# # DEBUG
# nopen3d()
# points3d(pt,size=7)
# text3d(pt[ind_xy_right[,1],], texts = paste(ind_xy_right[,2],ind_xy_right[,3],sep='|'), adj=1.3)
# points3d(pt[ind_xy_right[ind_xy_right[,3] %% 3 == 0,1],], size=10, col='magenta')

# ioa, left
dd <- ucl_rot_sm[order(i_match),][ind_left_lens,]
# dd <- sweep(dd, 1, sqrt(rowSums(dd^2)), '/')
nb_ind <- nb_ind_left

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd)) #arc dist
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(nb_ind[j,1],
                                 acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],])
  )
}

nb_dist_ucl_left <- nb_dist_ucl

# right
dd <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
nb_ind <- nb_ind_right

nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  # nb_dist_ucl[j,] <- c(nb_ind[j,1],
  #                      sqrt(rowSums(sweep(dd[nb_ind[j,-1], ], 2, dd[nb_ind[j,1],,drop=F])^2)) )
  nb_dist_ucl[nb_ind[j,1],] <- c(nb_ind[j,1],
                                 acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) )
}

nb_dist_ucl_right <- nb_dist_ucl

# SAVE
save(nb_ind_left, nb_ind_right, ind_xy_left, ind_xy_right,
     nb_dist_ucl_left, nb_dist_ucl_right,
     file = paste0("../microCT/2023_eyemap/20230926_nb.RData") )

# direction based on normals, only right side, 2023 ----------------------------------------------

f <- 2

load(paste0("../microCT/2023_eyemap/", fn[f], ".RData"))
load(paste0("../microCT/2023_eyemap/", fn[f], "_nb.RData"))

pts <- lens[!ind_left_lens,]
nb_ind <- nb_ind_right
ucl <- ucl_rot_sm[order(i_match),][!ind_left_lens,]

# - calculate normals
nucl <- matrix(nrow = nrow(pts), ncol = 3)
for (j in 1:nrow(nucl)) {
  if (sum(is.na(nb_ind[j,])) == 0) {
    pca <- prcomp(pts[nb_ind[j,],])
    vv <- pca$rotation[,3]
    if (vv %*% ucl[nb_ind[j,1],] > 0) {
      nucl[nb_ind[j,1], ] <- vv
    } else {
      nucl[nb_ind[j,1], ] <- -vv
    }
  }
}
save(nucl, file = paste("../microCT/2023_eyemap/", fn[f], "_normals.RData", sep='') )



# # direction based on normals, 2021 ----------------------------------------------
# # load("../microCT/data/eq_corr/12102019_female_new.RData")
# 
# f <- 1
# for (f in 1:length(fn)) {
#   # load(paste("../microCT/data/eq_corr/", fn[f], ".RData", sep=''))  
#   load(paste("../microCT/2023_eyemap/", fn[f], ".RData", sep=''))
#   
#   pts <- lens[!ind_left_lens,]
#   
#   # rotate to southern hemisphere
#   vcom <- colMeans(pts)
#   vcom <- vcom / sqrt(sum(vcom^2))
#   R_mat <- quaternion3D( cross3D(vcom, c(0,0,-1)), acos(vcom %*% c(0,0,-1)) /pi*180 )
#   pts_z <- as.matrix(pts) %*% t(R_mat)
#   
#   # stereographic projection
#   x <- pts_z[,1] / (1 - pts_z[,3])
#   y <- pts_z[,2] / (1 - pts_z[,3])
#   stereo <- cbind(x, y) 
#   stereo_df <- as.data.frame(stereo)
#   
#   # - deldir
#   # dtdt <- deldir(stereo[,1], stereo[,2], rw = c(-0.8,0.75,-1.1,0.8)) # mess up index ?
#   dtdt <- deldir(stereo[,1], stereo[,2])
#   
#   # -- deal with edge points, 
#   # if the boundary pt, then check the 2 circumcenters sharing this edge,
#   # remove the one with no corresponding circumcenter
#   # Alt,  deal with edge points wutg dummy pt. NOT tried ?
#   
#   ee_mat <- matrix(c(1,1,1,2,2,1,2,2), ncol = 2, byrow = T) #index matrix
#   
#   bp <- which(dtdt$dirsgs$bp1 | dtdt$dirsgs$bp2)
#   
#   delsgs_mod <- dtdt$delsgs # modify
#   for (j in bp) {
#     ind_edge <- dtdt$dirsgs[j, c('ind1', 'ind2')] #edge vertices 
#     ind1_tri <- c(delsgs_mod[delsgs_mod$ind1 %in% ind_edge[1], ]$ind2,
#                   delsgs_mod[delsgs_mod$ind2 %in% ind_edge[1], ]$ind1)
#     ind2_tri <- c(delsgs_mod[delsgs_mod$ind1 %in% ind_edge[2], ]$ind2,
#                   delsgs_mod[delsgs_mod$ind2 %in% ind_edge[2], ]$ind1)
#     ind_tri <- c(ind1_tri, ind2_tri)
#     ind_share <- ind_tri[duplicated(ind_tri)] # two vertices shared by the edge vertices
#     
#     if (length(ind_share) == 2) {
#       ee <- c() #exist a circumcenter
#       for (m in 1:2) {
#         for (n in 1:2) {
#           b1 <- ind_edge[m] %in% dtdt$dirsgs[dtdt$dirsgs$ind1 %in% ind_share[n], ]$ind2
#           b2 <- ind_edge[m] %in% dtdt$dirsgs[dtdt$dirsgs$ind2 %in% ind_share[n], ]$ind1
#           ee <- c(ee, b1 | b2)
#         }
#       }
#       if (sum(ee) < 4) {
#         for (k in which(!ee)) {
#           # to remove
#           ind1 <- delsgs_mod$ind1 %in% ind_edge[ee_mat[k,][1]] & delsgs_mod$ind2 %in% ind_share[ee_mat[k,][2]] 
#           ind2 <- delsgs_mod$ind2 %in% ind_edge[ee_mat[k,][1]] & delsgs_mod$ind1 %in% ind_share[ee_mat[k,][2]] 
#           ind_rm <- ind1 | ind2
#           delsgs_mod <- delsgs_mod[-which(ind_rm), ]
#         }
#       }
#     }
#   }
#   
#   # remove long edge
#   el <- pts[delsgs_mod$ind1,] - pts[delsgs_mod$ind2,]
#   el <- sqrt(rowSums(el^2))
#   delsgs_mod <- delsgs_mod[el < 30, ]
#   
#   # - calculate normals
#   nucl <- matrix(nrow = nrow(pts), ncol = 3)
#   for (j in 1:nrow(nucl)) {
#     jj <- c(delsgs_mod$ind2[delsgs_mod$ind1 == j], delsgs_mod$ind1[delsgs_mod$ind2 == j] )
#     nbs <- pts[c(j, jj),]
#     lm_plane <- lm(nbs[,2] ~ nbs[,3] + nbs[,1] )
#     vn <- c(lm_plane$coefficients[3], -1, lm_plane$coefficients[2]) #norm
#     vn <- vn/sqrt(sum(vn^2))
#     if (sum(vn * pts[j,]) < 0) {
#       vn <- -vn
#     }
#     nucl[j,] <- vn
#   }
#   
#   save(lens, cone, ind_Up, ind_Down, ind_left_lens, ind_left_cone,
#        eye_pca, vn, vf, vleft, ucl_rot_orig, L, i_match, ucl,
#        ind_Up2, ind_Down2, ucl_rot, nucl,
#        file = paste("../microCT/2023_eyemap/", fn[f], ".RData", sep='') )
#   
#   # save(lens, cone, ind_Up, ind_Down, ind_left_lens, ind_left_cone,
#   #      eye_pca, vn, vf, vleft, ucl_rot, L, i_match, ucl,
#   #      ind_Up2, ind_Down2, vn2, vf2, vleft2, ucl_rot2, nucl,
#   #      file = paste("../microCT/data/eq_corr/", fn[f], ".RData", sep='') )
#   
#   # save(lens, cone, ind_Up, ind_Down, ind_left_lens, ind_left_cone,
#   #      eye_pca, vn, vf, vleft, ucl_rot, L, i_match, ucl,
#   #      ind_Up2, ind_Down2, vn2, vf2, vleft2, ucl_rot2, nucl,
#   #      file = "../microCT/data/eq_corr/12102019_female_new.RData" )
# }
# 
# # #  PLOT
# # dev.new()
# # ggplot() +
# #   geom_segment(data = delsgs_mod, aes(x = x1, y = y1, xend = x2, yend = y2),
# #                size = 1, linetype = 1, color= "black") +
# #   geom_segment(data = dtdt$dirsgs, aes(x = x1, y = y1, xend = x2, yend = y2),
# #                size = 1, linetype = 1, color= "red") +
# #   theme_bw()
# # # coord_fixed(ratio = 1, xlim = c(-0.8, 0.75), ylim = c(-1.1, 0.8))
# 
# nopen3d()
# points3d(pts)
# linemat <- matrix(t(cbind(pts[delsgs_mod$ind1,],pts[delsgs_mod$ind2,])),
#                   ncol = 3, byrow = T)
# segments3d(linemat, color = "grey")
# 
# for (j in 1:nrow(ucl_rot)) {
#   arrow3d(pts[j,], pts[j,]+ucl_rot[j,]*10, theta=pi/18, n=8, s=0.3,width=0.3, col = "gray30", type = "rotation")
#   arrow3d(pts[j,], pts[j,]+nucl[j,]*10, theta=pi/18, n=8, s=0.3,width=0.3, col = "cyan", type = "rotation")
# }





# analysis ----------------------------------------------------------------

tb_fov <- matrix(ncol = 3, nrow = length(fn))
plt_pt_ls <- list()
plt_kde_ls <- list()
Zpred_lim <- matrix(ncol = 2, nrow = length(fn))
N_lens <- matrix(ncol = 2, nrow = length(fn))

for (f in 1:length(fn)) {
  # load(paste("../microCT/data/eq_corr/", fn[f], ".RData", sep=''))  
  load(paste("../microCT/2023_eyemap/", fn[f], ".RData", sep=''))  
  
  N_lens[f,] <- c(sum(ind_left_cone), sum(!ind_left_cone) )
  
  # - Mollweide
  ucl_rot_tp <- cart2sphZ(ucl_rot_sm)[,2:3] %>%
    as_tibble() %>%
    mutate(theta = theta/pi*180, phi = phi/pi*180) %>%
    mutate(phi = if_else(phi > 180, phi - 360, phi)) %>%
    mutate(phi = - phi) %>% #inside out
    as.matrix()
  ucl_rot_Mo <- Mollweide(ucl_rot_tp)
  ucl_rot_Mo <- cbind(ucl_rot_Mo, 1)
  ucl_rot_Mo[!ind_left_cone, 3] <- 2 # left == 1, right == 2
  colnames(ucl_rot_Mo) <- c('xM', 'yM', 'left')
  
  df_left <- ucl_rot_Mo %>% as_tibble() %>% filter(left == 1) %>% as.data.frame()
  df_right <- ucl_rot_Mo %>% as_tibble() %>% filter(left == 2) %>% as.data.frame()
  
  # - make ashape boundary
  xy_ashape <- ashape(ucl_rot_Mo[ind_left_cone, 1:2] +
                        matrix(runif(sum(ind_left_cone)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
  xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
  ahull_left <- xy_grid_ahull
  
  xy_ashape <- ashape(ucl_rot_Mo[!ind_left_cone, 1:2] +
                        matrix(runif(sum(!ind_left_cone)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
  xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
  ahull_right <- xy_grid_ahull
  
  # - point distr
  # df <- as.data.frame(ucl_rot_Mo[ind_left_cone, 1:3])
  df <- as.data.frame(ucl_rot_Mo)
  plt <- ggplot(df, aes(x = xM, y = yM)) +
    # ggplot(df, aes(x = xM, y = yM)) +
    geom_point(aes(colour = factor(left))) +
    scale_colour_manual(values = c("magenta", "green"), labels=c('left','right')) +
    # geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), 
    #              arrow = arrow(length = unit(0.005, "npc"),type = "closed"), size =1) +
    # annotate("text", x = 2.3, y = -1.1, label = "trans", size = 5, stroke = 3) +
    geom_path(data = ahull_left, aes(x1,y1), lwd = 1.5, colour = "magenta", alpha =0.7) +
    geom_path(data = ahull_right, aes(x1,y1), lwd = 1.5, colour = "green", alpha =0.7) +
    geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
    ylab("elevation") +
    xlab("azimuth") +
    theme_minimal() +
    scale_x_continuous(limits = c(-sqrt(8), sqrt(8)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-sqrt(2), sqrt(2)), breaks = c(-1.5,0,1.5), labels = c(-1.5,0,1.5), expand = c(0, 0)) + # set +y as above eq
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14) ) +
    labs(title = fn[f]) +
    # theme_void() +
    # theme(legend.position="none", panel.background = element_blank()) +
    coord_fixed(ratio = 1, xlim = c(-sqrt(8), sqrt(8)), ylim = c(-sqrt(2), sqrt(2)))
  
  plt_pt_ls[[f]] <- plt
  
  # windows(record = F, width = 16, height = 8)
  # print(plt)
  # file_name = paste(fn[f], "_pt.png", sep="")
  # png(file_name, width = 1600, height = 800, pointsize = 12)
  # dev.off()
  
  # # density plot with geom_density_2d
  # windows(record = F, width = 16, height = 8)
  # ggplot() +
  #   # geom_point() +
  #   geom_density_2d_filled(data=df_left, aes(xM, yM), alpha = 0.5, adjust=0.6) +
  #   # scale_fill_manual(values = pal_left, guide_legend("den"),
  #   # labels = breaks) +
  #   # scale_fill_viridis_d() +
  #   geom_density_2d_filled(data=df_right, aes(xM, yM), alpha = 0.5, adjust=0.6) +
  #   scale_fill_manual(values = pal_right, guide_legend("den"),
  #                     labels = breaks) +
  #   geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
  #   geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
  #   geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
  #   geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
  #   theme_minimal() +
  #   coord_fixed(ratio = 1)
  
  
  # - Gaussian density on S2
  grid_Gaussian_ls <- list()
  for (j in 1:2) { #left and right
    ii <- which(ucl_rot_Mo[,'left'] == j)
    grid_Gaussian <- as.data.frame(S2grid)
    grid_Gaussian$arcL = 0
    grid_Gaussian$Z = 0
    for (k in ii) {
      t0 <- ucl_rot_tp[k, 'theta'] / 180 * pi
      p0 <- ucl_rot_tp[k, 'phi'] / 180 * pi
      tp0 <- c(t0,p0)
      A <- 1
      # Let half-width at 10% = 5 deg ==> if -x^2/s^2, then s^2 = (5/180*pi)^2 / log(10)
      # need arc length here
      grid_Gaussian$arcL <- apply(grid_Gaussian, 1, function(x) arcLength(tp0, x[1:2]/180*pi) )
      grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[6] + 1*A*exp(-x[5]^2 / ((5/180*pi)^2 / log(10))) )
    }
    grid_Gaussian_ls[[j]] <- grid_Gaussian[, c('xM','yM','Z')]
  }
  # eye_den <- rbind(grid_Gaussian_ls[[1]], grid_Gaussian_ls[[2]])
  # eye_den$left <- c(rep(1,nrow(grid_Gaussian_ls[[1]])), rep(2,nrow(grid_Gaussian_ls[[2]])))


  # - pal
  Zpred_lim[f,] <- sapply(grid_Gaussian_ls, function(x) max(x$Z))

  n_lvl <- 10
  breaks <- seq(0,5,length.out = n_lvl+1)
  # getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
  getPalette <- colorRampPalette(brewer.pal(9, "Purples"))
  pal_left <- getPalette(n_lvl )
  # getPalette <- colorRampPalette(brewer.pal(9, "Purples"))
  getPalette <- colorRampPalette(brewer.pal(9, "Greens"))
  pal_right <- getPalette(n_lvl )

  # -- add two color scales
  pal_sum <- rep(0,n_lvl*(n_lvl+2))
  for (j in seq((n_lvl+1), n_lvl*(n_lvl+1), length.out = n_lvl)) {
    pal_sum[j] <- pal_left[j/(n_lvl+1)]
    c1 <- col2rgb(pal_left[j/(n_lvl+1)], alpha = F)
    for (k in 1:n_lvl) {
      pal_sum[k] <- pal_right[k]
      c2 <- col2rgb(pal_right[k], alpha = F)
      c3 <- (c1 + c2)/2
      pal_sum[j+k] <- rgb(t(c3), maxColorValue = 255)
    }
  }

  # -- grid for loess
  grid_M <- expand.grid(xM = seq(-sqrt(8), sqrt(8), length.out = 100),
                        yM = seq(-sqrt(2), sqrt(2), length.out = 50) )
  ii_inpoly <- sp::point.in.polygon(grid_M[,1], grid_M[,2], bkgd_chull[,1], bkgd_chull[,2])

  # -- loess
  fit_loess <- loess(Z ~ xM * yM, data = grid_Gaussian_ls[[1]], degree = 2, span =0.1,
                     control = loess.control(surface = "direct"))
  pred_loess <- predict(fit_loess, grid_M, se = T)
  df_pred <- grid_M
  df_pred$Z <- melt(pred_loess$fit)$value
  # df_pred$Z[df_pred$Z < 0] <- 0
  df_pred$equalSpace <- cut(df_pred$Z, breaks)
  df_pred_left <- df_pred[ii_inpoly == 1, ]

  fit_loess <- loess(Z ~ xM * yM, data = grid_Gaussian_ls[[2]], degree = 2, span =0.1,
                     control = loess.control(surface = "direct"))
  pred_loess <- predict(fit_loess, grid_M, se = T)
  df_pred <- grid_M
  df_pred$Z <- melt(pred_loess$fit)$value
  df_pred$equalSpace <- cut(df_pred$Z, breaks)
  df_pred_right <- df_pred[ii_inpoly == 1, ]

  df_pred_left$colSpace <- as.numeric(df_pred_left$equalSpace)
  df_pred_left$colSpace[is.na(df_pred_left$colSpace)] <- 0
  df_pred_left$colSpace[df_pred_left$colSpace == 1] <- 0
  df_pred_left$colSpace <- df_pred_left$colSpace * (n_lvl+1)

  df_pred_right$colSpace <- as.numeric(df_pred_right$equalSpace)
  df_pred_right$colSpace[is.na(df_pred_right$colSpace)] <- 0
  df_pred_right$colSpace[df_pred_right$colSpace == 1] <- 0

  df_pred_left$Z[df_pred_left$Z < 0] <- 0
  df_pred_right$Z[df_pred_right$Z < 0] <- 0

  df_pred <- df_pred_left[, c('xM','yM')]
  df_pred$Z <- df_pred_left$Z + df_pred_right$Z
  df_pred$colSpace <- df_pred_left$colSpace + df_pred_right$colSpace
  df_pred$colSpace[df_pred$colSpace == 0] <- NA



  # # -- hsv colorspace manipulations
  # library(colorspace)
  #
  # desat <- function(cols, sat=0.5) {
  #   X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  #   hsv(X[1,], X[2,], X[3,])
  # }
  #
  # pal_left <- desat(pal_left, sat=0.5)
  # pal_right <- desat(pal_right, sat=0.5)

  # pie(rep(1,length(pal_left)), col=pal_left)

  # library(ggnewscale)

  # -- plot
  plt <- ggplot() +
    geom_raster(data=df_pred, aes(xM, yM, fill = factor(colSpace)), interpolate = F) +
    # geom_point(data=df_pred, aes(xM, yM)) +
    scale_fill_manual(values = pal_sum, breaks = seq(1, n_lvl*(n_lvl+2), by = 1)) +
    # geom_raster(data=df_pred_left, aes(xM, yM, fill = equalSpace), interpolate = F) +
    # scale_fill_manual(values = pal_left) +
    # new_scale('fill') +
    # geom_raster(data=df_pred_right, aes(xM, yM, fill = equalSpace), interpolate = F) +
    # scale_fill_manual(values = pal_right) +
    # # geom_point(data=df_left, aes(xM,yM), colour='#74c476') +
    # geom_point(data=df_left, aes(xM,yM), colour='#fdd0a2') +
    # # geom_contour(data=df_pred_left, aes(x=xM,y=yM, z=Z), breaks = breaks[-c(1)], color = "black", alpha = 0.5) +
    # # geom_point(data=df_right, aes(xM,yM), colour='#fdd0a2') +
    # geom_point(data=df_right, aes(xM,yM), colour='gray') +
    # geom_contour(data=df_pred_right, aes(x=xM,y=yM, z=Z), breaks = breaks[-c(1)], color = "black", alpha = 0.5) +
  geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
    theme_void() +
    theme(legend.position="none", panel.background = element_blank()) +
    labs(title = fn[f]) +
    coord_fixed(ratio = 1, xlim = c(-sqrt(8), sqrt(8)), ylim = c(-sqrt(2), sqrt(2)))

  plt_kde_ls[[f]] <- plt

  # windows(width = 16, height = 8)
  # print(plt)
  # file_name = paste(fn[f], "_kde.png", sep="")
  # png(file_name, width = 1600, height = 800, pointsize = 12)
  # dev.off()
  
  
  # - fov area 
  
  #rotation for fov calculation, so left side is southern hemisphere
  vleft <- cross3D(vn, vf)
  mrot2 <- cbind(vf, -vn, vleft)
  ucl_rot <- as.matrix(ucl) %*% mrot2 
  
  # Mollweide
  ucl_rot2_tp <- cart2sphZ(ucl_rot)[,2:3] %>%
    as_tibble() %>%
    mutate(theta = theta/pi*180, phi = phi/pi*180) %>%
    mutate(phi = if_else(phi > 180, phi - 360, phi)) %>%
    mutate(phi = - phi) %>% #inside out
    as.matrix()
  ucl_rot2_Mo <- Mollweide(ucl_rot2_tp)
  ucl_rot2_Mo <- cbind(ucl_rot2_Mo, 1)
  ucl_rot2_Mo[!ind_left_cone, 3] <- 2 # left == 1, right == 2
  colnames(ucl_rot2_Mo) <- c('xM', 'yM', 'left')
  
  
  # # PLOT
  # df <- as.data.frame(ucl_rot2_Mo)
  # windows(record = F, width = 16, height = 8)
  # ggplot(df, aes(x = xM, y = yM)) +
  #   geom_point(aes(colour = factor(left))) +
  #   scale_colour_manual(values = c("magenta", "green"), labels=c('left','right')) +
  #   geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
  #   geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
  #   geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
  #   geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
  #   theme_void() +
  #   theme(legend.position="none", panel.background = element_blank()) +
  #   coord_fixed(ratio = 1, xlim = c(-sqrt(8), sqrt(8)), ylim = c(-sqrt(2), sqrt(2)))
  
  # # chull
  # hpts <- chull(ucl_rot2_Mo[ucl_rot2_Mo[,'left'] == 1, 1:2])
  # hpts <- c(hpts, hpts[1])
  # hpts_xy <- ucl_rot2_Mo[ucl_rot2_Mo[,'left'] == 1, 1:2][hpts, ] # hull edge points
  # polygon(hpts_xy)
  
  # ahull
  xy_ashape <- ashape(ucl_rot2_Mo[ucl_rot2_Mo[,'left'] == 1, 1:2] + 
                        matrix(runif(sum(ucl_rot2_Mo[,'left'] == 1)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
  
  # xy_grid_ahull_ls[[1]] TODO
  
  # dev.new()
  # plot(ucl_rot2_Mo[,1:2])
  # points(xy_grid_ahull, col='red', cex=2)
  # polygon(xy_grid_ahull)
  
  # ii_inpoly <- sp::point.in.polygon(S2grid_a[,'xM'], S2grid_a[,'yM'],
  #                                   xy_grid_ahull[,1], xy_grid_ahull[,2])
  # points(S2grid_a[ii_inpoly == 1, c('xM','yM')], col='blue', cex=1, pch=16)
  # sum(ii_inpoly) / nrow(S2grid_a)
  
  # fov
  ii_inpoly_left <- sp::point.in.polygon(S2grid[,'xM'], S2grid[,'yM'],
                                         xy_grid_ahull[,1], xy_grid_ahull[,2])
  # points(S2grid[ii_inpoly == 1, c('xM','yM')], col='purple', cex=1, pch=16)
  
  
  a_left <- sum(ii_inpoly_left) / nrow(S2grid) *4*pi
  
  # right
  xy_ashape <- ashape(ucl_rot2_Mo[ucl_rot2_Mo[,'left'] == 2, 1:2], alpha = 0.3)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
  ii_inpoly_right <- sp::point.in.polygon(S2grid[,'xM'], S2grid[,'yM'],
                                          xy_grid_ahull[,1], xy_grid_ahull[,2])
  a_right <- sum(ii_inpoly_right) / nrow(S2grid) *4*pi
  
  a_bi <- sum(ii_inpoly_right & ii_inpoly_left) / nrow(S2grid) *4*pi
  
  tb_fov[f,] <- c(a_left, a_right, a_bi)
  
}


for (j in 1:length(fn)) {
  windows(width = 16, height = 8)
  print(plt_pt_ls[[j]])
  # print(plt_kde_ls[[j]])
}


# # SAVE
# save(N_lens, tb_fov, Zpred_lim, plt_pt_ls, plt_kde_ls,
#      file = "data/plt_ptkde_.RData")



# curvature and sphericity, 2024 ------------------------------------------------

for (f in 1:length(fn)) {
  load(paste("../microCT/2023_eyemap/", fn[f], ".RData", sep=''))  
  
  for (lr in c(TRUE, FALSE)) {
    if (lr) { #left
      pts <- lens[ind_left_lens,]
    } else { #right
      pts <- lens[!ind_left_lens,]
    }
    
    # hex coord
    ls <- reghex(pts, lefteye = lr)
    ind_nb <- ls[[1]]
    ind_xy <- ls[[2]]
    
    ind_roc <- matrix(ncol = 6, nrow = nrow(pts)) #[ind, p,v,q,h,sph]
    ind_dia <- matrix(ncol = 2, nrow = nrow(pts)) #[ind, d]
    ind_roc[,1] <- ind_xy[,1]
    ind_dia[,1] <- ind_xy[,1]
    
    # - avg dia, iterate over all pts, compute distnace to all 6 nbs, 
    for (j in 1:dim(ind_dia)[1]) {
      if (sum(is.na(ind_nb[j,])) == 0) {
        ind_dia[j,2] <- 
          sweep(pts[ind_nb[j,-1],], 2, pts[ind_nb[j,1],,drop=F])^2 %>%
          rowSums() %>%
          sqrt()%>%
          mean()
      }
    }
    
    # - iterate over p\v\q\h=const and compute radius of cur.
    N <- 2 #no. of nbs on one side, total = 2N+1
    
    # p=const
    for (cc in seq(min(ind_xy[,2]),max(ind_xy[,2]))) {
      ixy <- ind_xy[ind_xy[,2] == cc,,drop=F]
      if (nrow(ixy) >= 2*N+1) {
        ixy <- ixy[order(ixy[,3]),] #order by q
        # iterate over point on the line with sufficient nbs
        for (j in (1+N):(nrow(ixy)-N)) {
          xyz <- pts[ixy[j + (-N:N),1], ]
          ind_roc[match(ixy[j,1],ind_xy[,1]), 4] <- roc_fit(xyz) #along q
        }
      }
    }
    
    # q=const
    for (cc in seq(min(ind_xy[,3]),max(ind_xy[,3]))) {
      ixy <- ind_xy[ind_xy[,3] == cc,,drop=F]
      if (nrow(ixy) >= 2*N+1) {
        ixy <- ixy[order(ixy[,2]),] #order by p
        # iterate over point on the line with sufficient nbs
        for (j in (1+N):(nrow(ixy)-N)) {
          xyz <- pts[ixy[j + (-N:N),1], ]
          ind_roc[match(ixy[j,1],ind_xy[,1]), 2] <- roc_fit(xyz) #along p
        }
      }
    }
    
    # v=const, use N-1
    for (cc in seq(min(ind_xy[,2]+ind_xy[,3]), max(ind_xy[,2]+ind_xy[,3]))) {
      ixy <- ind_xy[ind_xy[,2]+ind_xy[,3] == cc,,drop=F]
      if (nrow(ixy) >= 2*(N-1)+1) {
        ixy <- ixy[order(ixy[,3]-ixy[,2]),] #order by h
        # iterate over point on the line with sufficient nbs
        for (j in (1+(N-1)):(nrow(ixy)-(N-1))) {
          xyz <- pts[ixy[j + (-(N-1):(N-1)),1], ]
          ind_roc[match(ixy[j,1],ind_xy[,1]), 5] <- roc_fit(xyz) #along h
        }
      }
    }
    
    # h=const
    for (cc in seq(min(ind_xy[,3]-ind_xy[,2]),max(ind_xy[,3]-ind_xy[,2]))) {
      ixy <- ind_xy[ind_xy[,3]-ind_xy[,2] == cc,,drop=F]
      if (nrow(ixy) >= 2*N+1) {
        ixy <- ixy[order(ixy[,2]+ixy[,3]),] #order by v
        # iterate over point on the line with sufficient nbs
        for (j in (1+N):(nrow(ixy)-N)) {
          xyz <- pts[ixy[j + (-N:N),1], ]
          ind_roc[match(ixy[j,1],ind_xy[,1]), 3] <- roc_fit(xyz) #along v
        }
      }
    }
    
    # sphere
    NS <- 2 #avg over no. of layers of hex
    for (j in 1:nrow(ind_xy)) {
      ii <- ind_xy[j,1]
      for (k in 1:NS) {
        ii <- unique(c(ind_nb[match(ii,ind_nb[,1]),]))
        ii <- na.omit(ii)
      }
      if (length(ii) >= (19-3)) { #missing one side
        xyz_data <- pts[ii, ]
        colnames(xyz_data) <- c('x','y','z')
        xyz_data2 <- xyz_data^2
        Y <- rowSums(xyz_data2)
        X <- cbind(2*xyz_data, rep(1,dim(xyz_data)[1]))
        sph <- lsfit(X,Y,intercept = FALSE)
        r_fit <- sqrt(sph[[1]][4]+sum(sph[[1]][1:3]^2)) #radius
        # cosph <- sph[["coefficients"]][c('x','y','z')] #center of the sphere
        # xyzr[j,] <- c(cosph, r_fit)
        ind_roc[j, 6] <- r_fit
      }
    }
    
    if (lr) { #left
      ind_roc_left <- ind_roc
      pts_left <- pts
      ind_xy_left <- ind_xy
      ind_nb_left <- ind_nb
      ind_dia_left <- ind_dia
    } else { #right
      ind_roc_right <- ind_roc
      pts_right <- pts
      ind_xy_right <- ind_xy
      ind_nb_right <- ind_nb
      ind_dia_right <- ind_dia
    }
  }
  # SAVE
  save(ind_roc_left, ind_roc_right, ind_xy_left, ind_xy_right,
       file = paste0("../microCT/2023_eyemap/", fn[f], "_roc.RData"))
  save(ind_dia_left, ind_dia_right, ind_nb_left, ind_nb_right,
       file = paste0("../microCT/2023_eyemap/", fn[f], "_dia.RData"))
}#loop

  

# some plots, cf. Figure_3_uCT.R for paper plots
# 3D
# getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
# pal_tar <- getPalette(n_lvl - 1)
nopen3d()
points3d(pts)
points3d(pts, size=10,
         col=plotrix::color.scale(
           ind_roc[,6],
           c(0,1),0,1,color.spec="rgb"
           # c(200,300),90,80,color.spec="hcl"
           )
         )
title3d(fn[f])


# hist
dev.new()
# pdf(paste("hist_e", letters[ind_e], "_T4", letters[LL], ".pdf", sep = ""))
hh <- hist(ind_roc[,4], breaks = seq(50, 500, by=10), na.rm=T, plot = F)
plot(hh$mids, hh$counts, type='l', bty='n', xlim = c(100, 400), xaxt='n',yaxt='n',
     xlab ="radius [um]", ylab='counts', main = "radius of curvature")
axis(1, at = seq(100, 400, by =50), labels = seq(100, 400, by =50) )
# dev.off()


# PLOT
nopen3d()
points3d(pts)
text3d(pts[!is.na(ind_roc_right[,2]),], texts = round(ind_roc_right[,2],0), adj = 1)
# text3d(pt[ind_r_h !=0,], texts = round(ind_r_h[ind_r_h !=0],0), adj = 1)

dev.new()
plot(1:nrow(pts), ind_roc_right[,2], ylim=c(0, 500))
points(ind_roc_right[,3], pch=17, col='orange')
points(ind_roc_right[,4], pch=16, col='forestgreen')

dev.new()
# pdf(paste("hist_e", letters[ind_e], "_T4", letters[LL], ".pdf", sep = ""))
rr <- ind_roc_right[,2]
rr[is.na(rr)] <- 0
hh <- hist(rr, breaks = seq(0, 600, by=10), plot = F)
plot(hh$mids, hh$density, type='l', bty='n', xlim = c(50, 500), xaxt='n',yaxt='n',
     xlab ="radius [um]", ylab='counts', main = "radius of curvature")
axis(1, at = seq(00, 500, by =50), labels = seq(00, 500, by =50) )
hh <- hist(ind_roc_right[,3], breaks = seq(0, 600, by=10), plot = F)
points(hh$mids, hh$density, type='l', col='orange')
hh <- hist(ind_roc_right[,4], breaks = seq(0, 600, by=10), plot = F)
points(hh$mids, hh$density, type='l', col='forestgreen')

legend(x = "topright", legend = c("sphere", "v-row", "h-row"), col = c("black","orange","forestgreen"), lwd = 2)    


# Below are miscellaneous stuff -------------------------------------------

# some plot, need specify data --------------------------------------------

# bar PLOT
# counts
df <- data.frame(count = rowSums(N_lens))
df$sex <- if_else(grepl("female", fn), "F", "M")
df$strain <- fn
df$strain <- factor(df$strain)
df$count <- as.numeric(df$count)
windows(record = F, width = 8, height = 8)
ggplot(df, aes(x=strain, y=count, fill=sex)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold", size=10, angle = 30),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text =  element_text(size = 14),
        plot.title = element_text(size=18) ) +
  # ylab("FOV %") +
  ggtitle("ommatidium count") 

# fov
tb <- tb_fov / 4/pi
df <- data.frame(FOV = tb[,1]+tb[,2]-tb[,3])
df$sex <- if_else(grepl("female", fn), "F", "M")
df$strain <- fn
df$strain <- factor(df$strain)
df$FOV <- as.numeric(df$FOV)
windows(record = F, width = 8, height = 8)
ggplot(df, aes(x=strain, y=FOV, fill=sex)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold", size=10, angle = 30),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text =  element_text(size = 14),
        plot.title = element_text(size=18)) +
  ylab("FOV") +
  ggtitle("% FOV") 

# binocular %4pi
tb <- tb_fov
df <- data.frame(count = tb[,3]/4/pi)
df$sex <- if_else(grepl("female", fn), "F", "M")
df$strain <- fn
df$strain <- factor(df$strain)
df$count <- as.numeric(df$count)
windows(record = F, width = 8, height = 8)
ggplot(df, aes(x=strain, y=count, fill=sex)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold", size=10, angle = 30),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text =  element_text(size = 14),
        plot.title = element_text(size=18) ) +
  ylab("overlap") +
  ggtitle("% binocular overlap")


# binocular fov
tb <- tb_fov
df <- data.frame(count = tb[,3] / (tb[,1]+tb[,2]-tb[,3]))
df$sex <- if_else(grepl("female", fn), "F", "M")
df$strain <- fn
df$strain <- factor(df$strain)
df$count <- as.numeric(df$count)
windows(record = F, width = 8, height = 8)
ggplot(df, aes(x=strain, y=count, fill=sex)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold", size=10, angle = 30),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text =  element_text(size = 14),
        plot.title = element_text(size=18) ) +
  ylab("overlap") +
  ggtitle("% binocular overlap of FOV")



# binocular angle near equator
azim_bi <- matrix(nrow = length(fn), ncol = 4)
for (f in 1:length(fn)) {
  load(paste("../microCT/data/", fn[f], ".RData", sep=''))
  
  ucl_rot_tp <- cart2sphZ(ucl_rot)[,2:3] %>%
    as_tibble() %>%
    mutate(theta = theta/pi*180, phi = phi/pi*180) %>%
    mutate(phi = if_else(phi > 180, phi - 360, phi)) %>%
    mutate(phi = - phi) %>% #inside out
    as.matrix()
  
  tp <- ucl_rot_tp[ind_left_cone, ]
  azim_bi[f,1:2] <- range(tp[tp[,1] < 110 & tp[,1] > 70, 2])
  
  tp <- ucl_rot_tp[!ind_left_cone, ]
  azim_bi[f,3:4] <- range(tp[tp[,1] < 110 & tp[,1] > 70, 2])
}

tb <- azim_bi[, 2:3]
df <- data.frame(count = (tb[,1]-tb[,2])/2)
df$sex <- if_else(grepl("female", fn), "F", "M")
df$strain <- fn
df$strain <- factor(df$strain)
df$count <- as.numeric(df$count)
windows(record = F, width = 8, height = 8)
ggplot(df, aes(x=strain, y=count, fill=sex)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(face="bold", size=10, angle = 30),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text =  element_text(size = 14),
        plot.title = element_text(size=18) ) +
  ylab("azimuth [degree]") +
  ggtitle("binocular overlap near equator")

# # lens area ---------------------------------------------------------------
# 
# library(alphashape3d)
# 
# tb_area <- matrix(ncol = 2, nrow = 9)
# 
# for (f in 1:9) {
#   load(paste("data/proc_", fn[f], ".RData", sep='')) 
#   
#   pts <- lens[ind_left_lens, ]
#   
#   # lens_msh <- ashape3d(as.matrix(pts), alpha = 100) %>% as.mesh3d()
#   
#   vor <- deldir::deldir(df$xM, df$yM)
#   
# }


# # -- lens area 
# 
# nopen3d()
# points3d(pts)
# linemat <- matrix(t(cbind(pts[dtdt$delsgs$ind1,],pts[dtdt$delsgs$ind2,])),
#                   ncol = 3, byrow = T)
# segments3d(linemat, color = "grey")
# 
# 
# # df <- cbind(lens[ind_left_cone,1:2], ioa) %>% as.data.frame()
# 
# 
# 
# ggplot() +
#   geom_point(data = stereo_df[101:105,], aes(x,y), size = 5, colour='red') +
#   geom_segment( aes(x = x1, y = y1, xend = x2, yend = y2),
#                 size = 1,
#                 data = dtdt$dirsgs,
#                 linetype = 1,
#                 color= "black") 
# 
#   # geom_polygon( aes(x = x1, y = y1, xend = x2, yend = y2),
#   #               size = 2,
#   #               data = dtdt$dirsgs,
#   #               linetype = 1,
#   #               color= "grey") +
#   # geom_point(aes(colour = layer), size = 3) +
#   # scale_colour_manual(values = pal_tar, guide_legend("layer"),
#   #                     labels = c('empty', '1', '2','3','4')) +
#   # # scale_color_brewer(palette = "RdYlBu", direction = 1) +
#   # geom_path(data = poly_xy, aes(x,y), lwd = 1) +
#   # xlim(-1, pi) +
#   # ylim(-pi/2, pi/2) +
#   # scale_x_continuous(expand = c(0, 0)) +
#   # scale_y_reverse(expand = c(0, 0)) +
#   # theme_bw() +
#   # labs(title = paste(j, anno_neu[j,'name'], sep = ',')) +
#   # # labs(title = paste(round(LPTC_nodenum[j,k], digits = 2), LPTC_area[j,k], sep = ',')) +
#   # coord_fixed(ratio = 1)
# 
# 
# # PLOT
# 
# df <- cbind(ucl_rot_Mo[ind_left_cone,1:2], dtdt$summary$del.area) %>% as.data.frame()
# windows(record = F, width = 16, height = 8)
# ggplot(df, aes(x = xM, y = yM, z=V3)) +
#   geom_contour_filled(bins = 5)
# # ggplot(df, aes(x = xM, y = yM)) +
# geom_point(aes(colour = factor(left))) +
#   scale_colour_manual(values = c("magenta", "green"), labels=c('left','right')) +
#   # geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), 
#   #              arrow = arrow(length = unit(0.005, "npc"),type = "closed"), size =1) +
#   # geom_path(data = poly_xy, aes(x,y), lwd = 1) +
#   # annotate("text", x = 2.3, y = -1.1, label = "trans", size = 5, stroke = 3) +
#   geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
#   theme_void() +
#   theme(legend.position="none", panel.background = element_blank()) +
#   coord_fixed(ratio = 1, xlim = c(-sqrt(8), sqrt(8)), ylim = c(-sqrt(2), sqrt(2)))


linemat <- matrix(t(cbind(pts[delsgs_mod$ind1,],pts[delsgs_mod$ind2,])),
                  ncol = 3, byrow = T)
segments3d(linemat, color = "grey")


# inter-ommatidium angle, hex -------------------------------------------



# inter-ommatidium angle, Delaunay ----------------------------------------------

# Delaunay triangulation
plt_cont_ls <- list()
plt_cont_iod_ls <- list()
plt_outline_ls <- list()

for (f in 1:length(fn)) {
  load(paste("../microCT/data/eq_corr/", fn[f], ".RData", sep=''))  
  
  # - Mollweide
  ucl_rot_tp <- cart2sphZ(ucl_rot)[,2:3] %>%
    as_tibble() %>%
    mutate(theta = theta/pi*180, phi = phi/pi*180) %>%
    mutate(phi = if_else(phi > 180, phi - 360, phi)) %>%
    mutate(phi = - phi) %>% #inside out
    as.matrix()
  ucl_rot_Mo <- Mollweide(ucl_rot_tp)
  ucl_rot_Mo <- cbind(ucl_rot_Mo, 1)
  ucl_rot_Mo[!ind_left_cone, 3] <- 2 # left == 1, right == 2
  colnames(ucl_rot_Mo) <- c('xM', 'yM', 'left')
  
  # df_left <- ucl_rot_Mo %>% as_tibble() %>% filter(left == 1) %>% as.data.frame()
  # df_right <- ucl_rot_Mo %>% as_tibble() %>% filter(left == 2) %>% as.data.frame()
  
  # inter-ommatidium angle in Mollweide
  ioa_Mo <- list()
  iod_Mo <- list()
  for (lr in 1:2) {
    if (lr == 1) {
      pts <- ucl_rot[ind_left_cone, ]
    } else {
      pts <- ucl_rot[!ind_left_cone, ]
    }
    
    # rotate to southern hemisphere
    vcom <- colMeans(pts)
    vcom <- vcom / sqrt(sum(vcom^2))
    R_mat <- quaternion3D( cross3D(vcom, c(0,0,-1)), acos(vcom %*% c(0,0,-1)) /pi*180 )
    pts_z <- as.matrix(pts) %*% t(R_mat)
    
    # stereographic projection
    x <- pts_z[,1] / (1 - pts_z[,3])
    y <- pts_z[,2] / (1 - pts_z[,3])
    stereo <- cbind(x, y) 
    stereo_df <- as.data.frame(stereo)
    
    # # - geometry
    # library(geometry)
    # 
    # td <- delaunayn(stereo)
    # 
    # dev.new()
    # trimesh(td, stereo)
    # points(stereo[101:105,], col='red', cex = 2, pch = 16)
    
    
    # - deldir
    dtdt <- deldir(stereo[,1], stereo[,2])
    
    # -- deal with edge points, 
    # if the boundary pt, then check the 2 circumcenters sharing this edge,
    # remove the one with no corresponding circumcenter
    # Alt,  deal with edge points wutg dummy pt. NOT tried ?
    
    ee_mat <- matrix(c(1,1,1,2,2,1,2,2), ncol = 2, byrow = T) #index matrix
    
    bp <- which(dtdt$dirsgs$bp1 | dtdt$dirsgs$bp2)
    
    delsgs_mod <- dtdt$delsgs # modify
    for (j in bp) {
      ind_edge <- dtdt$dirsgs[j, c('ind1', 'ind2')] #edge vertices 
      ind1_tri <- c(delsgs_mod[delsgs_mod$ind1 %in% ind_edge[1], ]$ind2,
                    delsgs_mod[delsgs_mod$ind2 %in% ind_edge[1], ]$ind1)
      ind2_tri <- c(delsgs_mod[delsgs_mod$ind1 %in% ind_edge[2], ]$ind2,
                    delsgs_mod[delsgs_mod$ind2 %in% ind_edge[2], ]$ind1)
      ind_tri <- c(ind1_tri, ind2_tri)
      ind_share <- ind_tri[duplicated(ind_tri)] # two vertices shared by the edge vertices
      
      if (length(ind_share) == 2) {
        ee <- c() #exist a circumcenter
        for (m in 1:2) {
          for (n in 1:2) {
            b1 <- ind_edge[m] %in% dtdt$dirsgs[dtdt$dirsgs$ind1 %in% ind_share[n], ]$ind2
            b2 <- ind_edge[m] %in% dtdt$dirsgs[dtdt$dirsgs$ind2 %in% ind_share[n], ]$ind1
            ee <- c(ee, b1 | b2)
          }
        }
        # ind_edge[1] %in% dtdt$dirsgs[dtdt$dirsgs$ind1 %in% ind_share[2], ]$ind2
        # ind_edge[1] %in% dtdt$dirsgs[dtdt$dirsgs$ind2 %in% ind_share[2], ]$ind1
        # 
        # ind_edge[2] %in% dtdt$dirsgs[dtdt$dirsgs$ind1 %in% ind_share[2], ]$ind2
        # ind_edge[2] %in% dtdt$dirsgs[dtdt$dirsgs$ind2 %in% ind_share[2], ]$ind1
        
        if (sum(ee) < 4) {
          for (k in which(!ee)) {
            # to remove
            ind1 <- delsgs_mod$ind1 %in% ind_edge[ee_mat[k,][1]] & delsgs_mod$ind2 %in% ind_share[ee_mat[k,][2]] 
            ind2 <- delsgs_mod$ind2 %in% ind_edge[ee_mat[k,][1]] & delsgs_mod$ind1 %in% ind_share[ee_mat[k,][2]] 
            ind_rm <- ind1 | ind2
            delsgs_mod <- delsgs_mod[-which(ind_rm), ]
          }
        }
      }
    }
    
    # #  PLOT
    # dev.new()
    # ggplot() +
    #   geom_point(data = stereo_df[c(7,18),], aes(x,y), size = 5, colour='cyan') +
    #   geom_point(data = stereo_df[c(16),], aes(x,y), size = 5, colour='blue') +
    #   # geom_segment(data = dtdt$delsgs, aes(x = x1, y = y1, xend = x2, yend = y2),
    #   #              size = 1, linetype = 1, color= "black") +
    #   geom_segment(data = delsgs_mod, aes(x = x1, y = y1, xend = x2, yend = y2),
    #                size = 1, linetype = 1, color= "black") +
    #   geom_segment(data = dtdt$dirsgs, aes(x = x1, y = y1, xend = x2, yend = y2),
    #                size = 1, linetype = 1, color= "red") +
    #   theme_bw()
    
    # -- inter-omma angle
    ioa <- c()
    iod <- c() # distances between lenses
    for (j in 1:nrow(pts)) {
      dd <- c()
      ddlens <- c()
      ii <- which(delsgs_mod$ind1 %in% j)
      if (length(ii) > 0) {
        iinb <- delsgs_mod$ind2[ii]
        for (k in iinb) {
          dd <- c(dd, arcLength(pts[j,], pts[k,]))
          ddlens <- c(ddlens, sqrt(sum((lens[j,] - lens[k,])^2)) )
        }
      }
      ii <- which(delsgs_mod$ind2 %in% j)
      if (length(ii) > 0) {
        iinb <- delsgs_mod$ind1[ii]
        for (k in iinb) {
          dd <- c(dd, arcLength(pts[j,], pts[k,]))
          ddlens <- c(ddlens, sqrt(sum((lens[j,] - lens[k,])^2)) )
        }
      }
      ioa <- c(ioa, mean(dd)/pi*180)
      iod <- c(iod, mean(ddlens))
    }
    
    if (lr == 1) {
      ioa_Mo[[lr]] <- cbind(ucl_rot_Mo[ind_left_cone,1:2], ioa) %>% as.data.frame()
      iod_Mo[[lr]] <- cbind(ucl_rot_Mo[ind_left_cone,1:2], iod) %>% as.data.frame()
    } else {
      ioa_Mo[[lr]] <- cbind(ucl_rot_Mo[!ind_left_cone,1:2], ioa) %>% as.data.frame()
      iod_Mo[[lr]] <- cbind(ucl_rot_Mo[!ind_left_cone,1:2], iod) %>% as.data.frame()
    }
  }
  
  # # voronoi
  # dev.new()
  # ggplot() +
  #   geom_segment(data = dtdt$dirsgs, aes(x = x1, y = y1, xend = x2, yend = y2),
  #                 size = 1, linetype = 1, color= "black") +
  #   theme_bw()
  # 
  # # PLOT
  # nopen3d()
  # points3d(pts)
  # linemat <- matrix(t(cbind(pts[dtdt$delsgs$ind1,],pts[dtdt$delsgs$ind2,])),
  #                   ncol = 3, byrow = T)
  # segments3d(linemat, color = "grey")
  
  
  # # PLOT
  # df <- ioa_Mo[[2]]
  # df$equalSpace <- cut(df$ioa, seq(1,11,by = 2))
  # 
  # windows(record = F, width = 16, height = 8)
  # ggplot(df, aes(x = xM, y = yM)) +
  #   # scale_colour_viridis_b() 
  #   # geom_contour()
  #   geom_tile(aes(fill = equalSpace), width=0.2) +
  #   # geom_point(aes(colour=ioa)) 
  
  
  # -- pal and breaks
  n_lvl <- 8
  breaks_ioa <- seq(1,11,length.out = n_lvl+1)
  getPalette <- colorRampPalette(brewer.pal(9, "Purples"))
  pal_left <- getPalette(n_lvl) %>% rev()
  getPalette <- colorRampPalette(brewer.pal(9, "Greens"))
  pal_right <- getPalette(n_lvl) %>% rev()

  # -- add two color scales
  pal_sum <- rep(0,n_lvl*(n_lvl+2))
  for (j in seq((n_lvl+1), n_lvl*(n_lvl+1), length.out = n_lvl)) {
    pal_sum[j] <- pal_left[j/(n_lvl+1)]
    c1 <- col2rgb(pal_left[j/(n_lvl+1)], alpha = F)
    for (k in 1:n_lvl) {
      pal_sum[k] <- pal_right[k]
      c2 <- col2rgb(pal_right[k], alpha = F)
      c3 <- (c1 + c2)/2
      pal_sum[j+k] <- rgb(t(c3), maxColorValue = 255)
    }
  }
  
  # -- grid for loessm ioa
  grid_M <- expand.grid(xM = seq(-sqrt(8), sqrt(8), length.out = 100),
                        yM = seq(-sqrt(2), sqrt(2), length.out = 50) )
  ii_inpoly <- sp::point.in.polygon(grid_M[,1], grid_M[,2], bkgd_chull[,1], bkgd_chull[,2])
  
  # -- loess
  fit_loess <- loess(ioa ~ xM * yM, data = ioa_Mo[[1]], degree = 2, span = 0.1,
                     control = loess.control(surface = "direct"))
  pred_loess <- predict(fit_loess, grid_M, se = T)
  df_pred <- grid_M
  df_pred$Z <- melt(pred_loess$fit)$value
  df_pred$equalSpace <- cut(df_pred$Z, breaks_ioa)
  df_pred_left <- df_pred[ii_inpoly == 1, ]
  df_pred_left$colSpace <- as.numeric(df_pred_left$equalSpace)
  df_pred_left$colSpace[is.na(df_pred_left$colSpace)] <- 0
  # df_pred_left$colSpace[df_pred_left$colSpace == 1] <- 0
  # df_pred_left$colSpace <- df_pred_left$colSpace * (n_lvl+1)
  
  xy_ashape <- ashape(ucl_rot_Mo[ind_left_cone, 1:2] +
                        matrix(runif(sum(ind_left_cone)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
  xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
  ahull_left <- xy_grid_ahull
  ii_inpoly_left <- sp::point.in.polygon(grid_M[ii_inpoly == 1,1], grid_M[ii_inpoly == 1,2], xy_grid_ahull[,1], xy_grid_ahull[,2])
  df_pred_left$Z[!ii_inpoly_left] <- 20 # outside set to large angles
  df_pred_left$colSpace[!ii_inpoly_left] <- 0
  
  # right
  fit_loess <- loess(ioa ~ xM * yM, data = ioa_Mo[[2]], degree = 2, span = 0.1,
                     control = loess.control(surface = "direct"))
  pred_loess <- predict(fit_loess, grid_M, se = T)
  df_pred <- grid_M
  df_pred$Z <- melt(pred_loess$fit)$value
  df_pred$equalSpace <- cut(df_pred$Z, breaks_ioa)
  df_pred_right <- df_pred[ii_inpoly == 1, ]
  df_pred_right$colSpace <- as.numeric(df_pred_right$equalSpace)
  df_pred_right$colSpace[is.na(df_pred_right$colSpace)] <- 0
  
  xy_ashape <- ashape(ucl_rot_Mo[!ind_left_cone, 1:2] +
                        matrix(runif(sum(!ind_left_cone)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
  xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
  ahull_right <- xy_grid_ahull
  ii_inpoly_right <- sp::point.in.polygon(grid_M[ii_inpoly == 1,1], grid_M[ii_inpoly == 1,2], xy_grid_ahull[,1], xy_grid_ahull[,2])
  df_pred_right$Z[!ii_inpoly_right] <- 20
  df_pred_right$colSpace[!ii_inpoly_right] <- 0
  
  # df_pred <- df_pred_left[, c('xM','yM')]
  # df_pred$Z <- df_pred_left$Z + df_pred_right$Z
  # df_pred$colSpace <- df_pred_left$colSpace + df_pred_right$colSpace
  # df_pred$colSpace[df_pred$colSpace == 0] <- NA
  
  # df_pred_left$Z[!ii_inpoly_left] <- NA
  # df_pred_right$Z[!ii_inpoly_right] <- NA
  
  
  
  breaks_contour <- c(0, 2.5, 3, 4, 5)
  getPalette <- colorRampPalette(c('gray30', 'gray90'))
  pal_contour <- getPalette(length(breaks_contour)+1)
  
  # # PLOT contour
  # windows(record = F, width = 16, height = 8)
  
  # plt <- ggplot() +
  # ggplot() +
  # --- raster
  # geom_raster(data=df_pred, aes(xM, yM, fill = factor(colSpace)), alpha=0.2, interpolate = T) +
  # scale_fill_manual(values = pal_sum, breaks = seq(1, n_lvl*(n_lvl+2), by = 1)) +
  # --- filled contour
  # geom_contour_filled(data=df_pred_left, aes(x=xM,y=yM, z=Z), breaks=breaks_contour,
  #                     color=pal_left[3], alpha=0.5, lwd=2) +
  # scale_fill_manual(values = pal_contour, breaks=breaks_contour, labels = breaks_contour) +
  # geom_point(data=as.data.frame(ucl_rot_Mo),aes(x=xM, y=yM)) +
  # geom_contour_filled(data=df_pred_right, aes(x=xM,y=yM, z=Z), breaks=breaks_contour,
  #                     color=pal_right[3], alpha=0.5, lwd=2) +
  # scale_fill_manual(values = pal_contour, breaks=breaks_contour, labels=breaks_contour) +
  # --- contour 
  plt <- ggplot() +
    geom_contour(data=df_pred_left, aes(x=xM,y=yM, z=Z), breaks=breaks_contour,
                 colour=pal_left[1], alpha=0.5, lwd=2) +
    # stat_contour(data=df_pred_right, aes(x=xM,y=yM, z=Z, colour=..level..), colour=pal_right[1],
    #              breaks=breaks_contour, alpha=0.5, lwd=2) +
    geom_contour(data=df_pred_right, aes(x=xM,y=yM, z=Z), breaks=breaks_contour,
                 colour=pal_right[1], alpha=0.5, lwd=2) +
    geom_point(data = as.data.frame(ucl_rot_Mo), aes(xM, yM, colour = factor(left))) +
    scale_colour_manual(name = "side", values = c("magenta", "green"), labels=c('left','right'), guide = F) +
    geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
    theme_void() +
    # theme_bw() +
    theme(legend.position="none", panel.background = element_blank()) +
    labs(title = fn[f]) +
    coord_fixed(ratio = 1, xlim = c(-sqrt(8), sqrt(8)), ylim = c(-sqrt(2), sqrt(2)))
  
  # plt <- direct.label(plt, list("bottom.pieces", cex = 2))
  
  # USE # library(ggnewscale)
  # ggplot() +
  #   geom_raster(data=df_pred_left, aes(xM, yM, fill = Z),alpha=1, interpolate = T) +
  #   scale_fill_gradientn(colours = pal_left, limits = c(1,8), na.value = NA) +
  #   new_scale('fill') +
  #   geom_raster(data=df_pred_right, aes(xM, yM, fill = Z),alpha=0.5, interpolate = T) +
  #   scale_fill_gradientn(colours = pal_right, limits = c(1,8), na.value = NA) +
  
  plt_cont_ls[[f]] <- plt
  
  # file_name = paste(fn[f], "_contour.png", sep="")
  # png(file_name, width = 1600, height = 800, pointsize = 12)
  # print(plt)
  # dev.off()
  
  
  # -- grid for loessm iod
  grid_M <- expand.grid(xM = seq(-sqrt(8), sqrt(8), length.out = 100),
                        yM = seq(-sqrt(2), sqrt(2), length.out = 50) )
  ii_inpoly <- sp::point.in.polygon(grid_M[,1], grid_M[,2], bkgd_chull[,1], bkgd_chull[,2])
  
  # -- loess
  fit_loess <- loess(iod ~ xM * yM, data = iod_Mo[[1]], degree = 2, span = 0.1,
                     control = loess.control(surface = "direct"))
  pred_loess <- predict(fit_loess, grid_M, se = T)
  df_pred <- grid_M
  df_pred$Z <- melt(pred_loess$fit)$value
  df_pred$equalSpace <- cut(df_pred$Z, breaks_ioa)
  df_pred_left <- df_pred[ii_inpoly == 1, ]
  df_pred_left$colSpace <- as.numeric(df_pred_left$equalSpace)
  df_pred_left$colSpace[is.na(df_pred_left$colSpace)] <- 0
  
  xy_ashape <- ashape(ucl_rot_Mo[ind_left_cone, 1:2] +
                        matrix(runif(sum(ind_left_cone)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
  xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
  ahull_left <- xy_grid_ahull
  ii_inpoly_left <- sp::point.in.polygon(grid_M[ii_inpoly == 1,1], grid_M[ii_inpoly == 1,2], xy_grid_ahull[,1], xy_grid_ahull[,2])
  df_pred_left$Z[!ii_inpoly_left] <- 20 # outside set to large angles
  df_pred_left$colSpace[!ii_inpoly_left] <- 0
  
  # right
  fit_loess <- loess(iod ~ xM * yM, data = iod_Mo[[2]], degree = 2, span = 0.1,
                     control = loess.control(surface = "direct"))
  pred_loess <- predict(fit_loess, grid_M, se = T)
  df_pred <- grid_M
  df_pred$Z <- melt(pred_loess$fit)$value
  df_pred$equalSpace <- cut(df_pred$Z, breaks_ioa)
  df_pred_right <- df_pred[ii_inpoly == 1, ]
  df_pred_right$colSpace <- as.numeric(df_pred_right$equalSpace)
  df_pred_right$colSpace[is.na(df_pred_right$colSpace)] <- 0
  
  xy_ashape <- ashape(ucl_rot_Mo[!ind_left_cone, 1:2] +
                        matrix(runif(sum(!ind_left_cone)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
  xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
  ahull_right <- xy_grid_ahull
  ii_inpoly_right <- sp::point.in.polygon(grid_M[ii_inpoly == 1,1], grid_M[ii_inpoly == 1,2], xy_grid_ahull[,1], xy_grid_ahull[,2])
  df_pred_right$Z[!ii_inpoly_right] <- 20
  df_pred_right$colSpace[!ii_inpoly_right] <- 0

  breaks_contour <- seq(from = range(df_pred_right$Z)[1], to = range(df_pred_right$Z)[2], length.out = 5)
  getPalette <- colorRampPalette(c('gray30', 'gray90'))
  pal_contour <- getPalette(length(breaks_contour)+1)
  
  # --- contour 
  plt <- ggplot() +
    geom_contour(data=df_pred_left, aes(x=xM,y=yM, z=Z), breaks=breaks_contour,
                 colour=pal_left[1], alpha=0.5, lwd=2) +
    geom_contour(data=df_pred_right, aes(x=xM,y=yM, z=Z), breaks=breaks_contour,
                 colour=pal_right[1], alpha=0.5, lwd=2) +
    geom_point(data = as.data.frame(ucl_rot_Mo), aes(xM, yM, colour = factor(left))) +
    scale_colour_manual(name = "side", values = c("magenta", "green"), labels=c('left','right'), guide = F) +
    geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
    theme_void() +
    theme(legend.position="none", panel.background = element_blank()) +
    labs(title = paste("iod", fn[f])) +
    coord_fixed(ratio = 1, xlim = c(-sqrt(8), sqrt(8)), ylim = c(-sqrt(2), sqrt(2)))

  plt_cont_iod_ls[[f]] <- plt
  
  # - plot outlines
  # windows(record = F, width = 16, height = 8)
  plt <- ggplot() +
    geom_point(data = as.data.frame(ucl_rot_Mo), aes(xM, yM, colour = factor(left))) +
    scale_colour_manual(name = "side", values = c("magenta", "green"), labels=c('left','right'), guide = F) +
    geom_path(data = ahull_left, aes(x1,y1), lwd = 1.5, colour = "magenta", alpha =0.7) +
    geom_path(data = ahull_right, aes(x1,y1), lwd = 1.5, colour = "green", alpha =0.7) +
    geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
    theme_void() +
    # theme_bw() +
    theme(legend.position="none", panel.background = element_blank()) +
    labs(title = fn[f]) +
    coord_fixed(ratio = 1, xlim = c(-sqrt(8), sqrt(8)), ylim = c(-sqrt(2), sqrt(2)))
  
  plt_outline_ls[[f]] <- plt
  
  # file_name = paste(fn[f], "_outline.png", sep="")
  # png(file_name, width = 1600, height = 800, pointsize = 12)
  # print(plt)
  # dev.off()
}

for (j in 1:length(fn)) {
  windows(width = 16, height = 8)
  # print(plt_cont_ls[[j]])
  print(plt_cont_iod_ls[[j]])
  # print(plt_outline_ls[[j]])
}

# SAVE
save(plt_cont_ls, plt_cont_iod_ls, plt_outline_ls, file = "data/plt_contour_outline_ls.RData")




# match to ref data set ---------------------------------------------------

load(paste("../microCT/data/", fn[1], ".RData", sep=''))
ucl_rot_ref <- ucl_rot
ind_left_cone_ref <- ind_left_cone

ind_match_ls <- list()

for (f in 2:length(fn)) {
  load(paste("../microCT/data/", fn[f], ".RData", sep=''))
  
  pt1 <- ucl_rot_ref[!ind_left_cone_ref,]
  pt2 <- sweep(ucl_rot[!ind_left_cone,], 2, c(-0.5,0.5,0))
  
  # Hungarian
  L1 <- nrow(pt1)
  L2 <- nrow(pt2)
  
  Lmax <- max(L1, L2)
  Lmin <- min(L1, L2)
  mDist <- matrix(ncol = Lmax, nrow = Lmin)
  if (Lmin == L1) {
    dr <- pt1
    dc <- pt2
  } else {
    dr <- pt2
    dc <- pt1
  }
  for (j in 1:Lmin) {
    mDist[j,] <- as.matrix(dist(rbind(dr[j,], dc)))[-1,1]
  }
  i_match <- solve_LSAP(mDist)
  
  if (Lmin == L1) {
    ind_match <- cbind(seq(1,length(i_match)), i_match)
    # ind_match[!(ind_match[,2] %in% seq(1, L2)), 2] <- NA # def unmatched as NA
  } else {
    ii <- cbind(i_match, seq(1,length(i_match)) )
    ind_match <- ii[order(ii[,1]),]
    # ind_match[!(ind_match[,1] %in% seq(1, L1)), 1] <- NA
  }
  
  ind_match_ls[[f]] <- ind_match
  
  nopen3d()
  points3d(pt1, size = 5, col = 'gray30')
  points3d(pt2, size = 5, col = 'magenta')
  linemat <- matrix(t(cbind(pt1[ind_match[,1],], pt2[ind_match[,2],] )), ncol = 3, byrow = T)
  segments3d(linemat, color = "grey50")
  title3d(fn[f])
}

# save(fn, ind_match_ls, file = "data/uCT_match.RData")

# fov on a sphere of radius r ---------------------------------------------
# https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

r_env <- 500*10^3  # [um]

vleft <- cross3D(vn, vf)
mrot <- cbind(vf,vleft,vn)
lens_rot <- sweep(as.matrix(lens),2,eye_pca$center)  %*% mrot 
cone_rot <- sweep(as.matrix(cone),2,eye_pca$center) %*% mrot 

# # PLOT
# nopen3d()
# points3d(lens_rot[ind_left_lens,], size = 5, color = 'magenta')
# points3d(lens_rot[!ind_left_lens,], size = 5, color = 'green')
# points3d(cone_rot, size = 5, color = 'blue')
# arrow3d(c(0,0,0), c(0,0,1)*100, theta = pi/6,n = 4, col="blue", type = "rotation")
# planes3d(0,0,1,0, alpha = 0.5)
# arrow3d(c(0,0,0), c(1,0,0)*100, theta = pi/6,n = 4, col="brown", type = "rotation")
# linemat <- matrix(t(cbind(cone_rot, lens_rot[i_match,])), ncol = 3, byrow = T)
# segments3d(linemat, color = "grey")
# spheres3d(0,0,0,r_env, col='grey', alpha=0.1, lit=F)


omc <- sweep(cone_rot, 2, c(0,0,0)) # o - c term
d <- - rowSums(ucl_rot * omc) + sqrt(rowSums(ucl_rot * omc)^2 - (rowSums(omc^2) - r_env^2))
p <- cone_rot + sweep(ucl_rot, 1, d, FUN = '*')

nopen3d()
# points3d(lens_rot[ind_left_lens,], size = 5, color = 'magenta')
# points3d(lens_rot[!ind_left_lens,], size = 5, color = 'green')
points3d(cone_rot, size = 5, color = 'blue')
# linemat <- matrix(t(cbind(cone_rot, lens_rot[i_match,])), ncol = 3, byrow = T)
# segments3d(linemat, color = "grey")
points3d(p[ind_left_cone,], size = 6, color = 'magenta')
points3d(p[!ind_left_cone,], size = 6, color = 'green')
spheres3d(0,0,0,r_env, col='grey', alpha=0.1, lit=F)
planes3d(0,1,0, 0, alpha = 0.1)
planes3d(0,0,1, 0, alpha = 0.1)




# calculate the longitude FOV
p_xyz <- sweep(p[ii_left,], 2, cosph) #left or right
# p_xyz <- sweep(p, 2, cosph) #left or right
# p_xyz <- p_xyz %*% matrix(c(0,-1,0,1,0,0,0,0,1), ncol = 3) # front = y-axis
# colnames(p_xyz) <- c('x','y','z')
p_utp <- sweep(p_xyz, 1, sqrt(rowSums(p_xyz^2)), '/')
p_utp %<>% as_tibble() %>%
  mutate(theta = acos(z)) %>%
  mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
  mutate(phi = phi - 2*pi*(phi > pi) ) %>%
  as.data.frame()

range(p_utp[p_utp$theta < (pi/2 + 10 / 180 * pi) & p_utp$theta > (pi/2 - 10 / 180 * pi), ]$phi) / pi * 180

dev.new()
plot(p_utp[,c('phi', 'theta')])


# confocal data --------------------------------------------------------------
# ucl are very noisy. Data segmentation problem ?

# # Aaron's
# ends_data <- read.csv("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_eyemap/Lens_Medulla/data/ends_Detailed.csv", header=FALSE)
# colnames(ends_data) <- c("x", "y", "z")
# rownames(ends_data) <- seq(1,dim(ends_data)[1])
# tips_data <- read.csv("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_eyemap/Lens_Medulla/data/tips_Detailed.csv", header=FALSE)
# colnames(tips_data) <- c("x", "y", "z")
# rownames(tips_data) <- seq(1,dim(tips_data)[1])
# # medcol_data <- read.csv("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_eyemap/Lens_Medulla/data/eye nc82 rfp gfp 25x 2nd CT1 BUNDLES_Position.csv", header = TRUE)
# # medcol_data <- medcol_data[,1:3]
# # colnames(medcol_data) <- c("x", "y", "z")

lens_data <- read.csv("data/Eye nc82 RFP GFP bisect 25x 8th CONE_tips, ends, equator_lens_Detailed.csv", header=T)
lens_data <- lens_data[,c('Position.X','Position.Y','Position.Z')]
colnames(lens_data) <- c("x", "y", "z")

tip_data <- read.csv("data/Eye nc82 RFP GFP bisect 25x 8th CONE_tips, ends, equator_PRtips_Detailed.csv", header=T)
tip_data <- tip_data[,c('Position.X','Position.Y','Position.Z',"Default.Labels")]
colnames(tip_data) <- c("x", "y", "z", "label")
ind_Down <- which(tip_data$label == "Down")
ind_Down <- ind_Down[-match(c(747, 746, 745, 744),ind_Down)]
ind_Up <- which(tip_data$label == "Up")
tip_data <- tip_data[,c("x", "y", "z")]

lens <- lens_data
cone <- tip_data

# - alignment
# equator
v_up <- colMeans(cone[ind_Up,]) - colMeans(cone[ind_Down,])
cone_equator <- cone[c(ind_Up, ind_Down),]
pt_theta90 <- cone_equator
lm_plane <- lm(pt_theta90$x ~ pt_theta90$y + pt_theta90$z)
vn <- c(-1, lm_plane$coefficients[2], lm_plane$coefficients[3]) #norm
# vn <- prcomp(pt_theta90)$rotation[,3]
vn <- vn/sqrt(sum(vn^2))
# which side is up
if (sum(vn * v_up) < 0) {
  vn <- -vn
}
# pt_center <- colMeans(cone)

# front via pca, use both eyes
eye_pca <- prcomp(cone)
vf <- eye_pca$rotation[,3]
vf <- vf - c(vf %*% vn) * vn
vf <- vf / sqrt(sum(vf^2))
if (vf %*% (colMeans(lens) - eye_pca$center) < 0) {
  vf <- - vf
}
vf <- vf %*% t(quaternion3D(vn, 45))

vleft <- cross3D(vn, vf)

mrot <- cbind(t(vf),vleft,as.matrix(vn))
lens <- sweep(as.matrix(lens),2,eye_pca$center,'-') %*% mrot
cone <- sweep(as.matrix(cone),2,eye_pca$center,'-') %*% mrot


# - matching, Hungarian
# pt1 <- tips_data
# pt2 <- ends_data
pt1 <- lens
pt2 <- cone

L1 <- nrow(pt1)
L2 <- nrow(pt2)

Lmax <- max(L1, L2)
Lmin <- min(L1, L2)
mDist <- matrix(ncol = Lmax, nrow = Lmin)
if (Lmin == L1) {
  dr <- pt1
  dc <- pt2
} else {
  dr <- pt2
  dc <- pt1
}
for (j in 1:Lmin) {
  mDist[j,] <- as.matrix(dist(rbind(dr[j,], dc)))[-1,1]
}
i_match <- solve_LSAP(mDist)

nopen3d()
points3d(pt1, size = 5, col = 'gray30')
points3d(pt2, size = 5, col = 'magenta')
linemat <- matrix(t(cbind(dr, dc[i_match,] )), ncol = 3, byrow = T)
segments3d(linemat, color = "grey50")

# ucl
LM_ucl <- lens - cone[i_match,]
LM_ucl <- sweep(LM_ucl, 1, sqrt(rowSums(LM_ucl^2)), '/')

nopen3d()
points3d(LM_ucl, size=10)
spheres3d(0,0,0,0.99, alpha=0.1, lit=F)

# acos(as.numeric(LM_ucl[616,]) %*% as.numeric(LM_ucl[475,])) /pi*180 # 157 deg

# save(lens, cone, ind_Up, ind_Down, vn, vf, vleft, LM_ucl, i_match,
#      file = paste("../LM/aaron", ".RData", sep='') )

# more LM data, round 4 ------------------------------------------------------------

lens <- read.csv("../LM/CT1_round_4/fem 25x 1st weak GFP-Lenscutout-w-spots11.5um_Detailed.csv", header = T)
lens <- read.csv("../LM/CT1_round_4/fem 25x 2nd_Final- lense-w-11.5umSpots and equator_Detailed.csv", header = T)
lens <- read.csv("../LM/CT1_round_4/fem 25x 3rd weak GFP- lens cutout-w11.5umspots-wEquator_Detailed.csv", header = T)
lens <- read.csv("../LM/CT1_round_4/fem 25x 4th weak GFP_1-lenscutout-w 11.5um spots- w Equator_Detailed.csv", header = T)



# LM VS uCT ---------------------------------------------------------------

load('../microCT/data/eq_corr/12102019_female.RData')
lens_ref <- lens[!ind_left_lens,]

# least squares fit sphere
xyz_fit <- lens_ref
xyz_fit2 <- xyz_fit^2
Y <- rowSums(xyz_fit2)
X <- cbind(2*xyz_fit, rep(1,dim(xyz_fit)[1]))
sph <- lsfit(X,Y,intercept = FALSE)
r <- sqrt(sph[[1]][4]+sum(sph[[1]][1:3]^2))

# # DEBUG
# nopen3d()
# points3d(lens_ref)
# spheres3d(sph[[1]]['x'],sph[[1]]['y'],sph[[1]]['z'],r, col='gray', alpha=0.2)

lens_ref <- sweep(lens_ref, 2, c(sph[[1]]['x'],sph[[1]]['y'],sph[[1]]['z']))
lens_ref <- lens_ref / r

nopen3d()
points3d(lens_ref)
# spheres3d(sph[[1]]['x'], 0, sph[[1]]['z'],r, col='gray', alpha=0.2)

# LM
load('../LM/aaron.RData')
colnames(lens) <- c('x','y','z')

# least squares fit spher
xyz_fit <- lens
xyz_fit2 <- xyz_fit^2
Y <- rowSums(xyz_fit2)
X <- cbind(2*xyz_fit, rep(1,dim(xyz_fit)[1]))
sph <- lsfit(X,Y,intercept = FALSE)
r <- sqrt(sph[[1]][4]+sum(sph[[1]][1:3]^2))

lens <- sweep(lens, 2, c(sph[[1]]['x'],sph[[1]]['y'],sph[[1]]['z']))
lens <- lens / r

nopen3d()
# spheres3d(sph[[1]]['x'], 0, sph[[1]]['z'],r, col='gray', alpha=0.2)
spheres3d(0,0,0,1, col='gray', alpha=0.2)
points3d(lens, col='red',size=10)


# Hungarian
pt1 <- lens_ref
pt2 <- sweep(lens, 2, c(0.1,-0.1,0))

L1 <- nrow(pt1)
L2 <- nrow(pt2)

Lmax <- max(L1, L2)
Lmin <- min(L1, L2)
mDist <- matrix(ncol = Lmax, nrow = Lmin)
if (Lmin == L1) {
  dr <- pt1
  dc <- pt2
} else {
  dr <- pt2
  dc <- pt1
}
for (j in 1:Lmin) {
  mDist[j,] <- as.matrix(dist(rbind(dr[j,], dc)))[-1,1]
}
i_match <- solve_LSAP(mDist)

if (Lmin == L1) {
  ind_match <- cbind(seq(1,length(i_match)), i_match)
  # ind_match[!(ind_match[,2] %in% seq(1, L2)), 2] <- NA # def unmatched as NA
} else {
  ii <- cbind(i_match, seq(1,length(i_match)) )
  ind_match <- ii[order(ii[,1]),]
  # ind_match[!(ind_match[,1] %in% seq(1, L1)), 1] <- NA
}

nopen3d()
points3d(pt1, size = 5, col = 'gray30')
points3d(pt2, size = 5, col = 'magenta')
linemat <- matrix(t(cbind(pt1[ind_match[,1],], pt2[ind_match[,2],] )), ncol = 3, byrow = T)
segments3d(linemat, color = "grey50")


# function that maps lens to hex grid -----------------------------------------------
# [+v -v +q - q -p + p]

# uCT
load('../microCT/data/eq_corr/12102019_female.RData')
lens_ref <- lens[!ind_left_lens,]

load('../microCT/data/eq_corr/033021_female.RData')
pt <- lens[!ind_left_lens,]

# LM
load('../LM/aaron.RData')
colnames(lens) <- c('x','y','z')

lens <- read.csv("../LM/CT1_round_4/fem 25x 1st weak GFP-Lenscutout-w-spots11.5um_Detailed.csv", header = T)
lens <- lens[, c("Position.X", "Position.Y", "Position.Z")] %>% as.matrix()
colnames(lens) <- c('x','y','z')


# # loop
fn <- c('12102019_female', # match EM
        '12062019_female',
        '12052019_female',
        '033021_male',  # no overlap
        '033021_male_2', # ok
        '033021_female', # bad edge
        '033021_female_2', # large overlap
        '033121_female') # bad edge
# fn <- c("fem 25x 1st weak GFP-Lenscutout-w-spots11.5um_Detailed",
#         "fem 25x 2nd_Final- lense-w-11.5umSpots and equator_Detailed",
#         "fem 25x 3rd weak GFP- lens cutout-w11.5umspots-wEquator_Detailed",
#         "fem 25x 4th weak GFP_1-lenscutout-w 11.5um spots- w Equator_Detailed")

for (f in 1:length(fn)) {
  load(paste("../microCT/data/eq_corr/", fn[f], ".RData", sep=''))
  pt <- lens[!ind_left_lens,]
  
  # lens <- read.csv(paste("../LM/CT1_round_4/", fn[f], ".csv",sep = ''), header = T)
  # lens <- lens[, c("Position.X", "Position.Y", "Position.Z")] %>% as.matrix()
  # colnames(lens) <- c('x','y','z')
  # pt <- lens
  
  # pt <- lens_ref
  
  # data var
  ind_nb <- matrix(ncol = 7, nrow = nrow(pt)) # index of self and 6 nb in [p,v,q,-p,-v,-q]
  ind_xy <- matrix(ncol = 3, nrow = nrow(pt)) # [ind, x, y], use integer coordinates
  dv <- c(1, 1) # delta v
  dp <- c(1, 0)
  dq <- c(0, 1)
  dhex <- rbind(dp,dv,dq,-dp,-dv,-dq)
  
  ind_done <- c() # done checking
  ind_cur <- c() # need to check
  ind_new <- c() # save for next round
  ind_all <- seq(1, nrow(pt))
  
  # first hex, close to the center-of-mass
  # j0 <- 161 #lens LM 80 #lens_ref 391
  j0 <- rowSums(sweep(pt,2,colMeans(pt))^2) %>% which.min()
  zout <- pt[j0,] - colMeans(pt)
  j7 <- rowSums(sweep(pt, 2, pt[j0,])^2) %>% order() %>% head(7)
  pt7 <- pt[j7,]
  
  ## ## -- 2022-2-24
  # origin
  ind_nb[1,1] <- j7[1]
  ind_xy[1,] <- c(j7[1], 0, 0)
  # use vn = c(0,0,1) and zout to assign the other 6 pts
  vv <- sweep(pt7[-1,], 2, pt7[1,])
  # vv %*% vn
  iimax <- which.max(vv[,3])
  ind_nb[1,1+2] <- j7[iimax+1]
  ind_xy[1+2,] <- c(j7[iimax+1], dhex[2,])
  iip <- vv[,3] > 0
  iip[iimax] <- FALSE
  iip <- which(iip) + 1 #ind in j7
  if (length(iip) != 2) {
    print("iip error")
  }
  if (zout %*% cross3D(vv[iimax,], vv[iip[1]-1,]) > 0) {
    ind_nb[1,1+1] <- j7[iip[1]]
    ind_xy[1+1,] <- c(j7[iip[1]], dhex[1,])
    ind_nb[1,1+3] <- j7[iip[2]]
    ind_xy[1+3,] <- c(j7[iip[2]], dhex[3,])
  } else {
    ind_nb[1,1+1] <- j7[iip[2]]
    ind_xy[1+1,] <- c(j7[iip[2]], dhex[1,])
    ind_nb[1,1+3] <- j7[iip[1]]
    ind_xy[1+3,] <- c(j7[iip[1]], dhex[3,])
  }
  
  iimin <- which.min(vv[,3])
  ind_nb[1,1+5] <- j7[iimin+1]
  ind_xy[1+5,] <- c(j7[iimin+1], dhex[5,])
  iim <- vv[,3] < 0
  iim[iimin] <- FALSE
  iim <- which(iim) + 1
  if (length(iim) != 2) {
    print("iim error")
  }
  if (zout %*% cross3D(vv[iimin,], vv[iim[1]-1,]) > 0) {
    ind_nb[1,1+4] <- j7[iim[1]]
    ind_xy[1+4,] <- c(j7[iim[1]], dhex[4,])
    ind_nb[1,1+6] <- j7[iim[2]]
    ind_xy[1+6,] <- c(j7[iim[2]], dhex[6,])
  } else {
    ind_nb[1,1+4] <- j7[iim[2]]
    ind_xy[1+4,] <- c(j7[iim[2]], dhex[4,])
    ind_nb[1,1+6] <- j7[iim[1]]
    ind_xy[1+6,] <- c(j7[iim[1]], dhex[6,])
  }
  
  # ## ##
  # pc7 <- prcomp(pt7)
  # 
  # # +z outward
  # if ( pc7$rotation[,3] %*% zout > 0 ) {
  #   pc7$rotation[,3] <- -pc7$rotation[,3]
  # }
  # if (t(cross3D(pc7$rotation[,3],pc7$rotation[,1])) %*% pc7$rotation[,2] < 0 ) {
  #   pc7$rotation[,2] <- - pc7$rotation[,2]
  # }
  # pt7 <- sweep(pt7, 2, pc7$center) %*% pc7$rotation
  # 
  # # origin
  # ind_nb[1,1] <- j7[1]
  # ind_xy[1,] <- c(j7[1], 0, 0)
  # # right 3 pt
  # ii <- order(pt7[,1]) %>% tail(3)
  # ii2 <- pt7[ii, 2] %>% order()
  # ii3 <- j7[ii][ii2]
  # ind_nb[1,2:4] <- ii3
  # ind_xy[2:4,] <- cbind(ii3, dhex[1:3,])
  # # left 3 pt
  # ii <- order(pt7[,1]) %>% head(3)
  # ii2 <- pt7[ii, 2] %>% order() %>% rev()
  # ii3 <- j7[ii][ii2]
  # ind_nb[1,5:7] <- ii3
  # ind_xy[5:7,] <- cbind(ii3, dhex[4:6,])
  
  
  ind_done <- j7[1]
  ind_cur <- ind_nb[1,-1]
  rnb <- 2 # row number for ind_nb and ind_xy
  rxy <- 8
  while (length(ind_cur) > 0) {
    for (j in ind_cur) {
      ixy <- ind_xy[match(j, ind_xy[,1]), ]
      ind_nb[rnb, 1] <- ixy[1] # home ind
      for (k in 1:3) { #check pairs of nbs in order of p,v,q
        xy <- ixy[2:3] + dhex[k,] #get expected xy
        ii1 <- which(ind_xy[,2] %in% xy[1] & ind_xy[,3] %in% xy[2]) #?exist nb
        xy <- ixy[2:3] - dhex[k,]
        ii2 <- which(ind_xy[,2] %in% xy[1] & ind_xy[,3] %in% xy[2])
        if (length(ii1) + length(ii2) == 2) { #already in ind_xy
          ind_nb[rnb, 2*k+(0:1)] <- c(ind_xy[ii1,1], ind_xy[ii2,1]) # add to ind_nb
        } 
        if (length(ii1) + length(ii2) == 1) { #one is missing
          ii <- ifelse(length(ii1) == 1, ii1, ii2)
          vv <- pt[ixy[1],] - pt[ind_xy[ii,1],] 
          thr <- sqrt(sum(vv^2)) #dist
          xyz <- pt[ixy[1],] + vv #expected position
          dd <- sweep(pt, 2, xyz)^2 %>% rowSums() %>% sqrt()
          ii3 <- which.min(dd) #nearest
          if (min(dd) < thr*0.4 & !(ii3 %in% ind_xy[,1])) { # not too far off AND new 
            if (identical(ii, ii1)) {
              ind_nb[rnb, 2*k+(0:1)] <- c(ind_xy[ii1,1], ii3)
              ind_xy[rxy,] <- c(ii3, ixy[2:3] - dhex[k,])
            } else {
              ind_nb[rnb, 2*k+(0:1)] <- c(ii3, ind_xy[ii2,1])
              ind_xy[rxy,] <- c(ii3, ixy[2:3] + dhex[k,])
            }
            rxy <- rxy + 1
            ind_new <- c(ind_new, ii3)
          } else { #only register the existing nb
            if (identical(ii, ii1)) {
              ind_nb[rnb, 2*k+0] <- ind_xy[ii1,1]
            } else {
              ind_nb[rnb, 2*k+1] <- ind_xy[ii2,1]
            }
          }
        }
      }
      rnb <- rnb + 1
    }
    ind_done <- c(ind_done, ind_cur)
    ind_cur <- ind_new
    ind_new <- c()
  }
  
  # DEBUG
  nopen3d()
  points3d(pt)
  text3d(pt[ind_xy[,1],], texts = paste(ind_xy[,2],ind_xy[,3],sep=','), adj=1.3)
  points3d(pt[ind_xy[ind_xy[,2] %% 3 == 0,1],], size=10, col='magenta')
  

# cont. curvature 2021 ---------------------------------------------------------------

# # - sphere fit
# xyz_data <- pt
# colnames(xyz_data) <- c('x','y','z')
# 
# xyz_data2 <- xyz_data^2
# Y <- rowSums(xyz_data2)
# X <- cbind(2*xyz_data, rep(1,dim(xyz_data)[1]))
# sph <- lsfit(X,Y,intercept = FALSE)
# r_fit <- sqrt(sph[[1]][4]+sum(sph[[1]][1:3]^2)) #radius
# cosph <- sph[["coefficients"]][c('x','y','z')] #center of the sphere
# 
# nopen3d()
# points3d(xyz_data)
# spheres3d(cosph, radius = r_fit, alpha=0.2)


  # - local spheres
  N <- 3 #avg over no. of layers
  xyzr <- matrix(ncol = 4, nrow = nrow(pt)) # center of sphere and radius
  for (j in 1:nrow(pt)) {
    ii <- ind_nb[match(j,ind_nb[,1]),1]
    if (!is.na(ii)) {
      for (k in 1:N) {
        ii2 <- match(ii, ind_nb[,1]) # row ind in ind_nb
        ii <- unique(c(ind_nb[ii2,])) # row inf in pt
      }
      xyz_data <- pt[na.omit(ii), ]
      colnames(xyz_data) <- c('x','y','z')
      xyz_data2 <- xyz_data^2
      Y <- rowSums(xyz_data2)
      X <- cbind(2*xyz_data, rep(1,dim(xyz_data)[1]))
      sph <- lsfit(X,Y,intercept = FALSE)
      r_fit <- sqrt(sph[[1]][4]+sum(sph[[1]][1:3]^2)) #radius
      cosph <- sph[["coefficients"]][c('x','y','z')] #center of the sphere
      # if (length(ii) == 6*N*(N+1)/2+1) {
      xyzr[j,] <- c(cosph, r_fit)
      # }
    }
  }
  sum(is.na(ind_nb[,1]))
  
  # PLOT
  # getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
  # pal_tar <- getPalette(n_lvl - 1)
  nopen3d()
  points3d(pt, size=10, col=plotrix::color.scale(xyzr[,4],c(220,300),90,80,color.spec="hcl"))
  title3d(fn[f])
  
  # hh <- hist(xyzr[,4], breaks = seq(50, 400, by=10), plot = F)
  # points(hh$mids, hh$counts, type='l', col=brewer.pal(5, "Reds")[f+1])
} #loop

dev.new()
# pdf(paste("hist_e", letters[ind_e], "_T4", letters[LL], ".pdf", sep = ""))
hh <- hist(xyzr[,4], breaks = seq(100, 400, by=10), plot = F)
plot(hh$mids, hh$counts, type='l', bty='n', xlim = c(100, 400), xaxt='n',yaxt='n',
     xlab ="radius [um]", ylab='counts', main = "radius of curvature")
axis(1, at = seq(100, 400, by =50), labels = seq(100, 400, by =50) )
# dev.off()


# - circle fit
# choose rows
# ii <- ind_xy[,2] == ind_xy[,3] #v-row
ii <- ind_xy[,2] == -ind_xy[,3] #h-row 

# --- loop
# # v-row
# ind_r_v <- vector(mode="numeric", length=nrow(pt))
# for (k in -16:14) { #uCT, Aaron
#   ii <- ind_xy[,2] == ind_xy[,3] + k
# h-row
ind_r_h <- vector(mode="numeric", length=nrow(pt))
for (k in -28:29) { #uCT
# for (k in -27:28) { #Aaron
  ii <- ind_xy[,2] == -ind_xy[,3] + k

ii2 <- ind_xy[ii, ]
ii2 <- ii2[order(ii2[,2]),]
xyz <- pt[ii2[,1], ]
pca <- prcomp(xyz)
pc3 <- pca$rotation[,3]
com <- pca$center

xyz_pc <- sweep(xyz,2,com) %*% pca$rotation
xy_pc <- xyz_pc[, 1:2]

# # -- circumradius for 3 pts
# R <- vector(mode="numeric", length=nrow(xy_pc))
# N <- 3 # avg over 2*N+1
# for (j in (1+N):(nrow(xy_pc)-N)) {
#   xy <- xy_pc[j + (-N:N), ]
#   # points(xy, pch=16, cex=1+0.1*j, col=j)
#   s1 <- diff(xy[1:2,]) %>% .^2 %>% sum() %>% sqrt()
#   s2 <- diff(xy[2:3,]) %>% .^2 %>% sum() %>% sqrt()
#   s3 <- diff(xy[c(3,1),]) %>% .^2 %>% sum() %>% sqrt()
#   R[j] <- s1*s2*s3 / sqrt((s1+s2+s3)*(-s1+s2+s3)*(s1-s2+s3)*(s1+s2-s3))
# }


# -- circle fit
N <- 3 # avg over 2*N+1
# xyr <- matrix(ncol = 3, nrow = nrow(xy_pc))
for (j in (1+N):(nrow(xy_pc)-N)) {
  xyz_data <- xy_pc[j + (-N:N), ]
  colnames(xyz_data) <- c('x','y')
  xyz_data2 <- xyz_data^2
  Y <- rowSums(xyz_data2)
  X <- cbind(2*xyz_data, rep(1,dim(xyz_data)[1]))
  cir <- lsfit(X,Y,intercept = FALSE)
  r_fit <- sqrt(cir[[1]][3]+sum(cir[[1]][1:2]^2)) #radius cir[[1]] == cir[["coefficients"]]
  cocir <- cir[["coefficients"]][c('x','y')] #center of the sphere
  # xyr[j,] <- c(cocir, r_fit)
  # ind_r_v[ii2[j,1]] <- r_fit
  ind_r_h[ii2[j,1]] <- r_fit
}

}#loop

# PLOT
nopen3d()
points3d(pt)
text3d(pt[ind_r_v !=0,], texts = round(ind_r_v[ind_r_v !=0],0), adj = 1)
# text3d(pt[ind_r_h !=0,], texts = round(ind_r_h[ind_r_h !=0],0), adj = 1)

dev.new()
plot(1:nrow(pt), xyzr[,4], ylim=c(0, 500))
points(ind_r_v, pch=17, col='orange')
points(ind_r_h, pch=16, col='forestgreen')

dev.new()
# pdf(paste("hist_e", letters[ind_e], "_T4", letters[LL], ".pdf", sep = ""))
rr <- xyzr[,4]
rr[is.na(rr)] <- 0
hh <- hist(rr, breaks = seq(0, 500, by=10), plot = F)
plot(hh$mids, hh$density, type='l', bty='n', xlim = c(50, 500), xaxt='n',yaxt='n',
     xlab ="radius [um]", ylab='counts', main = "radius of curvature")
axis(1, at = seq(00, 500, by =50), labels = seq(00, 500, by =50) )
hh <- hist(ind_r_v, breaks = seq(0, 500, by=10), plot = F)
points(hh$mids, hh$density, type='l', col='orange')
hh <- hist(ind_r_h, breaks = seq(0, 500, by=10), plot = F)
points(hh$mids, hh$density, type='l', col='forestgreen')

legend(x = "topright", legend = c("sphere", "v-row", "h-row"), col = c("black","orange","forestgreen"), lwd = 2)    

# DEBUG
nopen3d()
points3d(xyz)
planes3d(pc3[1], pc3[2], pc3[3], -com %*% pc3, alpha=0.4)
arrow3d(com, com+pc3*100, theta= pi/9, n= 4, col= "cyan", type= "rotation")


