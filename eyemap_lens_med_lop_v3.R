# map lens to medulla to LOP 

# copy from lens_med_v7.R
# map Mi1s to uCT lens directly, then use ucl_rot (lens-tip) for eyemap
# eyemap = [Mi1 rownames, ucl_rot rownames]
# extend col structure to LOP via T4b
# 2021-05-05, no longer use CT1
# 2021-05-05, change RF_CT1 to RF_med
# 2021-05-17, rescale RF vec by quantile, centered on com
# 2021-05-26, resample T4 dendrite 
# 2021-12-05, recerter uCT and regenerate eyemap
# ucl_rot_sm smoothed by delaunay, med_xyz by nb_ind !


library(natverse)
library(tidyverse)
library(RColorBrewer)
library(np) # kernel regression to find the base pt on the eye
library(clue) # Hungarian matching
library(alphahull)
# library(ggplot2)
library(deldir)# for smoothing uCT data
library(sp) #in.poly


# setwd("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_eyemap/Lens_Medulla")

# clean everythign up.
rm(list=ls())
#close any open rgl windows
while (rgl.cur() > 0) { close3d() }
while (dev.cur() > 1) { dev.off() }

# # set up for 3d plots based on rgl package
# rgl::setupKnitr()

source("eyemap_func.R")


# load data ---------------------------------------------------------------

# uCT data
# load("data/uCT_12102019.RData")
# load(paste("../microCT/data/eq_corr/", '12102019_female', ".RData", sep=''))  
# load(paste("../microCT/data/eq_corr/", '12102019_female_new', ".RData", sep=''))  
# load(paste("../microCT/data/eq_corr/", '12102019_female', ".RData", sep='')) #new processed, should have '_new' correction

# load(paste0("../microCT/2023_eyemap/", '20230926', '.RData'))
# load(paste0("../microCT/2023_eyemap/", '20231107', '.RData'))
# load(paste0("../microCT/2023_eyemap/", '20240206', '.RData'))
# load(paste0("../microCT/2023_eyemap/", '20240510', '.RData'))
# load(paste0("../microCT/2023_eyemap/", '20240513', '.RData'))
# load(paste0("../microCT/2023_eyemap/", '20240520', '.RData'))
# load(paste0("../microCT/2023_eyemap/", '20240522', '.RData'))
# load(paste0("../microCT/2023_eyemap/", '20240524', '.RData'))
# load(paste0("../microCT/2023_eyemap/", '20240530', '.RData'))
# load(paste0("../microCT/2023_eyemap/", '20240530_dried', '.RData'))
# load(paste0("../microCT/2023_eyemap/", '20240612', '.RData'))
load(paste0("../microCT/2023_eyemap/", '20240701', '.RData'))

# lens_left <- lens
lens_left <- lens[ind_left_lens,]
# lens_mir <- lens_left
# lens_mir[,2] <- -lens_mir[,2]
lens <- lens[!ind_left_lens,] # match lens to EM
ind_Up_lens <- na.omit(match(i_match[ind_Up], rownames(lens)))
# ind_Up_lens_left <- na.omit(match(i_match[ind_Up], rownames(lens_left)))
ind_Down_lens <- na.omit(match(i_match[ind_Down], rownames(lens)))
cone_lens <- i_match[!ind_left_cone]

# re-order ucl_rot rows st same as lens and re-define rownames because mapping is betw Mi1 and lens
# ucl_rot_left <- ucl_rot[order(i_match),][ind_left_lens,]
ucl_rot_sm0 <- ucl_rot_sm[order(i_match),][!ind_left_lens,] # same order as lens now
ucl_rot_sm0_left <- ucl_rot_sm[order(i_match),][ind_left_lens,] # same order as lens now
rownames(ucl_rot_sm0) <- seq(1, nrow(ucl_rot_sm0))
rownames(ucl_rot_sm0_left) <- seq(1, nrow(ucl_rot_sm0_left))
colnames(ucl_rot_sm0) <- c('x','y','z')
colnames(ucl_rot_sm0_left) <- c('x','y','z')
rm(ucl_rot_sm)
rm(ucl)
rm(i_match)


# # CT1 to Mi1 mapping
# load("data/CT1_Mi1.RData")

# # - CT1
# load("data/neu_CT1.RData")

# Mi1
load("data/neu_Mi1.RData")

# # Mi1 numbering
# load('data/ab.RData')

# JFRC2010 mesh
load("data/JFRC2NP.surf.fafb.rda") # cf.ReiserGroup\p_R7R8
ME_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="ME_R")
LOP_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LOP_R")



# chiasm, Only need a reflection  ---------------------------------------------

load('data/chiasm.RData')

# # Mi1 connected to R7/8 which identify chiasm
# 
# anno_chi <- catmaid_query_by_annotation("^Chiasm adjacent$")
# anno_Mi1 <- catmaid_query_by_annotation("^Putative Mi1$")
# chiMi1_skid <- anno_Mi1$skid[anno_Mi1$skid %in% anno_chi$skid]
# chiR_skid <- anno_chi$skid[!(anno_chi$skid %in% chiMi1_skid)]
# chiMi1 <-  read.neurons.catmaid(chiMi1_skid, .progress='text')
# chiR <- read.neurons.catmaid(chiR_skid, .progress='text')
# 
# chiMi1_ep <- list()
# ind_chi <- c()
# for (j in 1:length(chiMi1)) {
#   chiMi1_ep[[j]] <- xyzmatrix(chiMi1[[j]]$d[chiMi1[[j]]$EndPoints,c("X","Y","Z")])
#   nb_dist_chi <- c()
#   for (k in 1:dim(med_xyz)[1]) {
#     nb_dist_chi[k] <- sqrt(min(rowSums(sweep(chiMi1_ep[[j]], 2, med_xyz[k,])^2)))
#   }
#   ind_chi[j] <- which.min(nb_dist_chi)
# }
# # 
# # # PLOT
# # nopen3d()
# # plot3d(chiMi1, WithNodes = F)
# # plot3d(chiR, WithNodes = F, col='gray', lwd=3)
# # points3d(med_xyz, col = 'gray', size = 8)
# # points3d(med_xyz[ind_chi,], col = 'green', size = 13)
# # points3d(med_xyz[match(vaxis_gen(0),eyemap[,2]),], size=16, col='blue')
# # # for (j in 1:length(chiMi1)) {
# # #   points3d(chiMi1_ep[[j]], col = 'blue', size = 5)
# # # }
# # 
# # anno_chi_nb <- catmaid_query_by_annotation("^Chiasm adjacent nbgp3$")
# # # anno_chi_nb <- anno_chi_nb[c(1:3,8), ]
# # chiR_nbgp3 <- read.neurons.catmaid(anno_chi_nb$skid, .progress='text')
# 
# # # SAVE
# # save(chiMi1, chiR, anno_chi, chiR_skid, chiMi1_skid, chiR_nbgp3, file = "data/chiasm.RData")

# project onto a sphere  ------------------------------------------------------------------------------------------

# --- eye lens from microCT
# xyz_data <- lens_uct_data
xyz_data <- lens
colnames(xyz_data) <- c('x','y','z')

xyz_data2 <- xyz_data^2
Y <- rowSums(xyz_data2)
X <- cbind(2*xyz_data, rep(1,dim(xyz_data)[1]))
sph <- lsfit(X,Y,intercept = FALSE)
r_fit <- sqrt(sph[[1]][4]+sum(sph[[1]][1:3]^2)) #radius
cosph <- sph[["coefficients"]][c('x','y','z')] #center of the sphere

# put on a unit sphere
utp <- cbind(xyz_data, sqrt(rowSums(sweep(xyz_data, 2, cosph)^2)))
# names(utp)[4] <- 'r'
colnames(utp) <- c('x','y','z','r')
utp %<>% as_tibble() %>%
  mutate(theta = acos((z-cosph['z'])/r)) %>%
  mutate(phi = 2*pi*((y-cosph['y']) < 0) + (-1)^((y-cosph['y']) < 0)*acos((x-cosph['x'])/r/sin(theta))) %>%
  mutate(ux = sin(theta)*cos(phi), uy = sin(theta)*sin(phi), uz = cos(theta)) %>%
  as.data.frame()

utp_lens <- utp

nopen3d()
points3d(utp_lens[, c('ux','uy','uz')], size = 6)
spheres3d(0,0,0,1, col='grey', alpha=0.2)
planes3d(0,0,1, 0, alpha = 0.2)
planes3d(0,1,0, 0, alpha = 0.2)


# -- EM Mi1
xyz_data <- Mi1_M10_xyz
colnames(xyz_data) <- c('x','y','z')

xyz_data2 <- xyz_data^2
Y <- rowSums(xyz_data2)
X <- cbind(2*xyz_data, rep(1,dim(xyz_data)[1]))
sph <- lsfit(X,Y,intercept = FALSE)
r_fit <- sqrt(sph[[1]][4]+sum(sph[[1]][1:3]^2))
cosph <- sph[["coefficients"]][c('x','y','z')]
# put on a unit sphere
utp <- cbind(xyz_data, sqrt(rowSums(sweep(xyz_data, 2, cosph)^2)))
utp <- as.data.frame(utp)
colnames(utp) <- c('x','y','z','r')
utp %<>% as_tibble() %>%
  mutate(theta = acos((z-cosph['z'])/r)) %>%
  mutate(phi = 2*pi*((y-cosph['y']) < 0) + (-1)^((y-cosph['y']) < 0)*acos((x-cosph['x'])/r/sin(theta))) %>%
  # mutate(theta = pi - theta) %>%
  # mutate(phi = phi + pi) %>%
  mutate(ux = sin(theta)*cos(phi), uy = sin(theta)*sin(phi), uz = cos(theta)) %>%
  as.data.frame()

utp_Mi1 <- utp # replace CT1


# # # plot
# nopen3d()
# spheres3d(cosph['x'],cosph['y'],cosph['z'],r_fit, col='grey', alpha=0.5)
# points3d(xyz_data, col = 'blue', size = 5)
# plot3d(tar, col= 'grey', soma=T, WithNodes = F)
# axes3d(c('x','y','z'))
# title3d('','','x','y','z')

# # # plot
# nopen3d()
# spheres3d(0,0,0,1, col='grey', alpha=0.2)
# points3d(utp[,c('ux','uy','uz')], col = 'blue', size = 5)
# axes3d(c('x','y','z'))
# title3d('','','x','y','z')


# align point sets, manual set --------------------------------------------------------------------------------------

# --- EM Mi1

# lens 33 x 36 col, centero of lens is pt_100
# x-y plane cross 
# pt_100 <- 35 # see 'choose_center_Mi1_cone.png'
# pt_100 <- 333
# pt_100 <- 581 
# pt_100 <- 8 # could be 12 as well, 2021
pt_100 <- 10 #upper eq # or 334 lower eq # 2023
pt_100_Mi1 <- pt_100
# pt_up <- 171
pt_up <- 315
pt_up_Mi1 <- pt_up

# use PR num
ind_up_Mi1 <- c(599,  39, 207, 606, 388,  10,  12, 157, 361, 287, 503, 105)
ind_down_Mi1 <- c(401, 334, 209, 210, 720, 667, 590, 628, 489, 289, 393)

pt_theta90_ID <- c(ind_up_Mi1, ind_down_Mi1)
pt_theta90_ID_Mi1 <- pt_theta90_ID
pt_theta90 <- rbind(c('ux'=0, 'uy'=0, 'uz'=0), utp_Mi1[pt_theta90_ID,c('ux','uy','uz')])
lm_plane <- lm(pt_theta90$ux ~ pt_theta90$uy + pt_theta90$uz + 0)
vn <- c(-1, lm_plane$coefficients[1], lm_plane$coefficients[2]) #norm
vn <- vn/sqrt(sum(vn^2))
if (sum(vn * utp_Mi1[pt_up, c('ux','uy','uz')]) < 0) {
  vn <- -vn
}
# pt_pca <- prcomp(pt_theta90)
# vn <- pt_pca$rotation[,3]
# if (sum(vn * utp_Mi1[pt_up, c('ux','uy','uz')]) < 0) {
#   vn <- -vn
# }

# plot
nopen3d()
points3d(utp_Mi1[,c('ux','uy','uz')], col='blue')
points3d(utp_Mi1[Mi1_ind_PR[[8]],c('ux','uy','uz')], size = 10, col='blue')
points3d(utp_Mi1[Mi1_ind_PR[[7]],c('ux','uy','uz')], size = 10, col='light blue')
points3d(utp_Mi1[unlist(Mi1_ind_PR[c(6,5)]),c('ux','uy','uz')], size = 10, col='gold2')
points3d(utp_Mi1[pt_100, c('ux','uy','uz')], col='cyan', size = 16)
points3d(utp_Mi1[pt_up, c('ux','uy','uz')], col='red', size = 16)
points3d(pt_theta90, col='red', size = 13)
arrow3d(c(0,0,0), vn*1, theta = pi/6,n = 4, col="blue", type = "rotation")
planes3d(vn[1],vn[2],vn[3],0)
spheres3d(0,0,0,1, col='grey', alpha=0.2)
# identify3d(utp_Mi1[,c('ux','uy','uz')], adj=c(-.5,-.5))

# plot3d(Mi1[ii]), use 
# ii <- sapply(ind_up_Mi1, function(x) x + sum(sub_Mi1_ind <= x))

# rotation to the new coordinate
vx <- utp_Mi1[pt_100, c('ux','uy','uz')]
utp_Mi1_rot <- rot_zx0(utp_Mi1[,c('ux','uy','uz')], c(vn), c(vx))
pt_theta90_rot <- rot_zx0(pt_theta90, c(vn), c(vx))

# # plot
nopen3d()
spheres3d(0,0,0,1, col='grey', alpha=0.2)
points3d(pt_theta90, col='red')
# arrow3d(c(0,0,0), vn*1, theta = pi/6,n = 4, col="blue", type = "rotation")
# planes3d(vn[1],vn[2],vn[3],0, col = "grey")
# arrow3d(c(0,0,0), t(vn1)*1, theta = pi/6,n = 4, col="brown", type = "rotation")
# arrow3d(c(0,0,0), t(vn2)*1, theta = pi/6,n = 4, col="red", type = "rotation")
# arrow3d(c(0,0,0), t(cpt3)*1, theta = pi/6,n = 4, col="green", type = "rotation")
points3d(utp_Mi1_rot, col='red', size = 5)
points3d(pt_theta90_rot, col='green', size = 8)
points3d(utp_Mi1_rot[pt_up, , drop=F], col='black', size = 15)
points3d(matrix(pt_theta90_rot[as.character(pt_100),], ncol = 3), col='black', size = 10)
axes3d(c('x','y','z')); title3d('','','x','y','z')


# -- chiasm , just need a reflection (enough if a single plane )
utp_Mi1_rot_chiasm <- utp_Mi1_rot
utp_Mi1_rot_chiasm[,'Y'] <- -utp_Mi1_rot_chiasm[,'Y']
# 
# utp_CT1_rot_chiasm <- utp_CT1_rot
# utp_CT1_rot_chiasm[,'Y'] <- -utp_CT1_rot_chiasm[,'Y']


# --- lens from micro CT

# # choose 100 and up
nopen3d()
points3d(lens)
points3d(lens[ind_Up_lens,], size=10, col='red')
points3d(lens[ind_Down_lens,], size=10, col='blue')
identify3d(lens)


# pt_100 <- 388 # 12102019
# pt_up <- 512 # 12102019

# pt_100 <- 642 # 20230926, pre-84
# pt_up <- 813
# pt_100 <- 647 # 20230926, 84
# pt_up <- 785

# pt_100 <- 417 # 20231107
# pt_up <- 565

# pt_100 <- 434 # 20240206
# pt_up <- 688

# pt_100 <- 366 # 20240510
# pt_up <- 624

# pt_100 <- 397 # 20240513
# pt_up <- 597

# pt_100 <- 403 # 20240520
# pt_up <- 667

# pt_100 <- 377 # 20240522
# pt_up <- 611

# pt_100 <- 411 # 20240524
# pt_up <- 673

# pt_100 <- 373 # 20240530
# pt_up <- 638

# pt_100 <- 389 # 20240530_dried
# pt_up <- 622

# pt_100 <- 398 # 20240612
# pt_up <- 669

pt_100 <- 395 # 20240701
pt_up <- 659

pt_100_lens <- pt_100
pt_up_lens <- pt_up

#  a little off c(551, 338, 313)
# ind_90 <- ind_down_uct_lens | ind_up_uct_lens
# ind_90[c(558, 333, 308)] <- FALSE
# pt_theta90 <- rbind(c('ux'=0, 'uy'=0, 'uz'=0), utp_lens[ind_90,c('ux','uy','uz')])


pt_theta90_ID_lens <- c(ind_Up_lens, ind_Down_lens)
pt_theta90 <- rbind(c('ux'=0, 'uy'=0, 'uz'=0), utp_lens[pt_theta90_ID_lens,c('ux','uy','uz')])
# lm_plane <- lm(pt_theta90$ux ~ pt_theta90$uy + pt_theta90$uz + 0)
# vn <- c(-1, lm_plane$coefficients[1], lm_plane$coefficients[2]) #norm
lm_plane <- lm(pt_theta90$uz ~ pt_theta90$ux + pt_theta90$uy + 0)
vn <- c(lm_plane$coefficients[1], lm_plane$coefficients[2], -1) #norm
vn <- vn/sqrt(sum(vn^2))
# which side is up
if (sum(vn * utp_lens[pt_up, c('ux','uy','uz')]) < 0) {
  vn <- -vn
}
# pt_pca <- prcomp(pt_theta90)
# vn <- pt_pca$rotation[,3]
# if (sum(vn * utp_lens[pt_up, c('ux','uy','uz')]) < 0) {
#   vn <- -vn
# }

# plot, find center and up
nopen3d()
points3d(utp_lens[,c('ux','uy','uz')], col='grey', size = 5)
# points3d(utp_lens[ind_Up_lens,c('ux','uy','uz')], col='blue', size = 8)
# points3d(utp_lens[ind_Down_lens,c('ux','uy','uz')], col='red', size = 8)
# points3d(utp_lens[pt_theta90_ID_lens, c('ux','uy','uz')], col='green', size = 15)
# identify3d(utp_lens[,c('ux','uy','uz')])
spheres3d(0,0,0,1, col='grey', alpha=0.2)
points3d(pt_theta90, col='cyan', size = 13)
points3d(utp_lens[pt_up, c('ux','uy','uz')], col='red', size = 13)
arrow3d(c(0,0,0), vn*1, theta = pi/6,n = 4, col="blue", type = "rotation")
planes3d(vn[1],vn[2],vn[3],0, color = "grey")
axes3d(c('x','y','z'))
title3d('','','x','y','z')

# rotation to the new coordinate
vx <- utp_lens[pt_100, c('ux','uy','uz')]
utp_lens_rot <- rot_zx0(utp_lens[,c('ux','uy','uz')], c(vn), c(vx))
pt_theta90_rot <- rot_zx0(pt_theta90, c(vn), c(vx))


# plot
nopen3d()
spheres3d(0,0,0,1, col='grey', alpha=0.2)
points3d(pt_theta90, col='red')
# arrow3d(c(0,0,0), vn*1, theta = pi/6,n = 4, col="blue", type = "rotation")
# planes3d(vn[1],vn[2],vn[3],0, col = "grey")
# arrow3d(c(0,0,0), t(vn1)*1, theta = pi/6,n = 4, col="brown", type = "rotation")
# arrow3d(c(0,0,0), t(vn2)*1, theta = pi/6,n = 4, col="red", type = "rotation")
# arrow3d(c(0,0,0), t(cpt3)*1, theta = pi/6,n = 4, col="green", type = "rotation")
points3d(utp_lens_rot, col='blue', size = 5)
points3d(pt_theta90_rot, col='green', size = 10)
points3d(matrix(pt_theta90_rot[paste(pt_100),], ncol = 3), col='black', size = 15)
points3d(utp_lens_rot[pt_up, , drop=F], col='black', size = 15)
planes3d(0,0,1,0, col = "grey", alpha = 0.2)
planes3d(1,0,0,0, col = "grey", alpha = 0.2)
axes3d(c('x','y','z')); title3d('','','x','y','z')
# identify3d(utp_lens_rot)


# # check alignment, EM CT1 and lens
# nopen3d()
# spheres3d(0,0,0,1, col='grey', alpha=0.2)
# points3d(utp_Mi1_rot_chiasm, size = 6)
# points3d(utp_Mi1_rot_chiasm[ind_up_Mi1,], col='blue', size = 14)
# points3d(matrix(utp_Mi1_rot_chiasm[pt_100_Mi1,],ncol=3), col='cyan', size = 20)
# points3d(utp_lens_rot, col='brown', size = 6)
# points3d(utp_lens_rot[ind_Up_lens,], col='pink', size = 14)
# points3d(matrix(utp_lens_rot[pt_100_lens,],ncol=3), col='red', size = 20)
# axes3d(c('x','y','z')); title3d('','','x','y','z')


# # parameter for generating hex lattice
# pteq <- pt_theta90_rot[-1,]
# pteq <- pteq[order(pteq[,2]),]
# hdx <- mean(sqrt(rowSums(diff(as.matrix(pteq))^2))) #this x is y-axis in the 3d plot
# hdy <- hdx / sqrt(3)
# hdia <- dim(pteq)[1]
# pteqN <- dim(pteq)[1]
# hdax <- (atan(pteq[pteqN,2]/pteq[pteqN,1]) - atan(pteq[1,2]/pteq[1,1]))/(pteqN - 1) #x is y-axis in the 3d plot
# hday <- hdax / sqrt(3)



# establish hex in lens -------------------------------------------------------------------

pt <- lens

# ls <- reghex(pt, j7 = c(pt_100_lens, 689, 693, 643, 596, 591, 636) ) # 20230926, pre-84
# ls <- reghex(pt, j7 = c(pt_100_lens, 693, 697, 648, 603, 600, 646) ) # 20230926
# ls <- reghex(pt, j7 = c(pt_100_lens,409,444,450,424,388,382) ) # 20231107
# ls <- reghex(pt, j7 = c(pt_100_lens,457, 468, 446, 410, 401, 421) ) # 202402026
# ls <- reghex(pt, j7 = c(pt_100_lens,378, 400,387, 353, 331, 344) ) # 20240510
# ls <- reghex(pt, j7 = c(pt_100_lens,417,430,412,377,362,381) ) # 20240513
# ls <- reghex(pt, j7 = c(pt_100_lens,416, 436, 423, 388, 368, 382 ) ) # 20240520
# ls <- reghex(pt, j7 = c(pt_100_lens,394, 412, 396, 364, 345, 360) ) # 20240522
# ls <- reghex(pt, j7 = c(pt_100_lens,425, 447, 435, 400, 377, 389 ) ) # 20240524
# ls <- reghex(pt, j7 = c(pt_100_lens,394, 407, 385, 353, 339, 364 ) ) # 20240530
# ls <- reghex(pt, j7 = c(pt_100_lens,394, 422, 420, 388, 356, 361 ) ) # 20240530_dried
# ls <- reghex(pt, j7 = c(pt_100_lens,422, 431, 409, 375, 365, 388 ) ) # 20240612
ls <- reghex(pt, j7 = c(pt_100_lens,419, 430, 406, 372, 362, 385 ) ) # 2024

nb_ind <- ls[[1]]
ind_xy <- ls[[2]]


# # DEBUG
# nopen3d()
# points3d(pt,size=7)
# text3d(pt[ind_xy[,1],], texts = paste(ind_xy[,2],ind_xy[,3],sep='|'), adj=1.3)
# points3d(pt[ind_xy[ind_xy[,3] %% 3 == 0,1],], size=10, col='magenta')


# manual  matching, new, recente ---------------------------------------------------------------------------------
# just EM Mi1 to uCT lens

# # 2023
# med_ixy <- sweep(med_ixy, 2, c(0,1,0), '+')
# save(med_ixy, file = "data/med_ixy.RData")

load('data/med_ixy.RData')

# shifting data
utp_lens_rot_shift <- sweep(utp_lens_rot,2,c(+1,0,0)) #NOTE the small area one should be on the convex side


# # 2021
# x <- as.data.frame(ind_xy)
# colnames(x) <- c('idx', 'x','y')
# y <- as.data.frame(med_ixy)
# colnames(y) <- c('idy', 'x','y')
# xy <- merge(x, y , by=c('x','y'))
# 
# lens_Mi1 <- xy[,3:4] %>% as.matrix()

# # 2023
lens_Mi1 <- matrix(ncol = 2, nrow = nrow(med_ixy))  # [lens index, Mi1 index ]
for (j in 1:nrow(med_ixy)) {
  ii <- which(med_ixy[j,2] == ind_xy[,2] & med_ixy[j,3] == ind_xy[,3])
  if (length(ii) == 1) {
    lens_Mi1[j,] <- c(ind_xy[ii,1], med_ixy[j,1])
  } else {
    lens_Mi1[j,] <- c(NA, med_ixy[j,1])
  }
}
# append un-matched Mi1
lens_Mi1_all <- rbind(lens_Mi1, cbind(ind_xy[!(ind_xy[,1] %in% na.omit(lens_Mi1[,1])), 1], NA))
# only matched
lens_Mi1 <- lens_Mi1[!is.na(lens_Mi1[,1]), ]


# -- PLOT matching
nopen3d()
points3d(utp_Mi1_rot_chiasm, size = 4, col = 'black')
# points3d(utp_Mi1_rot, size = 4, col = 'black')
points3d(utp_lens_rot_shift, size = 4, col = 'red')
linemat <- matrix(t(cbind(utp_Mi1_rot_chiasm[lens_Mi1[,2],],
                          utp_lens_rot_shift[lens_Mi1[,1],])), ncol = 3, byrow = T)
segments3d(linemat, color = "grey", lwd=2)

# exclu
# points3d(utp_lens_rot_shift[!(seq(1,nrow(utp_lens_rot_shift)) %in% lens_Mi1[,1]),], size = 15, col = 'cyan')
# points3d(utp_Mi1_rot_chiasm[!(seq(1,nrow(utp_Mi1_rot_chiasm)) %in% lens_Mi1[,2]),], size = 15, col = 'magenta')
points3d(utp_lens_rot_shift[lens_Mi1_all[is.na(lens_Mi1_all[,2]),1],], size = 15, col = 'cyan')
points3d(utp_Mi1_rot_chiasm[lens_Mi1_all[is.na(lens_Mi1_all[,1]),2],], size = 15, col = 'magenta')

# # check, 2021
# points3d(utp_lens_rot_shift[ya[[1]],], size = 8, col = 'cyan')
# points3d(utp_Mi1_rot_chiasm[yb[[1]],], size = 8, col = 'cyan')
# 
# points3d(utp_lens_rot_shift[xa_ls[[1]][[10]],], size = 10, col = 'blue')
# points3d(utp_Mi1_rot_chiasm[xb_ls[[1]][[10]],], size = 10, col = 'gold2')
# 
# 
# points3d(utp_lens_rot_shift[eyemap[match(xa_ls[[4]][[10]], eyemap[,2]),2],], size = 15, col = rainbow(5))
# points3d(utp_Mi1_rot_chiasm[eyemap[match(xb_ls[[4]][[10]], eyemap[,1]),1],], size = 15, col = rainbow(5))


# -- eyemap, [med, lens (cone or ucl indeed)]
# eyemap <- cbind(match(lens_Mi1[,2], CT1_Mi1), lens_Mi1[,1])
eyemap <- cbind(lens_Mi1[,2], lens_Mi1[,1])



# cont. re-order st. eyemap entry is neuron index, eyemap index is pts --------

med_xyz <- Mi1_M10_xyz[eyemap[,1],] 
rownames(med_xyz) <- eyemap[,1]
# ucl_rot <- ucl_rot[eyemap[,2], ] 
ucl_rot_sm <- ucl_rot_sm0[eyemap[,2], ] 
ind_Up_ucl <- match(ind_Up_lens, eyemap[,2])
ind_Down_ucl <- match(ind_Down_lens, eyemap[,2])

nb_ind <- matrix(match(nb_ind, eyemap[,2]), ncol = 7, byrow = F) # now it's eyemap index
nb_ind <- nb_ind[!is.na(nb_ind[,1]), ]

# save an un-ordered copy, including unmatched points
lens_ixy <- ind_xy

ind_xy[,1] <- match(ind_xy[,1], eyemap[,2])
ind_xy <- ind_xy[!is.na(ind_xy[,1]), ]

# Note 2023-11-01
ord <- order(nb_ind[,1])
nb_ind <- nb_ind[ord,] 
ind_xy <- ind_xy[ord,]

# ioa
# right
dd <- ucl_rot_sm
nb_dist_ucl <- matrix(ncol = 7, nrow = nrow(dd))
for (j in 1:nrow(dd)) {
  nb_dist_ucl[nb_ind[j,1],] <- c(
    nb_ind[j,1], acos(dd[nb_ind[j,-1], ] %*% dd[nb_ind[j,1],]) )
}

# lens_ixy[,1] <- match(ind_xy[,1], eyemap[,2])
# lens_ixy <- lens_ixy[!is.na(lens_ixy[,1]),]


# # DEBUG 2023
# nopen3d()
# points3d(ucl_rot_sm,size=7)
# text3d(ucl_rot_sm[lens_ixy[,1],], texts = paste(lens_ixy[,2],lens_ixy[,3],sep='|'), adj=1.3)
# points3d(ucl_rot_sm[lens_ixy[lens_ixy[,3] %% 3 == 0,1],], size=10, col='magenta')
# 
# nopen3d()
# points3d(ucl_rot_sm0,size=7)
# text3d(ucl_rot_sm0[ind_xy[,1],], texts = paste(ind_xy[,2],ind_xy[,3],sep='|'), adj=1.3)
# points3d(ucl_rot_sm0[ind_xy[ind_xy[,3] %% 3 == 0,1],], size=10, col='magenta')
# points3d(ucl_rot_sm0[ind_xy[ind_xy[,2] %% 3 == 0,1],], size=10, col='cyan')

# nb dist med [center, nbco]  -------------------------------------------

# 2023
nb_dist_med <- matrix(ncol = 7, nrow = nrow(eyemap))
for (j in 1:nrow(eyemap)) {
  nb_dist_med[nb_ind[j,1],] <- c(
    nb_ind[j,1],
    sqrt(rowSums(sweep(med_xyz[nb_ind[j,-1], ], 2, med_xyz[nb_ind[j,1],,drop=F])^2))
  )
}

nbco <- rbind(c(+1, 0), # p-axis
              c(+1,+1), # v-axis
              c(0, +1), # q-axis
              c(-1, 0),
              c(-1,-1),
              c(0, -1)
)

# SAVE
save(nb_ind, ind_xy, nb_dist_med, nb_dist_ucl, nbco, file = "data/hexnb_ind_dist.RData")


#  add aux points beyond the boundary for regression ----------------------

pts <- ucl_rot_sm

# rotate to southern hemisphere
vcom <- colMeans(pts)
vcom <- vcom / sqrt(sum(vcom^2))
R_mat <- quaternion3D( cross3D(vcom, c(0,0,-1)), acos(vcom %*% c(0,0,-1)) /pi*180 )
pts_z <- as.matrix(pts) %*% t(R_mat)

# stereographic projection
x <- pts_z[,1] / (1 - pts_z[,3])
y <- pts_z[,2] / (1 - pts_z[,3])
stereo <- cbind(x, y)

dtdt <- deldir(stereo[,1], stereo[,2])

# -- deal with edge points,
# if the boundary pt, then check the 2 circumcenters sharing this edge,
# remove the one with no corresponding circumcenter
# Alt,  deal with edge points with dummy pt. NOT tried ?

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

# remove long edge
el <- pts[delsgs_mod$ind1,] - pts[delsgs_mod$ind2,]
el <- sqrt(rowSums(el^2))
delsgs_mod <- delsgs_mod[el < 0.27, ]

stereo_as <- ashape(stereo, alpha = 0.2)
xy_edge <- stereo_as$edges[,1:6]
stereo_poly <- mkpoly(xy_edge)[[1]]

ind_bd <- stereo_poly$ind1 # boundary pt index

stereo_bd <- stereo_poly[,3:4]
stereo_bd <- rbind(stereo_bd, stereo_bd[1,])

# -- add aux pt
pt_aux <- matrix(ncol = 3, nrow = 0)
pt_med_aux <- matrix(ncol = 3, nrow = 0)
for (iibd in ind_bd) {
  # neighbors of this point
  ind_nb <- c(delsgs_mod$ind2[delsgs_mod$ind1 %in% iibd], delsgs_mod$ind1[delsgs_mod$ind2 %in% iibd])
  # aux pt for each neighbor
  for (iinb in ind_nb) {
    # ii_med <- match(rownames(pts)[iibd],eyemap[,2])
    # # ii_med <- lens_Mi1[match(rownames(pts)[iibd],lens_Mi1[,1]), 2]
    # if (!is.na(ii_med)) {
      xyz <- ucl_rot_sm[iibd,] + ucl_rot_sm[iibd,] - ucl_rot_sm[iinb,] 
      pt_aux <- rbind(pt_aux, matrix(xyz/sqrt(sum(xyz^2)),ncol = 3))
      
      # xyz <- med_xyz[ii_med,] + med_xyz[ii_med,] - med_xyz[ii_med,] 
      xyz <- med_xyz[iibd,] + med_xyz[iibd,] - med_xyz[iinb,] 
      pt_med_aux <- rbind(pt_med_aux, xyz)
    # }
  }
}

# remove pts within the polygon
pt_aux_z <- as.matrix(pt_aux) %*% t(R_mat)
# stereographic projection
x <- pt_aux_z[,1] / (1 - pt_aux_z[,3])
y <- pt_aux_z[,2] / (1 - pt_aux_z[,3])
stereo_aux <- cbind(x, y) 

ii_inpoly <- point.in.polygon(stereo_aux[,1], stereo_aux[,2], stereo_bd[,1], stereo_bd[,2])

pt_aux <- pt_aux[!ii_inpoly,]
pt_med_aux <- pt_med_aux[!ii_inpoly,]

# remove pts close to ucl_rot_sm
ii_dd <- c()
for (ii in 1:nrow(pt_aux)) {
  ii_dd <- c(ii_dd, min(sqrt(rowSums(sweep(ucl_rot_sm, 2, pt_aux[ii,])^2))) > 0.08)
}
pt_aux <- pt_aux[ii_dd,]
pt_med_aux <- pt_med_aux[ii_dd,]

# remove clustered pt
ii_o <- order(sqrt(rowSums(sweep(pt_aux, 2, colMeans(ucl_rot_sm))^2)))
pt_aux <- pt_aux[ii_o,]
pt_med_aux <- pt_med_aux[ii_o,]
ii <- 1
while (ii < nrow(pt_aux)) {
  dd <- sqrt(rowSums(sweep(pt_aux, 2, pt_aux[ii,])^2))
  ii_dd <- dd > 0.1
  ii_dd[ii] <- TRUE
  pt_aux <- pt_aux[ii_dd,]
  pt_med_aux <- pt_med_aux[ii_dd,]
  ii <- ii + 1
}
pt_aux <- sweep(pt_aux,1,sqrt(rowSums(pt_aux^2)),'/') #normalize

# remove pts too far from med_xyz
ii_dd <- c()
for (ii in 1:nrow(pt_med_aux)) {
  ii_dd <- c(ii_dd, min(sqrt(rowSums(sweep(med_xyz, 2, pt_med_aux[ii,])^2))) < 1e4)
}
pt_aux <- pt_aux[ii_dd,]
pt_med_aux <- pt_med_aux[ii_dd,]

nopen3d()
points3d(ucl_rot_sm)
points3d(pt_aux, size = 8, col='red')
text3d(pt_aux, texts = seq(1,nrow(pt_aux)), adj = -1)
# points3d(ucl_rot_sm[match(xa_ls[[1]][[1]], eyemap[,2]),], size = 10, col ='green')
# points3d(ucl_rot_sm[match(xa_ls[[2]][[1]], eyemap[,2]),], size = 10, col ='green')

nopen3d()
points3d(med_xyz)
points3d(pt_med_aux, size = 8, col='red')
text3d(pt_med_aux, texts = seq(1,nrow(pt_med_aux)), adj = -1.5)
# points3d(med_xyz[match(xb_ls[[1]][[1]], eyemap[,1]),], size = 10, col ='green')
# points3d(med_xyz[match(xb_ls[[2]][[1]], eyemap[,1]),], size = 10, col ='green')

# --- re-define
Npt <- nrow(eyemap)
rownames(pt_aux) <- seq(Npt + 1, Npt + nrow(pt_aux))
rownames(pt_med_aux) <- seq(Npt + 1, Npt + nrow(pt_med_aux))

# ### aux
ucl_rot_aux <- rbind(ucl_rot_sm, pt_aux)
med_xyz_aux <- rbind(med_xyz, pt_med_aux)
eyemap_aux <- rbind(eyemap, cbind(seq(nrow(eyemap)+1, nrow(eyemap)+nrow(pt_aux)),
                              seq(nrow(eyemap)+1, nrow(eyemap)+nrow(pt_aux))) )

# ### NO aux
# ucl_rot_aux <- ucl_rot_sm
# med_xyz_aux <- med_xyz
# eyemap_aux <- eyemap


# SAVE eyemap -------------------------------------------------------------

save(eyemap, med_xyz, delsgs_mod, ucl_rot_sm, lens_Mi1, lens_ixy,
     utp_lens_rot, utp_lens_rot_shift, utp_Mi1_rot, utp_Mi1_rot_chiasm, ind_Up_ucl, ind_Down_ucl,
     ucl_rot_aux, med_xyz_aux, eyemap_aux, Npt,
     file = "data/eyemap.RData")

# overlapping T4 subtypes to confirm local orthogonality ------------------------

library(stringr)

anno_overlap_M <- catmaid_query_by_annotation("^T4 overlap M$")
anno_overlap_LD <- catmaid_query_by_annotation("^T4 overlap LD$") #misnomer 
anno_overlap_LV <- catmaid_query_by_annotation("^T4 overlap LV$") #misnomer 
anno_overlap_VP <- catmaid_query_by_annotation("^T4 overlap VP$")
anno_overlap <- rbind(anno_overlap_M, anno_overlap_LD, anno_overlap_LV, anno_overlap_VP)
T4_overlap_pos <- c('M','AD','AV','LV')

ii <- regexpr("T4", anno_overlap$name) #start index
anno_overlap$ti <- match(substring(anno_overlap$name, ii+2, ii+2), letters)

T4_overlap <-  read.neurons.catmaid(anno_overlap[,"skid"], .progress='text')
T4_overlap <- kinkCT1(T4_overlap)# correct kink CT1

T4_overlap_N <- c(nrow(anno_overlap_M), nrow(anno_overlap_LD), nrow(anno_overlap_LV), nrow(anno_overlap_VP) )

# SAVE
# save(anno_overlap, T4_overlap_pos, T4_overlap_N, T4_overlap, file = "data/T4_hs.RData")


# # - exemple
# neu <- T4_overlap_M
# anno <- anno_overlap_M
# dend_com <- matrix(ncol = 3, nrow = length(neu))
# for (j in 1:length(neu)) {
#   tar <- neu[[j]]
#   ind_D = match(tar$tags$"dendrite start" , tar$d$PointNo)
#   mat_D = tar$d[ind_D,]
#   node_dend <- child_node(tar, mat_D)
#   dend_com[j,] <- colMeans(xyzmatrix(node_dend))
# }
# nopen3d()
# text3d(dend_com, texts = paste(seq(1,length(neu)), word(anno$name, 2)), adj = -0.5)
# 
# # 7   9880571  Putative T4a 9880571 BG neuron  9880570
# # 10 11364712 Putative T4b 11364712 AT neuron 11364711
# # 8  11301435 Putative T4c 11301435 AT neuron 11301434
# # 9  11347885 Putative T4d 11347885 AT neuron 11347884


T4_col <- brewer.pal(7, "RdYlBu")[c(1,3,6,7)]
pals <- brewer.pal(9, "Spectral")
# col_col <- pals[c(1,1,1,7,7,7,7,7,9,9,9,9,9,9,9)]
col_col <- pals[c(rep(1,20), rep(9,40))]

# # PLOT
# nopen3d()
# shade3d(ME_msh, alpha=0.1)
# shade3d(LOP_msh, alpha=0.1)
# plot3d(T4_overlap_M, lwd = 3)
# plot3d(T4_overlap_LD, lwd=3)
# plot3d(T4_overlap_LV, lwd = 3)
# plot3d(T4_overlap_VP, lwd = 3)
# # plot3d(CT1, col = 'yellow3')
# points3d(med_xyz, col = 'pink3', size = 3)
# # points3d(med_xyz[ind_chi,], col = 'green', size = 15)


# - method 5, Strahler order

neu <- T4_overlap
anno <- anno_overlap

# nbhd <- 20000 # radius of inclusion
nbhd_N <- 1+8+16 #1+6 +12 # num of nbs
ls_nb <- list() # ID of nb column labeled by CT1
RF_med <- matrix(nrow = length(neu), ncol = 6) # vector [base,  center of mass]
RF_lens <- matrix(nrow = length(neu), ncol = 6) # vector [base,  center of mass]
RF_var <- c()
count_multi_root <- c()
for (j in 1:length(neu)) {
  tar <- neu[[j]]
  # tar <- resample(tar, stepsize = 100) #lose tags
  ind_D = match(tar$tags$"dendrite start" , tar$d$PointNo)
  mat_D = tar$d[ind_D,]
  xyz_D <- xyzmatrix(mat_D)
  if (nrow(xyz_D) > 1) {
    count_multi_root <- count_multi_root + 1
  }
  xyz_D_avg <- colMeans(xyz_D) #for multiple dendrite tags
  node_dend <- child_node(tar, mat_D)
  node_xyz <- xyzmatrix(node_dend)
  
  # - find nearby columns, use dendrite start position
  ls_nb[[j]] <- sweep(med_xyz, 2, xyz_D_avg, '-')^2 %>% rowSums() %>% order() %>% head(nbhd_N)

  # - orientation of CT1 via Strahler order
  L <- sqrt(sum((colMeans(node_xyz) - xyz_D_avg)^2))
  RF_med[j,] <- neu_dir_SO(xneu = tar, xroot = "dendrite start", so = c(2,3),
                           sumType = 'avg', scaling = L)
  # UPDATE 20210517, rescale by quantile, centered on com
  v <- RF_med[j,4:6]-RF_med[j,1:3]
  v <- v/sqrt(sum(v^2))
  pv <- sweep(node_xyz,2,RF_med[j,1:3],'-') %*% v
  nq <- quantile(pv, c(0.25,0.75)) %>% diff()
  RF_med[j,] <- c(colMeans(node_xyz)- v*nq/2, colMeans(node_xyz)+ v*nq/2)

  # map to ucl_rot_sm
  xyz_vtail <- matrix(ncol = 3)  #base of vec
  xyz_vhead <- matrix(ncol = 3) #head of vec
  xyz_vtail_eval <- data.frame(mc.x = RF_med[j,1], mc.y = RF_med[j,2], mc.z = RF_med[j,3])
  xyz_vhead_eval <- data.frame(mc.x = RF_med[j,4], mc.y = RF_med[j,5], mc.z = RF_med[j,6])
  
  # -- kernel regression
  for (k in 1:3) {
    npdata <- data.frame(mc = med_xyz[ls_nb[[j]],], ec = ucl_rot_sm[ls_nb[[j]],k]) # with re-defined ucl_rot_sm and med_xyz
    # bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'generalized_nn', regtype= 'll')
    # bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'generalized_nn', regtype= 'lc')
    bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'fixed', regtype= 'll')
    # bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
    model_np <- npreg(bw)
    # model_np <- npreg(ec ~ mc.x + mc.y + mc.z, data = npdata)
    xyz_vtail[k] <- predict(model_np, newdata = xyz_vtail_eval)
    xyz_vhead[k] <- predict(model_np, newdata = xyz_vhead_eval)
  }
  utp_vtail_norm <- sqrt(sum(xyz_vtail^2))
  RF_lens[j,] <- c(xyz_vtail, xyz_vhead) / utp_vtail_norm  #normalize to vtail
  v1 <- RF_lens[j,1:3]
  v2 <- RF_lens[j,4:6]
  RF_lens[j,4:6] <- v2 - c((v2 - v1) %*% v1) * v1 #take tangent
}

RF_lens_overlap <- RF_lens
RF_med_overlap <- RF_med


# PLOT RF on lens, vec dir = optic flow
RF_lens <- RF_lens_overlap
anno <- anno_overlap

nopen3d()
points3d(ucl_rot_sm)
points3d(ucl_rot_sm[vaxis_gen(1, ind_xy),], size=12, col=pal_axes[3])
for (j in 1:13) {
  if (j %% 2 == 1) {
    points3d(ucl_rot_sm[ind_xy[ind_xy[,2] == j,1], ], size=12, col=pal_axes[2])
    points3d(ucl_rot_sm[ind_xy[ind_xy[,2] == -j,1], ], size=12, col=pal_axes[2])}
}

# points3d(ucl_rot, col = col_col, size = 5)
for (j in 1:nrow(RF_lens)) {
  # arrow3d(RF_lens[j,4:6], RF_lens[j,1:3], theta = pi / 8, n = 4, col = 'blue', alpha = 0.3, type = "rotation")
  # specify type
  ii <- regexpr("T4", anno[j,"name"]) #start index
  type_ii <- match(substring(anno[j,"name"], ii+2, ii+2), letters)
  arrow3d(RF_lens[j,4:6], RF_lens[j,1:3], theta = pi / 9, n = 4, col = T4_col[type_ii], alpha = 0.9, type = "rotation", lit=F)
  # pch3d(matrix(RF_lens[j,1:3], ncol = 3), col = "brown", cex = 2, pch = 96+type_ii)
}
title3d(main = "Strahler")
# axes3d(c('x','y','z')); title3d('','','x','y','z')
points3d(ucl_rot_sm[match(ind_Up_lens, eyemap[,2]),], col = "gold2", size = 15)
points3d(ucl_rot_sm[match(ind_Down_lens, eyemap[,2]),], col = "cyan", size = 15)


# PLOT, med
RF_med <- RF_med_overlap
anno <- anno_overlap

nopen3d()
points3d(med_xyz)
points3d(med_xyz[vaxis_gen(1, ind_xy),], size=12, col=pal_axes[3])
for (j in 1:13) {
  if (j %% 2 == 1) {
    points3d(med_xyz[ind_xy[ind_xy[,2] == j,1], ], size=12, col=pal_axes[2])
    points3d(med_xyz[ind_xy[ind_xy[,2] == -j,1], ], size=12, col=pal_axes[2])}
}

for (j in 1:nrow(RF_med)) {
  ii <- regexpr("T4", anno[j,"name"]) #start index
  type_ii <- match(substring(anno[j,"name"], ii+2, ii+2), letters)
  arrow3d(RF_med[j,4:6], RF_med[j,1:3], theta = pi / 12, n = 4, col = T4_col[type_ii], alpha = 0.9, type = "rotation", lit=F)
}
title3d(main = "Strahler")
points3d(med_xyz[match(ind_Up_lens, eyemap[,2]),], col = "gold2", size = 15)
points3d(med_xyz[match(ind_Down_lens, eyemap[,2]),], col = "cyan", size = 15)

# quantify alignment ------------------------------------------------------

# - PLOT 3D arrows
# center
n0 <- 1
n1 <- dim(anno_overlap_M)[1]
vM <- RF_lens[n0:n1,]
# vM <- RF_med[n0:n1,]
nopen3d()
for (j in n0:n1) {
  ii <- regexpr("T4", anno_overlap_M[j,"name"]) #start index
  type_ii <- match(substring(anno_overlap_M[j,"name"], ii+2, ii+2), letters)
  # arrow3d(c(0,0,0), vM[j,4:6]-vM[j,1:3], theta = pi / 18, n = 4, col = type_ii, type = "rotation")
  # # pch3d(matrix(RF_med[j,4:6], ncol = 3), col = "red", cex = 2, pch = 96+type_ii)
  # pch3d(matrix(vM[j,4:6]-vM[j,1:3], ncol = 3), col = "blue", cex = 2, pch = 96+type_ii)
  arrow3d(c(0,0,0), vM[j,4:6]-vM[j,1:3], theta = pi / 18, n = 4, col = T4_col[type_ii], type = "rotation")
  pch3d(matrix(vM[j,4:6]-vM[j,1:3], ncol = 3), col = "black", cex = 2, pch = 96+type_ii)
}

# dorsal
n0 <- dim(anno_overlap_M)[1]+1
n1 <- dim(anno_overlap_M)[1]+dim(anno_overlap_LD)[1]
vD <- RF_lens[n0:n1,]
nopen3d()
for (j in 1:dim(vD)[1]) {
  ii <- regexpr("T4", anno_overlap_LD[j,"name"]) #start index
  type_ii <- match(substring(anno_overlap_LD[j,"name"], ii+2, ii+2), letters)
  arrow3d(c(0,0,0), vD[j,4:6]-vD[j,1:3], theta = pi / 18, n = 4, col = 'red', type = "rotation")
  pch3d(matrix(vD[j,4:6]-vD[j,1:3], ncol = 3), col = "black", cex = 2, pch = 96+type_ii)
}

# ventral
n0 <- dim(anno_overlap_M)[1]+dim(anno_overlap_LD)[1]+1
n1 <- dim(anno_overlap_M)[1]+dim(anno_overlap_LD)[1]+dim(anno_overlap_LV)[1]
vV <- RF_lens[n0:n1,]
nopen3d()
for (j in 1:dim(vV)[1]) {
  ii <- regexpr("T4", anno_overlap_LV[j,"name"]) #start index
  type_ii <- match(substring(anno_overlap_LV[j,"name"], ii+2, ii+2), letters)
  arrow3d(c(0,0,0), vV[j,4:6]-vV[j,1:3], theta = pi / 18, n = 4, col = T4_col[type_ii], type = "rotation")
  pch3d(matrix(vV[j,4:6]-vV[j,1:3], ncol = 3), col = "black", cex = 2, pch = 96+type_ii)
}

# ventral posterior
n0 <- dim(anno_overlap_M)[1]+dim(anno_overlap_LD)[1]+dim(anno_overlap_LV)[1]+1
n1 <- dim(anno_overlap)[1]
vVP <- RF_lens[n0:n1,]
nopen3d()
for (j in 1:dim(vVP)[1]) {
  ii <- regexpr("T4", anno_overlap_VP[j,"name"]) #start index
  type_ii <- match(substring(anno_overlap_VP[j,"name"], ii+2, ii+2), letters)
  arrow3d(c(0,0,0), vVP[j,4:6]-vVP[j,1:3], theta = pi / 18, n = 4, col = T4_col[type_ii], type = "rotation")
  pch3d(matrix(vVP[j,4:6]-vVP[j,1:3], ncol = 3), col = "black", cex = 2, pch = 96+type_ii)
}


# - stats
vS <- vM[,4:6]-vM[,1:3]
vpc <- prcomp(vS)
vS <- vpc$x[,1:2]
dev.new()
plot(vpc$x[,1:2])
ii <- regexpr("T4", anno_overlap_M[j,"name"])
text(vpc$x[,1:2], labels =  match(substring(anno_overlap_M[,"name"], ii+2, ii+2), letters) )
angS <- acos(vpc$x[,1]/sqrt(rowSums(vpc$x[,1:2]^2)))*(-1)^(vpc$x[,2]<0) + 2*pi*(vpc$x[,2]<0)
angS <- angS/pi*180
anno_stats <- anno_overlap_M
pair_ang <- matrix(ncol = 3, nrow = choose(dim(vS)[1],2))
rr <- 1
for (j in 1:(dim(vS)[1]-1)) {
  for (k in (j+1):dim(vS)[1]) {
    jj <- regexpr("T4", anno_stats[j,"name"]) #start index
    kk <- regexpr("T4", anno_stats[k,"name"]) #start index
    type_jj <- match(substring(anno_stats[j,"name"], jj+2, jj+2), letters)
    type_kk <- match(substring(anno_stats[k,"name"], kk+2, kk+2), letters)
    if ((type_jj + type_kk ==3) | (type_jj + type_kk == 7)) {
      pair_ang[rr,1] <- abs(angS[j] - angS[k])
    } else {
      pair_ang[rr,1] <- acos(vpc$x[j,1:2] %*% vpc$x[k,1:2] / sqrt(sum(vpc$x[j,1:2]^2)) / sqrt(sum(vpc$x[k,1:2]^2))) /pi*180
    }
    pair_ang[rr,2] <- type_jj
    pair_ang[rr,3] <- type_kk
    rr <- rr+1
  }
}
pair_ang_M <- pair_ang


# pair_ang <- matrix(ncol = 3, nrow = choose(dim(vO)[1],2))
# rr <- 1
# for (j in 1:(dim(vO)[1]-1)) {
#   for (k in (j+1):dim(vO)[1]) {
#     jj <- regexpr("T4", anno_stats[j,"name"]) #start index
#     kk <- regexpr("T4", anno_stats[k,"name"]) #start index
#     type_jj <- match(substring(anno_stats[j,"name"], jj+2, jj+2), letters)
#     type_kk <- match(substring(anno_stats[k,"name"], kk+2, kk+2), letters)
#     pair_ang[rr,1] <- acos(vO[j,] %*% vO[k,] / sqrt(sum(vO[j,]^2)) / sqrt(sum(vO[k,]^2))) / pi * 180
#     pair_ang[rr,2] <- type_jj
#     pair_ang[rr,3] <- type_kk
#     rr <- rr+1
#   }
# }
# pair_ang_M <- pair_ang

vS <- vV[,4:6]-vV[,1:3]
vpc <- prcomp(vS)
angS <- acos(vpc$x[,1]/sqrt(rowSums(vpc$x[,1:2]^2)))*(-1)^(vpc$x[,2]<0) + 2*pi*(vpc$x[,2]<0)
angS <- angS/pi*180
anno_stats <- anno_overlap_LV
pair_ang <- matrix(ncol = 3, nrow = choose(dim(vS)[1],2))
rr <- 1
for (j in 1:(dim(vS)[1]-1)) {
  for (k in (j+1):dim(vS)[1]) {
    jj <- regexpr("T4", anno_stats[j,"name"]) #start index
    kk <- regexpr("T4", anno_stats[k,"name"]) #start index
    type_jj <- match(substring(anno_stats[j,"name"], jj+2, jj+2), letters)
    type_kk <- match(substring(anno_stats[k,"name"], kk+2, kk+2), letters)
    if ((type_jj + type_kk ==3) | (type_jj + type_kk == 7)) {
      pair_ang[rr,1] <- abs(angS[j] - angS[k])
    } else {
      pair_ang[rr,1] <- acos(vpc$x[j,1:2] %*% vpc$x[k,1:2] / sqrt(sum(vpc$x[j,1:2]^2)) / sqrt(sum(vpc$x[k,1:2]^2))) /pi*180
    }
    pair_ang[rr,2] <- type_jj
    pair_ang[rr,3] <- type_kk
    rr <- rr+1
  }
}
pair_ang_V <- pair_ang

vS <- vVP[,4:6]-vVP[,1:3]
vpc <- prcomp(vS)
angS <- acos(vpc$x[,1]/sqrt(rowSums(vpc$x[,1:2]^2)))*(-1)^(vpc$x[,2]<0) + 2*pi*(vpc$x[,2]<0)
angS <- angS/pi*180
anno_stats <- anno_overlap_VP
pair_ang <- matrix(ncol = 3, nrow = choose(dim(vS)[1],2))
rr <- 1
for (j in 1:(dim(vS)[1]-1)) {
  for (k in (j+1):dim(vS)[1]) {
    jj <- regexpr("T4", anno_stats[j,"name"]) #start index
    kk <- regexpr("T4", anno_stats[k,"name"]) #start index
    type_jj <- match(substring(anno_stats[j,"name"], jj+2, jj+2), letters)
    type_kk <- match(substring(anno_stats[k,"name"], kk+2, kk+2), letters)
    if ((type_jj + type_kk ==3) | (type_jj + type_kk == 7)) {
      pair_ang[rr,1] <- abs(angS[j] - angS[k])
    } else {
      pair_ang[rr,1] <- acos(vpc$x[j,1:2] %*% vpc$x[k,1:2] / sqrt(sum(vpc$x[j,1:2]^2)) / sqrt(sum(vpc$x[k,1:2]^2))) /pi*180
    }
    pair_ang[rr,2] <- type_jj
    pair_ang[rr,3] <- type_kk
    rr <- rr+1
  }
}
pair_ang_VP <- pair_ang

pair_ang <- rbind(pair_ang_M, pair_ang_V)
pair_ang <- rbind(pair_ang, pair_ang_VP)
pair_ang <- as.data.frame(pair_ang)
colnames(pair_ang) <- c('ang', 't1', 't2')
pair_ang["pos"] <- c(rep("M",dim(pair_ang_M)[1]), rep("VA",dim(pair_ang_V)[1]), rep("VP",dim(pair_ang_VP)[1]))

pair_ang_2 <- pair_ang
pair_ang_2["pair"] <- NA

pair_ang_2[pair_ang_2["t1"] == 1 & pair_ang_2["t2"] == 2, "pair"] <- "a-b"
pair_ang_2[pair_ang_2["t1"] == 3 & pair_ang_2["t2"] == 4, "pair"] <- "c-d"

pair_ang_2[pair_ang_2["t1"] == 1 & pair_ang_2["t2"] == 4, "pair"] <- "a-d"
pair_ang_2[pair_ang_2["t1"] == 2 & pair_ang_2["t2"] == 3, "pair"] <- "b-c"
pair_ang_2[pair_ang_2["t1"] == 1 & pair_ang_2["t2"] == 3, "pair"] <- "a-c"
pair_ang_2[pair_ang_2["t1"] == 2 & pair_ang_2["t2"] == 4, "pair"] <- "b-d"

pair_ang_2 <- pair_ang_2[!is.na(pair_ang_2["pair"]),]
pair_ang_2["pair"] <- factor(pair_ang_2$pair, levels = c("a-b", "c-d", "a-d","b-c","a-c","b-d"), ordered = T)
pair_ang_2["pos"] <- factor(pair_ang_2$pos, levels = c("VA", "M", "VP"), ordered = T)

pair_ang_ab <- pair_ang_2[pair_ang_2$pair %in% c("a-b", "c-d"), ]
pair_ang_ac <- pair_ang_2[pair_ang_2$pair %in% c("a-c", "b-d"), ]
pair_ang_ad <- pair_ang_2[pair_ang_2$pair %in% c("a-d", "c-d"), ]


# t1 <- 1
# t2 <- 4
# pair_ang_bool <- (pair_ang[,2]==t1 & pair_ang[,3]==t2) | (pair_ang[,2]==t2 & pair_ang[,3] == t1)
# ang_bd <- pair_ang[pair_ang_bool, 1] #b and d
# mean(ang_bd)
# sd(ang_bd)
# # pair_ang[,'t1'] == 1 & pair_ang[,'t2'] == 2


# plot ggplot
windows(record = F, width = 8, height = 8)
ggplot(pair_ang_2, aes(pair, ang)) +
  geom_boxplot(aes(fill = pos)) +
  scale_fill_manual(values = c("royalblue", "seagreen", "orange")) +
  # geom_point(aes(colour = pos)) +
  theme_bw() +
  geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 180), size = 1, linetype = 2, colour = "brown") +
  # geom_segment(aes(x = 4.5, y = 0, xend = 4.5, yend = 180), size = 1, linetype = 2, colour = "gold") +
  geom_segment(aes(x = 0, y = 90, xend = 7, yend = 90), size = 1, linetype = 2, colour = "gold2") 
  # geom_segment(aes(x = 0, y = 180, xend = 7, yend = 180), size = 1, linetype = 2, colour = "brown") 

# seperate plots, HHMI meeting
windows(record = F, width = 8, height = 8)
ggplot(pair_ang_ab, aes(pair, ang)) +
  geom_boxplot(aes(fill = pos)) +
  scale_fill_manual(values = c("royalblue", "seagreen", "orange")) +
  # geom_point(aes(colour = pos)) +
  theme_bw() +
  xlab("T4 pair") +
  ylab("Angle [degree]") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12,face="bold"),
        legend.title = element_text(size=14,face="bold")) +
  # geom_segment(aes(x = 2.5, y = 0, xend = 2.5, yend = 180), size = 1, linetype = 2, colour = "brown") +
  # geom_segment(aes(x = 4.5, y = 0, xend = 4.5, yend = 180), size = 1, linetype = 2, colour = "gold") +
  geom_segment(aes(x = 0.5, y = 180, xend = 2.5, yend = 180), size = 1, linetype = 2, colour = "brown") 


windows(record = F, width = 8, height = 8)
ggplot(pair_ang_ac, aes(pair, ang)) +
  geom_boxplot(aes(fill = pos)) +
  scale_fill_manual(values = c("royalblue", "seagreen", "orange")) +
  theme_bw() +
  xlab("T4 pair") +
  ylab("Angle [degree]") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12,face="bold"),
        legend.title = element_text(size=14,face="bold")) +
  geom_segment(aes(x = 0.5, y = 90, xend = 2.5, yend = 90), size = 1, linetype = 2, colour = "brown") 

windows(record = F, width = 8, height = 8)
ggplot(pair_ang_ad, aes(pair, ang)) +
  geom_boxplot(aes(fill = pos)) +
  scale_fill_manual(values = c("royalblue", "seagreen", "orange")) +
  theme_bw() +
  xlab("T4 pair") +
  ylab("Angle [degree]") +
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12,face="bold"),
        legend.title = element_text(size=14,face="bold")) +
  geom_segment(aes(x = 0.5, y = 90, xend = 2.5, yend = 90), size = 1, linetype = 2, colour = "brown") 

# plot
dev.new()
dx <- - 0.2 # =/- 0.2
cc <- 6 # 6 - pink, 4 - blue
t1 <- 1; t2 <- 2; xx <- 1+dx-1
pair_ang_bool <- (pair_ang[,2]==t1 & pair_ang[,3]==t2) | (pair_ang[,2]==t2 & pair_ang[,3] == t1)
ang_p <- pair_ang[pair_ang_bool, 1]

plot(rep(xx,length(ang_p)), ang_p, xlim = c(-1,7), ylim = c(0, 180), axes = F, xlab = 'T4 pair', ylab = 'Angle', col = cc)
points(rep(xx,length(ang_p)), ang_p, col = cc)
lines(c(xx-0.2, xx+0.2), c(mean(ang_p),mean(ang_p)), lwd = 3, col = cc)
t1 <- 3; t2 <- 4; xx <- 2+dx-1
pair_ang_bool <- (pair_ang[,2]==t1 & pair_ang[,3]==t2) | (pair_ang[,2]==t2 & pair_ang[,3] == t1)
ang_p <- pair_ang[pair_ang_bool, 1]
points(rep(xx,length(ang_p)), ang_p, col = cc)
lines(c(xx-0.2, xx+0.2), c(mean(ang_p),mean(ang_p)), lwd = 3, col = cc)
t1 <- 1; t2 <- 4; xx <- 3+dx
pair_ang_bool <- (pair_ang[,2]==t1 & pair_ang[,3]==t2) | (pair_ang[,2]==t2 & pair_ang[,3] == t1)
ang_p <- pair_ang[pair_ang_bool, 1]
points(rep(xx,length(ang_p)), ang_p, col = cc)
lines(c(xx-0.2, xx+0.2), c(mean(ang_p),mean(ang_p)), lwd = 3, col = cc)
t1 <- 2; t2 <- 3; xx <- 4+dx
pair_ang_bool <- (pair_ang[,2]==t1 & pair_ang[,3]==t2) | (pair_ang[,2]==t2 & pair_ang[,3] == t1)
ang_p <- pair_ang[pair_ang_bool, 1]
points(rep(xx,length(ang_p)), ang_p, col = cc)
lines(c(xx-0.2, xx+0.2), c(mean(ang_p),mean(ang_p)), lwd = 3, col = cc)
t1 <- 1; t2 <- 3; xx <- 5+dx
pair_ang_bool <- (pair_ang[,2]==t1 & pair_ang[,3]==t2) | (pair_ang[,2]==t2 & pair_ang[,3] == t1)
ang_p <- pair_ang[pair_ang_bool, 1]
points(rep(xx,length(ang_p)), ang_p, col = cc)
lines(c(xx-0.2, xx+0.2), c(mean(ang_p),mean(ang_p)), lwd = 3, col = cc)
t1 <- 2; t2 <- 4; xx <- 6+dx
pair_ang_bool <- (pair_ang[,2]==t1 & pair_ang[,3]==t2) | (pair_ang[,2]==t2 & pair_ang[,3] == t1)
ang_p <- pair_ang[pair_ang_bool, 1]
points(rep(xx,length(ang_p)), ang_p, col = cc)
lines(c(xx-0.2, xx+0.2), c(mean(ang_p),mean(ang_p)), lwd = 3, col = cc)

axis(side = 1, at = c(0,1,3,4,5,6), labels = c('a-b', 'c-d', 'a-d', 'b-c', 'a-c', 'b-d'))
axis(side = 2, at = seq(0,180,by = 30))
lines(c(2,2), c(0,180), lty = 5, lwd = 3, col = "black")
lines(c(-1,7), c(90,90), lty = 3, lwd = 3, col = "green")
lines(c(-1,7), c(60,60), lty = 4, lwd = 3, col = "green")


# - TALK connectome berlin
# nopen3d()
# tar <- T4_overlap_M[[2]] #a
# plot3d(tar, col = 'red', lwd = 1, alpha = 0.5) #a
# targ <- as.ngraph(tar)
# ii_root <- match(tar$tags$`T-bars start`, tar$d[,'PointNo']) # use a depth first search
# sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
# sub_tree <- subset(tar, sub_points)
# plot3d(sub_tree, col = 'red', lwd = 3) #a
# arrow3d(RF_med[2,4:6], RF_med[2,1:3],  theta = pi / 8, n = 4, col = 'red', type = "rotation")
# 
# tar <- T4_overlap_M[[4]] #b
# plot3d(tar, col = 'green', lwd = 1, alpha = 0.5) #b
# targ <- as.ngraph(tar)
# ii_root <- match(tar$tags$`T-bars start`, tar$d[,'PointNo']) # use a depth first search
# sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
# sub_tree <- subset(tar, sub_points)
# plot3d(sub_tree, col = 'green', lwd = 3) 
# arrow3d( RF_med[4,4:6], RF_med[4,1:3],theta = pi / 8, n = 4, col = 'green', type = "rotation")
# 
# 
# tar <- T4_overlap_M[[5]] #c
# plot3d(tar, col = 'brown', lwd = 1, alpha = 0.5) 
# targ <- as.ngraph(tar)
# ii_root <- match(tar$tags$`T-bars start`, tar$d[,'PointNo']) # use a depth first search
# sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
# sub_tree <- subset(tar, sub_points)
# plot3d(sub_tree, col = 'brown', lwd = 3) 
# arrow3d( RF_med[5,4:6], RF_med[5,1:3],theta = pi / 8, n = 4, col = 'brown', type = "rotation")
# 
# tar <- T4_overlap_M[[3]] #d
# plot3d(tar, col = 'blue', lwd = 1, alpha = 0.5) 
# targ <- as.ngraph(tar)
# ii_root <- match(tar$tags$`T-bars start`, tar$d[,'PointNo']) # use a depth first search
# sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
# sub_tree <- subset(tar, sub_points)
# plot3d(sub_tree, col = 'blue', lwd = 3) 
# arrow3d( RF_med[3,4:6], RF_med[3,1:3],theta = pi / 8, n = 4, col = 'blue', type = "rotation")
# 
# # T4a b awith arrows
# nopen3d()
# tar <- T4_overlap_M[[2]] #a
# plot3d(tar, col = 'red', lwd = 1, alpha = 0.5) #a
# targ <- as.ngraph(tar)
# ii_root <- match(tar$tags$`dendrite start`, tar$d[,'PointNo']) # use a depth first search
# sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
# sub_tree <- subset(tar, sub_points)
# plot3d(sub_tree, col = 'red', lwd = 3, alpha = 0.5) #a
# arrow3d(RF_med[2,4:6], RF_med[2,1:3],  theta = pi / 8, n = 4, col = 'red', type = "rotation")
# 
# tar <- T4_overlap_M[[4]] #b
# plot3d(tar, col = 'green', lwd = 1, alpha = 0.5) #b
# targ <- as.ngraph(tar)
# ii_root <- match(tar$tags$`dendrite start`, tar$d[,'PointNo']) # use a depth first search
# sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
# sub_tree <- subset(tar, sub_points)
# plot3d(sub_tree, col = 'green', lwd = 3, alpha = 0.5) #a
# arrow3d( RF_med[4,4:6], RF_med[4,1:3],theta = pi / 8, n = 4, col = 'green', type = "rotation")
# 
# # T4c d with arrows
# nopen3d()
# tar <- T4_overlap_M[[5]] #c
# plot3d(tar, col = 'black', lwd = 1, alpha = 0.5) 
# targ <- as.ngraph(tar)
# ii_root <- match(tar$tags$`dendrite start`, tar$d[,'PointNo']) # use a depth first search
# sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
# sub_tree <- subset(tar, sub_points)
# # plot3d(sub_tree, col = 'black', lwd = 3) #a
# arrow3d( RF_med[5,4:6], RF_med[5,1:3],theta = pi / 8, n = 4, col = 'black', type = "rotation")
# 
# tar <- T4_overlap_M[[3]] #d
# plot3d(tar, col = 'blue', lwd = 1, alpha = 0.5) 
# targ <- as.ngraph(tar)
# ii_root <- match(tar$tags$`dendrite start`, tar$d[,'PointNo']) # use a depth first search
# sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
# sub_tree <- subset(tar, sub_points)
# # plot3d(sub_tree, col = 'blue', lwd = 3) 
# arrow3d( RF_med[3,4:6], RF_med[3,1:3],theta = pi / 8, n = 4, col = 'blue', type = "rotation")




# add T4 ------------------------------------------------------------------

# status anno
LL <- 4
anno_str_axon <- paste("T4",letters[LL], " - axon", sep = "")
anno_str_dend <- paste("T4",letters[LL], " - dendrites", sep = "")
anno_str_comp <- paste("T4",letters[LL], " - complete", sep = "")
anno_axon <- catmaid_query_by_annotation(anno_str_axon)
anno_dend <- catmaid_query_by_annotation(anno_str_dend)
anno_comp <- catmaid_query_by_annotation(anno_str_comp)


# -- calculated RF
RF_lens <- RF_lens_T4d
RF_med <- RF_med_T4d

# # PLOT RF on lens
nopen3d()
points3d(ucl_rot, col = "grey", size = 5)
for (j in 1:dim(na.omit(RF_lens))[1]) {
  arrow3d(RF_lens[j,1:3], RF_lens[j,1:3]+(RF_lens[j,4:6]-RF_lens[j,1:3])*1, theta = pi / 4, n = 4, col = "red", type = "rotation")
  text3d(RF_lens[j,1:3], texts = paste(j), adj = -1)
}
planes3d(0,0,1, 0, alpha = 0.1);planes3d(1,0,0, 0, alpha = 0.1);planes3d(0,1,0, 0, alpha = 0.1)

# # PLOT, RF on CT1
nopen3d()
points3d(med_xyz, col = "grey", size = 5)
for (j in 1:dim(na.omit(RF_med))[1]) {
  arrow3d(RF_med[j,1:3], RF_med[j,1:3]+(RF_med[j,4:6]-RF_med[j,1:3])*1, theta = pi / 4, n = 4, col = "red", type = "rotation")
  text3d(RF_med[j,1:3], texts = paste(j), adj = -1)
}

# -- predicated ori
m <- 4
nopen3d()
points3d(ucl_rot, col = "grey", size = 5)
for (j in 1:dim(RF_lens_T4_pred)[1]) {
  arrow3d(RF_lens_T4_pred[j,(3*(m-1)+1):(3*(m-1)+3)], ucl_rot[j,], theta = pi / 4, n = 4, col = "red", type = "rotation")
}
planes3d(0,0,1, 0, alpha = 0.1)
planes3d(1,0,0, 0, alpha = 0.1)
planes3d(0,1,0, 0, alpha = 0.1)


# medulla to LOP, for LPT, map area to area,  ---------------------------------------------------------------------

# eyemap
load('data/eyemap.RData')

# lob layer
# load("data/LayerData_T4_20201212_7e7_1um.RData")
load('data/LayerData_T4_20210729_7e7_1um.RData')

layer_ls <- xyz_layer_T4

# combine grid points from 4 layers
len_layer <- c() #length of each grid
II <- 0
for (j in 1:4) {
  len_layer[j] <- dim(layer_ls[[j]])[1]
  row.names(layer_ls[[j]]) <- seq(II+1, II+len_layer[j]) # change names to canonical ordering, 1,2,3,etc
  II <- II+len_layer[j]
}
xyz_4layer <- do.call(rbind, layer_ls)

# load T4
load("data/neu_T4b.RData")

# LL <- 2
# anno_str <- paste("T4",letters[LL], " - complete", sep = "")
# anno_T4 <- catmaid_query_by_annotation(anno_str)
# neu_skid_comp <- anno_T4[,"skid"]
# neu_T4_comp <-  read.neurons.catmaid(neu_skid_comp, .progress = 'text')
# T4b <- neu_T4_comp
# anno_T4b <- anno_T4
# 
# # save(anno_T4b, T4b, file = "data/neu_T4b.RData")


# calculate T4 com
T4b_com <- matrix(nrow = length(T4b), ncol = 6) #com, [dendrite, axon]
for (j in 1:length(T4b)) {
  tar <- T4b[[j]]
  ind_T = match(tar$tags$"T-bars start" , tar$d$PointNo)
  mat_T = tar$d[ind_T,]
  xyz_T <- xyzmatrix(mat_T)
  ind_D = match(tar$tags$"dendrite start" , tar$d$PointNo)
  mat_D = tar$d[ind_D,]
  xyz_D <- xyzmatrix(mat_D)
  
  dp_T <- child_node(tar, mat_T) #downstream nodes of mat_T
  dp_T <- xyzmatrix(dp_T)
  dp_D <- child_node(tar, mat_D) 
  dp_D <- xyzmatrix(dp_D)
  
  T4b_com[j,] <- c(colMeans(dp_D), colMeans(dp_T))
}

# - BETTER
# map T4 onto L2 grid
T4b_L2_xyz <- matrix(ncol = 3, nrow = length(T4b))
for (j in 1:length(T4b)) {
  ii_L2 <- which.min(rowSums(sweep(xyz_layer_T4[[2]], 2, T4b_com[j, 4:6])^2))
  T4b_L2_xyz[j,] <- unlist(xyz_layer_T4[[2]][ii_L2,])
}

# kernel regression global
np_eval <- med_xyz
lop_pred <- matrix(nrow = nrow(np_eval), ncol = 3)
# ind_nb <- sweep(T4b_com[,1:3],2,np_eval[j,],'-')^2 %>% rowSums() %>% order(decreasing=F) %>% head(30)
for (k in 1:3) {
  # npdata <- data.frame(mc = T4b_com[ind_nb,1:3], ec = T4b_L2_xyz[ind_nb,k])
  npdata <- data.frame(mc = T4b_com[,1:3], ec = T4b_L2_xyz[,k])
  # bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'fixed', regtype= 'll') # lc are bad
  bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'generalized_nn', regtype= 'll') # gnn+ll promising, need more T4?
  # bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
  model_np <- npreg(bw)
  # model_np <- npreg(ec ~ mc.1 + mc.2 + mc.3, data = npdata)
  for (j in 1:nrow(np_eval)) {
    # med col
    np_eval_one <- data.frame(mc.1 = np_eval[j,1], mc.2 = np_eval[j,2], mc.3 = np_eval[j,3])
    lop_pred[j,k] <- predict(model_np, newdata = np_eval_one)
  }
}
lop_pred_global <- lop_pred


# # kernel regression local
# lop_pred <- matrix(nrow = nrow(med_xyz), ncol = 3)
# for (j in 1:nrow(med_xyz)) {
#   np_eval <- med_xyz[j,]
#   ind_nb <- sweep(T4b_com[,1:3],2,np_eval,'-')^2 %>% rowSums() %>% order(decreasing=F) %>% head(1+8+16)
#   for (k in 1:3) {
#     npdata <- data.frame(mc = T4b_com[ind_nb,1:3], ec = T4b_com[ind_nb,3+k])
#     bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'fixed', regtype= 'll')
#     # bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'generalized_nn', regtype= 'll')
#     model_np <- npreg(bw)
#     np_eval_one <- data.frame(mc.1 = np_eval[1], mc.2 = np_eval[2], mc.3 = np_eval[3])
#     lop_pred[j,k] <- predict(model_np, newdata = np_eval_one)
# 
#     # # 1D
#     # npdata <- data.frame(mc = T4b_com[ind_nb,k], ec = T4b_com[ind_nb,3+k])
#     # model_np <- npreg(ec ~ mc, data = npdata)
#     # np_eval_one <- data.frame(mc = np_eval[k])
#     # lop_pred[j,k] <- predict(model_np, newdata = np_eval_one)
#   }
# }
# lop_pred_local <- lop_pred
# # lop_pred_local_1D <- lop_pred


# # loess, has to be local
# T4b_com_df <- as.data.frame(T4b_com)
# colnames(T4b_com_df) <- c('m1','m2','m3','l1','l2','l3')
# med_xyz_df <- as.data.frame(med_xyz)
# colnames(med_xyz_df) <- c('m1','m2','m3')
# 
# lop_pred_loess <- matrix(ncol = 3, nrow = nrow(med_xyz))
# for (j in 1:nrow(med_xyz)) {
#   ind_nb <- sweep(T4b_com_df[,1:3],2,med_xyz[j,],'-')^2 %>% rowSums() %>% order(decreasing=F) %>% head(20)
#   df <- T4b_com_df[ind_nb, ]
# 
#   model_loess <- loess(l1 ~ m1 * m2 * m3, data = df, degree = 2, span = 1,
#                        control = loess.control(surface = "direct"))
#   lop_pred_loess[j,1] <- predict(model_loess, med_xyz_df[j,], se = T)$fit
#   model_loess <- loess(l2 ~ m1 * m2 * m3, data = df, degree = 2, span = 1,
#                        control = loess.control(surface = "direct"))
#   lop_pred_loess[j,2] <- predict(model_loess, med_xyz_df[j,], se = T)$fit
#   model_loess <- loess(l3 ~ m1 * m2 * m3, data = df, degree = 2, span = 1,
#                        control = loess.control(surface = "direct"))
#   lop_pred_loess[j,3] <- predict(model_loess, med_xyz_df[j,], se = T)$fit
# }


# - combine pt in layer 2 and associate with ii_stack
L2_med <- cbind(lop_pred, med_xyz) #pt in LOP, mapped by T4 from med col
rownames(L2_med) <- rownames(med_xyz)

# # SAVE
# save(L2_med, T4b_com, file = "data/L2_med.RData") 


# check mapping 
# T4b med to lop
nopen3d()
points3d(T4b_com[, 1:3], col = "pink", size = 6)
points3d(T4b_com[, 4:6], col = "gray", size = 12)
linemat <- matrix(t(cbind(T4b_com[, 1:3], T4b_com[, 4:6])), ncol = 3, byrow = T)
segments3d(linemat, color = "grey")

# plot, med lop map
nopen3d()
points3d(med_xyz, col = 'red', size = 5)
# points3d(sweep(L2_med, 2, c(50000,0,50000), FUN = "+"), col = 'blue', size = 5)
points3d(sweep(L2_med[,1:3], 2, c(50000,0,50000), FUN = "+"), col = 'blue', size = 5)
linemat <- matrix(t(cbind(med_xyz,
                          sweep(L2_med[,1:3], 2, c(50000,0,50000), FUN = "+")
)), ncol = 3, byrow = T)
segments3d(linemat, color = "grey")

# plot, med lop map
nopen3d()
points3d(med_xyz, col = 'red', size = 5)
points3d(sweep(L2_med, 2, c(50000,0,50000), FUN = "+"), col = 'blue', size = 5)
for (j in 1:dim(med_xyz)[1]) {
  lines3d(rbind(sweep(L2_med, 2, c(50000,0,50000), FUN = "+")[j,], med_xyz[as.integer(rownames(L2_med)[j]),]), col = "grey")
}
plot3d(CT1, col = 'green', lwd = 1, alpha = 1)

# label pt with med rowname
nopen3d()
points3d(med_xyz, col = 'red', size = 5, alpha = 0.3)
for (j in 1:dim(med_xyz)[1]) {
  text3d(med_xyz[j,], texts = rownames(med_xyz)[j])
}
nopen3d()
points3d(L2_med, col = 'blue', size = 5, alpha = 0.3)
for (j in 1:dim(med_xyz)[1]) {
  text3d(L2_med[j,], texts = rownames(L2_med)[j])
}

# check mapping -----------------------------------------------------------

# CT1 vs Mi1 lop_pred
nopen3d()
points3d(med_xyz_Mi1, size = 7, col = 'black')
points3d(med_xyz, size = 7, col = 'gold2')
linemat <- matrix(t(cbind(med_xyz_Mi1, med_xyz)), ncol = 3, byrow = T)
segments3d(linemat, color = "grey")

nopen3d()
points3d(lop_pred, size = 5, col = 'black')
points3d(lop_pred_CT1, size = 5, col = 'cyan')
linemat <- matrix(t(cbind(lop_pred, lop_pred_CT1)), ncol = 3, byrow = T)
segments3d(linemat, color = "grey")

# - lens vs Mi1
nopen3d()
points3d(utp_Mi1_rot_chiasm, size = 5, col = 'black')
points3d(utp_lens_rot_shift, size = 5, col = 'blue')
linemat <- matrix(t(cbind(utp_Mi1_rot_chiasm[lens_Mi1[,2],],
                          utp_lens_rot_shift[lens_Mi1[,1],])), ncol = 3, byrow = T)
segments3d(linemat, color = "grey")

points3d(utp_lens_rot_shift[!(seq(1,nrow(utp_lens_rot_shift)) %in% lens_Mi1[,1]),], size = 15, col = 'cyan')
points3d(utp_Mi1_rot_chiasm[!(seq(1,nrow(utp_Mi1_rot_chiasm)) %in% lens_Mi1[,2]),], size = 15, col = 'magenta')

# 
rgl.set(47)
points3d(utp_lens_rot_shift[a[[3]],], size = 10, col = 'pink')
points3d(head(utp_lens_rot_shift[a0,]), size = 20, col = 'green')

rgl.set(46)
points3d(utp_Mi1_rot_chiasm[b[[3]],], size = 10, col = 'red')

points3d(utp_lens_rot_shift[am[[14]],], size = 10, col = 'pink')
points3d(utp_Mi1_rot_chiasm[bm[[14]],], size = 10, col = 'blue')


# - Mi1 vs CT1
nopen3d()
points3d(Mi1_xyz, col = 'blue', size = 5)
# points3d(add_Mi1, col = 'cyan', size = 15, alpha = 0.7) # manually added
points3d(CT1_xyz, col = 'black', size = 5)

linemat <- matrix(t(cbind(CT1_xyz, Mi1_xyz[CT1_Mi1,])), ncol = 3, byrow = T)
segments3d(linemat, color = "grey10")

# 
rgl.set(49)
points3d(Mi1_xyz[b[[6]],], size = 10, col = 'red')
points3d(head(Mi1_xyz[b0,]), size = 20, col = 'green')
# 
rgl.set(48)
points3d(CT1_xyz[match(b[[6]], CT1_Mi1),], size = 10, col = 'red')
points3d(head(CT1_xyz[match(b0, CT1_Mi1),]), size = 20, col = 'green')

# - CT1 med vs lop
pca <- prcomp(med_xyz)
if (pca$rotation[,3] %*% c(1,0,1) < 0) {
  pca$rotation <- - pca$rotation
}
med_xyz_pca <- sweep(med_xyz, 2, pca$center) %*% pca$rotation

pca <- prcomp(lop_pred)
if (pca$rotation[,3] %*% c(1,0,1) < 0) {
  pca$rotation <- - pca$rotation
}
lop_pred_pca <- sweep(lop_pred, 2, pca$center) %*% pca$rotation
lop_pred_pca <- sweep(lop_pred_pca, 2, c(0,0,60000)) 

nopen3d()
points3d(med_xyz_pca, size=7, col='blue')
points3d(lop_pred_pca, col = "gray20", size = 5)
linemat <- matrix(t(cbind(med_xyz_pca, lop_pred_pca)), ncol = 3, byrow = T)
segments3d(linemat, color = "grey")

ii <- eyemap[match(a[[12]], eyemap[,2]),1] %>% match(rownames(med_xyz_pca)) %>% na.omit()

# 
rgl.set(50)
points3d(lop_pred_pca[ii,], size = 12, col='red')
points3d(head(lop_pred_pca[ii,]), size = 20, col='green')

points3d(med_xyz_pca[ii,], size = 12, col='red')

# T4
nopen3d()
points3d(lop_pred, size = 5)
points3d(T4b_com[,4:6], size = 7, col='green')
points3d(lop_pred[ii,], size = 12, col='red')

nopen3d()
points3d(med_xyz, size = 5)
points3d(T4b_com[,1:3], size = 7, col='green')
points3d(med_xyz[ii,], size = 12, col='red')

rgl.pop()



# Mi1
nopen3d()
points3d(utp_Mi1_rot_chiasm, size = 5, col = 'black')
points3d(utp_Mi1_rot_chiasm[bb0,], size = 15, col = 1:length(bb0))

points3d(utp_lens_rot_shift, size = 5, col = 'red')


# TmY5a to LO --------------------------------------------------------------

load('data/neu_Mi1.RData')
load('data/eyemap.RData')
load('data/Mi1_ind.RData')

# - JFRC2010 mesh
load("data/JFRC2NP.surf.fafb.rda") # from Greg Jefferis, cf.p_R7R8
ME_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="ME_R")
LOP_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LOP_R")
LO_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LO_R")
shade3d(ME_msh, alpha=0.1)

# - TmY5a

# anno_TmY5a <- catmaid_query_by_annotation("^putative TmY5a$", type = "neuron")
# neu <-  read.neurons.catmaid(anno_TmY5a$skid, .progress = 'text')
# TmY5a <- kinkCT1(neu) # corr kink

tag_xyz <- matrix(ncol = 9, nrow = length(TmY5a))
for (j in 1:length(TmY5a)) {
  tar <- TmY5a[[j]]
  # tar <- resample(tar, stepsize = 100) #lose tags
  M10 = match(tar$tags$"M10", tar$d$PointNo)
  # LO4 = match(tar$tags$"LO4", tar$d$PointNo)
  LO4 = match(tar$tags$"LO_LC4", tar$d$PointNo)
  LoP2 = match(tar$tags$"LoP2", tar$d$PointNo)
  tag_xyz[j,] <- c(t(xyzmatrix(tar$d[c(M10, LO4, LoP2),])))
}

# matching Mi1 (for length(TmY5a)==63 only)
TmY5a_Mi1 <- c()
for (j in 1:length(TmY5a)) {
  TmY5a_Mi1 <- c(TmY5a_Mi1, sweep(med_xyz, 2, tag_xyz[j,1:3])^2 %>% rowSums() %>% which.min())
}
TmY5a_Mi1[23] <- 326 # manual corr
TmY5a_Mi1[15] <- 650

# - kernel regression global
np_eval <- med_xyz
lo_pred <- matrix(nrow = nrow(np_eval), ncol = 3)
for (k in 1:3) {
  npdata <- data.frame(mc = tag_xyz[,1:3], ec = tag_xyz[,3+k])
  # bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'fixed', regtype= 'll') # lc are bad
  bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'generalized_nn', regtype= 'll') # gnn+ll promising, need more T4?
  # bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
  model_np <- npreg(bw)
  # model_np <- npreg(ec ~ mc.1 + mc.2 + mc.3, data = npdata)
  for (j in 1:nrow(np_eval)) {
    # lo col
    np_eval_one <- data.frame(mc.1 = np_eval[j,1], mc.2 = np_eval[j,2], mc.3 = np_eval[j,3])
    lo_pred[j,k] <- predict(model_np, newdata = np_eval_one)
  }
}

nopen3d()
# plot3d(TmY5a_me, col='gray50')
# plot3d(TmY5a[c(10,12,13,14)], col='gray50')
points3d(tag_xyz[,1:3], col ='magenta', size = 10)
points3d(tag_xyz[,4:6], col ='cyan', size = 10)
points3d(tag_xyz[,7:9], col ='dark green', size = 10)

# points3d(Mi1_M10_xyz, size=10)
points3d(med_xyz, size=7)
# points3d(Mi1_M10_xyz[Mi1_ind_PR[[8]],], size = 15, col='blue')
# points3d(Mi1_M10_xyz[Mi1_ind_PR[[7]],], size = 15, col='light blue')

points3d(lo_pred, size=8, col='blue')

linemat <- matrix(t(cbind(med_xyz[TmY5a_Mi1,], tag_xyz[,1:3])), ncol = 3, byrow = T)
segments3d(linemat, color = "grey")

# SAVE , overwrite
# save(lo_pred, tag_xyz, TmY5a, anno_TmY5a, file = "data/neu_TmY5a.RData")


# DEBUG, add more TmY5a
# #check completion
# ind_tag <- matrix(ncol = 4, nrow = length(TmY5a))
# for (j in 1:length(TmY5a)) {
#   tar <- TmY5a[[j]]
#   # M9 = match(tar$tags$"M9", tar$d$PointNo)
#   M10 = match(tar$tags$"M10", tar$d$PointNo)
#   LO4 = match(tar$tags$"LO4", tar$d$PointNo)
#   LOLC4 = match(tar$tags$"LO_LC4", tar$d$PointNo)
#   LoP2 = match(tar$tags$"LoP2", tar$d$PointNo)
#   ind_tag[j,] <- c(length(M10), length(LO4), length(LOLC4), length(LoP2))
# }
# 
# TmY5a <- TmY5a[ind_tag[,3] == 1]


# M10_msh <- catmaid_get_volume("v14.M10_thin", rval = 'mesh3d') %>% as.mesh3d()
# 
# tar <- TmY5a[[20]]
# 
# nopen3d()
# plot3d(tar)
# shade3d(M10_msh, alpha=0.1)
# 
# ii <- identify3d(xyzmatrix(tar$d))
# points3d(xyzmatrix(tar$d[ii,]), size = 20)
# 
# catmaid_get_labels(treenodes = tar$d[ii, 'PointNo'])
# catmaid_set_labels(node = tar$d[ii, 'PointNo'], labels = 'M10')
#
# TmY5a_me <- nlapply(TmY5a, function(x) subset(x, pointsinside(x, ME_msh, rval="distance") > - 2000) ) #me portion 
# TmY5a_lp <- nlapply(TmY5a, function(x) subset(x, pointsinside(x, LOP_msh, rval="distance") > 2000) ) 
# TmY5a_lo <- nlapply(TmY5a, function(x) subset(x, pointsinside(x, LO_msh, rval="distance") > 2000) ) 
# 
# 
# # LC14
# skid <- catmaid_query_by_name("Putative LC14", type = "neuron")$skid
# neu <-  read.neurons.catmaid(skid, .progress = 'text')
# LC14 <- kinkCT1(neu) # corr kink
# 
# skid <- catmaid_query_by_name("Putative LC14a", type = "neuron")$skid
# neu <-  read.neurons.catmaid(skid, .progress = 'text')
# LC14a <- kinkCT1(neu) # corr kink
# 
# # # - make ME layer 10
# # M10 <- ashape3d(Mi1_M10_xyz, alpha = 50000) %>% as.mesh3d()
# # catmaid_add_volume(M10, title="v14.M10_thin", comment="by AZ")
# 
# # - add more
# nopen3d()
# # plot3d(TmY5a_me, col='gray50')
# # plot3d(TmY5a[c(10,12,13,14)], col='gray50')
# points3d(tag_xyz[,1:3], col ='magenta', size = 10)
# points3d(tag_xyz[,4:6], col ='cyan', size = 10)
# # points3d(tag_xyz[,7:9], col ='dark green', size = 10)
# 
# points3d(Mi1_M10_xyz)
# shade3d(M10_msh, alpha=0.2)
# points3d(Mi1_M10_xyz[Mi1_ind_PR[[8]],], size = 12, col='blue')
# points3d(Mi1_M10_xyz[Mi1_ind_PR[[7]],], size = 12, col='light blue')
# 
# points3d(Mi1_M10_xyz[xb_ls[[1]][[1]],], col='green',size =15)
# points3d(Mi1_M10_xyz[xb_ls[[2]][[1]],], col='green',size =15)
# points3d(Mi1_M10_xyz[ybm[[1]],], col='cyan',size =15)
# 
# plot3d(LC14$"9079380", col='blue', lwd=3)
# plot3d(LC14$"15780602", col='gray30', lwd=3)
# # SAVE
# # htmlwidgets::saveWidget(rglwidget(width = 1200, height = 800), "LC14.html")
# 
# points3d(med_xyz[match(vaxis_gen(-1), eyemap[,2]),], col=pal_axes[3],size =15)
# points3d(med_xyz[match(haxis_gen(0), eyemap[,2]),], col=pal_axes[4],size =15)
# 
# # - edit lo tags
# # load('data/layer_LC6.RData')
# # surface3d(xx,yy,m_surf, color = "#fdae61", alpha = 0.8)
# 
# load('data/layer_LC4.RData')
# # shade3d(msh, alpha=.2)
# points3d(xyz_layer, color = "#fdae61", size = 2)
# 
# # # -- add tags
# # ii <- identify3d(tag_xyz[,4:6])
# # tar <- TmY5a[[ii]]
# # plot3d(tar, lwd=4, col='black')
# # 
# # ii <- identify3d(xyzmatrix(tar$d))
# # points3d(xyzmatrix(tar$d[ii,]), size = 15, col='blue')
# # 
# # catmaid_get_labels(treenodes = tar$d[ii, 'PointNo'])
# # # catmaid_remove_labels(node = tar$d[ii, 'PointNo'], labels = 'Mi1 column marker')
# # catmaid_set_labels(node = tar$d[ii, 'PointNo'], labels = 'LO_LC4')



# # map one LPTC onto eye --------------------------------------------------
# 
# a = -3; b = 2; c = 1; d = 380000 # cut out LP portion for all left LP neurons, for VH/SH
# neu_skid = c(1111992, 1088678) # layer 2 LPTC and H@,
# LPTC2 <-  read.neurons.catmaid(neu_skid, .progress='text')
# 
# tar <- LPTC2[[1]]
# # extract LPTC vertices
# df_ep = tar$d[tar$EndPoints, ]
# xyz_ep = xyzmatrix(df_ep)
# xyz_ep_LP <-  xyz_ep[a * xyz_ep[, 1] + b * xyz_ep[, 2] + c * xyz_ep[, 3] + d > 0,] # for VS/HS, '>' for LOP
# # xyz_ep_LP <- xyz_ep[rowSums((sweep(xyz_ep, 2, cc)) ^ 2) < rr ^ 2, ] # for T4, use '<' for LOP
# # ind_ax <-  tail(tar$SegList[[1]], 1) # a special bp on the axon
# # xyz_ax <- tar$d[ind_ax, 3:5]
# 
# xyz_dend_pj <- xyz_ep_LP
# 
# # project dendrite pts to the fitted grid by shortest distance
# ind_min <- c()
# # index in xyz_4layer with min distance
# ind_min <- apply(xyz_dend_pj, 1, function(pt) {which.min(rowSums(sweep(lop_pred, 2, pt) ^ 2))})
# # ind_min <- apply(xyz_dend, 1, function(pt){which.min(rowSums(sweep(xyz_4layer,2,pt)^2))}) # index in xyz_4layer with min distance
# ind_min <- unique(ind_min) #index of grid points
# ind_pj <- rownames(np_eval)[ind_min]
# 
# # PLOT, RF on lens
# nopen3d()
# points3d(utp_ends_rot, col = "grey", size = 5)
# points3d(matrix(utp_ends_rot[ii_stack[ii_stack[as.numeric(ind_pj),2],3],], ncol = 3), col = "red", size = 8)
# axes3d(c('x','y','z')); title3d('','','x','y','z')
# 
# 
# # reconstruct (x,y) projection based on (x,y)
# ii <- sort(as.integer(ind_pj))
# xy_tmp <- list()
# xy_ashape <- list()
# xy_edge <- vector("list", 4)
# II <- 0
# for (k in 1:4) {
#   xy_tmp[[k]] <- xyz_4layer[ii[ii>II+1 & ii<II+len_layer[k]], 1:2]
#   if (dim(xy_tmp[[k]]) >= 3) {
#     xy_ashape[[k]] <- ashape(xy_tmp[[k]], alpha = 20000)
#     xy_edge[[k]] <- xy_ashape[[k]]$edges[,1:6]
# 
#   }
#   II <- II+len_layer[k]
# }
# 
# # projection of grid on to xy plane
# xy_edge_grid <- list()
# xy_ashape_grid <- list()
# for (j in 1:4) {
#   xy_ashape_grid[[j]] <- ashape(layer_ls[[j]][,1:2], alpha = 10000)
#   xy_edge_grid[[j]] <- xy_ashape_grid[[j]]$edge[,3:6]
# }
# 
# # DEBUG
# dev.new()
# x_limit <- c(400000, 250000)
# y_limit <- c(350000, 180000)
# # plot(xy_edge[[k]][,1:2],type = "n")
# # segments(c(xy_edge[[k]][,1]), c(xy_edge[[k]][,2]), c(xy_edge[[k]][,3]), c(xy_edge[[k]][,4]))
# # j=39
# k=2
# plot(xy_edge_grid[[k]][,1:2],type = "n", xlim = x_limit, ylim = y_limit, xlab = "", ylab = "", main = paste("Layer", as.character(k), sep=" "))
# segments(c(xy_edge_grid[[k]][,1]), c(xy_edge_grid[[k]][,2]), c(xy_edge_grid[[k]][,3]), c(xy_edge_grid[[k]][,4]))
# # plot(xy_edge[[k]][,1:2],type = "p")
# # segments(c(xy_edge[[k]][,1]), c(xy_edge[[k]][,2]), c(xy_edge[[k]][,3]), c(xy_edge[[k]][,4]))
# points(xy_tmp[[k]], col = 'black',  cex = 1,  pch = 1)
# 
# 
# ### --- PLOT, 2D ---
# # assign utp_ends to 2D disk
# utp_ends_disk <- utp_ends_rot[,2:3]
# ptcir <- seq(0,2*pi,by = 0.1)
# ptcir_xy <- data.frame(x=cos(ptcir), y=sin(ptcir))
# 
# dev.new()
# plot(ptcir_xy,type = "l",xlab = "", ylab = "")
# points(utp_ends_disk, col = 'black',  pch = 1, cex = 0.2)
# 
# rfcol <- utp_ends_rot[ii_stack[ii_stack[as.numeric(ind_pj),2],3],2:3]
# colnames(rfcol) <- c('x','y')
# arrsz <- 0.05
# for (j in 1:dim(rfcol)[1]) {
#   arrows(rfcol[j,1], rfcol[j,2], rfcol[j,1]+arrsz, rfcol[j,2], length = 0.05, col = 2, angle = 30, lwd = 2)
# }


# symbols()




