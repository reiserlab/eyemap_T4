# T4 morphology and PDs

library(natverse)
library(tidyverse)
library(igraph)
library(RColorBrewer)
library(sf) #intersection
library(deldir)
library(np) 

# clean everything up.
rm(list=ls())

source("eyemap_func.R")

# load data ---------------------------------------------------------------

load("data/eyemap.RData")
load('data/hexnb_ind_dist.RData') 

# load neurons -------------------------------------------------------------------------------------------------

load("data/neu_T4_dend.RData")

# # query from VFB
# T4_dend <- list()
# anno_T4_dend <- list()
# for (LL in 1:4) {
#   anno_str_dend <- paste("T4",letters[LL], " - dendrites", sep = "")
#   anno_dend <- catmaid_query_by_annotation(anno_str_dend)
#   anno_T4_dend[[LL]] <- anno_dend
#   neu_skid_dend <- anno_dend[,"skid"]
#   T4 <- read.neurons.catmaid(neu_skid_dend, .progress = 'text')
#   T4_dend[[LL]] <- kinkCT1(T4)
# }
# # SAVE
# save(T4_dend, anno_T4_dend, file = "data/neu_T4_dend.RData")

# T4 analysis -------------------------------------------------------------
grid_x <- seq(-4e4, 4e4, by = 1e2)
grid_y <- seq(-4e4, 4e4, by = 1e2)
bkgd_grid_col <- expand.grid(grid_x, grid_y) # 0.1 um

ucl_rot_gal <- ucl_rot_sm

dir_type <- list() 
lens_type <- list()
dir_pc_type <- list()
so_type <- list() #SO summary
seg_summ_type <- list() # tree segment summary
seg_summ_pc_type <- list()
ell_type <- list()

for (LL in 1:4) {
  neu <- T4_dend[[LL]]
  anno <- anno_T4_dend[[LL]]
  
  # init
  root_xyz <- matrix(ncol = 3, nrow = length(neu))
  dend_com <- matrix(ncol = 3, nrow = length(neu))
  T4_dir <- matrix(ncol = 3, nrow = length(neu)) # SO direction
  T4_dir_rs <- matrix(ncol = 6, nrow = length(neu)) # SO dir re-scaled
  T4_od_rs <- matrix(ncol = 6, nrow = length(neu)) # dir orthogonal to PD
  
  root_xyz_pc <- matrix(ncol = 3, nrow = length(neu)) # in pca coord
  T4_dir_pc <- matrix(ncol = 3, nrow = length(neu)) 
  T4_dir_rs_pc <- matrix(ncol = 6, nrow = length(neu)) # SO dir re-scaled
  
  T4_dir_lens <- matrix(ncol = 6, nrow = length(neu)) # SO dir on lens
  T4_com_lens <- matrix(ncol = 3, nrow = length(neu)) # com on lens
  
  ang_dir_pc1 <- c() # angle betw dir and T4 pc1
  
  # -- so stats, so_summ = cbind(so_max, so_N x4, so_vL x4, so_cL x4)
  so_max <- c() # max SO
  so_N <- matrix(0, ncol = 4, nrow = length(neu)) # num of SO branches, order 1 to 4
  so_vL <- matrix(0, ncol = 4, nrow = length(neu)) # vector length of SO branches
  so_cL <- matrix(0, ncol = 4, nrow = length(neu)) # cable length of SO branches
  
  # -- 2D quantity after pca projection
  # --- ellipsoid model, ell_summ = cbind(...)
  ell_abc <- matrix(ncol = 3, nrow = length(neu)) # 3 axes of ellipsoid along T4 dir
  ell_area_pt <- c() #num of points
  ell_area_piab <- c() # pi*ea*eb
  ell_vol_piabc <- c() # 4/3*pi*ea*eb*ec
  
  ell_ab_eye <- matrix(ncol = 3+1+1, nrow = length(neu)) # mapped onto eye, [center, ea, eb]
  
  # -- segment stats, a list of matrices
  seg_summ_neu <- list() # seg starting point, so, length, ang wrt dir
  seg_summ_pc_neu <- list()

  ind_rtcom <- c() # com and root coincide
  for (j in 1:length(neu)) {
    tar <- neu[[j]]
    ind_D = match(tar$tags$`dendrite start` , tar$d$PointNo)
    df_D <-  tar$d[ind_D,]
    root_xyz[j,] <- xyzmatrix(df_D)
    
    # - subtree with root = dendrite start 
    targ <- as.ngraph(tar)
    ii_root <- ind_D
    # subtree and Strahler order
    sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
    subtree <- subset(tar, sub_points) 
    subtree <- resample(subtree, stepsize = 100) #resample
    
    subtree_g <- as.ngraph(subtree, weights = T)
    subtree_so <- strahler_order(subtree) # Strahler order
    so_max <- c(so_max, max(subtree_so$segments))
    
    dend_com[j,] <- colMeans(xyzmatrix(subtree$d))

    # seg starting point, dir, so, path length, ang, path length to root
    seg_summ <- matrix(0, ncol = 3+3+1+1+1+1+1+1+1, nrow = length(subtree$SegList)) 
    for (ii_seg in 1:length(subtree$SegList)) {
      node_ii <- subtree$SegList[[ii_seg]]
      seg_v <- subtree$d[tail(node_ii,1), c("X","Y","Z")] - subtree$d[node_ii[1], c("X","Y","Z")]
      L <- igraph::distances(subtree_g, v = node_ii[1], to = tail(node_ii, 1) )
      seg_summ[ii_seg,] <- unlist(c(subtree$d[node_ii[1], c("X","Y","Z")], seg_v,
                                    subtree_so$segments[ii_seg], round(L,0), 0, 0, 0, 0,
                                    igraph::distances(subtree_g, v= subtree$StartPoint, to= node_ii[1])) )
    }
    seg_summ <- data.frame(seg_summ)
    colnames(seg_summ) <- c('x0','y0','z0','xd','yd','zd','so','L','ang','angV','angH','angSO4','pathL')
    
    so_v_mat <- matrix(0, ncol = 3, nrow = 4)
    for (so in 1:4) {
      seg_so <- seg_summ[seg_summ$so == so, ]
      so_N[j, so] <- nrow(seg_so)
      vsum <- colSums(seg_so[, c('xd','yd','zd')])
      so_vL[j, so] <- sqrt(sum(vsum^2))
      so_cL[j, so] <- sum(seg_so$L)
      so_v_mat[so, ] <- vsum 
    }
    
    T4_dir[j,] <- so_v_mat[2,] + so_v_mat[3,] 
    T4_dir[j,] <- T4_dir[j,] / sqrt(sum(T4_dir[j,]^2)) # normalize
    
    # dir based on so=2,3,4
    d4 <- so_v_mat[2,] + so_v_mat[3,] + so_v_mat[4,]
    d4 <- d4 / sqrt(sum(d4^2))
    
    # re-scale dir
    pv <- sweep(xyzmatrix(subtree$d),2, dend_com[j,],'-') %*% T4_dir[j,]
    nq <- quantile(pv, c(0.01, 0.99)) 
    T4_dir_rs[j,] <- c(dend_com[j,] + T4_dir[j,]*nq[1], dend_com[j,]+ T4_dir[j,]*nq[2])
    
    # - dir on lens
    # use nb_ind
    ii <- sweep(med_xyz, 2, dend_com[j,], '-')^2 %>% rowSums() %>% which.min()
    nb_ii <- c(nb_ind[nb_ind[ii,], ]) %>% unique()%>% na.omit() %>% c()
    
    # +v-axis
    vaxis <- med_xyz[nb_ind[nb_ii[1], c(6,3)], ] %>% diff() 
    vaxis <- c(vaxis / sqrt(sum(vaxis^2)))
    # +h-axis
    haxis <- colMeans(med_xyz[nb_ind[nb_ii[1], c(4,5)], ]) - colMeans(med_xyz[nb_ind[nb_ii[1], c(2,7)], ]) 
    haxis <- c(haxis / sqrt(sum(haxis^2)))
    
    # -- kernel regression
    # pca to flatten med neighborhood
    nb_xyz <- med_xyz[nb_ii, ] 
    nb_pca <- prcomp(nb_xyz)
    
    nb_xyz_pc <- sweep(nb_xyz, 2,nb_pca$center) %*% nb_pca$rotation
    colnames(nb_xyz_pc) <- c('x','y','z')
    
    vec_xyz <- matrix(T4_dir_rs[j,], ncol=3, byrow = T)
    vec_xyz_pc <- sweep(vec_xyz, 2,nb_pca$center) %*% nb_pca$rotation
    colnames(vec_xyz_pc) <- c('x','y','z')
    
    com_xyz <- matrix(dend_com[j,], ncol = 3)
    com_xyz_pc <- (com_xyz - nb_pca$center) %*% nb_pca$rotation
    colnames(com_xyz_pc) <- c('x','y','z')
    
    # -- orthogonal dir
    od_pc <- diff(vec_xyz_pc)
    od_pc <- c(od_pc[2], -od_pc[1],0)
    od <- od_pc %*% t(nb_pca$rotation)
    od <- od / sqrt(sum(od^2))
    pv <- sweep(xyzmatrix(subtree$d),2, dend_com[j,],'-') %*% t(od)
    nq <- quantile(pv, c(0.01, 0.99)) 
    T4_od_rs[j,] <- c(dend_com[j,] + od*nq[1], dend_com[j,]+ od*nq[2])
    
    xyz_vtail <- matrix(ncol = 3)  #base of vec
    xyz_vhead <- matrix(ncol = 3) #head of vec
    xyz_com <- c()  # com
    xyz_vtail_eval <- data.frame(mc.x = vec_xyz_pc[1,1], mc.y = vec_xyz_pc[1,2], mc.z = 0)
    xyz_vhead_eval <- data.frame(mc.x = vec_xyz_pc[2,1], mc.y = vec_xyz_pc[2,2], mc.z = 0)
    xyz_com_eval <- data.frame(mc.x = com_xyz_pc[1], mc.y = com_xyz_pc[2], mc.z = 0)
    for (k in 1:3) {
      npdata <- data.frame(mc = nb_xyz_pc, ec = ucl_rot_sm[nb_ii,k])
      bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
      model_np <- npreg(bw)
      xyz_vtail[k] <- predict(model_np, newdata = xyz_vtail_eval)
      xyz_vhead[k] <- predict(model_np, newdata = xyz_vhead_eval)
      xyz_com[k] <- predict(model_np, newdata = xyz_com_eval)
    }
    T4_com_lens[j,] <- xyz_com / sqrt(sum(xyz_com^2))
    T4_dir_lens[j,] <- c(xyz_vtail, xyz_vhead) / sqrt(sum(xyz_com^2))  #normalize to com, 2022-02-07
    
    # - angle betw seg and T4_dir
    for (ii_seg in 1:length(subtree$SegList)) {
      vv <- as.numeric(seg_summ[ii_seg,c('xd','yd','zd')]) 
      seg_summ[ii_seg,'ang'] <- round(acos(vv %*% T4_dir[j,] / sqrt(sum(vv^2))) /pi*180, 1)
      seg_summ[ii_seg,'angV'] <- round(acos(vv %*% vaxis / sqrt(sum(vv^2))) /pi*180, 1)
      seg_summ[ii_seg,'angH'] <- round(acos(vv %*% haxis / sqrt(sum(vv^2))) /pi*180, 1)
      seg_summ[ii_seg,'angSO4'] <- round(acos(vv %*% d4 / sqrt(sum(vv^2))) /pi*180, 1)
    }
    
    # - pca subtree
    tar_pca <- prcomp(xyzmatrix(subtree$d)) # in rotated coord, x = 1st pc
    pc1 <- tar_pca$rotation[,1]
    pc2 <- tar_pca$rotation[,2]
    ang_dir_pc1 <- c(ang_dir_pc1, acos(abs(T4_dir[j,] %*% pc1)))
    Dref <- T4_dir[j,] # use PD as ref direction
    if (LL <= 2) {
      tar_pca$rotation[,1] <- pc2
      tar_pca$rotation[,2] <- pc1
    }
    if ( tar_pca$rotation[,1] %*% Dref < 0 ) {
      tar_pca$rotation[,1] <- - tar_pca$rotation[,1]
    }
    if (t(cross3D(tar_pca$rotation[,1],tar_pca$rotation[,2])) %*% tar_pca$rotation[,3] < 0 ) {
      tar_pca$rotation[,3] <- - tar_pca$rotation[,3]
    }
    
    # - choose 2D origin
    ori_xyz <- tar_pca$center
    
    # - transform to pc coord
    root_xyz_pc[j,] <- (root_xyz[j,] - ori_xyz) %*% tar_pca$rotation
    T4_dir_pc[j,] <- T4_dir[j,] %*% tar_pca$rotation
    T4_dir_rs_pc[j,1:3] <- (T4_dir_rs[j,1:3] - ori_xyz) %*% tar_pca$rotation
    T4_dir_rs_pc[j,4:6] <- (T4_dir_rs[j,4:6] - ori_xyz) %*% tar_pca$rotation
    
    seg_summ_pc <- seg_summ
    seg_summ_pc[, c('x0','y0','z0')] <- sweep(as.matrix(seg_summ_pc[, c('x0','y0','z0')]),2,ori_xyz,'-') %*% tar_pca$rotation
    seg_summ_pc[, c('xd','yd','zd')] <- as.matrix(seg_summ_pc[, c('xd','yd','zd')]) %*% tar_pca$rotation
    
    seg_summ_neu[[j]] <- seg_summ
    seg_summ_pc_neu[[j]] <- seg_summ_pc
    
    # - ellipsoid
    subtree_pc <- subtree
    subtree_pc$d[, c("X","Y","Z")] <- sweep(xyzmatrix(subtree$d), 2, ori_xyz, '-') %*% tar_pca$rotation
    xyz <- xyzmatrix(subtree_pc$d)
    ell_abc[j,1] <- diff(quantile(xyz[,1], c(0.01, 0.99))) / 2
    ell_abc[j,2] <- diff(quantile(xyz[,2], c(0.01, 0.99))) / 2
    ell_abc[j,3] <- diff(quantile(xyz[,3], c(0.01, 0.99))) / 2
   
    # - find the nb in med and eye
    nb_N <- 1 + 8 + 8*2
    nb_N_col <- 1 + 8 + 8*2 + 8*3
    nbd <- rowSums(sweep(med_xyz, 2, ori_xyz, '-')^2)
    ind_nb_N <- order(nbd)[1:nb_N]
    ind_nb_N_col <- order(nbd)[1:nb_N_col]
    nb <- med_xyz[ind_nb_N, ]
    nb_col <- med_xyz[ind_nb_N_col, ]
    
    nb_pc <- sweep(nb, 2, ori_xyz, '-') %*% tar_pca$rotation
    nb_col_pc <- sweep(nb_col, 2, ori_xyz, '-') %*% tar_pca$rotation
    
    # -  count grid points
    hpts_ind <- chull(xyzmatrix(subtree_pc$d)[,1:2]) #chull from [x,y] of subtree
    hpts_ind <- c(hpts_ind, hpts_ind[1])
    hpts_xy <- xyzmatrix(subtree_pc$d)[hpts_ind, 1:2] # hull edge points
    
    hpts_ind_col <- chull(nb_col_pc[, 1:2])
    hpts_ind_col <- c(hpts_ind_col, hpts_ind_col[1])
    hpts_xy_col <- nb_col_pc[hpts_ind_col, 1:2] # hull edge points
    poly_st_col <- st_polygon(list(hpts_xy_col))
    in_bkgd_col <- sp::point.in.polygon(bkgd_grid_col[,1], bkgd_grid_col[,2], hpts_xy_col[,1], hpts_xy_col[,2])
    in_bkgd <- sp::point.in.polygon(bkgd_grid_col[,1], bkgd_grid_col[,2], hpts_xy[,1], hpts_xy[,2])
    
    # - T4 ellipsoid
    # in_bkgd <- sp::point.in.polygon(bkgd_grid[,1], bkgd_grid[,2], hpts_xy[,1], hpts_xy[,2])
    ell_area_pt <- c(ell_area_pt, round(sum(in_bkgd)/100, digits = 1) )
    ell_area_piab <- c(ell_area_piab, round(pi*ell_abc[j,1]*ell_abc[j,2]/1e6, digits = 1) )
    ell_vol_piabc <- c(ell_vol_piabc, round(4/3*pi*ell_abc[j,1]*ell_abc[j,2]*ell_abc[j,3], digits=1))
    
    # - map onto eye
    ucl_nb <- ucl_rot_gal[ind_nb_N, ]
    
    # np
    xyz <- rbind(dend_com[j,],
                 dend_com[j,]-tar_pca$rotation[,1]*ell_abc[j,1],
                 dend_com[j,]+tar_pca$rotation[,1]*ell_abc[j,1],
                 dend_com[j,]-tar_pca$rotation[,2]*ell_abc[j,2],
                 dend_com[j,]+tar_pca$rotation[,2]*ell_abc[j,2] )
    xyz_eval <- data.frame(mc.x = xyz[,1], mc.y = xyz[,2], mc.z = xyz[,3])
    np_ab <- matrix(nrow = nrow(xyz_eval), ncol = ncol(xyz_eval))
    # -- kernel regression
    for (k in 1:3) {
      npdata <- data.frame(mc = med_xyz[ind_nb_N,], ec = ucl_rot_gal[ind_nb_N,k]) 
      model_np <- npreg(ec ~ mc.x + mc.y + mc.z, data = npdata)
      for (m in 1:nrow(xyz_eval)) {
        np_ab[m, k] <- predict(model_np, newdata = xyz_eval[m,])
      }
    }
    np_ab <- sweep(np_ab, 1, sqrt(rowSums(np_ab^2)), '/')
    
    ab_tp <- cart2sphZ(np_ab[2:5,])[,c('theta','phi')] 
    # ellipse, [center, ea, eb]
    ell_ab_eye[j,] <- c(np_ab[1,],
                        arcLength(ab_tp[1,], ab_tp[2,])/2,
                        arcLength(ab_tp[3,], ab_tp[4,])/2 )
  }
  
  tmp <- as.data.frame(cbind(root_xyz, dend_com, T4_dir, T4_dir_rs, ang_dir_pc1, T4_od_rs))
  colnames(tmp) <- c('rtx','rty','rtz','comx','comy','comz','xd','yd','zd','rsx0','rsy0','rsz0','rsxd','rsyd','rszd','ang','odx0','ody0','odz0','odxd','odyd','odzd')
  dir_type[[LL]] <- tmp
  
  tmp <- as.data.frame(cbind(T4_dir_lens, T4_com_lens))
  colnames(tmp) <- c('x0','y0','z0','xd','yd','zd','comx','comy','comz')
  lens_type[[LL]] <- tmp
  
  tmp <- as.data.frame(cbind(root_xyz_pc, T4_dir_pc, T4_dir_rs_pc))
  colnames(tmp) <- c('pcrtx','pcrty','pcrtz','pcxd','pcyd','pczd','pcrsx0','pcrsy0','pcrsz0','pcrsxd','pcrsyd','pcrszd')
  dir_pc_type[[LL]] <- tmp
  
  tmp <- as.data.frame(cbind(so_max, so_N, so_vL, so_cL))
  colnames(tmp) <- c('so_max', paste('so_N',1:4,sep=''), paste('so_vL',1:4,sep=''),paste('so_cL',1:4,sep=''))
  so_type[[LL]] <- tmp
  
  tmp <- as.data.frame(cbind(ell_abc, ell_vol_piabc, ell_area_piab, ell_area_pt, ell_ab_eye))
  colnames(tmp) <- c('ea','eb','ec', 'V_piabc', 'A_piab','A_pt','comx','comy','comz','eye_ea','eye_eb')
  ell_type[[LL]] <- tmp
  
  seg_summ_type[[LL]] <- seg_summ_neu
  seg_summ_pc_type[[LL]] <- seg_summ_pc_neu
}

# # SAVE
# save(dir_type, lens_type, dir_pc_type, so_type, seg_summ_type, seg_summ_pc_type, ell_type,
#      file = "data/T4_gallery.RData")


# npreg vector field on lens, T4b, use global np with fixed bw  ------------------------

LL <- 2

RF_lens <- lens_type[[LL]][, c('x0','y0','z0','xd','yd','zd')] %>% as.matrix()
colnames(RF_lens) <- NULL
RF_com <- lens_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
colnames(RF_com) <- NULL

np_eval <- ucl_rot_sm
vf_pred <- matrix(nrow = dim(np_eval)[1], ncol = 3)
for (k in 1:3) {
  npdata <- data.frame(mc = RF_com, ec = RF_lens[,3+k]) #use com
  bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'fixed', regtype= 'll')
  model_np <- npreg(bw)
  for (j in 1:nrow(np_eval)) {
    np_eval_one <- data.frame(mc.1 = np_eval[j,1], mc.2 = np_eval[j,2], mc.3 = np_eval[j,3])
    vf_pred[j,k] <- predict(model_np, newdata = np_eval_one)
  }
}
RF_lens_T4b_pred <- 2*(vf_pred - ucl_rot_sm) + ucl_rot_sm # restore length

# use nb_ind to smooth
vv <- RF_lens_T4b_pred - ucl_rot_sm
vvnew <- vv
for (m in 1:nrow(nb_ind)) {
  if (!any(is.na(nb_ind[m,]))) {
    vvnew[nb_ind[m,1],] <- vv[nb_ind[m,1],]/2 + colMeans(vv[nb_ind[m,-1],])/2
  }
}
RF_lens_T4b_pred_sm <- ucl_rot_sm + vvnew
RF_lens_T4b_pred_sm <- sweep(RF_lens_T4b_pred_sm,1,sqrt(rowSums(RF_lens_T4b_pred_sm^2)),'/') #normalize

# npreg vector field on lens, T4d ------------------------

LL <- 4

RF_lens <- lens_type[[LL]][, c('x0','y0','z0','xd','yd','zd')] %>% as.matrix()
colnames(RF_lens) <- NULL
RF_com <- lens_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
colnames(RF_com) <- NULL

# - np
np_eval <- ucl_rot_sm
vf_pred <- matrix(nrow = dim(np_eval)[1], ncol = 3)
for (k in 1:3) {
  npdata <- data.frame(mc = RF_com, ec = RF_lens[,3+k] ) #use com
  bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'fixed', regtype= 'll')
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

# all T4 direction assuming antiparallelism  -----------------------------------

# This is RF or preferred direction (PD), opposite of predicated orientation of T4 dendrites
RF_lens_T4_pred <- cbind(RF_lens_T4b_pred,
                         ucl_rot_sm*2 - RF_lens_T4b_pred,
                         RF_lens_T4d_pred,
                         ucl_rot_sm*2 - RF_lens_T4d_pred )

RF_lens_T4_pred_sm <- cbind(RF_lens_T4b_pred_sm,
                         ucl_rot_sm*2 - RF_lens_T4b_pred_sm,
                         RF_lens_T4d_pred_sm,
                         ucl_rot_sm*2 - RF_lens_T4d_pred_sm )

# # SAVE
# save(RF_lens_T4_pred, RF_lens_T4_pred_sm, file = "data/T4_RF_pred.RData")


