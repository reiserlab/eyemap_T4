# T4 morphology 
# cf  setwd("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_T4T5")

library(natverse)
library(tidyverse)
# library(ggplot2)
library(igraph)
library(rjson)
library(RColorBrewer)
library(sf) #intersection
library(deldir)
library(gridExtra) #save ggplot to pdf

library(np) 

library(ggraph) #morphology analysis

# clean everythign up.
rm(list=ls())

source("eyemap_func.R")

pal_so <- c('gray', brewer.pal(5,"RdBu")[c(2,1,5,4)])

# pal_T4 <- brewer.pal(5, "BrBG")[c(2,1,4,5)]
pal_T4 <- c("turquoise", "#A6611A",  "royalblue", "plum")

# load data ---------------------------------------------------------------

# load("data/CT1_Mi1.RData")
# load("data/neu_CT1.RData")
# load("data/neu_Mi1.RData")

load("data/eyemap.RData")

# update 2022-01-26
load('data/hexnb_ind_dist.RData') #  in 'eyemap_lens_med_lop.R'


# load neurons -------------------------------------------------------------------------------------------------

load("data/neu_T4_dend.RData")


# # - T4
# T4_dend <- list()
# anno_T4_dend <- list()
# for (LL in 1:4) {
#   # anno_str_axon <- paste("T4",letters[LL], " - axon", sep = "")
#   anno_str_dend <- paste("T4",letters[LL], " - dendrites", sep = "")
#   # anno_str_comp <- paste("T4",letters[LL], " - complete", sep = "")
#   # anno_axon <- catmaid_query_by_annotation(anno_str_axon)
#   anno_dend <- catmaid_query_by_annotation(anno_str_dend)
#   # anno_comp <- catmaid_query_by_annotation(anno_str_comp)
#   anno_T4_dend[[LL]] <- anno_dend
#   # neu_skid_axon <- anno_axon[,"skid"]
#   neu_skid_dend <- anno_dend[,"skid"]
#   # neu_skid_comp <- anno_comp[,"skid"]
#   T4 <- read.neurons.catmaid(neu_skid_dend, .progress = 'text')
#   T4_dend[[LL]] <- kinkCT1(T4)
#   # T4[[LL]] <- resample(T4_LL, stepsize = 100) # NO -> lose tags
# }
# 
# # # SAVE
# # save(T4_dend, anno_T4_dend, file = "data/neu_T4_dend.RData")
 
# # LL <- 2
# # nopen3d()
# # for (j in 1:length(T4[[LL]])) {
# #   tar <- T4[[LL]][[j]]
# #   targ <- as.ngraph(tar)
# #   ii_root <- match(tar$tags$`dendrite start`, tar$d[,'PointNo']) # use a depth first search
# #   sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
# #   sub_tree <- subset(tar, sub_points, invert = F)
# #   plot3d(sub_tree, col = j, lwd = 2, WithNodes = F)
# # }

# # T4 tracing --------------------------------------------------------------
# 
# # check for multi-marker
# neu <- T4_dend[[LL]]
# ind_cor <- c()
# ind_br <- c()
# for (j in 1:length(neu)) {
#   tar <- neu[[j]]
#   ind_D = match(tar$tags$`dendrite start` , tar$d$PointNo)
#   if (length(ind_D)==0) {
#     ind_cor <- c(ind_cor,j)
#   }
#   if (length(ind_D) > 1) {
#     ind_br <- c(ind_br,j)
#   }
# }



# branching angles --------------------------------------------------------

bp3 <- matrix(ncol = 3, nrow = 0) # non-3 branching point
seg_ang_type <- list()
for (LL in 1:4) {
  neu <- T4_dend[[LL]]
  anno <- anno_T4_dend[[LL]]
  
  seg_ang_neu <- list()
  for (j in 1:length(neu)) {
    tar <- neu[[j]]
    ind_D = match(tar$tags$`dendrite start` , tar$d$PointNo)
    # df_D <-  tar$d[ind_D,]
    # root_xyz[j,] <- xyzmatrix(df_D)
    
    # - subtree with root = dendrite start 
    targ <- as.ngraph(tar)
    ii_root <- ind_D
    # subtree and Strahler order
    sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
    subtree <- subset(tar, sub_points) 
    subtree <- resample(subtree, stepsize = 100) #resample
    
    subtree_g <- as.ngraph(subtree, weights = T)
    subtree_so <- strahler_order(subtree) # Strahler order
    
    xyz <- xyzmatrix(subtree$d)
    seg_ang <- matrix(ncol = 6, nrow = 0) # [ang, LL, j, ind parent, so parent, so daughter]
    for (bp in subtree$BranchPoints) {
      ind_seg <- sapply(subtree$SegList, function(x) bp %in% x) %>% which()
      # check if 3 branches
      if (length(ind_seg) != 3) {
        bp3 <- rbind(bp3, c(LL, j, bp))
        next
      }
      # dir vec
      vda <- matrix(ncol = 3, nrow = 0)
      ind_vda <- c()
      for (iiseg in ind_seg) {
        seg <- subtree$SegList[[iiseg]]
        if (bp == head(seg,1)) {
          vda <- rbind(vda, diff(xyz[head(seg,2),])) # daughter vector
          ind_vda <- c(ind_vda, iiseg)
        } else if (bp == tail(seg,1)) {
          vpa <- diff(xyz[tail(seg,2),]) #parent vec
          ind_vpa <- iiseg
        }
      }
      # angles
      for (iia in 1:2) {
        ca <- vpa %*% vda[iia,] / sqrt(sum(vpa^2)) / sqrt(sum(vda[iia,]^2))
        ca <- if_else(ca > 1, 1, ca)
        ca <- if_else(ca < -1, -1, ca)
        ang <- acos(ca) /pi*180
        seg_ang <- rbind(seg_ang, c(round(ang,2), LL, j, ind_vpa, subtree_so$segments[ind_vpa], subtree_so$segments[ind_vda[iia]]) )
      }
    }
    seg_ang_neu[[j]] <- seg_ang
  }
  seg_ang_type[[LL]] <- seg_ang_neu
}
    

# - distr angles
LL <- 4
ang_so <- do.call(rbind, seg_ang_type[[LL]])
colnames(ang_so) <- c('ang', 'type', 'neu', 'indseg','soPa','soDa')
ang_so <- ang_so[!is.na(ang_so[,1]),]
# df <- ang_so %>%
#   as_tibble() %>%
#   filter(soPa == 2) %>%
#   group_by(indseg) %>%
#   # mutate(n = n()) %>%
#   head()

# ang diff betw daughters
ang <- ang_so[which(ang_so[,'soPa'] == 2),]
angdiff <- matrix(ncol = ncol(ang), nrow = 0)
for (j in 1:(nrow(ang)/2)) {
  if (ang[2*j-1,'indseg'] != ang[2*j, 'indseg']) {
    print(j)
    stop("check")
  }
  angdiff <- rbind(angdiff, c(abs(diff(ang[(1:2)+2*j-2, 'ang'])), ang[2*j-1, -1]) )
}

# T4 analysis -------------------------------------------------------------
# dir differs from 'lens_med_lop.R' by child_node() vs resample()

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
  # ell_e <- c()
  
  ell_ab_eye <- matrix(ncol = 3+1+1, nrow = length(neu)) # mapped onto eye, [center, ea, eb]

  # -- segment stats, a list of matrices
  seg_summ_neu <- list() # seg starting point, so, length, ang wrt dir
  seg_summ_pc_neu <- list()
  
  # T4_colnum <- c() # num of col covered by T4
  # T4_colnum_para <- c()# num of col along preferred dir
  # T4_colnum_perp <- c()
  # L_inter_col <- c() #inter-col distance
  # col_area <- c() # local col area, assuming hex
  
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
    
    # # DEBUG
    # nopen3d()
    # # plot3d(subtree, col = 'grey')
    # for (so in 1:max(subtree_so$points)) {
    # # for (so in 1:3) {
    #   pt_p <- subtree$d$PointNo %in% unique(subtree$d[subtree_so$points==so,"Parent"])
    #   pt_so <- pt_p | subtree_so$points==so
    #   plot3d(subset(subtree, pt_so), col=pal_so[so], add = T, boundingbox = boundingbox(subtree), lwd= 3, WithNodes = F)
    # }
    # points3d(xyzmatrix(df_D), size=10)
    # # plot3d(tar, col='gray', lwd=3, soma=T)
    
    # >>>
    # define orientation based on strahler order
    # then calculate pca and adjust according to so dir
    # then get the seg angle in the pca ref frame
    
    # so_v_neu <- list() #list of so of individual branch direction, for histo
    
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
    # coord_D[j, ] <- as.numeric(subtree$d[subtree$StartPoint, c("X","Y","Z")])
    seg_summ <- data.frame(seg_summ)
    colnames(seg_summ) <- c('x0','y0','z0','xd','yd','zd','so','L','ang','angV','angH','angSO4','pathL')
    
    so_v_mat <- matrix(0, ncol = 3, nrow = 4)
    for (so in 1:4) {
      seg_so <- seg_summ[seg_summ$so == so, ]
      so_N[j, so] <- nrow(seg_so)
      vsum <- colSums(seg_so[, c('xd','yd','zd')])
      so_vL[j, so] <- sqrt(sum(vsum^2))
      so_cL[j, so] <- sum(seg_so$L)
      # so_v_mat[so, ] <- vsum / so_vL[j, so] # normalize per so, Diff from neu_dir_SO
      so_v_mat[so, ] <- vsum # Same as in neu_dir_SO
    }
    
    # T4_dir[j,] <- (so_v_mat[2,] + so_v_mat[3,]) / 2 #weight so = 2,3 for normalize per so
    T4_dir[j,] <- so_v_mat[2,] + so_v_mat[3,] # difference between this and eyemap_T4_RF is sumType='avg' for the latter
    T4_dir[j,] <- T4_dir[j,] / sqrt(sum(T4_dir[j,]^2)) # normalize
    
    # dir based on so=2,3,4
    # d4 <- so_v_mat[2,] + so_v_mat[3,] + so_v_mat[1,] + so_v_mat[4,]
    d4 <- so_v_mat[2,] + so_v_mat[3,] + so_v_mat[4,]
    d4 <- d4 / sqrt(sum(d4^2))
    
    # re-scale dir
    pv <- sweep(xyzmatrix(subtree$d),2, dend_com[j,],'-') %*% T4_dir[j,]
    # nq <- quantile(pv, c(0.01,0.99)) %>% diff()
    # T4_dir_rs[j,] <- c(dend_com[j,]- T4_dir[j,]*nq/2, dend_com[j,]+ T4_dir[j,]*nq/2)
    nq <- quantile(pv, c(0.01, 0.99)) # 2021-11-11
    T4_dir_rs[j,] <- c(dend_com[j,] + T4_dir[j,]*nq[1], dend_com[j,]+ T4_dir[j,]*nq[2])
    
    # - dir on lens
    # ## ## use nearest
    # nbhd_N <- 1+8+16 # num of nbs
    # nb_ii <- sweep(med_xyz, 2, dend_com[j,], '-')^2 %>% rowSums() %>% order() %>% head(nbhd_N)
    ## ## use nb_ind
    ii <- sweep(med_xyz, 2, dend_com[j,], '-')^2 %>% rowSums() %>% which.min()
    nb_ii <- c(nb_ind[nb_ind[ii,], ]) %>% unique()%>% na.omit() %>% c()
    # nb_ii <- c(nb_ind[nb_ind[match(ii,nb_ind[,1]),], ]) %>% unique()%>% na.omit() %>% c() #2023
    
    # +v-axis
    # vaxis <- med_xyz[nb_ind[nb_ii[1], c(7,6)], ] %>% diff() #2021
    vaxis <- med_xyz[nb_ind[nb_ii[1], c(6,3)], ] %>% diff() #2023
    vaxis <- c(vaxis / sqrt(sum(vaxis^2)))
    # +h-axis
    # haxis <- colMeans(med_xyz[nb_ind[nb_ii[1], c(2,3)], ]) - colMeans(med_xyz[nb_ind[nb_ii[1], c(4,5)], ]) #2021
    haxis <- colMeans(med_xyz[nb_ind[nb_ii[1], c(4,5)], ]) - colMeans(med_xyz[nb_ind[nb_ii[1], c(2,7)], ]) #2023
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
    # arrow3d(dend_com[j,], dend_com[j,]+T4_dir[j,]*1e4, theta = pi / 9, n = 8, col = "red", type = "rotation")
    # arrow3d(dend_com[j,], dend_com[j,]+od*1e4, theta = pi / 9, n = 8, col = "blue", type = "rotation")
    pv <- sweep(xyzmatrix(subtree$d),2, dend_com[j,],'-') %*% t(od)
    # od_L <- c(od_L, diff(quantile(pv, c(0.01,0.99))) )
    # nq <- quantile(pv, c(0.01,0.99)) %>% diff()
    # T4_od_rs[j,] <- c(dend_com[j,]- od*nq/2, dend_com[j,]+ od*nq/2)
    nq <- quantile(pv, c(0.01, 0.99)) # 2021-11-11
    T4_od_rs[j,] <- c(dend_com[j,] + od*nq[1], dend_com[j,]+ od*nq[2])
    
    xyz_vtail <- matrix(ncol = 3)  #base of vec
    xyz_vhead <- matrix(ncol = 3) #head of vec
    xyz_com <- c()  # com
    xyz_vtail_eval <- data.frame(mc.x = vec_xyz_pc[1,1], mc.y = vec_xyz_pc[1,2], mc.z = 0)
    xyz_vhead_eval <- data.frame(mc.x = vec_xyz_pc[2,1], mc.y = vec_xyz_pc[2,2], mc.z = 0)
    xyz_com_eval <- data.frame(mc.x = com_xyz_pc[1], mc.y = com_xyz_pc[2], mc.z = 0)
    # xyz_vtail_eval <- data.frame(mc.x = T4_dir_rs[j,1], mc.y = T4_dir_rs[j,2], mc.z = T4_dir_rs[j,3])
    # xyz_vhead_eval <- data.frame(mc.x = T4_dir_rs[j,4], mc.y = T4_dir_rs[j,5], mc.z = T4_dir_rs[j,6])
    # xyz_com_eval <- data.frame(mc.x = dend_com[j,1], mc.y = dend_com[j,2], mc.z = dend_com[j,3])
    for (k in 1:3) {
      # npdata <- data.frame(mc = med_xyz[nb_ii,], ec = ucl_rot_sm[nb_ii,k]) # with re-defined ucl_rot_sm and med_xyz
      npdata <- data.frame(mc = nb_xyz_pc, ec = ucl_rot_sm[nb_ii,k]) # with re-defined ucl_rot_sm and med_xyz
      
      # bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'generalized_nn', regtype= 'll')# with 1+8+16, all lc works well
      # bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'generalized_nn', regtype= 'lc') # lc seems to generate tiny arrows
      # bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'fixed', regtype= 'll')
      bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
      model_np <- npreg(bw)
      xyz_vtail[k] <- predict(model_np, newdata = xyz_vtail_eval)
      xyz_vhead[k] <- predict(model_np, newdata = xyz_vhead_eval)
      xyz_com[k] <- predict(model_np, newdata = xyz_com_eval)
    }
    # utp_vtail_norm <- sqrt(sum(xyz_vtail^2))
    # T4_dir_lens[j,] <- c(xyz_vtail, xyz_vhead) / utp_vtail_norm  #normalize to vtail
    # v1 <- T4_dir_lens[j,1:3]
    # v2 <- T4_dir_lens[j,4:6]
    # T4_dir_lens[j,4:6] <- v2 - c((v2 - v1) %*% v1) * v1 #take tangent
    T4_com_lens[j,] <- xyz_com / sqrt(sum(xyz_com^2))
    T4_dir_lens[j,] <- c(xyz_vtail, xyz_vhead) / sqrt(sum(xyz_com^2))  #normalize to com, 2022-02-07
    
    
    
    # - angle betw seg and T4_dir
    for (ii_seg in 1:length(subtree$SegList)) {
      # dir_rot <- as.numeric(seg_summ[ii_seg,c('xd','yd','zd')]) %*% tar_pca$rotation
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
    # Dref <- dend_com - xyz_D_avg # use com-mean as ref dir
    Dref <- T4_dir[j,] # use PD as ref direction
    if (LL <= 2) {
      # if ( abs(pc1 %*% Dref) > abs(pc2 %*% Dref)) {
        tar_pca$rotation[,1] <- pc2
        tar_pca$rotation[,2] <- pc1
      # }
    }
    # if (LL >= 3) {
    #   if ( abs(pc1 %*% Dref) < abs(pc2 %*% Dref)) {
    #     tar_pca$rotation[,1] <- pc2
    #     tar_pca$rotation[,2] <- pc1
    #   }
    # }
    if ( tar_pca$rotation[,1] %*% Dref < 0 ) {
      tar_pca$rotation[,1] <- - tar_pca$rotation[,1]
    }
    if (t(cross3D(tar_pca$rotation[,1],tar_pca$rotation[,2])) %*% tar_pca$rotation[,3] < 0 ) {
      tar_pca$rotation[,3] <- - tar_pca$rotation[,3]
    }
    
    # - choose 2D origin
    # ori_xyz <- xyz_D_avg
    # ori_xyz <- dend_com
    ori_xyz <- tar_pca$center
    
    # - transform to pc coord
    root_xyz_pc[j,] <- (root_xyz[j,] - ori_xyz) %*% tar_pca$rotation
    # so_v_mat_pc <- so_v_mat %*% tar_pca$rotation
    # so_v_mat_raw <- so_v_mat_raw %*% tar_pca$rotation
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
    
    
    # # ---  DEBUG
    # nopen3d()
    # for (k in 1:max(subtree_so$points)) {
    #   pt_p <- subtree$d[,"PointNo"] %in% unique(subtree$d[subtree_so$points == k,"Parent"])
    #   pt_so <- pt_p | subtree_so$points==k
    #   plot3d(subset(subtree, pt_so), col=pal_so[k], add = k!=1, boundingbox = boundingbox(tar), lwd = 3)
    #   # plot3d(subset(subtree_pc, pt_so), col=pal_so[k], add = k!=1, boundingbox = boundingbox(tar), lwd = 3)
    # }
    # points3d(matrix(root_xyz[j,], ncol = 3), col = 'grey', size = 25)
    # for (ii_so in 1:4) {
    #   arrow3d(root_xyz[j,], root_xyz[j,]+so_v_mat[ii_so,]*1e3, theta = pi / 9, n = 8, col = pal_so[ii_so], type = "rotation")
    # }
    # 
    # nopen3d()
    # for (k in 1:max(subtree_so$points)) {
    #   pt_p <- subtree$d[,"PointNo"] %in% unique(subtree$d[subtree_so$points == k,"Parent"])
    #   pt_so <- pt_p | subtree_so$points==k
    #   plot3d(subset(subtree_pc, pt_so), col=pal_so[k], add = k!=1, boundingbox = boundingbox(subtree_pc), lwd = 3)
    #   # plot3d(subset(subtree_pc, pt_so), col=pal_so[k], add = k!=1, boundingbox = boundingbox(tar), lwd = 3)
    # }
    # points3d(matrix(root_xyz_pc[j,], ncol = 3), col = 'grey30', size = 25)
    # points3d(matrix(c(0,0,0), ncol = 3), col = 'grey', size = 25)
    # for (ii_so in 1:4) {
    #   arrow3d(root_xyz_pc[j,], root_xyz_pc[j,]+so_v_mat_pc[ii_so,]*1e3, theta = pi / 9, n = 8, col = pal_so[ii_so], type = "rotation")
    # }
    # arrow3d(root_xyz_pc[j,], root_xyz_pc[j,]+T4_dir_pc[j,]*1e4, theta = pi / 9, n = 8, col = 'dark green', type = "rotation")
    # nopen3d()
    # plot3d(subtree, WithNodes = F, lwd=2, col='gray')
    # points3d(root_xyz,size =30)
    # arrow3d(root_xyz[j,], root_xyz[j,]+T4_dir[j,]*5e3, theta = pi / 9, n = 8, col = 'dark green', type = "rotation")
    # arrow3d(tar_pca$center, tar_pca$center+tar_pca$rotation[,1]*5e3, theta = pi / 9, n = 8, col = 'blue', type = "rotation")
    # arrow3d(tar_pca$center, tar_pca$center+tar_pca$rotation[,2]*5e3, theta = pi / 9, n = 8, col = 'red', type = "rotation")
    
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
    
    # ### ### ellipticity via polygon
    # hpts_xy_med <- hpts_xy
    # poly_st <- st_polygon(list(hpts_xy))
    # # line_st_pc1 <- st_linestring(matrix(c(-2e4, dend_com_pc[2], 2e4, dend_com_pc[2]), nrow = 2, byrow = T)) # x-axis
    # line_st_pc1 <- st_linestring(matrix(c(dend_com_pc[1]-T4_dir[j,1]*2e5, dend_com_pc[2]-T4_dir[j,2]*2e5, dend_com_pc[1]+T4_dir[j,1]*2e5, dend_com_pc[2]+T4_dir[j,2]*2e5), nrow = 2, byrow = T)) # x-axis
    # int_st_pc1 = st_intersection(line_st_pc1, poly_st) # intersections
    # # line_st_pc2 <- st_linestring(matrix(c(dend_com_pc[1], -2e4, dend_com_pc[1], 2e4), nrow = 2, byrow = T))
    # line_st_pc2 <- st_linestring(matrix(c(dend_com_pc[1]+T4_dir[j,2]*2e5, dend_com_pc[2]-T4_dir[j,1]*2e5, dend_com_pc[1]-T4_dir[j,2]*2e5, dend_com_pc[2]+T4_dir[j,1]*2e5), nrow = 2, byrow = T)) # x-axis
    # int_st_pc2 = st_intersection(line_st_pc2, poly_st) #  intersection
    # 
    # 
    # # intersections of major and minor axes
    # ea <- sqrt(sum(as.numeric(diff(int_st_pc1))^2)) / 2 / 1000 #[um]
    # eb <- sqrt(sum(as.numeric(diff(int_st_pc2))^2)) / 2 / 1000
    # 
    # ell_abc[j,1] <- ea
    # ell_abc[j,2] <- eb
    
    
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
save(dir_type, lens_type, dir_pc_type, so_type, seg_summ_type, seg_summ_pc_type, ell_type,
     file = "data/T4_gallery.RData")


# generate T4b gallery -----------------------------------------------------------------------------------------------------
# choose nb via nb_ind
# align to meridian by mapping ref points to med via np

load('data/T4_gallery.RData')
load('data/hexnb_ind_dist.RData') #  in 'eyemap_lens_med_lop.R'

nb_coord <- rbind(c(0,0),
                  c(-5/180*pi,0),
                  c(+5/180*pi,0),
                  c(0, -5/180*pi),
                  c(0, +5/180*pi) )

# Mollweide guidelines
Mollweide_ori <- c(-16e3, 7e3)
Mollweide_mul <- 1000

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


# main
for (LL in 1:4) {
  pdf(file = paste("plots/gallery_v10_T4", letters[LL], ".pdf", sep = ''))
  
  # LL <- 2 # choose type
  # eye
  com_xyz_eye <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
  T4_utp <- cart2sph2tp(com_xyz_eye) 
  # med
  com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
  v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix()
  v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
  v3 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
  v4 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()
  # eab <- ell_type[[LL]][,c('ea','eb')] %>% as.matrix()
  
  for (j in 1:length(T4_dend[[LL]])) {
  # for (j in 1:3) {
    # 1+6+12 hex nb
    ii <- sweep(ucl_rot_sm, 2, com_xyz_eye[j, ] )^2 %>% rowSums() %>% which.min()
    ii2 <- c(nb_ind[nb_ind[ii,],]) %>% unique()
    
    # ref nb
    pt <- ucl_rot_sm[ii, ] %>% matrix(ncol=3)
    tp <- cart2sphZ(pt)[,2:3]
    rtp <- cbind(1, sweep(nb_coord,2,tp,'+'))
    xyz_ref <- sph2cartZ(rtp)
    
    # np to med
    xyz_ref_med <- matrix(nrow = nrow(xyz_ref), ncol = 3) 
    xyz_eval <- data.frame(mc.x = xyz_ref[,1], mc.y = xyz_ref[,2], mc.z = xyz_ref[,3])
    for (k in 1:3) {
      npdata <- data.frame(mc = ucl_rot_sm[ii2,], ec = med_xyz[ii2,k])
      bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
      model_np <- npreg(bw)
      xyz_ref_med[,k] <- predict(model_np, newdata = xyz_eval)
    }
    
    # pc
    # zz <- cross3D((xyz_ref_med[2,]-xyz_ref_med[1,]), (xyz_ref_med[5,]-xyz_ref_med[1,])) #pointing inwards
    zz <- -(cross3D((xyz_ref_med[2,]-xyz_ref_med[1,]), (xyz_ref_med[5,]-xyz_ref_med[1,]))) #pointing outwards
    pc <- prcomp(xyz_ref_med)
    if (pc$rotation[,3] %*% zz < 0) {
      pc$rotation <- -pc$rotation
    }
    if (c(cross3D(pc$rotation[,3], pc$rotation[,1])) %*% pc$rotation[,2] < 0) {
      pc$rotation[,2] <- -pc$rotation[,2]
    }
    
    xyz_ref_med_pc <- sweep(xyz_ref_med, 2, xyz_ref_med[1,]) %*% pc$rotation
    y_ref <- diff(xyz_ref_med_pc[c(3,2),1:2])
    ang <- acos(y_ref %*% c(0,1) / sqrt(sum(y_ref^2)))
    ang <- if_else(y_ref[1] < 0, -ang, ang)
    rotM <- matrix(c(cos(ang), sin(ang),0,
                     -sin(ang), cos(ang), 0,
                     0, 0, 1), ncol = 3, byrow = T)
    xyz_ref_med_pc_rot <- xyz_ref_med_pc %*% rotM
    
    # T4
    tar <- T4_dend[[LL]][[j]]
    ind_D = match(tar$tags$`dendrite start` , tar$d$PointNo)
    df_D <-  tar$d[ind_D,]
    root_xyz <- tar$d[ind_D, c("X","Y","Z")]
    
    # - find the subtree with root = dendrite start
    targ <- as.ngraph(tar)
    ii_root <- ind_D
    # subtree and Strahler order
    sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
    subtree <- subset(tar, sub_points) 
    subtree_g <- as.ngraph(subtree, weights = T)
    subtree_so <- strahler_order(subtree) # Strahler order
    max(subtree_so$segments)
    
    PD <- rbind(com_xyz[j,], v1[j,], v0[j,]) #center and arrow
    OD <- rbind(v3[j,], v4[j,]) #OD
    
    # transform to pc coord, xyz_ref_med[1,] as origin
    subtree_pc <- subtree
    subtree_pc$d[, c("X","Y","Z")] <- sweep(xyzmatrix(subtree$d), 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM
    root_xyz_pc <- as.numeric(root_xyz - xyz_ref_med[1,]) %*% pc$rotation %*% rotM 
    PD_pc <- sweep(PD, 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM #center and arrow
    OD_pc <- sweep(OD, 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM #center and arrow
    nb_pc <- sweep(med_xyz[ii2,], 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM # nb col
    
    
    # PLOT
    # windows(width = 9, height = 9)
    # png(filename = paste("SO_T4", letters[LL], "_", j, ".png", sep = ''))
    # pdf(paste("SO_T4", letters[LL], "_", j, ".pdf", sep = ''), width = 9, height = 12)
    plot(nb_pc[1:7, 1:2], col = 'gray50', pch=1, lwd=2, cex= 10, asp= 1, xlab='', ylab ='',xaxt="n", yaxt="n",
         main = "", xlim = c(-20e3, 10e3), ylim = c(-15e3, 15e3))
    axis(1, at = seq(-20e3, 10e3, by = 5e3), labels = seq(-20,10,by=5), las=1)
    # axis(1, at = seq(-20e3, 10e3, by = 5e3), labels = c('-20','-15','-10','-5','0 μm', '5','10'), las=1)
    axis(2, at = seq(-15e3, 15e3, by = 5e3), labels = seq(-15,15,by=5), las=1)
    par("mai"=c(0.2,0.2,0.2,0.2), "omi"=c(0.2,0.2,0.2,0.2))
    for (k in 1:max(subtree_so$points)) {
      pt_p <- subtree_pc$d[,"PointNo"] %in% unique(subtree_pc$d[subtree_so$points == k,"Parent"])
      pt_so <- pt_p | subtree_so$points==k
      if (sum(pt_so) > 1) {
        plot(subset(subtree_pc, pt_so), col=pal_so[k], add = T, lwd = 0.8*k+ 1, WithNodes = F)
      }
    }
    points(root_xyz_pc[1], root_xyz_pc[2], pch=16, cex=2, col='black')
    arrows(PD_pc[2,1], PD_pc[2,2], PD_pc[3,1], PD_pc[3,2], lwd= 3,)
    arrows(OD_pc[1,1], OD_pc[1,2], OD_pc[2,1], OD_pc[2,2], lwd= 3, code = 3, angle = 90, length = 0.2)
    # title(paste("SO_T4", letters[LL], "_#", j, sep = ''))
    if (j == 1) {
      # segments(0,-11000,5000,-11000,lwd = 2) # 10 um
      # text(x=2500, y= -11000, cex=0.7, adj=c(0.5,-0.4), labels = "5 μm")
      points(-20e3, 15e3, cex=1, col='magenta', pch=16)
      text(x=-18e3, y= 15e3, cex=0.7, adj=0, labels = "position in eye coordinates")
      points(-20e3, 14e3, pch=16, cex=2, col='black')
      text(x=-18e3, y= 14e3, cex=0.7, adj=0, labels = "root of dendrite")
    }
  
    x0 <- -20e3
    # text(x=x0 - 1e3, y= 15e3, cex=0.7, adj=0, labels = "[μm]")
    text(x=x0, y= 10e3, cex=0.7, adj=0, labels = paste0("T4", letters[LL]," #", j, " skid= ",tar$skid))
    # text(x=x0, y= 4e3, cex=0.7, adj=0, labels = "eye coord")
    text(x=x0, y= 4e3, cex=0.7, adj=0, labels = paste0("elev = ", round(90 - T4_utp[j,'t'], 1),"°"))
    text(x=x0, y= 3e3, cex=0.7, adj=0, labels = paste0("azim = ", round(T4_utp[j,'p'],1),"°"))
    vPD <- diff(PD_pc[2:3,1:2])
    # aPD <- acos(vPD[2] / sqrt(sum(vPD^2))) /pi*180
    # aPD <- if_else(vPD[1] > 0, -aPD, aPD)
    aPD <- 90 - atan2(vPD[2], vPD[1])/pi*180 #clockwise st in eye is counterclockwise
    text(x=x0, y= 2e3, cex=0.7, adj=0, labels = paste0("meridian angle = ", round(aPD, 2),"°"))
    
    text(x=x0, y= 0e3, cex=0.7, adj=0, labels = paste0("|PD| = ", round(sqrt(sum(vPD^2))/1e3, 1)))
    text(x=x0, y= -1e3, cex=0.7, adj=0, labels = paste0("|OD| = ", round(sqrt(sum((v4[j,]-v3[j,])^2))/1e3, 1)))
    text(x=x0, y= -2e3, cex=0.7, adj=0, labels = paste0("|OD| / |PD| = ", round(sqrt(sum((v4[j,]-v3[j,])^2)) / sqrt(sum(vPD^2)), 2)))
    # text(x=x0, y= -4e3, cex=0.7, adj=0, labels = paste("pts area = ", ell_type[[LL]][j, 'A_pt'], sep = ''))
    # text(x=x0, y= -5e3, cex=0.7, adj=0, labels = paste("pi*a*b area = ", ell_type[[LL]][j, 'A_piab'], sep = ''))
  
    text(x=x0, y= -7e3, cex=0.7, adj=0, labels = paste0("SN max = ", so_type[[LL]][j,'so_max']))
    text(x=x0, y= -8e3, cex=0.7, adj=0, labels = "For segments of SN = 1, 2, 3, 4:")
    text(x=x0, y= -9e3, cex=0.7, adj=0, labels = paste0("No. of segments = ", toString(so_type[[LL]][j,c('so_N1','so_N2','so_N3','so_N4')])))
    text(x=x0, y= -10e3, cex=0.7, adj=0, labels = paste0("Path length[μm] = ", toString(round(so_type[[LL]][j,c('so_cL1','so_cL2','so_cL3','so_cL4')]/1e3, 1))))
    text(x=x0, y= -11e3, cex=0.7, adj=0, labels = paste0("Vector sum length [μm] = ", toString(round(so_type[[LL]][j,c('so_vL1','so_vL2','so_vL3','so_vL4')]/1e3, 1))))
    
    # coord on eye
    polygon(bkgd_mer)
    lines(bkgd_eq)
    lines(bkgd_eq_p45)
    lines(bkgd_eq_m45)
    thetaphi <- T4_utp[j, c('t','p')]
    # thetaphi['t'] <- 90 - thetaphi['t'] # back to [0 180]
    points(Mollweide(thetaphi)*Mollweide_mul + Mollweide_ori, cex=1, col='magenta', pch=16)
    
    # if (j == 1) {
    #   text(x= -18e3, y= 13e3, cex=0.7, adj=0, labels = "blue dot = med col, \n
    #      green diamond = center of mass \n
    #      x, y axes = pc1, 2 in nm \n
    #      (0,0) = dendrite start tag")
    # }
  }
  dev.off()  
}




# extra stuff below -----------------------------------------------------

# angular histo of SO vector ----------------------------------------------

for (LL in 1:1) {
  pdf(file = paste("plots/anghist_v2/T4", letters[LL], "_anghist.pdf", sep = ''))
  neu <- T4_dend[[LL]]
  seg_summ <- seg_summ_type[[LL]]
  for (j in 1:length(neu)) {
    ssn <- seg_summ[[j]] # seg_summ
    df <- matrix(ncol = 3, nrow = 0)
    for (ii_h in 1:4) {
      if (nrow(ssn[ssn$so == ii_h,]) > 0 ) {
        hh <- hist(ssn[ssn$so == ii_h,]$ang, breaks = seq(0,180,10), plot = F)
        df <- rbind(df, cbind(hh$mids, hh$density, ii_h))
      }
    }
    df <- data.frame(df)
    colnames(df) <- c('ang', 'freq', 'so')
    plt <- ggplot(df, aes(x=ang, y=freq)) +
      geom_bar(stat='identity', width = 5) +
      # coord_cartesian(ylim = c(0, 1)) +
      scale_x_continuous(breaks = seq(0,180,45), labels = seq(0,180,45)) +
      facet_grid(rows = vars(so)) +
      theme_minimal() 
    # do.call("grid.arrange", plt) 
    print(plt)
  }
  dev.off()
}



# - summary stats
dev.new()
ggplot(aa, aes(x=so, y=L)) +
  geom_point(position = position_nudge())

ggplot(df, aes(x=ang, y=freq)) +
  geom_bar(stat='identity', width = 5) +
  # coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0,180,45), labels = seq(0,180,45)) +
  facet_grid(rows = vars(so)) +
  theme_minimal() 



# T4 seg_summ analysis, seg dir wrt to SO specific PD or v-axis -------------------------
# gallery, distr of seg
# x-axis = PD, normal to med pointing inward = on central meridian northern hemisphere

load('data/T4_gallery.RData')
load('data/hexnb_ind_dist.RData') #  in 'eyemap_lens_med_lop.R'

library(gridExtra)

# main
for (LL in 1:4) {
  # LL <- 2 # choose type
  # pdf(file = paste("plots/segSummary_segPD_T4", letters[LL], ".pdf", sep = ''), onefile = T, width = 12, height = 7)
  pdf(file = paste("plots/segSummary_vaxis_T4", letters[LL], ".pdf", sep = ''), onefile = T, width = 12, height = 7)
  # eye
  com_xyz_eye <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
  T4_utp <- cart2sph2tp(com_xyz_eye) 
  # med
  com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
  v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix()
  v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
  v3 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
  v4 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()
  # eab <- ell_type[[LL]][,c('ea','eb')] %>% as.matrix()
  # seg
  seg_summ_neu <- seg_summ_type[[LL]]
  
  for (j in 1:length(T4_dend[[LL]])) {
    ii <- sweep(ucl_rot_sm, 2, com_xyz_eye[j, ] )^2 %>% rowSums() %>% which.min()
    if (sum(is.na(nb_ind[ii,])) == 0) {
      tar <- T4_dend[[LL]][[j]]
      
      ii19 <- c(nb_ind[nb_ind[ii,],]) %>% unique()
      # v_in <- colMeans(med_xyz) - med_xyz[391,,drop=F]
      # v_in <- t(v_in / sqrt(sum(v_in^2)))
      v_in <- cross3D(diff(med_xyz[nb_ind[ii,c(2,1)],]), diff(med_xyz[nb_ind[ii,c(3,1)],])) 
      v_in <- t(v_in / sqrt(sum(v_in^2)))
      
      # med normal
      pc <- prcomp(med_xyz[na.omit(ii19),])
      if (sum(pc$rotation[,3] * v_in) > 0) {
        pc$rotation <- -pc$rotation
      }
      v_normal <- pc$rotation[,3] # going outwards
      
      # seg
      seg_summ <- seg_summ_neu[[j]] %>%
        select(xd,yd,zd, so) %>%
        mutate(L = sqrt(xd^2+yd^2+zd^2)) %>% #vec length
        mutate(xn= xd/L, yn= yd/L, zn= zd/L) %>%
        as.data.frame()
      # ## ## PD per so
      # PDso <- matrix(ncol = 4, nrow = 4)
      # for (iiso in 1:4) {
      #   ss <- seg_summ[seg_summ$so == iiso,]
      #   if (nrow(ss) > 0) {
      #     # PD
      #     pdseg <- colSums(ss[, c('xd','yd','zd')])
      #     PDso[iiso, ] <- c(pdseg / sqrt(sum(pdseg^2)), sqrt(sum(pdseg^2)))
      #   }
      # }
      ## ## v-axis, update 20220127
      vaxis <- med_xyz[nb_ind[ii, c(7,6)], ] %>% diff()
      vaxis <- c(vaxis / sqrt(sum(vaxis^2)), sqrt(sum(vaxis^2)))
      PDso <- rbind(vaxis, vaxis, vaxis, vaxis)
      
      # rotation and Mollweide
      ls_gg <- list()
      for (iiso in 1:4) {
        v_x <- PDso[iiso, 1:3]
        v_x <- v_x / sqrt(sum(v_x^2))
        v_z <- v_normal - c(v_normal %*% v_x) * v_x
        v_z <- v_z / sqrt(sum(v_z^2))
        
        ss <- seg_summ[seg_summ$so == iiso,]
        if (nrow(ss) > 0) {
          segn <- ss[, c('xn','yn','zn')]
          
          seg_rot <- rot_zx0(segn, zv = v_z, xv = v_x)
          v_normal_rot <- rot_zx0(matrix(v_normal,ncol = 3), zv = v_z, xv = v_x)
          
          # plot Mollweide
          plt_df <- Mollweide(cart2sph2tp(seg_rot)[,c('t','p')])
          plt_df <- as.data.frame(cbind(plt_df, ss$L / PDso[iiso, 4] * 100))
          colnames(plt_df) <- c('xM','yM','vL')
          plt_df2 <- Mollweide(cart2sph2tp(v_normal_rot)[,c('t','p')])
          plt_df2 <- as.data.frame(plt_df2)
          colnames(plt_df2) <- c('xM','yM')
          ls_gg[[iiso]] <- plt_Mo + 
            geom_point(data=plt_df, aes(x=xM,y=yM, size=vL), shape=1) +
            geom_point(data=plt_df2, aes(x=xM,y=yM), shape=2, size=5) +
            labs(title = paste("so= ",iiso,sep=''), x='',y='')
        }
      }
      ls_gg[[1]] <- ls_gg[[1]] + labs(title= paste("T4", letters[LL]," #", j, " skid= ", tar$skid,", so= 1", sep=''))
      
      # PLOT
      # windows(width = 16, height = 9)
      # ggdraw() +
      #   draw_plot(ls_gg[[1]] + theme(legend.justification = "bottom"), x = 0, y = 0.5, width = 0.5, height = 0.45) +
      #   draw_plot(ls_gg[[2]] + theme(axis.title = element_text()), 0.5, 0.5, 0.5, 0.45) +
      #   draw_plot(ls_gg[[3]] + theme(axis.title = element_text()), 0, 0, 0.5, 0.45) +
      #   draw_plot(ls_gg[[4]] + theme(axis.title = element_text()), 0.5, 0, 0.5, 0.45) 
      
      do.call("grid.arrange", ls_gg)
    }
  }
  dev.off()  
}



# T4 seg_summ analysis, per type, seg dir/pos wrt to neuron PD or v-axis -------------------------
# x-axis = PD or v-axis, 
# z-axis = local normal to med pointing distally into medulla
# choose ang or pos

load('data/T4_gallery.RData')
load('data/hexnb_ind_dist.RData') #  in 'eyemap_lens_med_lop.R'

library(gridExtra)


plt_df_type <- list()
plt_df2_type <- list()
nb_counter <- c(0,0,0,0)
for (LL in 1:4) {
  # LL <- 2 # choose type
   # eye
  com_xyz_eye <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
  T4_utp <- cart2sph2tp(com_xyz_eye) 
  # med
  com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
  v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix()
  v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
  # v3 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
  # v4 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()
  # eab <- ell_type[[LL]][,c('ea','eb')] %>% as.matrix()
  # seg
  seg_summ_neu <- seg_summ_type[[LL]]
  
  plt_df_neu <- list()
  plt_df2_neu <- list()
  ang_vaxis_pca <- vector(mode='numeric', length = length(T4_dend[[LL]]))
  
  seg_rot_neu <- list()  
  # dalpha_new <- c()
  
  for (j in 1:length(T4_dend[[LL]])) {
    # tar <- T4_dend[[LL]][[j]]
    ii <- sweep(ucl_rot_sm, 2, com_xyz_eye[j, ] )^2 %>% rowSums() %>% which.min() # hex nb
    # # max diff in local skew angles
    # dalpha <- nb_ind[ii,] %>% nb_ind[.,] %>% c() %>% unique() %>% nb_ang_med[.,2] %>% na.omit() %>% dist() %>% max()
    # dalpha_new <- c(dalpha_new, dalpha)
    
    # if (sum(is.na(nb_ind[ii,])) == 0 & dalpha <15) {
    if (sum(is.na(nb_ind[ii,])) == 0 &
        cart2sph2tp(com_xyz_eye)[j, 't'] < 900 &
        cart2sph2tp(com_xyz_eye)[j, 'p'] > -60) {
    # if (sum(is.na(nb_ind[ii,])) == 0) {
      # distal direction
      v_out <- cross3D(diff(med_xyz[nb_ind[ii,c(1,2)],]), diff(med_xyz[nb_ind[ii,c(1,3)],])) 
      v_out <- t(v_out / sqrt(sum(v_out^2)))
      
      # med normal
      ii19 <- c(nb_ind[nb_ind[ii,],]) %>% unique()
      pc <- prcomp(med_xyz[na.omit(ii19),])
      if (sum(pc$rotation[,3] * v_out) > 0) {
        pc$rotation <- -pc$rotation
      }
      v_normal <- pc$rotation[,3] # going outwards
      
      ### ### seg dir
      seg_summ <- seg_summ_neu[[j]] %>%
        select(xd,yd,zd, so) %>%
        mutate(L = sqrt(xd^2+yd^2+zd^2)) %>% #vec length
        mutate(xn= xd/L, yn= yd/L, zn= zd/L) %>%
        as.data.frame()
      
      # ### ### seg position
      # seg_summ <- seg_summ_neu[[j]] %>% select(x0,y0,z0, so) #keep position + SO
      # ii_root <- which(seg_summ_neu[[j]]$so == max(seg_summ_neu[[j]]$so) & seg_summ_neu[[j]]$pathL == 0) # pathL=0 and max so
      # seg_summ <- sweep(seg_summ, 2, c(unlist(seg_summ[ii_root[1],1:3]),0) ) #subtract the root
      # seg_summ %<>% mutate(L = sqrt(x0^2+y0^2+z0^2)) %>%
      #   mutate(xn= x0/L, yn= y0/L, zn= z0/L) %>%
      #   as.data.frame()
      # seg_summ <- seg_summ[-which(seg_summ_neu[[j]]$pathL == 0), ]
      
      ## ## PD
      PD <- v1[j,] - v0[j,]
      PD <- c(PD/sqrt(sum(PD^2)), sqrt(sum(PD^2)))
      v_ref <- PD[1:3]
      
      ## ## v-axis, update 20220127
      vaxis <- med_xyz[nb_ind[ii, c(7,6)], ] %>% diff()
      vaxis <- vaxis / sqrt(sum(vaxis^2))
      v_ref <- vaxis
      
      v_x <- v_ref
      v_x <- v_x / sqrt(sum(v_x^2))
      v_z <- v_normal - sum(v_normal * v_x) * v_x
      v_z <- v_z / sqrt(sum(v_z^2))
      
      ang_vaxis_pca[j] <- acos(sum(v_normal * v_x))/pi*180
      
      # rotation and Mollweide
      for (iiso in 1:4) {
        # ss <- seg_summ[seg_summ$so == iiso,]
        segn <- seg_summ[seg_summ$so == iiso, c('xn','yn','zn')]
        if (nrow(segn) > 0) {
          seg_rot <- rot_zx0(segn, zv = v_z, xv = v_x)
          v_normal_rot <- rot_zx0(matrix(v_normal,ncol = 3), zv = v_z, xv = v_x)
          
          # offset the y flip in cart2sph2tp
          seg_rot[,2] <- -seg_rot[,2]
          v_normal_rot[2] <- -v_normal_rot[2]
          
          # plot Mollweide
          plt_df <- Mollweide(cart2sph2tp(seg_rot)[,c('t','p')])
          plt_df <- as.data.frame(cbind(plt_df, seg_summ[seg_summ$so == iiso, "L"] / PD[4] * 100))
          colnames(plt_df) <- c('xM','yM','vL')
          plt_df2 <- Mollweide(cart2sph2tp(v_normal_rot)[,c('t','p')])
          plt_df2 <- as.data.frame(plt_df2)
          colnames(plt_df2) <- c('xM','yM')
          
          if (length(plt_df_neu) < iiso) {
            plt_df_neu[[iiso]] <- plt_df
            plt_df2_neu[[iiso]] <- plt_df2
          } else {
            plt_df_neu[[iiso]] <- rbind(plt_df_neu[[iiso]], plt_df)
            plt_df2_neu[[iiso]] <- rbind(plt_df2_neu[[iiso]], plt_df2)
          }
          
          # # all seg_rot in 3D
          # if (length(seg_rot_neu) < iiso) {
          #   seg_rot_neu[[iiso]] <- seg_rot
          # } else {
          #   seg_rot_neu[[iiso]] <- rbind(seg_rot_neu[[iiso]], seg_rot)
          # }
        }
      }
    } else {
      nb_counter[LL] <- nb_counter[LL] + 1
    }
  }
  plt_df_type[[LL]] <- plt_df_neu
  plt_df2_type[[LL]] <- plt_df2_neu
}

# save
# save(plt_df_type, plt_df2_type, file = "data/seg_pos_plt.rda")
# save(plt_df_type, plt_df2_type, file = "data/seg_ang_plt.rda")

# # DEBUG
# windows(width = 12, height = 7)
# plt + geom_point(data=plt_df, aes(x=xM,y=yM), size=1, shape=1, col=j)
# 
# plt <- plt_Mo + 
#   # geom_point(data=plt_df_type[[LL]][[iiso]], aes(x=xM,y=yM, size=vL), shape=1) +
#   geom_point(data=plt_df, aes(x=xM,y=yM), size=1, shape=1) +
#   # stat_bin_2d(data=plt_df_type[[LL]][[iiso]][,1:2], aes(x=xM,y=yM), binwidth = c(0.1, 0.1)) +
#   # geom_raster(data=plt_df_type[[LL]][[iiso]], aes(x=xM,y=yM, fill = density)) +
#   geom_point(data=plt_df2, aes(x=xM,y=yM), shape=2, size=3) +
#   labs(title = paste("so= ",iiso,sep=''), x='',y='')


# PLOT
for (LL in 1:4) {
  ls_gg <- list()
  for (iiso in 1:4) {
    ls_gg[[iiso]] <- plt_Mo + 
      # geom_point(data=plt_df_type[[LL]][[iiso]], aes(x=xM,y=yM, size=vL), shape=1) +
      geom_point(data=plt_df_type[[LL]][[iiso]], aes(x=xM,y=yM), size=1, shape=1) +
      # stat_bin_2d(data=plt_df_type[[LL]][[iiso]][,1:2], aes(x=xM,y=yM), binwidth = c(0.1, 0.1)) +
      # geom_raster(data=plt_df_type[[LL]][[iiso]], aes(x=xM,y=yM, fill = density)) +
      geom_point(data=plt_df2_type[[LL]][[iiso]], aes(x=xM,y=yM), shape=2, size=3) +
      labs(title = paste("so= ",iiso,sep=''), x='',y='')
  }
  
  ls_gg[[1]] <- ls_gg[[1]] + labs(title= paste("T4", letters[LL], ",so= 1", sep=''))
  # ls_gg[[1]] <- ls_gg[[1]] + labs(title= paste("T4", letters[LL], ",so= 1", ", branch position", sep=''))
  
  windows(width = 12, height = 7)
  # pdf(file = paste("plots/segAngle_neuronPD_T4", letters[LL], ".pdf", sep = ''), onefile = T, width = 12, height = 7)
  # pdf(file = paste("plots/segAngle_vaxis_T4", letters[LL], ".pdf", sep = ''), onefile = T, width = 12, height = 7)
  # pdf(file = paste("plots/segPos_vaxis_T4", letters[LL], ".pdf", sep = ''), onefile = T, width = 12, height = 7)
  do.call("grid.arrange", ls_gg)
  # dev.off()
}


# T4 seg_summ analysis, branching angle, new order param?,  ---------------

load('data/T4_gallery.RData')
load('data/hexnb_ind_dist.RData') #  in 'eyemap_lens_med_lop.R'

library(gridExtra)

plt_df_type <- list()
plt_df2_type <- list()
nb_counter <- c(0,0,0,0)
for (LL in 1:4) {
  # LL <- 2 # choose type
  # eye
  com_xyz_eye <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
  T4_utp <- cart2sph2tp(com_xyz_eye) 
  # med
  com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
  v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix()
  v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
  v3 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
  v4 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()
  # eab <- ell_type[[LL]][,c('ea','eb')] %>% as.matrix()
  # seg
  seg_summ_neu <- seg_summ_type[[LL]]
  
  plt_df_neu <- list()
  plt_df2_neu <- list()
  for (j in 1:length(T4_dend[[LL]])) {
    # tar <- T4_dend[[LL]][[j]]
    ii <- sweep(ucl_rot_sm, 2, com_xyz_eye[j, ] )^2 %>% rowSums() %>% which.min() # hex nb
    # if (sum(is.na(nb_ind[ii,])) == 0 & com_xyz_eye[j,3] > cos(60/180*pi)) {
    # if (sum(is.na(nb_ind[ii,])) == 0 & cart2sph2tp(com_xyz_eye)[j, 'p'] < 30) {
    if (sum(is.na(nb_ind[ii,])) == 0) {
      ii19 <- c(nb_ind[nb_ind[ii,],]) %>% unique()
      v_in <- cross3D(diff(med_xyz[nb_ind[ii,c(2,1)],]), diff(med_xyz[nb_ind[ii,c(3,1)],])) 
      v_in <- t(v_in / sqrt(sum(v_in^2)))
      
      # med normal
      pc <- prcomp(med_xyz[na.omit(ii19),])
      if (sum(pc$rotation[,3] * v_in) > 0) {
        pc$rotation <- -pc$rotation
      }
      v_normal <- pc$rotation[,3] # going outwards
      
      ### ### seg dir
      seg_summ <- seg_summ_neu[[j]] %>%
        select(xd,yd,zd, so) %>%
        mutate(L = sqrt(xd^2+yd^2+zd^2)) %>% #vec length
        mutate(xn= xd/L, yn= yd/L, zn= zd/L) %>%
        as.data.frame()
      # ### ### seg position
      # seg_summ <- seg_summ_neu[[j]] %>% select(x0,y0,z0, so)
      # ii_root <- which(seg_summ_neu[[j]]$so == max(seg_summ_neu[[j]]$so) & seg_summ_neu[[j]]$pathL == 0) # pathL=0 and max so
      # seg_summ <- sweep(seg_summ, 2, c(unlist(seg_summ[ii_root[1],1:3]),0) )
      # seg_summ %<>% mutate(L = sqrt(x0^2+y0^2+z0^2)) %>% 
      #   mutate(xn= x0/L, yn= y0/L, zn= z0/L) %>%
      #   as.data.frame()
      # seg_summ <- seg_summ[-which(seg_summ_neu[[j]]$pathL == 0), ]
      
      # ## ## PD
      # PD <- v1[j,] - v0[j,]
      # PD <- c(PD/sqrt(sum(PD^2)), sqrt(sum(PD^2)))
      ## ## v-axis, update 20220127
      vaxis <- med_xyz[nb_ind[ii, c(7,6)], ] %>% diff()
      vaxis <- c(vaxis / sqrt(sum(vaxis^2)), sqrt(sum(vaxis^2)))
      PD <- vaxis
      
      v_x <- PD[1:3]
      v_z <- v_normal - c(v_normal %*% v_x) * v_x
      v_z <- v_z / sqrt(sum(v_z^2))
      v_x <- v_x / sqrt(sum(v_x^2))
      # rotation and Mollweide
      for (iiso in 1:4) {
        # ss <- seg_summ[seg_summ$so == iiso,]
        segn <- seg_summ[seg_summ$so == iiso, c('xn','yn','zn')]
        if (nrow(segn) > 0) {
          seg_rot <- rot_zx0(segn, zv = v_z, xv = v_x)
          v_normal_rot <- rot_zx0(matrix(v_normal,ncol = 3), zv = v_z, xv = v_x)
          
          # plot Mollweide
          plt_df <- Mollweide(cart2sph2tp(seg_rot)[,c('t','p')])
          plt_df <- as.data.frame(cbind(plt_df, seg_summ[seg_summ$so == iiso, "L"] / PD[4] * 100))
          colnames(plt_df) <- c('xM','yM','vL')
          plt_df2 <- Mollweide(cart2sph2tp(v_normal_rot)[,c('t','p')])
          plt_df2 <- as.data.frame(plt_df2)
          colnames(plt_df2) <- c('xM','yM')
          
          if (length(plt_df_neu) < iiso) {
            plt_df_neu[[iiso]] <- plt_df
            plt_df2_neu[[iiso]] <- plt_df2
          } else {
            plt_df_neu[[iiso]] <- rbind(plt_df_neu[[iiso]], plt_df)
            plt_df2_neu[[iiso]] <- rbind(plt_df2_neu[[iiso]], plt_df2)
          }
        }
      }
    } else {
      nb_counter[LL] <- nb_counter[LL] + 1
    }
  }
  plt_df_type[[LL]] <- plt_df_neu
  plt_df2_type[[LL]] <- plt_df2_neu
}

# PLOT
for (LL in 1:4) {
  ls_gg <- list()
  for (iiso in 1:4) {
    ls_gg[[iiso]] <- plt_Mo + 
      # geom_point(data=plt_df_type[[LL]][[iiso]], aes(x=xM,y=yM, size=vL), shape=1) +
      geom_point(data=plt_df_type[[LL]][[iiso]], aes(x=xM,y=yM), size=1, shape=1) +
      # stat_bin_2d(data=plt_df_type[[LL]][[iiso]][,1:2], aes(x=xM,y=yM), binwidth = c(0.1, 0.1)) +
      # geom_raster(data=plt_df_type[[LL]][[iiso]], aes(x=xM,y=yM, fill = density)) +
      geom_point(data=plt_df2_type[[LL]][[iiso]], aes(x=xM,y=yM), shape=2, size=3) +
      labs(title = paste("so= ",iiso,sep=''), x='',y='')
  }
  
  ls_gg[[1]] <- ls_gg[[1]] + labs(title= paste("T4", letters[LL], ",so= 1", sep=''))
  # ls_gg[[1]] <- ls_gg[[1]] + labs(title= paste("T4", letters[LL], ",so= 1", ", branch position", sep=''))
  
  # windows(width = 12, height = 7)
  pdf(file = paste("plots/segAngle_vaxis_T4", letters[LL], ".pdf", sep = ''), onefile = T, width = 12, height = 7)
  do.call("grid.arrange", ls_gg)
  dev.off()
}

# ellipsoid -----------------------------------------------------------------

ell <- bind_rows(ell_type[1:4], .id = "type")

dev.new()
plot(ell$ea/ell$ec, ell$eb/ell$ec, col=pal_T4[as.numeric(ell$type)], cex=1, pch=16)

dev.new()
plot(ell$type, ell$V_piabc, col=pal_T4[as.numeric(ell$type)], cex=1, pch=16)

dev.new()
plot(ell$type, ell$ec, col=pal_T4[as.numeric(ell$type)], cex=1, pch=16)

dev.new()
plot(ell$type, ell$ea/ell$eb,col=pal_T4[as.numeric(ell$type)], cex=1, pch=16)
dev.new()
plot(ell$type, ell$eye_ea/ell$eye_eb,col=pal_T4[as.numeric(ell$type)], cex=1, pch=16)

dev.new()
plot(ell$ea/ell$eb, ell$eye_ea/ell$eye_eb, col=pal_T4[as.numeric(ell$type)], cex=1, pch=16)

# T4 medulla depth --------------------------------------------------------

T4_dep_type <- list()
seg_dep_type <- list()

for (LL in 1:4) {
  T4_dep_neu <- list()
  seg_dep_neu <- list()
  
  neu <- T4_dend[[LL]]
  anno <- anno_T4_dend[[LL]]
  
  for (j in 1:length(neu)) {
    tar <- neu[[j]]
    ind_D = match(tar$tags$`dendrite start` , tar$d$PointNo)
    df_D <-  tar$d[ind_D,]
    
    # - find the subtree with root = dendrite start
    targ <- as.ngraph(tar)
    ii_root <- ind_D
    # subtree and Strahler order
    sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
    subtree <- subset(tar, sub_points) 
    subtree <- resample(subtree, stepsize = 100) #resample
    subtree_g <- as.ngraph(subtree, weights = T)
    subtree_so <- strahler_order(subtree) # Strahler order
    
    nb_N <- 1 + 6 + 12
    ori_xyz <- colMeans(xyzmatrix(subtree$d))
    nbd <- rowSums(sweep(med_xyz, 2, ori_xyz, '-')^2)
    ind_nb_N <- order(nbd)[1:nb_N]
    nb <- med_xyz[ind_nb_N, ]
    
    # - med pca for depth, TODO depth plot ??
    nb_pca <- prcomp(nb) # in rotated coord, x = 1st pc
    pt_prox <- c(350000, 250000, 300000)
    # +z goes distally
    if ( nb_pca$rotation[,3] %*% pt_prox > 0 ) {
      nb_pca$rotation[,3] <- - nb_pca$rotation[,3]
    }
    if (t(cross3D(nb_pca$rotation[,3],nb_pca$rotation[,1])) %*% nb_pca$rotation[,2] < 0 ) {
      # nb_pca$rotation[,1] <- - nb_pca$rotation[,1] # ok too
      nb_pca$rotation[,2] <- - nb_pca$rotation[,2]
    }
    
    # # DEBUG
    # nopen3d()
    # plot3d(subtree, WithNodes = F, lwd=2)
    # arrow3d(ori_xyz, ori_xyz+nb_pca$rotation[,3]*5e3, col ='blue',type = "rotation")
    # points3d(med_xyz)
    
    T4_dep_neu[[j]] <- sweep(xyzmatrix(subtree$d),2,nb_pca$center,'-') %*% nb_pca$rotation[,3]
    
    seg_xyz <- seg_summ_type[[LL]][[j]][,c('x0','y0','z0')] %>% as.matrix()
    seg_dep_neu[[j]] <- sweep(seg_xyz, 2, nb_pca$center,'-') %*% nb_pca$rotation[,3]
    # root_xyz <- dir_type[[LL]][j, c('rtx','rty','rtz')] %>% as.matrix()
    # seg_dep_neu[[j]] <- sweep(seg_xyz, 2, root_xyz, '-') %*% nb_pca$rotation[,3]
  }
  T4_dep_type[[LL]] <- T4_dep_neu
  seg_dep_type[[LL]] <- seg_dep_neu
}

# PLOT --> not much structure 
dev.new()
plot(seg_summ_type[[LL]][[j]]$pathL, seg_dep_type[[LL]][[j]])
seg_xyz <- seg_summ_type[[LL]][[j]][,c('x0','y0','z0')] %>% as.matrix()
root_xyz <- dir_type[[LL]][j, c('rtx','rty','rtz')] %>% as.matrix()
vecL <- sweep(seg_xyz, 2, root_xyz, '-')^2 %>% rowSums() %>% sqrt()
points(vecL, seg_dep_type[[LL]][[j]], pch = 3, col='red')

# x <- c()
# y <- c()
for (LL in 1:4) {
  anno <- anno_T4_dend[[LL]]
  x <- c()
  y <- c()
  z <- c()
  for (j in 1:nrow(anno)) {
    x <- c(x, seg_summ_type[[LL]][[j]]$pathL)
    y <- c(y, seg_dep_type[[LL]][[j]])
    seg_xyz <- seg_summ_type[[LL]][[j]][,c('x0','y0','z0')] %>% as.matrix()
    root_xyz <- dir_type[[LL]][j, c('rtx','rty','rtz')] %>% as.matrix()
    z <- c(z, sweep(seg_xyz, 2, root_xyz, '-')^2 %>% rowSums() %>% sqrt())
  }
  # points(x, y, pch=3, col=LL)
  points(z, y, pch=1, col=LL)
}
plot(z,y,type = 'n')



# T4 tree morphology -----------------------------------------------------------


# - compare tree structure

load("data/neu_T4_overlap.RData")
T4_overlap_loc <- c('M', 'LD', 'LV', 'VP')
T4_overlap_ind <- rbind(c(1, T4_overlap_N[1]),
                        c(T4_overlap_N[1] + 1, sum(T4_overlap_N[1:2])),
                        c(sum(T4_overlap_N[1:2]) + 1, sum(T4_overlap_N[1:3])),
                        c(sum(T4_overlap_N[1:3]) + 1, sum(T4_overlap_N[1:4])) )
ii_overlap <- 2
anno <- anno_overlap[T4_overlap_ind[ii_overlap,1]:T4_overlap_ind[ii_overlap,2],]
neu <- T4_overlap[T4_overlap_ind[ii_overlap,1]:T4_overlap_ind[ii_overlap,2]]
  

# LL <- 2
# neu <- T4_dend[[LL]]
# anno <- anno_T4_dend[[LL]]


j <- 2

tar <- neu[[j]]
ind_D = match(tar$tags$`dendrite start` , tar$d$PointNo)

df_D = tar$d[ind_D,]
xyz_D <- xyzmatrix(df_D)
# xyz_D_avg <- colMeans(xyz_D) #for multiple dendrite tags
# node_dend <- child_node(tar, df_D)
# dend_xyz <- xyzmatrix(node_dend)
# dend_com <- colMeans(dend_xyz)

# - find the subtree with root = dendrite start (NB use its parent to avoid plot segments with single node)
targ <- as.ngraph(tar)
ii_root <- ind_D

# ii_root <- 1145

# sub_points <- igraph::graph.dfs(targ, root = ii_root[br], unreachable=FALSE, neimode='out')$order
sub_points <- igraph::graph.bfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
subtree <- subset(tar, sub_points)
subtree <- resample(subtree, stepsize = 400)
subtree_so <- strahler_order(subtree) # Strahler order
subtree_g <- as.ngraph(subtree, weights = T)

# add so to edge
E(subtree_g)$so <- NA
for (ii_e in 1:length(E(subtree_g))) {
  evn <- attributes(E(subtree_g)[ii_e])$vnames
  ii <- regexpr("\\|",  evn)[1]
  pn <- substring(evn, ii+1, nchar(evn))
  ii <- match(pn, subtree$d$PointNo)
  E(subtree_g)$so[ii_e] <- subtree_so$points[ii]
}

# remove non-vertex points
subtree_g2 <- subtree_g 
for (ii_seg in 1:length(subtree$SegList)) {
  seg_ind <- subtree$SegList[[ii_seg]]
  if (length(seg_ind) > 2) {
    segend_ind <- seg_ind[c(1, length(seg_ind))]
    seg <- subtree$d$PointNo[seg_ind] %>% as.character() # PointNo
    segend <- seg[c(1, length(seg))] #ends
    segmid <- seg[-c(1,length(seg))] #middle
    #cable dist
    seg_len <- igraph::distances(subtree_g, v= segend_ind[1], to= segend_ind[2]) 
    subtree_g2 <- subtree_g2 + edge(segend, weight= seg_len, so= subtree_so$segments[ii_seg]) - vertices(segmid)
  }
}
E(subtree_g2)$weight <- round(E(subtree_g2)$weight, 0)

# E(subtree_g2)[[]]


# nodes <- xyzmatrix(subtree$d)
g1 <- subtree_g
windows(width = 10, height = 12)
# ggraph(g1, layout = 'dendrogram', circular = F) +
ggraph(g1, layout = 'tree') +
  # geom_edge_diagonal(aes(label = seq(1,ecount(g1)), colour = as.factor(so)), edge_width = 1) +
  geom_edge_diagonal(aes(colour = as.factor(so)), edge_width = 2) +
  scale_edge_colour_manual(name = "SO", values = pal_so ) +
  geom_node_point() +
  # geom_node_text(aes(label = seq(1,vcount(g1))), nudge_x = 0.1 ) +
  theme_void()


# igraph::distances(g1, v = 18, to = 33 )
# g2 <- g1 + edge(segend, weight = seg_len) - vertices(segmid)
g2 <- subtree_g2
windows(width = 10, height = 8)
# ggraph(g2, layout = 'tree') +
#   geom_edge_diagonal(aes(colour = as.factor(so)), edge_width = 2) +
ggraph(g2, layout = 'dendrogram',  circular = F) +
  # geom_edge_elbow(aes(colour = weight), edge_width = 2) +
  # geom_edge_elbow(aes(label = seq(1,ecount(g2)), colour = as.factor(so)), edge_width = 2) +
  geom_edge_elbow(aes(colour = as.factor(so)), edge_width = 2) +
  scale_edge_colour_manual(name = "SO", values = pal_so ) +
  geom_node_point() +
  # geom_node_text(aes(label = seq(1,vcount(g2))), nudge_x = 0.1 ) +
  theme_void()

# # DEBUG
# nopen3d()
# # plot3d(subtree, col = 'grey')
# for (so in 1:max(subtree_so$points)) {
#   pt_p <- subtree$d$PointNo %in% unique(subtree$d[subtree_so$points==so,"Parent"])
#   pt_so <- pt_p | subtree_so$points==so
#   plot3d(subset(subtree, pt_so), col=pal_so[so], add = T, boundingbox = boundingbox(subtree), lwd= 3, WithNodes = F)
# }

# # find longest path across graph
# d=get.diameter(gw)
# # make a new neuron using the longest path
# gw_spine=as.neuron(induced.subgraph(gw, d))
# # make a new neuron containing all nodes except those in longest path
# gw_antispine=as.neuron(delete.vertices(gw, d))


# summary -------------------------------------------------------------------------------------------------------
so_plt <- lacroix

# df <- cbind(T4_utp[, c('t','p')], ell_area_pt, T4_colnum, T4_colnum_para, ell_e, so_N, so_vL, so_cL)
# df[is.na(df)] <- 0
# colnames(df) <- c('t','p', 'A', 'colNum', 'colNumPara', 'e',
#                   'N1','N2','N3','N4',
#                   'L1','L2','L3','L4',
#                   'C1','C2','C3','C4')
# dfsort <- df[order(df$t), ]


df <- cbind(T4_utp[, c('t','p')], ell_area_pt, ell_abc[,1:2],ell_e, ell_ab_eye[,1:2], ell_e_eye, so_N, so_vL, so_cL)
df[is.na(df)] <- 0
colnames(df) <- c('t','p', 'A', 'ea', 'eb', 'ellip', 'ea_eye', 'eb_eye', 'ellip_eye',
                  'N1','N2','N3','N4',
                  'L1','L2','L3','L4',
                  'C1','C2','C3','C4')
dfsort <- df[order(df$t), ]


windows(width = 20, height = 5)
ggplot(dfsort) + 
  geom_ribbon(aes(x = t, ymin = 0, ymax = L1, fill = "SO 1")) +
  geom_ribbon(aes(x = t, ymin = L1, ymax = L1+L2, fill = "SO 2")) +
  geom_ribbon(aes(x = t, ymin = L1+L2, ymax = L1+L2+L3, fill = "SO 3")) +
  geom_ribbon(aes(x = t, ymin = L1+L2+L3, ymax = L1+L2+L3+L4, fill = "SO 4")) +
  scale_color_manual(name = "SO", aesthetics = "fill", values = unname(so_plt[1:4]) ) +
  xlab('elev [deg]') +
  ylab('SO length sum [nm]') +
  labs(title = paste("summary SO vec length", sep = ',')) 

windows(width = 20, height = 5)
ggplot(dfsort) + 
  geom_ribbon(aes(x = t, ymin = 0, ymax = C1, fill = "SO 1")) +
  geom_ribbon(aes(x = t, ymin = C1, ymax = C1+C2, fill = "SO 2")) +
  geom_ribbon(aes(x = t, ymin = C1+C2, ymax = C1+C2+C3, fill = "SO 3")) +
  geom_ribbon(aes(x = t, ymin = C1+C2+C3, ymax = C1+C2+C3+C4, fill = "SO 4")) +
  scale_color_manual(name = "SO", aesthetics = "fill", values = unname(so_plt[1:4]) ) +
  xlab('elev [deg]') +
  ylab('SO length sum [nm]') +
  labs(title = paste("summary cable ength", sep = ',')) 

windows(width = 20, height = 5)
ggplot(dfsort, aes(x = t)) + 
  geom_ribbon(aes(ymin = 0, ymax = N1, fill = "SO 1")) +
  geom_ribbon(aes(ymin = N1, ymax = N1+N2, fill = "SO 2")) +
  geom_ribbon(aes(ymin = N1+N2, ymax = N1+N2+N3, fill = "SO 3")) +
  geom_ribbon(aes(ymin = N1+N2+N3, ymax = N1+N2+N3+N4, fill = "SO 4")) +
  scale_color_manual(name = "SO", aesthetics = "fill", values = unname(so_plt[1:4]) ) +
  xlab('elev [deg]') +
  ylab('SO number') +
  labs(title = paste("summary SO number", sep = ',')) 
  
dev.new()
ggplot(dfsort) +
  geom_point(aes(x = t, y = A), color = 'red') +
  xlab('elev [deg]') +
  ylab('SO number') +
  labs(title = paste("summary area", sep = ',')) 
  
dev.new()
ggplot(dfsort) +
  geom_point(aes(x = t, y = ellip), color = 'blue') +
  xlab('elev [deg]') +
  ylab('SO number') +
  labs(title = paste("summary ellipticity", sep = ',')) 

dev.new()
ggplot(dfsort) +
  geom_point(aes(x = t, y = ellip_eye), color = 'blue') +
  xlab('elev [deg]') +
  ylab('SO number') +
  labs(title = paste("summary eye ellipticity", sep = ',')) 

dev.new()
ggplot(dfsort) +
  geom_point(aes(x = t, y = colNum), color = 'red') +
  xlab('elev [deg]') +
  ylab('med col number') +
  labs(title = paste("med col num", sep = ',')) 

dev.new()
ggplot(dfsort) +
  geom_point(aes(x = t, y = colNumPara), color = 'red') +
  xlab('elev [deg]') +
  ylab('med col number') +
  labs(title = paste("med col num Parallel", sep = ',')) 


dfsort <- df[order(df$p), ]
dev.new()
ggplot(dfsort) +
  geom_point(aes(x = p, y = colNumPara), color = 'red') +
  xlab('azim [deg]') +
  ylab('med col number') +
  labs(title = paste("med col num Parallel", sep = ',')) 


# -- heatmap for num of col along preferred dir
# library(alphahull)
# 
# # transform st. looking out the right eye, [t=-45,p=45] = top right
# load("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_eyemap/Lens_Medulla/eyemap_20200319.RData")
# colnames(ucl_rot) <- c('x','y','z')
# ucl_rot %<>% as_tibble() %>%  
#   mutate(y = -y) %>%
#   mutate(theta = acos(z)) %>%
#   mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
#   mutate(t = theta / pi * 180, p = phi/pi*180) %>% 
#   mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
#   as.data.frame()
# 
# xyMollweide <- Mollweide(ucl_rot[,c('t', 'p')])
# colnames(xyMollweide) <- c('xM','yM')
# ucl_rot <- cbind(ucl_rot, xyMollweide)
# 
# ucl_rot$t <- ucl_rot$t - 90
# 
# ii_excl <- c(754, 745, 733)
# # points(ucl_rot[-ii_excl,c('xM','yM')])
# ucl_rot_excl <- ucl_rot[-ii_excl,]
# 
# ucl_ashape <- ashape(ucl_rot_excl[,c('xM','yM')] + matrix(runif(dim(ucl_rot_excl)[1]*2, 1e-9, 2e-9), ncol = 2), alpha = 0.5)
# poly_xy <- mkpoly(ucl_ashape$edges)[[1]]
# poly_xy <- poly_xy[,c(3,4)]
# poly_xy <- rbind(poly_xy, poly_xy[1,])
# poly_xy <- as.data.frame(poly_xy)
# colnames(poly_xy) <- c('x','y')
# 
# 
# bd_grid <- expand.grid(seq(), bd_theta) # equirectangular
# 
# # PLOT,  Gaussian around each com with synap# as height,  cp data via binning
# ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,'phi_deg'], xy_bd_chull[,'theta_deg'])


mg <- df[, c('t','p','colNumPara')]
# gplt <- ggplot(data = mg, aes(x = p, y = t)) +
#   geom_contour(aes(z = colNumPara)) +
#   # scale_fill_distiller(palette = "Blues", direction = 1) +
#   scale_y_reverse(expand = c(0, 0)) +
#   coord_fixed(ratio = 1) +
#   xlab('azim [deg]') +
#   ylab('elev [deg]') +
#   labs(title = "num of medulla col along the preferred direction")

gplt <- ggplot(data = mg, aes(x = p, y = t)) +
  geom_point() +
  scale_y_reverse() +
  coord_fixed(ratio = 1) +
  xlab('azim [deg]') +
  ylab('elev [deg]') +
  labs(title = "num of medulla col along the preferred direction")
for (j in 1:dim(mg)[1]) {
  gplt <- gplt + 
    annotate("text", x = df[j,'p'] + 2, y = df[j,'t'] +2 , label = df[j,'colNumPara'])
}

dev.new()
gplt


# gallery local T4 overlap --------------------------------------------------


load("data/neu_T4_overlap.RData")

rownames(Mi1_xyz) <- seq(1, nrow(Mi1_xyz))
med_xyz <- Mi1_xyz[CT1_Mi1[eyemap[,1]],]

ucl_rot_gal <- ucl_rot_sm
colnames(ucl_rot_gal) <- c('x','y','z')
ucl_rot_gal %<>% as_tibble() %>%  
  mutate(y = -y) %>%
  mutate(theta = acos(z)) %>%
  mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
  mutate(t = theta / pi * 180 - 90, p = phi/pi*180) %>% 
  mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
  as.data.frame()

pal_so <- brewer.pal(6, "RdYlBu")[c(5,6,3,2,1)]

T4_overlap_loc <- c('M', 'LD', 'LV', 'VP')
T4_overlap_ind <- rbind(c(1, T4_overlap_N[1]),
                        c(T4_overlap_N[1] + 1, sum(T4_overlap_N[1:2])),
                        c(sum(T4_overlap_N[1:2]) + 1, sum(T4_overlap_N[1:3])),
                        c(sum(T4_overlap_N[1:3]) + 1, sum(T4_overlap_N[1:4])) )

for (ii_overlap in 1:4) {
  anno <- anno_overlap[T4_overlap_ind[ii_overlap,1]:T4_overlap_ind[ii_overlap,2],]
  neu <- T4_overlap[T4_overlap_ind[ii_overlap,1]:T4_overlap_ind[ii_overlap,2]]
  N <- length(neu)
  
  # backgroud cols
  ii <- sweep(med_xyz, 2, colMeans(xyzmatrix(neu)))^2 %>% rowSums() %>% order() %>% head(1+6+12)
  bkgd_col <- med_xyz[ii,]
  
  # pca 
  # arrow3d(colMeans(med_xyz), colMeans(med_xyz) + c(-0.84,  0.20, -0.49)*1e5, theta = pi/6,n = 4, col="blue", type = "rotation")
  med_pca <- prcomp(bkgd_col)
  if (med_pca$rotation[,1] %*% c(-0.84,  0.20, -0.49) < 0) {
    med_pca$rotation <- - med_pca$rotation
  }
  if (t(cross3D(med_pca$rotation[,1],med_pca$rotation[,2])) %*% med_pca$rotation[,3] < 0 ) {
    med_pca$rotation[,3] <- - med_pca$rotation[,3]
  }
  bkgd_col_pc <- sweep(bkgd_col, 2, med_pca$center) %*% med_pca$rotation
  
  # # DEBUG
  # nopen3d()
  # points3d(med_xyz[1:Npt,], col = "grey", size = 5)
  # points3d(bkgd_col, col = "blue", size = 10)
  # for (j in 1:N) {
  #   ii <- regexpr("T4", anno[j,"name"]) #start index
  #   type_ii <- match(substring(anno[j,"name"], ii+2, ii+2), letters)
  #   arrow3d(RF_CT1[j,4:6], RF_CT1[j,1:3], theta = pi / 12, n = 4, col = T4_col[type_ii],alpha = 0.9, type = "rotation", lit=F) 
  #   text3d(RF_CT1[j,1:3], texts = paste(j), adj = -1)
  # }
  # arrow3d(med_pca$center, med_pca$center + med_pca$rotation[,1]*1e5, theta = pi/6,n = 4, col="blue", type = "rotation")
  # # CT1
  # # points3d(med_xyz[CT1_Mi1[as.numeric(rownames(med_xyz))] %in% bm[[1]],], size = 10, col = 'cyan')
  # # points3d(med_xyz[CT1_Mi1[as.numeric(rownames(med_xyz))] %in% bb0,], size = 10, col = 'cyan')
  # # Mi1
  # points3d(med_xyz[rownames(med_xyz) %in% bm[[1]],], size = 10, col = 'cyan')
  # points3d(med_xyz[rownames(med_xyz) %in% bm[[3]],], size = 10, col = 'cyan')
  # points3d(med_xyz[rownames(med_xyz) %in% bb0,], size = 10, col = 'cyan')
  
  
  # grid for area calculation, physical unit
  grid_x <- seq(-1.5e4, 1.5e4, by = 1e3)
  grid_y <- seq(-1.5e4, 1.5e4, by = 1e3)
  bkgd_grid <- expand.grid(grid_x, grid_y)
  # grid for area calculation, for counting columns
  grid_x <- seq(-4e4, 4e4, by = 1e2)
  grid_y <- seq(-4e4, 4e4, by = 1e2)
  bkgd_grid_col <- expand.grid(grid_x, grid_y) # 0.1 um
  
  # - loop and save
  
  # Np <- 1 # neuron per page
  # pdf(file = paste("plots/T4_overlap_", T4_overlap_loc[ii_overlap], "_SO23.pdf", sep = '') )
  
  so_N <- matrix(0, ncol = 4, nrow = length(neu)) # num of SO branches, order 1 to 4
  so_vL <- matrix(0, ncol = 4, nrow = length(neu)) # vector length of SO branches
  so_cL <- matrix(0, ncol = 4, nrow = length(neu)) # cable length of SO branches
  so_v <- list() # list of T4 dir based on SO
  so_v_raw <- list() # list of T4 dir based on SO, no normalization
  so_max <- c() # max SO
  
  ell_abc <- matrix(ncol = 3, nrow = length(neu)) # major and minor axes along pca, and subtype
  ell_area_pt <- c() #num of points
  ell_area_piab <- c() # pi*ea*eb
  ell_e <- c()
  T4_dir <- matrix(ncol = 3, nrow = length(neu)) # SO direction
  # T4_utp <- as.data.frame(matrix(ncol = 5, nrow = length(neu))) # ucl_rot cloest to dendrite com
  # colnames(T4_utp) <- c('x','y','z','t','p')
  
  ell_ab_eye <- matrix(ncol = 3, nrow = length(neu)) # major and minor axes along pca, and subtype
  ell_e_eye <- c()
  
  T4_chull_med <- list()
  T4_chull_lens <- list()
  
  for (j in 1:length(neu)) {
    tar <- neu[[j]]
    ind_D = match(tar$tags$`dendrite start`, tar$d$PointNo)
    df_D = tar$d[ind_D,]
    xyz_D <- xyzmatrix(df_D)
    xyz_D_avg <- colMeans(xyz_D) #for multiple dendrite tags
    node_dend <- child_node(tar, df_D)
    dend_xyz <- xyzmatrix(node_dend)
    dend_com <- colMeans(dend_xyz)
    
    # find the subtree with root = dendrite start (NB use its parent to avoid plot segments with single node)
    targ <- as.ngraph(tar)
    ii_root <- ind_D
    
    # in case of multiple branches
    sub_points <- c()
    so_max_br <- c()
    subtree_br_so <- list() # branch SO list
    subtree_br <- list() #branch subtree list
    subtree_br_g <- list()
    for (br in 1:length(ii_root)) {
      sub_points_br <- igraph::graph.dfs(targ, root = ii_root[br], unreachable=FALSE, neimode='out')$order
      sub_points <- c(sub_points, sub_points_br) # all sub nodes
      
      subtree_br[[br]] <- subset(tar, sub_points_br) 
      subtree_br_so[[br]] <- strahler_order(subtree_br[[br]]) # Strahler order
      subtree_br_g[[br]] <- as.ngraph(subtree_br[[br]], weights = T)
      so_max_br <- c(so_max_br, max(subtree_br_so[[br]]$points) )
    }
    subtree <- subset(tar, sub_points)
    subtree_g <- as.ngraph(subtree, weights = T)
    so_max <- c(so_max, max(so_max_br) )
    
    # -- define orientation based on strahler order
    so_v_mat <- matrix(0, ncol = 3, nrow = 4)
    so_v_mat_raw <- matrix(0, ncol = 3, nrow = 4) #no normalization
    for (so in 1:4) {
      segnum <- 0
      seg_v <- matrix(ncol = 3, nrow = 0)
      for (br in 1:length(ii_root)) {
        seg_ls <- subtree_br[[br]]$SegList[subtree_br_so[[br]]$segments == so]
        if (length(seg_ls) > 0) {
          # add up so segments
          for (k in 1:length(seg_ls)) {
            node_ii <- seg_ls[[k]]
            if (length(node_ii) > 1) {
              seg_v <- rbind(seg_v, subtree_br[[br]]$d[tail(node_ii,1), c("X","Y","Z")] - subtree_br[[br]]$d[head(node_ii,1), c("X","Y","Z")])
              so_cL[j, so] <- so_cL[j, so] + igraph::distances(subtree_br_g[[br]], v = node_ii[1], to = tail(node_ii, 1) )
            }
          }
        } 
        segnum <- segnum + length(seg_ls)
      }
      
      so_N[j, so] <- segnum
      so_vL[j, so] <- sqrt(sum(colSums(seg_v)^2))
      so_v_mat[so, ] <- as.numeric(colSums(seg_v)) / sqrt(sum(colSums(seg_v)^2))
      so_v_mat_raw[so, ] <- as.numeric(colSums(seg_v))
    }
    
    # -- translation (xyz_D or dend_com as origin) then rotation 
    # ori_xyz <- xyz_D_avg
    # ori_xyz <- dend_com
    ori_xyz <- med_pca$center
    nodes <- xyzmatrix(subtree$d)
    nodes_pc <-  sweep(nodes, 2, med_pca$center, '-') %*% med_pca$rotation
    subtree_pc <- subtree
    subtree_pc$d[, c("X","Y","Z")] <- nodes_pc
    subtree_br_pc <- subtree_br #for branches
    for (br in 1:length(ii_root)) {
      nodes <- xyzmatrix(subtree_br[[br]]$d)
      nodes_pc <-  sweep(nodes, 2, med_pca$center, '-') %*% med_pca$rotation
      subtree_br_pc[[br]]$d[, c("X","Y","Z")] <- nodes_pc
    }
    xyz_D_pc <- sweep(xyz_D, 2, ori_xyz, '-') %*% med_pca$rotation
    dend_com_pc <- (dend_com - ori_xyz) %*% med_pca$rotation
    so_v_mat <- so_v_mat %*% med_pca$rotation
    so_v_mat_raw <- so_v_mat_raw %*% med_pca$rotation
    
    so_v[[j]] <- so_v_mat
    so_v_raw[[j]] <- so_v_mat_raw
    T4_dir[j,] <- (so_v_mat[2,] + so_v_mat[3,]) / 2 #weight so = 2,3 equally
    
    # -- find the position on the eye
    nb_pc <- bkgd_col_pc
    nb_col_pc <- nb_pc
    
    # # find the eye lens
    # ind_nnb <- as.numeric(rownames(nb)[1])
    # ind_nnb_col <- as.numeric(rownames(nb_col)[1])
    # T4_utp[j, ] <- as.numeric(ucl_rot[eyemap[ind_nnb,2], c('x','y','z', 't', 'p')])
    
    # -- ellipticity 
    # ellipticity
    hpts_ind <- chull(xyzmatrix(subtree_pc$d)[,1:2])
    hpts_ind <- c(hpts_ind, hpts_ind[1])
    hpts_xy <- xyzmatrix(subtree_pc$d)[hpts_ind, 1:2] # hull edge points
    hpts_xy_med <- hpts_xy
    poly_st <- st_polygon(list(hpts_xy))
    # line_st_pc1 <- st_linestring(matrix(c(-2e4, dend_com_pc[2], 2e4, dend_com_pc[2]), nrow = 2, byrow = T)) # x-axis
    line_st_pc1 <- st_linestring(matrix(c(dend_com_pc[1]-T4_dir[j,1]*2e5, dend_com_pc[2]-T4_dir[j,2]*2e5, dend_com_pc[1]+T4_dir[j,1]*2e5, dend_com_pc[2]+T4_dir[j,2]*2e5), nrow = 2, byrow = T)) # x-axis
    int_st_pc1 = st_intersection(line_st_pc1, poly_st) # intersections
    # line_st_pc2 <- st_linestring(matrix(c(dend_com_pc[1], -2e4, dend_com_pc[1], 2e4), nrow = 2, byrow = T))
    line_st_pc2 <- st_linestring(matrix(c(dend_com_pc[1]+T4_dir[j,2]*2e5, dend_com_pc[2]-T4_dir[j,1]*2e5, dend_com_pc[1]-T4_dir[j,2]*2e5, dend_com_pc[2]+T4_dir[j,1]*2e5), nrow = 2, byrow = T)) # x-axis
    int_st_pc2 = st_intersection(line_st_pc2, poly_st) #  intersection
    
    
    # intersections of major and minor axes
    ea <- sqrt(sum(as.numeric(diff(int_st_pc1))^2)) / 2 / 1000 #[um]
    eb <- sqrt(sum(as.numeric(diff(int_st_pc2))^2)) / 2 / 1000
    
    # subtype
    ii <- regexpr("T4", anno[j,"name"]) #start index
    type_ii <- match(substring(anno[j,"name"], ii+2, ii+2), letters)
    
    # ellipticity 
    ell_abc[j,] <- c(ea,eb,type_ii)
    ell_e <- c(ell_e, (ea-eb)/ea)
    
    # - count grid points
    hpts_ind_col <- chull(nb_col_pc[, 1:2])
    hpts_ind_col <- c(hpts_ind_col, hpts_ind_col[1])
    hpts_xy_col <- nb_col_pc[hpts_ind_col, 1:2] # hull edge points
    poly_st_col <- st_polygon(list(hpts_xy_col))
    in_bkgd_col <- sp::point.in.polygon(bkgd_grid_col[,1], bkgd_grid_col[,2], hpts_xy_col[,1], hpts_xy_col[,2])
    in_bkgd <- sp::point.in.polygon(bkgd_grid_col[,1], bkgd_grid_col[,2], hpts_xy[,1], hpts_xy[,2])
    
    # - T4 area
    # in_bkgd <- sp::point.in.polygon(bkgd_grid[,1], bkgd_grid[,2], hpts_xy[,1], hpts_xy[,2])
    ell_area_pt <- c(ell_area_pt, round(sum(in_bkgd)/100, digits = 1) )
    ell_area_piab <- c(ell_area_piab, round(pi*ea*eb, digits = 1) )
    
    # - map onto eye
    # find nearby columns
    nbhd <- 20000
    ii <- sqrt(rowSums(sweep(med_xyz, 2, dend_com)^2)) <= nbhd
    nb_ID <- which(ii)
    
    #axes pts
    axpt_pc <- rbind(cbind(as.matrix(int_st_pc1)*1.5,0), cbind(as.matrix(int_st_pc2)*1.5,0)) 
    axpt <- sweep(axpt_pc %*% t(med_pca$rotation), 2, med_pca$center, '+')
    
    # np
    xyz <- rbind(axpt, xyzmatrix(subtree$d)[hpts_ind, ])
    T4_chull_med[[j]] <- xyz
    xyz_chull <- matrix(ncol = 3, nrow = nrow(xyz))
    xyz_chull_eval <- data.frame(mc.x = xyz[,1], mc.y = xyz[,2], mc.z = xyz[,3])
    # -- kernel regression
    for (k in 1:3) {
      npdata <- data.frame(mc = med_xyz[nb_ID,], ec = ucl_rot_gal[nb_ID,k]) 
      model_np <- npreg(ec ~ mc.x + mc.y + mc.z, data = npdata)
      for (m in 1:nrow(xyz_chull_eval)) {
        xyz_chull[m, k] <- predict(model_np, newdata = xyz_chull_eval[m,])
      }
    }
    T4_chull_lens[[j]] <- xyz_chull
    
    # ea and eb on eye
    eye_pca <- prcomp(ucl_rot_gal[nb_ID, 1:3])
    if (eye_pca$rotation[,3] %*% c(1,1,1) < 0) {
      eye_pca$rotation <- - eye_pca$rotation
    }
    if (t(cross3D(eye_pca$rotation[,3],eye_pca$rotation[,1])) %*% eye_pca$rotation[,2] < 0 ) {
      eye_pca$rotation[,2] <- - eye_pca$rotation[,2]
    }
    
    xyz_chull_pc <- sweep(xyz_chull, 2, eye_pca$center) %*% eye_pca$rotation
    hpts_ind <- chull(xyz_chull_pc[-c(1,2,3,4),1:2])
    hpts_ind <- c(hpts_ind, hpts_ind[1])
    hpts_xy <- xyz_chull_pc[-c(1,2,3,4),1:2][hpts_ind, ] # hull edge points
    poly_st <- st_polygon(list(hpts_xy))
    line_st_pc1 <- st_linestring(xyz_chull_pc[1:2,1:2]) # ea-axis
    int_st_pc1 = st_intersection(line_st_pc1, poly_st) # intersections
    line_st_pc2 <- st_linestring(xyz_chull_pc[3:4,1:2])  # eb
    int_st_pc2 = st_intersection(line_st_pc2, poly_st) #  intersection
    
    # intersections of major and minor axes on unit circle eye
    ea_eye <- sqrt(sum(as.numeric(diff(int_st_pc1))^2)) / 2
    eb_eye <- sqrt(sum(as.numeric(diff(int_st_pc2))^2)) / 2
    # ellipticity 
    ell_ab_eye[j,] <- c(ea_eye,eb_eye,type_ii)
    ell_e_eye <- c(ell_e_eye, (ea_eye - eb_eye)/ea_eye)
    
#   }
# }

    # PLOT
    # dev.new()
    # plot(nb_pc[, 1:2], col = 'blue', pch = 1, cex = 3)
    # for (k in 1:max(subtree_br_so$points)) {
    #   pt_p <- subtree_pc$d[,"PointNo"] %in% unique(subtree_pc$d[subtree_br_so$points == k,"Parent"])
    #   pt_so <- pt_p | subtree_br_so$points==k
    #   # plot(subset(subtree_pc, pt_so), col=pal_so[k], add = k!=1, lwd = 2)
    #   plot(subset(subtree_pc, pt_so), col=pal_so[k], add = T, lwd = 2)
    # }
    # points(matrix(dend_com_pc[1:2], ncol = 2), col = 'green', pch = 18, cex = 2)
    # points(matrix(xyz_D_pc[,1:2], ncol = 2), col = 'black', pch = 16, cex = 1)
    # # points(matrix(T4_dir[j, 1:2]*5000, ncol = 2), col = 'cyan', pch = 16, cex = 2)
    # # lines(rbind(xyz_D_pc[,1:2], T4_dir[j, 1:2]*10000), lwd = 2)
    # arrows(0,0, T4_dir[j,1]*1.5e4, T4_dir[j,2]*1.5e4, lwd = 2)
    # polygon(hpts_xy)
    
    # pdf PLOT
    # if (j %% Np == 1) {
    plot(nb_pc[, 1:2], col = 'blue', pch = 16, cex = 1, asp = 1, xlab = '', ylab = '',
         main = "", xlim = c(-20e3, 15e3), ylim = c(-15e3, 15e3*(2*Np-1)))
    for (br in 1:length(ii_root)) {
      subtree_plt <- subtree_br_pc[[br]]
      for (k in 1:max(subtree_br_so[[br]]$points)) {
        pt_p <- subtree_plt$d[,"PointNo"] %in% unique(subtree_plt$d[subtree_br_so[[br]]$points == k,"Parent"])
        pt_so <- pt_p | subtree_br_so[[br]]$points==k
        if (sum(pt_so) > 1) {
          plot(subset(subtree_plt, pt_so), col=pal_so[k], add = T, lwd = 2, WithNodes = F)
        }
      }
    }
    
    points(matrix(dend_com_pc[1:2], ncol = 2), col = 'green', pch = 18, cex = 2)
    points(matrix(xyz_D_pc[,1:2], ncol = 2), col = 'brown', pch = 16, cex = 1)
    # single arrow
    arrows(0,0, T4_dir[j,1]*2e4, T4_dir[j,2]*2e4, lwd = 1)
    # # all SO arrows
    # for (ii_so in 1:4) {
    #   arrows(0,0, so_v_raw[[j]][ii_so,1], so_v_raw[[j]][ii_so,2], lwd = 1, col = pal_so[ii_so])
    # }
    polygon(hpts_xy_med)
    # } else {
    #   subtree_plt$d[,"Y"] <- subtree_plt$d[,"Y"] + 40e3*((j-1)%%Np)
    #   for (k in 1:max(subtree_br_so$points)) {
    #     pt_p <- subtree_plt$d[,"PointNo"] %in% unique(subtree_plt$d[subtree_br_so$points == k,"Parent"])
    #     pt_so <- pt_p | subtree_br_so$points==k
    #     plot(subset(subtree_plt, pt_so), col=pal_so[k], add = T, lwd = 2, axes = F, xlab = '', ylab = '')
    #   }
    # }
    
    ii <- regexpr("T4", anno[j,"name"]) #start index
    type_ii <- substring(anno[j,"name"], ii+2, ii+2)
    
    x0 <- -20e3
    text(x= x0, y= 40e3*((j-1)%%Np)+ 8e3, cex=0.7, adj=0,
         labels = paste("T4", type_ii, ", ", tar$NeuronName, ", ", j,sep = ''))
    text(x= x0, y= 40e3*((j-1)%%Np)+ 4e3, cex=0.7, adj=0, labels = "eye coord")
    # text(x= x0, y= 40e3*((j-1)%%Np)+ 3e3, cex=0.7, adj=0, labels = paste("elev = ", round(T4_utp[j,'t'],1), sep = ''))
    # text(x= x0, y= 40e3*((j-1)%%Np)+ 2e3, cex=0.7, adj=0, labels = paste("azim = ", round(T4_utp[j,'p'],1), sep = ''))
    
    # text(x= x0, y= 40e3*((j-1)%%Np)+ 0e3, cex=0.7, adj=0, labels = paste("local col area = ", col_area[j], sep = ''))
    # text(x= x0, y= 40e3*((j-1)%%Np)+ -1e3, cex=0.7, adj=0, labels = paste("pts area = ", ell_area_pt[j], sep = ''))
    text(x= x0, y= 40e3*((j-1)%%Np)+ -2e3, cex=0.7, adj=0, labels = paste("pi*a*b area = ", ell_area_piab[j], sep = ''))
    
    # text(x= x0, y= 40e3*((j-1)%%Np)+ -4e3, cex=0.7, adj=0, labels = paste("col num = ", T4_colnum[j], sep = ''))
    # text(x= x0, y= 40e3*((j-1)%%Np)+ -5e3, cex=0.7, adj=0, labels = paste("col num para = ", T4_colnum_para[j], sep = ''))
    
    text(x= x0, y= 40e3*((j-1)%%Np)+ -5e3, cex=0.7, adj=0, labels = paste("ea = ", round(ea,2), sep = ''))
    text(x= x0, y= 40e3*((j-1)%%Np)+ -6e3, cex=0.7, adj=0, labels = paste("eb = ", round(eb,2), sep = ''))
    text(x= x0, y= 40e3*((j-1)%%Np)+ -7e3, cex=0.7, adj=0, labels = paste("ellip = ", round(ell_e[j],2), sep = ''))
    # text(x= x0, y= 40e3*((j-1)%%Np)+ -8e3, cex=0.7, adj=0, labels = paste("inter col dist = ", L_inter_col[j], sep = ''))
    
    text(x= x0, y= 40e3*((j-1)%%Np)+ -8e3, cex=0.7, adj=0, labels = paste("ea eye = ", round(ea_eye,2), sep = ''))
    text(x= x0, y= 40e3*((j-1)%%Np)+ -9e3, cex=0.7, adj=0, labels = paste("eb eye = ", round(eb_eye,2), sep = ''))
    text(x= x0, y= 40e3*((j-1)%%Np)+ -10e3, cex=0.7, adj=0, labels = paste("ellip eye = ", round(ell_e_eye[j],2), sep = ''))
    
    
    text(x= x0, y= 40e3*((j-1)%%Np)+ -12e3, cex=0.7, adj=0, labels = paste("SO max = ", so_max[j], sep = ''))
    text(x= x0, y= 40e3*((j-1)%%Np)+ -13e3, cex=0.7, adj=0, labels = paste("SO num = ", toString(so_N[j,]), sep = ''))
    text(x= x0, y= 40e3*((j-1)%%Np)+ -14e3, cex=0.7, adj=0, labels = paste("vector length = ", toString(round(so_vL[j,]/1e3, 1)), sep = ''))
    text(x= x0, y= 40e3*((j-1)%%Np)+ -15e3, cex=0.7, adj=0, labels = paste("cable length = ", toString(round(so_cL[j,]/1e3, 1)), sep = ''))
    
    
    if (j == 1) {
      text(x= -18e3, y= 13e3, cex=0.7, adj=0, labels = "blue dot = med col, \n
         green diamond = center of mass \n
         x, y axes = pc1, 2 in nm \n
         (0,0) = dendrite start tag")
    }
  }
  dev.off()

} #ii_overlap