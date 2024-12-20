# generate T4b gallery and save to pdf ----------------------------------------
# choose nb via nb_ind
# align to meridian by mapping ref points to med via np

library(natverse)
library(tidyverse)

source("eyemap_func.R")

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
  pdf(file = paste("plots/gallery_T4", letters[LL], ".pdf", sep = ''))
  
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
    plot(nb_pc[1:7, 1:2], col = 'gray50', pch=1, lwd=2, cex= 10, asp= 1, xlab='', ylab ='',xaxt="n", yaxt="n",
         main = "", xlim = c(-20e3, 10e3), ylim = c(-15e3, 15e3))
    axis(1, at = seq(-20e3, 10e3, by = 5e3), labels = paste0(seq(-20,10,by=5), 'μm'), las=1)
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
    if (j == 1) {
      points(-20e3, 15e3, cex=1, col='magenta', pch=16)
      text(x=-18e3, y= 15e3, cex=0.7, adj=0, labels = "position in eye coordinates")
      points(-20e3, 14e3, pch=16, cex=2, col='black')
      text(x=-18e3, y= 14e3, cex=0.7, adj=0, labels = "root of dendrite")
    }
    
    x0 <- -20e3
    text(x=x0, y= 10e3, cex=0.7, adj=0, labels = paste0("T4", letters[LL]," #", j, " skid= ",tar$skid))
    text(x=x0, y= 4e3, cex=0.7, adj=0, labels = paste0("elev = ", round(90 - T4_utp[j,'t'], 1),"°"))
    text(x=x0, y= 3e3, cex=0.7, adj=0, labels = paste0("azim = ", round(T4_utp[j,'p'],1),"°"))
    vPD <- diff(PD_pc[2:3,1:2])
    aPD <- 90 - atan2(vPD[2], vPD[1])/pi*180 #clockwise st in eye is counterclockwise
    text(x=x0, y= 2e3, cex=0.7, adj=0, labels = paste0("meridian angle = ", round(aPD, 2),"°"))
    
    text(x=x0, y= 0e3, cex=0.7, adj=0, labels = paste0("|PD| = ", round(sqrt(sum(vPD^2))/1e3, 1)))
    text(x=x0, y= -1e3, cex=0.7, adj=0, labels = paste0("|width| = ", round(sqrt(sum((v4[j,]-v3[j,])^2))/1e3, 1)))
    text(x=x0, y= -2e3, cex=0.7, adj=0, labels = paste0("|width| / |PD| = ", round(sqrt(sum((v4[j,]-v3[j,])^2)) / sqrt(sum(vPD^2)), 2)))
    text(x=x0, y= -5e3, cex=0.7, adj=0, labels = paste("pi*a*b area = ", ell_type[[LL]][j, 'A_piab'], sep = ''))
    
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
    points(Mollweide(thetaphi)*Mollweide_mul + Mollweide_ori, cex=1, col='magenta', pch=16)
  }
  dev.off()  
}


