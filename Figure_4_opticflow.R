# Fig.4A, visual -- neuronal ------------------------------------------------------

# # find eg
# com <-  sapply(T4_dend[[2]], function(x) colMeans(xyzmatrix(x))) %>% t()
# nopen3d()
# points3d(med_xyz)
# points3d(med_xyz[vaxis_gen(clv,ixy = ind_xy),], col=pal_axes[3],size =12)
# points3d(med_xyz[ind_xy[ind_xy[,2] == clp,1], ], col=pal_axes[2],size =12)
# points3d(med_xyz[ind_xy[ind_xy[,3] == clq,1], ], col=pal_axes[1],size =12)
# points3d(com, size=12, col='red')

LL <- 2
ind_eg2 <- c(29, 142, 163,  139, 169, 121,   65, 143,  170) #posterior
ind_eg2 <- c(86,  168,   3) #anterior
ind_eg2 <- c(162, 87, 171, 143) #south

mrot <- matrix(c(cos(60/180*pi), sin(60/180*pi), 0,
                 -sin(60/180*pi), cos(60/180*pi), 0,
                 0, 0, 1), ncol = 3, byrow = T)
ind_axis <- vaxis_gen(clv) 

# eye
ucl_rot_rhs_shift <- ucl_rot_right %*% mrot 
ucl_rot_rhs_shift <- sweep(ucl_rot_rhs_shift, 2, c(-10,0,0))
ucl_rot_sm_shift <- ucl_rot_sm %*% mrot 
ucl_rot_sm_shift <- sweep(ucl_rot_sm_shift, 2, c(-10,0,0))

# med
xyz_data <- Mi1_M10_xyz
colnames(xyz_data) <- c('x','y','z')
xyz_data2 <- xyz_data^2
Y <- rowSums(xyz_data2)
X <- cbind(2*xyz_data, rep(1,dim(xyz_data)[1]))
sph <- lsfit(X,Y,intercept = FALSE)
r_fit <- sqrt(sph[[1]][4]+sum(sph[[1]][1:3]^2)) *0.9 #radius
cosph <- sph[["coefficients"]][c('x','y','z')] #center of the sphere

med_M10_shift <- sweep(Mi1_M10_xyz,2,cosph)/r_fit
med_xyz_shift <- sweep(med_xyz,2,cosph)/r_fit

pt_theta90 <- rbind(c(0,0,0), med_M10_shift[Mi1_ind_PR[[8]],])
lm_plane <- lm(pt_theta90[,1] ~ pt_theta90[,2] + pt_theta90[,3] + 0)
vn <- c(-1, lm_plane$coefficients[1], lm_plane$coefficients[2]) #norm
vn <- vn/sqrt(sum(vn^2))
if (sum(vn * med_M10_shift[pt_up_Mi1, ]) < 0) {
  vn <- -vn
}
vx <- med_M10_shift[pt_100_Mi1, ]
med_M10_shift <- rot_zx0(med_M10_shift, vn, vx)
med_xyz_shift <- rot_zx0(med_xyz_shift, vn, vx)


# PLOT
bb <- rbind(c(0.06303674, -1.11337304,-0.94910240),
            c(11.43733883,   1.05451393,1.61241937) )
colnames(bb) <- c('x','y','z')

nopen3d()
# par3d('windowRect' = c(100,100,2100,1500))
par3d('windowRect' = c(50,50,2100,1200))

# par3d('userMatrix')
M <- matrix(c(0.19098200, -0.98094434, 0.03567285, -0.080707046,
              0.00933612,  0.03815538, 0.99922776,  0.290239177,
              -0.98154825, -0.19050151, 0.01644530,  0.004202841,
              0.00000000,  0.00000000, 0.00000000,  1.000000000),
            ncol = 4, byrow = T)

plot3d(boundingbox(bb), alpha=0)
rgl.viewpoint(fov=0,zoom=0.1762879, userMatrix= M ) #inside-out


# - med col
points3d(med_xyz_shift[vaxis_gen(clv,ixy = ind_xy),], col=pal_axes[3],size =12)
# points3d(med_M10_shift[eyemap[match(ind_axis, eyemap[,2]),1],], col=pal_axes[3],size =12) #alt
points3d(med_xyz_shift[ind_xy[ind_xy[,2] == clp,1], ], col=pal_axes[2],size =12)
points3d(med_xyz_shift[ind_xy[ind_xy[,3] == clq,1], ], col=pal_axes[1],size =12)
points3d(med_xyz_shift[haxis_gen(clh, ind_xy),], col=pal_axes[4],size =12)
points3d(med_M10_shift[seq(1,nrow(utp_Mi1_rot_chiasm)) %in% lens_Mi1[,2],], size = 8, col = 'gray')

# - lens xyz
points3d(ucl_rot_rhs_shift[ind_axis,], size=12, col=pal_axes[3])
points3d(ucl_rot_rhs_shift[lens_ixy[lens_ixy[,2] == 1,1], ], size=12, col=pal_axes[2])
points3d(ucl_rot_rhs_shift[lens_ixy[lens_ixy[,3] == 0,1], ], size=12, col=pal_axes[1])
points3d(ucl_rot_rhs_shift[haxis_gen(1, lens_ixy),], size=12, col=pal_axes[4])
points3d(ucl_rot_rhs_shift[seq(1,nrow(ucl_rot_rhs_shift)) %in% lens_Mi1[,1],], size = 8, col = 'gray')

# - T4
# tar <- T4_dend[[2]][[ind_eg2[6]]] #posterior
tar <- T4_dend[[2]][[ind_eg2[2]]] #anterior
ind_T = match(tar$tags$"dendrite start" , tar$d$PointNo)
targ <- as.ngraph(tar)
ii_root <- ind_T
# subtree and resample
sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order

sk <- tar
sk <- subset(tar, sub_points)

xyz <- xyzmatrix(sk$d)
xyz <- sweep(xyz,2,cosph)/r_fit
xyz <- rot_zx0(xyz, vn, vx)
sk$d[, c("X","Y","Z")] <- xyz
# 
# plot3d(sk, lwd=3, col=pal_T4[2], soma = F, WithNodes = F) # dendrite

# PD
# for (j in ind_eg2[c(3,6,9)]) {
# for (j in ind_eg2[c(6)]) {
for (j in ind_eg2[1]) {
  # -- med
  vcom <- as.matrix(dir_type[[LL]][j, c("comx","comy","comz")])
  v0 <- as.matrix(dir_type[[LL]][j, c("rsx0","rsy0","rsz0")])
  vd <- as.matrix(dir_type[[LL]][j, c("rsxd","rsyd","rszd")])
  # nbhd_N <- 1+8 +16 # num of nbs
  nbhd_N <- 1+6 +12 # num of nbs
  # ii_nb <- sweep(Mi1_M10_xyz, 2, vcom, '-')^2 %>% rowSums() %>% order() %>% head(nbhd_N)
  # pch3d(med_M10_shift[ii_nb,],col=pal_T4[LL], pch=1, cex=1)
  ii_nb <- sweep(med_xyz, 2, vcom, '-')^2 %>% rowSums() %>% order() %>% head(nbhd_N)
  pch3d(med_xyz_shift[ii_nb,],col=pal_T4[LL], pch=1, cex=1)
  
  v_shift <- rbind(v0, vd) 
  v_shift <- sweep(v_shift,2,cosph)/r_fit
  v_shift <- rot_zx0(v_shift, vn, vx )
  arrow3d(v_shift[2,], v_shift[1,], theta = pi/6, n = 8, col =pal_T4[LL], type = "rotation", lit=F)
  
  # -- eye
  vcom <- as.matrix(lens_type[[LL]][j, c("comx","comy","comz")])
  v0 <- as.matrix(lens_type[[LL]][j, c("x0","y0","z0")])
  vd <- as.matrix(lens_type[[LL]][j, c("xd","yd","zd")])
  
  vcom_shift <- vcom %*% mrot +c(10,0,0)
  v0_shift <- v0 %*% mrot +c(10,0,0)
  vd_shift <- vd %*% mrot +c(10,0,0)

  pch3d(ucl_rot_sm_shift[ii_nb,],col=pal_T4[LL], pch=1, cex= 0.8)
  arrow3d(vd_shift , v0_shift, theta = pi/12, n = 8, col =pal_T4[LL], type = "rotation", lit=F)
}


# rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(-60/180*pi,1,0,0) %*% rotationMatrix(-35/180*pi,0,0,1)) #outside
rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(-90/180*pi,1,0,0) %*% rotationMatrix(125/180*pi,0,0,1)) #inside-out

# - connecting lines
linemat <- matrix(t(cbind(med_M10_shift[lens_Mi1[match(ind_axis, lens_Mi1[,1]),2],],
                          ucl_rot_rhs_shift[ind_axis,])), ncol = 3, byrow = T)
segments3d(linemat, color = 'gray', lwd=3, alpha=0.5)

# - exclu
pch3d(ucl_rot_rhs_shift[!(seq(1,nrow(ucl_rot_rhs_shift)) %in% lens_Mi1[,1]),], pch=1, cex=0.1, lwd=1.5, col ='gray30')
pch3d(med_M10_shift[!(seq(1,nrow(utp_Mi1_rot_chiasm)) %in% lens_Mi1[,2]),,drop=F], pch=1, cex=0.25, lwd=1.5, col ='gray30')

aa <- med_M10_shift[!(seq(1,nrow(utp_Mi1_rot_chiasm)) %in% lens_Mi1[,2]),,drop=F]
# pch3d(aa, pch=3, cex=0.25, lwd=1.5, col =c('gray30','white')) # ??some problem with pch3d
pch3d(rbind(aa,aa*5), pch=3, cex=0.25, lwd=1.5, col =c('gray30','white')) #this works...

# rgl.snapshot(filename = paste("T4b_med_eye_excl", ".png", sep = ''))
# rgl.snapshot(filename = paste("T4b_med_eye", ".png", sep = ''))

# check np regression -----------------------------------------------------

# see eyemap_np_test.R

# Fig.4B,  Mollweide, real T4 --------------------------------------------------------------

## ## chose type, 2 or 4
LL <- 2 

# eye
com_xyz <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
v0 <- as.matrix(lens_type[[LL]][, c("x0","y0","z0")])
v1 <- as.matrix(lens_type[[LL]][, c("xd","yd","zd")])

v0_Mo <- sweep(v0, 1, sqrt(rowSums(v0^2)), '/')
v0_Mo <- cart2sph2tp(v0_Mo)
v0_Mo <- Mollweide(v0_Mo[,c('t', 'p')])
colnames(v0_Mo) <- c('xM','yM')
v1_Mo <- sweep(v1, 1, sqrt(rowSums(v1^2)), '/')
v1_Mo <- cart2sph2tp(v1_Mo)
v1_Mo <- Mollweide(v1_Mo[,c('t', 'p')])
colnames(v1_Mo) <- c('xM','yM')

# PLOT
df_pos <- data.frame(ucl_rot_Mo)
df_arrow <- data.frame(cbind(v0_Mo[,c('xM','yM')], v1_Mo[,c('xM','yM')]) )
colnames(df_arrow) <- c('x','y','xend','yend')

plt <- plt_Momin +
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  geom_point(data=df_pos, aes(x = xM, y = yM), colour='gray')
if (LL == 2) {
  plt <- plt +
    geom_segment(data = df_arrow[ind_eg2[c(1)],], aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =3) +
    geom_segment(data = df_arrow[-ind_eg2[c(1)],], aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =1)
} else {
  # geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[2],arrow = arrow(length = unit(0.008, "npc"), type = "closed"), size =1) +
  plt <- plt +
    geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =1)
}
  
# add v-axis
# ii_ax <- -14:14
# for (j in ii_ax ) {
#   ind_axis <- vaxis_gen(j)
#   vl <- df_pos[match(ind_axis, eyemap[,2]),]
#   vl <- vl[order(vl$yM),]
#   plt <- plt +
#     geom_path(data = vl, aes(x=xM, y=yM), colour = 'grey50', lwd=0.5)
# }
windows(width = 12, height = 8)
# pdf(paste("T4", letters[LL], "_RF_Mollweide_real.pdf", sep = ''), width = 9.5, height = 6.5)
plt
# dev.off()
# ggsave("T4d_RF_Mollweide.pdf", width = 8.5, height = 4.5)

# Fig.4C, Mollweide with indi and avg H2 with edge ----------------------------------

## ## choose T4 type
LL = 2

vv <- ucl_rot_sm + 0.5*(RF_lens_T4_pred[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)
# vv <- ucl_rot_sm + 0.5*(RF_lens_T4_pred_sm[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)
vv_Mo <- sweep(vv, 1, sqrt(rowSums(vv^2)), '/')
vv_Mo <- cart2sph2tp(vv_Mo)
vv_Mo <- Mollweide(vv_Mo[,c('t', 'p')])
colnames(vv_Mo) <- c('xM','yM')
rownames(vv_Mo) <- rownames(ucl_rot_sm)

# H2 
# combine arenaAng and headAnd
uxy <- unique(tb[, c('stimPosX', 'stimPosY')])

df <- matrix(ncol = 8, nrow = 0)
df_indi <- matrix(ncol = 8, nrow = 0)
for (j in 1:nrow(uxy)) {
  for (k in c(0,1)) {
    ii <- tb$stimPosX == uxy$stimPosX[j] &
      tb$stimPosY == uxy$stimPosY[j] &
      tb$edgeVal == k
    
    thetaphi <- tb_v[ii, 1:2] #base
    
    if (sum(is.na(thetaphi[,1])) < nrow(thetaphi)) {
      thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
      thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
      thetaphi <- thetaphi / 180*pi
      thetaphi <- cbind(1, thetaphi)
      xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
      xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
      
      pxyz <- xyz
      
      ## ## spk
      thetaphi <- tb_v[ii,3:4] 
      # ## ## sub-threshold
      # thetaphi <- tb_v[ii,5:6] 
      
      thetaphi[, 1] <- 90 - thetaphi[, 1]
      thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
      thetaphi <- thetaphi / 180*pi
      thetaphi <- cbind(1, thetaphi)
      xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
      xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
      
      dxyz <- xyz
      
      df_tmp2 <- cbind(pxyz, dxyz, j, k)
      
      df_tmp <- c(colMeans(pxyz),
                  colMeans(dxyz), j, k)
    }
    df <- rbind(df, df_tmp)
    df_indi <- rbind(df_indi, df_tmp2)
  }
}
colnames(df) <- c('x','y','z','xend','yend','zend','pos', 'edge')
colnames(df_indi) <- c('x','y','z','xend','yend','zend','pos', 'edge')
 
# average over birght and dark edge
df_avg <- as.data.frame(df) %>%
  group_by(pos) %>%
  summarise(x=mean(x),y=mean(y),z=mean(z),
            xend=mean(xend),yend=mean(yend), zend=mean(zend)) %>%
  ungroup() %>%
  as.data.frame()

# dark (=0) or bright (=1)
# pxyz <- df[df[,8]== 1, 1:3] 
pxyz <- df_avg[,c('x','y','z')] #combine
H2p_Merc <- cart2Mercator(pxyz) #Mercator
H2p_Mo <- cart2sph2tp(pxyz) # mollweide
H2p_Mo <- Mollweide(H2p_Mo[,c('t', 'p')])
colnames(H2p_Mo) <- c('xM','yM')

# dxyz <- df[df[,8]== 1, 4:6]
dxyz <- df_avg[, c('xend','yend','zend')]
H2d_Merc <- cart2Mercator(dxyz)
H2d_Mo <- cart2sph2tp(dxyz) 
H2d_Mo <- Mollweide(H2d_Mo[,c('t', 'p')])
colnames(H2d_Mo) <- c('xM','yM')

# PLOT
# Mollweide
df_H2 <- cbind(H2p_Mo, H2d_Mo) %>% as.data.frame()
df_H2 <- df_H2[!is.na(df_H2[,1]),]
colnames(df_H2) <- c('x','y','xend','yend')
df_H2 <- as.data.frame(df_H2)

df_pos <- data.frame(ucl_rot_Mo)
df_arrow <- data.frame(cbind(ucl_rot_Mo[,c('xM','yM')], vv_Mo[,c('xM','yM')]) )
colnames(df_arrow) <- c('x','y','xend','yend')

plt <- plt_Momin +
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  # geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],arrow = arrow(length = unit(0.008, "npc"), type = "closed"), size =1) +
  geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =1) +
  annotate("text", x = df_H2$x+0.02, y = df_H2$y+0.02, label = 1:nrow(df_H2), size = 3) +
  geom_segment(data = df_H2, aes(x = x,y = y, xend = xend,yend = yend), colour='red',arrow = arrow(length = unit(0.02, "npc"), type = "closed"), size =1)

windows(width = 12, height = 8)
# pdf(paste0("T4", letters[LL], "_RF_H2_edge_Moll.pdf"), width = 9.5, height = 6.5)
plt
# dev.off()
# ggsave("T4d_RF_Mollweide.pdf", width = 8.5, height = 4.5)


# ED Fig.6A, Mercator projection ---------------------------------------------

LL = 2
vv <- ucl_rot_sm + 0.5*(RF_lens_T4_pred_sm[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)

vv_Mer <- cart2Mercator(vv)

# H2
df_H2 <- cbind(H2p_Merc, H2d_Merc) %>% as.data.frame()
df_H2 <- df_H2[!is.na(df_H2[,1]),]
colnames(df_H2) <- c('x','y','xend','yend')
df_H2 <- as.data.frame(df_H2)

# # add a v-axis
# ind_axis <- vaxis_gen(5)
# xy <- data.frame(ucl_rot_sm)[match(ind_axis, eyemap[,2]),]
# ax_Mer <- cart2Mercator(xy)

# PLOT
df_arrow <- data.frame(cbind(ucl_rot_Merc, vv_Mer) )
colnames(df_arrow) <- c('x','y','xend','yend')


plt <- plt_Mer +
  geom_path(data = cmer_Merc, aes(x=x, y=y), colour = pal_axes[3], lwd=1) +
  # geom_path(data = ax_Mer, aes(x=x, y=y), colour = 'blue', lwd=1) +
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =1) +
  geom_segment(data = df_H2, aes(x = x,y = y, xend = xend,yend = yend), colour='red', size=1) +
  # geom_segment(data = df_H2, aes(x = x,y = y, xend = xend,yend = yend), colour='red',arrow = arrow(length = unit(0.02, "npc"), type = "closed"), size =1)
  # geom_point(data=df_arrow[ind_eghex,], aes(x=x, y=y), size=3) + #eg hex position
  labs(title = "")

# add hex
# for (k in 1:3) {
#   # j <- c(510, 662, 366)[k]
#   # j <- c(500, 554, 322)[k]
#   j <- ind_eghex[k]
#   pt <- ucl_rot_Merc[nb_ind[j,], ]
#   # pt <- pt[c(2,3,6,4,5,7,2), ] #2021
#   pt <- pt[c(2,3,4,5,6,7,2), ] #2023
#   plt <- plt + geom_path(data=as.data.frame(pt), aes(x=x, y=y), colour='gray50',lwd=0.5)
# }

windows(width = 5, height = 6.5)
plt
# ggsave(paste("T4", letters[LL], "_RF_Mercator_PD.pdf",sep=''), width = 5, height = 6.5)


# cont. Fig.4C, summary of comparison and explanaition of bias -----------------------------------------------------------------

kk <- c()
ang <- matrix(ncol = 2, nrow = 6)
for (j in seq(1,6)) {
  m <- sweep(as.matrix(df_arrow[,c('x','y')]), 2, as.matrix(df_H2[j,c('x','y')]))
  k <- rowSums(m^2) %>% which.min() #closet T4
  kk <- c(kk, k)
  ang[j,] <- c(atan2(df_H2$yend[j] - df_H2$y[j], 
                     df_H2$xend[j] - df_H2$x[j]),
               atan2(df_arrow$yend[k] - df_arrow$y[k], 
                     df_arrow$xend[k] - df_arrow$x[k])
  )
}

ang <- ang[c(1,4,6,5,3,2), ]

# PLOT
# df <- data.frame(da = ang[,1] - ang[,2])
# df <- df /pi*180
# windows(width = 6, height = 3)
# # pdf(paste("H2_edge_DvsBvsG", '.pdf', sep = ''), width = 10, height = 5)
# ggplot(df) + 
#   geom_point(aes(x = seq(1,6), y=da), size = 1) +
#   coord_cartesian(ylim = c(-15, 15)) +
#   # scale_y_continuous(breaks= seq(-75,50,by=25), labels = paste0(seq(-75, 50,by=25), "°")) +
#   theme_minimal() +
#   labs(title= paste("H2 tuning", sep = ''),
#        x='position', y='Angle wrt +h [deg]') 


# --- plot both angs
df <- data.frame(
  rbind(
    cbind(ang[,1]/pi*180, 1), cbind(ang[,2]/pi*180,2)
  )
)
colnames(df) <- c('ang','type')
df$ang <- df$ang + 180

windows(width = 6, height = 3)
# pdf(paste("H2_vs_T4b", '.pdf', sep = ''), width = 6, height = 3)
ggplot(df ) + 
  geom_point(aes(x=c(seq(1,6),seq(1,6)), y=ang, colour= factor(type)), size = 1) +
  scale_discrete_manual(values = c('red', pal_T4[2]),aesthetics = "colour",guide= guide_legend(title="")) +
  coord_cartesian(ylim = c(0, 40)) +
  scale_y_continuous(breaks= seq(0,40,by=10), labels = paste0(seq(0,40,by=10), "°")) +
  scale_x_continuous(breaks= seq(0,6,by=1), labels = seq(0,6,by=1)) +
  theme_minimal() +
  labs(title= paste("H2 vs T4b", sep = ''),
       x='position', y='') 
# dev.off()

# -- indi
pxyz <- df_indi[,c('x','y','z')] 
H2p_Merc_indi <- cart2Mercator(pxyz) #Mercator
dxyz <- df_indi[, c('xend','yend','zend')]
H2d_Merc_indi <- cart2Mercator(dxyz)

df_H2_indi <- cbind(H2p_Merc_indi, H2d_Merc_indi, df_indi[,'pos']) %>% as.data.frame()
# df_H2_indi <- df_H2_indi[!is.na(df_H2_indi[,1]),]
colnames(df_H2_indi) <- c('x','y','xend','yend','pos')
df_H2_indi <- as.data.frame(df_H2_indi)

ang_indi <- matrix(ncol = 2, nrow = nrow(df_H2_indi))
for (j in 1:nrow(ang_indi)) {
  ang_indi[j,] <- c(atan2(df_H2_indi$yend[j] - df_H2_indi$y[j], 
                          df_H2_indi$xend[j] - df_H2_indi$x[j]),
                    df_H2_indi$pos[j]
  )
}

# --- plot both angs
df <- cbind(ang_indi, 1)
df <- as.data.frame(df)
colnames(df) <- c('ang','pos','type')
df$x <- 1
# reorder, c(1,4,6,5,3,2)
df$x[df$pos == 2] <- 6
df$x[df$pos == 3] <- 5
df$x[df$pos == 4] <- 2
df$x[df$pos == 5] <- 4
df$x[df$pos == 6] <- 3

df$ang <- df$ang /pi*180 + 180
df$ang <- if_else(df$ang > 180, df$ang-360, df$ang)

df3 <- df %>%
  group_by(pos, x) %>%
  summarise(mean = mean(ang),
            sd = sd(ang),
            n = n(),
            se = sd / sqrt(n)
            # lower = mean(ang)-se(ang),
            # higher = mean(ang)+se(ang)
            ) %>%
  ungroup() %>%
  data.frame()

# T4b
df2 <- cbind(ang[,2], seq(1,6)-0.1, 2)
df2 <- as.data.frame(df2)
colnames(df2) <- c('ang','pos','type')
df2$ang <- df2$ang /pi*180 + 180

windows(width = 6, height = 3)
# pdf(paste("H2_vs_T4b_se", '.pdf', sep = ''), width = 2, height = 2)
ggplot( ) + 
  geom_point(data=df3, aes(x=x, y=mean), colour='red', size = 1) +
  geom_errorbar(data=df3, aes(x=x, ymin=mean-se, ymax=mean+se),colour='red', width=.2) +
  # geom_errorbar(data=df3, aes(x=x, ymin=lower, ymax=higher), width=.2) +
  geom_point(data=df2, aes(x=pos, y=ang), colour= pal_T4[2], size = 1) +
  # scale_discrete_manual(values = c('red', pal_T4[2]),aesthetics = "colour",guide= guide_legend(title="")) +
  coord_cartesian(ylim = c(-5, 35)) +
  scale_y_continuous(breaks= seq(0,40,by=10), labels = paste0(seq(0,40,by=10), "°")) +
  scale_x_continuous(breaks= seq(0,6,by=1), labels = seq(0,6,by=1)) +
  theme_minimal() +
  labs(title= paste("H2 vs T4b", sep = ''),
       x='position', y='') 
# dev.off()


# # - source of bias --> prob not this alone, only about 3 degrees. 
# # lens positions alignment at 90 azimuth, compare with Eyal head angle pictures
# nopen3d()
# points3d(lens)
# points3d(lens[ind_Up_lens,], size=10)
# planes3d(1,0,0, 0, alpha = 0.2)
# planes3d(0,0,1, 0, alpha = 0.2)
# 
# dxy = lens[c(435,409), c('x','z')] %>% diff()
# atan2(dxy[2], dxy[1]) /pi*180
# 
# xy1 <- colMeans(lens[c(432,398),c('x','z')])
# xy2 <- colMeans(lens[c(445, 411),c('x','z')])
# dxy <- xy1 - xy2


# cont. ED Fig.7G, compare to Henning et al 2022 ----------------------------

# load data
T4b_pxy_H <- read.csv("data/data_Henning/T4b_pxy.csv")
T4b_vxy_H <- read.csv("data/data_Henning/T4b_vxy.csv")
T4d_pxy_H <- read.csv("data/data_Henning/T4d_pxy.csv")
T4d_vxy_H <- read.csv("data/data_Henning/T4d_vxy.csv")
T4b_12_H <- read.csv("data/data_Henning/T4b_12.csv")

# choose between T4b subgroup 1 or 2, or T4d. 
# The next section depends on this choice too

## ## T4b
LL <- 2
vv <- ucl_rot_sm + 0.5*(RF_lens_T4_pred_sm[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)
vv_Mer <- cart2Mercator(vv)
## ## subgroup, 1 or 2
sg <- 2 
pxy0 <- T4b_pxy_H[T4b_12_H == sg,]
vxy0 <- T4b_vxy_H[T4b_12_H == sg,]

# ## ## T4d
# LL <- 4 
# vv <- ucl_rot_sm + 0.5*(RF_lens_T4_pred_sm[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)
# vv_Mer <- cart2Mercator(vv)
# pxy0 <- T4d_pxy_H
# vxy0 <- T4d_vxy_H

# remove na
nona <- !is.na(pxy0[,1])
pxy0 <- pxy0[nona,]
vxy0 <- vxy0[nona,]
pxy <- pxy0
vxy <- pxy0 + vxy0

# convert to Mercator
elaz <- pxy[, c(2,1)]
elaz[,2] <- -elaz[,2]
pxy <- cart2Mercator(sph2cartZ(elaz2sph(elaz)))
elaz <- vxy[, c(2,1)]
elaz[,2] <- -elaz[,2]
vxy <- cart2Mercator(sph2cartZ(elaz2sph(elaz)))

df_arrow_H <- data.frame(cbind(pxy, vxy ) )
colnames(df_arrow_H) <- c('x','y','xend','yend')

df_arrow <- data.frame(cbind(ucl_rot_Merc, vv_Mer) )
colnames(df_arrow) <- c('x','y','xend','yend')

plt <- plt_Mer +
  # geom_path(data = cmer_Merc, aes(x=x, y=y), colour = pal_axes[3], lwd=1) +
  geom_segment(data=df_arrow_H, aes(x = x,y = y, xend = xend,yend = yend), colour='black',size =1) +
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =1) +
  # geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],arrow = arrow(length = unit(0.008, "npc"), type = "closed"), size =1) +
  coord_fixed(ratio = 1) +
  scale_y_continuous(limits = log(tan(pi/4 + c(-60,60)/180*pi/2)),
                     breaks = log(tan(pi/4 + seq(-60,60,by=10)/180*pi/2)),
                     labels= seq(-60,60,by=10), expand = c(0.02, 0))+
  scale_x_continuous(limits = c(-20,120)/180*pi, breaks = seq(-20,120,by=10)/180*pi, 
                     labels = seq(-20,120,by=10), expand = c(0.02, 0)) +
  labs(title = "")

windows(width = 8, height = 8)
plt
# ggsave(paste("T4", letters[LL], "_group_", sg, "_Henning_Merc.pdf",sep=''), width = 8, height = 8)


# cont. 10-dge bin avg ----------------------------------------------------------
# [-10, 40], [-10 80]
divxdeg <- seq(-10,80,by=10)
divx <- divxdeg/ 180*pi
divydeg <- seq(-10,40,by=10)
divy <- log(tan(pi/4 + divydeg /180*pi/2))

# -- Henning
nr <- nrow(pxy0)
angs_ref <- matrix(rep(c(0,1), nr), ncol = 2, byrow = T)
angs_H <- angcos_vf(angs_ref, vxy0)
angs_H <- angs_H + 90

angavg_x_H <- vector(mode = "numeric", length = length(divx)-1)
for (j in 1:(length(divx)-1)) {
  # in x,y bins
  ii <- (pxy0[,2] > divydeg[1]) & (pxy0[,2] < divydeg[length(divydeg)])
  ii2 <- (pxy0[,1] > divxdeg[j]) & (pxy0[,1] < divxdeg[j+1])
  ii <- ii & ii2
  angavg_x_H[j] <- mean(angs_H[ii])
}

angavg_y_H <- vector(mode = "numeric", length = length(divy)-1)
for (j in 1:(length(divy)-1)) {
  # in x,y bins
  ii <- (pxy0[,1] > divxdeg[1]) & (pxy0[,1] < divxdeg[length(divxdeg)])
  ii2 <- (pxy0[,2] > divydeg[j]) & (pxy0[,2] < divydeg[j+1])
  ii <- ii & ii2
  angavg_y_H[j] <- mean(angs_H[ii])
}

# -- T4 reg
pp <- ucl_rot_Merc
vv <- vv_Mer - ucl_rot_Merc
nr <- nrow(vv)
angs_ref <- matrix(rep(c(0,1), nr), ncol = 2, byrow = T)
angs <- angcos_vf(angs_ref, vv)
angs <- angs + 90

angavg_x <- vector(mode = "numeric", length = length(divx)-1)
for (j in 1:(length(divx)-1)) {
  # in x,y bins
  ii <- (pp[,2] > divy[1]) & (pp[,2] < divy[length(divy)])
  ii2 <- (pp[,1] > divx[j]) & (pp[,1] < divx[j+1])
  ii <- ii & ii2
  angavg_x[j] <- mean(angs[ii])
}

angavg_y <- vector(mode = "numeric", length = length(divy)-1)
for (j in 1:(length(divy)-1)) {
  # in x,y bins
  ii <- (pp[,1] > divx[1]) & (pp[,1] < divx[length(divx)])
  ii2 <- (pp[,2] > divy[j]) & (pp[,2] < divy[j+1])
  ii <- ii & ii2
  angavg_y[j] <- mean(angs[ii])
}

# PLOT
## ## xbin
pp <- cbind(divxdeg[-1] - 5, 0)
vv <- cbind(cos(angavg_x/180*pi), sin(angavg_x/180*pi))
vv_H <- cbind(cos(angavg_x_H/180*pi), sin(angavg_x_H/180*pi))
# ## ## ybin
# pp <- cbind(0, divydeg[-1] - 5)
# vv <- cbind(cos(angavg_y/180*pi), sin(angavg_y/180*pi))
# vv_H <- cbind(cos(angavg_y_H/180*pi), sin(angavg_y_H/180*pi))

df_arrow <- data.frame(cbind(pp, pp+vv*3 ) )
colnames(df_arrow) <- c('x','y','xend','yend')
df_arrow_H <- data.frame(cbind(pp, pp+vv_H*3 ) )
colnames(df_arrow_H) <- c('x','y','xend','yend')

plt <- ggplot() +
  geom_segment(data=df_arrow_H, aes(x = x,y = y, xend = xend,yend = yend), colour='black',size =1) +
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =1) +
  coord_fixed(ratio = 1) +
  scale_y_continuous(limits = range(divydeg), breaks = divydeg, 
                     labels= divydeg, expand = c(0.02, 0))+
  scale_x_continuous(limits = range(divxdeg), breaks = divxdeg, 
                     labels = divxdeg, expand = c(0.02, 0)) +
  labs(title = "")

windows(width = 8, height = 5)
plt

# ggsave(paste0("avg_ang_xbins_T4", letters[LL], "_group_", sg, ".pdf"), width = 7, height = 5)
# ggsave(paste0("avg_ang_ybins_T4", letters[LL], "_group_", sg, ".pdf"), width = 7, height = 5)


# ED Fig.6B, ED Fig.7D, size, T4b and T4d on eye, use T4 gallery, -------------------------------------

LL <- 2 # choose type

# - eye
com_xyz <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
v0 <- as.matrix(lens_type[[LL]][, c("x0","y0","z0")])
v1 <- as.matrix(lens_type[[LL]][, c("xd","yd","zd")])

PD_eye <- matrix(ncol = 1, nrow = nrow(v0)) # PD length normalized
for (j in 1:nrow(v0)) {
  # PD_eye[j,] <- acos((sum(v0[j,]^2)+sum(v1[j,]^2)-sum((v1[j,]-v0[j,])^2))/2/sqrt(sum(v0[j,]^2))/sqrt(sum(v1[j,]^2))) /pi*180 # in [deg]
  # PD_eye[j,] <- acos((1+sum(v1[j,]^2)-sum((v1[j,]-v0[j,])^2))/2/1/sqrt(sum(v1[j,]^2))) /pi*180 # in [deg]
  PD_eye[j,] <- angcos(v1[j,], v0[j,])
} 


# -- size with interp
np_eval <- ucl_rot_sm
np_pred <- matrix(nrow = nrow(np_eval), ncol = 1)

npdata <- data.frame(mc = com_xyz,ec = PD_eye)
colnames(npdata) <- c('mc.1','mc.2','mc.3','ec')
bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'fixed', regtype= 'll')
model_np <- npreg(bw)
for (j in 1:nrow(np_eval)) {
  np_eval_one <- data.frame(mc.1 = np_eval[j,1], mc.2 = np_eval[j,2], mc.3 = np_eval[j,3])
  np_pred[j] <- predict(model_np, newdata = np_eval_one)
}


# - size from 3d
vv <- (RF_lens_T4_pred_sm[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)
size_3d <- sqrt(rowSums(vv^2))

# # - histo
# range((np_pred))
# dev.new()
# # pdf(paste("PD_ang_eye_T4", letters[LL], ".pdf", sep = ''))
# hh <- hist(np_pred, breaks = seq(6, 36, by=1), plot = F)
# plot(hh$mids, hh$counts, type='l', bty='n', xlim = c(8, 24), xaxt='n',yaxt='n', xlab ="angle [deg]", ylab='counts')
# axis(1, at = seq(8,24,by=2), labels = paste(seq(8, 24, by =2), "°", sep = '') )
# # dev.off()


# - plot
df_pos <- data.frame(ucl_rot_Mo)
df_pos$val <- np_pred
# df_pos$val <- size_3d /pi*180

range((df_pos$val))
quantile(na.omit(df_pos$val), c(0, 0.05,0.25, 0.5, 0.75,0.95, 1))

# df_pos$valgp <- cut(df_pos$val,c(8, 12, 14, 16, 24)) #T4b
# df_pos$valgp <- cut(df_pos$val, quantile(na.omit(df_pos$val), c(0,0.25, 0.5, 0.75, 1))) #T4b
# df_pos$valgp <- cut(df_pos$val,c(8, 9, 10, 12, 18)) #T4d
# df_pos$valgp <- cut(df_pos$val,c(8, 10, 12, 16, 24)) #T4b & d, 2021
df_pos$valgp <- cut(df_pos$val,c(8, 10, 12, 15, 18)) #T4b & d, 2024

# -- 2d
plt <- ggplot() +
  geom_polygon(data = as.data.frame(bkgd_str_equa), aes(x=xM, y=yM), fill = 'grey90', alpha=1) +
  geom_polygon(data = as.data.frame(bkgd_str_meri), aes(x=xM, y=yM), fill = 'grey90', alpha=1) +
  geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  # geom_path(data=as.data.frame(bkgd_str_equa), aes(xM, yM), colour='gray', alpha=1,linetype=1, lwd=1) +
  # geom_path(data=as.data.frame(bkgd_str_meri), aes(xM, yM), colour='gray', alpha=1,linetype=1, lwd=1) +
  geom_point(data=df_pos, aes(x = xM, y = yM, colour = valgp), size = 2, na.rm = T) +
  scale_color_manual(values = pal_heat2,guide= guide_legend(title="ang[deg]"), na.translate=T) +
  scale_x_continuous(limits = c(- 2*sqrt(2)/2, 2*sqrt(2)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-sqrt(2), sqrt(2)), expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = c(.9, .9), panel.background = element_blank()) +
  coord_fixed(ratio = 1) +
  labs(title = paste("T4", letters[LL],"_PD [deg]", "_6nb_eye", sep = ""))
windows(width = 12, height = 8)
# pdf(paste("T4", letters[LL], "_PD_eye", '.pdf', sep = ""), width = 8.5, height = 4.5)
plt
# dev.off()


# - hist along equ
str_hw <- 15

# xlab ="azim[deg]", pch=15, ylab='PD ang [deg]'

windows(width = 8, height = 4)
# pdf(paste0("PD_T4", letters[LL], "_equ.pdf"), width = 8, height = 4)
thr <- cos((90- str_hw)/180*pi)
ii <- com_xyz[,3] < thr & com_xyz[,3] > -thr
com_rtp <- cart2sphZ(com_xyz)
x <- 360 - com_rtp[ii,3]/pi*180
x <- if_else(x>180, x-360, x)
plot(x, PD_eye[ii], xlab ="", pch=15, ylab='', main = paste("T4", letters[LL], " PD along equ", sep = ""),
     ylim = c(5, 25), xlim= c(-10, 150), xaxt='n',yaxt='n')
# plot(x, PD_eye[ii], xlab ="azim[deg]", pch=15, ylab='PD ang [deg]', main = paste("T4", letters[LL], " PD along equ", sep = ""),
#      ylim = c(5, 20), xlim= c(-30, 150), xaxt='n',yaxt='n')
axis(1, at = seq(-30, 150, by = 30), labels = seq(-30, 150, by =30) )
axis(2, at = seq(5, 25, by =5), labels = seq(5, 25, by =5) )

#interp
thr <- cos((90- str_hw)/180*pi)
ii <- ucl_rot_sm[,3] < thr & ucl_rot_sm[,3] > -thr
com_rtp <- cart2sphZ(ucl_rot_sm)
x <- 360 - com_rtp[ii,3]/pi*180 
x <- if_else(x>180, x-360, x)
points(x, np_pred[ii], col=pal_heat2[2], type='p', pch=16)
# dev.off()

# - hist along meri
windows(width = 8, height = 4)
# pdf(paste0("PD_T4", letters[LL], "_mer.pdf"), width = 8, height = 4)
xyz <- com_xyz
xyz[,2] <- -xyz[,2]
com_rtp <- cart2sphZ(xyz)
ii <- com_rtp[,3] > (45-str_hw)/180*pi & com_rtp[,3] < (45+str_hw)/180*pi
x <- com_rtp[ii,2]/pi*180
# plot(x, PD_eye[ii], pch=15, xlab ="", ylab='', main = paste("T4", letters[LL], " PD along central meridian", sep = ""),
#      ylim = c(10,25), xlim= c(10, 170), xaxt='n',yaxt='n')
plot(x, PD_eye[ii], pch=15, xlab ="", ylab='', main = paste("T4", letters[LL], " PD along central meridian", sep = ""),
     ylim = c(5,20), xlim= c(170, 10), xaxt='n',yaxt='n')
axis(1, at = seq(180, 0, by = -45), labels = seq(-90, 90, by = 45) )
axis(2, at = seq(5, 25, by =5), labels = seq(5, 25, by =5) )

#interp
xyz <- ucl_rot_sm
xyz[,2] <- -xyz[,2]
com_rtp <- cart2sphZ(xyz)
ii <- com_rtp[,3] > (45-str_hw)/180*pi & com_rtp[,3] < (45+str_hw)/180*pi
x <- com_rtp[ii,2]/pi*180 
points(x, np_pred[ii], col=pal_heat2[2], type='p', pch=16)
# dev.off()

# # - real T4
# com_xyz_Mo <- com_xyz
# colnames(com_xyz_Mo) <- c('x','y','z')
# com_xyz_Mo %<>% as_tibble() %>%  
#   mutate(y = -y) %>%
#   mutate(theta = acos(z)) %>%
#   mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
#   mutate(t = theta / pi * 180, p = phi/pi*180) %>%  
#   mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
#   as.data.frame()
# com_xyz_Mo <- Mollweide(com_xyz_Mo[,c('t', 'p')])
# colnames(com_xyz_Mo) <- c('xM','yM')
# 
# df_pos <- as.data.frame(com_xyz_Mo)
# df_pos$val <- PD_eye
# 
# range(na.omit(df_pos$val))
# quantile(na.omit(df_pos$val), c(0.01, 0.99))
# 
# # df_pos$valgp <- cut(df_pos$val,c(9, 14, 17, 20, 25)) # T4b
# df_pos$valgp <- cut(df_pos$val,c(6, 10, 14, 18, 22)) # T4d
# windows(width = 12, height = 8)
# ggplot() +
#   geom_polygon(data = as.data.frame(bkgd_str_equa), aes(x=xM, y=yM), fill = 'grey90', alpha=1) +
#   geom_polygon(data = as.data.frame(bkgd_str_meri), aes(x=xM, y=yM), fill = 'grey90', alpha=1) +
#   geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
#   geom_point(data=df_pos, aes(x = xM, y = yM, colour = valgp), size = 3) +
#   scale_color_manual(values = pal_heat2,guide= guide_legend(title="ang[deg]")) +
#   geom_text(data=df_pos, aes(x=xM, y=yM, label= round(val,1)), nudge_x = 0.1 )+
#   ylab("elevation") + xlab("azimuth") +
#   theme_minimal() +
#   scale_x_continuous(limits = c(-2*sqrt(2)/2, 2*sqrt(2)), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(-sqrt(2), sqrt(2)), breaks = c(-1.5,0,1.5), labels = c(-1.5,0,1.5), expand = c(0, 0)) + # set +y as above eq
#   theme(panel.grid = element_blank(), axis.text = element_blank()) +
#   labs(title = paste("T4", letters[LL],"_real_PD [deg]", "_6nb_eye", sep = "")) +
#   coord_fixed(ratio=1)

# Fig.4D, angles between eye v-axis and T4b, 2D heat map ------------------------------------------------

# np in medulla
LL <- 2
v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix()
v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
RF_med <- cbind(v0, v1)
colnames(RF_med) <- c()
np_eval <- med_xyz_aux
Npt <- nrow(eyemap)
vf_pred <- matrix(nrow = dim(np_eval)[1], ncol = 3)
for (k in 1:3) {
  npdata <- data.frame(mc = RF_med[,1:3], ec = RF_med[,3+k])
  bw <- npregbw(formula= ec~mc.1+mc.2+mc.3, data= npdata, bwtype= 'fixed', regtype= 'll')
  model_np <- npreg(bw)
  for (j in 1:nrow(np_eval)) {
    np_eval_one <- data.frame(mc.1 = np_eval[j,1], mc.2 = np_eval[j,2], mc.3 = np_eval[j,3])
    vf_pred[j,k] <- predict(model_np, newdata = np_eval_one)
  }
}
RF_med_T4b_pred <- vf_pred[1:Npt,]

# # DEBUG
# nopen3d()
# points3d(med_xyz)
# for (j in 1:Npt) {
#   arrow3d(med_xyz[j,], RF_med_T4b_pred[j,], col='blue' )
# }
# for (j in 1:nrow(RF_med_T4b)) {
#   arrow3d(RF_med_T4b[j,1:3], RF_med_T4b[j,4:6], col='red' )
# }

v_data <- (RF_lens_T4_pred_sm[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm) 
v_data_med <- -(RF_med_T4b_pred - med_xyz) # note chiasm

# nopen3d()
# for (j in 1:nrow(v_data)) {
#   arrow3d(ucl_rot_sm[j,1:3], ucl_rot_sm[j,1:3]+v_data[j,], col='red' )
# }

# -- vertical axis
ang_T4_vaxis_eye <- matrix(NA, ncol = 1, nrow = nrow(eyemap))
for (j in 1:nrow(eyemap)) {
  # for (j in 1:12) {
  # 2023
  if (sum(is.na(nb_ind[j,c(3,6)])) == 0) {
    xyz2 <- ucl_rot_sm[nb_ind[j,1],,drop=F]
    vv12 <- ucl_rot_sm[nb_ind[j,3],] - ucl_rot_sm[nb_ind[j,1],] #vector goes bottom --> top, along +h axis
    vv12 <- vv12 - sum(vv12 * xyz2) * xyz2 # tangent vector from one lens to another
    vv23 <- ucl_rot_sm[nb_ind[j,1],] - ucl_rot_sm[nb_ind[j,6],]
    vv23 <- vv23 - sum(vv23 * xyz2) * xyz2 # tangent vector from one lens to another
    vvT4 <- v_data[nb_ind[j,1], ,drop=F ]
    ang12 <- acos(vv12 %*% t(vvT4) / sqrt(sum(vv12^2)) / sqrt(sum(vvT4^2)) )
    ang12 <- ang12 * sign( xyz2 %*% cross3D(vvT4, vv12) )
    ang23 <- acos(vv23 %*% t(vvT4) / sqrt(sum(vv23^2)) / sqrt(sum(vvT4^2)) )
    ang23 <- ang23 * sign( xyz2 %*% cross3D(vvT4, vv23) )
    ang_T4_vaxis_eye[j] <- (ang12+ang23)/2/pi*180
  }
} 
ang_T4b_vaxis_eye <- ang_T4_vaxis_eye


range(ang_T4_vaxis_eye, na.rm = T)
quantile(na.omit(ang_T4b_vaxis_eye), c(0.05,0.1, 0.25, 0.5, 0.75, 0.9, 0.95))


# - hist

# -- eye, combine with hex hist -- cf Figure_3_lattice.R, 
dev.new()
# pdf("hist T4b and v-axis on eye.pdf")
# hh <- hist(ang_T4_vaxis_eye, breaks = seq(25, 165, by=5), plot = F) #2021
# plot(hh$mids, hh$counts, type='l', bty='n', col=pal_T4[2], xlim = c(30, 150), xaxt='n',yaxt='n', xlab ="", ylab="", main = "T4b and v-axis on eye")
hh <- hist(ang_T4_vaxis_eye, breaks = seq(30, 150, by=5), plot = F) #2023
plot(hh$mids, hh$counts, type='l', bty='n', col=pal_T4[2], xlim = c(60, 120), xaxt='n',yaxt='n', xlab ="", ylab="", main = "T4b and v-axis on eye")

nb_ang_ucl <- matrix(ncol = 2, nrow = nrow(eyemap))
for (j in 1:nrow(nb_ind)) {
  if (sum(complete.cases(nb_ind[j,])) == 7) {
    # eye
    # 2023
    bt <- ucl_rot_sm[nb_ind[j,3], ] - ucl_rot_sm[nb_ind[j,6], ]
    bf <- colMeans(ucl_rot_sm[nb_ind[j,c(4,5)],]) - colMeans(ucl_rot_sm[nb_ind[j,c(2,7)],])
    ang <- acos(bt %*% bf / sqrt(sum(bt^2)) / sqrt(sum(bf^2)) ) /pi*180 
    nb_ang_ucl[j,] <- c(j, ang)
  }
}
hh <- hist(nb_ang_ucl[,2], breaks = seq(30, 150, by=5), plot = F)
points(hh$mids, hh$counts, type='l', bty='n', xlim = c(30, 150), xaxt='n',yaxt='n', xlab ="angle [deg]", ylab='counts')

axis(1, at = seq(30, 150, by =30), labels =  paste(seq(30, 150, by =30), "°", sep = '') )
# dev.off()


# - density
windows(width = 5, height = 4)
# pdf(paste("den_T4b_vaxis.pdf", sep = ""), width = 5, height = 4)
dd <- density(ang_T4_vaxis_eye, from= 30, to= 150, bw='SJ', na.rm = T)
plot(dd$x, dd$y, type='l', bty='n', col='black', lwd=3, xlim = c(30, 150), xaxt='n',yaxt='n',
     xlab ="", ylab='',main = '')
axis(1, at = seq(30, 150, by =30), labels = paste(seq(30, 150, by =30), "°", sep = ''), cex.axis = 1.5 )
axis(2, at = c(0, 1/512/diff(dd$x[1:2]) ), labels = c('0', '1'), cex.axis = 1.5 )
# dev.off()

# # -- h-axis
# dev.new()
# hh <- hist(ang_T4_haxis_eye, breaks = seq(-40, 40, by=5), plot = F)
# plot(hh$mids, hh$counts, type='l', bty='n', col=pal_T4[2], xlim = c(-40, 40), yaxt='n', xlab ="", ylab="", main = "T4b and h-axis on eye")


# # -- med
# dev.new()
# # pdf("hist T4b and v-axis in med.pdf")
# hh <- hist(ang_T4_vaxis_med, breaks = seq(30, 155, by=5), plot = F)
# plot(hh$mids, hh$counts, type='l', bty='n', xlim = c(30, 150), xaxt='n',yaxt='n', xlab ="", ylab="", main = "T4b and v-axis in med")
# axis(1, at = seq(30, 150, by =60), labels = paste(seq(30, 150, by =60), "°", sep = '') )
# # dev.off()


# - PLOT 2d Mollweide, angle
df_pos <- data.frame(ucl_rot_Mo)
df_pos$quan <- ang_T4b_vaxis_eye
quantile(ang_T4b_vaxis_eye, c(0.01,0.05,0.25,0.5,0.75,0.95,0.99), na.rm = T)
rg <- c(60, 90, 120) 

plt <- plt_Momin + 
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  geom_point(data=df_pos, aes(x = xM, y = yM, colour = quan), size = 2) +
  scale_color_gradientn(colours = pal_heat1, values = scales::rescale(rg), limits=range(rg), oob=scales::squish,
                        breaks= rg, labels= rg, guide = guide_colorbar(title = "ang") ) +
  # geom_point(data=df_na, aes(x=xM, y=yM), colour='gray', shape=1, size=2, stroke=1) +
  theme(legend.position = c(.9, .9) ) +
  # geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],arrow = arrow(length = unit(0.008, "npc"), type = "closed"), size =1) +
  labs(title = "ang betw T4b and v-axis on eye")
for (k in 1:3) {
  j <- ind_eghex_2[k]
  pt <- ucl_rot_Mo[nb_ind[j,], ]
  pt <- pt[c(2,3,4,5,6,7,2), ] 
  plt <- plt + geom_path(data=as.data.frame(pt), aes(x=xM, y=yM), colour='gray50',lwd=0.5)
}

windows(width = 9, height = 6)
# pdf("ang_T4b_v_eye.pdf", width = 8.5, height = 4.5)
plt
# dev.off()


#  example T4 angle, local meridian aligned ----------------------------------------------------------
nb_coord <- rbind(c(-5/180*pi,0),
                  c(+5/180*pi,0),
                  c(0, -5/180*pi),
                  c(0, +5/180*pi) )

# PD
LL <- 2
T4com <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
T4dir <- dir_type[[LL]][, c('rsxd','rsyd','rszd', 'rsx0','rsy0','rsz0')] %>% as.matrix() # reverse st optic flow dir
T4comeye <- lens_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
T4direye <- lens_type[[LL]][, c('xd','yd','zd', 'x0','y0','z0')] %>% as.matrix() # reverse st optic flow dir

# PD field
PD_eye <- 0.5*(RF_lens_T4_pred[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)
vv_eye <- ucl_rot_sm + PD_eye
# vv <- ucl_rot_sm + 0.5*(RF_lens_T4_pred_sm[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)
vv_Mo <- sweep(vv_eye, 1, sqrt(rowSums(vv_eye^2)), '/')
vv_Mo <- cart2sph2tp(vv_Mo)
vv_Mo <- Mollweide(vv_Mo[,c('t', 'p')])
colnames(vv_Mo) <- c('xM','yM')
rownames(vv_Mo) <- rownames(ucl_rot_sm)

# cardinal translation
vt_car <- t_gen(c(-1,0,0), as.matrix(ucl_rot_sm))*0.1 + ucl_rot_sm
vt_car <- sweep(vt_car, 1, sqrt(rowSums(vt_car^2)), '/')

vR_car <- R_gen(c(0,0,-1), as.matrix(ucl_rot_sm))*0.1 + ucl_rot_sm
vt_car <- sweep(vR_car, 1, sqrt(rowSums(vR_car^2)), '/')

hex_eg <- list()
vv_eg <- list()
vt_car_eg <- list()
vR_car_eg <- list()
for (k in 1:3) {
  j <- ind_eghex_2[k]
  # ind_T4 <- rowSums(sweep(T4com, 2, med_xyz[j,])^2) %>% which.min()
  # ind_T4 <- c(108,139,170)[k] #picked in figure 2
  # vv <- matrix(T4direye[ind_T4,], nrow = 2, byrow=T) # real T4
  vv <- rbind(ucl_rot_sm[j,], PD_eye[j,]) #regression T4
  vtcar <- vt_car[j,,drop=F]
  vRcar <- vR_car[j,,drop=F]
  
  # eye
  pt <- ucl_rot_sm[nb_ind[j,], ] 
  tp <- cart2sphZ(pt[1,,drop=F])[,2:3]
  rtp <- cbind(1, sweep(nb_coord,2,tp,'+'))
  xyz <- sph2cartZ(rtp)
  
  ## ## med
  # ind_nb <- rowSums(sweep(med_xyz, 2, as.numeric(T4com[ind_T4, ]))^2) %>% order() %>% head(1+6+12)
  # pt <- med_xyz[nb_ind[j,], ] # med
  # # np
  # xyz <- eyemap_np(ucl_rot_sm[ind_nb,], med_xyz[ind_nb,], xyz)
  # vv <- matrix(T4dir[ind_T4,], nrow = 2, byrow=T)
  
  zz <- cross3D((xyz[1,]-pt[1,]), (xyz[4,]-pt[1,])) #pointing inwards
  pc <- prcomp(xyz)
  if (pc$rotation[,3] %*% zz < 0) {
    pc$rotation <- -pc$rotation
  }
  if (c(cross3D(pc$rotation[,3], pc$rotation[,1])) %*% pc$rotation[,2] < 0) {
    pc$rotation[,2] <- -pc$rotation[,2]
  }
  
  xyz <- sweep(xyz, 2, pc$center) %*% pc$rotation
  pt <- sweep(pt, 2, pc$center) %*% pc$rotation
  
  vv <- sweep(vv, 2, pc$center) %*% pc$rotation
  vtcar <- sweep(vtcar, 2, pc$center) %*% pc$rotation
  vRcar <- sweep(vRcar, 2, pc$center) %*% pc$rotation
  
  xyz[,3] <- 0
  pt[,3] <- 0
  vv[,3] <- 0
  vtcar[,3] <- 0
  vRcar[,3] <- 0
  
  # align to v-axis
  va <- (xyz[1,] - xyz[2,]) #vertical axis of hex
  ang <- acos(va %*% c(0,1,0) / sqrt(sum(va^2)))
  if (c(cross3D(va, c(0,1,0))) %*% c(0,0,1) < 0) {
    ang <- -ang
  }
  pt <- pt %*% matrix(c(cos(ang), sin(ang), 0,
                        -sin(ang), cos(ang), 0,
                        0,0,1), ncol=3, byrow=T)
  vv <- vv %*% matrix(c(cos(ang), sin(ang), 0,
                        -sin(ang), cos(ang), 0,
                        0,0,1), ncol=3, byrow=T)
  vtcar <- vtcar %*% matrix(c(cos(ang), sin(ang), 0,
                        -sin(ang), cos(ang), 0,
                        0,0,1), ncol=3, byrow=T)
  vRcar <- vRcar %*% matrix(c(cos(ang), sin(ang), 0,
                              -sin(ang), cos(ang), 0,
                              0,0,1), ncol=3, byrow=T)
  
  hex_eg[[k]] <- pt
  vv_eg[[k]] <- vv
  vt_car_eg[[k]] <- vtcar
  vR_car_eg[[k]] <- vRcar
}

# PLOT
j <- 3 # j=1,2,or3
pt <- hex_eg[[j]]
vv <- vv_eg[[j]]/2
windows(width = 8, height = 8)
# pdf(file = paste("hex_eg_T4_reg_", c("D","C","V")[j], ".pdf", sep='')) # from regression
plot(pt, pch=16, axes=F, ann=F, cex=2, asp=1)
# plot(pt, xlim=c(-0.1,0.4), ylim=c(-0.1,0.5))

lines(pt[c(2,3,4,5,6,7,2),1:2], lwd=3)
lines(pt[c(3,6), 1:2], lty = "dotted", lwd=2)
a1 <- colMeans(pt[c(2,7),]) #right midpoint
a2 <- colMeans(pt[c(4,5),]) # left
va <- pt[3,] - pt[6,] # vertical vec upwards

vT4 <- c(diff(vv))
# arrows(a1[1], a1[2], a2[1], a2[2], lwd=2)
arrows(vv[1,1], vv[1,2], vv[2,1], vv[2,2], lwd=2, col=pal_T4[LL])
# arrows(pt[1,1], pt[1,2], vt_car_eg[[j]][1], vt_car_eg[[j]][2], lwd=2, col=pal_TR[1])
# arrows(pt[1,1], pt[1,2], vR_car_eg[[j]][1], vR_car_eg[[j]][2], lwd=2, col=pal_TR[6])
title(paste(c("D","C","V")[j],
            " alpha= ",
            round(acos(vT4 %*%  va / sqrt(sum(va^2)) / sqrt(sum(vT4^2))) /pi*180, 1),
            sep = '')
      )
# dev.off()

# comparison: from RF_pred_RData ------------------------------------------------------

# generate smapling axes
maxis_tpxyz <- axis_sampling(da = 2)

# ED Fig.6E, t,R prediction, T4b  -------------------------------------------------

# position
d <- ucl_rot_sm
colnames(d) <- c("X","Y","Z")

# - optimal t, R
LL <- 2
v_data <- RF_lens_T4_pred[,(3*(LL-1)+1):(3*(LL-1)+3)]- ucl_rot_sm 
# v_data <- -v_data  #T4a
v_data <- v_data - sweep(ucl_rot_sm, 1, rowSums(v_data * ucl_rot_sm), '*') #tangent
ii <- ucl_rot_sm[,2] < -0.98
p <- v_data / mean(sqrt(rowSums(v_data[ii,]^2)))
vf_amp <- mean(sqrt(rowSums(v_data[ii,]^2)))


# start_time <- Sys.time()
# main loop
R_cErr <- c()
t_cErr <- c()
# nd <- sum(!is.na(p[,1]))
ii <- ucl_rot_sm[,2] < -0.98
for (k in 1:dim(maxis_tpxyz)[1]) {
  # rot
  v <- R_gen(maxis_tpxyz[k,3:5], d)
  v <- v / mean(sqrt(rowSums(v[ii,]^2)))
  
  # Err <- sqrt(rowSums((v - p)^2)) ## L2
  Err <- angcos_vf(v, p) ## ang
  R_cErr <- c(R_cErr, mean(Err,na.rm = T))
  
  # trans
  # v <- t_gen(maxis_tpxyz[k,3:5], d) * max(vf_amp)
  v <- t_gen(maxis_tpxyz[k,3:5], d)
  v <- v / mean(sqrt(rowSums(v[ii,]^2)))
  
  # Err <- sqrt(rowSums((v - p)^2)) ## L2
  Err <- angcos_vf(v, p) ## ang
  t_cErr <- c(t_cErr, mean(Err, na.rm = T))
}
# end_time <- Sys.time()
# end_time - start_time

# mean square error
R_errLand <- data.frame(phi = maxis_tpxyz[,2], theta = maxis_tpxyz[,1], error = R_cErr)
t_errLand <- data.frame(phi = maxis_tpxyz[,2], theta = maxis_tpxyz[,1], error = t_cErr)

# tR prediction, [theta, phi, x, y ,z]
pred_R <- maxis_tpxyz[which.min(R_cErr),]
pred_t <- maxis_tpxyz[which.min(t_cErr),]

pred_t_T4b <- pred_t[3:5]
pred_R_T4b <- pred_R[3:5]


# - h-axis
hex_hvf_eye <- matrix(ncol = 3, nrow = nrow(eyemap))
hex_hvf_med <- matrix(ncol = 3, nrow = nrow(eyemap))
for (j in 1:nrow(nb_ind)) {
  if (sum(complete.cases(nb_ind[j,])) == 7) {
    # eye
    bt <- ucl_rot_sm[nb_ind[j,3], ] - ucl_rot_sm[nb_ind[j,6], ]
    bf <- colMeans(ucl_rot_sm[nb_ind[j,4:5],]) - colMeans(ucl_rot_sm[nb_ind[j,c(2,7)],])
    hex_hvf_eye[j,] <- bf
    # med
    bt <- med_xyz[nb_ind[j,3], ] - med_xyz[nb_ind[j,6], ]
    bf <- colMeans(med_xyz[nb_ind[j,4:5],]) - colMeans(med_xyz[nb_ind[j,c(2,7)],])
    hex_hvf_med[j,] <- bf
  }
}
v_data <- hex_hvf_eye
p <- v_data / max(na.omit(sqrt(rowSums(v_data^2))))

# main loop
R_cErr <- c()
t_cErr <- c()
ii <- ucl_rot_sm[,2] < -0.98
for (k in 1:dim(maxis_tpxyz)[1]) {
  # rot
  v <- R_gen(maxis_tpxyz[k,3:5], d)
  v <- v / mean(sqrt(rowSums(v[ii,]^2)))
  
  # Err <- sqrt(rowSums((v - p)^2)) ## L2
  Err <- angcos_vf(v, p) ## ang
  R_cErr <- c(R_cErr, mean(Err,na.rm = T))
  
  # trans
  # v <- t_gen(maxis_tpxyz[k,3:5], d) * max(vf_amp)
  v <- t_gen(maxis_tpxyz[k,3:5], d)
  v <- v / mean(sqrt(rowSums(v[ii,]^2)))
  
  # Err <- sqrt(rowSums((v - p)^2)) ## L2
  Err <- angcos_vf(v, p) ## ang
  t_cErr <- c(t_cErr, mean(Err, na.rm = T))
}
# end_time <- Sys.time()
# end_time - start_time

# mean square error
R_errLand <- data.frame(phi = maxis_tpxyz[,2], theta = maxis_tpxyz[,1], error = R_cErr)
t_errLand <- data.frame(phi = maxis_tpxyz[,2], theta = maxis_tpxyz[,1], error = t_cErr)

# tR prediction, [theta, phi, x, y ,z]
pred_R <- maxis_tpxyz[which.min(R_cErr),]
pred_t <- maxis_tpxyz[which.min(t_cErr),]

pred_t_hex <- pred_t[3:5]
pred_R_hex <- pred_R[3:5]


# PLOT
windows(record = F, width = 16, height = 16)

plt_df <- data.frame(azim = t_errLand$phi, elev = t_errLand$theta, 
                     err = t_errLand$error) #need to move azim to [-180,180], elev can be chnaged via plotting axis
plt_df$azim <- if_else(plt_df$azim > 180, plt_df$azim - 360, plt_df$azim)
star <- c(plt_df$azim[which.min(plt_df$err)], plt_df$elev[which.min(plt_df$err)])
gg1 <- ggplot(plt_df, aes(x=azim, y=elev, fill = err)) + 
  geom_tile(aes(width = 5, height = 5), color = "white") +
  scale_fill_gradient(low = "#fee8c8", high = "#e34a33", limit = range(t_errLand$error)) +
  annotate("text", x = star[1], y = star[2], label = "*", size = 6) + 
  # scale_x_continuous(breaks = seq(-180,180,by = 90), labels = seq(-180,180,by = 90) ) +
  scale_x_continuous(breaks = seq(-180,180,by = 90), labels = seq(-180,180,by = 90) ) +
  scale_y_reverse(breaks = seq(0,180,by = 90), labels = seq(90,-90,by = -90)) + 
  ggtitle(paste('translation ', paste0(round(pred_t_T4b,2),collapse = '/'), sep='')) +
  theme_minimal() +
  coord_fixed(ratio=1)

plt_df <- data.frame(x=seq(1,dim(maxis_tpxyz)[1]), error=sort(t_errLand[,3]))
gg2 <- ggplot(plt_df, aes(x,y=error )) + 
  geom_point() +
  theme_minimal()

plt_df <- data.frame(azim = R_errLand$phi, elev = R_errLand$theta, 
                     error = R_errLand$error)
plt_df$azim <- if_else(plt_df$azim > 180, plt_df$azim - 360, plt_df$azim)
star <- c(plt_df$azim[which.min(plt_df$err)], plt_df$elev[which.min(plt_df$err)])
gg3 <- ggplot(plt_df, aes(x=azim, y=elev, fill = error)) +
  geom_tile(aes(width = 5, height = 5), color = "white") +
  scale_fill_gradient(low = "#fee8c8", high = "#e34a33", limit = range(R_errLand[,3])) +
  annotate("text", x = star[1], y = star[2], label = "*", size = 6) +
  # scale_x_continuous(breaks = seq(-180,180,by = 90), labels = seq(-180,180,by = 90) ) +
  scale_x_continuous(breaks = seq(-180,180,by = 90), labels = seq(-180,180,by = 90) ) +
  scale_y_reverse(breaks = seq(0,180,by = 90), labels = seq(90,-90,by = -90)) +
  ggtitle(paste('rotation ', paste0(round(pred_R_T4b,2),collapse = '/'), sep='')) +
  theme_minimal() +
  coord_fixed(ratio=1)

plt_df <- data.frame(x=seq(1,dim(maxis_tpxyz)[1]), error=sort(R_errLand[,3]))
gg4 <- ggplot(plt_df, aes(x,y=error )) + 
  geom_point() +
  theme_minimal() 

ggdraw() +
  draw_plot(gg1 + theme(legend.justification = "bottom"), x = 0, y = 0.5, width = 0.7, height = 0.45) +
  draw_plot(gg2 + theme(axis.title = element_text()), 0.75, 0.6, 0.25, 0.25) +
  draw_plot(gg3 + theme(axis.title = element_text()), 0, 0, 0.7, 0.45) +
  draw_plot(gg4 + theme(axis.title = element_text()), 0.75, 0.1, 0.25, 0.25) 

# ggsave("pred_motion.pdf", width = 8.5, height = 4.5) # big!



# -- Mollweide 
star <- Mollweide(cart2sph2tp(matrix(pred_t_T4b,ncol=3))[, c('t','p')])
star_op <- Mollweide(cart2sph2tp(matrix(-pred_t_T4b,ncol=3))[, c('t','p')])
star2 <- Mollweide(cart2sph2tp(matrix(pred_R_T4b,ncol=3))[, c('t','p')])
star2_op <- Mollweide(cart2sph2tp(matrix(-pred_R_T4b,ncol=3))[, c('t','p')])

# --- trans
err_xyz <- data.frame(x = maxis_tpxyz[,3], y = maxis_tpxyz[,4], z = maxis_tpxyz[,5])
err_Mo <- Mollweide(cart2sph2tp(err_xyz)[, c('t','p')]) %>% as.data.frame()
colnames(err_Mo) <- c('xM','yM')
err_Mo$quan <- t_cErr
err_Mo <- err_Mo[err_Mo$xM != 0, ]

plt <- plt_Mo +
  # geom_point(data = err_Mo, aes(x=xM, y=yM, colour = quan), size=3 ) +
  # scale_color_gradient(low = "gray85", high= "red", limits=c(0, 180), oob=scales::squish, 
  #                      breaks= c(0,90, 180), labels= c(0,90,180), guide = guide_colorbar(title = "error [deg]"), na.value = 'yellow')
  geom_tile(data=err_Mo, aes(x=xM, y=yM, fill = quan), width = 1/sqrt(2)/10, height = 1/sqrt(8)/5, color = "white", na.rm = T) +
  scale_fill_gradient(low = "#fee8c8", high = "#e34a33", limit = range(R_errLand[,3]), guide = "colorbar", breaks=seq(30,150,by=60)) +
  # scale_x_continuous(limits = c(-sqrt(8), sqrt(8)), expand = c(0, 1)) +
  annotate("text", x = star[1], y = star[2], label = "+", size = 20) +
  annotate("text", x = star_op[1], y = star_op[2], label = "+", size = 20)+
  theme(legend.position = c(.9, .9))

windows(width = sqrt(8)*10, height = sqrt(2)*10)
plt
# ggsave("errLand_t.png", dpi = 600)
# ggsave("errLand_legend_t.pdf", width = sqrt(8)*5, height = sqrt(2)*5) # big!
# ggsave("errLand_legend_t.pdf", width = sqrt(8)*2, height = sqrt(2)*2) 

# --- rot
err_xyz <- data.frame(x = maxis_tpxyz[,3], y = maxis_tpxyz[,4], z = maxis_tpxyz[,5])
err_Mo <- Mollweide(cart2sph2tp(err_xyz)[, c('t','p')]) %>% as.data.frame()
colnames(err_Mo) <- c('xM','yM')
err_Mo$quan <- R_cErr
err_Mo <- err_Mo[err_Mo$xM != 0, ]

plt <- plt_Mo +
  geom_tile(data=err_Mo, aes(x=xM, y=yM, fill = quan), width = 1/sqrt(2)/10, height = 1/sqrt(8)/5, color = "white", na.rm = T) +
  scale_fill_gradient(low = "#fee8c8", high = "#e34a33", limit = range(R_errLand[,3])) +
  annotate("text", x = star2[1], y = star2[2], label = "x", size = 20) +
  annotate("text", x = star2_op[1], y = star2_op[2], label = "x", size = 20) 

windows(width = sqrt(8)*10, height = sqrt(2)*10)
plt
# ggsave("errLand_R.png", dpi = 600)
# ggsave("errLand_legend_R.pdf", width = sqrt(8)*2, height = sqrt(2)*2)

# --- outline
windows(width = sqrt(8)*10, height = sqrt(2)*10)
plt_Mo
# ggsave("errLand_ref.pdf", width = sqrt(8)*2, height = sqrt(2)*2)

# cont. Fig.4E, Ideal optic flow field, Mollweide ----------------------------------------

LL = 2
vv <- ucl_rot_sm + 0.5*(RF_lens_T4_pred[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)
vv_Mo <- sweep(vv, 1, sqrt(rowSums(vv^2)), '/')
vv_Mo <- cart2sph2tp(vv_Mo)
vv_Mo <- Mollweide(vv_Mo[,c('t', 'p')])
colnames(vv_Mo) <- c('xM','yM')
rownames(vv_Mo) <- rownames(ucl_rot_sm)

# cardinal directions
## ##
vt_car <- t_gen(c(-1, 0, 0), as.matrix(ucl_rot_sm))*0.2 + ucl_rot_sm
vt_car_Mo <- sweep(vt_car, 1, sqrt(rowSums(vt_car^2)), '/')
vt_car_Mo <- cart2sph2tp(vt_car_Mo)
vt_car_Mo <- Mollweide(vt_car_Mo[,c('t', 'p')])
colnames(vt_car_Mo) <- c('xM','yM')
rownames(vt_car_Mo) <- rownames(ucl_rot_sm)

vR_car <- R_gen(c(0,0,-1), as.matrix(ucl_rot_sm))*0.2 + ucl_rot_sm
vR_car_Mo <- sweep(vR_car, 1, sqrt(rowSums(vR_car^2)), '/')
vR_car_Mo <- cart2sph2tp(vR_car_Mo)
vR_car_Mo <- Mollweide(vR_car_Mo[,c('t', 'p')])
colnames(vR_car_Mo) <- c('xM','yM')
rownames(vR_car_Mo) <- rownames(ucl_rot_sm)

# T4 data
vt_opt <- t_gen(pred_t_T4b, as.matrix(ucl_rot_sm))*0.1 + ucl_rot_sm
vt_opt_Mo <- sweep(vt_opt, 1, sqrt(rowSums(vt_opt^2)), '/')
vt_opt_Mo <- cart2sph2tp(vt_opt_Mo)
vt_opt_Mo <- Mollweide(vt_opt_Mo[,c('t', 'p')])
colnames(vt_opt_Mo) <- c('xM','yM')
rownames(vt_opt_Mo) <- rownames(ucl_rot_sm)
star <- Mollweide(cart2sph2tp(matrix(-pred_t_T4b,ncol=3))[, c('t','p')])

# PLOT
df_pos <- data.frame(ucl_rot_Mo)
df_arrow <- data.frame(cbind(ucl_rot_Mo[,c('xM','yM')], vv_Mo[,c('xM','yM')]) )
colnames(df_arrow) <- c('x','y','xend','yend')

df_arrow2 <- data.frame(cbind(ucl_rot_Mo[,c('xM','yM')], vt_car_Mo[,c('xM','yM')]) )
colnames(df_arrow2) <- c('x','y','xend','yend')

df_arrow3 <- data.frame(cbind(ucl_rot_Mo[,c('xM','yM')], vR_car_Mo[,c('xM','yM')]) )
colnames(df_arrow3) <- c('x','y','xend','yend')

df_arrow4 <- data.frame(cbind(ucl_rot_Mo[,c('xM','yM')], vt_opt_Mo[,c('xM','yM')]) )
colnames(df_arrow4) <- c('x','y','xend','yend')


# subsampling
xy_sub <- seq(-20,20,by=3)
indlens_sub <- lens_ixy[lens_ixy[,2] %in% xy_sub & lens_ixy[,3] %in% xy_sub, 1]
ind_sub <- match(indlens_sub, eyemap[,2])


plt <- plt_Momin +
  # geom_path(data = chull_fov_Mo, aes(x=xM,y=yM), lwd = 0.1, colour = "black", alpha =0.5) +
  # plt_Momin +
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  # geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],arrow = arrow(length = unit(0.008, "npc"), type = "closed"), size =1) +
  # geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =1) +
  # geom_segment(data = df_arrow2[ind_sub,], aes(x = x,y = y, xend = xend,yend = yend), colour=pal_TR[1], size =1) +
  # geom_segment(data = df_arrow2[ind_sub,], aes(x = x,y = y, xend = xend,yend = yend), colour=pal_TR[2], size =1) +
  geom_segment(data = df_arrow3[ind_sub,], aes(x = x,y = y, xend = xend,yend = yend), colour=pal_TR[6], size =1)
  # geom_segment(data = df_arrow4, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_TR[1], size =1) +
  annotate("text", x = star[1], y = star[2], label = "", size = 9)

windows(width = 12, height = 8)
# pdf(paste("T4", letters[LL], "_RF_Mollweide.pdf", sep = ''), width = 9.5, height = 6.5)
# pdf(paste("vt_car_mx_Moll.pdf", sep = ''), width = 9.5, height = 6.5)
# pdf(paste("vt_car_my_Moll.pdf", sep = ''), width = 9.5, height = 6.5)
# pdf(paste("vR_car_mz_Moll.pdf", sep = ''), width = 9.5, height = 6.5)
# pdf(paste("vt_opt_Mollweide.pdf", sep = ''), width = 9.5, height = 6.5)
plt
# dev.off()
# ggsave("T4b_RF_Mollweide.pdf", width = 8.5, height = 4.5)


# Ideal flow fields, subsampled as T4 ------------------------------

# subsampling as T4 field
# xy_sub <- seq(-18,18,by=3) # LPT
xy_sub <- seq(-20,20,by=3) # as above
# indlens_sub <- lens_ixy[lens_ixy[,2] %in% seq(-16,17,by=3) & lens_ixy[,3] %in% seq(-14,16,by=3), 1]
indlens_sub <- lens_ixy[lens_ixy[,2] %in% xy_sub & lens_ixy[,3] %in% xy_sub, 1]
ind_sub <- match(indlens_sub, eyemap[,2])


pos3d <- ucl_rot_sm # index order should be the same as utp_ends_rot and RF_lens_T4_pred
colnames(pos3d) <- c("X","Y","Z")

# vf <- vf_medi_norm_ls[[j]]
# 
# # threshold small values, NOT in use
# vf_amp <- sqrt(rowSums(vf^2))
# vf_thre <- vf_amp > 0 #(max(vf_amp) * 0.01)
# 
# ii <- which(vf_thre)
# ii <- ii[ii %in% ind_sub]

ii <- na.omit(ind_sub)
# position
d <- pos3d[ii,]
nd <- dim(d)[1]

# # vector field
# p <- vf[ii,]
# p <- p - sweep(d, MARGIN = 1, STATS = rowSums(p*d), FUN = '*') # perpendicular comp
# p <- p / max(sqrt(rowSums(p^2))) # norm
# 
# # - flow field
# vv <- d + p
# vv_Mo <- sweep(vv, 1, sqrt(rowSums(vv^2)), '/')
# vv_Mo <- cart2sph2tp(vv_Mo)
# vv_Mo <- Mollweide(vv_Mo[,c('t', 'p')])
# colnames(vv_Mo) <- c('xM','yM')

df_pos <- data.frame(ucl_rot_Mo)[ii, ]

# df_arrow <- data.frame(cbind(df_pos - 0*(vv_Mo-df_pos), df_pos + 0.1*(vv_Mo-df_pos) ))
# colnames(df_arrow) <- c('x','y','xend','yend')

## ## rot
samples_xyz <- rbind(pred_R_T4b,-pred_R_T4b)
## ## trans
samples_xyz <- rbind(pred_t_T4b,-pred_t_T4b)

samples_Mo <- Mollweide(cart2sph2tp(samples_xyz)[, c('t','p')]) %>% as.data.frame()
colnames(samples_Mo) <- c('xM','yM')

# -- chull in Mollweid
ii <- chull(ucl_rot_Mo)
ii <- c(ii, ii[1])
ind_chull <- ii
chull_fov_Mo <- ucl_rot_Mo[ind_chull,] %>% as.data.frame()


plt_ls <- list()
# mang <- c()
for (k in 1:dim(samples_xyz)[1]) {
  # ## ## rot
  # v <- R_gen(samples_xyz[k, ], d) #already normalized
  ## ## trans
  v <- t_gen(samples_xyz[k, ], d) #already normalized
  
  vv <- d + v /100
  vv_Mo <- sweep(vv, 1, max(sqrt(rowSums(vv^2))), '/')
  # vv_Mo <- vv
  vv_Mo <- cart2sph2tp(vv_Mo)
  vv_Mo <- Mollweide(vv_Mo[,c('t', 'p')])
  colnames(vv_Mo) <- c('xM','yM')
  
  df_arrow_2 <- data.frame(cbind(df_pos - 0*(vv_Mo-df_pos), df_pos + 10*(vv_Mo-df_pos) ))
  colnames(df_arrow_2) <- c('x','y','xend','yend')
  
  plt_ls[[k]] <-  plt_Mo +
    # plt_Momin +
    geom_path(data = chull_fov_Mo, aes(x=xM,y=yM), lwd = 0.1, colour = "black", alpha =0.5) +
    geom_segment(data = df_arrow_2, aes(x = x,y = y, xend = xend,yend = yend), colour='black', size =0.1,alpha=1)+
    # arrow= arrow(angle=15, length = unit(0.01, "npc"), type = "closed", ends = "last")) +
    geom_point(data= samples_Mo[k,,drop=F], aes(x=xM, y=yM)) +
    labs(title = "")
}

# # save
# for (k in 1:dim(samples_xyz)[1]) {
#   # pdf(file = paste0("of_optimal_trans", k, ".pdf"), width = 5, height = 2.5 )
#   pdf(file = paste0("of_optimal_rot", k, ".pdf"), width = 9.5, height = 6.5 )
#   print(plt_ls[[k]])
#   dev.off()
# }

# cont. Fig.4F, angular error wrt cardinal/optimized motion VS lattice --------------

# - lattice h-axis
ii_6nb <- rep(0,nrow(nb_ind))
v_hex <- matrix(ncol = 6, nrow = nrow(nb_ind))

for (j in 1:nrow(nb_ind)) {
  if (sum(complete.cases(nb_ind[j,])) == 7) {
    ii_6nb[j] <- 1
    # eye
    # bt <- ucl_rot_sm[nb_ind[j,6], ] - ucl_rot_sm[nb_ind[j,7], ]
    bf <- colMeans(ucl_rot_sm[nb_ind[j,4:5],]) - colMeans(ucl_rot_sm[nb_ind[j,c(2,7)],])
    v_hex[j,] <- c(ucl_rot_sm[nb_ind[j,1], ], bf)
  }
}

# - T4b
LL <- 2
v_data <- RF_lens_T4_pred[,(3*(LL-1)+1):(3*(LL-1)+3)]- ucl_rot_sm 
v_data <- v_data - sweep(ucl_rot_sm, 1, rowSums(v_data * ucl_rot_sm), '*') #tangent
ii <- ucl_rot_sm[,2] < -0.98
v_data <- v_data / mean(sqrt(rowSums(v_data[ii,]^2)))

# - artificial motoin
# use lens dir data
pt0 <- as.matrix(ucl_rot_sm) #pts at t=0 on screen

# of
vt_car <- t_gen(c(-1,0,0), pt0)*0.1
vR_car <- R_gen(c(0,0,-1), pt0)*0.1
vt_opt <- t_gen(pred_t_T4b, pt0)*0.1
vR_opt <- R_gen(pred_R_T4b, pt0)*0.1
# vt_hex <- t_gen(pred_t_hex, pt0)*0.1
# vR_hex <- R_gen(pred_R_hex, pt0)*0.1


vt_car_y <- t_gen(c(0,-1,0), pt0)*0.1


# - diff in ang
ang_motion <- matrix(ncol = 6, nrow = nrow(ucl_rot_sm))
for (j in 1:nrow(ucl_rot_sm)) {
  if (ii_6nb[j] == 1) {
    ang_motion[j,] <- c(angcos(v_data[j,], v_hex[j,4:6]),
                        angcos(v_data[j,], vR_car[j,]),
                        angcos(v_data[j,], vt_car[j,]),
                        angcos(v_data[j,], vt_car_y[j,]),
                        angcos(v_data[j,], vR_opt[j,]),
                        angcos(v_data[j,], vt_opt[j,])   )
  }
}

df <- ang_motion 
colnames(df) <- c('PD_hex', 'vR_car','vt_car','vt_car_y', 'vR_opt', 'vt_opt')

dfm <- reshape2::melt(df)
dfm <- na.omit(dfm)

windows(width = 6, height = 4)
# png("comparison_T4b_PD.png", width = 1200, height = 900, pointsize=24)
ggplot(dfm, aes(factor(Var2), value) ) + 
  geom_jitter(colour='black',size = 1,height = 0, width = 0.1) +
  # scale_colour_manual(values = col_4_65, guide=FALSE) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.3, y = 0),
               colour = 'black', size = 1, width = .2, na.rm = T ) +
  stat_summary(fun = median, geom = "crossbar",  colour = 'black', width = 0.2, lwd=0.4, position = position_nudge(x = 0.3, y = 0), na.rm = T) +
  coord_cartesian() +
  theme_minimal_grid() +
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, by=45), labels = paste(seq(0, 180, by=45),"°", sep = ''), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  # theme(axis.text.y = element_text(size = 20),axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank(),axis.text.x = element_blank()) +
  labs(title = "", x='', y='') 
# ggsave(filename = "comparison_T4b_all.png")
# dev.off()

colMeans(ang_motion, na.rm = T) %>% round(1)
colMeans(ang_hex, na.rm = T) %>% round(1)
apply(ang_motion, 2, median, na.rm=T)
# apply(ang_hex, 2, median, na.rm=T)

# t,R prediction, T4d ---------------------------------------------------------

# position
d <- ucl_rot_sm
colnames(d) <- c("X","Y","Z")

# -- optimal
LL <- 4
v_data <- RF_lens_T4_pred[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm
v_data <- v_data - sweep(ucl_rot_sm, 1, rowSums(v_data * ucl_rot_sm), '*') #tangent
ii <- ucl_rot_sm[,2] < -0.98
p <- v_data / mean(sqrt(rowSums(v_data[ii,]^2)))
vf_amp <- mean(sqrt(rowSums(v_data[ii,]^2)))

# main loop
R_cErr <- c()
t_cErr <- c()
for (k in 1:dim(maxis_tpxyz)[1]) {
  # rot
  v <- R_gen(maxis_tpxyz[k,3:5], d)
  v <- v / mean(sqrt(rowSums(v[ii,]^2)))
  Err <- angcos_vf(v, p) ## ang
  R_cErr <- c(R_cErr, mean(Err,na.rm = T))

  # trans
  v <- t_gen(maxis_tpxyz[k,3:5], d)
  v <- v / mean(sqrt(rowSums(v[ii,]^2)))
  Err <- angcos_vf(v, p) ## ang
  t_cErr <- c(t_cErr, mean(Err, na.rm = T))
}
# mean square error
R_errLand <- data.frame(phi = maxis_tpxyz[,2], theta = maxis_tpxyz[,1], error = R_cErr)
t_errLand <- data.frame(phi = maxis_tpxyz[,2], theta = maxis_tpxyz[,1], error = t_cErr)

# tR prediction, [theta, phi, x, y ,z]
pred_R <- maxis_tpxyz[which.min(R_cErr),]
pred_t <- maxis_tpxyz[which.min(t_cErr),]

pred_t_T4d <- pred_t[3:5]
pred_R_T4d <- pred_R[3:5]


# -- v-axis
hex_hvf_eye <- matrix(ncol = 3, nrow = nrow(eyemap))
for (j in 1:nrow(nb_ind)) {
  if (sum(complete.cases(nb_ind[j,])) == 7) {
    # eye
    bt <- ucl_rot_sm[nb_ind[j,3], ] - ucl_rot_sm[nb_ind[j,6], ]
    hex_hvf_eye[j,] <- -bt
  }
}
v_data <- hex_hvf_eye
p <- v_data / max(na.omit(sqrt(rowSums(v_data^2))))

# main loop
R_cErr <- c()
t_cErr <- c()
for (k in 1:dim(maxis_tpxyz)[1]) {
  # rot
  v <- R_gen(maxis_tpxyz[k,3:5], d)
  v <- v / mean(sqrt(rowSums(v[ii,]^2)))
  Err <- angcos_vf(v, p) ## ang
  R_cErr <- c(R_cErr, mean(Err,na.rm = T))
  
  # trans
  v <- t_gen(maxis_tpxyz[k,3:5], d)
  v <- v / mean(sqrt(rowSums(v[ii,]^2)))
  Err <- angcos_vf(v, p) ## ang
  t_cErr <- c(t_cErr, mean(Err, na.rm = T))
}
# mean square error
R_errLand <- data.frame(phi = maxis_tpxyz[,2], theta = maxis_tpxyz[,1], error = R_cErr)
t_errLand <- data.frame(phi = maxis_tpxyz[,2], theta = maxis_tpxyz[,1], error = t_cErr)

# tR prediction, [theta, phi, x, y ,z]
pred_R <- maxis_tpxyz[which.min(R_cErr),]
pred_t <- maxis_tpxyz[which.min(t_cErr),]

pred_t_hex <- pred_t[3:5]
pred_R_hex <- pred_R[3:5]

# cont. ED Fig.7E, angular error wrt cardinal/optimized motion VS lattice---------------------------

# - lattice v-axis
ii_6nb <- rep(0,nrow(nb_ind))
v_hex <- matrix(ncol = 6, nrow = nrow(nb_ind))

for (j in 1:nrow(nb_ind)) {
  if (sum(complete.cases(nb_ind[j,])) == 7) {
    ii_6nb[j] <- 1
    # eye
    bt <- ucl_rot_sm[nb_ind[j,3], ] - ucl_rot_sm[nb_ind[j,6], ] 
    v_hex[j,] <- c(ucl_rot_sm[nb_ind[j,1], ], -bt)
  }
}

# - T4d
LL <- 4
v_data <- RF_lens_T4_pred[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm
v_data <- v_data - sweep(ucl_rot_sm, 1, rowSums(v_data * ucl_rot_sm), '*') #tangent
ii <- ucl_rot_sm[,2] < -0.98
v_data <- v_data / mean(sqrt(rowSums(v_data[ii,]^2)))

# - optic flow
pt0 <- as.matrix(ucl_rot_sm) #pts at t=0 on screen

# of
vt_car <- t_gen(c(0,0,1), pt0)*0.1
vt_car_my <- t_gen(c(0,-1,0), pt0)*0.1
vR_car <- R_gen(c(-1,0,0), pt0)*0.1
vt_opt <- t_gen(pred_t_T4d, pt0)*0.1
vR_opt <- R_gen(pred_R_T4d, pt0)*0.1
vt_hex <- t_gen(pred_t_hex, pt0)*0.1
vR_hex <- R_gen(pred_R_hex, pt0)*0.1

# # - diff in ang
ang_motion <- matrix(ncol = 6, nrow = nrow(ucl_rot_sm))
for (j in 1:nrow(ucl_rot_sm)) {
  if (ii_6nb[j] == 1) {
    ang_motion[j,] <- c(angcos(v_data[j,], v_hex[j,4:6]),
                        angcos(v_data[j,], vR_car[j,]),
                        angcos(v_data[j,], vt_car[j,]),
                        angcos(v_data[j,], vt_car_my[j,]),
                        angcos(v_data[j,], vR_opt[j,]),
                        angcos(v_data[j,], vt_opt[j,])   )
  }
}

# PLOT
df <- ang_motion 
colnames(df) <- c('PD_hex', 'vR_car','vt_car','vt_car_my', 'vR_opt', 'vt_opt')

dfm <- reshape2::melt(df)
dfm <- na.omit(dfm)

# windows(width = 3, height = 4)
windows(width = 6, height = 4)
# png("comparison_T4d_PD.png", width = 1200, height = 900, pointsize=24)
ggplot(dfm, aes(factor(Var2), value) ) + 
  geom_jitter(colour='black',size = 1,height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.3, y = 0),
               colour = 'black', size = 1, width = .2, na.rm = T ) +
  stat_summary(fun = median, geom = "crossbar",  colour = 'black', width = 0.2, lwd=0.4, position = position_nudge(x = 0.3, y = 0), na.rm = T) +
  coord_cartesian() +
  theme_minimal_grid() +
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, by=45), labels = paste(seq(0, 180, by=45),"°", sep = ''), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  # theme(axis.text.y = element_text(size = 20),axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank(),axis.text.x = element_blank()) +
  labs(title = "", x='', y='') 
# ggsave(filename = "comparison_T4d_PD.png")
# dev.off()

colMeans(ang_motion, na.rm = T) %>% round(1)
colMeans(ang_hex, na.rm = T) %>% round(1)
apply(ang_motion, 2, median, na.rm=T)
apply(ang_hex, 2, median, na.rm=T)


# Fig.4G, ED Fig.7F, partition by best fit to canonical/cardinal  motion ----------------------------------

# position
d <- ucl_rot_sm
colnames(d) <- c("X","Y","Z")

## ##
LL <- 2
v_data <- RF_lens_T4_pred[,(3*(LL-1)+1):(3*(LL-1)+3)]- ucl_rot_sm 
v_data <- v_data - sweep(ucl_rot_sm, 1, rowSums(v_data * ucl_rot_sm), '*') #tangent
ii <- ucl_rot_sm[,2] < -0.98
p <- v_data / mean(sqrt(rowSums(v_data[ii,]^2))) # normalize to 
vf_amp <- mean(sqrt(rowSums(v_data[ii,]^2)))

# # all cardinal motions
# ls_car <- list(t_gen(c(+1,0,0), d),
#                t_gen(c(-1,0,0), d),
#                t_gen(c(0,+1,0), d),
#                t_gen(c(0,-1,0), d),
#                t_gen(c(0,0,+1), d),
#                t_gen(c(0,0,-1), d),
#                R_gen(c(+1,0,0), d),
#                R_gen(c(-1,0,0), d),
#                R_gen(c(0,+1,0), d),
#                R_gen(c(0,-1,0), d),
#                R_gen(c(0,0,+1), d),
#                R_gen(c(0,0,-1), d) )

## ## relevant to T4b
pal_3 <- pal_TR[c(1,2,6)] 
ls_car <- list(#t_gen(c(+1,0,0), d),
               t_gen(c(-1,0,0), d),
               # t_gen(c(0,+1,0), d),
               t_gen(c(0,-1,0), d),
               # t_gen(c(0,0,+1), d),
               # t_gen(c(0,0,-1), d),
               # R_gen(c(+1,0,0), d),
               # R_gen(c(-1,0,0), d),
               # R_gen(c(0,+1,0), d),
               # R_gen(c(0,-1,0), d),
               # R_gen(c(0,0,+1), d),
               R_gen(c(0,0,-1), d) )


## ## relevant to T4d
# pal_3 <- pal_TR[c(3,4,5)]
# ls_car <- list(
#   # t_gen(c(+1,0,0), d),
#   # t_gen(c(-1,0,0), d),
#   # t_gen(c(0,+1,0), d),
#   # t_gen(c(0,-1,0), d),
#   t_gen(c(0,0,+1), d),
#   # t_gen(c(0,0,-1), d),
#   # R_gen(c(+1,0,0), d),
#   R_gen(c(-1,0,0), d),
#   # R_gen(c(0,+1,0), d),
#   R_gen(c(0,-1,0), d)
#   # R_gen(c(0,0,+1), d),
#   # R_gen(c(0,0,-1), d) 
#   )


angdiff <- matrix(ncol = length(ls_car), nrow = nrow(d))
for (j in 1:length(ls_car)) {
  # dd <- sqrt(rowSums((v - p)^2))
  angdiff[,j] <- angcos_vf(p, ls_car[[j]]) 
}

# -- eyal plot
star <- Mollweide(cart2sph2tp(matrix(pred_t_T4b,ncol=3))[, c('t','p')])
star2 <- Mollweide(cart2sph2tp(matrix(pred_R_T4b,ncol=3))[, c('t','p')])
# 
# subsampling
xy_sub <- seq(-20,20,by=3)
indlens_sub <- lens_ixy[lens_ixy[,2] %in% xy_sub & lens_ixy[,3] %in% xy_sub, 1]
ind_sub <- match(indlens_sub, eyemap[,2])

# plot
ang <- (180-angdiff)
# ang <- angdiff
amp <- 0.1
df_1 <- data.frame(cbind(ucl_rot_Mo[,c('xM','yM')],
                         ucl_rot_Mo[,c('xM')] + ang[,1]/180 * amp * cos(-30/180*pi),
                         ucl_rot_Mo[,c('yM')] + ang[,1]/180 * amp * sin(-30/180*pi)))
colnames(df_1) <- c('x','y','xend','yend')
df_2 <- data.frame(cbind(ucl_rot_Mo[,c('xM','yM')],
                         ucl_rot_Mo[,c('xM')] + ang[,2]/180 * amp * cos(210/180*pi),
                         ucl_rot_Mo[,c('yM')] + ang[,2]/180 * amp * sin(210/180*pi)))
colnames(df_2) <- c('x','y','xend','yend')
df_3 <- data.frame(cbind(ucl_rot_Mo[,c('xM','yM')],
                         ucl_rot_Mo[,c('xM')] + ang[,3]/180 * amp * cos(90/180*pi),
                         ucl_rot_Mo[,c('yM')] + ang[,3]/180 * amp * sin(90/180*pi)))
colnames(df_3) <- c('x','y','xend','yend')
df_scale <- data.frame(cbind(-.5, -c(.2,.4,.6,.8), -.5+seq(45,180,by=45)/180 *amp, -c(.2,.4,.6,.8)))
colnames(df_scale) <- c('x','y','xend','yend')

plt <- plt_Mo34 +
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  geom_segment(data = df_1[ind_sub,], aes(x = x,y = y, xend = xend,yend = yend), colour=pal_3[1],size =2, lineend = "round") +  
  geom_segment(data = df_2[ind_sub,], aes(x = x,y = y, xend = xend,yend = yend), colour=pal_3[2],size =2, lineend = "round") +  
  geom_segment(data = df_3[ind_sub,], aes(x = x,y = y, xend = xend,yend = yend), colour=pal_3[3],size =2, lineend = "round") +  
  geom_segment(data = df_scale, aes(x = x,y = y, xend = xend,yend = yend), colour="black",size =2, lineend = "round") +  
  # scale_color_manual(values = pal_3,guide= guide_legend(title="ang"), na.value="gray") +
  annotate("text", x = star[1], y = star[2], label = "+", size = 9) +
  annotate("text", x = star2[1], y = star2[2], label = "x", size = 9) +
  labs(title = "ang diff color")
windows(width = 12, height = 8)
# pdf(paste("T4", letters[LL], "_partition_eyal.pdf", sep = ''), width = 8.5, height = 4.5)
plt
# dev.off()

# ## ## -- Mollweide
# LL = m
# vv <- ucl_rot_sm + 0.5*(RF_lens_T4_pred[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)
# vv_Mo <- sweep(vv, 1, sqrt(rowSums(vv^2)), '/')
# vv_Mo <- cart2sph2tp(vv_Mo)
# vv_Mo <- Mollweide(vv_Mo[,c('t', 'p')])
# colnames(vv_Mo) <- c('xM','yM')
# rownames(vv_Mo) <- rownames(ucl_rot_sm)
# 
# # -- PLOT
# df_pos <- data.frame(ucl_rot_Mo)
# df_arrow <- data.frame(cbind(ucl_rot_Mo[,c('xM','yM')], vv_Mo[,c('xM','yM')]) )
# colnames(df_arrow) <- c('x','y','xend','yend')
# 
# # plyr::count(df_arrow$col)
# df_arrow$col <- apply(angdiff, 1, which.min) %>% factor()
# 
# star <- Mollweide(cart2sph2tp(matrix(pred_t_T4b,ncol=3))[, c('t','p')])
# star2 <- Mollweide(cart2sph2tp(matrix(pred_R_T4b,ncol=3))[, c('t','p')])
# 
# plt <- plt_Mo34 +
#   geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
#   geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=col),size =1) +  
#   scale_color_manual(values = pal_3,guide= guide_legend(title="ang"), na.value="gray") +
#   annotate("text", x = star[1], y = star[2], label = "+", size = 9) +
#   annotate("text", x = star2[1], y = star2[2], label = "x", size = 9) +
#   labs(title = "ang diff color")
# windows(width = 12, height = 8)
# # pdf(paste("T4", letters[LL], "_partition.pdf", sep = ''), width = 8.5, height = 4.5)
# plt
# # dev.off()

# ## ## -- error for each motion
# k <- 2
# df_pos <- data.frame(ucl_rot_Mo)
# df_pos$size <- angdiff[, k]
# 
# plt <- plt_Momin +
#   geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
#   geom_point(data=df_pos, aes(x=xM, y=yM, size = size),colour = pal_3[k]) +
#   scale_size(name= 'ang diff', breaks = seq(30, 180, by=30), labels = seq(30, 180, by=30)) +
#   theme(legend.position = c(.9, .9) ) +
#   labs(title = "angdiff")
# windows(width = 12, height = 8)
# # pdf(paste("T4", letters[LL], "_angdiff_", k, ".pdf", sep = ''), width = 8.5, height = 4.5)
# plt
# # dev.off()


# ## ## -- error for each motion, heat map
# k <- 3
# df_pos <- data.frame(ucl_rot_Mo)
# df_pos$quan <- angdiff[, k]
# 
# rg <- c(0.01, 180)
# 
# plt <- plt_Mo34 +
#   geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
#   geom_point(data=df_pos, aes(x = xM, y = yM, colour = quan), size = 2) +
#   # scale_color_gradientn(colours = pal_heat1, values = scales::rescale(rg), limits=range(rg), oob=scales::squish,
#   #                       breaks= rg, labels= rg, guide = guide_colorbar(title = "ang") ) +
#   scale_color_gradient(low = "blue", high= "gray90", limits=rg, oob=scales::squish, trans='log',
#                        breaks= rg, labels= rg, guide = guide_colorbar(title = "angle")) +
#   theme(legend.position = c(.9, .9) ) +
#   labs(title = "angdiff")
# windows(width = 12, height = 8)
# # pdf(paste("T4", letters[LL], "_angdiff_heat_", k, ".pdf", sep = ''), width = 8.5, height = 4.5)
# plt
# # dev.off()



# Fig.4H,  3d optimal axes ------------------------------------------------------------
# image tracer + line art
# atan2(pred_t_T4b[2], pred_t_T4b[1]) /pi*180 + 180

eq_xyz <- cbind(cos(seq(0,360,by=5)/180*pi),
               sin(seq(0,360,by=5)/180*pi),
               0)
para_n <- cbind(sin(pi/4)*cos(seq(0,360,by=5)/180*pi),
                sin(pi/4)*sin(seq(0,360,by=5)/180*pi),
                cos(pi/4))
para_s <- cbind(sin(3*pi/4)*cos(seq(0,360,by=5)/180*pi),
                sin(3*pi/4)*sin(seq(0,360,by=5)/180*pi),
                cos(3*pi/4))
fr_xyz <- cbind(cos(seq(-90,270,by=5)/180*pi),
                0,
                sin(seq(-90,270,by=5)/180*pi) )
fr_45 <- fr_xyz %*% matrix(c(cos(pi/4),sin(pi/4),0,
                  -sin(pi/4), cos(pi/4), 0,
                  0,0,1), ncol = 3, byrow = T) 
fr_m45 <- fr_xyz %*% matrix(c(cos(pi/4),-sin(pi/4),0,
                  sin(pi/4), cos(pi/4), 0,
                  0,0,1), ncol = 3, byrow = T)
side_xyz <- cbind(0,
                  sin(seq(-90,270,by=5)/180*pi),
                  cos(seq(-90,270,by=5)/180*pi)
                   )

t_a <- pred_t_T4b
t_a[2] <- -t_a[2]
r_a <- pred_R_T4b
r_a[2] <- -r_a[2]
t_c <- pred_t_T4d
t_c[2] <- -t_c[2]
r_c <- pred_R_T4d
r_c[2] <- -r_c[2]

nopen3d()
par3d('windowRect' = c(50,50,1550,1550))
# par3d('windowRect' = c(50,50,550,550))
spheres3d(0,0,0,0.999, col='grey100', alpha=1, lit=F)
arrow3d(pred_t_T4b*1.1, -pred_t_T4b*1.7, theta = pi / 18, n = 8, s=0.05, col = pal_T4[2], type = "rotation", lwd=2, lit=F)
arrow3d(pred_R_T4b*1.1, -pred_R_T4b*1.75, theta = pi / 9, n = 8, s=0.05, col = pal_T4[2], type = "rotation", lwd=2, lit=F)
arrow3d(pred_t_T4d*1.1, -pred_t_T4d*1.7, theta = pi / 18, n = 8, s=0.05, col = pal_T4[4], type = "rotation", lwd=2, lit=F)
arrow3d(pred_R_T4d*1.1, -pred_R_T4d*1.75, theta = pi / 9, n = 8, s=0.05, col = pal_T4[4], type = "rotation", lwd=2, lit=F)

# arrow3d(t_a*1.6, -t_a*1.7, theta = pi / 18, n = 8, s=0.05, col = pal_T4[1], type = "rotation", lwd=2, lit=F)
# arrow3d(r_a*1.6, -r_a*1.75, theta = pi / 9, n = 8, s=0.05, col = pal_T4[1], type = "rotation", lwd=2, lit=F)
# arrow3d(t_c*1.6, -t_c*1.7, theta = pi / 18, n = 8, s=0.05, col = pal_T4[3], type = "rotation", lwd=2, lit=F)
# arrow3d(r_c*1.6, -r_c*1.75, theta = pi / 9, n = 8, s=0.05, col = pal_T4[3], type = "rotation", lwd=2, lit=F)
# lines3d(rbind(-pred_t_T4b*1.6, pred_t_T4b*1.7), lwd=2)
# lines3d(rbind(-pred_R_T4b*1.6, pred_R_T4b*1.75), lwd=2)
lines3d(eq_xyz, lwd=2)
lines3d(para_n, lwd=2)
lines3d(para_s, lwd=2)
lines3d(fr_xyz, lwd=2)
lines3d(fr_45, lwd=2)
lines3d(fr_m45, lwd=2)
lines3d(side_xyz, lwd=2)
# lines3d(rbind(c(1,0,0), c(1.6, 0, 0)), lwd=2)
# lines3d(rbind(c(0,0,1), c(0, 0, 1.5)), lwd=2)
# lines3d(rbind(c(0,-1,0), c(0, -1.8, 0)), lwd=2)
# points3d(ucl_rot_sm[rownames(ucl_rot_sm) %in% lens_ixy[lens_ixy[,3] == 0, 1], ], size=12, col=pal_axes[1])
# points3d(ucl_rot_sm[rownames(ucl_rot_sm) %in% lens_ixy[lens_ixy[,2] == -1, 1], ], size=12, col=pal_axes[2])
# points3d(ucl_rot_sm[rownames(ucl_rot_sm) %in% vaxis_gen(-1), ], size=12, col=pal_axes[3])
# points3d(ucl_rot_sm[c(ind_Up_ucl, ind_Down_ucl), ], size=12, col=pal_axes[4])
# points3d(ucl_rot_sm, size=7, col='gray50')
# arrow3d(c(1,0,0), c(1.6, 0, 0), theta = pi / 18, n = 8, col = "gray30", type = "rotation", lwd=2)
segments3d(rbind(c(1,0,0), c(1.8, 0, 0)),col = "gray30",lwd=3)
# text3d(c(1.3, -0.1, -0.3), texts = 'front', cex=3)
# arrow3d(c(0,0,1), c(0, 0, 1.5), theta = pi / 12, n = 8, col = "gray30", type = "rotation", lwd=2)
segments3d(rbind(c(0,0,1), c(0, 0, 1.5)),col = "gray30",lwd=3)
# text3d(c(-0.2, -0.2, 1.5), texts = 'up', cex=3)
# arrow3d(c(0,-1,0), c(0, -1.8, 0), theta = pi / 18, n = 8, col = "gray30", type = "rotation", lwd=2)
segments3d(rbind(c(0,-1,0), c(0, -1.6, 0)),col = "gray30",lwd=3)
# text3d(c(0.2, -1.8, 0), texts = 'side', cex=3)
# planes3d(0,1,0, 0, alpha = 0.2, lit=F)
# planes3d(0,0,1, 0, alpha = 0.2, lit=F)
rgl.viewpoint(fov=0,zoom=0.7, userMatrix= rotationMatrix(-60/180*pi,1,0,0) %*% rotationMatrix(-55/180*pi,0,0,1))
# rgl.viewpoint(fov=0,zoom=0.7, userMatrix= rotationMatrix(-80/180*pi,1,0,0) %*% rotationMatrix(-47/180*pi,0,0,1))

# rgl.snapshot("3d_guide.png")
# rgl.snapshot("3d_opt_noarrow.png")
# rgl.snapshot("3d_opt.png")
# rgl.postscript("3d_opt.pdf", "pdf") #no