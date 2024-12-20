# Fig.4A, visual -- neuronal ------------------------------------------------------

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
par3d('windowRect' = c(50,50,1800,1000))
plot3d(boundingbox(bb), alpha=0)

# - med col
points3d(med_xyz_shift[vaxis_gen(clv,ixy = ind_xy),], col=pal_axes[3],size =12)
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

# PD
for (j in ind_eg2[1]) {
  # -- med
  vcom <- as.matrix(dir_type[[LL]][j, c("comx","comy","comz")])
  v0 <- as.matrix(dir_type[[LL]][j, c("rsx0","rsy0","rsz0")])
  vd <- as.matrix(dir_type[[LL]][j, c("rsxd","rsyd","rszd")])
  nbhd_N <- 1+6 +12 # num of nbs
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

# - connecting lines
ii <- lens_ixy[lens_ixy[,2] == 1,]
ii <- ii[order(ii[,3]),]
ii_ucl <- ii[,1]
ii <- ind_xy[ind_xy[,2] == 0,]
ii <- ii[order(ii[,3]),]
ii_med <- rownames(med_xyz[ii[,1],])
linemat <- matrix(t(cbind(
  med_M10_shift[ii_med, ],
  ucl_rot_rhs_shift[ii_ucl[-1], ])
), ncol = 3, byrow = T)
segments3d(linemat, color = 'gray', lwd=3, alpha=0.5)

ii <- lens_ixy[lens_ixy[,2] == 1,]
ii <- ii[order(ii[,3]),]
ii_ucl <- ii[,1]
ii <- ind_xy[ind_xy[,2] == 0,]
ii <- ii[order(ii[,3]),]
ii_med <- rownames(med_xyz[ii[,1],])

# - exclu
xyz_exclu <- rbind(
  ucl_rot_rhs_shift[!(seq(1,nrow(ucl_rot_rhs_shift)) %in% lens_Mi1[,1]),],
  matrix(med_M10_shift[!(seq(1,nrow(utp_Mi1_rot_chiasm)) %in% lens_Mi1[,2]),],nrow=1)
)
pch3d(xyz_exclu, pch=1, cex=0.03, lwd=1.5, col ='gray30')

rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(-90/180*pi,1,0,0) %*% rotationMatrix(125/180*pi,0,0,1)) #inside-out

# rgl.snapshot(filename = paste("T4b_med_eye_excl", ".png", sep = ''))


# Fig.4B, ED Fig.7A, Mollweide, real T4 --------------------------------------------------------------

## ## chose type, T4b=2 or T4d=4
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

# Fig.4C, ED Fig.7B, Mollweide  ------------------------------------------------

## ## chose type, T4b=2 or T4d=4
LL = 2 

vv <- ucl_rot_sm + 0.5*(RF_lens_T4_pred[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)
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


# cont. ED Fig.6A, Mercator projection ---------------------------------------------

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
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =1) +
  geom_segment(data = df_H2, aes(x = x,y = y, xend = xend,yend = yend), colour='red', size=1) +
  labs(title = "")

windows(width = 5, height = 6.5)
plt
# ggsave(paste("T4", letters[LL], "_RF_Mercator_PD.pdf",sep=''), width = 5, height = 6.5)


# cont. Fig.4C, summary of comparison and explanation of bias -----------------------------------------------------------------

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
ang <- ang[c(2,3,5,6,4,1), ]


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

# -- indivisual 
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
df$x <- 0
# reorder, c(2,3,5,6,4,1)
df$x[df$pos == 1] <- 6
df$x[df$pos == 2] <- 1
df$x[df$pos == 3] <- 2
df$x[df$pos == 4] <- 5
df$x[df$pos == 5] <- 3
df$x[df$pos == 6] <- 4

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

# ttest
j <- 3
t.test(df$ang[df$x==j], mu=df2$ang[j])


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

v_data <- (RF_lens_T4_pred_sm[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm) 
v_data_med <- -(RF_med_T4b_pred - med_xyz) # note chiasm

# -- vertical axis
ang_T4_vaxis_eye <- matrix(NA, ncol = 1, nrow = nrow(eyemap))
for (j in 1:nrow(eyemap)) {
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
  theme(legend.position = c(.9, .9) ) +
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
  vv <- rbind(ucl_rot_sm[j,], PD_eye[j,]) #regression T4
  vtcar <- vt_car[j,,drop=F]
  vRcar <- vR_car[j,,drop=F]
  
  # eye
  pt <- ucl_rot_sm[nb_ind[j,], ] 
  tp <- cart2sphZ(pt[1,,drop=F])[,2:3]
  rtp <- cbind(1, sweep(nb_coord,2,tp,'+'))
  xyz <- sph2cartZ(rtp)

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
arrows(vv[1,1], vv[1,2], vv[2,1], vv[2,2], lwd=2, col=pal_T4[LL])
title(paste(c("D","C","V")[j],
            " alpha= ",
            round(acos(vT4 %*%  va / sqrt(sum(va^2)) / sqrt(sum(vT4^2))) /pi*180, 1),
            sep = '')
      )
# dev.off()


# ED Fig.6E, t,R prediction, T4b  -------------------------------------------------

# generate sampling axes
maxis_tpxyz <- axis_sampling(da = 2)


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


# main loop
R_cErr <- c()
t_cErr <- c()
# nd <- sum(!is.na(p[,1]))
ii <- ucl_rot_sm[,2] < -0.98
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


ii <- na.omit(ind_sub)
# position
d <- pos3d[ii,]
nd <- dim(d)[1]

df_pos <- data.frame(ucl_rot_Mo)[ii, ]


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