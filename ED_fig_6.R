# see Fig_4.R for ED Fig.6A, ED Fig.6F

# ED Fig.6A, Mercator projection ---------------------------------------------

LL = 2
vv <- ucl_rot_sm + 0.5*(RF_lens_T4_pred_sm[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm)

vv_Mer <- cart2Mercator(vv)

# H2
df_H2 <- cbind(H2p_Merc, H2d_Merc) %>% as.data.frame()
df_H2 <- df_H2[!is.na(df_H2[,1]),]
colnames(df_H2) <- c('x','y','xend','yend')
df_H2 <- as.data.frame(df_H2)

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


# ED Fig.6B,  shear angle  vs T4  --------------------------------------------------

# ang betw (bot->top) and (back->front), that is, +v
nb_ang_ucl <- matrix(ncol = 2, nrow = nrow(eyemap))
hex_hvf_eye <- matrix(ncol = 3, nrow = nrow(eyemap))
hex_hvf_med <- matrix(ncol = 3, nrow = nrow(eyemap))
for (j in 1:nrow(nb_ind)) {
  if (sum(complete.cases(nb_ind[j,])) == 7) {
    bt <- ucl_rot_sm[nb_ind[j,3], ] - ucl_rot_sm[nb_ind[j,6], ]
    bf <- colMeans(ucl_rot_sm[nb_ind[j,4:5],]) - colMeans(ucl_rot_sm[nb_ind[j,c(2,7)],])
    hex_hvf_eye[nb_ind[j,1],] <- bf
    ang <- acos(bt %*% bf / sqrt(sum(bt^2)) / sqrt(sum(bf^2)) ) /pi*180 
    nb_ang_ucl[nb_ind[j,1],] <- c(j, ang)
  }
}

# angles between eye v-axis and T4b
v_data <- (RF_lens_T4_pred_sm[,(3*(LL-1)+1):(3*(LL-1)+3)] - ucl_rot_sm) 

# -- vertical axis
ang_T4_vaxis_eye <- matrix(NA, ncol = 1, nrow = nrow(eyemap))
for (j in 1:nrow(nb_ind)) {
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


# - PLOT 2d Mollweide, angle diff
df_pos <- data.frame(ucl_rot_Mo)
df_pos$quan <- ang_T4b_vaxis_eye - nb_ang_ucl[,2]

quantile(df_pos$quan, c(0.01,0.05,0.25,0.5,0.75,0.95,0.99), na.rm = T)
rg <- c(-15, 0, 15) 

plt <- plt_Momin + 
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  geom_point(data=df_pos, aes(x = xM, y = yM, colour = quan), size = 2) +
  scale_color_gradientn(colours = pal_heat1, values = scales::rescale(rg), limits=range(rg), oob=scales::squish,
                        breaks= rg, labels= rg, guide = guide_colorbar(title = "t4b-v-shear") ) +
  theme(legend.position = c(.9, .9) ) +
  labs(title = "T4b-h on eye")
windows(width = 9, height = 6)
# pdf("t4b_haxis_eye.pdf", width = 8.5, height = 4.5)
plt
# dev.off()

# ED Fig.6C, size, T4b -----------------------------------------------------------

LL <- 2 # choose type

com_xyz <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
v0 <- as.matrix(lens_type[[LL]][, c("x0","y0","z0")])
v1 <- as.matrix(lens_type[[LL]][, c("xd","yd","zd")])

PD_eye <- matrix(ncol = 1, nrow = nrow(v0)) # PD length normalized
for (j in 1:nrow(v0)) {
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

# - plot
df_pos <- data.frame(ucl_rot_Mo)
df_pos$val <- np_pred

range((df_pos$val))
quantile(na.omit(df_pos$val), c(0, 0.05,0.25, 0.5, 0.75,0.95, 1))
df_pos$valgp <- cut(df_pos$val,c(8, 10, 12, 15, 18))

# -- 2d
plt <- ggplot() +
  geom_polygon(data = as.data.frame(bkgd_str_equa), aes(x=xM, y=yM), fill = 'grey90', alpha=1) +
  geom_polygon(data = as.data.frame(bkgd_str_meri), aes(x=xM, y=yM), fill = 'grey90', alpha=1) +
  geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
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

windows(width = 8, height = 4)
# pdf(paste0("PD_T4", letters[LL], "_equ.pdf"), width = 8, height = 4)
thr <- cos((90- str_hw)/180*pi)
ii <- com_xyz[,3] < thr & com_xyz[,3] > -thr
com_rtp <- cart2sphZ(com_xyz)
x <- 360 - com_rtp[ii,3]/pi*180
x <- if_else(x>180, x-360, x)
plot(x, PD_eye[ii], xlab ="", pch=15, ylab='', main = paste("T4", letters[LL], " PD along equ", sep = ""),
     ylim = c(5, 25), xlim= c(-10, 150), xaxt='n',yaxt='n')
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


# ED Fig.6G, regression med_xyz to ucl_rot_sm -----------------------------------------------

load("data/eyemap.RData")
N <- nrow(med_xyz)
iis <- seq(1, N)

nnn <- c(1+8, 1+8+8*2, 1+8+8*2+8*3, 1+8+8*2+8*3+8*4)
nnn <- c(1+8, 1+8+8*2)

sy_fix_ll_ls <- list()
sy_fix_lc_ls <- list()
for (ind_nnn in 1:length(nnn)) {
  nbhd_N <- nnn[ind_nnn]
  sy_ll <- matrix(ncol = 3, nrow = length(iis))
  sy_lc <- matrix(ncol = 3, nrow = length(iis))
  for (j in 1:length(iis)) {
    nb_ind <- sweep(med_xyz, 2, med_xyz[iis[j],], '-')^2 %>% rowSums() %>% order() %>% head(nbhd_N)
    sx <- data.frame(mc.x = med_xyz[iis[j],1], mc.y = med_xyz[iis[j],2], mc.z = med_xyz[iis[j],3])
    # -- kernel regression
    for (k in 1:3) {
      npdata <- data.frame(mc = med_xyz[nb_ind,], ec = ucl_rot_sm[nb_ind,k]) 
      
      bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'fixed', regtype= 'll')
      model_np <- npreg(bw)
      sy_ll[j, k] <- predict(model_np, newdata = sx)
      
      bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'fixed', regtype= 'lc')
      model_np <- npreg(bw)
      sy_lc[j, k] <- predict(model_np, newdata = sx)
    }
  }
  sy_fix_ll_ls[[ind_nnn]] <- sy_ll
  sy_fix_lc_ls[[ind_nnn]] <- sy_lc
}

# - cp bwtype
nnn <- c(1+8, 1+8+8*2)
sy_gnn_ls <- list()
sy_ann_ls <- list()
for (ind_nnn in 1:length(nnn)) {
  nbhd_N <- nnn[ind_nnn]
  sy_ann <- matrix(ncol = 3, nrow = length(iis)) #adaptive
  sy_gnn <- matrix(ncol = 3, nrow = length(iis)) # general
  for (j in 1:length(iis)) {
    nb_ind <- sweep(med_xyz, 2, med_xyz[iis[j],], '-')^2 %>% rowSums() %>% order() %>% head(nbhd_N)
    sx <- data.frame(mc.x = med_xyz[iis[j],1], mc.y = med_xyz[iis[j],2], mc.z = med_xyz[iis[j],3])
    # -- kernel regression
    for (k in 1:3) {
      npdata <- data.frame(mc = med_xyz[nb_ind,], ec = ucl_rot_sm[nb_ind,k]) 
      
      bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'generalized_nn', regtype= 'll')
      model_np <- npreg(bw)
      sy_gnn[j, k] <- predict(model_np, newdata = sx)
      
      bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
      model_np <- npreg(bw)
      sy_ann[j, k] <- predict(model_np, newdata = sx)
    }
  }
  sy_gnn_ls[[ind_nnn]] <- sy_gnn
  sy_ann_ls[[ind_nnn]] <- sy_ann
}


# - with aux
med_xyz <- med_xyz_aux
ucl_rot_sm <- ucl_rot_aux

nnn <- c(1+8, 1+8+8*2)
sy_fix_ll_aux_ls <- list()
sy_fix_lc_aux_ls <- list()
for (ind_nnn in 1:length(nnn)) {
  nbhd_N <- nnn[ind_nnn]
  
  sy_ll_aux <- matrix(ncol = 3, nrow = length(iis))
  sy_lc_aux <- matrix(ncol = 3, nrow = length(iis))
  for (j in 1:length(iis)) {
    nb_ind <- sweep(med_xyz, 2, med_xyz[iis[j],], '-')^2 %>% rowSums() %>% order() %>% head(nbhd_N)
    sx <- data.frame(mc.x = med_xyz[iis[j],1], mc.y = med_xyz[iis[j],2], mc.z = med_xyz[iis[j],3])
    # -- kernel regression
    for (k in 1:3) {
      npdata <- data.frame(mc = med_xyz[nb_ind,], ec = ucl_rot_sm[nb_ind,k]) 
      
      bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'fixed', regtype= 'll')
      model_np <- npreg(bw)
      sy_ll_aux[j, k] <- predict(model_np, newdata = sx)
      
      bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'fixed', regtype= 'lc')
      model_np <- npreg(bw)
      sy_lc_aux[j, k] <- predict(model_np, newdata = sx)
    }
  }
  sy_fix_ll_aux_ls[[ind_nnn]] <- sy_ll_aux
  sy_fix_lc_aux_ls[[ind_nnn]] <- sy_lc_aux
}

# # SAVE
# save(nnn,
#      sy_fix_lc_ls, sy_fix_ll_ls,
#      sy_gnn_ls, sy_ann_ls,
#      sy_fix_lc_aux_ls, sy_fix_ll_aux_ls,
#      file = "data/sy_ls.RData")

load(file = "data/sy_ls.RData")

# PLOT
for (ind_nnn in 1:2) {
  nopen3d()
  # adaptive
  points3d(ucl_rot_sm[iis,], size=8)
  sy <- sy_ann_ls[[ind_nnn]]
  # push out a bit to see better
  sy <- sweep(sy, 1, sqrt(rowSums(sy^2)), '/')
  sy <- sweep(sy, 1, 1.005, '*')
  points3d(sy, col='red', size = 8)
  linemat <- matrix(t(cbind(ucl_rot_sm[iis,], sy)), ncol = 3, byrow = T)
  segments3d(linemat, color = "grey")
  title3d(paste("aux", nnn[ind_nnn]))
}

# SAVE
# rgl.snapshot("np_test_nnb9_ann.png")


# - hist
dev.new()
# pdf("np_compare_dist.pdf")
dt <- (sy_ann_ls[[1]] - ucl_rot_sm[iis,])^2 %>% rowSums() %>% sqrt()
range(dt)
hist(dt, breaks = seq(0, 0.5,by=0.001), xlim = c(0, 0.05),ylim=c(0,200),
     xlab ="pairwise distance", ylab='counts', plot = T,  xaxt = "n")
lines(c(5/180*pi/2,5/180*pi/2), c(0,100))
axis(1, at = seq(0, 3, by=0.5)/180*pi, labels = paste0(seq(0,3, by =0.5),"°"))
# dev.off()

# - compare
stats <- matrix(ncol = 2, nrow = 0)
sy_ls <- c(sy_fix_lc_ls, sy_fix_ll_ls, sy_gnn_ls, sy_ann_ls, sy_fix_lc_aux_ls, sy_fix_ll_aux_ls)
for (j in 1:length(sy_ls)) {
  dd <- sqrt(rowSums((ucl_rot_sm[iis,] - sy_ls[[j]])^2))
  stats <- rbind(stats, c(mean(dd), sd(dd)))
}

df <- as.data.frame(cbind(seq(1,length(sy_ls)), stats))
dev.new()
ggplot(df, aes(x=V1, y=V2)) + 
  geom_errorbar(aes(ymin=V2-V3, ymax=V2+V3), width=.2) +
  geom_line() +
  geom_point()
