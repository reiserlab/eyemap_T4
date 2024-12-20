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


# ED Fig.6B, size, T4b -----------------------------------------------------------

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



# regression med_xyz to ucl_rot_sm -----------------------------------------------

load("data/eyemap.RData")
N <- nrow(med_xyz)
# iis <- sample(N, 50)
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
      
      # bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'generalized_nn', regtype= 'll')
      # model_np <- npreg(bw)
      # sy_gnn[j, k] <- predict(model_np, newdata = sx)
      
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

# save.image("C:/Users/artxz/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_eyemap/R_eyemap/ws_np_test.RData")


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
  # sy <- sy_fix_ll_ls[[ind_nnn]]
  sy <- sy_ann_ls[[ind_nnn]]
  # sy <- sy_fix_ll_aux_ls[[ind_nnn]]
  # push out a bit
  sy <- sweep(sy, 1, sqrt(rowSums(sy^2)), '/')
  sy <- sweep(sy, 1, 1.005, '*')
  points3d(sy, col='red', size = 8)
  linemat <- matrix(t(cbind(ucl_rot_sm[iis,], sy)), ncol = 3, byrow = T)
  segments3d(linemat, color = "grey")
  
  # general
  # points3d(ucl_rot_sm)
  # # sy <- sy_fix_lc_ls[[ind_nnn]]
  # sy <- sy_gnn_ls[[ind_nnn]]
  # # sy <- sy_fix_lc_aux_ls[[ind_nnn]]
  # points3d(sy, col='blue', size = 7)
  # linemat <- matrix(t(cbind(ucl_rot_sm[iis,], sy)), ncol = 3, byrow = T)
  # segments3d(linemat, color = "grey")
  
  # title3d(paste("fixed", nnn[ind_nnn]))
  # title3d(paste("nn", nnn[ind_nnn]))
  title3d(paste("aux", nnn[ind_nnn]))
  # # title3d(round(sse,4))
}

# SAVE
# rgl.snapshot("np_test_nnb9_ann.png")


# - hist
dev.new()
# pdf("np_compare_dist.pdf")
dt <- (sy_ann_ls[[1]] - ucl_rot_sm[iis,])^2 %>% rowSums() %>% sqrt()
range(dt)
# hh <- hist(dt, breaks = seq(0, 0.5,by=0.02), plot = F)
# plot(hh$mids, hh$density, type='l', bty='n', xlim = c(0, 0.5),ylim=c(0,4), xaxt='n',yaxt='n', xlab ="asp. ratio", ylab='')
# axis(1, at = seq(0,0.15, by =0.05), labels = seq(0,0.15, by =0.05))
hist(dt, breaks = seq(0, 0.5,by=0.001), xlim = c(0, 0.05),ylim=c(0,200),
     xlab ="pairwise distance", ylab='counts', plot = T,  xaxt = "n")
lines(c(5/180*pi/2,5/180*pi/2), c(0,100))
axis(1, at = seq(0, 3, by=0.5)/180*pi, labels = paste0(seq(0,3, by =0.5),"°"))
# dev.off()

# - stats compare
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
