# see Fig_4.R for ED Fig.7A, ED Fig.7B, ED Fig.7E, ED Fig.7F


# ED Fig.7D, size, T4d on eye --------------------------------------------------

LL <- 4 # choose type

# - eye
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


# ED Fig.7G, compare to Henning et al 2022 ----------------------------

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

## ## T4d
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
  geom_segment(data=df_arrow_H, aes(x = x,y = y, xend = xend,yend = yend), colour='black',size =1) +
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =1) +
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

