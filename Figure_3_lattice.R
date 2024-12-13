
# T4 com etc --------------------------------------------------------------

LL <- 2
T4com <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
T4dir <- dir_type[[LL]][, c('rsxd','rsyd','rszd', 'rsx0','rsy0','rsz0')] %>% as.matrix() # reverse st optic flow dir
T4direye <- lens_type[[LL]][, c('xd','yd','zd', 'x0','y0','z0')] %>% as.matrix() # reverse st optic flow dir


# ED Fig.5C, ED Fig.5D, ED Fig.5E, ED Fig.5J, asp ratio------------------------------------

# - histo
# dt <- rowMeans(nb_dist_ucl[, 2:5]) / rowMeans(nb_dist_ucl[, 6:7]) #2021
dt <- rowMeans(nb_dist_ucl[, c(2,4,5,7)]) / rowMeans(nb_dist_ucl[, c(3,6)]) #2023
range(na.omit(dt))
quantile(na.omit(dt),0.98)
dev.new()
# pdf("hist_latt_ratio_eye.pdf")
hh <- hist(dt, breaks = seq(0.5, 4.5,by=0.1), plot = F)
plot(hh$mids, hh$counts, type='l', bty='n', xlim = c(0.5, 2), xaxt='n',yaxt='n', xlab ="asp. ratio", ylab='counts')
axis(1, at = seq(0.5,2, by =0.5), labels = seq(0.5,2, by =0.5))
# dev.off()

# -- med
# dt <- rowMeans(nb_dist_med[, 2:5]) / rowMeans(nb_dist_med[, 6:7]) #2021
dt <- rowMeans(nb_dist_med[, c(2,4,5,7)]) / rowMeans(nb_dist_med[, c(3,6)]) #2023
range(na.omit(dt))
dev.new()
# pdf("hist_latt_ratio_med.pdf", height=3)
hh <- hist(dt, breaks = seq(0.5, 1.5,by=0.1), plot=F)
plot(hh$mids, hh$counts, type='l', bty='n', lwd=3, xlim = c(0.5, 1.5), xaxt='n',yaxt='n', xlab ="asp. ratio", ylab='')
axis(1, at = seq(0.5,1.5, by =0.5), labels = seq(0.5,1.5, by =0.5))
# dev.off()

# -- compare
dev.new()
# pdf("hist_latt_ratio_eyeVSmed.pdf")
# dt <- rowMeans(nb_dist_ucl[, 2:5]) / rowMeans(nb_dist_ucl[, 6:7]) #2021
dt <- rowMeans(nb_dist_ucl[, c(2,4,5,7)]) / rowMeans(nb_dist_ucl[, c(3,6)]) #2023
hh <- hist(dt, breaks = seq(0.5, 3.5,by=0.1), plot = F)
plot(hh$mids, hh$density, type='l', bty='n', xlim = c(0.5, 2),ylim=c(0,4), xaxt='n',yaxt='n', xlab ="asp. ratio", ylab='')
axis(1, at = seq(0.5,2, by =0.5), labels = seq(0.5,2, by =0.5))

# dt <- rowMeans(nb_dist_med[, 2:5]) / rowMeans(nb_dist_med[, 6:7]) #2021
dt <- rowMeans(nb_dist_med[, c(2,4,5,7)]) / rowMeans(nb_dist_med[, c(3,6)]) #2023
hh <- hist(dt, breaks = seq(0.4, 1.5,by=0.1), plot=F)
lines(hh$mids, hh$density, lty=2, lwd=3 )
# dev.off()

# -- PLOT density
windows(width = 4, height = 4)
# pdf(paste("asp_comp.pdf", sep = ""), width = 4, height = 4)
# dt <- rowMeans(nb_dist_ucl[, 2:5]) / rowMeans(nb_dist_ucl[, 6:7]) #2021
dt <- rowMeans(nb_dist_ucl[, c(2,4,5,7)]) / rowMeans(nb_dist_ucl[, c(3,6)]) #2023
dd <- density(dt, from= 0.4, to= 3.1, bw='SJ', na.rm = T)
plot(dd$x, dd$y, type='l', bty='n', col='black', lwd=3, xlim = c(0.5, 2), ylim=c(0,5), xaxt='n',yaxt='n',
     xlab ="", ylab='',main = '')
# dt <- rowMeans(nb_dist_med[, 2:5]) / rowMeans(nb_dist_med[, 6:7]) #2021
dt <- rowMeans(nb_dist_med[, c(2,4,5,7)]) / rowMeans(nb_dist_med[, c(3,6)]) #2023
dd <- density(dt, from= 0.4, to= 3.1, bw='SJ', na.rm = T)
points(dd$x, dd$y, type='l', bty='n', col='blue', lwd=3)
axis(1, at = seq(0.5, 2, by =0.5), labels = seq(0.5, 2, by =0.5), cex.axis = 1.5 )
axis(2, at = c(0, 1/512/diff(dd$x[1:2])*((3.1-0.4)/(2-0.5)) ), labels = c('0', '1'), cex.axis = 1.5 )
# dev.off()


# - 2d, asp ratio
df_pos <- data.frame(ucl_rot_Mo)

# -- eye
rr <- rowMeans(nb_dist_ucl[, c(2,4,5,7)]) / rowMeans(nb_dist_ucl[, c(3,6)]) #2023
quantile(na.omit(rr), c(0, 0.05,0.5,0.95, 1))
scalecuts <- c(0.8, 1.6)
df_pos$ratio <- rr

# --- Mollweide
plt <- plt_Mo34 + 
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  geom_point(data=df_pos, aes(x = xM, y = yM, colour = ratio), size = 2) +
  scale_color_gradientn(colours = pal_heat1, values = scales::rescale(c(scalecuts[1], 1, scalecuts[2])), limits=scalecuts, oob=scales::squish,
                        breaks= scalecuts, labels= scalecuts, guide = guide_colorbar(title = "ratio") ) +
  # scale_color_gradientn(colours = c(munsell::mnsl("5PB 7/10"),"white",munsell::mnsl("5R 7/10")), values = scales::rescale(c(-0.3, 0, 0.5)),  limits=log10(c(0.5, 3.5)),
  #                       breaks= log10(c(0.5,1,3.5)), labels=c(-0.3,0,0.5), guide = guide_colorbar(title = "ratio [log10]") ) +
  # scale_color_gradientn(colours = c(munsell::mnsl("5PB 7/10"),munsell::mnsl("5Y 9/2"),munsell::mnsl("5R 7/10")), values = scales::rescale(c(0.5, 1, 3)),  limits=c(0.5,3.5),
  #                       breaks= c(0.5,1,2,3), labels= c(0.5,1,2,3), guide = guide_colorbar(title = "ratio") ) +
  theme(legend.position = c(.9, .9) ) +
  labs(title = "aspect ratio eye")
# for (k in 1:3) {
#   # 2023
#   j <- match(ind_eghex[k], nb_ind[,1])
#   pt <- ucl_rot_Mo[nb_ind[j,], ]
#   pt <- pt[c(2,3,4,5,6,7,2), ]
#   plt <- plt + geom_path(data=as.data.frame(pt), aes(x=xM, y=yM), colour='gray50',lwd=0.5)
# }

windows(width = 12, height = 8)
# pdf("asp_eye.pdf", width = 8.5, height = 4.5)
plt
# dev.off()

# --- Mercator
df <- ucl_rot_Merc
df$quan <- rr

plt <- plt_Mer +
  geom_path(data = cmer_Merc, aes(x=x, y=y), colour = pal_axes[3], lwd=1) +
  geom_point(data=df, aes(x=x, y=y, colour = quan), size = 1) +
  scale_color_gradientn(colours = pal_heat1, values = scales::rescale(c(scalecuts[1], 1, scalecuts[2])), limits=scalecuts, oob=scales::squish,
                        breaks= scalecuts, labels= scalecuts, guide = guide_colorbar(title = "ratio") ) +
  labs(title = "asp ratio eye")

windows(width = 8, height = 10)
# pdf("asp_eye_Merc.pdf",width = 8, height = 10)
plt
# dev.off()

# -- med
rr <- rowMeans(nb_dist_med[, c(2,4,5,7)]) / rowMeans(nb_dist_med[, c(3,6)]) #2023
quantile(na.omit(rr), c(0, 0.01,0.05,0.5,0.95,0.99, 1))
scalecuts <- c(0.8, 1.6)
df_pos$ratio <- rr

# Mollweide
plt <- plt_Mo34 + 
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  geom_point(data=df_pos, aes(x = xM, y = yM, colour = ratio), size = 2) +
  scale_color_gradientn(colours = pal_heat1, values = scales::rescale(c(scalecuts[1], 1, scalecuts[2])), limits=scalecuts, oob=scales::squish,
                        breaks= scalecuts, labels= scalecuts, guide = guide_colorbar(title = "ratio") ) +
  theme(legend.position = c(.9, .9) ) +
  labs(title = "aspect ratio med")
# for (k in 1:3) {
#   # 2023
#   j <- match(ind_eghex[k], nb_ind[,1])
#   pt <- ucl_rot_Mo[nb_ind[j,], ]
#   pt <- pt[c(2,3,4,5,6,7,2), ]
#   plt <- plt + geom_path(data=as.data.frame(pt), aes(x=xM, y=yM), colour='gray50',lwd=0.5)
# }

windows(width = 12, height = 8)
# pdf("asp_med.pdf", width = 8.5, height = 4.5)
plt
# dev.off()


# #  example asp ratio, v-axis aligned ---------------------------------
# 
# 
# # 226, 77, 558, ratio <1 =1 >2
# # c(470, 707, 364) %in% match(vaxis_gen(1), eyemap[,2])
# hex_eg <- list()
# for (k in 1:3) {
#   # 2021
#   # j <- c(470, 707, 364)[k]
#   # # pt <- ucl_rot_sm[nb_ind[j,], ]
#   # zz <- cross3D((pt[3,]-pt[1,]), (pt[2,]-pt[1,])) #pointing inwards
#   # 2023
#   j <- match(ind_eghex[k], nb_ind[,1])
#   pt <- ucl_rot_sm[nb_ind[j,], ]
#   # pt <- med_xyz[nb_ind[j,], ]
#   zz <- cross3D((pt[2,]-pt[1,]), (pt[3,]-pt[1,])) #pointing inwards
#   pc <- prcomp(pt)
# 
#   if (pc$rotation[,3] %*% zz < 0) {
#     pc$rotation <- -pc$rotation
#   }
#   if (c(cross3D(pc$rotation[,3], pc$rotation[,1])) %*% pc$rotation[,2] < 0) {
#     pc$rotation[,2] <- -pc$rotation[,2]
#   }
# 
#   pt <- sweep(pt, 2, pc$center) %*% pc$rotation
#   pt[,3] <- 0
#   va <- (pt[6,] - pt[7,]) #vertical axis of hex
#   ang <- acos(va %*% c(0,1,0) / sqrt(sum(va^2)))
#   if (c(cross3D(va, c(0,1,0))) %*% c(0,0,1) < 0) {
#     ang <- -ang
#   }
#   pt <- pt %*% matrix(c(cos(ang), sin(ang), 0,
#                         -sin(ang), cos(ang), 0,
#                         0,0,1), ncol=3, byrow=T)
#   hex_eg[[k]] <- pt
# }
# 
# # PLOT
# j <- 3
# pt <- hex_eg[[j]]
# windows(width = 8, height = 8)
# # pdf(file = paste("hex_eg_", c("DD","D","C")[j], ".pdf", sep=''))
# plot(pt, pch=16, axes=F, ann=F, cex=2, asp=1)
# #2021
# # lines(pt[c(2,3,6,4,5,7,2),1:2], lwd=3)
# # # lines(pt[c(6,7), 1:2], lty = "dotted", lwd=2)
# # a1 <- colMeans(pt[4:5,]) #right midpoint
# # a2 <- colMeans(pt[2:3,]) # left
# # va <- pt[6,] - pt[7,] # vertical vec upwards
# # 2023
# lines(pt[c(2,3,4,5,6,7,2),1:2], lwd=3) #2023
# a1 <- colMeans(pt[4:5,]) #right midpoint
# a2 <- colMeans(pt[c(2,7),]) # left
# va <- pt[3,] - pt[6,] # vertical vec upwards
# # arrows(a1[1], a1[2], a2[1], a2[2], lwd=2)
# title(paste(c("DD","D","C")[j],
#             " asp= ",
#             # round(mean(sqrt(rowSums(sweep(pt[2:5,],2,pt[1,])^2))) / mean(sqrt(rowSums(sweep(pt[6:7,],2,pt[1,])^2))), 2), #2021
#             round(mean(sqrt(rowSums(sweep(pt[c(2,4,5,7),],2,pt[1,])^2))) / mean(sqrt(rowSums(sweep(pt[c(3,6),],2,pt[1,])^2))), 2), #2023
#             sep = '')
# )
# # dev.off()


# ED Fig.5C, example asp ratio, local meridian aligned----------------------------------------------------------

nb_coord <- rbind(c(-5/180*pi,0),
                  c(+5/180*pi,0),
                  c(0, -5/180*pi),
                  c(0, +5/180*pi) )

hex_eg <- list()
for (k in 1:3) {
  j <- match(ind_eghex[k], nb_ind[,1]) # 2023
  # ind_T4 <- rowSums(sweep(T4com, 2, med_xyz[j,])^2) %>% which.min()
  
  # eye
  pt <- ucl_rot_sm[nb_ind[j,], ] 
  tp <- cart2sphZ(pt[1,,drop=F])[,2:3]
  rtp <- cbind(1, sweep(nb_coord,2,tp,'+'))
  xyz <- sph2cartZ(rtp)
  
  # # ### med
  # ind_nb <- rowSums(sweep(med_xyz, 2, as.numeric(T4com[ind_T4, ]))^2) %>% order() %>% head(1+6+12)
  # pt <- med_xyz[nb_ind[j,], ] # med
  # # np
  # xyz <- eyemap_np(ucl_rot_sm[ind_nb,], med_xyz[ind_nb,], xyz)
  
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
  xyz[,3] <- 0
  pt[,3] <- 0
  
  # align to meridian
  va <- (xyz[1,] - xyz[2,]) #vertical axis of hex
  ang <- acos(va %*% c(0,1,0) / sqrt(sum(va^2)))
  if (c(cross3D(va, c(0,1,0))) %*% c(0,0,1) < 0) {
    ang <- -ang
  }
  pt <- pt %*% matrix(c(cos(ang), sin(ang), 0,
                        -sin(ang), cos(ang), 0,
                        0,0,1), ncol=3, byrow=T)
  
  hex_eg[[k]] <- pt
}

# PLOT
j <- 3
pt <- hex_eg[[j]]
windows(width = 8, height = 8)
# pdf(file = paste("hex_eg_", c("DD","D","C")[j], ".pdf", sep=''))
plot(pt, pch=16, axes=F, ann=F, cex=2, asp=1)
lines(pt[c(2,3,4,5,6,7,2),1:2], lwd=3)
a1 <- colMeans(pt[4:5,]) #right midpoint
a2 <- colMeans(pt[c(2,7),]) # left
va <- pt[3,] - pt[6,] # vertical vec upwards
title(paste(c("DD","D","C")[j],
            " asp= ", 
            round(mean(sqrt(rowSums(sweep(pt[c(2,4,5,7),],2,pt[1,])^2))) / mean(sqrt(rowSums(sweep(pt[c(3,6),],2,pt[1,])^2))), 2), #2023
            sep = '')
)
# dev.off()


# Fig.3J, ED Fig.5H, shear angle -------------------------------------------------------------
# ang betw (bot->top) and (back->front), that is, +v

nb_ang_ucl <- matrix(ncol = 2, nrow = nrow(eyemap))
nb_ang_med <- matrix(ncol = 2, nrow = nrow(eyemap))
hex_hvf_eye <- matrix(ncol = 3, nrow = nrow(eyemap))
hex_hvf_med <- matrix(ncol = 3, nrow = nrow(eyemap))
for (j in 1:nrow(nb_ind)) {
  if (sum(complete.cases(nb_ind[j,])) == 7) {
    # eye
    # 2023
    bt <- ucl_rot_sm[nb_ind[j,3], ] - ucl_rot_sm[nb_ind[j,6], ]
    bf <- colMeans(ucl_rot_sm[nb_ind[j,4:5],]) - colMeans(ucl_rot_sm[nb_ind[j,c(2,7)],])
    hex_hvf_eye[nb_ind[j,1],] <- bf
    ang <- acos(bt %*% bf / sqrt(sum(bt^2)) / sqrt(sum(bf^2)) ) /pi*180 
    nb_ang_ucl[nb_ind[j,1],] <- c(j, ang)
    # med
    # 2023
    bt <- med_xyz[nb_ind[j,3], ] - med_xyz[nb_ind[j,6], ]
    bf <- colMeans(med_xyz[nb_ind[j,4:5],]) - colMeans(med_xyz[nb_ind[j,c(2,7)],])
    hex_hvf_med[nb_ind[j,1],] <- bf
    ang <- acos(bt %*% bf / sqrt(sum(bt^2)) / sqrt(sum(bf^2)) ) /pi*180 
    nb_ang_med[nb_ind[j,1],] <- c(j, ang)
  }
}

# # DEBUG med
# nopen3d()
# points3d(med_xyz)
# for (j in seq(-12,12, by=3)) {
#   ii <- vaxis_gen(j)
#   points3d(med_xyz[match(ii, eyemap[,2]),], col='orange', size = 9)
# }
# for (j in 1:nrow(nb_ind)) {
#   if (sum(complete.cases(nb_ind[j,])) == 7) {
#     # med
#     bt <- med_xyz[nb_ind[j,6], ] - med_xyz[nb_ind[j,7], ]
#     bf <- colMeans(med_xyz[nb_ind[j,2:3],]) - colMeans(med_xyz[nb_ind[j,4:5],])
#     arrow3d(med_xyz[nb_ind[j,1],], med_xyz[nb_ind[j,1],]+bf, col='blue',theta= pi/12, n=6, type= "rotation")
#   }
# }

# # - histo
# range(na.omit(nb_ang_med[,2]))
# dev.new()
# # pdf("hist_latt_ang_med.pdf")
# hh <- hist(nb_ang_med[,2], breaks = seq(30, 150, by=5), plot = F)
# plot(hh$mids, hh$counts, type='l', bty='n', xlim = c(30, 150), xaxt='n',yaxt='n', xlab ="angle [deg]", ylab='counts')
# axis(1, at = seq(30,150, by =30), labels = seq(30,150, by =30))
# # dev.off()
# 
# range(na.omit(nb_ang_ucl[,2]))
# dev.new()
# # pdf("hist_latt_ang_eye.pdf")
# hh <- hist(nb_ang_ucl[,2], breaks = seq(30, 150, by=5), plot = F)
# plot(hh$mids, hh$counts, type='l', bty='n', xlim = c(30, 150), xaxt='n',yaxt='n', xlab ="angle [deg]", ylab='counts')
# axis(1, at = seq(30, 150, by =30), labels =  paste(seq(30, 150, by =30), "°", sep = '') )
# # dev.off()

# - density
windows(width = 5, height = 4)
# pdf(paste("den_ang_eye.pdf", sep = ""), width = 5, height = 4)
dd <- density(nb_ang_ucl[,2], from= 30, to= 150, bw='SJ', na.rm = T)
plot(dd$x, dd$y, type='l', bty='n', col='black', lwd=3, xlim = c(30, 150), xaxt='n',yaxt='n',
     xlab ="", ylab='',main = '')
axis(1, at = seq(30, 150, by =30), labels = paste(seq(30, 150, by =30), "°", sep = ''), cex.axis = 1.5 )
axis(2, at = c(0, 1/512/diff(dd$x[1:2]) ), labels = c('0', '1'), cex.axis = 1.5 )
# dev.off()

# - 2d plot
df_pos <- data.frame(ucl_rot_Mo)
### ###
df_pos$quan <- nb_ang_ucl[,2]
quantile(na.omit(df_pos$quan), c(0, 0.01,0.5,0.99, 1))
### ###
df_pos$quan <- nb_ang_med[,2]

# -- Mollweide
# rg <- c(45, 90, 135)
rg <- c(60, 90, 120) # 20240701
plt <- plt_Mo34 + 
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  geom_point(data=df_pos, aes(x = xM, y = yM, colour = quan), size = 2) +
  scale_color_gradientn(colours = pal_heat1, values = scales::rescale(rg),
                        limits=range(rg), oob=scales::squish,breaks= rg,
                        labels= rg, guide = guide_colorbar(title = "skewness") ) +
  # scale_color_gradientn(colours = c('gray95', 'gray45','gray0'), values = scales::rescale(rg),
  #                       limits=range(rg), oob=scales::squish,breaks= rg,
  #                       labels= rg, guide = guide_colorbar(title = "skewness") ) +
  theme(legend.position = c(.9, .9) ) +
  # annotate("text", df_pos$xM, df_pos$yM, label=round(df_pos$quan,0)) +
  labs(title = "skewness eye")
for (k in 1:3) {
  # 2023
  j <- match(ind_eghex_2[k], nb_ind[,1])
  pt <- ucl_rot_Mo[nb_ind[j,], ]
  pt <- pt[c(2,3,4,5,6,7,2), ]
  plt <- plt + geom_path(data=as.data.frame(pt), aes(x=xM, y=yM), colour='gray50',lwd=0.5)
}

windows(width = 9, height = 6)
# pdf("skew_eye.pdf", width = 8.5, height = 4.5)
# pdf("skew_med.pdf", width = 8.5, height = 4.5)
plt
# dev.off()

# -- mercator
df <- ucl_rot_Merc
df$quan <- nb_ang_ucl[,2]
df$quan <- nb_ang_med[,2]

# rg <- c(45, 90, 135)
rg <- c(60, 90, 120)
plt <- plt_Mer +
  geom_path(data = cmer_Merc, aes(x=x, y=y), colour = pal_axes[3], lwd=1) +
  geom_point(data=df, aes(x=x, y=y, colour = quan), size = 1) +
  scale_color_gradientn(colours = pal_heat1, values = scales::rescale(rg), limits=range(rg), oob=scales::squish,
                        breaks= rg, labels= rg, guide = guide_colorbar(title = "skewness") ) +
  theme(legend.position = c(.9, .9) ) +
  labs(title = "skewness med")

windows(width = 8, height = 10)
# pdf("skew_eye_Merc.pdf",width = 8, height = 10)
# pdf("skew_med_Merc.pdf",width = 8, height = 10)
plt
# dev.off()

# Fig.3J,  example skew angle, local meridian aligned----------------------------------------------------------

nb_coord <- rbind(c(-5/180*pi,0),
                  c(+5/180*pi,0),
                  c(0, -5/180*pi),
                  c(0, +5/180*pi) )

hex_eg <- list()
for (k in 1:3) {
  # j <- ind_eghex_2[k] #2021
  j <- match(ind_eghex_2[k], nb_ind[,1]) # 2023
  
  # eye
  pt <- ucl_rot_sm[nb_ind[j,], ] 
  tp <- cart2sphZ(pt[1,,drop=F])[,2:3]
  rtp <- cbind(1, sweep(nb_coord,2,tp,'+'))
  xyz <- sph2cartZ(rtp)
  
  # # ### med
  # ind_nb <- rowSums(sweep(med_xyz, 2, as.numeric(T4com[ind_T4, ]))^2) %>% order() %>% head(1+6+12)
  # pt <- med_xyz[nb_ind[j,], ] # med
  # # np
  # xyz <- eyemap_np(ucl_rot_sm[ind_nb,], med_xyz[ind_nb,], xyz)
  
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
  xyz[,3] <- 0
  pt[,3] <- 0
  
  # align to meridian
  va <- (xyz[1,] - xyz[2,]) #vertical axis of hex
  ang <- acos(va %*% c(0,1,0) / sqrt(sum(va^2)))
  if (c(cross3D(va, c(0,1,0))) %*% c(0,0,1) < 0) {
    ang <- -ang
  }
  pt <- pt %*% matrix(c(cos(ang), sin(ang), 0,
                        -sin(ang), cos(ang), 0,
                        0,0,1), ncol=3, byrow=T)
  
  hex_eg[[k]] <- pt
}


# PLOT
j <- 3
pt <- hex_eg[[j]]
windows(width = 8, height = 8)
# pdf(file = paste("hex_eg_alpha_eye_", c("D","C","V")[j], ".pdf", sep=''))
# pdf(file = paste("hex_eg_alpha_med_", c("D","C","V")[j], ".pdf", sep=''))
# pdf(file = paste("hex_eg_alpha_med_2_", c("D","C","V")[j], ".pdf", sep=''))
plot(pt, pch=16, axes=F, ann=F, cex=2, asp=1)
# 2023
lines(pt[c(2,3,4,5,6,7,2),1:2], lwd=3)
lines(pt[c(3,6), 1:2], lty = "dotted", lwd=2)
a1 <- colMeans(pt[c(2,7),]) #right midpoint 
a2 <- colMeans(pt[c(4,5),]) # left
va <- pt[3,] - pt[6,] # vertical vec upwards
arrows(a1[1], a1[2], a2[1], a2[2], lwd=2)
title(paste(c("D","C","V")[j],
            " alpha= ",
            round(acos((a2-a1) %*%  va / sqrt(sum(va^2)) / sqrt(sum((a2-a1)^2))) /pi*180, 1),
            sep = '')
)
# dev.off()


# #  example skewness, v-axis aligned -----------------------------------------------
# 
# # c(87, 463, 513) %in% match(vaxis_gen(11), eyemap[,2]) #eye
# # c(219, 328, 582) %in% match(vaxis_gen(-6), eyemap[,2]) #med
# hex_eg <- list()
# for (k in 1:3) {
#   # ###
#   # j <- c(87, 463, 513)[k] #eye
#   # pt <- ucl_rot_sm[nb_ind[j,], ]
#   # ###
#   j <- c(219, 328, 582)[k] #med
#   pt <- med_xyz[nb_ind[j,], ]
# 
#   zz <- cross3D((pt[3,]-pt[1,]), (pt[2,]-pt[1,])) #pointing inwards
#   pc <- prcomp(pt)
# 
#   if (pc$rotation[,3] %*% zz < 0) {
#     pc$rotation <- -pc$rotation
#   }
#   if (c(cross3D(pc$rotation[,3], pc$rotation[,1])) %*% pc$rotation[,2] < 0) {
#     pc$rotation[,2] <- -pc$rotation[,2]
#   }
# 
#   pt <- sweep(pt, 2, pc$center) %*% pc$rotation
#   pt[,3] <- 0
#   va <- (pt[6,] - pt[7,]) #vertical axis of hex
#   ang <- acos(va %*% c(0,1,0) / sqrt(sum(va^2)))
#   if (c(cross3D(va, c(0,1,0))) %*% c(0,0,1) < 0) {
#     ang <- -ang
#   }
#   pt <- pt %*% matrix(c(cos(ang), sin(ang), 0,
#                         -sin(ang), cos(ang), 0,
#                         0,0,1), ncol=3, byrow=T)
#   hex_eg[[k]] <- pt
# }
# 
# # PLOT
# j <- 3
# pt <- hex_eg[[j]]
# windows(width = 8, height = 8)
# # pdf(file = paste("hex_eg_alpha_eye_", c("D","C","V")[j], ".pdf", sep=''))
# # pdf(file = paste("hex_eg_alpha_med_", c("D","C","V")[j], ".pdf", sep=''))
# # pdf(file = paste("hex_eg_alpha_med_2_", c("D","C","V")[j], ".pdf", sep=''))
# plot(pt, pch=16, axes=F, ann=F, cex=2, asp=1)
# lines(pt[c(2,3,6,4,5,7,2),1:2], lwd=3)
# lines(pt[c(6,7), 1:2], lty = "dotted", lwd=2)
# a1 <- colMeans(pt[4:5,]) #right midpoint
# a2 <- colMeans(pt[2:3,]) # left
# va <- pt[6,] - pt[7,] # vertical vec upwards
# arrows(a1[1], a1[2], a2[1], a2[2], lwd=2)
# title(paste(c("D","C","V")[j],
#             " alpha= ",
#             round(acos((a2-a1) %*%  va / sqrt(sum(va^2)) / sqrt(sum((a2-a1)^2))) /pi*180, 1),
#             sep = '')
# )
# # dev.off()



# # ang betw +h and parallels, +v and meridians -------------------------------------------------------------
# 
# nb_ang_para <- matrix(ncol = 2, nrow = nrow(eyemap))
# nb_ang_meri <- matrix(ncol = 2, nrow = nrow(eyemap))
# for (j in 1:nrow(nb_ind)) {
#   if (sum(complete.cases(nb_ind[j,])) == 7) {
#     # eye, meri
#     bt <- ucl_rot_sm[nb_ind[j,6], ] - ucl_rot_sm[nb_ind[j,7], ] #bottom -> top
#     bt <- bt - c(bt %*% ucl_rot_sm[nb_ind[j,1], ]) * ucl_rot_sm[nb_ind[j,1], ]
#     v_meri <- c(0,0,1) - ucl_rot_sm[nb_ind[j,1], ]
#     v_meri <- v_meri - c(v_meri %*% ucl_rot_sm[nb_ind[j,1], ]) * ucl_rot_sm[nb_ind[j,1], ]
#     
#     # para
#     bf <- colMeans(ucl_rot_sm[nb_ind[j,2:3],]) - colMeans(ucl_rot_sm[nb_ind[j,4:5],])
#     bf <- bf - c(bf %*% ucl_rot_sm[nb_ind[j,1], ]) * ucl_rot_sm[nb_ind[j,1], ]
#     v_para <- matrix(c(cos(0.1), -sin(0.1), 0,
#                        sin(0.1), cos(0.1), 0,
#                        0, 0, 1), ncol = 3, byrow = T) %*% matrix(ucl_rot_sm[nb_ind[j,1], ], ncol = 1)
#     v_para <- c(v_para)
#     v_para <- v_para - c(v_para %*% ucl_rot_sm[nb_ind[j,1], ]) * ucl_rot_sm[nb_ind[j,1], ]
#     
#     # meri
#     ang <- acos(bf %*% v_meri / sqrt(sum(bf^2)) / sqrt(sum(v_meri^2)) ) /pi*180 
#     nb_ang_meri[j,] <- c(j, ang-90)
#     # para
#     ang <- acos(bt %*% v_para / sqrt(sum(bt^2)) / sqrt(sum(v_para^2)) ) /pi*180 
#     nb_ang_para[j,] <- c(j, 90-ang)
#   }
# }
# 
# 
# # - histo
# range(na.omit(nb_ang_para[,2]))
# range(na.omit(nb_ang_meri[,2]))
# # dev.new()
# # # pdf("hist_latt_ang_med.pdf")
# # hh <- hist(nb_ang_med[,2], breaks = seq(30, 150, by=5), plot = F)
# # plot(hh$mids, hh$counts, type='l', bty='n', xlim = c(30, 150), xaxt='n',yaxt='n', xlab ="angle [deg]", ylab='counts')
# # axis(1, at = seq(30,150, by =30), labels = seq(30,150, by =30))
# # # dev.off()
# # 
# # range(na.omit(nb_ang_ucl[,2]))
# # dev.new()
# # # pdf("hist_latt_ang_eye.pdf")
# # hh <- hist(nb_ang_ucl[,2], breaks = seq(30, 150, by=5), plot = F)
# # plot(hh$mids, hh$counts, type='l', bty='n', xlim = c(30, 150), xaxt='n',yaxt='n', xlab ="angle [deg]", ylab='counts')
# # axis(1, at = seq(30, 150, by =30), labels =  paste(seq(30, 150, by =30), "°", sep = '') )
# # # dev.off()
# 
# # - 2d plot
# df_pos <- data.frame(ucl_rot_Mo)
# ## ##
# df_pos$quan <- nb_ang_para[,2]
# quantile(na.omit(df_pos$quan), c(0, 0.01,0.5,0.99, 1))
# ## ##
# df_pos$quan <- nb_ang_meri[,2]
# 
# # -- Mollweide
# rg <- c(-75, 0, 75)
# rg <- c(-35, 0, 35)
# plt <- plt_Mo34 + 
#   geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
#   geom_point(data=df_pos, aes(x = xM, y = yM, colour = quan), size = 2) +
#   scale_color_gradientn(colours = pal_heat1, values = scales::rescale(rg), limits=range(rg), oob=scales::squish,
#                         breaks= rg, labels= rg, guide = guide_colorbar(title = "angle") ) +
#   theme(legend.position = c(.9, .9) ) +
#   # labs(title = "eye +h vs para")
#   labs(title = "eye +v vs meri")
# # for (k in 1:3) {
# #   # j <- c(510, 662, 366)[k]
# #   j <- c(500, 554, 322)[k]
# #   pt <- ucl_rot_Mo[nb_ind[j,], ]
# #   pt <- pt[c(2,3,6,4,5,7,2), ]
# #   plt <- plt + geom_path(data=as.data.frame(pt), aes(x=xM, y=yM), colour='gray50',lwd=0.5)
# # }
# 
# windows(width = 9, height = 6)
# # pdf("skew_eye.pdf", width = 8.5, height = 4.5)
# # pdf("skew_med.pdf", width = 8.5, height = 4.5)
# plt
# # dev.off()

# # const size(rather, on a plane) on curved surface --------------------------------------------
# 
# c1 <- Mollweide(cbind(45, seq(-60, 90, by = 2)))
# colnames(c1) <- c('xM','yM')
# c2 <- Mollweide(cbind(65, seq(-60, 90, by = 2)))
# colnames(c2) <- c('xM','yM')
# c3 <- Mollweide(cbind(25, seq(-60, 90, by = 2)))
# colnames(c3) <- c('xM','yM')
# m1 <- Mollweide(cbind(seq(0, 90, by = 2), 0))
# colnames(m1) <- c('xM','yM')
# m2 <- Mollweide(cbind(seq(0, 90, by = 2), 20))
# colnames(m2) <- c('xM','yM')
# m3 <- Mollweide(cbind(seq(0, 90, by = 2), 40))
# colnames(m3) <- c('xM','yM')
# m4 <- Mollweide(cbind(seq(0, 90, by = 2), -20))
# colnames(m4) <- c('xM','yM')
# 
# # aL <- arcLength(c(45,0)/180*pi, c(45,20)/180*pi) /2*sqrt(3)
# aL <- 20 /180*pi * sin(45/180*pi)
# ang1 <- aL / sin(45/180*pi) /pi*180
# ang2 <- aL / sin(65/180*pi) /pi*180
# ang3 <- aL / sin(25/180*pi) /pi*180
# p1 <- Mollweide(cbind(45,c(-ang1,0,ang1,ang1*2,ang1*3)))
# colnames(p1) <- c('xM','yM')
# p2 <- Mollweide(cbind(65,c(-ang2,0,ang2,ang2*2)+ang2/2) )
# colnames(p2) <- c('xM','yM')
# p3 <- Mollweide(cbind(25,c(-ang3,0,ang3,ang3*2)+ang3/2) )
# colnames(p3) <- c('xM','yM')
# 
# df_pos <- data.frame(rbind(p1,p2,p3))
# 
# windows(width = 12, height = 8)
# # pdf("const_size_onS2.pdf", width = 9.5, height = 6.5, useDingbats = F)
# ggplot(df_pos, aes(x = xM, y = yM)) +
#   geom_point(size = 10) +
#   scale_color_gradientn(colours = c(munsell::mnsl("5PB 7/10"),"white",munsell::mnsl("5R 7/10")), values = scales::rescale(c(-90, 0, 60)), limits=c(-90, 60),
#                         breaks= seq(-90, 60,by=30), labels=seq(-90, 60,by=30), guide = guide_colorbar(title = "angle [deg]")) +
#   geom_path(data = as.data.frame(c1), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(c2), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(c3), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(m1), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(m2), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(m3), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(m4), aes(x=xM, y=yM), colour = 'grey50') +
#   ylab("elevation") +
#   xlab("azimuth") +
#   theme_minimal() +
#   scale_x_continuous(limits = c(- 2*sqrt(2)/2, 2*sqrt(2)), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0, sqrt(2)), breaks = c(0,1.5), labels = c(0,-1.5), expand = c(0, 0)) + # set +y as above eq
#   # theme(legend.position="none") +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14)) +
#   labs(title = "const size on sphere") +
#   coord_fixed(ratio=1)
# # dev.off()

# # eye-med difference ------------------------------------------------------
# 
# # -asp
# # - 2d, asp ratio
# df_pos <- data.frame(ucl_rot_Mo)
# 
# # eye - med
# df_pos$ratio <- log10(rowMeans(nb_dist_ucl[, 2:5]) / rowMeans(nb_dist_ucl[, 6:7])) - log10(rowMeans(nb_dist_med[, 2:5]) / rowMeans(nb_dist_med[, 6:7]))
# 
# windows(width = 12, height = 8)
# # pdf("asp_eyeVSmed.pdf", width = 9.5, height = 6.5)
# ggplot(df_pos, aes(x = xM, y = yM)) +
#   geom_point(aes(colour = ratio), size = 3) +
#   # scale_color_gradientn(colours = c(munsell::mnsl("5PB 7/10"),"white",munsell::mnsl("5Y 7/10")), values = scales::rescale(c(0.5, 1, 3.5)), guide = "colorbar", limits=c(0.5, 3.5)) +
#   scale_color_gradientn(colours = c(munsell::mnsl("5PB 7/10"),"white",munsell::mnsl("5R 7/10")), values = scales::rescale(c(-0.4, 0, 0.7)),  limits=log10(c(0.4, 5)),
#                         breaks= log10(c(0.4, 1, 5)), labels=c(-0.4, 0, 0.7), guide = guide_colorbar(title = "ratio [log10]") ) +
#   geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
#   # geom_path(data = mer10, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
#   ylab("elevation") +
#   xlab("azimuth") +
#   theme_minimal() +
#   scale_x_continuous(limits = c(- 2*sqrt(2)/2, 2*sqrt(2)), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(-sqrt(2), sqrt(2)), breaks = c(-1.5,0,1.5), labels = c(-1.5,0,1.5), expand = c(0, 0)) + # set +y as above eq
#   # theme(legend.position="none") +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14)) +
#   labs(title = "aspect ratio eye-med") +
#   coord_fixed(ratio=1)
# # dev.off()
# 
# 
# # -ang
# df_pos <- data.frame(ucl_rot_Mo)
# df_pos$quan <- nb_ang_ucl[,2] - nb_ang_med[,2]
# 
# windows(width = 12, height = 8)
# # pdf("ang_latt_eyeVSmed.pdf", width = 9.5, height = 6.5, useDingbats = F)
# ggplot(df_pos, aes(x = xM, y = yM)) +
#   geom_point(aes(colour = quan), size = 3) +
#   scale_color_gradientn(colours = c(munsell::mnsl("5PB 7/10"),"white",munsell::mnsl("5R 7/10")), values = scales::rescale(c(-90, 0, 60)), limits=c(-90, 60),
#                         breaks= seq(-90, 60,by=30), labels=seq(-90, 60,by=30), guide = guide_colorbar(title = "angle [deg]")) +
#   geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
#   geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
#   ylab("elevation") +
#   xlab("azimuth") +
#   theme_minimal() +
#   scale_x_continuous(limits = c(- 2*sqrt(2)/2, 2*sqrt(2)), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(-sqrt(2), sqrt(2)), breaks = c(-1.5,0,1.5), labels = c(-1.5,0,-1.5), expand = c(0, 0)) + # set +y as above eq
#   # theme(legend.position="none") +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14)) +
#   labs(title = "angle with vertrical axis eye-med") +
#   coord_fixed(ratio=1)
# # dev.off()
