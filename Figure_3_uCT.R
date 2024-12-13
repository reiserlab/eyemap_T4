# re-run Figure_0 after this script

# ED Fig.5A, ED Fig.5F, 3 axes ------------------------------------------------------------------

# 3D
nopen3d()
points3d(ucl_rot_sm)
points3d(ucl_rot_sm[match(lens_ixy[lens_ixy[,3] == clp,1],eyemap[,2]),], size=10, col=pal_axes[1])
points3d(ucl_rot_sm[match(lens_ixy[lens_ixy[,2] == clq,1],eyemap[,2]),], size=10, col=pal_axes[2])
# points3d(ucl_rot_sm[match(c(ind_Up_lens, ind_Down_lens),eyemap[,2]),], size=10, col=pal_axes[4])
points3d(ucl_rot_sm[vaxis_gen(clv, ind_xy), ], size=10, col=pal_axes[3])
points3d(ucl_rot_sm[haxis_gen(clh, ind_xy), ], size=10, col=pal_axes[4])

# eq
# nopen3d()
# points3d(ucl_rot_sm)
# points3d(ucl_rot_sm[ind_Up_ucl,], size=10, col='red')
# points3d(ucl_rot_sm[ind_Down_ucl,], size=10, col='green')
# planes3d(0,0,1, 0, alpha = 0.5)


# regular hex grid
ind_axis <- vaxis_gen(clv, lens_ixy) 

# - in hex lattice
lens_ixy_hex <- lens_ixy
unit_hex <- t(cbind(c(cos(30/180*pi), sin(30/180*pi)), c(-cos(30/180*pi), sin(30/180*pi))))
lens_ixy_hex[,2:3] <- lens_ixy_hex[,2:3] %*% unit_hex

windows(width = 8, height = 8)
# pdf("axes_lens_hex.pdf")
plot(lens_ixy_hex[, 2:3], xaxt="n", yaxt="n", pch=16, cex=1.25, col='gray') #2021
points(lens_ixy_hex[lens_ixy[,2] == clq,  2:3], pch=16, col=pal_axes[2], cex=1.5)
points(lens_ixy_hex[lens_ixy[,3] == clp, 2:3], pch=16, col=pal_axes[1], cex=1.5)
points(lens_ixy_hex[lens_ixy[,1] %in% haxis_gen(clh), 2:3], pch=16, col=pal_axes[4], cex=1.5)
points(lens_ixy_hex[lens_ixy[,1] %in% ind_axis, 2:3], pch=16, col=pal_axes[3], cex=1.5)
title("4 axes in hex")
# dev.off()


# # - Mollweide, same as "inter-ommatidiam angle..."
# df <- as.data.frame(ucl_rot_Mo)
# plt <- plt_Mo +
#   geom_point(data= as.data.frame(ucl_rot_Mo_right), aes(x = xM, y = yM), colour = 'gray', size=1) +
#   # 2023
#   geom_point(data = df[ind_xy[ind_xy[,2] == clp,1], ], aes(xM, yM), colour = pal_axes[2], size=2) +
#   geom_point(data = df[ind_xy[ind_xy[,3] == clq,1], ], aes(xM, yM), colour = pal_axes[1], size=2) +
#   # geom_point(data = df[c(ind_Up_ucl, ind_Down_ucl), ], aes(xM, yM), colour = pal_axes[4], size=2) +
#   geom_point(data = df[haxis_gen(clh, ind_xy), ], aes(xM, yM), colour = pal_axes[4], size=2) +
#   geom_point(data = df[vaxis_gen(clv, ind_xy), ], aes(xM, yM), colour = pal_axes[3], size=2) +
#   labs(title = "Mollweide 3 axes with eq")
# windows(width = 17, height = 9)
# # pdf("eye_right_Mollweide_ioa_hex_4axes.pdf", width = 8.5, height = 4.5)
# plt
# # dev.off()

# - on Mercator 
df <- as.data.frame(ucl_rot_Merc)
plt <- plt_Mer +
  geom_point(data= as.data.frame(ucl_rot_Merc_right), aes(x,y), colour = 'gray') +
  # 2023
  geom_point(data = df[ind_xy[ind_xy[,2] == clp,1], ], aes(x, y), colour = pal_axes[2], size=2) +
  geom_point(data = df[ind_xy[ind_xy[,3] == clq,1], ], aes(x, y), colour = pal_axes[1], size=2) +
  geom_point(data = df[haxis_gen(clh, ind_xy), ], aes(x, y), colour = pal_axes[4], size=2) +
  geom_point(data = df[vaxis_gen(clv, ind_xy), ], aes(x, y), colour = pal_axes[3], size=2) +
  labs(title = "Mercator 3 axes with eq")
windows(width = 8, height = 10)
# pdf("axes_Mercator.pdf",width = 8, height = 10)
plt
# dev.off()


# ED Fig.4D, Mollweide projection of 2 eyes --------------------------------

ucl_rot_tp <- cart2sphZ(ucl_rot_2eye)[,2:3] %>%
  as_tibble() %>%
  mutate(theta = theta/pi*180, phi = phi/pi*180) %>%
  mutate(phi = if_else(phi > 180, phi - 360, phi)) %>%
  mutate(phi = - phi) %>% #inside out
  as.matrix()
ucl_rot_Mo2 <- Mollweide(ucl_rot_tp)
ucl_rot_Mo2 <- cbind(ucl_rot_Mo2, 1)
ucl_rot_Mo2[!ind_left_cone, 3] <- 2 # left == 1, right == 2
colnames(ucl_rot_Mo2) <- c('xM', 'yM', 'left')

xy_ashape <- ashape(ucl_rot_Mo2[ind_left_cone, 1:2] +
                      matrix(runif(sum(ind_left_cone)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
xy_edge <- xy_ashape$edges[,1:6]
xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
ahull_left <- xy_grid_ahull

xy_ashape <- ashape(ucl_rot_Mo2[!ind_left_cone, 1:2] +
                      matrix(runif(sum(!ind_left_cone)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
xy_edge <- xy_ashape$edges[,1:6]
xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
ahull_right <- xy_grid_ahull

# point distr
df <- as.data.frame(ucl_rot_Mo2)
windows(width = 16, height = 8)
# pdf("viewing_angles_bi.pdf")
plt <- plt_Mo + 
  geom_point(data=df, aes(x=xM, y=yM, colour = factor(left)), size=1) +
  scale_colour_manual(values = pal_lr, labels=c('left','right')) +
  # geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), 
  #              arrow = arrow(length = unit(0.005, "npc"),type = "closed"), size =1) +
  geom_path(data = ahull_left, aes(x1,y1), lwd = 1, colour = pal_lr[1], alpha =0.7) +
  geom_path(data = ahull_right, aes(x1,y1), lwd = 1, colour = pal_lr[2], alpha =0.7) +
  theme(legend.position = c(.9, .9) ) +
  labs(title = '2 eyes')
plt
# dev.off()

# cont. Fig.3F, inter-ommatidia angle --------------------------

ioa_hex <- cbind(ucl_rot_Mo_right, rowMeans(nb_dist_ucl_right[,2:7], na.rm = T)/pi*180) # 2023
ioa_hex <- as.data.frame(ioa_hex)
colnames(ioa_hex) <- c('xM','yM','ioa')
ioa_hex_right <- ioa_hex

ioa_hex <- cbind(ucl_rot_Mo_left, rowMeans(nb_dist_ucl_left[,2:7], na.rm = T)/pi*180) # 2023
ioa_hex <- as.data.frame(ioa_hex)
colnames(ioa_hex) <- c('xM','yM','ioa')
ioa_hex_left <- ioa_hex

# -- pal and breaks
n_lvl <- 8
breaks_ioa <- seq(1,11,length.out = n_lvl+1)

# -- grid for loessm ioa
bkgd_chull <- rbind(bkgd_mer_ww, bkgd_mer_ee[seq(nrow(bkgd_mer_ee),1,by=-1),])
colnames(bkgd_chull) <- c('xM','yM')
grid_M <- expand.grid(xM = seq(-sqrt(8), sqrt(8), length.out = 100),
                      yM = seq(-sqrt(2), sqrt(2), length.out = 50) )
ii_inpoly <- sp::point.in.polygon(grid_M[,1], grid_M[,2], bkgd_chull[,1], bkgd_chull[,2])

# - loess
# -- left
# fit_loess <- loess(ioa ~ xM * yM, data = ioa_Mo[[1]], degree = 2, span = 0.1,
#                    control = loess.control(surface = "direct"))
fit_loess <- loess(ioa ~ xM * yM, data =ioa_hex_left, degree = 2, span = 0.1,
                   control = loess.control(surface = "direct"))
pred_loess <- predict(fit_loess, grid_M, se = T)
df_pred <- grid_M
df_pred$Z <- melt(pred_loess$fit)$value
df_pred$equalSpace <- cut(df_pred$Z, breaks_ioa)
df_pred_left <- df_pred[ii_inpoly == 1, ]
df_pred_left$colSpace <- as.numeric(df_pred_left$equalSpace)
df_pred_left$colSpace[is.na(df_pred_left$colSpace)] <- 0
# df_pred_left$colSpace[df_pred_left$colSpace == 1] <- 0
# df_pred_left$colSpace <- df_pred_left$colSpace * (n_lvl+1)

xy_ashape <- ashape(ucl_rot_Mo2[ind_left_cone, 1:2] +
                      matrix(runif(sum(ind_left_cone)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
xy_edge <- xy_ashape$edges[,1:6]
xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
ahull_left <- xy_grid_ahull
ii_inpoly_left <- sp::point.in.polygon(grid_M[ii_inpoly == 1,1], grid_M[ii_inpoly == 1,2], xy_grid_ahull[,1], xy_grid_ahull[,2])
df_pred_left$Z[!ii_inpoly_left] <- 20 # outside set to large angles
df_pred_left$colSpace[!ii_inpoly_left] <- 0

# -- right
# fit_loess <- loess(ioa ~ xM * yM, data = ioa_Mo[[2]], degree = 2, span = 0.1,
#                    control = loess.control(surface = "direct"))
fit_loess <- loess(ioa ~ xM * yM, data =ioa_hex_right, degree = 2, span = 0.1,
                   control = loess.control(surface = "direct"))
pred_loess <- predict(fit_loess, grid_M, se = T)
df_pred <- grid_M
df_pred$Z <- melt(pred_loess$fit)$value
df_pred$equalSpace <- cut(df_pred$Z, breaks_ioa)
df_pred_right <- df_pred[ii_inpoly == 1, ]
df_pred_right$colSpace <- as.numeric(df_pred_right$equalSpace)
df_pred_right$colSpace[is.na(df_pred_right$colSpace)] <- 0

xy_ashape <- ashape(ucl_rot_Mo2[!ind_left_cone, 1:2] +
                      matrix(runif(sum(!ind_left_cone)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
xy_edge <- xy_ashape$edges[,1:6]
xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
ahull_right <- xy_grid_ahull
ii_inpoly_right <- sp::point.in.polygon(grid_M[ii_inpoly == 1,1], grid_M[ii_inpoly == 1,2], xy_grid_ahull[,1], xy_grid_ahull[,2])
df_pred_right$Z[!ii_inpoly_right] <- 20
df_pred_right$colSpace[!ii_inpoly_right] <- 0

# -- plot
# range(df_pred_right$Z)
# breaks_contour <- c(0, 2.5, 3, 4, 5)
breaks_contour <- c(0, 4, 5, 6)
# breaks_contour <- c(0, 4, 5, 5.5)
getPalette <- colorRampPalette(c('gray30', 'gray90'))
pal_contour <- getPalette(length(breaks_contour)+1)

# PLOT contour, choose left or right eye
df <- as.data.frame(ucl_rot_Mo_right)
plt <- plt_Mo +
  geom_contour(data=df_pred_right, aes(x=xM,y=yM, z=Z), breaks=breaks_contour,
               colour=pal_lr[2], alpha=0.9, lwd=0.5) +
  geom_point(data = as.data.frame(ucl_rot_Mo2[!ind_left_cone,]), aes(xM, yM), colour ='gray70', size=1, alpha=1) +
  # geom_contour(data=df_pred_left, aes(x=xM,y=yM, z=Z), breaks=breaks_contour,
  #              colour=pal_lr[2], alpha=0.9, lwd=0.5) +
  geom_point(data = df[lens_ixy[lens_ixy[,3] == clq,1], ], aes(xM, yM), colour = pal_axes[1], size=2) +
  geom_point(data = df[c(ind_Up_lens, ind_Down_lens), ], aes(xM, yM), colour = pal_axes[4], size=2) +
  geom_point(data = df[lens_ixy[lens_ixy[,2] == clp,1], ], aes(xM, yM), colour = pal_axes[2], size=2) +
  geom_point(data = cmer_right, aes(x=xM, y=yM), colour = pal_axes[3], size=2) +
  labs(title = paste("ioa", fn))
windows(width = 16, height = 8)
# pdf("eye_right_Mollweide_ioa_hex_4axes.pdf", width = 8.5, height = 4.5)
plt
# dev.off()


# # v-ioa along v-direction, -14:14 ---------------------------------------------------
# # there are 4 combinations one can do.Here are 2 of them.
# 
# ind_v <- matrix(ncol = 1, nrow = 59) 
# ioa_v <- matrix(ncol = 2, nrow = 59)
# iic <- vaxis_gen(0, ind_xy) 
# for (j in -28:30) { 
#   ii <- haxis_gen(j, ind_xy) 
#   if (length(ii[ii %in% iic]) > 0) {
#     ind_v[j+29] <- ii[ii %in% iic]  # the central meridian pt 
#     ioa <- rowMeans(nb_dist_ucl[ii, c(3,6)], na.rm = T)/pi*180 
#     ioa_v[j+29,] <- c(mean(ioa), sd(ioa) )
#   }
# }
# 
# # PLOT
# windows(width = 8, height = 6)
# # pdf("ioa_v_along_v.pdf")
# xx <- cart2sph2tp(ucl_rot_sm)[ind_v, 't'] %>% rev()
# yy <- ioa_v[,1] %>% rev()
# yysd <- ioa_v[,2] %>% rev()
# plot(xx, yy, xlab ="", pch=16, ylab='', main = "v-ioa vs v-direction",
#      ylim = c(3, 8.5), xlim= c(0, 180), xaxt='n',yaxt='n')
# # arrows(xx, yy-yysd, xx, yy+yysd, length=0.05, angle=90, code=3)
# axis(1, at = seq(180, 0, by = -45), labels = paste(seq(-90, 90, by =45), "°", sep = '') )
# axis(2, at = seq(3, 8.5, by =.5), labels = paste(seq(3, 8.5, by =.5), "°", sep = '') )
# # dev.off()

# # h-ioa along h-axis, -15:14 ---------------------------------------------------
# 
# ind_h <- matrix(ncol = 1, nrow = 30)
# ioa_h <- matrix(ncol = 2, nrow = 30)
# iic <- c(haxis_gen(1, ind_xy), haxis_gen(0, ind_xy))
# for (j in -13:16) {
#   ii <- vaxis_gen(j, ind_xy)
#   if (length(ii[ii %in% iic]) > 0) {
#     ind_h[j+14] <- ii[ii %in% iic]  # the central meridian pt
#     ioa <- rowMeans(nb_dist_ucl[ii,c(2,4,5,7)], na.rm = T)/pi*180
#     ioa_h[j+14,] <- c(mean(ioa), sd(ioa) )
#   }
# }
# 
# # PLOT
# windows(width = 8, height = 6)
# # pdf("ioa_h_along_h.pdf")
# xx <- cart2sph2tp(ucl_rot_sm)[ind_h, 'p']
# yy <- ioa_h[,1]
# plot(xx, yy, xlab ="", pch=16, ylab='', main = "h-ioa along h-direction",
#      ylim = c(4, 8), xlim= c(-10, 140), xaxt='n',yaxt='n')
# axis(1, at = seq(-20, 140, by = 20), labels = paste(seq(-20, 140, by =20), "°", sep = '') )
# axis(2, at = seq(4, 8, by =.5), labels = paste(seq(4, 8, by =.5), "°", sep = '') )
# # dev.off()

# Fig.3H, ioa_h, ioa_v, in Mollweide ---------------------------------------------------

ioa_h <- matrix(NA, ncol = 1, nrow = nrow(ucl_rot_sm))
for (j in 1:nrow(ucl_rot_sm)) {
  if (sum(is.na(nb_ind[j,])) == 0) {
    dd <- vector(mode='numeric', length=2L)
    dd[1] <- acos(ucl_rot_sm[nb_ind[j,2], ] %*% ucl_rot_sm[nb_ind[j,4],])
    dd[2] <- acos(ucl_rot_sm[nb_ind[j,5], ] %*% ucl_rot_sm[nb_ind[j,7],])
    ioa_h[nb_ind[j,1]] <- mean(dd)/pi*180
  }
}

ioa_v <- matrix(NA, ncol = 1, nrow = nrow(ucl_rot_sm))
for (j in 1:nrow(ucl_rot_sm)) {
  if (sum(is.na(nb_ind[j,])) == 0) {
    ioa_v[nb_ind[j,1]] <- mean(nb_dist_ucl[nb_ind[j,1], c(3,6)], na.rm = T)/pi*180
  }
}

quantile(na.omit(ioa_h), c(0, 0.05,0.1, 0.25,0.5,0.75, 0.9,0.95,1))
quantile(na.omit(ioa_v), c(0, 0.05,0.1, 0.25,0.5,0.75, 0.9,0.95,1))

# - Mollweide
df_pos <- data.frame(ucl_rot_Mo)

### ###
df_pos$quan<- ioa_h
bin_lim <- c(7, 12) #5-95%
### ###
df_pos$quan<- ioa_v
bin_lim <- c(3, 6)

df <- df_pos[, c('xM','yM','quan')]
df$quan[is.na(df$quan)] <- 20
colnames(df) <- c('x','y','z')
grid <- with(df, interp::interp(x, y, z))
# griddf <- subset(data.frame(x = rep(grid$x, nrow(grid$z)),
#                             y = rep(grid$y, each = ncol(grid$z)),
#                             z = as.numeric(grid$z)),
#                  !is.na(z))

plt <- plt_Mo34 + 
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1.1) +
  geom_point(data=df_pos, aes(x = xM, y = yM, colour = quan), size = 2) +
  scale_color_gradient(low = "blue", high= "gray85", limits=bin_lim, oob=scales::squish, trans='log',
                       breaks=bin_lim, labels=bin_lim, guide = guide_colorbar(title = "angle"), na.value = 'yellow') +
  # geom_contour(data=griddf, aes(x=x,y=y, z=z), breaks=breaks_contour,
  #              colour=pal_lr[2], alpha=0.9, lwd=0.5) +
  theme(legend.position = c(.9, .9) ) +
  labs(title = "ioa eye along h-axis")

windows(width = 12, height = 8)
# pdf("ioa_h.pdf", width = 8.5, height = 4.5)
# pdf("ioa_v.pdf", width = 8.5, height = 4.5)
plt
# dev.off()


# ED Fig.5G, axes in Mercator --------------------------------------------------------

# -- vaxis
plt <- plt_Mer 
for (j in seq(-13,16)) { 
  ind_axis <- vaxis_gen(j, ind_xy) 
  xy <- data.frame(ucl_rot_Merc)[ind_axis,]
  mer_Merc <- xy[order(xy$y),]
  plt <- plt + geom_path(data = mer_Merc, aes(x=x, y=y), colour = pal_axes[3], lwd=1)
}
windows(width = 5, height = 6.5)
plt
# ggsave(paste("v-axis_Mercator.pdf",sep=''), width = 5, height = 6.5)

# -- h-axis
plt <- plt_Mer # + geom_point(data=ucl_rot_Merc, aes(x = x, y = y), size = 2)
for (j in seq(-28,32)) {
  ind_axis <- haxis_gen(j, ind_xy)
  xy <- data.frame(ucl_rot_Merc)[ind_axis,]
  mer_Merc <- xy[order(xy$x),]
  plt <- plt + geom_path(data = mer_Merc, aes(x=x, y=y), colour = pal_axes[3], lwd=1)
}
windows(width = 5, height = 6.5)
plt
# ggsave(paste("h-axis_Mercator.pdf",sep=''), width = 5, height = 6.5)



# ED Fig.4B, normal VS tip-lens ------------------------------------------------------

nucl_eyemap <- nucl[eyemap[,2], ] 
lens_eyemap <- lens[eyemap[,2], ] 

# nopen3d()
# points3d(lens_eyemap)
# for (j in 1:nrow(eyemap)) {
#   arrow3d(lens_eyemap[j,], lens_eyemap[j,]+ucl_rot_sm[j,]*10, theta=pi/18, n=8, s=0.3,width=0.3, col = "gray30", type = "rotation")
#   arrow3d(lens_eyemap[j,], lens_eyemap[j,]+nucl_eyemap[j,]*10, theta=pi/18, n=8, s=0.3,width=0.3, col = "cyan", type = "rotation")
# }

# full right
nucl_rot_Mo_right <- nucl
colnames(nucl_rot_Mo_right) <- c('x','y','z')
nucl_rot_Mo_right %<>% as_tibble() %>%  
  mutate(y = -y) %>%
  mutate(theta = acos(z)) %>%
  mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
  mutate(t = theta / pi * 180, p = phi/pi*180) %>%  
  mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
  as.data.frame()
nucl_rot_Mo_right <- Mollweide(nucl_rot_Mo_right[,c('t', 'p')])
colnames(nucl_rot_Mo_right) <- c('xM','yM')


df_pos <- data.frame(ucl_rot_Mo_right)
df_pos$quan <- acos(rowSums(ucl_rot_right * nucl))

df_arrow <- data.frame(cbind(ucl_rot_Mo_right, nucl_rot_Mo_right) )
colnames(df_arrow) <- c('x','y','xend','yend')
df_arrow <- df_arrow[vaxis_gen(c(-10,1,10)),]

plt <- plt_Mo34 +
  geom_point(data=df_pos, aes(x = xM, y = yM), colour ='gray', size = 1) +
  geom_point(data= as.data.frame(nucl_rot_Mo_right),aes(x = xM, y = yM),colour='salmon', size = 1) +
  geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour="gray60",size =0.5) +
  labs(title = "gray ucl, pink normals") 
windows(width = 12, height = 8)
# pdf("ucl_vs_normal.pdf", width = 9.5, height = 6.5)
plt
# dev.off()


# Fig.3E, just 3d ucl ------------------------------------------------------------

nopen3d()
par3d('windowRect' = c(100,100,1000,1000))
spheres3d(0,0,0,0.99, col='grey95', alpha=1, lit=F)
points3d(ucl_rot_sm[ind_xy[ind_xy[,2] == clp, 1], ], size=12, col=pal_axes[2])
points3d(ucl_rot_sm[ind_xy[ind_xy[,3] == clq, 1], ], size=12, col=pal_axes[1])
points3d(ucl_rot_sm[vaxis_gen(clv, ind_xy), ], size=12, col=pal_axes[3])
points3d(ucl_rot_sm[haxis_gen(clh, ind_xy), ], size=12, col=pal_axes[4])

points3d(ucl_rot_sm, size=7, col='gray50')
arrow3d(c(1,0,0), c(1.6, 0, 0), theta = pi / 18, n = 8, col = "gray30", type = "rotation")
# text3d(c(1.3, -0.1, -0.3), texts = 'front', cex=3)
arrow3d(c(0,0,1), c(0, 0, 1.5), theta = pi / 12, n = 8, col = "gray30", type = "rotation")
# text3d(c(-0.2, -0.2, 1.5), texts = 'up', cex=3)
arrow3d(c(0,-1,0), c(0, -1.8, 0), theta = pi / 18, n = 8, col = "gray30", type = "rotation")
# text3d(c(0.2, -1.8, 0), texts = 'side', cex=3)
# planes3d(0,1,0, 0, alpha = 0.2, lit=F)
# planes3d(0,0,1, 0, alpha = 0.2, lit=F)
rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(-60/180*pi,1,0,0) %*% rotationMatrix(-35/180*pi,0,0,1))
# rgl.viewpoint(fov=0,zoom=0.7, userMatrix= rotationMatrix(-80/180*pi,1,0,0) %*% rotationMatrix(-47/180*pi,0,0,1))

# rgl.snapshot("ucl_rot.png")


# ED Fig.4H, ED Fig.4G, curvature and dia, NEED re-run Figure_0.R here --------

# - PLOT sphere
## ## all points
pts <- lens
xyz_data <- pts
## ## central points
ii <- nb_ind[nb_ind[1,],] %>% unique() #central
pts <- lens[as.integer(rownames(ucl_rot_sm[ii,])),]
xyz_data <- pts

colnames(xyz_data) <- c('x','y','z')
xyz_data2 <- xyz_data^2
Y <- rowSums(xyz_data2)
X <- cbind(2*xyz_data, rep(1,dim(xyz_data)[1]))
sph <- lsfit(X,Y,intercept = FALSE)
r_fit <- sqrt(sph[[1]][4]+sum(sph[[1]][1:3]^2)) #radius

nopen3d()
par3d('windowRect' = c(100,100,1000,1000))
spheres3d(sph[[1]][1],
          sph[[1]][2],
          sph[[1]][3],
          r_fit, col='grey80', alpha=0.8, lit=F)
spheres3d(lens[,1],lens[,2],lens[,3], 5, lit=F)
spheres3d(pts[,1],pts[,2],pts[,3], 6, col='red', lit=F)
# planes3d(0,1,0, 0, alpha = 0.2)
# planes3d(0,0,1, 0, alpha = 0.2)
view3d(fov=0,zoom=0.8, userMatrix= rotationMatrix(-90/180*pi,1,0,0) %*% rotationMatrix(-90/180*pi,0,0,1))
spheres3d(lens_left[,1],lens_left[,2],lens_left[,3], 5, lit=F)
# rgl.snapshot("curvature_2eye.png")

# - radius of curvature and lens dia
# load(file = paste0("../microCT/2023_eyemap/", "20240701", "_roc.RData"))
# load(file = paste0("../microCT/2023_eyemap/", "20240701", "_dia.RData"))
# load(paste0("data/microCT/20240701", "_roc.RData"))
# load(paste0("data/microCT/20240701", "_dia.RData"))

ind_roc <- ind_roc_right # [ind, p,v,q,h,sph]
ind_dia <- ind_dia_right
lens_ixy_hex <- ind_xy_right
unit_hex <- t(cbind(c(cos(30/180*pi), sin(30/180*pi)), c(-cos(30/180*pi), sin(30/180*pi))))
lens_ixy_hex[,2:3] <- lens_ixy_hex[,2:3] %*% unit_hex


# - PLOT cur
quantile(ind_roc[,c(3,5,6)], c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T)
bin_lim <- c(120, 320)
## ## choose a direction, 3-v, 5-h, 6-sph
df <- cbind(lens_ixy_hex[,2:3], ind_roc[,3]) %>% as.data.frame()
colnames(df) <- c('x','y','r')

plt <- ggplot() + 
  geom_point(data=df, aes(x = x, y = y, colour = r), size = 2) +
  scale_color_gradient(
    low = "gray95", high= "magenta", limits=bin_lim, oob=scales::squish,# trans='log',
    breaks=bin_lim, labels=bin_lim, guide = guide_colorbar(title = "r[um]"), 
    na.value = 'yellow') +
  # scale_x_continuous(limits = c(-sqrt(8), sqrt(8)), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(-sqrt(2), sqrt(2)), breaks = c(-1.5,0,1.5), labels = c(-1.5,0,1.5), expand = c(0, 0)) + # set +y as above eq
  theme_void() +
  theme(legend.position = c(.9, .9) ) +
  # theme(legend.position="none", panel.background = element_blank()) +
  coord_fixed(ratio = 1) +
  labs(title = "radius of curvature along v-axis")
  # labs(title = "radius of curvature along h-axis")
  # labs(title = "radius of curvature sph")
windows(width = 8, height = 8)
# pdf("roc_p.pdf")
# pdf("roc_v.pdf")
# pdf("roc_q.pdf")
# pdf("roc_h.pdf")
plt
# dev.off()


# - example h=const , central 5 lenses, re-run Figure_0 afterwards
load(paste0("data/microCT/20240701", ".RData"))
pts <- lens[!ind_left_lens,]
ls <- reghex(pts)
ind_nb <- ls[[1]]
ind_xy <- ls[[2]]

ii <- c(nb_ind[1,c(1,3,6)],nb_ind[nb_ind[1,3],3], nb_ind[nb_ind[1,6],6])
cc <- 0
N <- 2
ixy <- ind_xy[ind_xy[,3]-ind_xy[,2] == cc,,drop=F]
ixy <- ixy[order(ixy[,2]+ixy[,3]),] #order by v
xyz <- pts[ixy[,1],]
xyz5 <- pts[ixy[15 + (-N:N),1], ]

xyz_pca <- prcomp(xyz)
xyz_pc <- sweep(xyz, 2,xyz_pca$center) %*% xyz_pca$rotation
xyz5_pc <- sweep(xyz5, 2,xyz_pca$center) %*% xyz_pca$rotation

xyz_data <- xyz5_pc[, c(2,1)]
colnames(xyz_data) <- c('x','y')
xyz_data2 <- xyz_data^2
Y <- rowSums(xyz_data2)
X <- cbind(2*xyz_data, rep(1,dim(xyz_data)[1]))
cir <- lsfit(X,Y,intercept = FALSE)
r_fit <- sqrt(cir[[1]][3]+sum(cir[[1]][1:2]^2)) #radius cir[[1]] == cir[["coefficients"]]

# plot
dev.new()
# pdf('cir_fit_eg.pdf', width = 6, height = 6)
plot(xyz_pc[,2], xyz_pc[,1], pch=16, 
     asp = 1, xlim = c(-100, 650), ylim = c(-400,400))
points(xyz5_pc[,2], xyz5_pc[,1], pch=16, cex=1.5, col='red')
draw.circle(cir[[1]][1], cir[[1]][2], r_fit, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
# dev.off()

# - PLOT dia
quantile(ind_dia[,2], c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T)
bin_lim <- c(15.5, 17.5)
df <- cbind(lens_ixy_hex[,2:3], ind_dia[,2]) %>% as.data.frame()
colnames(df) <- c('x','y','r')

plt <- ggplot() + 
  geom_point(data=df, aes(x = x, y = y, colour = r), size = 2) +
  scale_color_gradient(
    low = "gray95", high= "magenta", limits=bin_lim, oob=scales::squish,# trans='log',
    breaks=bin_lim, labels=bin_lim, guide = guide_colorbar(title = "r[um]"), 
    na.value = 'yellow') +
  theme_void() +
  theme(legend.position = c(.9, .9) ) +
  # theme(legend.position="none", panel.background = element_blank()) +
  coord_fixed(ratio = 1) +
  labs(title = "diameter")
windows(width = 8, height = 8)
# pdf("dia_6.pdf")
plt
# dev.off()
