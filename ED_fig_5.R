# ED Fig.5A, ------------------------------------------------------------------

# regular hex grid
ind_axis <- vaxis_gen(clv, lens_ixy) 

# - in hex lattice
lens_ixy_hex <- lens_ixy
unit_hex <- t(cbind(c(cos(30/180*pi), sin(30/180*pi)), c(-cos(30/180*pi), sin(30/180*pi))))
lens_ixy_hex[,2:3] <- lens_ixy_hex[,2:3] %*% unit_hex

windows(width = 8, height = 8)
# pdf("axes_lens_hex.pdf")
plot(lens_ixy_hex[, 2:3], xaxt="n", yaxt="n", pch=16, cex=1.25, col='gray') 
points(lens_ixy_hex[lens_ixy[,2] == clq,  2:3], pch=16, col=pal_axes[2], cex=1.5)
points(lens_ixy_hex[lens_ixy[,3] == clp, 2:3], pch=16, col=pal_axes[1], cex=1.5)
points(lens_ixy_hex[lens_ixy[,1] %in% haxis_gen(clh), 2:3], pch=16, col=pal_axes[4], cex=1.5)
points(lens_ixy_hex[lens_ixy[,1] %in% ind_axis, 2:3], pch=16, col=pal_axes[3], cex=1.5)
title("4 axes in hex")
# dev.off()


# ED Fig.5F ---------------------------------------------------------------

# - on Mercator 
df <- as.data.frame(ucl_rot_Merc)
plt <- plt_Mer +
  geom_point(data= as.data.frame(ucl_rot_Merc_right), aes(x,y), colour = 'gray') +
  geom_point(data = df[ind_xy[ind_xy[,2] == clp,1], ], aes(x, y), colour = pal_axes[2], size=2) +
  geom_point(data = df[ind_xy[ind_xy[,3] == clq,1], ], aes(x, y), colour = pal_axes[1], size=2) +
  geom_point(data = df[haxis_gen(clh, ind_xy), ], aes(x, y), colour = pal_axes[4], size=2) +
  geom_point(data = df[vaxis_gen(clv, ind_xy), ], aes(x, y), colour = pal_axes[3], size=2) +
  labs(title = "Mercator 3 axes with eq")
windows(width = 8, height = 10)
# pdf("axes_Mercator.pdf",width = 8, height = 10)
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


# ED Fig.5E, --------------------------------------------------------------

windows(width = 4, height = 4)
# pdf(paste("asp_comp.pdf", sep = ""), width = 4, height = 4)
dt <- rowMeans(nb_dist_ucl[, c(2,4,5,7)]) / rowMeans(nb_dist_ucl[, c(3,6)]) #2023
dd <- density(dt, from= 0.4, to= 3.1, bw='SJ', na.rm = T)
plot(dd$x, dd$y, type='l', bty='n', col='black', lwd=3, xlim = c(0.5, 2), ylim=c(0,5), xaxt='n',yaxt='n',
     xlab ="", ylab='',main = '')
dt <- rowMeans(nb_dist_med[, c(2,4,5,7)]) / rowMeans(nb_dist_med[, c(3,6)]) #2023
dd <- density(dt, from= 0.4, to= 3.1, bw='SJ', na.rm = T)
points(dd$x, dd$y, type='l', bty='n', col='blue', lwd=3)
axis(1, at = seq(0.5, 2, by =0.5), labels = seq(0.5, 2, by =0.5), cex.axis = 1.5 )
axis(2, at = c(0, 1/512/diff(dd$x[1:2])*((3.1-0.4)/(2-0.5)) ), labels = c('0', '1'), cex.axis = 1.5 )
# dev.off()


# ED Fig.5C, ED Fig.5D,  ED Fig.5J, asp ratio------------------------------------


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
  theme(legend.position = c(.9, .9) ) +
  labs(title = "aspect ratio eye")

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

windows(width = 12, height = 8)
# pdf("asp_med.pdf", width = 8.5, height = 4.5)
plt
# dev.off()


# ED Fig.5C, example asp ratio, local meridian aligned----------------------------------------------------------

nb_coord <- rbind(c(-5/180*pi,0),
                  c(+5/180*pi,0),
                  c(0, -5/180*pi),
                  c(0, +5/180*pi) )

hex_eg <- list()
for (k in 1:3) {
  j <- match(ind_eghex[k], nb_ind[,1]) # 2023
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