# T4 com etc --------------------------------------------------------------

LL <- 2
T4com <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
T4dir <- dir_type[[LL]][, c('rsxd','rsyd','rszd', 'rsx0','rsy0','rsz0')] %>% as.matrix() # reverse st optic flow dir
T4direye <- lens_type[[LL]][, c('xd','yd','zd', 'x0','y0','z0')] %>% as.matrix() # reverse st optic flow dir


# Fig.3J,  shear angle -------------------------------------------------------------
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

# - density
windows(width = 5, height = 4)
# pdf(paste("den_ang_eye.pdf", sep = ""), width = 5, height = 4)
dd <- density(nb_ang_ucl[,2], from= 30, to= 150, bw='SJ', na.rm = T)
plot(dd$x, dd$y, type='l', bty='n', col='black', lwd=3, xlim = c(30, 150), xaxt='n',yaxt='n',
     xlab ="", ylab='',main = '')
axis(1, at = seq(30, 150, by =30), labels = paste(seq(30, 150, by =30), "°", sep = ''), cex.axis = 1.5 )
axis(2, at = c(0, 1/512/diff(dd$x[1:2]) ), labels = c('0', '1'), cex.axis = 1.5 )
# dev.off()

# - 2D
df_pos <- data.frame(ucl_rot_Mo)
df_pos$quan <- nb_ang_ucl[,2]

# -- Mollweide
rg <- c(60, 90, 120) 
plt <- plt_Mo34 + 
  geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
  geom_point(data=df_pos, aes(x = xM, y = yM, colour = quan), size = 2) +
  scale_color_gradientn(colours = pal_heat1, values = scales::rescale(rg),
                        limits=range(rg), oob=scales::squish,breaks= rg,
                        labels= rg, guide = guide_colorbar(title = "shear") ) +
  theme(legend.position = c(.9, .9) ) +
  labs(title = "shear eye")
for (k in 1:3) {
  j <- match(ind_eghex_2[k], nb_ind[,1])
  pt <- ucl_rot_Mo[nb_ind[j,], ]
  pt <- pt[c(2,3,4,5,6,7,2), ]
  plt <- plt + geom_path(data=as.data.frame(pt), aes(x=xM, y=yM), colour='gray50',lwd=0.5)
}
windows(width = 9, height = 6)
# pdf("shear_eye.pdf", width = 8.5, height = 4.5)
plt
# dev.off()


# cont. ED Fig.5H, --------------------------------------------------------

# -- mercator
df <- ucl_rot_Merc
df$quan <- nb_ang_ucl[,2]

# rg <- c(45, 90, 135)
rg <- c(60, 90, 120)
plt <- plt_Mer +
  geom_path(data = cmer_Merc, aes(x=x, y=y), colour = pal_axes[3], lwd=1) +
  geom_point(data=df, aes(x=x, y=y, colour = quan), size = 1) +
  scale_color_gradientn(colours = pal_heat1, values = scales::rescale(rg), limits=range(rg), oob=scales::squish,
                        breaks= rg, labels= rg, guide = guide_colorbar(title = "shear") ) +
  theme(legend.position = c(.9, .9) ) +
  labs(title = "shear med")

windows(width = 8, height = 10)
# pdf("shear_eye_Merc.pdf",width = 8, height = 10)
plt
# dev.off()

# Fig.3J,  example shear angle, local meridian aligned----------------------------------------------------------

nb_coord <- rbind(c(-5/180*pi,0),
                  c(+5/180*pi,0),
                  c(0, -5/180*pi),
                  c(0, +5/180*pi) )

hex_eg <- list()
for (k in 1:3) {
  j <- match(ind_eghex_2[k], nb_ind[,1])
  
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
## ## choose a cell
j <- 3
pt <- hex_eg[[j]]
windows(width = 8, height = 8)
# pdf(file = paste("hex_eg_alpha_eye_", c("D","C","V")[j], ".pdf", sep=''))
plot(pt, pch=16, axes=F, ann=F, cex=2, asp=1)
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
