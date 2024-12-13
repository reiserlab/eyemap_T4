# change hex ordering 20240416
# "Angle in med ref"
# "RF in common med ref"

# Mi1 equator, use PR num -------------------------------------------
ind_up_Mi1 <- c(599,  39, 207, 606, 388,  10,  12, 157, 361, 287, 503, 105)
ind_down_Mi1 <- c(401, 334, 209, 210, 720, 667, 590, 628, 489, 289, 393)

# ME layer ----------------------------------------------------------------

# - coord transform via pca
tar <- T4_eg[[1]]
ind_D = match(tar$tags$"dendrite start" , tar$d$PointNo)
node_dend <- child_node(tar, tar$d[ind_D,])
cen <- colMeans(xyzmatrix(node_dend)) 
cen_ii <- rowSums(sweep(Mi1_M10_xyz, 2, cen)^2) %>% order(decreasing=F) %>% head(1+6+12+18)

node_xyz <- Mi1_M10_xyz[cen_ii,]
me_pca <- prcomp(node_xyz)
if (me_pca$rotation[,3] %*% c(326665, 259455, 235136) > 0) {
  me_pca$rotation[,3] <- - me_pca$rotation[,3]
}
if (me_pca$rotation[,1] %*% c(-19607, -44222, 12645) < 0) {
  me_pca$rotation[,1] <- - me_pca$rotation[,1]
}
if (t(cross3D(me_pca$rotation[,3],me_pca$rotation[,1])) %*% me_pca$rotation[,2] < 0 ) {
  me_pca$rotation[,2] <- - me_pca$rotation[,2]
}


# # DEBUG
# nopen3d()
# points3d(Mi1_M10_xyz, size = 7)
# points3d(Mi1_M10_xyz[pt_100_Mi1,,drop=F], size = 15, col ='green')
# points3d(Mi1_M10_xyz[ind_up_Mi1, ], size = 10, col ='gold2')
# points3d(Mi1_M10_xyz[ind_down_Mi1, ], size = 10, col ='magenta')
# 
# plot3d(Mi1[[pt_100_Mi1]], lwd=2)
# plot3d(T4_eg, soma = T)
# 
# points3d(Mi1_M10_xyz[cen_ii,], col='cyan', size = 17)
# 
# arrow3d(me_pca$center, me_pca$center+me_pca$rotation[,1]*5e4, col ='blue',type = "rotation")
# arrow3d(me_pca$center, me_pca$center+me_pca$rotation[,2]*5e4, col ='red',type = "rotation")
# arrow3d(me_pca$center, me_pca$center+me_pca$rotation[,3]*5e4, col ='pink',type = "rotation")


# - xform
xyz <- t(ME_msh$vb[1:3,])
xyz <- sweep(xyz, 2, me_pca$center) %*% me_pca$rotation
ME_msh_xform <- ashape3d(xyz, alpha = 50000) %>% as.mesh3d()

xyz <- t(LOP_msh_mod$vb[1:3,])
xyz <- sweep(xyz, 2, me_pca$center) %*% me_pca$rotation
LOP_msh_xform <- ashape3d(xyz, alpha = 50000) %>% as.mesh3d()

T4_eg_xform <- xEucl_neu(T4_eg, me_pca$rotation, me_pca$center)
Mi1_eg_xform <- xEucl_neu(Mi1[Mi1_eg_ind], me_pca$rotation, me_pca$center)
PR_xform <- xEucl_neu(PR, me_pca$rotation, me_pca$center)
L1_xform <- xEucl_neu(L1, me_pca$rotation, me_pca$center)
R7_xform <- xEucl_neu(R7, me_pca$rotation, me_pca$center)
TmY5a_xform <- xEucl_neu(TmY5a, me_pca$rotation, me_pca$center)
# Mi1_eg_xform <- xEucl_neu(Mi1[c(334,583,10)], me_pca$rotation, me_pca$center)
# Mi1_eg_xform <- xEucl_neu(Mi1[c(334,209, 401)], me_pca$rotation, me_pca$center)
Mi1_M10_xyz_xform <- sweep(Mi1_M10_xyz, 2, me_pca$center) %*% me_pca$rotation
Mi1_M5_xyz_xform <- sweep(Mi1_M5_xyz, 2, me_pca$center) %*% me_pca$rotation

# - make new local mesh 
xyz <- t(ME_msh_xform$vb[1:3,])
dd <- sqrt(rowSums(xyz[,1:2]^2))
xyz_msh <- xyz[dd < 0.75*max(dist(node_xyz)), ]
ME_msh_local <- ashape3d(xyz_msh, alpha = 60000) %>% as.mesh3d()

# # DEBUG
# nopen3d()
# # points3d(Mi1_M10_xyz_xform, size = 7)
# points3d(Mi1_M10_xyz_xform[pt_100_Mi1,,drop=F], size = 15, col ='green')
# # points3d(Mi1_M10_xyz_xform[ind_up_Mi1, ], size = 10, col ='gold2')
# # points3d(Mi1_M10_xyz_xform[ind_down_Mi1, ], size = 10, col ='magenta')
# points3d(Mi1_M10_xyz_xform[cen_ii,], col='cyan', size = 17)
# plot3d(Mi1_eg_xform, lwd=2, soma=T)
# plot3d(T4_eg_xform, soma = T)
# shade3d(ME_msh_local, alpha=0.1)
# rgl.viewpoint(fov=0,zoom=1, userMatrix= rotationMatrix(-90/180*pi,1,0,0) %*%
#                 rotationMatrix(90/180*pi,0,0,1) )
# # segments3d(layers_me, lwd=1)

# - ME boundary
# fit a surface to the mesh boundary and them sample the fit in the desired positions
xy <- data.frame(x= seq(-10000, 10000, length.out = 20), y = seq(10000, -10000, length.out = 20))

# top
xyz_top <- xyz_msh[xyz_msh[,3] > 0, ]
x <- xyz_top[,1]; y <- xyz_top[,2]; z <- xyz_top[,3]
fitlm <- lm(z ~ poly(x, y, degree = 2, raw = T))
valfit <- predict(fitlm, xy) #generate values from the fit
top_bd <- cbind(xy, valfit)
# bottom
xyz_bot <- xyz_msh[xyz_msh[,3] < 0, ]
x <- xyz_bot[,1]; y <- xyz_bot[,2]; z <- xyz_bot[,3]
fitlm <- lm(z ~ poly(x, y, degree = 2, raw = T))
valfit <- predict(fitlm, xy) #generate values from the fit
bot_bd <- cbind(xy, valfit)

# # 2D mesh as outline
# ME_msh_local_bd <- rbind(top_bd, apply(bot_bd, MARGIN = 2, rev)) 
# ME_msh_local_bd <- rbind(ME_msh_local_bd[-1,], ME_msh_local_bd[1,]) %>%
#   cbind(ME_msh_local_bd, .) %>%
#   t() %>%
#   matrix(., ncol = 3, byrow = T)


# --layer definition from 7-column paper
# https://www.pnas.org/content/112/44/13711?ijkey=682c42b2318fb1efb9fe6d7ad6da91752ce3fa22&keytype2=tf_ipsecsha
# bp7c <- c(0, 5.5,14.5,20,23,28,36,42,48,56,61) #layer 1 to 10
bp7c <- c(0, 5.5, 14.5, 20, 23, 28, 33, 39, 45, 56, 61) #layer 1 to 10, modify 6-7-8-9 by Mi1, C2, Tm5 and Tm20
bp7c <- max(bp7c) - bp7c %>% rev() #layer 10 to 1
bp7c_prob <- bp7c/max(bp7c)
z_qua <- apply(cbind(top_bd$valfit, bot_bd$valfit), MARGIN = 1,
               function(x) quantile(x, probs = bp7c_prob) )
layers <- matrix(ncol = 3, nrow = 0)
for (j in 1:length(bp7c_prob)) {
  cc <- cbind(xy, z_qua[j,])
  c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
    t() %>%
    matrix(., ncol = 3, byrow = T)
  layers <- rbind(layers, c2)
}
layers[,3] <- layers[,3] * 0.97 + 1100
layers_me <- layers


# Fig.2B, plot Mi1 eg with layers -------------------------------------------------

# - with neuropil
zAng <- 80
zAngRot <- zAng/180*pi
zRot <- matrix(c(cos(zAngRot), sin(zAngRot), 0,
                 -sin(zAngRot), cos(zAngRot), 0,
                 0, 0, 1), ncol = 3, byrow = T)
sbar <- matrix(c(0,0,0, 1e4, 0,0), ncol=3,byrow=T) %*% zRot

nopen3d()
par3d('windowRect' = c(100,100,1300,1300))
plot3d(PR_xform, col='brown', alpha=0.5, lwd=3, soma=T)
plot3d(L1_xform, col='blue', alpha=0.5, lwd=3, soma=T)
# plot3d(R7_xform, col='blue', alpha=0.5, lwd=3, soma=T)
# plot3d(TmY5a_xform[c(13)], col='green', alpha=0.8, lwd=3, soma=T)
plot3d(Mi1_eg_xform[[1]], col='gray30',lwd=3, soma=T, alpha=1)
plot3d(T4_eg_xform[[2]], col=pal_T4[[2]], lwd=3, soma=T)
points3d(sweep(med_xyz,2,me_pca$center) %*% me_pca$rotation, size=5, col='magenta')
shade3d(ME_msh_xform, alpha=0.2, col='gray')
shade3d(LOP_msh_xform, alpha=0.2, col='gray')
segments3d(sweep(sbar,2,c(0,-50000,-70000),'+'), lwd=2)
rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(-90/180*pi,1,0,0) %*%
                rotationMatrix(-80/180*pi,0,0,1) )
# rgl.snapshot(filename = paste("Mi1_T4_mesh.png", sep = ''))

rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(-90/180*pi,1,0,0) %*%
                rotationMatrix(-(80+90)/180*pi,0,0,1) )

# rgl.snapshot(filename = paste("Mi1_T4_mesh_sideview.png", sep = ''))
# rgl.snapshot(filename = paste("L1_Mi1_T4_mesh_sideview.png", sep = ''))

# - row of Mi1
zAng <- 135
zAngRot <- zAng/180*pi
zRot <- matrix(c(cos(zAngRot), sin(zAngRot), 0,
                 -sin(zAngRot), cos(zAngRot), 0,
                 0, 0, 1), ncol = 3, byrow = T)
sbar <- matrix(c(0,0,0, 1e4, 0,0), ncol=3,byrow=T) %*% zRot

nopen3d()
par3d('windowRect' = c(100,100,1300,1300))
plot3d(Mi1_eg_xform[[1]], col='black',lwd=2, soma=T, lit=F)
plot3d(Mi1_eg_xform[c(2,5,6)], col=c('gray60','gray80','gray40'),lwd=2, soma=T, lit=F)

segments3d(sweep(layers_me,2,c(0,0,-1000)), lwd=1)
# pch3d(Mi1_M10_xyz_xform[Mi1_eg_ind[1],,drop=F], pch=16, radius=7e3, alpha=0.5, col='magenta')
# pch3d(Mi1_M5_xyz_xform[Mi1_eg_ind[1],,drop=F], pch=16, radius=7e3, alpha=0.5, col='magenta')
# plot3d(T4_eg_xform, col=pal_T4, lwd=2, soma=T)
# segments3d(layers_LOP_mepca, lwd=1)
segments3d(sweep(sbar,2,c(-10000,0000,-10000),'+'), lwd=2)
# segments3d(sweep(sbar,2,c(0,-10000,-10000),'+'), lwd=2)
rgl.viewpoint(fov=0,zoom=1, userMatrix= rotationMatrix(-90/180*pi,1,0,0) %*%
                rotationMatrix(-135/180*pi,0,0,1) )

# rgl.postscript("Mi1_layer_med.pdf", fmt = 'pdf')
# rgl.snapshot(filename = paste("Mi1_eg_one.png", sep = ''))
# rgl.snapshot(filename = paste("Mi1_eg.png", sep = ''))
# rgl.snapshot(filename = paste("Mi1_T4_eg.png", sep = ''))


# Fig.2C, Mi1 med col with annotation (and T4b) -------------------------------------------------------------

LL <- 2
vcom <- as.matrix(dir_type[[LL]][, c("comx","comy","comz")])
ind_eg2 <- c(142, 163,  139, 169, 121, 108,  65, 143,  170)
ind_eg2_3 <- ind_eg2[c(1,4,7)]
  
# - M10
ps = 13 #point size
nopen3d()
par3d('windowRect' = c(100,100,1300,1300))

points3d(Mi1_M10_xyz[c(Mi1_ind_PR[[8]], c(35,8,755,334,10,12,590,361,89)),], size = ps, col='tan') # =7
# points3d(Mi1_M10_xyz[Mi1_ind_PR[[8]],], size = 10, col='tan')
# points3d(Mi1_M10_xyz[c( c(35,8,755,334,10,12,590,361,89)),], size = 10, col='blue')
points3d(Mi1_M10_xyz[Mi1_ind_PR[[7]],], size = ps, col='gold') # =8
# points3d(Mi1_M10_xyz[unlist(Mi1_ind_PR[c(5,6)]),], size = 12, col='turquoise')
points3d(Mi1_M10_xyz[Mi1_neu_ind[match(anno_Mi1_DRA$skid, anno_Mi1$skid)],], size = ps, col='#CE7FFF')
# points3d(Mi1_M10_xyz[Mi1_ind_noPR,], size = 10, col='#CE7FFF')
pch3d(Mi1_M10_xyz[Mi1_ind_noPR,], pch = 1, col='black', cex= 0.11)

points3d(Mi1_M10_xyz[med_ixy[med_ixy[,3] == 0, 1],], size= ps, col=pal_axes[1])
points3d(Mi1_M10_xyz[med_ixy[med_ixy[,2] == -1, 1][-c(14,27)],], size= ps, col=pal_axes[2])
points3d(Mi1_M10_xyz[med_ixy[med_ixy[,2] == med_ixy[,3]-1,1][1:29],], size= ps, col=pal_axes[3])
points3d(Mi1_M10_xyz[med_ixy[med_ixy[,2] == -med_ixy[,3]-1,1],], size=ps, col=pal_axes[4])


ii <- c(na.omit(unlist(Mi1_ind_PR[7:8])), c(35,8,755,334,10,12,590,361,89),
        Mi1_neu_ind[match(anno_Mi1_DRA$skid, anno_Mi1$skid)], Mi1_ind_noPR,
        med_ixy[med_ixy[,3] == 0, 1],med_ixy[med_ixy[,2] == -1, 1],med_ixy[med_ixy[,2] == med_ixy[,3]-1,1] )

ii <- c(med_ixy[med_ixy[,3] == 0, 1],
        med_ixy[med_ixy[,2] == -1, 1],
        med_ixy[med_ixy[,2] == med_ixy[,3]-1,1],
        med_ixy[med_ixy[,2] == -med_ixy[,3]-1,1],
        Mi1_ind_noPR)

points3d(Mi1_M10_xyz[-ii,], size = ps, col='gray')

# # points3d(Mi1_M10_xyz[unlist(Mi1_ind_PR[c(4,3,2)]),], size = 20, col='salmon')
# # text3d(Mi1_M10_xyz[1:length(Mi1_neu_ind), ], texts = Mi1_neu_ind, adj=-0.4)
# pch3d(Mi1_M10_xyz[c(xb_ls[[1]][[1]],xb_ls[[2]][[1]]),], pch=3, cex=0.5, lwd=2, col='pink')
# yb0 <- sapply(xb_ls[c(1,3)], function(x) sapply(x, function(x){x[1]})) %>% unlist()
# yb0 <- yb0[!(yb0 %in% c(773, 767, 763))]
# pch3d(Mi1_M10_xyz[yb0,], pch=3, cex=0.5, lwd=2, col='pink')

# pch3d(vcom[-ind_eg2_3,], pch=1, cex=0.5, lwd=1, color=pal_T4[LL])
# pch3d(vcom[ind_eg2_3,], pch=1, cex=1, lwd=3, color='black')

# rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(180/180*pi,1,0,0) %*%
#                 rotationMatrix(-60/180*pi,0,1,0) %*%
#                 rotationMatrix(10/180*pi,1,0,0))

rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(180/180*pi,1,0,0) %*%
                rotationMatrix(120/180*pi,0,1,0) %*%
                rotationMatrix(20/180*pi,1,0,0))
# scale bar
c(10000,0,0) %>% 
  rotate3d(180/180*pi,1,0,0) %>%
  rotate3d(120/180*pi,0,1,0) %>%
  rotate3d(20/180*pi,1,0,0) %>%
  rbind(c(0,0,0), .) %>%
  sweep(2, c(3e5, 3.5e5, 2.5e5), '+') %>%
  segments3d()


# rgl.snapshot(filename = paste("Mi1_", "equ.png", sep = ''))
# rgl.postscript(filename = "Mi1_equ.svg", fmt = 'svg', drawText = T)
# rgl.snapshot(filename = paste("Mi1_", "4axes.png", sep = ''))


# Fig.2D, 3 axes in hex or squre grid --------------------------------------------------

# - in hex lattice
med_ixy_hex <- med_ixy
unit_hex <- t(cbind(c(cos(30/180*pi), sin(30/180*pi)), c(-cos(30/180*pi), sin(30/180*pi))))
med_ixy_hex[,2:3] <- med_ixy_hex[,2:3] %*% unit_hex

windows(width = 8, height = 8)
# pdf("axes_med_hex.pdf")
plot(med_ixy_hex[, 2:3], xaxt="n", yaxt="n", pch=16, cex=1.25, col='gray')
points(med_ixy_hex[med_ixy[,3] == 0, 2:3], pch=16, col=pal_axes[1], cex=1.5)
points(med_ixy_hex[med_ixy[,2] == 0,  2:3], pch=16, col=pal_axes[2], cex=1.5)
points(med_ixy_hex[med_ixy[,1] %in% eyemap[match(haxis_gen(0), eyemap[,2]),1], 2:3], pch=16, col=pal_axes[4], cex=1.5)
points(med_ixy_hex[med_ixy[,1] %in% eyemap[match(vaxis_gen(0), eyemap[,2]),1], 2:3], pch=16, col=pal_axes[3], cex=1.5)
# text(med_ixy[,2:3], labels = med_ixy[,1], cex = 0.7, adj = -0.5)
title("4 axes in hex")
# dev.off()


# - in square lattice
med_ixy_square <- med_ixy_hex
med_ixy_square[,2] <- med_ixy_square[,2] / sqrt(3)

windows(width = 5, height = 8)
# pdf("axes_med_square.pdf", width = 5, height = 8)
plot(med_ixy_square[, 2:3], asp=1, xaxt="n", yaxt="n",ann=F, pch=16, cex=1.1, col='gray')
points(med_ixy_square[med_ixy[,3] == 0, 2:3], pch=16, col=pal_axes[1], cex=1.1)
points(med_ixy_square[med_ixy[,2] == 0,  2:3], pch=16, col=pal_axes[2], cex=1.1)
points(med_ixy_square[med_ixy[,1] %in% eyemap[match(haxis_gen(0), eyemap[,2]),1], 2:3], pch=16, col=pal_axes[4], cex=1.1)
points(med_ixy_square[med_ixy[,1] %in% eyemap[match(vaxis_gen(0), eyemap[,2]),1], 2:3], pch=16, col=pal_axes[3], cex=1.1)
# text(med_ixy[,2:3], labels = med_ixy[,1], cex = 0.7, adj = -0.5)
# title("4 axes in square")
# dev.off()

# Fig.2E, all T4b in med and lop -----------------------------------------------------------------

# T4b_me <- nlapply(T4b, function(x) subset(x, pointsinside(x,ME_msh,rval='distance') > -5000))
# T4b_lop <- nlapply(T4b, function(x) subset(x, pointsinside(x,LOP_msh,rval='distance') > 0))

T4b_me <- T4b
for (j in 1:length(T4b)) {
  tar <- T4b[[j]]
  ind_T = match(tar$tags$"dendrite start" , tar$d$PointNo)
  targ <- as.ngraph(tar)
  ii_root <- ind_T
  # subtree and resample
  sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
  T4b_me[[j]] <- subset(tar, sub_points) 
}

# ii <- c(139,108,170) # cf ind_eg2
ii <- c(139,170) # cf ind_eg2

nopen3d()
par3d('windowRect' = c(100,100,1300,1300))
plot3d(T4b_me[ii], lwd=4, col="black")
# plot3d(T4b_me[-ii], lwd=2, col=c("gray70","gray75","gray80","gray85"))
plot3d(T4b_me[-ii], lwd=2, 
       col=rep(c("#C69F9F","#B3C69F","#9FC6C6","#B29FC6","gray"),35)[1:length(T4b_me[-ii])]
       )

# shade3d(ME_msh, col='gray', alpha=0.1)
rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(180/180*pi,1,0,0) %*%
                rotationMatrix(120/180*pi,0,1,0) %*%
                rotationMatrix(20/180*pi,1,0,0))
# scale bar
c(10000,0,0) %>% 
  rotate3d(180/180*pi,1,0,0) %>%
  rotate3d(120/180*pi,0,1,0) %>%
  rotate3d(20/180*pi,1,0,0) %>%
  rbind(c(0,0,0), .) %>%
  sweep(2, c(3e5, 3.5e5, 2.5e5), '+') %>%
  segments3d(lwd=4)
# rgl.snapshot(filename = "T4b_me.png")

# nopen3d()
# plot3d(T4b_lop, lwd=2)
# shade3d(LOP_msh_mod, col='gray', alpha=0.2)
# # rgl.snapshot(filename = "T4b_lop.png")

# ED Fig.2G, all T4d in med  -----------------------------------------------------------------

T4d <- T4_dend[[4]]
for (j in 1:length(T4d)) {
  tar <- T4d[[j]]
  ind_T = match(tar$tags$"dendrite start" , tar$d$PointNo)
  targ <- as.ngraph(tar)
  ii_root <- ind_T
  # subtree and resample
  sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
  T4d[[j]] <- subset(tar, sub_points) 
}

ii <- c(109, 71, 52) # cf ind_eg2
ii <- c(76, 69) # cf ind_eg2

nopen3d()
par3d('windowRect' = c(100,100,1300,1300))
plot3d(T4d[ii], lwd=4, col="black")
# plot3d(T4d[-ii], lwd=2, col=c("gray70","gray75","gray80","gray85"))
plot3d(T4d[-ii], lwd=2, 
       col=rep(c("#C69F9F","#B3C69F","#9FC6C6","#B29FC6","gray"),35)[1:length(T4d[-ii])]
)

# shade3d(ME_msh, col='gray', alpha=0.1)
rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(180/180*pi,1,0,0) %*%
                rotationMatrix(120/180*pi,0,1,0) %*%
                rotationMatrix(20/180*pi,1,0,0))
# scale bar
c(10000,0,0) %>% 
  rotate3d(180/180*pi,1,0,0) %>%
  rotate3d(120/180*pi,0,1,0) %>%
  rotate3d(20/180*pi,1,0,0) %>%
  rbind(c(0,0,0), .) %>%
  sweep(2, c(3e5, 3.5e5, 2.5e5), '+') %>%
  segments3d(lwd=4)
# rgl.snapshot(filename = "T4d_me.png")

# ED Fig.2E, chiasm vs v-axis --------------------------------------------------------

# ind_chi <- c(280, 279, 378, 404, 684, 683)
ind_chi <- match(na.omit(match(anno_chi$skid, anno_Mi1$skid)), rownames(med_xyz))
 
# PLOT R7
nopen3d()
par3d('windowRect' = c(100,100,1300,1300))
plot3d(chiR_nbgp3[c(1:3,8)], lwd=3, col= brewer.pal(6, "Blues")[2:5])
# plot3d(chiMi1[c(10)], col="gray30", lwd=2, WithNodes = F, soma = T)
# plot3d(chiMi1[c(1)], col="brown", lwd=2, WithNodes = F, soma = T)
rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(180/180*pi,1,0,0) %*%
                rotationMatrix(120/180*pi,0,1,0) %*%
                rotationMatrix(20/180*pi,1,0,0))
# rgl.snapshot(filename = paste("chiasm_R7.png", sep = ''))

# PLOT
nopen3d()
par3d('windowRect' = c(100,100,1300,1300))
# ii <- as.integer(rownames(med_xyz)[match(vaxis_gen(0),eyemap[,2])])
# ii <- which(med_ixy[,1] %in% eyemap[match(vaxis_gen(0), eyemap[,2]),1])
ii <- na.omit(eyemap[match(vaxis_gen(0), eyemap[,2]),1])
points3d(Mi1_M10_xyz[ii,], size=10, col=pal_axes[3])
points3d(Mi1_M10_xyz[-ii,], size=10, col='gray')

# points3d(med_xyz[match(vaxis_gen(0),eyemap[,2]),], size=10, col=pal_axes[3])
# points3d(med_xyz[-match(vaxis_gen(0),eyemap[,2]),], col = 'gray', size = 10)
plot3d(chiMi1[c(2,1,6)], col="gray30", lwd=2, WithNodes = F, soma = T)
plot3d(chiMi1[c(3,10,7)], col="brown", lwd=2, WithNodes = F, soma = T)
points3d(med_xyz[ind_chi,], col = 'green', size = 17)
plot3d(chiR, WithNodes = F, col='gray', lwd=3)

rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(180/180*pi,1,0,0) %*%
                rotationMatrix(120/180*pi,0,1,0) %*%
                rotationMatrix(20/180*pi,1,0,0))
# scale bar
c(10000,0,0) %>% 
  rotate3d(180/180*pi,1,0,0) %>%
  rotate3d(120/180*pi,0,1,0) %>%
  rotate3d(20/180*pi,1,0,0) %>%
  rbind(c(0,0,0), .) %>%
  sweep(2, c(3e5, 3.5e5, 2.5e5), '+') %>%
  segments3d()

# rgl.snapshot(filename = paste("chiasm_vaxis.png", sep = ''))

# # more T4b eg with Strahler ----------------------------------------------------------------
# # choose nb via nb_ind
# # align to meridian by mapping ref points to med via np
# 
# nb_coord <- rbind(c(0,0),
#                   c(-5/180*pi,0),
#                   c(+5/180*pi,0),
#                   c(0, -5/180*pi),
#                   c(0, +5/180*pi) )
# 
# LL <- 2 # choose type
# # eye
# com_xyz_eye <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
# # med
# com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
# v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix()
# v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
# v3 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
# v4 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()
# # nopen3d()
# # points3d(ucl_rot_sm)
# # points3d(ucl_rot_sm[match(c(xb_ls[[1]][[1]],xb_ls[[2]][[1]]), eyemap[,1]),], col=pal_axes[1],size =10,alpha=0.5)
# # points3d(ucl_rot_sm[match(ybm[[1]], eyemap[,1]),], col=pal_axes[2],size =10,alpha=0.5)
# # points3d(ucl_rot_sm[match(vaxis_gen(-1), eyemap[,2]),], col=pal_axes[3],size =10,alpha=0.5)
# # points3d(ucl_rot_sm[c(ind_Up_ucl, ind_Down_ucl),], col=pal_axes[4],size =10, alpha=0.5)
# # points3d(com_xyz_eye, size=10, col='blue',alpha=0.9)
# # identify3d(com_xyz_eye, adj=-1)
# 
# # # 1 center + 4 equ + 4 meridian
# # ind_eg2 <- c(142, 163,  25, 139,  23,  80, 107, 133, 112) 
# # 3 + 3 + 3
# # ind_eg2 <- c(142, 163,  22, 169, 121, 156,  65, 143,  12)
# # # ind_eg2 <- c(2,22,150) # med_xyz =c(510, 662, 366)
# # ind_eg2 <- c(108, 139,170) # med_xyz = , ucl=c(500,554,273)
# ind_eg2 <- c(142, 163,  139, 169, 121, 108,  65, 143,  170)
# 
# T4b_eg <- T4_dend[[LL]][ind_eg2]
# 
# # nopen3d()
# # points3d(med_xyz)
# # text3d(com_xyz[ind_eg2,], texts = ind_eg2, adj = 1,col='red')
# 
# for (j in 1:length(T4b_eg)) {
#   # 1+6+12 hex nb
#   ii <- sweep(ucl_rot_sm, 2, com_xyz_eye[ind_eg2[j], ] )^2 %>% rowSums() %>% which.min()
#   ii2 <- c(nb_ind[nb_ind[ii,],]) %>% unique()
#   # points3d(ucl_rot_sm[ii2,], size=10, col='black')
#   
#   # ref nb
#   pt <- ucl_rot_sm[ii, ] %>% matrix(ncol=3)
#   tp <- cart2sphZ(pt)[,2:3]
#   rtp <- cbind(1, sweep(nb_coord,2,tp,'+'))
#   xyz_ref <- sph2cartZ(rtp)
#   
#   # np to med
#   xyz_ref_med <- matrix(nrow = nrow(xyz_ref), ncol = 3) 
#   xyz_eval <- data.frame(mc.x = xyz_ref[,1], mc.y = xyz_ref[,2], mc.z = xyz_ref[,3])
#   for (k in 1:3) {
#     npdata <- data.frame(mc = ucl_rot_sm[ii2,], ec = med_xyz[ii2,k])
#     bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
#     model_np <- npreg(bw)
#     xyz_ref_med[,k] <- predict(model_np, newdata = xyz_eval)
#   }
#   
#   # pc
#   zz <- cross3D((xyz_ref_med[2,]-xyz_ref_med[1,]), (xyz_ref_med[5,]-xyz_ref_med[1,])) #pointing inwards
#   pc <- prcomp(xyz_ref_med)
#   if (pc$rotation[,3] %*% zz < 0) {
#     pc$rotation <- -pc$rotation
#   }
#   if (c(cross3D(pc$rotation[,3], pc$rotation[,1])) %*% pc$rotation[,2] < 0) {
#     pc$rotation[,2] <- -pc$rotation[,2]
#   }
#   
#   xyz_ref_med_pc <- sweep(xyz_ref_med, 2, xyz_ref_med[1,]) %*% pc$rotation
#   y_ref <- diff(xyz_ref_med_pc[c(3,2),1:2])
#   ang <- acos(y_ref %*% c(0,1) / sqrt(sum(y_ref^2)))
#   ang <- if_else(y_ref[1] < 0, -ang, ang)
#   rotM <- matrix(c(cos(ang), sin(ang),0,
#                    -sin(ang), cos(ang), 0,
#                    0, 0, 1), ncol = 3, byrow = T)
#   xyz_ref_med_pc_rot <- xyz_ref_med_pc %*% rotM
#   
#   # T4
#   tar <- T4b_eg[[j]]
#   ind_D = match(tar$tags$`dendrite start` , tar$d$PointNo)
#   df_D <-  tar$d[ind_D,]
#   root_xyz <- tar$d[ind_D, c("X","Y","Z")]
#   
#   # - find the subtree with root = dendrite start
#   targ <- as.ngraph(tar)
#   ii_root <- ind_D
#   # subtree and Strahler order
#   sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
#   subtree <- subset(tar, sub_points) 
#   subtree_g <- as.ngraph(subtree, weights = T)
#   subtree_so <- strahler_order(subtree) # Strahler order
#   max(subtree_so$segments)
#   
#   PD <- rbind(com_xyz[ind_eg2[j],], v1[ind_eg2[j],], v0[ind_eg2[j],]) #center and arrow
#   OD <- rbind(v3[ind_eg2[j],], v4[ind_eg2[j],]) #OD
#   
#   # transform to pc coord, xyz_ref_med[1,] as origin
#   subtree_pc <- subtree
#   subtree_pc$d[, c("X","Y","Z")] <- sweep(xyzmatrix(subtree$d), 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM
#   root_xyz_pc <- as.numeric(root_xyz - xyz_ref_med[1,]) %*% pc$rotation %*% rotM 
#   PD_pc <- sweep(PD, 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM #center and arrow
#   OD_pc <- sweep(OD, 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM #center and arrow
#   nb_pc <- sweep(med_xyz[ii2,], 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM # nb col
#   
#   # # DEBUG
#   # nopen3d()
#   # plot3d(subtree_pc,lwd=2, col='gray')
#   # # points3d(nb_pc, size=20, col=pal_T4, alpha=0.5)
#   # points3d(xyz_ref_med_pc_rot, size=15, col=pal_T4)
#   # points3d(PD_pc, size=10, col=pal_T4)
#   # axes3d(c('x','y','z')); title3d('','','x','y','z')
#   
#   # PLOT
#   windows(width = 9, height = 12)
#   # png(filename = paste("SO_T4", letters[LL], "_", j, ".png", sep = ''))
#   # pdf(paste("SO_T4", letters[LL], "_", j, ".pdf", sep = ''), width = 9, height = 12)
#   plot(nb_pc[1:7, 1:2], xlim= range(nb_pc[1:7, 1])*1.5, ylim= range(nb_pc[1:7, 2])*1.5, col= scales::alpha('gray20', 1), pch=1, lwd=2, cex= 10, asp= 1, xlab='', ylab ='', xaxt="n", yaxt="n", bty='n')
#   for (k in 1:max(subtree_so$points)) {
#     pt_p <- subtree_pc$d[,"PointNo"] %in% unique(subtree_pc$d[subtree_so$points == k,"Parent"])
#     pt_so <- pt_p | subtree_so$points==k
#     if (sum(pt_so) > 1) {
#       plot(subset(subtree_pc, pt_so), col=pal_so[k], add = T, lwd = 4, WithNodes = F)
#     }
#   }
#   points(root_xyz_pc[1], root_xyz_pc[2], pch=16, cex=2, col='black')
#   # points(PD_pc[1,1], PD_pc[1,2], pch=16, cex=5, col='gray')
#   # points(PD_pc[2,1],PD_pc[2,2], pch=16, cex=4, col=pal_so[2])
#   # points(PD_pc[3,1],PD_pc[3,2], pch=16, cex=4, col=pal_so[4])
#   arrows(PD_pc[2,1], PD_pc[2,2], PD_pc[3,1], PD_pc[3,2], lwd= 4,)
#   arrows(OD_pc[1,1], OD_pc[1,2], OD_pc[2,1], OD_pc[2,2], lwd= 4, code = 3, angle = 90, length = 0.2)
#   title(paste("SO_T4", letters[LL], "_", ind_eg2[j], sep = ''))
#   # segments(OD_pc[1,1], OD_pc[1,2], OD_pc[2,1], OD_pc[2,2], lend="square",lwd = 2)
#   segments(-5000,-10000,5000,-10000,lwd = 2) # 10 um
#   # dev.off()
# }
# 
# 
# 
# # more T4d eg with Strahler ----------------------------------------------------------------
# # choose nb via nb_ind
# # align to meridian by mapping ref points to med via np
# 
# colnames(ucl_rot_sm) <- c('x','y','z')
# 
# nb_coord <- rbind(c(0,0),
#                   c(-5/180*pi,0),
#                   c(+5/180*pi,0),
#                   c(0, -5/180*pi),
#                   c(0, +5/180*pi) )
# 
# LL <- 4 # choose type
# # eye
# com_xyz_eye <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
# # med
# com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
# v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix()
# v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
# v3 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
# v4 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()
# # nopen3d()
# # points3d(ucl_rot_sm)
# # points3d(com_xyz_eye, size=12, col='red',alpha=0.9)
# # identify3d(com_xyz_eye, adj = 1)
# 
# # 3 + 3 + 3
# ind_eg4 <- c(86, 98, 109, 80, 84, 71, 51, 87, 52)
# T4d_eg <- T4_dend[[LL]][ind_eg4]
# 
# for (j in 1:length(T4d_eg)) {
#   # 1+6+12 hex nb
#   ii <- sweep(ucl_rot_sm, 2, com_xyz_eye[ind_eg4[j], ] )^2 %>% rowSums() %>% which.min()
#   ii2 <- c(nb_ind[nb_ind[ii,],]) %>% unique()
#   
#   # ref nb
#   pt <- ucl_rot_sm[ii, ] %>% matrix(ncol=3)
#   tp <- cart2sphZ(pt)[,2:3]
#   rtp <- cbind(1, sweep(nb_coord,2,tp,'+'))
#   xyz_ref <- sph2cartZ(rtp)
#   
#   # np to med
#   xyz_ref_med <- matrix(nrow = nrow(xyz_ref), ncol = 3) 
#   xyz_eval <- data.frame(mc.x = xyz_ref[,1], mc.y = xyz_ref[,2], mc.z = xyz_ref[,3])
#   for (k in 1:3) {
#     npdata <- data.frame(mc = ucl_rot_sm[ii2,], ec = med_xyz[ii2,k])
#     bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
#     model_np <- npreg(bw)
#     xyz_ref_med[,k] <- predict(model_np, newdata = xyz_eval)
#   }
#   
#   # pc
#   zz <- cross3D((xyz_ref_med[2,]-xyz_ref_med[1,]), (xyz_ref_med[5,]-xyz_ref_med[1,])) #pointing inwards
#   pc <- prcomp(xyz_ref_med)
#   if (pc$rotation[,3] %*% zz < 0) {
#     pc$rotation <- -pc$rotation
#   }
#   if (c(cross3D(pc$rotation[,3], pc$rotation[,1])) %*% pc$rotation[,2] < 0) {
#     pc$rotation[,2] <- -pc$rotation[,2]
#   }
#   
#   xyz_ref_med_pc <- sweep(xyz_ref_med, 2, xyz_ref_med[1,]) %*% pc$rotation
#   y_ref <- diff(xyz_ref_med_pc[c(3,2),1:2])
#   ang <- acos(y_ref %*% c(0,1) / sqrt(sum(y_ref^2)))
#   ang <- if_else(y_ref[1] < 0, -ang, ang)
#   rotM <- matrix(c(cos(ang), sin(ang),0,
#                    -sin(ang), cos(ang), 0,
#                    0, 0, 1), ncol = 3, byrow = T)
#   xyz_ref_med_pc_rot <- xyz_ref_med_pc %*% rotM
#   
#   # T4
#   tar <- T4d_eg[[j]]
#   ind_D = match(tar$tags$`dendrite start` , tar$d$PointNo)
#   df_D <-  tar$d[ind_D,]
#   root_xyz <- tar$d[ind_D, c("X","Y","Z")]
#   
#   # - find the subtree with root = dendrite start
#   targ <- as.ngraph(tar)
#   ii_root <- ind_D
#   # subtree and Strahler order
#   sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
#   subtree <- subset(tar, sub_points) 
#   subtree_g <- as.ngraph(subtree, weights = T)
#   subtree_so <- strahler_order(subtree) # Strahler order
#   max(subtree_so$segments)
#   
#   PD <- rbind(com_xyz[ind_eg4[j],], v1[ind_eg4[j],], v0[ind_eg4[j],]) #center and arrow
#   OD <- rbind(v3[ind_eg4[j],], v4[ind_eg4[j],]) #OD
#   
#   # transform to pc coord, xyz_ref_med[1,] as origin
#   subtree_pc <- subtree
#   subtree_pc$d[, c("X","Y","Z")] <- sweep(xyzmatrix(subtree$d), 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM
#   root_xyz_pc <- as.numeric(root_xyz - xyz_ref_med[1,]) %*% pc$rotation %*% rotM 
#   PD_pc <- sweep(PD, 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM #center and arrow
#   OD_pc <- sweep(OD, 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM #center and arrow
#   nb_pc <- sweep(med_xyz[ii2,], 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM # nb col
#   
#   # PLOT
#   windows(width = 9, height = 12)
#   # png(filename = paste("SO_T4", letters[LL], "_", j, ".png", sep = ''))
#   # pdf(paste("SO_T4", letters[LL], "_", j, ".pdf", sep = ''), width = 9, height = 12)
#   plot(nb_pc[1:7, 1:2], xlim= range(nb_pc[1:7, 1])*1.6, ylim= range(nb_pc[1:7, 2])*1.6, col= scales::alpha('gray20', 1), pch=1, lwd=2, cex= 10, asp= 1, xlab='', ylab ='', xaxt="n", yaxt="n", bty='n')
#   for (k in 1:max(subtree_so$points)) {
#     pt_p <- subtree_pc$d[,"PointNo"] %in% unique(subtree_pc$d[subtree_so$points == k,"Parent"])
#     pt_so <- pt_p | subtree_so$points==k
#     if (sum(pt_so) > 1) {
#       plot(subset(subtree_pc, pt_so), col=pal_so[k], add = T, lwd = 4, WithNodes = F)
#     }
#   }
#   points(root_xyz_pc[1], root_xyz_pc[2], pch=16, cex=2, col='black')
#   # points(PD_pc[1,1], PD_pc[1,2], pch=16, cex=5, col='gray')
#   # points(PD_pc[2,1],PD_pc[2,2], pch=16, cex=4, col=pal_so[2])
#   # points(PD_pc[3,1],PD_pc[3,2], pch=16, cex=4, col=pal_so[4])
#   arrows(PD_pc[2,1], PD_pc[2,2], PD_pc[3,1], PD_pc[3,2], lwd= 4,)
#   arrows(OD_pc[1,1], OD_pc[1,2], OD_pc[2,1], OD_pc[2,2], lwd= 4, code = 3, angle = 90, length = 0.2)
#   title(paste("SO_T4", letters[LL], "_", ind_eg4[j], sep = ''))
#   # segments(OD_pc[1,1], OD_pc[1,2], OD_pc[2,1], OD_pc[2,2], lend="square",lwd = 2)
#   segments(-5000,-10000,5000,-10000,lwd = 2) # 10 um
#   # dev.off()
#   
# }


# 2D Angle of PD vs hex in med ref, without np ---------------------------------------------------
# check exist 19 nb, use 7 nb to determine angles

LL <- 2 # T4b

com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
com_xyz_eye <- lens_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix() 
v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
# v6 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
# v7 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()

T4b_north <- com_xyz_eye[,3] > 0
T4b_south <- com_xyz_eye[,3] < 0

cc <- c()
PD_ang_T4b_med <- matrix(ncol = 1, nrow = nrow(v0))
for (j in 1:nrow(v0)) {
  ii <- sweep(med_xyz, 2, com_xyz[j,])^2 %>% rowSums() %>% which.min()
  ii7 <- nb_ind[ii,]
  # ii19 <- c(ii7, nb_ind[ii7[3],2:3],nb_ind[ii7[6],c(3,6,4)],nb_ind[ii7[4],4:5], nb_ind[ii7[5],5],nb_ind[ii7[7],c(5,7,2)],nb_ind[ii7[2],2])
  ii19 <- c(ii7, nb_ind[ii7[2],2:3],nb_ind[ii7[3],c(3,4)],nb_ind[ii7[4],4:5], nb_ind[ii7[5],c(5,6)],nb_ind[ii7[6],c(6,7)],nb_ind[ii7[7],c(7,2)]) #2024
  cc <- c(cc, sum(!is.na(ii19)))
  
  if (sum(is.na(ii19)) == 0) { # consistent with common ref
    # vh <- colMeans(med_xyz[ii7[2:3],]) - colMeans(med_xyz[ii7[4:5],]) #+h
    # vv <- med_xyz[ii7[6],] - med_xyz[ii7[7],] #+v
    vh <- colMeans(med_xyz[ii7[c(4,5)],]) - colMeans(med_xyz[ii7[c(2,7)],]) #+h, 2024
    vv <- med_xyz[ii7[3],] - med_xyz[ii7[6],] #+v
    vPD <- v1[j,] - v0[j,] # PD
    
    # # 3D angle
    # ang <- acos(sum(vh * vPD) / sqrt(sum(vh^2)) / sqrt(sum(vPD^2)) ) /pi*180
    # ang <- if_else(t(cross3D(vv, vh)) %*% cross3D(vh, vPD) > 0, ang, -ang) # v -> h is positive angle
    
    # 2D angle
    # project onto plane
    tar_pca <- prcomp(xyzmatrix(med_xyz[ii19,])) 
    if ( tar_pca$rotation[,3] %*% cross3D(vh, vv) < 0 ) {
      tar_pca$rotation <- - tar_pca$rotation
    }
    # if (t(cross3D(tar_pca$rotation[,3],tar_pca$rotation[,1])) %*% tar_pca$rotation[,2] < 0 ) {
    #   tar_pca$rotation[,2] <- - tar_pca$rotation[,2]
    # }
    vv2 <- vv %*% tar_pca$rotation
    vv2 <- vv2[,1:2]
    vPD2 <- vPD %*% tar_pca$rotation
    vPD2 <- vPD2[,1:2]
    ang <- acos( sum(vv2 * vPD2) / sqrt(sum(vv2^2)) / sqrt(sum(vPD2^2)) ) /pi*180
    
    PD_ang_T4b_med[j] <- ang
  }
}
cc_b <- cc


LL <- 4 # T4d

com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
com_xyz_eye <- lens_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix() 
v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
# v6 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
# v7 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()

T4d_north <- com_xyz_eye[,3] > 0
T4d_south <- com_xyz_eye[,3] < 0

cc <- c()
PD_ang_T4d_med <- matrix(ncol = 1, nrow = nrow(v0))
for (j in 1:nrow(v0)) {
  ii <- sweep(med_xyz, 2, com_xyz[j,])^2 %>% rowSums() %>% which.min()
  ii7 <- nb_ind[ii,]
  ii19 <- c(ii7, nb_ind[ii7[2],2:3],nb_ind[ii7[3],c(3,4)],nb_ind[ii7[4],4:5], nb_ind[ii7[5],c(5,6)],nb_ind[ii7[6],c(6,7)],nb_ind[ii7[7],c(7,2)])
  cc <- c(cc, sum(!is.na(ii19)))
  
  if (sum(is.na(ii19)) == 0) { # consistent with common ref
    vh <- colMeans(med_xyz[ii7[c(4,5)],]) - colMeans(med_xyz[ii7[c(2,7)],]) #+h, 2024
    vv <- med_xyz[ii7[3],] - med_xyz[ii7[6],] #+v
    vPD <- v1[j,] - v0[j,] # PD
    
    # 3D
    ang <- acos(sum( (-vv) * vPD) / sqrt(sum(vv^2)) / sqrt(sum(vPD^2)) ) /pi*180
    ang <- if_else(t(cross3D(vv, vh)) %*% cross3D( (-vv), vPD) > 0, ang, -ang) # v -> h is positive angle
    
    # 2D angle
    # project onto plane
    tar_pca <- prcomp(xyzmatrix(med_xyz[ii19,])) 
    if ( tar_pca$rotation[,3] %*% cross3D(vh, vv) < 0 ) {
      tar_pca$rotation <- - tar_pca$rotation
    }
    # if (t(cross3D(tar_pca$rotation[,3],tar_pca$rotation[,1])) %*% tar_pca$rotation[,2] < 0 ) {
    #   tar_pca$rotation[,2] <- - tar_pca$rotation[,2]
    # }
    vh2 <- vh %*% tar_pca$rotation
    vh2 <- vh2[,1:2]
    vPD2 <- vPD %*% tar_pca$rotation
    vPD2 <- vPD2[,1:2]
    ang <- acos( sum(vh2 * vPD2) / sqrt(sum(vh2^2)) / sqrt(sum(vPD2^2)) ) /pi*180
    
    PD_ang_T4d_med[j] <- ang
  }
}
cc_d <- cc


# -- PLOT ang density
windows(width = 8, height = 4)
# pdf(paste("PD_ang_med.pdf", sep = ""), width = 8, height = 4)

# dd <- density(PD_ang_T4b_med-90, from= -45, to= 45, bw='SJ', adjust=1, na.rm = T)
dd <- density(PD_ang_T4b_med[T4b_south,]-90, from= -45, to= 45, bw='SJ', adjust=1, na.rm = T)
plot(dd$x, dd$y, type='l', bty='n', col=pal_T4[2], lwd=3, 
     xlim = c(-45, 45), 
     ylim = c(0, 0.03), 
     xaxt='n',yaxt='n', xlab ="", ylab='',main = '')
# dd <- density(PD_ang_T4d_med-90, from= -45, to= 45, bw='SJ', adjust=1,na.rm = T)
dd <- density(PD_ang_T4d_med[T4d_north,]-90, from= -45, to= 45, bw='SJ', adjust=1,na.rm = T)
points(dd$x, dd$y, type='l', bty='n',col=pal_T4[4], lwd=3)
axis(1, at = seq(-45, 45, by =45), labels = paste(seq(-45, 45, by =45), "°", sep = ''), cex.axis = 1.5 )
axis(2, at = c(0, 1/512/diff(dd$x[1:2]) ), labels = c('0', '1'), cex.axis = 1.5 )
# dev.off()


sum(na.omit(PD_ang_T4d_med) < 0)
t.test(na.omit(PD_ang_T4b_med), mu = 0)
t.test(na.omit(PD_ang_T4d_med), mu = 0)

  
# Fig.2H, RF in a common med ref, T4b ---------------------------------------------------

# match nbco / nb_ind order
refhex_1 <- rbind(c(0,0),
                  c(-1,+1),
                  c(0, +2),
                  c(+1,+1),
                  c(+1,-1),
                  c(0, -2),
                  c(-1,-1) )

# this is viewed from distal side
refhex_2 <- rbind(c(0,0),
                  c(-1,+1),
                  c(0, +2),
                  c(+1,+1),
                  c(+1,-1),
                  c(0, -2),
                  c(-1,-1),
                  c(-2,+2),
                  c(-1,+3),
                  c(0,+4),
                  c(+1,+3),
                  c(+2,+2),
                  c(+2,0),
                  c(+2,-2),
                  c(+1,-3),
                  c(0,-4),
                  c(-1,-3),
                  c(-2,-2),
                  c(-2,0) )


dev.new()
plot(refhex_2, asp=1)
text(refhex_2, labels = seq(1,nrow(refhex_2)), adj = 2)


# rescaling factor, med -> refhex
rs <- mean(nb_dist_med, na.rm =T) / mean(sqrt(rowSums(refhex_1[-1,]^2)))
rs <- 1

LL <- 2 

com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
com_xyz_eye <- lens_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix() 
v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
v6 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
v7 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()

# - np
v3 <- matrix(ncol = 2, nrow = nrow(v0))  #base of vec of T4_dir (orientation not PD).
v4 <- matrix(ncol = 2, nrow = nrow(v1)) #head of vec
v5 <- matrix(ncol = 2, nrow = nrow(com_xyz)) # for com
v8 <- matrix(ncol = 2, nrow = nrow(v6))  #base of vec
v9 <- matrix(ncol = 2, nrow = nrow(v7)) #head of vec
cc <- c()
for (j in 1:nrow(v0)) {
  ii <- sweep(med_xyz, 2, com_xyz[j,])^2 %>% rowSums() %>% which.min()
  ii7 <- nb_ind[ii,]
  # ii19 <- c(ii7, nb_ind[ii7[3],2:3],nb_ind[ii7[6],c(3,6,4)],nb_ind[ii7[4],4:5], nb_ind[ii7[5],5],nb_ind[ii7[7],c(5,7,2)],nb_ind[ii7[2],2])
  ii19 <- c(ii7, nb_ind[ii7[2],2:3],nb_ind[ii7[3],c(3,4)],nb_ind[ii7[4],4:5], nb_ind[ii7[5],c(5,6)],nb_ind[ii7[6],c(6,7)],nb_ind[ii7[7],c(7,2)]) #2024
  cc <- c(cc, sum(!is.na(ii19)))
  
  # nbs
  iinna <- !is.na(ii19)
  nb_xyz <- med_xyz[ii19[iinna], ]
  
  # pca
  nb_pca <- prcomp(nb_xyz)
  nb_xyz_pc <- sweep(nb_xyz, 2,nb_pca$center) %*% nb_pca$rotation
  colnames(nb_xyz_pc) <- c('x','y','z')
  nb_xy <- nb_xyz_pc[, 1:2] / rs
  
  vec_xyz <- rbind(v0[j,], v1[j,])
  vec_xyz_pc <- sweep(vec_xyz, 2,nb_pca$center) %*% nb_pca$rotation
  colnames(vec_xyz_pc) <- c('x','y','z')
  vec_xy <- vec_xyz_pc[, 1:2] / rs
  
  vec_xyz_od <- rbind(v6[j,], v7[j,])
  vec_xyz_od_pc <- sweep(vec_xyz_od, 2,nb_pca$center) %*% nb_pca$rotation
  colnames(vec_xyz_od_pc) <- c('x','y','z')
  vec_xy_od <- vec_xyz_od_pc[, 1:2] / rs
  
  com_xyz_pc <- (com_xyz[j,] - nb_pca$center) %*% nb_pca$rotation
  colnames(com_xyz_pc) <- c('x','y','z')
  com_xy <- com_xyz_pc[1:2] / rs
  
  # -- np
  xyz_vtail_eval <- data.frame(mc.x = vec_xy[1,1], mc.y = vec_xy[1,2])
  xyz_vhead_eval <- data.frame(mc.x = vec_xy[2,1], mc.y = vec_xy[2,2])
  xyz_vtail_eval_od <- data.frame(mc.x = vec_xy_od[1,1], mc.y = vec_xy_od[1,2])
  xyz_vhead_eval_od <- data.frame(mc.x = vec_xy_od[2,1], mc.y = vec_xy_od[2,2])
  xyz_com_eval <- data.frame(mc.x = com_xy[1], mc.y = com_xy[2])
  for (k in 1:2) {
    npdata <- data.frame(mc = nb_xy, ec = refhex_2[iinna,k])
    bw <- npregbw(formula= ec~mc.x+mc.y, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
    model_np <- npreg(bw)
    v3[j,k] <- predict(model_np, newdata = xyz_vtail_eval)
    v4[j,k] <- predict(model_np, newdata = xyz_vhead_eval)
    v5[j,k] <- predict(model_np, newdata = xyz_com_eval)
    v8[j,k] <- predict(model_np, newdata = xyz_vtail_eval_od)
    v9[j,k] <- predict(model_np, newdata = xyz_vhead_eval_od)
  }
}

df_b <- rbind(cbind(v4,1), cbind(v3,2), cbind(v5,3) ) %>% as.data.frame() #reverse arrow st. it's PD
colnames(df_b) <- c('x','y', 'gp')
df_b$gp <- factor(df_b$gp)
PD_hex_b <- cbind(v4, v3)
PD_com_b <- v5
OD_hex_b <- cbind(v9, v8)
cc_b <- cc

# PLOT
df0 <- as.data.frame(refhex_2)
colnames(df0) <- c('x','y')

df <- df_b
cc <- cc_b
ii <- cc==19
ii3 <- rep(ii, 3)
df[!ii3, ] <- NA

df_arrow <- as.data.frame(rbind(cbind(PD_hex_b[ii,1:2], PD_com_b[ii,]),
                                  cbind(PD_com_b[ii,], PD_hex_b[ii,3:4])) )
colnames(df_arrow) <- c('x','y','xend','yend')

# -- plot
windows(width = 7, height = 10)
plt <- ggplot() +
  geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour='gray',size =0.2) +
  geom_point(data=df0, aes(x,y), size=30, pch=1) +
  geom_point(data=df, aes(x=x,y=y, col=gp), size=3, na.rm = T ) +
  scale_color_manual(values = c(scales::alpha('violetred', 0.7), scales::alpha('cyan3', 0.7), scales::alpha('gray20', 0.5)), guide= guide_legend(title="")) +
  # scale_shape_manual(values=c(2, 1, 3))+
  ylab("y") +   xlab("x") +
  theme_minimal() +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3,3), labels = seq(-3,3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4, 4), expand = c(0, 0), breaks = seq(-5,5)) + # set +y as above eq
  labs(title = paste("T4", letters[LL],  sep = "")) +
  theme(axis.title = element_blank(), panel.grid.major = element_line(colour = "gray70"), panel.grid.minor =  element_blank() )+
  coord_fixed(ratio=1)
# pdf(paste("hex_ref_", letters[LL], '.pdf', sep = ""), width = 7, height = 10)
plt 
# dev.off()


# cont. T4d ----------------------------------------------

LL <- 4 # choose type
com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
com_xyz_eye <- lens_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix() 
v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
v6 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
v7 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()

# - np
v3 <- matrix(ncol = 2, nrow = nrow(v0))  #base of vec of T4_dir (orientation not PD).
v4 <- matrix(ncol = 2, nrow = nrow(v1)) #head of vec
v5 <- matrix(ncol = 2, nrow = nrow(com_xyz)) # for com
v8 <- matrix(ncol = 2, nrow = nrow(v6))  #base of vec
v9 <- matrix(ncol = 2, nrow = nrow(v7)) #head of vec
cc <- c()
for (j in 1:nrow(v0)) {
  ii <- sweep(med_xyz, 2, com_xyz[j,])^2 %>% rowSums() %>% which.min()
  ii7 <- nb_ind[ii,]
  ii19 <- c(ii7, nb_ind[ii7[2],2:3],nb_ind[ii7[3],c(3,4)],nb_ind[ii7[4],4:5], nb_ind[ii7[5],c(5,6)],nb_ind[ii7[6],c(6,7)],nb_ind[ii7[7],c(7,2)]) #2024
  cc <- c(cc, sum(!is.na(ii19)))
  # find nbs
  iinna <- !is.na(ii19)
  nb_xyz <- med_xyz[ii19[iinna], ]
  # pca alignment
  nb_pca <- prcomp(nb_xyz)
  nb_xyz_pc <- sweep(nb_xyz, 2,nb_pca$center) %*% nb_pca$rotation
  colnames(nb_xyz_pc) <- c('x','y','z')
  nb_xy <- nb_xyz_pc[, 1:2] / rs
  
  vec_xyz <- rbind(v0[j,], v1[j,])
  vec_xyz_pc <- sweep(vec_xyz, 2,nb_pca$center) %*% nb_pca$rotation
  colnames(vec_xyz_pc) <- c('x','y','z')
  vec_xy <- vec_xyz_pc[, 1:2] / rs
  
  vec_xyz_od <- rbind(v6[j,], v7[j,])
  vec_xyz_od_pc <- sweep(vec_xyz_od, 2,nb_pca$center) %*% nb_pca$rotation
  colnames(vec_xyz_od_pc) <- c('x','y','z')
  vec_xy_od <- vec_xyz_od_pc[, 1:2] / rs
  
  com_xyz_pc <- (com_xyz[j,] - nb_pca$center) %*% nb_pca$rotation
  colnames(com_xyz_pc) <- c('x','y','z')
  com_xy <- com_xyz_pc[1:2] / rs
  
  # -- np
  xyz_vtail_eval <- data.frame(mc.x = vec_xy[1,1], mc.y = vec_xy[1,2])
  xyz_vhead_eval <- data.frame(mc.x = vec_xy[2,1], mc.y = vec_xy[2,2])
  xyz_vtail_eval_od <- data.frame(mc.x = vec_xy_od[1,1], mc.y = vec_xy_od[1,2])
  xyz_vhead_eval_od <- data.frame(mc.x = vec_xy_od[2,1], mc.y = vec_xy_od[2,2])
  xyz_com_eval <- data.frame(mc.x = com_xy[1], mc.y = com_xy[2])
  for (k in 1:2) {
    npdata <- data.frame(mc = nb_xy, ec = refhex_2[iinna,k])
    bw <- npregbw(formula= ec~mc.x+mc.y, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
    model_np <- npreg(bw)
    v3[j,k] <- predict(model_np, newdata = xyz_vtail_eval)
    v4[j,k] <- predict(model_np, newdata = xyz_vhead_eval)
    v5[j,k] <- predict(model_np, newdata = xyz_com_eval)
    v8[j,k] <- predict(model_np, newdata = xyz_vtail_eval_od)
    v9[j,k] <- predict(model_np, newdata = xyz_vhead_eval_od)
  }
}

df_d <- rbind(cbind(v4,1), cbind(v3,2), cbind(v5,3) ) %>% as.data.frame()
colnames(df_d) <- c('x','y', 'gp')
df_d$gp <- factor(df_d$gp)
PD_hex_d <- cbind(v4, v3)
PD_com_d <- v5
OD_hex_d <- cbind(v9, v8)
cc_d <- cc

# PLOT
df0 <- as.data.frame(refhex_2)
colnames(df0) <- c('x','y')

df <- df_d
cc <- cc_d
ii <- cc==19
ii3 <- rep(ii, 3)
df[!ii3, ] <- NA

df_arrow <- as.data.frame(rbind(cbind(PD_hex_d[ii,1:2], PD_com_d[ii,]),
                                cbind(PD_com_d[ii,], PD_hex_d[ii,3:4])) )
colnames(df_arrow) <- c('x','y','xend','yend')

# -- plot
windows(width = 7, height = 10)
plt <- ggplot() +
  geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour='gray',size =0.2) +
  geom_point(data=df0, aes(x,y), size=30, pch=1) +
  geom_point(data=df, aes(x=x,y=y, col=gp), size=3, na.rm = T ) +
  scale_color_manual(values = c(scales::alpha('violetred', 0.7), scales::alpha('cyan3', 0.7), scales::alpha('gray20', 0.5)), guide= guide_legend(title="")) +
  ylab("y") +   xlab("x") +
  theme_minimal() +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3,3), labels = seq(-3,3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4, 4), expand = c(0, 0), breaks = seq(-5,5)) + # set +y as above eq
  labs(title = paste("T4", letters[LL],  sep = "")) +
  theme(axis.title = element_blank(), panel.grid.major = element_line(colour = "gray70"), panel.grid.minor =  element_blank() )+
  coord_fixed(ratio=1)
# pdf(paste("hex_ref_", letters[LL], '.pdf', sep = ""), width = 7, height = 10)
plt 
# dev.off()

# cont., ED Fig.3D, Fig.2J,Fig.2K, PD amplitude and angle in common hex ref ---------------

# T4b
ii <- cc_b ==19
PD_hex <- PD_hex_b[ii, ]
OD_hex <- OD_hex_b[ii, ]
PDOD_T4b <- cbind(sqrt(rowSums((PD_hex[,1:2] - PD_hex[,3:4])^2)), sqrt(rowSums((OD_hex[,1:2] - OD_hex[,3:4])^2)) )
PD_ang_T4b <- atan2(PD_hex[,4] - PD_hex[,2], PD_hex[,3] - PD_hex[,1]) /pi*180
PD_ang_T4b <- 90 - PD_ang_T4b # wrt +v
# T4d
ii <- cc_d ==19
PD_hex <- PD_hex_d[ii, ]
OD_hex <- OD_hex_d[ii, ]
PDOD_T4d <- cbind(sqrt(rowSums((PD_hex[,1:2] - PD_hex[,3:4])^2)), sqrt(rowSums((OD_hex[,1:2] - OD_hex[,3:4])^2)) )
PD_ang_T4d <- atan2(PD_hex[,4] - PD_hex[,2], PD_hex[,3] - PD_hex[,1]) /pi*180
PD_ang_T4d <- 90 - PD_ang_T4d

# scatter plot
df <- rbind(cbind(PDOD_T4b, 2), cbind(PDOD_T4d, 4) )
df <- as.data.frame(df)
colnames(df) <- c('PD','OD','gp')
df$gp <- factor(df$gp)
plt <- ggplot(df) +
  geom_point(aes(x = PD, y = OD, colour = gp), size = 3, alpha=0.8) +
  scale_color_manual(values= pal_T4[c(2,4)],labels = c("T4b", "T4d"), guide= guide_legend(title="dist[norm]"), na.value="gray") +
  ylab("OD") + xlab("PD") +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid.minor = element_blank() )+
  scale_x_continuous(limits = c(2, 6), breaks = seq(2,6), labels = seq(2,6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 6), breaks = seq(2,6), labels = seq(2,6), expand = c(0, 0)) + # set +y as above eq
  labs(title = paste("OD vs PD_med", "_norm_unit", sep = "")) +
  coord_fixed(ratio=1)
windows(width = 8, height = 8)
# pdf("OD_PD_hex.pdf",width = 8, height = 8)
# pdf("OD_PD_hex_sepa_norm.pdf",width = 8, height = 8)
# plt
# ggMarginal(plt, margins = "y", size = 2, type = "boxplot", outlier.size =3, groupColour = TRUE, groupFill = TRUE)
ggMarginal(plt, margins = "both", size = 4, type = "density", groupColour = TRUE, groupFill = F, lwd=2)
# dev.off()

# # -- ang hist
# windows(width = 8, height = 4)
# # pdf(paste("PD_ang_hex.pdf", sep = ""), width = 8, height = 4)
# hh <- hist(PD_ang_T4b, breaks = seq(-45, 45, by=5)+90, plot = F)
# plot(hh$mids, hh$density, type='l', bty='n', col=pal_T4[2], lwd=3, xlim = c(-45, 45+90), xaxt='n',yaxt='n',
#      xlab ="", ylab='',main = '')
# hh <- hist(PD_ang_T4d, breaks = seq(-45, 45, by=5), plot = F)
# points(hh$mids, hh$density, type='l', bty='n',col=pal_T4[4], lwd=3)
# axis(1, at = seq(-45, 45+90, by =45), labels = paste(seq(-45, 45+90, by =45), "°", sep = ''), cex.axis = 1.5 )
# axis(2, at = c(0, 1/5/length(hh$mids) ), labels = c('0', '1/n'), cex.axis = 1.5 )
# # dev.off()

# -- ang density
windows(width = 8, height = 4)
# pdf(paste("PD_ang_hex.pdf", sep = ""), width = 8, height = 4)
dd <- density(PD_ang_T4b, from= -45+90, to= 45+90, bw='SJ')
plot(dd$x, dd$y, type='l', bty='n', col=pal_T4[2], lwd=3, 
     xlim = c(-45+90, 45+180), ylim = c(0,0.03),
     xaxt='n',yaxt='n',xlab ="", ylab='',main = '')
dd <- density(PD_ang_T4d, from= -45+180, to= 45+180, bw='SJ')
points(dd$x, dd$y, type='l', bty='n',col=pal_T4[4], lwd=3)
axis(1, at = seq(-45+90, 45+180, by =45), labels = paste(seq(-45+90, 45+180, by =45), "°", sep = ''), cex.axis = 1.5 )
axis(2, at = c(0, 1/512/diff(dd$x[1:2]) ), labels = c('0', '1'), cex.axis = 1.5 )
# dev.off()


sum(PD_ang_T4d > 180)
sum(PD_ang_T4b < 90)
t.test(PD_ang_T4b, mu = 90)
t.test(PD_ang_T4d, mu = 0)

# -- PD amp density
windows(width = 8, height = 4)
# pdf(paste("PD_amp_hex.pdf", sep = ""), width = 8, height = 4)
dd <- density(PDOD_T4b[,1], from= 2, to= 6, bw='SJ')
plot(dd$x, dd$y, type='l', bty='n', col=pal_T4[2], lwd=3, xlim = c(2, 6), xaxt='n',yaxt='n',
     xlab ="", ylab='',main = '')
dd <- density(PDOD_T4d[,1], from= 2, to= 6, bw='SJ')
points(dd$x, dd$y, type='l', bty='n',col=pal_T4[4], lwd=3)
axis(1, at =  seq(2, 6, by =1), labels =  seq(2, 6, by =1), cex.axis = 1.5 )
axis(2, at = c(0, 1/512/diff(dd$x[1:2]) ), labels = c('0', '1'), cex.axis = 1.5 )
# dev.off()

wilcox.test(PDOD_T4b[,1], PDOD_T4d[,1])
t.test(PDOD_T4b[,1], PDOD_T4d[,1])


# cont. mean and sd -------------------------------------------------------
# 
# xy <- df_b[df_b$gp == 1, 1:2] #purple, tail
# colMeans(xy)
# sd(xy[,1])/col_dia
# sd(xy[,2])/col_dia
# 
# xy <- df_b[df_b$gp == 2, 1:2] #cyan head
# colMeans(xy)
# sd(xy[,1])/col_dia
# sd(xy[,2])/col_dia
# 
# 
# xy <- df_d[df_d$gp == 1, 1:2] #purple, tail
# colMeans(xy)
# sd(xy[,1])/col_dia
# sd(xy[,2])/col_dia
# 
# xy <- df_d[df_d$gp == 2, 1:2] #cyan head
# colMeans(xy)
# sd(xy[,1])/col_dia
# sd(xy[,2])/col_dia


# # cont.  spacial spread ---------------------------------------------------------
# col_dia <- sqrt(2) #diameter/nb-distance in regular grid
# 
# ## ##
# mul <- 1
# mul <- sqrt(2*log(2)) * 2 #convert to FWHM
# 
# ## ##
# df <- df_b
# df <- df_d
#   
# PD_ptdistr <- matrix(ncol = 3, nrow = 4) # spatial spred
# colnames(PD_ptdistr) <- c('PDtail','PDhead','center')
# rownames(PD_ptdistr) <- c('h','v','dh','dv')
# # PD tail
# fnorm <- MASS::fitdistr(x= df$x[df$gp==1], "normal")
# mux <- fnorm$estimate['mean']
# wx <- fnorm$estimate['sd'] * mul
# fnorm <- MASS::fitdistr(x= df$y[df$gp==1], "normal")
# muy <- fnorm$estimate['mean']
# wy <- fnorm$estimate['sd'] * mul
# PD_ptdistr[,1] <- c(mux, muy, wx, wy)
# # PD head
# fnorm <- MASS::fitdistr(x= df$x[df$gp==2], "normal")
# mux <- fnorm$estimate['mean']
# wx <- fnorm$estimate['sd'] * mul
# fnorm <- MASS::fitdistr(x= df$y[df$gp==2], "normal")
# muy <- fnorm$estimate['mean']
# wy <- fnorm$estimate['sd'] * mul
# PD_ptdistr[,2] <- c(mux, muy, wx, wy)
# # PD com
# fnorm <- MASS::fitdistr(x= df$x[df$gp==3], "normal")
# mux <- fnorm$estimate['mean']
# wx <- fnorm$estimate['sd'] * mul
# fnorm <- MASS::fitdistr(x= df$y[df$gp==3], "normal")
# muy <- fnorm$estimate['mean']
# wy <- fnorm$estimate['sd'] * mul
# PD_ptdistr[,3] <- c(mux, muy, wx, wy)
# 
# 
# # plot
# h <- MASS::kde2d(x= df_b$x[df_b$gp==1], y=df_b$y[df_b$gp==1],
#                  n =c(20,30), lims = c(-2, -1, -1, 1.5))
# dev.new()
# image(h)
# 
# # # - fit_gaussian_2D NOT working
# # df <- cbind(expand.grid(h$x, h$y), h$z)
# # colnames(df) <- c('X_values', 'Y_values', 'response')
# # gauss_fit_cir <- gaussplotR::fit_gaussian_2D(df, method = "circular")
# 
# # # - test 
# # sigma <- 10
# # w <- sqrt(2*log(2)) * sigma * 2 #FWHM
# # mu <- 100
# # x <- rnorm(100, mean = mu, sd = sigma)
# # # bw.nrd(x) #NO
# # MASS::fitdistr(x, "normal")
# # 
# # # https://stackoverflow.com/questions/37660950/how-to-create-2d-data-set-from-gaussian-distribution-in-r
# # # p <- 10  ## we want p-dimensional multivariate normal
# # # set.seed(0); X <- matrix(runif(p * p), p, p)  ## this random matrix has full rank
# # # 
# # # COV <- crossprod(X)  ## t(X) %*% X but about 2 times faster
# # # mu <- rep(0, p)
# # # 
# # # x <- MASS::mvrnorm(1000, mu, COV)  ## mvrnorm(sample.size, mean, covariance)
# 
# # - PLOT
# # todo, add distr plot along h/v axis for b/d


# ED Fig.3B, 3x3 compartments --------------------------------------------------------

com_xyz_eye_b <- lens_type[[2]][, c('comx','comy','comz')] %>% as.matrix()
com_xyz_eye_d <- lens_type[[4]][, c('comx','comy','comz')] %>% as.matrix()

# fnnb <- 13
fnnb <- 19

# - T4b
xyz <- com_xyz_eye_b
xyz[,2] <- -xyz[,2]
tp <- cart2sphZ(xyz)[,2:3] /pi*180
tp[,2] <- if_else(tp[,2] > 180, tp[,2] -360, tp[,2])

angdiv <- list()
inddiv <- list()
# tdiv <- c(0, 75, 105, 160)
# pdiv <- c(-15, 25, 65, 150)
tdiv <- c(0, 75, 100, 160)
pdiv <- c(-15, 35, 70, 150)
for (j in 1:3) {
  for (k in 1:3) {
    ii <- cc_b >= fnnb &
      tp[,1] > tdiv[j] & tp[,1] < tdiv[j+1] &
      tp[,2] > pdiv[k] & tp[,2] < pdiv[k+1]
    inddiv <- c(inddiv, list(which(ii)))
    PD_hex <- PD_hex_b[ii, ]
    ang <- atan2(PD_hex[,4] - PD_hex[,2], PD_hex[,3] - PD_hex[,1]) /pi*180
    ang <- 90 - ang #wrt +v
    angdiv <- c(angdiv, list(ang))
  }
}

# # DEBUG
# lapply(angdiv, length) %>% unlist()
# 
# nopen3d()
# points3d(ucl_rot_sm)
# for (j in 1:length(inddiv)) {
#   pt <- xyz[inddiv[[j]],]
#   pt[,2] <- -pt[,2]
#   points3d(pt, col=pal_9[j], size=10)
# }

# nopen3d()
# points3d(med_xyz)
# for (j in 1:nrow(v0)) {
#   if (ii[j]) {
#     arrow3d(v0[j,], v1[j,], theta = pi/6, n = 8, col =pal_T4[LL], type = "rotation", lit=F)
#   }
# }
# for (j in seq(-12, 12, by=2)) {
#   points3d(med_xyz[match(haxis_gen(j), eyemap[,2]),], col='red', size=15)
#   # points3d(med_xyz[match(vaxis_gen(j), eyemap[,2]),], col='red', size=15)
# }

# -- 3x3 plot
df0 <- as.data.frame(refhex_2)
colnames(df0) <- c('x','y')

windows(width = 9, height = 9)
# pdf(paste("hex_ref_b", '_3x3.pdf', sep = ""), width = 9, height = 9)
layout(matrix(seq(1,9),nrow = 3, ncol = 3, byrow = TRUE))
for (k in 1:length(inddiv)) {
  # plot(rbind(PD_hex_b[cc_b >= fnnb,1:2], PD_hex_b[cc_b >= fnnb,3:4]),
  #      ylim=c(-2.5,2.5), xlim=c(-3.5,3.5), mar=c(1,1,1,1),
  #      xlab ="", ylab='',main = '', asp = 1)
  plot(rbind(PD_hex_b[inddiv[[k]],1:2], PD_hex_b[inddiv[[k]],3:4]),
       ylim=c(-3.2,3.2), xlim=c(-3,3), mar=c(0.1,0.1,0.1,0.1),pch=16,col=pal_9[k],
       xlab ="", ylab='', xaxt="n",yaxt="n", frame.plot = F, main = '', asp = 1)
  points(df0, pch=1, cex=5.5)
  for (j in 1:nrow(PD_hex_b)) {
    if (j %in% inddiv[[k]]) {
      segments(PD_hex_b[j,1], PD_hex_b[j,2], PD_hex_b[j,3], PD_hex_b[j,4], col=pal_9[k])
    }
  }
  segments(0,0, v1[k,1]*3, v1[k,2]*3, col=pal_9[k], lwd=2)
}
# dev.off()


# -- avg vector
v1 <- matrix(ncol = 2, nrow = 0)
for (j in 1:length(inddiv)) {
  vv <- data.frame(PD_hex_b[inddiv[[j]],])
  colnames(vv) <- c('x0','y0','xd','yd')
  vv %<>% as_tibble() %>%
    mutate(ang = atan2(yd-y0, xd-x0)) %>%
    as.data.frame()
  v1 <- rbind(v1, c(cos(mean(vv$ang)), sin(mean(vv$ang))))
}

# PLOT
windows(width = 6, height = 6)
# pdf(paste("hex_ref_b", '_3x3_vectors.pdf', sep = ""), width = 6, height = 6)
plot(c(0,0), type="n",
     ylim=c(-1,1), xlim=c(-1,1), mar=c(0.1,0.1,0.1,0.1),
     xlab ="", ylab='', xaxt="n",yaxt="n", frame.plot = F, main = '', asp = 1)
for (j in 1:9) {
  segments(0,0, v1[j,1], v1[j,2], col=pal_9[j])
}
# dev.off()


# -- hist
dev.new()
plot(c(0,0), type = 'n', xlim = c(40, 150), ylim = c(0,0.05))
for (j in 1:length(angdiv)) {
  hh <- hist(angdiv[[j]], breaks = seq(45, 145, by=10), plot = F)
  lines(hh$mids, hh$density, lty=1, bty='n',col= pal_9[j], lwd=3)
}

# # 3x3 T4b in med ----------------------------------------------------------
# 
# T4b <- T4_dend[[2]]
# T4b_me <- T4_dend[[2]]
# for (j in 1:length(T4b)) {
#   tar <- T4b[[j]]
#   ind_T = match(tar$tags$"dendrite start" , tar$d$PointNo)
#   targ <- as.ngraph(tar)
#   ii_root <- ind_T
#   # subtree and resample
#   sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
#   T4b_me[[j]] <- subset(tar, sub_points) 
# }
# 
# nopen3d()
# par3d('windowRect' = c(100,100,1200,1200))
# for (k in 1:length(inddiv)) {
#   plot3d(T4b_me[inddiv[[k]]], lwd=2, col=pal_9[k])
# }
# plot3d(T4b_me[-unlist(inddiv)], col='gray')
# # shade3d(ME_msh, col='gray', alpha=0.1)
# rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(180/180*pi,1,0,0) %*%
#                 rotationMatrix(120/180*pi,0,1,0) %*%
#                 rotationMatrix(20/180*pi,1,0,0))
# # scale bar
# c(10000,0,0) %>% 
#   rotate3d(180/180*pi,1,0,0) %>%
#   rotate3d(120/180*pi,0,1,0) %>%
#   rotate3d(20/180*pi,1,0,0) %>%
#   rbind(c(0,0,0), .) %>%
#   sweep(2, c(3e5, 3.5e5, 2.5e5), '+') %>%
#   segments3d(lwd=4)
# # rgl.snapshot(filename = "T4b_me_3x3.png")


# # 3x3 Mercator in eye coord, real T4  --------------------------------------------------------
# 
# LL <- 2 # choose type
# # eye
# com_xyz <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
# v0 <- as.matrix(lens_type[[LL]][, c("x0","y0","z0")])
# v1 <- as.matrix(lens_type[[LL]][, c("xd","yd","zd")])
# 
# v0_Mer <- sweep(v0, 1, sqrt(rowSums(v0^2)), '/')
# v0_Mer <- cart2Mercator(v0_Mer)
# colnames(v0_Mer) <- c('x','y')
# v1_Mer <- sweep(v1, 1, sqrt(rowSums(v1^2)), '/')
# v1_Mer <- cart2Mercator(v1_Mer)
# colnames(v1_Mer) <- c('x','y')
# 
# # add a v-axis
# ind_axis <- vaxis_gen(5)
# xy <- data.frame(ucl_rot_sm)[match(ind_axis, eyemap[,2]),]
# ax_Mer <- cart2Mercator(xy)
# 
# # PLOT
# df_arrow <- data.frame(cbind(v0_Mer, v1_Mer) )
# colnames(df_arrow) <- c('x','y','xend','yend')
# 
# 
# plt <- plt_Mer +
#   geom_path(data = cmer_Merc, aes(x=x, y=y), colour = pal_axes[3], lwd=1) +
#   labs(title = "")
# 
# for (j in 1:length(inddiv)) {
#   plt <- plt +
#     geom_segment(data=df_arrow[inddiv[[j]],], aes(x = x,y = y, xend = xend,yend = yend),
#                  colour=pal_9[j],size =1)
# }
# 
# windows(width = 7, height = 8.5)
# plt
# 
# # - ommatidia in Mercator with axes
# plt <- plt_Mer +
#   geom_point(data=ucl_rot_Merc, aes(x=x, y=y))
# 
# for (j in seq(-15,14, by=3)) {
#   ind_axis <- vaxis_gen(j)
#   xy <- data.frame(ucl_rot_Merc)[match(ind_axis, eyemap[,2]),]
#   mer_Merc <- xy[order(xy$y),]
#   plt <- plt + geom_path(data = mer_Merc, aes(x=x, y=y), colour = 'gray30', lwd=0.5)
# }
# for (j in seq(-29,30, by=4)) {
#   ind_axis <- haxis_gen(j)
#   xy <- data.frame(ucl_rot_Merc)[match(ind_axis, eyemap[,2]),]
#   mer_Merc <- xy[order(xy$x),]
#   plt <- plt + geom_path(data = mer_Merc, aes(x=x, y=y), colour = 'gray30', lwd=0.5)
# }
# 
# windows(width = 7, height = 8.5)
# plt
# 
# # ggsave(paste("v-axis_Mercator.pdf",sep=''), width = 5, height = 6.5)

# # cont. T4d ---------------------------------------------------------------
# 
# xyz <- com_xyz_eye_d
# xyz[,2] <- -xyz[,2]
# tp <- cart2sphZ(xyz)[,2:3] /pi*180
# tp[,2] <- if_else(tp[,2] > 180, tp[,2] -360, tp[,2])
# 
# angdiv <- list()
# inddiv <- list()
# tdiv <- c(0, 75, 105, 160)
# pdiv <- c(-15, 25, 65, 150)
# for (j in 1:3) {
#   for (k in 1:3) {
#     ii <- cc_d >= fnnb &
#       tp[,1] > tdiv[j] & tp[,1] < tdiv[j+1] &
#       tp[,2] > pdiv[k] & tp[,2] < pdiv[k+1]
#     inddiv <- c(inddiv, list(which(ii)))
#     PD_hex <- PD_hex_d[ii,,drop=F ]
#     ang <- atan2(PD_hex[,4] - PD_hex[,2], PD_hex[,3] - PD_hex[,1]) /pi*180
#     ang <- 90 - ang #wrt +v
#     angdiv <- c(angdiv, list(ang))
#   }
# }
# 
# # DEBUG
# lapply(angdiv, length) %>% unlist()
# 
# # -- 3x3 plot
# windows(width = 9, height = 9)
# # pdf(paste("hex_ref_d", '_3x3.pdf', sep = ""), width = 9, height = 9)
# layout(matrix(seq(1,9),nrow = 3, ncol = 3, byrow = TRUE))
# for (k in 1:length(inddiv)) {
#   plot(rbind(PD_hex_d[cc_d >= fnnb,1:2], PD_hex_d[cc_d >= fnnb,3:4]),
#        ylim=c(-4.5,4.5), asp = 1,xlab ="", ylab='',main = '')
#   for (j in 1:nrow(PD_hex_d)) {
#     if (j %in% inddiv[[k]]) {
#       segments(PD_hex_d[j,1], PD_hex_d[j,2], PD_hex_d[j,3], PD_hex_d[j,4], col=pal_9[k])
#     }
#   }
# }
# # dev.off()
# 
# # -- hist
# dev.new()
# plot(c(0,0), type = 'n', xlim = c(40, 150), ylim = c(0,0.05))
# for (j in 1:length(angdiv)) {
#   hh <- hist(angdiv[[j]] - 90, breaks = seq(40, 130, by=10), plot = F)
#   lines(hh$mids, hh$density, lty=1, bty='n',col= pal_9[j], lwd=3)
# }


# Fig.2M, ED Fig.3D, OD vs PD, scatter + hist, normalized by column edge distance --------------------

ii <- cc_b ==19 
PD_hex <- PD_hex_b[ii, ]
OD_hex <- OD_hex_b[ii, ]
# PDOD_T4b <- cbind(sqrt(rowSums((PD_hex[,1:2] - PD_hex[,3:4])^2)), sqrt(rowSums((OD_hex[,1:2] - OD_hex[,3:4])^2)) ) /(2+sqrt(2))
PDOD_T4b <- cbind(sqrt(rowSums((PD_hex[,1:2] - PD_hex[,3:4])^2)) /(2+sqrt(2)),
                  sqrt(rowSums((OD_hex[,1:2] - OD_hex[,3:4])^2)) /(4+sqrt(2)) ) 

ii <- cc_d ==19 
PD_hex <- PD_hex_d[ii, ]
OD_hex <- OD_hex_d[ii, ]
# PDOD_T4d <- cbind(sqrt(rowSums((PD_hex[,1:2] - PD_hex[,3:4])^2)), sqrt(rowSums((OD_hex[,1:2] - OD_hex[,3:4])^2)) ) /(4+sqrt(2))
PDOD_T4d <- cbind(sqrt(rowSums((PD_hex[,1:2] - PD_hex[,3:4])^2)) /(4+sqrt(2)),
                  sqrt(rowSums((OD_hex[,1:2] - OD_hex[,3:4])^2)) /(2+sqrt(2)) ) 

# scatter plot
df <- rbind(cbind(PDOD_T4b, 2), cbind(PDOD_T4d, 4) )
df <- as.data.frame(df)
colnames(df) <- c('PD','OD','gp')
df$gp <- factor(df$gp)
plt <- ggplot(df) +
  geom_point(aes(x = PD, y = OD, colour = gp), size = 3, alpha=0.8) +
  scale_color_manual(values= pal_T4[c(2,4)],labels = c("T4b", "T4d"), guide= guide_legend(title="dist[norm]"), na.value="gray") +
  ylab("OD") + xlab("PD") +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid.minor = element_blank() )+
  scale_x_continuous(limits = c(0.4, 1.4), breaks = seq(0,1.4,by=0.2), labels = seq(0,1.4,by=0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 1.4), breaks = seq(0,1.4,by=0.2), labels = seq(0,1.4,by=0.2), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(0.2, 1.8), breaks = seq(0.2,1.8,by=0.2), labels = seq(0.2,1.8,by=0.2), expand = c(0, 0)) + 
  labs(title = paste("OD vs PD_med", "_separate_unit", sep = "")) +
   coord_fixed(ratio=1)
windows(width = 8, height = 8)
# pdf("OD_PD_hex_sepa_norm.pdf",width = 8, height = 8)
# plt
ggMarginal(plt, margins = "both", size = 4, type = "density", binwidth= 0.2, boundary= 0.4, groupColour = TRUE, groupFill = F, lwd=2)
# dev.off()

wilcox.test(PDOD_T4b[,1], PDOD_T4d[,1])
t.test(PDOD_T4b[,1], PDOD_T4d[,1])

sd(PDOD_T4b[,1]) / mean(PDOD_T4b[,1])
sd(PDOD_T4b[,2]) / mean(PDOD_T4b[,2])
sd(PDOD_T4d[,1]) / mean(PDOD_T4d[,1])
sd(PDOD_T4d[,2]) / mean(PDOD_T4d[,2])

# # -- hist of PD
# windows(width = 8, height = 4)
# # pdf(paste("PD_normDist.pdf", sep = ""), width = 8, height = 4)
# hh <- hist(PDOD_T4b[,1], breaks = seq(0.4, 1.4, by=0.1), plot = F)
# plot(hh$mids, hh$density, type='l', bty='n', col=pal_T4[2], lwd=3, xlim = c(0.4, 1.4), xaxt='n',yaxt='n',
#      xlab ="", ylab='',main = '')
# hh <- hist(PDOD_T4d[,1], breaks = seq(0.4, 1.4, by=0.1), plot = F)
# lines(hh$mids, hh$density, lty=1, bty='n',col=scales::alpha(pal_T4[4], 0.7), lwd=3)
# axis(1, at = seq(0.4, 1.4, by =0.2), labels = seq(0.4, 1.4, by =0.2), cex.axis = 1.5 )
# axis(2, at = c(0, 10/length(hh$mids) ), labels = c('0', '1/n'), cex.axis = 1.5 )
# # dev.off()


# --  density
windows(width = 8, height = 4)
# pdf(paste("PD_normDist.pdf", sep = ""), width = 8, height = 4)
dd <- density(PDOD_T4b[,1], from= 0.4, to= 1.4, bw='SJ')
plot(dd$x, dd$y, type='l', bty='n', col=pal_T4[2], lwd=3, xlim = c(0.4, 1.4), xaxt='n',yaxt='n',
     xlab ="", ylab='',main = '')
dd <- density(PDOD_T4d[,1], from= 0.4, to= 1.4, bw='SJ')
points(dd$x, dd$y, type='l', bty='n',col=pal_T4[4], lwd=3)
axis(1, at =  seq(0.4, 1.4, by =0.2), labels =  seq(0.4, 1.4, by =0.2), cex.axis = 1.5 )
axis(2, at = c(0, 1/512/diff(dd$x[1:2]) ), labels = c('0', '1'), cex.axis = 1.5 )
# dev.off()

t.test(PDOD_T4b[,1], PDOD_T4d[,1])

# -- density with ggplot
dd <- density(PDOD_T4b[,1], from= 0.4, to= 1.4)

df <- data.frame(rbind(cbind(PDOD_T4b[,1], 2), cbind(PDOD_T4d[,1], 4)))
colnames(df) <- c('amp','gp')
df$gp <- as.factor(df$gp)

mu <- plyr::ddply(df, "gp", summarise, grp.mean=mean(amp))

windows(width = 8, height = 4)
# pdf(paste("PD_normDist.pdf", sep = ""), width = 8, height = 4)
ggplot() +
  geom_density(data= df, aes(amp, color=gp), position= 'identity', bw= 'SJ', lwd=3) +
  scale_color_manual(values= pal_T4[c(2,4)],labels = c("T4b", "T4d"), guide= guide_legend(title="dist[norm]")) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=gp), linetype="dashed") +
  scale_x_continuous(limits = c(0.4, 1.4), breaks = seq(0.4,1.4,by=0.2), labels = seq(0.4,1.4,by=0.2), expand = c(0, 0.05)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1), labels = c(0, 1), expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "black"))
# dev.off()

t.test(PDOD_T4b[,1], PDOD_T4d[,1])
t.test(PDOD_T4b[,2], PDOD_T4d[,2])

# ED Fig.3C, OD vs PD, scatter in um ------------------------------------------

LL <- 2 # choose type
ii <- cc_b ==19 

v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix() 
v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
v6 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
v7 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()

PDOD_um_T4b <- cbind( sqrt(rowSums((v1 - v0)^2)) /1000, sqrt(rowSums((v7 - v6)^2)) /1000 ) #in um
PDOD_um_T4b <- PDOD_um_T4b[ii, ]


LL <- 4 # choose type
ii <- cc_d ==19 

v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix() 
v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
v6 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
v7 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()

PDOD_um_T4d <- cbind( sqrt(rowSums((v1 - v0)^2)) /1000, sqrt(rowSums((v7 - v6)^2)) /1000 ) #in um
PDOD_um_T4d <- PDOD_um_T4d[ii, ]


# scatter plot
df <- rbind(cbind(PDOD_um_T4b, 2), cbind(PDOD_um_T4d, 4) )
df <- as.data.frame(df)
colnames(df) <- c('PD','OD','gp')
df$gp <- factor(df$gp)
plt <- ggplot(df) +
  geom_point(aes(x = PD, y = OD, colour = gp), size = 3, alpha=0.8) +
  scale_color_manual(values= pal_T4[c(2,4)],labels = c("T4b", "T4d"), guide= guide_legend(title="dist[norm]"), na.value="gray") +
  ylab("OD") + xlab("PD") +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid.minor = element_blank() )+
  scale_x_continuous(limits = c(5, 25), breaks = seq(0,25,by=5), labels = seq(0,25,by=5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(5, 25), breaks = seq(0,25,by=5), labels = seq(0,25,by=5), expand = c(0, 0)) +
  labs(title = paste("OD vs PD_med", "_um", sep = "")) +
  coord_fixed(ratio=1)
windows(width = 8, height = 8)
# pdf("OD_PD_um.pdf",width = 8, height = 8)
# plt
ggMarginal(plt, margins = "both", size = 4, type = "density", groupColour = TRUE, groupFill = F, lwd=2)
# dev.off()


# COV
sd(PDOD_um_T4b[,1]) / mean(PDOD_um_T4b[,1])
sd(PDOD_um_T4b[,2]) / mean(PDOD_um_T4b[,2])
sd(PDOD_um_T4d[,1]) / mean(PDOD_um_T4d[,1])
sd(PDOD_um_T4d[,2]) / mean(PDOD_um_T4d[,2])


# # spatial profile, T4b and T4d in med, use T4 gallery, Fig.2  ---------------------------
# 
# PDOD_T4b <- PDOD_T4(LL=2)
# PDOD_T4d <- PDOD_T4(LL=4)
# 
# df <- rbind(cbind(PDOD_T4b[,5:6], 2),
#             cbind(PDOD_T4d[,5:6], 4) )
# df <- as.data.frame(df)
# colnames(df) <- c('PD','OD','gp')
# df$gp <- factor(df$gp)
# 
# 
# plt <- ggplot(df) +
#   geom_point(aes(x = PD, y = OD, colour = gp), size = 3, alpha=0.8) +
#   scale_color_manual(values= pal_T4[c(2,4)],labels = c("T4b", "T4d"), guide= guide_legend(title="dist[norm]"), na.value="gray") +
#   ylab("OD") + xlab("PD") +
#   theme_minimal() +
#   scale_x_continuous(limits = c(1, 3.5), breaks = c(1,2,3), labels = c(1,2,3), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(1, 3.5), breaks = c(1,2,3), labels = c(1,2,3), expand = c(0, 0)) + # set +y as above eq
#   labs(title = paste("OD vs PD_med", "_norm_unit", sep = "")) +
#   coord_fixed(ratio=1)
# windows(width = 8, height = 8)
# # pdf("OD_PD.pdf",width = 8, height = 10)
# # plt
# # ggMarginal(plt, margins = "y", size = 2, type = "boxplot", outlier.size =3, groupColour = TRUE, groupFill = TRUE)
# ggMarginal(plt, margins = "both", size = 4, type = "density", groupColour = TRUE, groupFill = F, lwd=2)
# # dev.off()

# ED Fig.2F, extension by TmY5a -------------------------------------------------------

# # DEBUG
# nopen3d()
# # plot3d(TmY5a_me, col='gray50')
# # plot3d(TmY5a[c(10,12,13,14)], col='gray50')
# points3d(tag_xyz[,1:3], col ='magenta', size = 10)
# points3d(tag_xyz[,4:6], col ='cyan', size = 10)
# points3d(tag_xyz[,7:9], col ='dark green', size = 10)
# 
# # points3d(Mi1_M10_xyz, size=10)
# points3d(med_xyz, size=7)
# # points3d(Mi1_M10_xyz[Mi1_ind_PR[[8]],], size = 15, col='blue')
# # points3d(Mi1_M10_xyz[Mi1_ind_PR[[7]],], size = 15, col='light blue')
# 
# points3d(lo_pred, size=8, col='blue')


# PLOT, lo with col
nopen3d()
shade3d(LO_msh, col='gray', alpha=0.1)

ii <- med_ixy[med_ixy[,3] == 0, 1]
points3d(lo_pred[match(ii, eyemap[,1]),] , size=10, col=pal_axes[1])
# points3d(lo_pred[ii, ] , size=10, col=pal_axes[1])
ii <- med_ixy[med_ixy[,2] == -1, 1]
points3d(lo_pred[match(ii, eyemap[,1]), ], size=10, col=pal_axes[2])
ii <- med_ixy[med_ixy[,2] == med_ixy[,3]-1, 1]
points3d(lo_pred[match(ii, eyemap[,1]), ], size=10, col=pal_axes[3])
ii <- med_ixy[med_ixy[,2] == -med_ixy[,3]-1, 1]
points3d(lo_pred[match(ii, eyemap[,1]), ], size=10, col=pal_axes[4])
points3d(lo_pred, size=9, col='gray')

rgl.viewpoint(fov=0,zoom=0.8, userMatrix= rotationMatrix(30/180*pi,1,0,0) %*%
                rotationMatrix(60/180*pi,0,1,0) %*%
                rotationMatrix(0/180*pi,1,0,0))
# scale bar
c(10000,0,0) %>% 
  rotate3d(30/180*pi,1,0,0) %>%
  rotate3d(60/180*pi,0,1,0) %>%
  rotate3d(0/180*pi,1,0,0) %>%
  rbind(c(0,0,0), .) %>%
  sweep(2, c(3e5, 2e5, 2.8e5), '+') %>%
  segments3d()

# rgl.snapshot("LO_col.png")


# TmY5a
# nopen3d()
# plot3d(TmY5a[[10]], col='black', lwd=3, soma = T)
# rgl.snapshot(filename = "TmY5a_10.png")