# ED Fig.2E, chiasm vs v-axis --------------------------------------------------------

ind_chi <- match(na.omit(match(anno_chi$skid, anno_Mi1$skid)), rownames(med_xyz))

# PLOT
nopen3d()
par3d('windowRect' = c(100,100,1300,1300))
ii <- na.omit(eyemap[match(vaxis_gen(0), eyemap[,2]),1])
points3d(Mi1_M10_xyz[ii,], size=10, col=pal_axes[3])
points3d(Mi1_M10_xyz[-ii,], size=10, col='gray')

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

# ED Fig.2F, extension by TmY5a -------------------------------------------------------

# PLOT, lo with col
nopen3d()
shade3d(LO_msh, col='gray', alpha=0.1)

ii <- med_ixy[med_ixy[,3] == 0, 1]
points3d(lo_pred[match(ii, eyemap[,1]),] , size=10, col=pal_axes[1])
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


ii <- c(76, 69) # eg skids

nopen3d()
par3d('windowRect' = c(100,100,1300,1300))
plot3d(T4d[ii], lwd=4, col="black")
plot3d(T4d[-ii], lwd=2, 
       col=rep(c("#C69F9F","#B3C69F","#9FC6C6","#B29FC6","gray"),35)[1:length(T4d[-ii])]
)

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