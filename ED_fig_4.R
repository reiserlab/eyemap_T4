# see Fig_3_uCT.R for ED Fig.4D, ED Fig.4E, ED Fig.4G

# ED Fig.4B, normal VS tip-lens ------------------------------------------------------

nucl_eyemap <- nucl[eyemap[,2], ] 
lens_eyemap <- lens[eyemap[,2], ] 

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

# ED Fig.4H, ED Fig.4J, curvature and dia, -----------------------------------

# plot sphere

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
view3d(fov=0,zoom=0.8, userMatrix= rotationMatrix(-90/180*pi,1,0,0) %*% rotationMatrix(-90/180*pi,0,0,1))
spheres3d(lens_left[,1],lens_left[,2],lens_left[,3], 5, lit=F)
# rgl.snapshot("curvature_2eye.png")

ind_roc <- ind_roc_right # [ind, p,v,q,h,sph]
ind_dia <- ind_dia_right
lens_ixy_hex <- ind_xy_right
unit_hex <- t(cbind(c(cos(30/180*pi), sin(30/180*pi)), c(-cos(30/180*pi), sin(30/180*pi))))
lens_ixy_hex[,2:3] <- lens_ixy_hex[,2:3] %*% unit_hex


# - PLOT cur
quantile(ind_roc[,c(3,5,6)], c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T)
bin_lim <- c(120, 320)

## ## choose a direction, 3-v, 5-h, 6-sph
ii <- 3
df <- cbind(lens_ixy_hex[,2:3], ind_roc[,ii]) %>% as.data.frame()
colnames(df) <- c('x','y','r')

plt <- ggplot() + 
  geom_point(data=df, aes(x = x, y = y, colour = r), size = 2) +
  scale_color_gradient(
    low = "gray95", high= "magenta", limits=bin_lim, oob=scales::squish,# trans='log',
    breaks=bin_lim, labels=bin_lim, guide = guide_colorbar(title = "r[um]"), 
    na.value = 'yellow') +
  theme_void() +
  theme(legend.position = c(.9, .9) ) +
  coord_fixed(ratio = 1) +
  labs(title = "radius of curvature along v-axis")
# labs(title = "radius of curvature along h-axis")
# labs(title = "radius of curvature sph")
windows(width = 8, height = 8)
# pdf("roc_v.pdf")
# pdf("roc_h.pdf")
# pdf("roc_sph.pdf")
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

# PLOT circle fit
dev.new()
# pdf('cir_fit_eg.pdf', width = 6, height = 6)
plot(xyz_pc[,2], xyz_pc[,1], pch=16, 
     asp = 1, xlim = c(-100, 650), ylim = c(-400,400))
points(xyz5_pc[,2], xyz5_pc[,1], pch=16, cex=1.5, col='red')
draw.circle(cir[[1]][1], cir[[1]][2], r_fit, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
# dev.off()

# - PLOT diameters
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
  coord_fixed(ratio = 1) +
  labs(title = "diameter")
windows(width = 8, height = 8)
# pdf("dia_6.pdf")
plt
# dev.off()