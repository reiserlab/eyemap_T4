
# Fig.2G, Strahler new ------------------------------------------------------------
# cf. eyemap_T4_v2.R
# align with meridian lines in visual space

nb_coord <- rbind(c(0,0),
                  c(-5/180*pi,0),
                  c(+5/180*pi,0),
                  c(0, -5/180*pi),
                  c(0, +5/180*pi) )

# choose type 
LL <- 2
# and one of these 4 examples in Fig.2G
j <- 139
# j <- 170
# j <- 76
# j <- 69

com_xyz_eye <- as.matrix(lens_type[[LL]][, c("comx","comy","comz")])
T4_utp <- cart2sph2tp(com_xyz_eye) 
# med
com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix()
v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
v3 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
v4 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()


ii <- sweep(ucl_rot_sm, 2, com_xyz_eye[j, ] )^2 %>% rowSums() %>% which.min()
ii2 <- c(nb_ind[nb_ind[ii,],]) %>% unique()

# ref nb
pt <- ucl_rot_sm[ii, ] %>% matrix(ncol=3)
tp <- cart2sphZ(pt)[,2:3]
rtp <- cbind(1, sweep(nb_coord,2,tp,'+'))
xyz_ref <- sph2cartZ(rtp)

# np to med
xyz_ref_med <- matrix(nrow = nrow(xyz_ref), ncol = 3) 
xyz_eval <- data.frame(mc.x = xyz_ref[,1], mc.y = xyz_ref[,2], mc.z = xyz_ref[,3])
for (k in 1:3) {
  npdata <- data.frame(mc = ucl_rot_sm[ii2,], ec = med_xyz[ii2,k])
  bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
  model_np <- npreg(bw)
  xyz_ref_med[,k] <- predict(model_np, newdata = xyz_eval)
}

# pc
zz <- cross3D((xyz_ref_med[2,]-xyz_ref_med[1,]), (xyz_ref_med[5,]-xyz_ref_med[1,])) #pointing inwards
pc <- prcomp(xyz_ref_med)
if (pc$rotation[,3] %*% zz < 0) {
  pc$rotation <- -pc$rotation
}
if (c(cross3D(pc$rotation[,3], pc$rotation[,1])) %*% pc$rotation[,2] < 0) {
  pc$rotation[,2] <- -pc$rotation[,2]
}

xyz_ref_med_pc <- sweep(xyz_ref_med, 2, xyz_ref_med[1,]) %*% pc$rotation
y_ref <- diff(xyz_ref_med_pc[c(3,2),1:2])
ang <- acos(y_ref %*% c(0,1) / sqrt(sum(y_ref^2)))
ang <- if_else(y_ref[1] < 0, -ang, ang)
rotM <- matrix(c(cos(ang), sin(ang),0,
                 -sin(ang), cos(ang), 0,
                 0, 0, 1), ncol = 3, byrow = T)
xyz_ref_med_pc_rot <- xyz_ref_med_pc %*% rotM

# T4
tar <- T4_dend[[LL]][[j]]
ind_D = match(tar$tags$`dendrite start` , tar$d$PointNo)
df_D <-  tar$d[ind_D,]
root_xyz <- tar$d[ind_D, c("X","Y","Z")]

# - find the subtree with root = dendrite start
targ <- as.ngraph(tar)
ii_root <- ind_D
# subtree and Strahler order
sub_points <- igraph::graph.dfs(targ, root = ii_root, unreachable=FALSE, neimode='out')$order
subtree <- subset(tar, sub_points) 
subtree_g <- as.ngraph(subtree, weights = T)
subtree_so <- strahler_order(subtree) # Strahler order
max(subtree_so$segments)

PD <- rbind(com_xyz[j,], v1[j,], v0[j,]) #center and arrow
OD <- rbind(v3[j,], v4[j,]) #OD

# transform to pc coord, xyz_ref_med[1,] as origin
subtree_pc <- subtree
subtree_pc$d[, c("X","Y","Z")] <- sweep(xyzmatrix(subtree$d), 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM
root_xyz_pc <- as.numeric(root_xyz - xyz_ref_med[1,]) %*% pc$rotation %*% rotM 
PD_pc <- sweep(PD, 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM #center and arrow
OD_pc <- sweep(OD, 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM #center and arrow
nb_pc <- sweep(med_xyz[ii2,], 2, xyz_ref_med[1,]) %*% pc$rotation %*% rotM # nb col


# PLOT
windows(width = 12, height = 12)
# pdf(paste("T4_eg_SO_", LL, "_", j, ".pdf", sep = ''), width = 12, height = 12)
plot(nb_pc[1:7, 1:2], col = 'gray20', pch = 1, cex = 25, asp = 1, lwd=1,
     xlab = '', ylab = '', xaxt="n", yaxt="n", xlim = c(-10e3, 10e3), ylim = c(-10e3, 10e3) )
# par("mai"=c(0.2,0.2,0.2,0.2), "omi"=c(0.2,0.2,0.2,0.2))

for (k in 1:max(subtree_so$points)) {
  # k <- 1
  pt_p <- subtree_pc$d[,"PointNo"] %in% unique(subtree_pc$d[subtree_so$points == k,"Parent"])
  pt_so <- pt_p | subtree_so$points==k
  if (sum(pt_so) > 1) {
    plot(subset(subtree_pc, pt_so), col=pal_so[k], add = T, lwd = 1.5*k + 1, WithNodes = F)
  }
}
arrows(PD_pc[2,1], PD_pc[2,2], PD_pc[3,1], PD_pc[3,2], lwd= 3,)
arrows(OD_pc[1,1], OD_pc[1,2], OD_pc[2,1], OD_pc[2,2], lwd= 3, code = 3, angle = 90, length = 0.2)

# dev.off()

# PLOT, vec for each seg
kk_seg <- which(subtree_so$segments == k)
for (m in kk_seg) {
  ii <- subtree_pc$SegList[[m]]
  if (length(ii) >= 2 & runif(1,0,1) > 0) {
    head <- subtree_pc$d[head(ii,1), c("X","Y","Z")] %>% as.matrix()
    tail <- subtree_pc$d[tail(ii,1), c("X","Y","Z")] %>% as.matrix()
    # arrows(head[1], head[2], tail[1], tail[2], lwd= 1.5, length= 0.05)
    segments(head[1], head[2], tail[1], tail[2], lwd= 1.5)
  }
}

# ED Fig.3A, SO histogram ------------------------------------------------------------

# annulus  area
# 5 deg elev, solid angles
dtheta <- 10
N <- 180 / dtheta
theta <- seq(0,pi, length.out = (N+1))
h <- diff(cos(theta)) %>% abs()
A <- 2*pi*1*h
h_norm <- round(h / min(h), 2)


# - histo wrt PD
LL <- 3
ssn <- dplyr::bind_rows(seg_summ_type[[LL]])
# ssn <- seg_summ_type[[LL]][[1]]

df <- matrix(ncol = 3, nrow = 0)
for (ii_h in 1:4) {
  if (nrow(ssn[ssn$so == ii_h,]) > 0 ) {
    # hh <- hist(ssn[ssn$so == ii_h,]$ang, breaks = seq(0,180,dtheta), plot = F)
    # hh <- hist(ssn[ssn$so == ii_h,]$angV, breaks = seq(0,180,dtheta), plot = F)
    # hh <- hist(ssn[ssn$so == ii_h,]$angH, breaks = seq(0,180,dtheta), plot = F)
    hh <- hist(ssn[ssn$so == ii_h,]$angSO4, breaks = seq(0,180,dtheta), plot = F)
    # df <- rbind(df, cbind(hh$mids, hh$density, ii_h))
    df <- rbind(df, cbind(hh$mids, hh$density/h_norm, ii_h))
  }
}
df <- as.data.frame(df)
colnames(df) <- c('ang', 'freq', 'so')

windows(width = 5, height = 16)
# pdf(file = paste("SO_hist_angSO4_areaNorm_T4", letters[LL], ".pdf", sep = ''), width = 5, height = 16)
ggplot(df, aes(x=ang, y=freq)) +
  geom_bar(stat='identity') +
  scale_x_continuous(breaks = seq(0,180,45), labels = paste(seq(0, 180, by=45),"°", sep = '') ) +
  # coord_polar(theta = "x") +
  facet_grid(rows = vars(so) ) +
  theme_minimal() +
  theme( axis.title = element_blank(),
         axis.text.y = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.y = element_blank()) +
  # labs(title = paste("T4", letters[LL], "_SO_histo", sep = '') )
  # labs(title = paste("T4", letters[LL], "_SO_histo_areaNorm_vaxis", sep = '') )
  labs(title = paste("T4", letters[LL], "_SO_histo_areaNorm_haxis", sep = '') )
# dev.off()

# 
# seg_summ <- seg_summ_type[[LL]]
# for (j in 1:length(neu)) {
#   ssn <- seg_summ[[j]] # seg_summ
#   df <- matrix(ncol = 3, nrow = 0)
#   for (ii_h in 1:4) {
#     if (nrow(ssn[ssn$so == ii_h,]) > 0 ) {
#       hh <- hist(ssn[ssn$so == ii_h,]$ang, breaks = seq(0,180,10), plot = F)
#       df <- rbind(df, cbind(hh$mids, hh$density, ii_h))
#     }
#   }
#   df <- data.frame(df)
#   colnames(df) <- c('ang', 'freq', 'so')
#   plt <- ggplot(df, aes(x=ang, y=freq)) +
#     geom_bar(stat='identity', width = 5) +
#     # coord_cartesian(ylim = c(0, 1)) +
#     scale_x_continuous(breaks = seq(0,180,45), labels = seq(0,180,45)) +
#     facet_grid(rows = vars(so)) +
#     theme_minimal()
#   # do.call("grid.arrange", plt)
#   print(plt)
# }


# # hotspot 3d ----------------------------------------------------------------
# 
# T4_overlap_pos <- c('M','AD','AV','LV')
# 
# hs_ind <- c(0, cumsum(T4_overlap_N))
# 
# hs_dir <- matrix(ncol = 6, nrow = nrow(anno_overlap))
# hs_dir_lens <- matrix(ncol = 6, nrow = nrow(anno_overlap))
# hs_com <- matrix(ncol = 3, nrow = nrow(anno_overlap))
# hs_com_lens <- matrix(ncol = 3, nrow = nrow(anno_overlap))
# for (j in 1:nrow(anno_overlap)) {
#   anno <- anno_overlap[j,]
#   ii <- regexpr("T4", anno$name) #start index
#   type_ii <- match(substring(anno$name, ii+2, ii+2), letters)
#   i2 <- match(anno$skid, anno_T4_dend[[type_ii]]$skid) # T4 index
#   hs_dir[j,] <- unlist(dir_type[[type_ii]][i2, c('rsx0','rsy0','rsz0','rsxd','rsyd','rszd')])
#   hs_dir_lens[j,] <- unlist(lens_type[[type_ii]][i2, c('x0','y0','z0','xd','yd','zd')])
#   hs_com[j,] <- unlist(dir_type[[type_ii]][i2, c('comx','comy','comz')])
#   hs_com_lens[j,] <- unlist(lens_type[[type_ii]][i2, c('comx','comy','comz')])
# }
# 
# # lens
# nopen3d()
# points3d(ucl_rot_sm, size = 5)
# # points3d(hs_com_lens, size =7, col='red')
# for (j in 1:nrow(anno_overlap)) {
#   type_ii <- anno_overlap$ti[j]
#   arrow3d(hs_dir_lens[j,1:3], hs_dir_lens[j,4:6], theta = pi / 9, n = 4, col =pal_T4[type_ii], type = "rotation", lit=F)
#   # arrow3d(RF_lens[j,1:3], RF_lens[j,4:6], theta = pi / 18, n = 4, col =pal_T4[type_ii], type = "rotation", lit=F)
# }
# 
# # PLOT, hot spot on eye
# nopen3d()
# par3d('windowRect' = c(100,100,1000,1000))
# points3d(ucl_rot_sm, size=7)
# # points3d(ucl_rot_sm[rownames(ucl_rot_sm) %in% vaxis_gen(-1), ], size=12, col=pal_axes[3])
# # points3d(ucl_rot_sm[c(ind_Up_ucl, ind_Down_ucl), ], size=12, col=pal_axes[4])
# # points3d(hs_com_lens, size = 10, col='red')
# spheres3d(0,0,0,0.99, col='grey80', alpha=1, lit=F)
# arrow3d(c(1,0,0), c(1.6, 0, 0), theta = pi / 18, n = 8, col = "gray30", type = "rotation")
# text3d(c(1.7, -0.2, -0.3), texts = 'front', cex=3)
# arrow3d(c(0,0,1), c(0, 0, 1.2), theta = pi / 9, n = 8, col = "gray30", type = "rotation")
# text3d(c(0, -0.2, 1.3), texts = 'up', cex=3)
# planes3d(0,1,0, 0, alpha = 0.2)
# planes3d(0,0,1, 0, alpha = 0.2)
# for (j in 1:4) {
#   xyz <- hs_com_lens[hs_ind[j]:hs_ind[j+1],] %>% colMeans()
#   xyz <- xyz / sqrt(sum(xyz^2))
#   points3d(matrix(xyz,ncol = 3), size = 30, col='light green')
# }
# rgl.viewpoint(fov=0,zoom=0.7, userMatrix= rotationMatrix(-80/180*pi,1,0,0) %*% rotationMatrix(-47/180*pi,0,0,1))
# 
# # rgl.snapshot(filename = paste("hotspot_eye.png", sep = ''))


# ED Fig.6C, hotspot Mollweide --------------------------------------------------

# cf. "hot spot  Mollweide" in Figure_2_T4.R
hs_Mo <- matrix(ncol = 3, nrow = 4)
for (j in 1:4) {
  xyz <- hs_com_lens[hs_ind[j]:hs_ind[j+1],] %>% colMeans()
  xyz <- xyz / sqrt(sum(xyz^2))
  ii <- sweep(ucl_rot_sm, 2, xyz)^2 %>% rowSums() %>% sqrt() %>% which.min()
  hs_Mo[j,] <- ucl_rot_sm[ii,]
}
colnames(hs_Mo) <- c('x','y','z')
hs_Mo %<>% as_tibble() %>%  
  mutate(y = -y) %>%
  mutate(theta = acos(z)) %>%
  mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
  mutate(t = theta / pi * 180, p = phi/pi*180) %>%  
  mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
  as.data.frame()
hs_Mo <- Mollweide(hs_Mo[,c('t', 'p')])
colnames(hs_Mo) <- c('xM','yM')

df <- as.data.frame(ucl_rot_Mo)
plt <- plt_Mo +
  geom_point(data=df, aes(x = xM, y = yM), colour = 'gray', size=1) +
  # geom_contour(data=df_pred_right, aes(x=xM,y=yM, z=Z), breaks=breaks_contour,
  #              colour='black', alpha=0.9, lwd=0.5) +
  geom_point(data = df[ind_xy[ind_xy[,2] == clp,1], ], aes(xM, yM), colour = pal_axes[2], size=2) +
  geom_point(data = df[ind_xy[ind_xy[,3] == clq,1], ], aes(xM, yM), colour = pal_axes[1], size=2) +
  # geom_point(data = df[c(ind_Up_ucl, ind_Down_ucl), ], aes(xM, yM), colour = pal_axes[4], size=2) +
  geom_point(data = df[haxis_gen(clh, ind_xy), ], aes(xM, yM), colour = pal_axes[4], size=2) +
  geom_point(data = df[vaxis_gen(clv, ind_xy), ], aes(xM, yM), colour = pal_axes[3], size=2) +
  geom_point(data = as.data.frame(hs_Mo), aes(xM, yM), colour = 'gray', shape = 19, size=25, alpha=0.5 ) +
  labs(title = "Mollweide 3 axes with eq, hot spots")
windows(width = 17, height = 9)
# pdf("Mollweide_ioa_hex_4axes_hotspots.pdf", width = 8.5, height = 4.5)
plt
# dev.off()


# coord transform, NEW, use ref grid --------------------------------------

# add home col -> too variable to use
anno_overlap_hc <- anno_overlap
anno_overlap_hc$hc <- NA
for (j in 1:nrow(hs_com)) {
  anno_overlap_hc$hc[j] <- sweep(med_xyz, 2, hs_com[j,])^2 %>% rowSums() %>% which.min
}


# - np
# 2024
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
# change view, now viewed from proximal side, ie, looking outwards
refhex_2[, 1] <- -refhex_2[, 1]

# rescaling factor, med -> refhex
rs <- mean(nb_dist_med, na.rm =T) / mean(sqrt(rowSums(refhex_1[-1,]^2)))
rs <- 1

# hotspot dir
com_xyz <- hs_com
v0 <- hs_dir[, 1:3]
v1 <- hs_dir[, 4:6]

# - use np
v3 <- matrix(ncol = 2, nrow = nrow(v0))  #base of vec of T4_dir (orientation not PD).
v4 <- matrix(ncol = 2, nrow = nrow(v1)) #head of vec
v5 <- matrix(ncol = 2, nrow = nrow(com_xyz)) # for com
cc <- c()
for (j in 1:nrow(v0)) {
  ii <- sweep(med_xyz, 2, com_xyz[j,])^2 %>% rowSums() %>% which.min()
  ii7 <- nb_ind[ii,]
  # ii19 <- c(ii7, nb_ind[ii7[3],2:3],nb_ind[ii7[6],c(3,6,4)],nb_ind[ii7[4],4:5], nb_ind[ii7[5],5],nb_ind[ii7[7],c(5,7,2)],nb_ind[ii7[2],2]) #2021
  ii19 <- c(ii7, nb_ind[ii7[2],2:3],nb_ind[ii7[3],c(3,4)],nb_ind[ii7[4],4:5], nb_ind[ii7[5],c(5,6)],nb_ind[ii7[6],c(6,7)],nb_ind[ii7[7],c(7,2)]) #2024
  cc <- c(cc, sum(!is.na(ii19)))

  iinna <- !is.na(ii19)
  nb_xyz <- med_xyz[ii19[iinna], ]

  nb_pca <- prcomp(nb_xyz)
  nb_xyz_pc <- sweep(nb_xyz, 2,nb_pca$center) %*% nb_pca$rotation
  colnames(nb_xyz_pc) <- c('x','y','z')
  nb_xy <- nb_xyz_pc[, 1:2] / rs

  vec_xyz <- rbind(v0[j,], v1[j,])
  vec_xyz_pc <- sweep(vec_xyz, 2,nb_pca$center) %*% nb_pca$rotation
  colnames(vec_xyz_pc) <- c('x','y','z')
  vec_xy <- vec_xyz_pc[, 1:2] / rs

  com_xyz_pc <- (com_xyz[j,] - nb_pca$center) %*% nb_pca$rotation
  colnames(com_xyz_pc) <- c('x','y','z')
  com_xy <- com_xyz_pc[1:2] / rs

  # -- np
  xyz_vtail_eval <- data.frame(mc.x = vec_xy[1,1], mc.y = vec_xy[1,2])
  xyz_vhead_eval <- data.frame(mc.x = vec_xy[2,1], mc.y = vec_xy[2,2])
  xyz_com_eval <- data.frame(mc.x = com_xy[1], mc.y = com_xy[2])
  for (k in 1:2) {
    npdata <- data.frame(mc = nb_xy, ec = refhex_2[iinna,k])
    bw <- npregbw(formula= ec~mc.x+mc.y, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
    model_np <- npreg(bw)
    v3[j,k] <- predict(model_np, newdata = xyz_vtail_eval)
    v4[j,k] <- predict(model_np, newdata = xyz_vhead_eval)
    v5[j,k] <- predict(model_np, newdata = xyz_com_eval)
  }
}


# PLOT
df0 <- as.data.frame(refhex_2)
colnames(df0) <- c('x','y')

df <- rbind(cbind(v4,1), cbind(v3,2), cbind(v5,3) ) %>% as.data.frame() #reverse arrow st. it's PD
colnames(df) <- c('x','y', 'gp')
df$gp <- factor(df$gp)

df_arrow <- as.data.frame(rbind(cbind(v4, v5),
                                cbind(v5, v3)) )
colnames(df_arrow) <- c('x','y','xend','yend')

# -- 2d
windows(width = 7, height = 10)
plt <- ggplot() +
  geom_segment(data = df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour='gray',size =0.2) +
  geom_point(data=df0, aes(x,y), size=30, pch=1) +
  geom_point(data=df, aes(x=x,y=y, col=gp), size=3, na.rm = T ) +
  # geom_point(data=df, aes(x=x,y=y, col=gp, shape=gp), size=3, na.rm = T ) +
  # geom_density_2d()
  # scale_color_manual(values = pal_T4, guide= guide_legend(title="")) +
  scale_color_manual(values = c(scales::alpha('violetred', 0.7), scales::alpha('cyan3', 0.7), scales::alpha('gray20', 0.5)), guide= guide_legend(title="")) +
  # scale_shape_manual(values=c(2, 1, 3))+
  ylab("y") +   xlab("x") +
  theme_minimal() +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3,3), labels = seq(-3,3), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(-2.5, 2.5), expand = c(0, 0), breaks = seq(-2,2), labels = seq(-2,2)) + # set +y as above eq
  scale_y_continuous(limits = c(-4, 4), expand = c(0, 0), breaks = seq(-5,5)) + # set +y as above eq
  labs(title = paste("T4", letters[LL],  sep = "")) +
  theme(axis.title = element_blank(), panel.grid.major = element_line(colour = "gray70"), panel.grid.minor =  element_blank() )+
  coord_fixed(ratio=1)
# pdf(paste("hex_ref_", letters[LL], '_N.pdf', sep = ""), width = 7, height = 10)
plt

# - calculate angles
PD_ang <- data.frame(v3 - v4)
colnames(PD_ang) <- c('x','y')
PD_ang$ang <- apply(PD_ang, 1, function(x) atan2(x[2],x[1])/pi*180)


# cont. ED Fig.6D, hot spot stats, NEW -------------------------------------------------------------------

# angles betw 
TT <- matrix(c(1,2, 3,4,
               1,3, 2,4,
               1,4, 2,3), ncol = 2, byrow = T)
df <- matrix(ncol = 3, nrow = 0)
for (ii_t in 1:nrow(TT)) {
  t1 <- TT[ii_t,1]
  t2 <- TT[ii_t,2]
  for (j in 1:4) {
    anno <- anno_overlap[(hs_ind[j]+1):hs_ind[j+1],]
    ang <- PD_ang$ang[(hs_ind[j]+1):hs_ind[j+1]]
    ang_ls <- list()
    for (k in 1:4) {
      ang_ls[[k]] <- ang[anno$ti == k]
    }
    angtt <- c()
    for (m in 1:length(ang_ls[[t1]])) {
      for (n in 1:length(ang_ls[[t2]])) {
        angtt <- c(angtt, abs(ang_ls[[t1]][m] - ang_ls[[t2]][n]))
      }
    }
    
    df <- rbind(df, cbind(angtt, j, 10*ii_t))
  }
}
df <- as.data.frame(df)
colnames(df) <- c('ang', 'hs', 'pair')
# df$ang <- if_else(df$ang > 180 & df$pair >= 30, 360 - df$ang, df$ang)
hs_ang <- df


# PLOT
ii <- hs_ang$pair <= 20 # ab and cd
df <- data.frame(x = hs_ang$pair[ii] + hs_ang$hs[ii]*2, ang = hs_ang$ang[ii])
# df$ang <- if_else(df$ang > 180 & hs_ang$pair >= 30, 360 - df$ang, df$ang)

windows(width = 12, height = 6)
# pdf(paste("hs_stat_anti", '.pdf', sep = ''), width = 12, height = 6)
ggplot(df, aes(factor(x), ang) ) + 
  # geom_jitter(aes(colour=neu),size = 4,height = 0, width = 0.1) +
  geom_jitter(colour= 'gray', size = 1,height = 0, width = 0.1) +
  # scale_colour_manual(values = col_4_65, guide=FALSE) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 1, position = position_nudge(x = 0.2, y = 0)) +
  # coord_cartesian(ylim = c(0, 180)) +
  # scale_y_continuous(breaks= seq(0,180,by=45)) +
  coord_cartesian(ylim = c(120, 240)) +
  scale_y_continuous(breaks= seq(120, 240,by=60), labels = paste0(seq(120, 240,by=60), "°")) +
  scale_x_discrete(labels=c("12"="M", "14"="AD", "16"="AV", "18"="LV",
                            "22"="M", "24"="AD", "26"="AV", "28"="LV")) +
  theme_minimal() +
  labs(title= paste("angles between anti T4", sep = ''),
       x='Hotspot position', y='Angle [deg]') 
# dev.off()


# aa <- hs_ang[hs_ang$hs != 1 & hs_ang$pair %in% c(30,40, 50,60), ]
# aa$pair[aa$pair %in% c(30,40)] <- 10
# aa$pair[aa$pair %in% c(50,60)] <- 100
# df <- aa
# df$group <- df$hs * df$pair
# 
# windows(width = 6, height = 3)
# # pdf(paste("hs_ang_90", '.pdf', sep = ''), width = 5, height = 3)
# ggplot(df, aes(factor(group), ang) ) + 
#   # geom_jitter(aes(colour=neu),size = 4,height = 0, width = 0.1) +
#   geom_jitter(colour= 'gray', size = 1,height = 0, width = 0.1) +
#   # scale_colour_manual(values = col_4_65, guide=FALSE) +
#   stat_summary(fun.min = function(z) { quantile(z,0.25) },
#                fun.max = function(z) { quantile(z,0.75) },
#                geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
#                colour = 'black', size = 1, width = .2 ) +
#   stat_summary(fun = median, geom = "point",  colour = 'black', size = 1, position = position_nudge(x = 0.2, y = 0)) +
#   coord_cartesian(ylim = c(0, 180)) +
#   scale_y_continuous(breaks= seq(0,180,by=45)) +
#   scale_x_discrete(labels=c("20"="AD", "30"="AV", "40"="LV", "200"="AD", "300"="AV", "400"="LV")) +
#   theme_minimal() +
#   labs(title= paste("angles between T4s", sep = ''), x='Hotspot', y='Angle [degree]')
# # dev.off()

# # ellipsoid on eye --------------------------------------------------------
# 
# # # - hot spot
# # load("data/neu_T4_overlap.RData")
# 
# hs_ind <- c(0, cumsum(T4_overlap_N))
# 
# hs_ell <- matrix(ncol = 3, nrow = nrow(anno_overlap))
# hs_ell_eye <- matrix(ncol = 2, nrow = nrow(anno_overlap))
# for (j in 1:nrow(anno_overlap)) {
#   anno <- anno_overlap[j,]
#   ii <- regexpr("T4", anno$name) #start index
#   type_ii <- match(substring(anno$name, ii+2, ii+2), letters)
#   i2 <- match(anno$skid, anno_T4_dend[[type_ii]]$skid) # T4 index
#   hs_ell[j,] <- unlist(ell_type[[type_ii]][i2, c('ea','eb','ec')])
#   hs_ell_eye[j,] <- unlist(ell_type[[type_ii]][i2, c('eye_ea','eye_eb')])
# }
# colnames(hs_ell) <- c('ea','eb','ec')
# colnames(hs_ell_eye) <- c('ea','eb')
# 
# hs_ell[anno_overlap$ti < 2.5, 1]
# 
# hs_ell[,1] / hs_ell[,2] * hs_ell_eye[,1] / hs_ell_eye[,2] 
# 
# eeratio <- hs_ell[,1] / hs_ell[,2] / (hs_ell_eye[,1] / hs_ell_eye[,2] )
# eeratio <- if_else(anno_overlap$ti > 2.5, 1/eeratio, eeratio)
# mean(eeratio)
# sd(eeratio)
# 
# dev.new()
# # plot(eeratio, col='blue', pch=16)
# hist(eeratio)
# 
# # - all T4
# T4_N <- sapply(ell_type, nrow)
# # TODO
# 
# # transform st avg T4a -> x-axis, both in med and on eye, plot in orca -------
# 
# # # ###  2 choose 1
# # # overlap_dir <- hs_dir; pt_pc3 <- c(350000,250000,250000) # med
# # overlap_dir <- hs_dir_lens; pt_pc3 <- c(0,0,0)  #lens
# # 
# # avg_sd_hs <- list() # 4 hot spots
# # vrn_hs <- list()
# # ang_hs <- list()
# # for (j in 1:4) {
# #   v0 <- overlap_dir[(hs_ind[j]+1):hs_ind[j+1], 1:3]
# #   vd <- overlap_dir[(hs_ind[j]+1):hs_ind[j+1], 4:6]
# #   anno <- anno_overlap[(hs_ind[j]+1):hs_ind[j+1],]
# #   pc <- prcomp(rbind(v0,vd))
# #   if (pc$rotation[,3] %*% (colMeans(rbind(v0,vd)) - pt_pc3) < 0 ) {
# #     pc$rotation[,3] <- - pc$rotation[,3]
# #   }
# #   if (t(cross3D(pc$rotation[,3], pc$rotation[,1])) %*% pc$rotation[,2] < 0 ) {
# #     pc$rotation[,2] <- - pc$rotation[,2]
# #   }
# #   
# #   v0_pc <- sweep(v0, 2, pc$center) %*% pc$rotation %>% .[,1:2]
# #   vd_pc <- sweep(vd, 2, pc$center) %*% pc$rotation %>% .[,1:2]
# #   # rotate and normalize
# #   v1 <- (vd_pc - v0_pc)[anno$ti == 1,]
# #   amp <- sqrt(sum(colMeans(v1)^2))
# #   # ang1 <- apply((vd_pc - v0_pc)[anno$ti == 1,], 1, function(x) atan2(x[2], x[1]))
# #   # ang1 <- if_else(ang1 < 0, 2*pi+ang1, ang1)
# #   # angr <-  mean(ang1)
# #   v1n <- sweep(v1, 1, sqrt(rowSums(v1^2)), '/')
# #   angr <- atan2(colMeans(v1n)[2], colMeans(v1n)[1])
# #   angr <- angr + pi #T4a -> -x
# #   rot <- matrix(c(cos(-angr), sin(-angr),
# #                   -sin(-angr), cos(-angr)), ncol = 2, byrow = T)
# #   vrn <- (vd_pc - v0_pc) %*% rot / amp
# #   avg_sd <- matrix(ncol = 2, nrow = 4)
# #   ang_cell <- c()
# #   for (k in 1:4) {
# #     ang <- apply(vrn[anno$ti == k,], 1, function(x) atan2(x[2], x[1])) 
# #     if (k == 1) {
# #       ang <- if_else(ang < 0, 2*pi+ang, ang)
# #     }
# #     ang_cell <- c(ang_cell, ang)
# #     avg_sd[k,1] <- round(mean(ang) /pi*180, 1)
# #     avg_sd[k,2] <- round(sd(ang) /pi*180, 1) # /(n-1)
# #   }
# #   ang_cell[order(anno$ti)] <- ang_cell
# #   ang_hs[[j]] <- ang_cell
# #   avg_sd_hs[[j]] <- avg_sd
# #   vrn_hs[[j]] <- vrn
# # }
# # 
# # 
# # # PLOT polar plot
# # rmax <- sapply(vrn_hs, function(x) sqrt(rowSums(x^2))) %>% unlist %>% max()
# # rmax <- 2.5
# # fig <- list()
# # for (j in 1:4) {
# #   fig[[j]] <- plot_ly(type = 'scatterpolar', mode = "markers") %>%
# #     add_trace(mode = "markers",
# #               r = sqrt(rowSums(vrn_hs[[j]]^2)),
# #               theta = apply(vrn_hs[[j]], 1, function(x) atan2(x[2], x[1])) /pi*180,
# #               marker = list(color = pal_T4[anno_overlap[(hs_ind[j]+1):hs_ind[j+1],'ti']],
# #                             symbol = 'circle',
# #                             size = 10
# #               ),
# #               name = "T4",
# #               showlegend = FALSE
# #     ) 
# #   for (k in 1:4) {
# #     fig[[j]] <- fig[[j]] %>%
# #       add_trace(mode = 'lines', # add wedge
# #                 r = c(0,1,1,1),
# #                 theta = c(0,
# #                           avg_sd_hs[[j]][k,1] - avg_sd_hs[[j]][k,2],
# #                           avg_sd_hs[[j]][k,1],
# #                           avg_sd_hs[[j]][k,1] + avg_sd_hs[[j]][k,2]),
# #                 fill = 'toself',
# #                 fillcolor = adjustcolor(pal_T4[k], alpha.f = 0.3),
# #                 line = list(
# #                   width = 0
# #                 ),
# #                 name = paste("T4", letters[k], sep = '')
# #       ) 
# #   }
# #   fig[[j]] <- fig[[j]] %>%
# #     layout(
# #       title = list(
# #         text = paste("position:", T4_overlap_pos[j]),
# #         x = 0.1,
# #         y = 0.9
# #       ),
# #       font = list(
# #         family = 'Arial',
# #         size = 16,
# #         color = '#000'
# #       ),
# #       polar = list(
# #         radialaxis = list(
# #           visible = T,
# #           showline = F,
# #           range = c(0,rmax),
# #           tick0 = 0,
# #           dtick = 1,
# #           showticklabels = F,
# #           ticks = "",
# #           angle = 45
# #         ),
# #         angularaxis = list(
# #           thetaunit= "degree",
# #           dtick = 45,
# #           showline = F,
# #           size = 16
# #         )
# #       ),
# #       margin = list(l = 35, r = 15, t = 25, b = 10),
# #       showlegend = T
# #     )
# #   
# #   # SAVE
# #   orca(fig[[j]], file = paste("overlap_med_", T4_overlap_pos[j], ".png", sep=''))
# #   # orca(fig[[j]], file = paste("overlap_eye_", T4_overlap_pos[j], ".png", sep='') )
# #   # orca(fig[[j]],  file = paste("overlap_eye_", T4_overlap_pos[j], ".pdf", sep=''))
# # }