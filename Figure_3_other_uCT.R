# re-run Figure_0 afterwards

# ED Fig.4F, ED Fig.4G, ioa and skewness, 2023, redo nb_ind  -------------------

# fn <- c(
#   # '20230926',
#   # '20231107',
#   # '20240206'
#   '20240701',
#   '20231107',
#   '20240530'
# )

fn <- c(
  '20240701',
  '20231107',
  '20240530'
)

for (f in 1:length(fn)) {
  # load(paste0("../microCT/2023_eyemap/", fn[f], '.RData'))
  load(paste0("data/microCT/", fn[f], ".RData"))

  lens_left <- lens[ind_left_lens,]
  lens_mir <- lens_left
  lens_mir[,2] <- -lens_mir[,2]
  lens_2eye <- lens
  lens_right <- lens[!ind_left_lens,]
  ind_Up_lens_right <- na.omit(match(i_match[ind_Up], rownames(lens_right)))
  ind_Up_lens_left <- na.omit(match(i_match[ind_Up], rownames(lens_left)))
  ind_Down_lens_right <- na.omit(match(i_match[ind_Down], rownames(lens_right)))
  ind_Down_lens_left <- na.omit(match(i_match[ind_Down], rownames(lens_left)))
  # cone_lens <- i_match[!ind_left_cone]
  # re-order ucl_rot rows and re-define rownames because mapping is betw Mi1 and lens
  ucl_rot_left <- ucl_rot_sm[order(i_match),][ind_left_lens,]
  colnames(ucl_rot_left) <- c('x','y','z')
  ucl_rot_right <- ucl_rot_sm[order(i_match),][!ind_left_lens,]
  colnames(ucl_rot_right) <- c('x','y','z')
  
  for (lr in c(1,2)) {
    if (lr == 1) {
      pt <- lens_left
      ind_Up <- ind_Up_lens_left
      ind_Down <- ind_Down_lens_left
      ucl_rot <- ucl_rot_left
    } else {
      pt <- lens_right
      ind_Up <- ind_Up_lens_right
      ind_Down <- ind_Down_lens_right
      ucl_rot <- ucl_rot_right
    }
    
    hexcoord <- reghex(pt)
    ind_nb <- hexcoord[[1]]
    ind_xy <- hexcoord[[2]]
    
    
    # # PLOT hex grid
    # nopen3d()
    # points3d(pt)
    # text3d(pt[ind_xy[,1],], texts = paste(ind_xy[,2],ind_xy[,3],sep=','), adj=1.3)
    # points3d(pt[ind_xy[ind_xy[,2] %% 3 == 0,1],], size=10, col='magenta')
    # 
    # # PLOT regular grid
    # lens_ixy_hex <- ind_xy
    # unit_hex <- t(cbind(c(cos(30/180*pi), sin(30/180*pi)), c(-cos(30/180*pi), sin(30/180*pi))))
    # lens_ixy_hex[,2:3] <- lens_ixy_hex[,2:3] %*% unit_hex
    # 
    # windows(width = 8, height = 8)
    # # pdf("axes_lens_hex.pdf")
    # # plot(lens_ixy_hex[, 2:3], xaxt="n", yaxt="n", pch=16, cex=1.25, col='gray')
    # plot(lens_ixy_hex[, 2:3], pch=16, cex=1.25, col='gray')
    # points(lens_ixy_hex[ind_xy[,2] == 0,  2:3], pch=16, col=pal_axes[2], cex=1.5)
    # points(lens_ixy_hex[ind_xy[,3] == 0, 2:3], pch=16, col=pal_axes[1], cex=1.5)
    # points(lens_ixy_hex[ind_xy[,2] == ind_xy[,3] , 2:3], pch=16, col=pal_axes[3], cex=1.5)
    # points(lens_ixy_hex[ind_xy[,2] == -ind_xy[,3], 2:3], pch=16, col=pal_axes[4], cex=1.5)
    # # text(ind_xy[,2:3], labels = ind_xy[,1], cex = 0.7, adj = -0.5)
    # title(paste("hex", fn[f], c('left','right')[lr], 'N=', nrow(ind_xy), sep='_'))
    # # dev.off()
    
    
    # - Mollweide
    ucl_rot_tp <- cart2sphZ(ucl_rot)[,2:3] %>%
      as_tibble() %>%
      mutate(theta = theta/pi*180, phi = phi/pi*180) %>%
      mutate(phi = if_else(phi > 180, phi - 360, phi)) %>%
      mutate(phi = - phi) %>% #inside out
      as.matrix()
    ucl_rot_Mo <- Mollweide(ucl_rot_tp)
    colnames(ucl_rot_Mo) <- c('xM','yM')
    rownames(ucl_rot_Mo) <- rownames(ucl_rot)
    
    
    # - use ind_nb to determine nb and ioa, assume  Euclidean dist ~ arclength
    nb_dist_ucl <- matrix(nrow = nrow(ucl_rot), ncol = 7)
    for (j in 1:nrow(ucl_rot)) {
      nb_dist_ucl[ind_nb[j,1],] <- c(ind_nb[j,1],
                                     acos(ucl_rot[ind_nb[j,-1], ] %*% ucl_rot[ind_nb[j,1],]) )
    }
    ioa_hex <- cbind(ucl_rot_Mo, rowMeans(nb_dist_ucl[,2:7], na.rm = T)/pi*180)
    ioa_hex <- as.data.frame(ioa_hex)
    colnames(ioa_hex) <- c('xM','yM','ioa')
    
    
    # - loessm ioa
    # -- grid
    bkgd_chull <- rbind(bkgd_mer_ww, bkgd_mer_ee[seq(nrow(bkgd_mer_ee),1,by=-1),])
    colnames(bkgd_chull) <- c('xM','yM')
    grid_M <- expand.grid(xM = seq(-sqrt(8), sqrt(8), length.out = 100),
                          yM = seq(-sqrt(2), sqrt(2), length.out = 50) )
    ii_inpoly <- sp::point.in.polygon(grid_M[,1], grid_M[,2], bkgd_chull[,1], bkgd_chull[,2])
    
    # -- loess
    fit_loess <- loess(ioa ~ xM * yM, data = ioa_hex, degree = 2, span = 0.1,
                       control = loess.control(surface = "direct"))
    pred_loess <- predict(fit_loess, grid_M, se = T)
    df_pred <- grid_M
    df_pred$Z <- melt(pred_loess$fit)$value
    df_pred <- df_pred[ii_inpoly == 1, ]
    
    xy_ashape <- ashape(ucl_rot_Mo[, 1:2] +
                          matrix(runif(nrow(ucl_rot)*2, 1e-9, 2e-9), ncol=2), alpha = 0.3)
    xy_edge <- xy_ashape$edges[,1:6]
    xy_grid_ahull <- mkpoly(xy_edge)[[1]][,3:4] # hull edge points
    xy_grid_ahull <- rbind(xy_grid_ahull, xy_grid_ahull[1,])
    ii_inpoly <- sp::point.in.polygon(grid_M[ii_inpoly == 1,1], grid_M[ii_inpoly == 1,2], xy_grid_ahull[,1], xy_grid_ahull[,2])
    df_pred$Z[!ii_inpoly] <- 20 # outside set to large angles
    
    
    
    # -- plot
    # range(df_pred$Z)
    breaks_contour <- c(0, 4, 5, 6)
    # getPalette <- colorRampPalette(c('gray30', 'gray90'))
    # pal_contour <- getPalette(length(breaks_contour)+1)
    
    df <- as.data.frame(ioa_hex)
    plt <- plt_Mo +
      geom_contour(data=df_pred, aes(x=xM,y=yM, z=Z), breaks=breaks_contour,
                   colour='gray20', alpha=0.9, lwd=0.5) +
      geom_point(data = df, aes(xM, yM), colour =pal_lr[lr], size=1, alpha=1) +
      # geom_text(data = df, aes(xM, yM, label=round(ioa,1)), size= 3) +
      # geom_point(data = df[ind_xy[ind_xy[,2] == 0, 1], ], aes(xM, yM), colour = pal_axes[2], size=2) +
      # geom_point(data = df[ind_xy[ind_xy[,3] == 0, 1], ], aes(xM, yM), colour = pal_axes[1], size=2) +
      # geom_point(data = df[ind_xy[ind_xy[,2] == ind_xy[,3], 1], ], aes(xM, yM), colour = pal_axes[3], size=2) +
      # geom_point(data = df[ind_xy[ind_xy[,2] == -ind_xy[,3], 1], ], aes(xM, yM), colour = pal_axes[4], size=2) +
      theme(legend.position = c(.9, .9) ) +
      labs(title = paste("ioa ", c("left","right")[lr], " eye_", fn[f], sep=''))
    # windows(width = 16, height = 8); plt
    pdf(paste("ioa_", c("left","right")[lr], " eye_", fn[f], ".pdf", sep=''), width = 8.5, height = 4.5)
    print(plt)
    dev.off()
    
    
    # - skewness
    nb_ang_ucl <- matrix(ncol = 2, nrow = nrow(ucl_rot))
    hex_hvf_eye <- matrix(ncol = 3, nrow = nrow(ucl_rot))
    for (j in 1:nrow(ind_nb)) {
      if (sum(complete.cases(ind_nb[j,])) == 7) {
        # eye
        bt <- ucl_rot[ind_nb[j,3], ] - ucl_rot[ind_nb[j,6], ]
        bf <- colMeans(ucl_rot[ind_nb[j,4:5],]) - colMeans(ucl_rot[ind_nb[j,c(2,7)],])
        hex_hvf_eye[ind_nb[j,1],] <- bf
        ang <- acos(bt %*% bf / sqrt(sum(bt^2)) / sqrt(sum(bf^2)) ) /pi*180 
        nb_ang_ucl[ind_nb[j,1],] <- c(j, ang)
      }
    }
    
    # # - histo
    # range(na.omit(nb_ang_ucl[,2]))
    # dev.new()
    # # pdf("hist_latt_ang_eye.pdf")
    # hh <- hist(nb_ang_ucl[,2], breaks = seq(0, 180, by=5), plot = F)
    # plot(hh$mids, hh$counts, type='l', bty='n', xlim = c(30, 150), xaxt='n',yaxt='n', xlab ="angle [deg]", ylab='counts')
    # axis(1, at = seq(30, 150, by =30), labels =  paste(seq(30, 150, by =30), "°", sep = '') )
    # # dev.off()
    
    df_pos <- data.frame(ucl_rot_Mo)
    df_pos$quan <- nb_ang_ucl[,2]
    # df_pos$quan <- NA
    # df_pos$quan[1:5] <- nb_ang_ucl[1:5,2]
    quantile(na.omit(df_pos$quan), c(0, 0.01, 0.1, 0.5, 0.9, 0.99, 1))
    
    
    # -- Mollweide
    # rg <- c(45, 90, 135)
    rg <- c(60, 90, 120)
    plt <- plt_Mo + 
      # geom_path(data = cmer, aes(x=xM, y=yM), colour = pal_axes[3], lwd=1) +
      geom_point(data=df_pos, aes(x = xM, y = yM, colour = quan), size = 2) +
      scale_color_gradientn(colours = pal_heat1, values = scales::rescale(rg), limits=range(rg), oob=scales::squish,
                            breaks= rg, labels= rg, guide = guide_colorbar(title = "skewness") ) +
      # geom_text(data=df_pos, aes(xM, yM, label=round(quan,1)), size=3, nudge_x = 0.03) +
      # geom_text(data=df_pos, aes(xM, yM, label= seq(1,nrow(df_pos))), size=3, nudge_y = 0.02) +
      theme(legend.position = c(.9, .9) ) +
      labs(title = "skew angle")
    # windows(width = 9, height = 6)
    pdf(paste0("skew_eye_", fn[f], "_", c('left','right')[lr], ".pdf"), width = 8.5, height = 4.5)
    print(plt)
    dev.off()
    
  }
}



# ED Fig.4E, ioa along the equator -------------------------------------------

# for (f in 1:length(fn)) {
  # for (f in 1:2) {
ioa_Mo <- list()
  f <- 1
  # load(paste0("../microCT/2023_eyemap/", fn[f], '.RData'))
  load(paste0("data/microCT/", fn[f], ".RData"))
  
  lens_2eye <- lens
  lens_left <- lens[ind_left_lens,]
  lens_right <- lens[!ind_left_lens,]
  ind_Up_lens_right <- na.omit(match(i_match[ind_Up], rownames(lens_right)))
  ind_Up_lens_left <- na.omit(match(i_match[ind_Up], rownames(lens_left)))
  ind_Down_lens_right <- na.omit(match(i_match[ind_Down], rownames(lens_right)))
  ind_Down_lens_left <- na.omit(match(i_match[ind_Down], rownames(lens_left)))
  # cone_lens <- i_match[!ind_left_cone]
  # re-order ucl_rot rows and re-define rownames because mapping is betw Mi1 and lens
  ucl_rot_2eye <- ucl_rot
  ucl_rot_left <- ucl_rot[order(i_match),][ind_left_lens,]
  rownames(ucl_rot_left) <- seq(1, nrow(ucl_rot_left))
  colnames(ucl_rot_left) <- c('x','y','z')
  ucl_rot_right <- ucl_rot[order(i_match),][!ind_left_lens,]
  rownames(ucl_rot_right) <- seq(1, nrow(ucl_rot_right))
  colnames(ucl_rot_right) <- c('x','y','z')

  
  for (lr in c(1,2)) {
    if (lr == 1) {
      pt <- lens_left
      ind_Up <- ind_Up_lens_left
      ind_Down <- ind_Down_lens_left
      ucl_rot <- ucl_rot_left
      hexcoord <- reghex(pt, lefteye = T)
    } else {
      pt <- lens_right
      ind_Up <- ind_Up_lens_right
      ind_Down <- ind_Down_lens_right
      ucl_rot <- ucl_rot_right
      hexcoord <- reghex(pt, lefteye = F)
    }
    
    # hexcoord <- reghex(pt)
    ind_nb <- hexcoord[[1]]
    ind_xy <- hexcoord[[2]]
  
    # - Mollweide
    ucl_rot_Mo <- ucl_rot
    colnames(ucl_rot_Mo) <- c('x','y','z')
    ucl_rot_Mo %<>% as_tibble() %>%
      mutate(y = -y) %>%
      mutate(theta = acos(z)) %>%
      mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
      mutate(t = theta / pi * 180, p = phi/pi*180) %>%
      mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
      as.data.frame()
    ucl_rot_Mo <- Mollweide(ucl_rot_Mo[,c('t', 'p')])
    colnames(ucl_rot_Mo) <- c('xM','yM')
    rownames(ucl_rot_Mo) <- rownames(ucl_rot)
    
    
    # - use ind_nb to determine nb and ioa, assume  Euclidean dist ~ arclength
    nb_dist_ucl <- matrix(nrow = nrow(ucl_rot), ncol = 7)
    for (j in 1:nrow(ucl_rot)) {
      # nb_dist_ucl[ind_nb[j,1],] <- sweep(ucl_rot[ind_nb[j,],], 2, ucl_rot[ind_nb[j,1],])^2 %>% rowSums() %>% sqrt()
      nb_dist_ucl[ind_nb[j,1],] <- c(
        ind_nb[j,1],acos(ucl_rot[ind_nb[j,-1], ] %*% ucl_rot[ind_nb[j,1],])
      )
    }
    ioa_hex <- cbind(ucl_rot_Mo, rowMeans(nb_dist_ucl[,2:7], na.rm = T)/pi*180)
    ioa_hex <- as.data.frame(ioa_hex)
    colnames(ioa_hex) <- c('xM','yM','ioa')
    
    ioa_Mo[[lr]] <- ioa_hex
  }
    
  
delev <- 15 # +/-15 deg in Mollweide
# Mollweide(rbind(c(70, 0), c(90,0), c(110,0)))
dMx <-  Mollweide(rbind(c(90,0), c(90,5)))[2,1] # 5 deg in Mollweide
eqbw <-  Mollweide(rbind(c(90,0), c(90,delev)))[2,1] 

# # avg over 10 degree azim
# xazim <- seq(0, sqrt(8), by = dMx/5)
# xazim <- c(rev(xazim[-1]), xazim)

lr <- 1
xya <- ioa_Mo[[lr]]
xya <- xya[abs(xya$yM) < eqbw, ] 
# xioa <- c()
# for (j in 1:length(xazim)) {
#   xioa[j] <- mean(xya$ioa[xya$xM > xazim[j]- dMx & xya$xM < xazim[j]+ dMx])
# }
# xioa[is.na(xioa)] <- 0
# xioa_l <- 
df_l <- data.frame(x=xya$xM, y=xya$ioa)

lr <- 2
xya <- ioa_Mo[[lr]]
xya <- xya[abs(xya$yM) < eqbw, ]
# xioa <- c()
# for (j in 1:length(xazim)) {
#   xioa[j] <- mean(xya$ioa[xya$xM > xazim[j]- dMx & xya$xM < xazim[j]+ dMx])
# }
# xioa[is.na(xioa)] <- 0
# xioa_r <- xioa
df_r <- data.frame(x=xya$xM, y=xya$ioa)

# df_l <- data.frame(x=seq(-180,180), y=xioa_l)
# df_r <- data.frame(x=seq(-180,180), y=xioa_r)

xazim <- seq(0, sqrt(8), by = dMx*12)
xazim <- c(-rev(xazim[-1]), xazim)

windows(width =7.5, height = 6)
# pdf(paste0("ioa_equ_", delev, "deg.pdf"))
ggplot() +
  geom_point(data=df_l, aes(x,y),size=2, col=pal_lr[1]) +
  geom_smooth(data=df_l, aes(x,y), se=F, col='blue', lwd=1) +
  geom_point(data=df_r, aes(x,y),size=2, col=pal_lr[2]) +
  geom_smooth(data=df_r, aes(x,y), se=F, col='black', lwd=1) +
  # geom_line(data=df_r, aes(x,y),lwd=2, col=pal_lr[2]) +
  scale_x_continuous(limits = c(-sqrt(8), sqrt(8)), breaks =xazim, labels =paste0(xazim/dMx/6*30,"°"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(3, 9), breaks = seq(3,10), labels =paste0(seq(3,10),"°"), expand = c(0, 0)) +
  theme_minimal() +
  ylab("inter-ommatidium angle [deg]") +
  xlab("azimuth [deg]") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14) ) +
  # theme(legend.position="none", panel.background = element_blank()) +
  labs(title = paste("ioa +/-", delev, fn)) +
  coord_fixed(ratio = 0.5)
# dev.off()
