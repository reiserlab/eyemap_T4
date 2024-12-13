

# calculation, 2022  ---------------------------------------------------------------------

# tb <- read_excel("data_Eyal/H2_2021/movGrtTabwEZ_ThetaCorr.xlsx", sheet = 1, col_names = T ) %>% data.frame()
# tb_v <- Eyal_arena_2023(tb)

# A <- (4+10/32 + 4+20/32)/2  #radius
# 
# Dx <- (5+6/32) - A  # area center wrt to fly
# Dy <- (4+20/32) - A
# Dz <- - (1+ 13/32)
# 
# # LED size
# DL <- A * 1.25/180*pi
# 
# # response amp factor
# resfac <- 20
# 
# # - loop
# tb_v <- matrix(ncol = 2+2+2, nrow = nrow(tb)) # [elev azim]
# 
# for (j in 1:nrow(tb)) {
#   # - set para
#   t1 <- tb$meanThetaSpk[j] #mean theta spike
#   t2 <- tb$meanThetaSub[j] #subthr
#   # t12 <- (t1+t2)/2
#   tt <- c(0, t1, t2) #base, head(spk), head(sub)
#   # tt <- -tt #clockwise
#   DSI1 <- tb$meanDSISpk[j] # spike
#   DSI2 <- tb$meanDSISub[j] #subthr
#   DD <- c(0, DSI1, DSI2)
#   nx <- tb$stimPosX[j] #azim angle of LED
#   ny <- tb$stimPosY[j] #elev of LED
#   arot <- tb$arenaAng[j] /180*pi #arena rotation ang
#   afly <- tb$headAng[j] /180*pi #eqator below holder
# 
#   if (!is.na(afly)) {
#     for (k in 1:3) {
#       # in arena 2D coord
#       # ?? use t1/2/t12 * (DL or 0)
#       # DL * resfac is arbitrary
#       # 96 pixel * 1.25 = 120 deg per side
#       ax <- (120 - nx*1.25) /180*pi - cos(tt[k])* DL /A * DD[k] * resfac * (tt[k] != 0) #angle from midline
#       ay <- DL * ny + sin(tt[k]) * DL  * DD[k] * resfac * (tt[k] != 0)   #dist from bottom
# 
#       # go to lab coord via a translation
#       # origin at arena center, xy-plane is level
#       x <- A * cos(ax)
#       y <- A * sin(ax)
#       z <- ay - (4 + 24/32)/2
# 
# 
#       xyz <- matrix(c(x,y,z), ncol=3)
# 
#       # arena rotation,
#       # arot is angle below ground, use -arot
#       xyz <- matrix(c(cos(arot), 0, -sin(-arot),
#                       0, 1, 0,
#                       sin(-arot), 0, cos(arot)), ncol=3, byrow = T) %*% t(xyz)
#       xyz <- t(xyz)
# 
#       # go to fly's coord
#       xyz <- xyz + c(Dx, Dy, Dz)
# 
#       # fly rotation
#       # afly is angle below ground, use afly since it's ref that's rotating
#       xyz <- matrix(c(cos(afly), 0, -sin(afly),
#                       0, 1, 0,
#                       sin(afly), 0, cos(afly)), ncol = 3, byrow = T) %*% t(xyz)
#       xyz <- t(xyz)
# 
#       # to spherical coord
#       thetaphi <- cart2sphZ(xyz) /pi*180
# 
#       tb_v[j, (k-1)*2+(1:2)] <- round(c(90 - thetaphi[2], thetaphi[3]), 2)
#     }
#   }
# }
# 
# # tb_v [theta(base), phi(base), spk, spk, sub, sub]

# write.csv(tb_v, file = "data/H2_tuning.csv", row.names = F)
# save(tb_v, tb, file = 'data/H2_tuning.rda')


# plot  --------------------------------------------------------------------

load('data/H2_tuning.rda' )

# pick a cell
j <- 7
ii <- tb$ï..cellNum == j

thetaphi <- tb_v[ii, 1:2] #spk

thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
thetaphi <- thetaphi / 180*pi
thetaphi <- cbind(1, thetaphi)
xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])

pxyz <- xyz

thetaphi <- tb_v[ii,3:4]

thetaphi[, 1] <- 90 - thetaphi[, 1]
thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
thetaphi <- thetaphi / 180*pi
thetaphi <- cbind(1, thetaphi)
xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])

dxyz <- xyz

df <- cbind(cart2Mercator(pxyz), cart2Mercator(dxyz))
df <- df[!is.na(df[,1]),]
colnames(df) <- c('x','y','xend','yend')
df_arrow <- as.data.frame(df)

# mercator
plt <- plt_Mer +
  # geom_path(data = cmer_Merc, aes(x=x, y=y), colour = pal_axes[3], lwd=1) +
  # geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour=pal_T4[LL],size =1, arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend), colour='black',size =1) +
  scale_x_continuous(limits = c(-180,180)/180*pi,breaks = seq(-180,180,by=45)/180*pi, labels = seq(-180,180,by=45),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(-75,75)/180*pi,breaks = log(tan(pi/4 + seq(-75,75,by=15)/180*pi/2)), labels = seq(-75,75,by=15),  expand = c(0, 0)) +
  # scale_x_continuous(limits = c(-20,80)/180*pi,breaks = seq(-20,80,by=10)/180*pi, labels = seq(-20,80,by=10),  expand = c(0, 0)) +
  # scale_y_continuous(limits = c(-40,60)/180*pi,breaks = log(tan(pi/4 + seq(-40,60,by=20)/180*pi/2)), labels = seq(-40,60,by=20),  expand = c(0, 0)) +
  # annotate("text", x = df_arrow$x+0.02, y = df_arrow$y+0.02, label = 1:nrow(df_arrow), size = 3) + 
  labs(title = "")
windows(width = 14, height = 6.5)
plt
# ggsave(paste("H2_Mercator_cell_", j, ".pdf",sep=''), width = 14, height = 6.5)


# plot all 5 cells --------------------------------------------------------

load('data/H2_tuning.rda' )

df <- matrix(ncol = 5, nrow = 0)
for (j in 1:7) {

  ii <- tb$ï..cellNum == j
  thetaphi <- tb_v[ii, 1:2] #base
  
  if (sum(is.na(thetaphi[,1])) < nrow(thetaphi)) {
    thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    pxyz <- xyz
    
    # thetaphi <- tb_v[ii,3:4] #spk
    thetaphi <- tb_v[ii,5:6] #sub
    
    thetaphi[, 1] <- 90 - thetaphi[, 1]
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    dxyz <- xyz
    
    df_tmp <- cbind(cart2Mercator(pxyz), cart2Mercator(dxyz), j)
  }
  
  df <- rbind(df, df_tmp[!is.na(df_tmp[,1]),])
}

colnames(df) <- c('x','y','xend','yend','cell')
df_arrow <- as.data.frame(df)
df_arrow$cell <- as.factor(df_arrow$cell)

# mercator
plt <- plt_Mer +
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=cell),size =1)+
  # geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=cell),size =1, arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
  scale_colour_manual(values = brewer.pal(7,"Set1"), breaks = c("1", "2", "3","4","5","6","7")) +
  scale_x_continuous(limits = c(-45,90)/180*pi,breaks = seq(-45,90,by=45)/180*pi, labels = paste0(seq(-45,90,by=45),"°"),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(-15,75)/180*pi,breaks = log(tan(pi/4 + seq(-15,75,by=15)/180*pi/2)), labels = paste0(seq(-15,75,by=15),"°"),  expand = c(0, 0)) +
  # annotate("text", x = df_arrow$x+0.02, y = df_arrow$y+0.02, label = 1:nrow(df_arrow), size = 3) + 
  labs(title = "")
windows(width = 9, height = 8)
plt
# ggsave(paste("H2_Mercator_spk.pdf",sep=''), width = 9, height = 8)
# ggsave(paste("H2_Mercator_sub.pdf",sep=''), width = 9, height = 8)


# process new data 2023 edge -----------------------------------------------------------

# edges
# tb <- read_excel("data_Eyal/H2_2023/movEdgeTabwPY.xlsx", sheet = 1, col_names = T ) %>% data.frame()
# tb <- read_excel("data_Eyal/H2_2023/movEdgeTabwPY_ThetaCorr.xlsx", sheet = 1, col_names = T ) %>% data.frame()
# tb <- read_excel("data_Eyal/H2_2023/movEdgeTabwPY_ThetaCorr_fullResp.xlsx", sheet = 1, col_names = T ) %>% data.frame()
tb <- read_excel("data/movEdgeTabwPY_ThetaCorr_fullResp.xlsx", sheet = 1, col_names = T ) %>% data.frame()
names(tb)[names(tb) == "headAngPitch"] <- "headAng"

tb_v <- Eyal_arena_2023(tb)

# # tb_v [theta(base), phi(base), spk, spk, sub, sub, dark=0/bright=1]
# save(tb_v, tb, file = 'data/H2_tuning_2023.rda')

# plot all 7 cells, edge ---------------------------------------------------

load('data/H2_tuning_2023.rda' )

df <- matrix(ncol = 5, nrow = 0)
for (j in 1:7) {
  for (k in c(0,1)) {
    
    ii <- tb$cellNum == j & tb$edgeVal == k
    
    thetaphi <- tb_v[ii, 1:2] #base
    
    if (sum(is.na(thetaphi[,1])) < nrow(thetaphi)) {
      thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
      thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
      thetaphi <- thetaphi / 180*pi
      thetaphi <- cbind(1, thetaphi)
      xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
      xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
      
      pxyz <- xyz
      
      ### ###
      thetaphi <- tb_v[ii,3:4] #spk
      ### ###
      # thetaphi <- tb_v[ii,5:6] #sub
      
      thetaphi[, 1] <- 90 - thetaphi[, 1]
      thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
      thetaphi <- thetaphi / 180*pi
      thetaphi <- cbind(1, thetaphi)
      xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
      xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
      
      dxyz <- xyz
      
      df_tmp <- cbind(cart2Mercator(pxyz), cart2Mercator(dxyz), j, k)
    }
    df <- rbind(df, df_tmp[!is.na(df_tmp[,1]),])
  }
}
colnames(df) <- c('x','y','xend','yend','cell', 'edge')


## ### dark
df_arrow <- as.data.frame(df[df$edge == 0,])
# ### ### bright
# df_arrow <- as.data.frame(df[df$edge == 1,])
# ### ### birght x1
# df_arrow <- as.data.frame(df[df$edge == 1 & df$cell == 1,])

df_arrow$cell <- as.factor(df_arrow$cell)
# mercator
plt <- plt_Mer +
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=cell),size =1)+
  # geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=cell),size =1, arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
  scale_colour_manual(values = brewer.pal(7,"Set1"), breaks = c("1", "2", "3","4","5","6","7")) +
  scale_x_continuous(limits = c(-45,90)/180*pi,breaks = seq(-45,90,by=45)/180*pi, labels = paste0(seq(-45,90,by=45),"°"),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(-15,75)/180*pi,breaks = log(tan(pi/4 + seq(-15,75,by=15)/180*pi/2)), labels = paste0(seq(-15,75,by=15),"°"),  expand = c(0, 0)) +
  # annotate("text", x = df_arrow$x+0.02, y = df_arrow$y+0.02, label = 1:nrow(df_arrow), size = 3) +
  labs(title = "edge - spk - dark")
windows(width = 9, height = 8)
plt
# ggsave(paste("H2_Mercator_edge_spk_bright.pdf",sep=''), width = 9, height = 8)
# ggsave(paste("H2_Mercator_edge_spk_dark.pdf",sep=''), width = 9, height = 8)



# process new data 2023 gratings ----------------------------------------------------------------

# tb <- read_excel("data_Eyal/H2_2023/movGrtTabwPY_ThetaCorr.xlsx", sheet = 1, col_names = T ) %>% data.frame()
tb <- read_excel("data/movGrtTabwPY_ThetaCorr.xlsx", sheet = 1, col_names = T ) %>% data.frame()
names(tb)[names(tb) == "headAngPitch"] <- "headAng"

tb_v <- Eyal_arena_2023(tb)

# save(tb_v, tb, file = 'data/H2_tuning_2023_grating.rda')

# plot all 7 cells, grating ----------------------------------------------------

load('data/H2_tuning_2023_grating.rda' )

df <- matrix(ncol = 5, nrow = 0)
for (j in 1:7) {
  
  ii <- tb$cellNum == j 
  
  thetaphi <- tb_v[ii, 1:2] #base
  
  if (sum(is.na(thetaphi[,1])) < nrow(thetaphi)) {
    thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    pxyz <- xyz
    
    thetaphi <- tb_v[ii,3:4] #spk
    # thetaphi <- tb_v[ii,5:6] #sub
    
    thetaphi[, 1] <- 90 - thetaphi[, 1]
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    dxyz <- xyz
    
    df_tmp <- cbind(cart2Mercator(pxyz), cart2Mercator(dxyz), j)
  }
  df <- rbind(df, df_tmp[!is.na(df_tmp[,1]),])
}
colnames(df) <- c('x','y','xend','yend','cell')

df_arrow <- as.data.frame(df)
df_arrow$cell <- as.factor(df_arrow$cell)

# mercator
plt <- plt_Mer +
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=cell),size =1)+
  # geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=cell),size =1, arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
  scale_colour_manual(values = brewer.pal(7,"Set1"), breaks = c("1", "2", "3","4","5","6","7")) +
  scale_x_continuous(limits = c(-45,90)/180*pi,breaks = seq(-45,90,by=45)/180*pi, labels = paste0(seq(-45,90,by=45),"°"),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(-15,75)/180*pi,breaks = log(tan(pi/4 + seq(-15,75,by=15)/180*pi/2)), labels = paste0(seq(-15,75,by=15),"°"),  expand = c(0, 0)) +
  # annotate("text", x = df_arrow$x+0.02, y = df_arrow$y+0.02, label = 1:nrow(df_arrow), size = 3) + 
  labs(title = "")
windows(width = 9, height = 8)
plt
# ggsave(paste("H2_Mercator_spk.pdf",sep=''), width = 9, height = 8)
# ggsave(paste("H2_Mercator_sub.pdf",sep=''), width = 9, height = 8)


# combine all grating data---------------------------------------------------

load('data/H2_tuning.rda' )

df <- matrix(ncol = 5, nrow = 0)
for (j in 1:7) {
  
  ii <- tb$cellNum == j
  thetaphi <- tb_v[ii, 1:2, drop=F] #base
  
  if (sum(is.na(thetaphi[,1])) < nrow(thetaphi)) {
    thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    pxyz <- xyz
    
    thetaphi <- tb_v[ii,3:4] #spk
    # thetaphi <- tb_v[ii,5:6] #sub
    
    thetaphi[, 1] <- 90 - thetaphi[, 1]
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    dxyz <- xyz
    
    df_tmp <- cbind(cart2Mercator(pxyz), cart2Mercator(dxyz), j)
  }
  
  df <- rbind(df, df_tmp[!is.na(df_tmp[,1]),])
}

colnames(df) <- c('x','y','xend','yend','cell')
df5 <- df


# - new data
load('data/H2_tuning_2023_grating.rda' )

df <- matrix(ncol = 5, nrow = 0)
for (j in 1:7) {
  
  ii <- tb$cellNum == j 
  
  thetaphi <- tb_v[ii, 1:2] #base
  
  if (sum(is.na(thetaphi[,1])) < nrow(thetaphi)) {
    thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    pxyz <- xyz
    
    thetaphi <- tb_v[ii,3:4] #spk
    # thetaphi <- tb_v[ii,5:6] #sub
    
    thetaphi[, 1] <- 90 - thetaphi[, 1]
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    dxyz <- xyz
    
    df_tmp <- cbind(cart2Mercator(pxyz), cart2Mercator(dxyz), j)
  }
  df <- rbind(df, df_tmp[!is.na(df_tmp[,1]),])
}
colnames(df) <- c('x','y','xend','yend','cell')

df5$cell <- df5$cell + 7
df <- rbind(df, df5)

# PLOT
df_arrow <- as.data.frame(df)
df_arrow$cell <- as.factor(df_arrow$cell)

# mercator
plt <- plt_Mer +
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=cell),size =1)+
  # geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=cell),size =1, arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
  scale_colour_manual(values = c( brewer.pal(7,"Set1")[c(1,2,5,6,7)], rep('black',5)) ) +
  scale_x_continuous(limits = c(-45,90)/180*pi,breaks = seq(-45,90,by=45)/180*pi, labels = paste0(seq(-45,90,by=45),"°"),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(-15,75)/180*pi,breaks = log(tan(pi/4 + seq(-15,75,by=15)/180*pi/2)), labels = paste0(seq(-15,75,by=15),"°"),  expand = c(0, 0)) +
  # annotate("text", x = df_arrow$x+0.02, y = df_arrow$y+0.02, label = 1:nrow(df_arrow), size = 3) + 
  labs(title = "grating - spk")
windows(width = 9, height = 8)
plt


# ED Fig.1D, combine all grating data, avg ------------------------------------

load('data/H2_tuning.rda' )

# change from 87 to 96 for averaging
# tb$stimPosX[tb$stimPosX == 87] <- 96

# remove for avg
tb_v <- tb_v[!tb$stimPosX == 87,]
tb <- tb[!tb$stimPosX == 87,]
tb_v <- tb_v[!tb$stimPosY == 24,]
tb <- tb[!tb$stimPosY == 24,]

uxy <- unique(tb[, c('stimPosX', 'stimPosY')])

df <- matrix(ncol = 5, nrow = 0)
for (j in 1:nrow(uxy)) {
  ii <- tb$stimPosX == uxy$stimPosX[j] &
    tb$stimPosY == uxy$stimPosY[j] &
    !is.na(tb$headAng)
  
  thetaphi <- tb_v[ii, 1:2, drop=F] #base
  
  if (sum(is.na(thetaphi[,1])) < nrow(thetaphi)) {
    thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    pxyz <- xyz
    
    thetaphi <- tb_v[ii,3:4] #spk
    # thetaphi <- tb_v[ii,5:6] #sub
    
    thetaphi[, 1] <- 90 - thetaphi[, 1]
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    dxyz <- xyz
    
    df_tmp <- c(colMeans(cart2Mercator(pxyz)), 
                    colMeans(cart2Mercator(dxyz)),
                    as.integer(paste0(uxy$stimPosX[j],uxy$stimPosY[j]))
                    )
  }
  
  df <- rbind(df, df_tmp)
}

colnames(df) <- c('x','y','xend','yend','pos')
df5 <- df


# - new data
load('data/H2_tuning_2023_grating.rda' )

uxy <- unique(tb[, c('stimPosX', 'stimPosY')])

df <- matrix(ncol = 5, nrow = 0)
for (j in 1:nrow(uxy)) {
  ii <- tb$stimPosX == uxy$stimPosX[j] &
    tb$stimPosY == uxy$stimPosY[j]
  
  thetaphi <- tb_v[ii, 1:2, drop=F] #base
  
  if (sum(is.na(thetaphi[,1])) < nrow(thetaphi)) {
    thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    pxyz <- xyz
    
    thetaphi <- tb_v[ii,3:4] #spk
    # thetaphi <- tb_v[ii,5:6] #sub
    
    thetaphi[, 1] <- 90 - thetaphi[, 1]
    thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
    thetaphi <- thetaphi / 180*pi
    thetaphi <- cbind(1, thetaphi)
    xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
    xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
    
    dxyz <- xyz
    
    df_tmp <- c(colMeans(cart2Mercator(pxyz)),
                    colMeans(cart2Mercator(dxyz)),
                    as.integer(paste0(uxy$stimPosX[j],uxy$stimPosY[j]))
                    )
  }
  df <- rbind(df, df_tmp)
}
colnames(df) <- c('x','y','xend','yend','pos')


# PLOT
df <- rbind(df, df5)
# combine same position
df_arrow <- as.data.frame(df) %>%
  group_by(pos) %>%
  summarise(x=mean(x),
            y=mean(y),
            xend = mean(xend),
            yend = mean(yend),
            pos = mean(pos)) %>%
  as.data.frame()
df_arrow$pos <- as.factor(df_arrow$pos)

# mercator
plt <- plt_Mer +
  geom_segment(data=df_arrow, 
               aes(x = x,y = y, xend = xend,yend = yend,colour=pos),
               # colour='black', 
               size =1)+
  # scale_colour_manual(values = c( brewer.pal(7,"Set1")[c(1,2,5,6,7)], rep('black',5)) ) +
  scale_colour_manual(values = brewer.pal(10,"Paired"), guide="none") +
  scale_x_continuous(limits = c(-15,95)/180*pi,breaks = seq(-180,180,by=45)/180*pi, labels = paste0(seq(-180,180,by=45), "°"),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(-20,70)/180*pi,breaks = log(tan(pi/4 + seq(-75,75,by=15)/180*pi/2)), labels = paste0(seq(-75,75,by=15),"°"),  expand = c(0, 0)) +
  # annotate("text", x = df_arrow$x+0.02, y = df_arrow$y+0.02, label = 1:nrow(df_arrow), size = 3) + 
  labs(title = "grating - spk - all")
# windows(width = 9, height = 8)
windows(width = 6.5, height = 6.5)
plt
# ggsave(paste("H2_all_grating_avg.pdf",sep=''), width = 9, height = 8)

# average edge data --------------------------------------------------------

load('data/H2_tuning_2023.rda' )

# position uniquely determined by  c('stimPosX', 'stimPosY')
# combine arenaAng and headAnd,
uxy <- unique(tb[, c('stimPosX', 'stimPosY')])

df <- matrix(ncol = 6, nrow = 0)
for (j in 1:nrow(uxy)) {
  for (k in c(0,1)) {
    ii <- tb$stimPosX == uxy$stimPosX[j] &
      tb$stimPosY == uxy$stimPosY[j] &
      tb$edgeVal == k
    
    thetaphi <- tb_v[ii, 1:2] #base
    
    if (sum(is.na(thetaphi[,1])) < nrow(thetaphi)) {
      thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
      thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
      thetaphi <- thetaphi / 180*pi
      thetaphi <- cbind(1, thetaphi)
      xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
      xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
      
      pxyz <- xyz
      
      ### spk
      thetaphi <- tb_v[ii,3:4] #spk
      ### sub
      # thetaphi <- tb_v[ii,5:6] #sub
      
      thetaphi[, 1] <- 90 - thetaphi[, 1]
      thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
      thetaphi <- thetaphi / 180*pi
      thetaphi <- cbind(1, thetaphi)
      xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
      xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
      
      dxyz <- xyz
      
      # mean position
      pxyz_mean <- colMeans(pxyz)
      pxyz_mean <- pxyz_mean / sqrt(sum(pxyz_mean^2))
      pxyz_mean <- matrix(pxyz_mean, ncol = 3)
      # mean arrow head position
      dxyz_mean <- colMeans(dxyz - pxyz) + pxyz_mean
      dxyz_mean <- dxyz_mean / sqrt(sum(dxyz_mean^2))
      dxyz_mean <- matrix(dxyz_mean, ncol = 3)
      
      df_tmp <- c(unlist(cart2Mercator(pxyz_mean)),
                  unlist(cart2Mercator(dxyz_mean)),
                  as.integer(paste0(uxy$stimPosX[j],uxy$stimPosY[j])),
                  k)
      
      df <- rbind(df, df_tmp)
    }
  }
}
df <- data.frame(df)
colnames(df) <- c('x','y','xend','yend','pos', 'edge')

### ### combine bright/dark
df_arrow <- df %>%
  group_by(pos) %>%
  summarise(x=first(x),y=first(y),xend=mean(xend),yend=mean(yend),edge=max(edge)) %>%
  ungroup() %>%
  as.data.frame()

### ### not combine
# df_arrow <- as.data.frame(df)

df_arrow$edge <- as.factor(df_arrow$edge)

# # longer arrows
# df_arrow$xend <- 1 * (df_arrow$xend - df_arrow$x) + df_arrow$x
# df_arrow$yend <- 1 * (df_arrow$yend - df_arrow$y) + df_arrow$y

# mercator
plt <- plt_Mer +
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=edge),size =1)+
  # geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=cell),size =1, arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
  # scale_colour_manual(values = brewer.pal(7,"Set1"), breaks = c("1", "2", "3","4","5","6","7")) +
  # scale_colour_manual(values = c('black', 'red'),
  #                     breaks = c("0", "1"),
  #                     labels = c('dark', 'bright')) +
  scale_colour_discrete(guide="none") +
  scale_x_continuous(limits = c(-15,95)/180*pi,breaks = seq(-180,180,by=45)/180*pi, labels = paste0(seq(-180,180,by=45), "°"),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(-20,70)/180*pi,breaks = log(tan(pi/4 + seq(-75,75,by=15)/180*pi/2)), labels = paste0(seq(-75,75,by=15),"°"),  expand = c(0, 0)) +
  annotate("text", x = df_arrow$x+0.02, y = df_arrow$y+0.02, label = df_arrow$pos, size = 3) +
  labs(title = "edge - spk ")
# windows(width = 9, height = 8)
windows(width = 6.5, height = 6.5)
plt
# ggsave(paste("H2_new_edge_spk.pdf",sep=''), width = 9, height = 8)
# ggsave(paste("H2_new_edge_sub.pdf",sep=''), width = 9, height = 8)



# ED Fig.1E, dark vs. bright vs. grating at each location ---------------------

load('data/H2_tuning_2023.rda' )

# position uniquely determined by  c('stimPosX', 'stimPosY')
# make a factor for plotting
uxy <- unique(tb[, c('edgeVal', 'stimPosX', 'stimPosY')]) %>%
  # mutate(stimPosX = as.character(stimPosX),
  #        stimPosY = as.character(stimPosY) ) %>%
  mutate(stimType = paste0(edgeVal,stimPosX, stimPosY),
         stimTypePos = paste0(stimPosX, stimPosY))
# uxy <- lapply(uxy, factor)
# recode(uxy, "c('A', 'B')='A+B';c('D', 'E') = 'D+E'")


#base
thetaphi <- tb_v[,1:2] 
thetaphi[, 1] <- 90 - thetaphi[, 1]
thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
thetaphi <- thetaphi / 180*pi
thetaphi <- cbind(1, thetaphi)
xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
# xyz[!is.na(thetaphi[,3]), ] <- sph2cartZ(thetaphi[!is.na(thetaphi[,3]), ])
xyz <- sph2cartZ(thetaphi) %>% cart2Mercator()
dxyz <- xyz

#spk
thetaphi <- tb_v[,3:4] 
thetaphi[, 1] <- 90 - thetaphi[, 1]
thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
thetaphi <- thetaphi / 180*pi
thetaphi <- cbind(1, thetaphi)
xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
xyz <- sph2cartZ(thetaphi) %>% cart2Mercator()
pxyz <- xyz

#sub
# thetaphi <- tb_v[ii,5:6] 

# angle with x-axis
df_ang <- cbind(dxyz, pxyz)
colnames(df_ang) <- c('x0','y0','xsp','ysp')
df_ang <- as.data.frame(df_ang) %>%
  transmute(ang = atan2(ysp-y0, -(xsp-x0)) /pi*180)

# # at each location
# df <- cbind(df_ang, tb) %>%
#   as_tibble() %>%
#   group_by(stimPosX,stimPosY,edgeVal) %>%
#   summarise(
#     ang_mean = mean(ang),
#     ang_sd = sd(ang)
#   ) %>%
#   as.data.frame()

df <- cbind(df_ang, tb[,1:5])
df <- merge(df, uxy, by=c('edgeVal','stimPosX', 'stimPosY'))
# df$stimType <- factor(
#   df$stimType,
#   levels=c('02913', '12913', '05828', '15828', '08736', '18736',
#            '08721', '18721', '08729', '18729', '08713', '18713')
#   )


# - add gratings, 2023
load('data/H2_tuning_2023_grating.rda' )

uxy <- unique(tb[, c('stimPosX','stimPosY')]) %>%
  mutate(stimType = paste0(2,stimPosX, stimPosY),
         stimTypePos = paste0(stimPosX, stimPosY))

#base
thetaphi <- tb_v[,1:2] 
thetaphi[, 1] <- 90 - thetaphi[, 1]
thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
thetaphi <- thetaphi / 180*pi
thetaphi <- cbind(1, thetaphi)
xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
xyz <- sph2cartZ(thetaphi) %>% cart2Mercator()
dxyz <- xyz

#spk
thetaphi <- tb_v[,3:4] 
thetaphi[, 1] <- 90 - thetaphi[, 1]
thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
thetaphi <- thetaphi / 180*pi
thetaphi <- cbind(1, thetaphi)
xyz <- matrix(ncol = 3, nrow = nrow(thetaphi))
xyz <- sph2cartZ(thetaphi) %>% cart2Mercator()
pxyz <- xyz

#sub
# thetaphi <- tb_v[ii,5:6] 

# angle with x-axis
df_ang <- cbind(dxyz, pxyz)
colnames(df_ang) <- c('x0','y0','xsp','ysp')
df_ang <- as.data.frame(df_ang) %>%
  transmute(ang = atan2(ysp-y0, -(xsp-x0)) /pi*180)

df2 <- cbind(df_ang, tb[,1:4])
df2 <- merge(df2, uxy, by=c('stimPosX', 'stimPosY'))
# df2$stimType <- factor(
#   df2$stimType,
#   levels=c('22913', '25828', '28736','28721', '28729', '28713')
# )
df2 <- cbind(2, df2)
colnames(df2)[1] <- "edgeVal"

# - combine
df3 <- rbind(df, df2)
df3$stimType <- factor(
  df3$stimType,
  levels=c('02913', '12913', '22913', "s1", '05828', '15828', '25828', "s2",
           '08736', '18736', '28736', "s3", '08721', '18721', '28721', "s4",
           '08729', '18729', '28729', "s5", '08713', '18713', '28713')
)

# for MR, edgeVal dark=0, bright=1, gratings=2
# write.csv(df3, file = "Eyal_exp2.csv", row.names = F)

# PLOT
windows(width = 10, height = 5)
# pdf(paste("H2_edge_DvsBvsG", '.pdf', sep = ''), width = 10, height = 5)
ggplot(df3, aes(stimType, ang, colour= factor(edgeVal)) ) + 
  geom_jitter(size = 1,height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               aes(colour = factor(edgeVal)),
               size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  aes(colour = factor(edgeVal)), size = 1, position = position_nudge(x = 0.2, y = 0)) +
  scale_discrete_manual(values = c('#762A83','#1B7837','#000000'),aesthetics = "colour",guide= guide_legend(title="")) +
  coord_cartesian(ylim = c(-75, 50)) +
  scale_y_continuous(breaks= seq(-75,50,by=25), labels = paste0(seq(-75, 50,by=25), "°")) +
  scale_x_discrete(drop = F,
                   labels=c('02913'='1D', '12913'='1B', '22913'='1G', 's1'='',
                            '05828'='2D', '15828'='2B', '25828'='2G', 's2'='',
                            '08736'='3D', '18736'='3B', '28736'='3G', 's3'='',
                            '08721'='4D', '18721'='4B', '28721'='4G', 's4'='',
                            '08729'='5D', '18729'='5B', '28729'='5G', 's5'='',
                            '08713'='6D', '18713'='6B', '28713'='6G')) +
  theme_minimal() +
  labs(title= paste("H2 tuning", sep = ''),
       x='position', y='Angle wrt +h [deg]') 
# dev.off()


# -- 2-way anova, no interaction
df4 <- df3
df4$stimTypePos <- factor(df4$stimTypePos)
df4$edgeVal <- factor(df4$edgeVal)
resp <- aov(ang ~ stimTypePos + edgeVal, data = df4)
summary(resp)

# dev.new()
# par(mfrow=c(2,2))
# plot(resp)
# 
# tukey_resp <- TukeyHSD(resp)
# tukey_resp$stimType[order(tukey_resp$stimType[,4]), ]
# 
# dev.new()
# plot(tukey_resp, las = 1)

# -- wilcox test, avg bright/dark bars
df_test <- df3[df3$edgeVal < 2, ] 
test_wilc <- matrix(nrow = 0, ncol = 2)
for (j in unique(df_test$stimTypePos)) {
  data <- df_test[df_test$stimTypePos == j,]
  test_wilc <- rbind(
    test_wilc,
    c(j,
      wilcox.test(data[,'ang'])$p.value
      )
  )
}


# ED Fig.1F, deformation, mod from "new data 2023 edge" ------------------------

# make template, center plus 8 directions
# 22.5/1.25=18
NN <- 8.5 # num of pixel as radius

NXY <- data.frame(dx = c(0,NN,0,-NN,0,NN/sqrt(2),-NN/sqrt(2),-NN/sqrt(2),NN/sqrt(2)),
                  dy = c(0,0,NN,0,-NN,NN/sqrt(2),NN/sqrt(2),-NN/sqrt(2),-NN/sqrt(2)))

# stim center
tb <- read_excel("data_Eyal/H2_2023/movEdgeTabwPY_ThetaCorr_fullResp.xlsx", sheet = 1, col_names = T ) %>% data.frame()
names(tb)[names(tb) == "headAngPitch"] <- "headAng"

stimpos <- tb %>% 
  mutate(pos = as.integer(paste0(stimPosX,stimPosY))) %>%
  group_by(pos) %>%
  summarise(stimPosX = first(stimPosX),
            stimPosY = first(stimPosY),
            arenaAng = mean(arenaAng),
            headAng = mean(headAng)) %>%
  as.data.frame()
# stimpos <- data.frame(stimPosX = c(29,87,87, 58,87,87),
#                  stimPosY = c(13,13,29, 28,21,36),
#                  arenaAng = c(49,49,49, 19,19,19),
#                  headAng = rep(64,6))

pos <- merge(NXY, stimpos, all=T) %>%
  transmute(arenaAng = arenaAng, headAng = headAng,
            px = stimPosX + dx, py = stimPosY + dy)

# area params
A <- (4+10/32 + 4+20/32)/2  #radius
Dx <- (5+6/32) - A  # area center wrt to fly
Dy <- (4+20/32) - A
Dz <- - (1+ 13/32)

# LED size, assume 1.25 at equator
DL <- A * 1.25/180*pi

# response amp factor
resfac <- 1

# - loop
stims <- matrix(ncol = 2, nrow = nrow(pos)) # [elev azim]

for (j in 1:nrow(stims)) {
  arot <- pos$arenaAng[j] /180*pi #arena rotation ang
  afly <- pos$headAng[j] /180*pi #eqator below holder
  
  # cylindrical arena, 2d
  #angle from sagittal plane
  ax <- (120 - pos$px[j]*1.25) /180*pi
  #dist from bottom
  ay <- DL * pos$py[j]
  
  
  # go to lab coord via a translation
  # origin at arena center, xy-plane is level
  x <- A * cos(ax)
  y <- A * sin(ax)
  z <- ay - (4 + 24/32)/2
  
  
  xyz <- matrix(c(x,y,z), ncol=3)
  
  # arena rotation, 
  # arot is angle below ground, use -arot
  xyz <- matrix(c(cos(-arot), 0, -sin(-arot),
                  0, 1, 0,
                  sin(-arot), 0, cos(-arot)), ncol=3, byrow = T) %*% t(xyz)
  xyz <- t(xyz)
  
  # go to fly's coord
  xyz <- xyz + c(Dx, Dy, Dz)
  
  # fly rotation 
  # afly is angle below ground, use afly since it's ref that's rotating
  xyz <- matrix(c(cos(afly), 0, -sin(afly),
                  0, 1, 0, 
                  sin(afly), 0, cos(afly)), ncol = 3, byrow = T) %*% t(xyz)
  xyz <- t(xyz)
  
  # to spherical coord
  thetaphi <- cart2sphZ(xyz) /pi*180
  
  stims[j, ] <- round(c(90 - thetaphi[2], thetaphi[3]), 2)
}


# - Mercator projection
thetaphi <- stims

thetaphi[, 1] <- 90 - thetaphi[, 1] # elev to theta
thetaphi[, 2] <- 360 - thetaphi[, 2] # left to right, need it cos cart2Merc invert y
thetaphi <- thetaphi / 180*pi
thetaphi <- cbind(1, thetaphi)

xyz <- sph2cartZ(thetaphi)
xyz_mer <- cart2Mercator(xyz)


# PLOT
df <- as.data.frame(cbind(xyz_mer,
                          factor(rep(seq(1,6), each=9)),
                          rep(c(1,2,rep(1,7)), 6)
                          ))
colnames(df) <- c('x','y','cell','size')

df_arrow <- matrix(ncol = 4, nrow = 0)
for (j in 0:5) {
  df_arrow <- rbind(df_arrow,
              cbind(
                matrix(xyz_mer[9*j+1,], nrow=8, ncol = 2, byrow = T),
                xyz_mer[9*j+2:9,]
              ))
}
df_arrow <- cbind(df_arrow, rep(seq(1,6), each=8))
df_arrow <- matrix(unlist(df_arrow), ncol = 5) 
df_arrow <- data.frame(df_arrow)
colnames(df_arrow) <- c('x','y','xend','yend','cell')

# mercator
plt <- plt_Mer +
  geom_point(data=df, aes(x = x,y = y,colour=cell,size =size))+
  geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour= as.factor(cell)))+
  # geom_segment(data=df_arrow, aes(x = x,y = y, xend = xend,yend = yend, colour=cell),size =1, arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
  # scale_colour_manual(values = brewer.pal(7,"Set1"), breaks = c("1", "2", "3","4","5","6","7")) +
  scale_colour_manual(values = pal_9[c(5,4,3,6,2,1)],
                      breaks = as.character(seq(1,6))) +
  # scale_x_continuous(limits = c(0,90)/180*pi,breaks = seq(-45,90,by=45)/180*pi, labels = paste0(seq(-45,90,by=45),"°"),  expand = c(0, 0)) +
  # scale_y_continuous(limits = c(-20,75)/180*pi,breaks = log(tan(pi/4 + seq(-15,75,by=15)/180*pi/2)), labels = paste0(seq(-15,75,by=15),"°"),  expand = c(0, 0)) +
  scale_x_continuous(limits = c(-15,95)/180*pi,breaks = seq(-180,180,by=45)/180*pi, labels = paste0(seq(-180,180,by=45), "°"),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(-20,70)/180*pi,breaks = log(tan(pi/4 + seq(-75,75,by=15)/180*pi/2)), labels = paste0(seq(-75,75,by=15),"°"),  expand = c(0, 0)) +
  # annotate("text", x = df_arrow$x+0.02, y = df_arrow$y+0.02, label = df$pos, size = 3) +
  theme(legend.position="none")+
  labs(title = "templates ")
# windows(width = 8, height = 8)
windows(width = 6.5, height = 6.5)
plt
# ggsave("arena_deform.pdf", width = 8, height = 8)


# - avg amplitudes 1 vs 4
p1 <- sweep(xyz[2:9,], 2, xyz[1,,drop=F], '-')^2 %>% rowSums() %>% mean() %>% sqrt()
p4 <- sweep(xyz[3*9+2:9,], 2, xyz[3*9+1,,drop=F], '-')^2 %>% rowSums() %>% mean() %>% sqrt()
abs(p1-p4)/p1


