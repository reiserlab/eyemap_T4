

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

df <- matrix(ncol = 3, nrow = 0)
for (ii_h in 1:4) {
  if (nrow(ssn[ssn$so == ii_h,]) > 0 ) {
    hh <- hist(ssn[ssn$so == ii_h,]$angSO4, breaks = seq(0,180,dtheta), plot = F)
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
  facet_grid(rows = vars(so) ) +
  theme_minimal() +
  theme( axis.title = element_blank(),
         axis.text.y = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.y = element_blank()) +
  labs(title = paste("T4", letters[LL], "_SO_histo_areaNorm_haxis", sep = '') )
# dev.off()

# ED Fig.3B, 3x3 compartments --------------------------------------------------------

com_xyz_eye_b <- lens_type[[2]][, c('comx','comy','comz')] %>% as.matrix()
com_xyz_eye_d <- lens_type[[4]][, c('comx','comy','comz')] %>% as.matrix()


fnnb <- 19
# - T4b
xyz <- com_xyz_eye_b
xyz[,2] <- -xyz[,2]
tp <- cart2sphZ(xyz)[,2:3] /pi*180
tp[,2] <- if_else(tp[,2] > 180, tp[,2] -360, tp[,2])

angdiv <- list()
inddiv <- list()
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

# -- 3x3 plot
df0 <- as.data.frame(refhex_2)
colnames(df0) <- c('x','y')

windows(width = 9, height = 9)
# pdf(paste("hex_ref_b", '_3x3.pdf', sep = ""), width = 9, height = 9)
layout(matrix(seq(1,9),nrow = 3, ncol = 3, byrow = TRUE))
for (k in 1:length(inddiv)) {
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

