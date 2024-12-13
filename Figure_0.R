# Overview ---------------------------------------------------------------

# This script sets the stage (library/func/data/color+plot conventions/etc) while other "Figure_*.R" scripts plot figures.
# It loads processed data. To re-process, see:
# eyemap_uCT
# eyemap_lens_med_lop_v3.R
# eyemap_T4_v2
# eyemap_T4_RF

# Notes on some variables
# eyemap=[rownames of med_xyz, rownames of ucl_rot_sm]. Note that all these 3 variables have the same row order/index.
# Mi1_ixy=[rownames of ucl_rot_sm, x, y]

# Notation and convention
# "## ##" = need to specify a choice, say, between T4b and T4d
# "Mollweide" and "Mercator" = different projections of the same plot

# Plotting:
# By default, the plot is shown in a plotting device, usually with "nopen()" or "windows()".
# To save to disk, look for commented commands:
# "pdf()" paired with "dev.off()", "ggsave()", or "rgl.snapshot()".


# load library ------------------------------------------------------------

library(natverse)
library(tidyverse)
library(RColorBrewer) #palette
library(readxl)
library(alphashape3d) # ashape3d
library(np) # kernel regression to find the base pt on the eye
library(ggExtra) #ggMarginal
library(alphahull) # ashape
library(reshape2) #melt
library(cowplot)#ggdraw

# library("plotrix")

# # library(deldir)# for smoothing uCT data
# # library(sp) #in.poly
# # library(plotly) # polar plot
# # library(orca)
# # library(processx)
# # library(clue) # Hungarian matching
# # library(colorspace)

# # clean everythign up.
# rm(list=ls())
# #close any open rgl windows
# while (rgl.cur() > 0) { close3d() }
# # close opened dev
# while (dev.cur() > 1) { dev.off() }

source("eyemap_func.R")

# server ------------------------------------------------------------------

# catmaid_login(server= catmaid.server,
#               authname= catmaid.authname,
#               authpassword= catmaid.authpassword,
#               token= catmaid.token)

vfbcatmaid("fafb") # https://catmaid.virtualflybrain.org/

print(catmaid_login())
# expected output: 
# Connection to catmaid server:
# https://fafb.catmaid.virtualflybrain.org
# Login active since: Mon, 25 Nov 2024 21:08:52 GMT

# for MD only ----------------------------------------------------------

# PR, skid=17640612
# L1, skid=17572546 
# central vertical line of Mi1
# 14858811 14926430 13780764 13774456 11250042 14908550 13870804
# T4, same as T4_ind_eg
# 9880570, 11364711, 11301434, 11347884


# load uCT data ---------------------------------------------------------

# cf. eyemap_uCT.R

# old
# load(paste("../microCT/data/eq_corr/", '12102019_female', ".RData", sep=''))

# 2023
# fn <- "2023_eyemap/20230926"
# fn <- "2023_eyemap/20231107" #gg
# fn <- "2023_eyemap/20240206"
# fn <- "2023_eyemap/20240510"
# fn <- "2023_eyemap/20240513"
# fn <- "2023_eyemap/20240520"
# fn <- "2023_eyemap/20240522"
# fn <- "2023_eyemap/20240524"
# fn <- "2023_eyemap/20240530" #gg
# fn <- "2023_eyemap/20240530_dried"
# fn <- "2023_eyemap/20240612"

# fn <- "2023_eyemap/20240701" #gg
# load(paste0("../microCT/", fn, '.RData'))
# load(paste0("../microCT/", fn, '_nb.RData'))
# load(paste0("../microCT/", fn, "_normals.RData"))

load(paste0("data/microCT/20240701", ".RData"))
load(paste0("data/microCT/20240701", "_nb.RData"))
load(paste0("data/microCT/20240701", "_normals.RData"))
load(paste0("data/microCT/20240701", "_roc.RData"))
load(paste0("data/microCT/20240701", "_dia.RData"))

# separate left and right, re-indexing
lens_left <- lens[ind_left_lens,]
lens_mir <- lens_left
lens_mir[,2] <- -lens_mir[,2]
lens_2eye <- lens
lens <- lens[!ind_left_lens,]
ind_Up_lens <- na.omit(match(i_match[ind_Up], rownames(lens)))
ind_Up_lens_left <- na.omit(match(i_match[ind_Up], rownames(lens_left)))
ind_Down_lens <- na.omit(match(i_match[ind_Down], rownames(lens)))
ind_Down_lens_left <- na.omit(match(i_match[ind_Down], rownames(lens_left)))
cone_lens <- i_match[!ind_left_cone]
# re-order ucl_rot rows and re-define rownames cos mapping is betw Mi1 and lens
# ucl_rot_2eye <- ucl_rot #2021
ucl_rot_2eye <- ucl_rot_sm #2023

ucl_rot_left <- ucl_rot_2eye[order(i_match),][ind_left_lens,]
colnames(ucl_rot_left) <- c('x','y','z')
ucl_rot_right <- ucl_rot_2eye[order(i_match),][!ind_left_lens,]
colnames(ucl_rot_right) <- c('x','y','z')
rm(ucl)
i_match_uct <- i_match
rm(i_match)

# for Emil
# save(ucl_rot_left, lens_left, file = "data/ucl_rot_left.rda") 

# for Brad
# write.csv(ucl_rot_right, file = "for_Brad/20240701_ucl_right.csv", row.names = F)
# write.csv(ucl_rot_left, file = "for_Brad/20240701_ucl_left.csv", row.names = F)

# for male-cns, see for_malecns/
# write.csv(ucl_rot_sm, file = "for_malecns/20230926/ucl_rot_sm.csv", row.names = F)
# write.csv(nb_ind, file = "for_malecns/20230926/nb_ind.csv", row.names = F)
# write.csv(ind_xy, file = "for_malecns/20230926/ind_xy.csv", row.names = F)

# for Max Josch
# tb <- cbind(ucl_rot_2eye, as.integer(ind_left_cone))
# colnames(tb) <- c('x','y','z','left')
# write.csv(tb, file = "ommatidia_dir_202309.csv", row.names = F)
# write.csv(tb, file = "ommatidia_dir_2019.csv", row.names = F)

# # DEBUG
# nopen3d()
# points3d(lens_2eye, size = 7, color = 'pink', alpha=1)
# points3d(cone, size = 5, color = 'blue')
# linemat <- matrix(t(cbind(lens_2eye[i_match_uct,], cone)), ncol = 3, byrow = T)
# points3d(lens, size = 7, color = 'pink', alpha=1)
# points3d(cone[!ind_left_cone,], size = 5, color = 'blue')
# linemat <- matrix(t(cbind(lens[match(i_match_uct[!ind_left_cone],rownames(lens)),], cone[!ind_left_cone,])), ncol = 3, byrow = T)
# segments3d(linemat, color = "grey")

# H2 data -----------------------------------------------------------------

# load('data/H2_tuning.rda' )

# tb_v [theta(base), phi(base), spk, spk, sub, sub, dark=0/bright=1]
load('data/H2_tuning_2023.rda' )

# find unit hex eg.------------------------------------------------------------

# dev.new()
# plot(ucl_rot_Mo)
# points(ucl_rot_Mo[c(129, 464, 743),], cex=1.3, pch=2,col='cyan')
# points(cmer, col='orange', pch=16)
# points(bkgd_eq, type='l')
# points(bkgd_mer_e, type='l')
# points(bkgd_mer, type='l')
# # points(ucl_rot_Mo[match(haxis_gen(0), eyemap[,2]),], col='black', pch=16)
# points(ucl_rot_Mo[nb_ind[match(129, nb_ind[,1]),],], col='magenta', pch=16)

# ii <- identify(ucl_rot_Mo)


# asp ratio, [DD D M]
# ind_eghex <-  c(185, 72, 424) # 2024 0206, 
ind_eghex <-  c(409, 223, 107) # 2024 0701

# skewness, [D M V]
# ind_eghex_2 <- c(185, 424, 745) # 2024 0206
ind_eghex_2 <- c(210, 348,  559) # 2024 0701
# ind_eghex_2 <- c(211, 347, 592) 
# ind_eghex_2 <- c(177, 305,  660)  


# load EM data ---------------------------------------------------------------

# [index, p, q] of lens/medulla points
load('data/lens_ixy.RData')
load('data/med_ixy.RData')

# index of special Mi1: _DRA, _PR, _noPR
load('data/Mi1_ind.RData') 

# neighboring points index and distance, from 'eyemap_lens_med_lop.R'
load('data/hexnb_ind_dist.RData') 

# load EM neuron H2
load('data/neu_H2.RData')
# Alt. load from database
# anno_H2 <- catmaid_query_by_name("^H2 neuron \\(right", type = 'neuron')
# H2 <- read.neuron.catmaid(anno_H2$skid, .progress = 'text')
# save(anno_H2, H2, file = 'data/neu_H2.RData')

# PR
load('data/neu_PR.RData')

# L1
load('data/neu_L1.RData')

# R7
load('data/neu_R7.RData')

# Tm5a and lo xyz 
load('data/neu_TmY5a.RData')  

# Mi1
load("data/neu_Mi1.RData")

# T4
load('data/neu_T4b.RData')
load("data/neu_T4_dend.RData")
load('data/T4_hs.RData')

# T4 RF
load('data/T4_RF_pred.RData')

# T4 branching patten
load('data/T4_gallery.RData')

# eyemap, NB, load this after loading uCT data
load("data/eyemap.RData")

# chiasm
load('data/chiasm.RData')

# LOP xyz
load('data/L2_med.RData') 

# LOP layers, xyz_layer_T4, from p_layerLOP
load("data/LayerData_T4_20210729_7e7_1um.RData")

# # - JFRC2010 mesh, from Greg Jefferis, cf.p_R7R8
load("data/JFRC2NP.surf.fafb.rda")
ME_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="ME_R")
# LOP_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LOP_R")
LO_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LO_R")
# # fix LOP_msh
# xyz <- t(LOP_msh$vb)[-ii,1:3]
# msh.a <- ashape3d(xyz, alpha = 20000) # 20000 look ok
# LOP_msh_mod <- as.mesh3d(msh.a)
# save(LOP_msh_mod, file = "data/LOP_msh_mod.RData")

load('data/LOP_msh_mod.RData')

# # load a single LOP mesh, missing a corner
# load("../../p_layerLOP/LP_mesh.RData") # LP


# special pts in medulla ------------------------
# cf "eyemap_lens_med_lop"

# pt_100_med <- 10 # could be 12 as well, 2021
# pt_100_med <- 8 # 2023

# 2024
pt_100 <- 10 
pt_100_Mi1 <- pt_100
pt_up <- 315
pt_up_Mi1 <- pt_up

# color palette --------------------------------------------------------------

pal_9 <- brewer.pal(9,"Paired")
pal_lr <- c("light blue", "gray50") #uCT left and right eye
pal_axes <- c('#227722', '#89DCB3','black', 'orange') # +p/q/v/h axes
pal_heat1 <- c("blue","#e7d4e8","red")
pal_heat2 <- rev(brewer.pal(8,"PuBu")[c(2,4,6,8)])
pal_so <- c('gray70', '#fdbb84', '#d7301f', 'gray35', 'gray35')
pal_T4 <- c("turquoise", "#A6611A",  "royalblue", "plum")
# -tx -ty tz -rx ry -rz
pal_TR <- c("#B48ADC", "#FFC125","#278384","#3D29A3","#A37B29","#66CD00")

# exemple T4 ------------------------------------------------------------

T4_ind_eg <- c(9880570, 11364711, 11301434, 11347884) # match(T4_ind_eg, anno_overlap$skid) = c(7, 10, 8 ,9)

T4_eg <- neuronlist()
for (LL in 1:4) {
  skid <- sapply(T4_dend[[LL]], function(x) x$skid)
  T4_eg <- c(T4_eg, T4_dend[[LL]][na.omit(match(T4_ind_eg, skid))])
}


# example Mi1 -----------------------------------------------------------

Mi1_eg_skid <- c(14858811, 10888497, 10850334, 13726215, 15137473, 14567424, 13975332)
Mi1_eg_ind <- match(Mi1_eg_skid, anno_Mi1$skid)

# catmaid_query_connected(Mi1_eg_skid[1])$outgoing %>%
#   mutate(name = catmaid_get_neuronnames(partner))

# nopen3d()
# points3d(med_xyz)
# points3d(med_xyz[match(Mi1_neu_ind[c(334,583,10)], rownames(med_xyz)),],size=20,col='black')

# # load neurons ------------------------------------------------------------
# 
# # -- T5, few annotated
# T5 <- list()
# for (LL in 1:4) {
#   anno_str_axon <- paste("T5",letters[LL], " - axon", sep = "")
#   anno_str_dend <- paste("T5",letters[LL], " - dendrites", sep = "")
#   # anno_str_comp <- paste("T5",letters[LL], " - complete", sep = "")
#   anno_axon <- catmaid_query_by_annotation(anno_str_axon)
#   anno_dend <- catmaid_query_by_annotation(anno_str_dend)
#   # anno_comp <- catmaid_query_by_annotation(anno_str_comp)
#   neu_skid_axon <- anno_axon[,"skid"]
#   neu_skid_dend <- anno_dend[,"skid"]
#   # neu_skid_comp <- anno_comp[,"skid"]
#   T5_LL <- read.neurons.catmaid(neu_skid_dend, .progress = 'text')
#   T5[[LL]] <- kinkCT1(T5_LL)
# }
# 
# # anno_T5 <- catmaid_query_by_name("^Putative T5 ", type = 'neuron')
# # T5 <- read.neurons.catmaid(anno_T5$skid, .progress = 'text')


# Mollweide guidelines and ucl_rot_sm -----------------------------------------

# - ommatidia directions
ucl_rot_Mo <- ucl_rot_sm
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
rownames(ucl_rot_Mo) <- rownames(ucl_rot_sm)

# full left
ucl_rot_Mo_left <- ucl_rot_left
colnames(ucl_rot_Mo_left) <- c('x','y','z')
ucl_rot_Mo_left %<>% as_tibble() %>%  
  mutate(y = -y) %>%
  mutate(theta = acos(z)) %>%
  mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
  mutate(t = theta / pi * 180, p = phi/pi*180) %>%  
  mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
  as.data.frame()
ucl_rot_Mo_left <- Mollweide(ucl_rot_Mo_left[,c('t', 'p')])
colnames(ucl_rot_Mo_left) <- c('xM','yM')
rownames(ucl_rot_Mo_left) <- rownames(ucl_rot_left)

# full right
ucl_rot_Mo_right <- ucl_rot_right
colnames(ucl_rot_Mo_right) <- c('x','y','z')
ucl_rot_Mo_right %<>% as_tibble() %>%  
  mutate(y = -y) %>%
  mutate(theta = acos(z)) %>%
  mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
  mutate(t = theta / pi * 180, p = phi/pi*180) %>%  
  mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
  as.data.frame()
ucl_rot_Mo_right <- Mollweide(ucl_rot_Mo_right[,c('t', 'p')])
colnames(ucl_rot_Mo_right) <- c('xM','yM')
rownames(ucl_rot_Mo_right) <- rownames(ucl_rot_right)

# - Mollweide guidelines
Mollweide_ori <- c(0,0)
Mollweide_mul <- 1
bydeg <- 5

bkgd_eq <- Mollweide(cbind(90, seq(-180, 180, by = bydeg)))
bkgd_eq <- sweep(bkgd_eq*Mollweide_mul, 2, Mollweide_ori, '+')
colnames(bkgd_eq) <- c('xM','yM')
bkgd_eq_p45 <- Mollweide(cbind(45, seq(-180, 180, by = bydeg)))
bkgd_eq_p45 <- sweep(bkgd_eq_p45*Mollweide_mul, 2, Mollweide_ori, '+')
colnames(bkgd_eq_p45) <- c('xM','yM')
bkgd_eq_m45 <- Mollweide(cbind(135, seq(-180, 180, by = bydeg)))
bkgd_eq_m45 <- sweep(bkgd_eq_m45*Mollweide_mul, 2, Mollweide_ori, '+')
colnames(bkgd_eq_m45) <- c('xM','yM')
bkgd_mer_ww <- Mollweide(cbind(seq(0, 180, by = bydeg), rep(-180, 180/bydeg+1)))
bkgd_mer_w <- Mollweide(cbind(seq(180, 0, by = -bydeg), rep(-90, 180/bydeg+1)))
bkgd_mer_c <- Mollweide(cbind(seq(0, 180, by = bydeg), rep(0,180/bydeg+1)))
bkgd_mer_e <- Mollweide(cbind(seq(180, 0, by = -bydeg), rep(90,180/bydeg+1)))
bkgd_mer_ee <- Mollweide(cbind(seq(0, 180, by = bydeg), rep(180,180/bydeg+1)))
bkgd_mer <- rbind(bkgd_mer_ww, bkgd_mer_w, bkgd_mer_c,bkgd_mer_e,bkgd_mer_ee)
bkgd_mer <- sweep(bkgd_mer*Mollweide_mul, 2, Mollweide_ori, '+')
colnames(bkgd_mer) <- c('xM','yM')

# -- whole
plt_Mo <- ggplot() +
  geom_path(data= as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour= 'grey50') +
  geom_path(data= as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour= 'grey50') +
  geom_path(data= as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour= 'grey50') +
  geom_path(data= as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
  scale_x_continuous(limits = c(-sqrt(8), sqrt(8)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-sqrt(2), sqrt(2)), breaks = c(-1.5,0,1.5), labels = c(-1.5,0,1.5), expand = c(0, 0)) + # set +y as above eq
  theme_void() +
  theme(legend.position="none", panel.background = element_blank()) +
  coord_fixed(ratio = 1)

# -- 3/4
bkgd_eq <- Mollweide(cbind(90, seq(-90, 180, by = bydeg)))
colnames(bkgd_eq) <- c('xM','yM')
bkgd_eq_p45 <- Mollweide(cbind(45, seq(-90, 180, by = bydeg)))
colnames(bkgd_eq_p45) <- c('xM','yM')
bkgd_eq_m45 <- Mollweide(cbind(135, seq(-90, 180, by = bydeg)))
colnames(bkgd_eq_m45) <- c('xM','yM')
bkgd_mer <- rbind(bkgd_mer_w, bkgd_mer_c,bkgd_mer_e,bkgd_mer_ee)
colnames(bkgd_mer) <- c('xM','yM')

plt_Mo34 <- ggplot() +
  geom_path(data= as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data= as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data= as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data= as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
  theme_void() +
  theme(legend.position="none", panel.background = element_blank()) +
  scale_x_continuous(limits = c(- 2*sqrt(2)/2, 2*sqrt(2)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-sqrt(2), sqrt(2)), expand = c(0, 0)) +
  coord_fixed(ratio = 1)

# -- min
bkgd_eq <- Mollweide(cbind(90, seq(-30, 150, by = bydeg)))
colnames(bkgd_eq) <- c('xM','yM')
bkgd_mer <- rbind(bkgd_mer_c,bkgd_mer_e)
colnames(bkgd_mer) <- c('xM','yM')

plt_Momin <- ggplot() +
  geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
  geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
  theme_void() +
  theme(legend.position="none", panel.background = element_blank()) +
  scale_x_continuous(limits = c(- 2*sqrt(2)/2, 2*sqrt(2)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-sqrt(2), sqrt(2)), expand = c(0, 0)) +
  coord_fixed(ratio = 1)
  
# theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
#       axis.title.x = element_text(size = 14),
#       axis.title.y = element_text(size = 14, angle = 90),
#       panel.grid = element_blank(), legend.position = c(.9, .9) )


# equator and meridian strips -----------------------------------------------
str_hw <- 15
bkgd_str_equa <- Mollweide(
  rbind(
    cbind(rep(90- str_hw, 10), seq(-45, 160, length.out=10)),
    cbind(seq(90- str_hw, 90+ str_hw, length.out=5), rep(160, 5)),
    cbind(rep(90+ str_hw, 10), seq(160, -45, length.out=10)),
    cbind(seq(90+ str_hw, 90- str_hw, length.out=5), rep(-45, 5)) ) )
colnames(bkgd_str_equa) <- c('xM','yM')
bkgd_str_meri <- Mollweide(
  rbind(
    cbind(seq(0, 180, length.out=20), rep(45- str_hw,20)),
    cbind(seq(180, 0, length.out=20), rep(45+ str_hw,20)) ) )
colnames(bkgd_str_meri) <- c('xM','yM')

# Mercator ----------------------------------------------------------------

ucl_rot_Merc_right <- cart2Mercator(ucl_rot_right) %>% as.data.frame()
rownames(ucl_rot_Merc_right) <- rownames(ucl_rot_right)

ucl_rot_Merc <- cart2Mercator(ucl_rot_sm) %>% as.data.frame()
rownames(ucl_rot_Merc) <- rownames(ucl_rot_sm)

xlat <- seq(-45,180,by=45)
xtick <- xlat/180*pi
ylat <- c(seq(-75,75,by=15), 85)
ytick <- log(tan(pi/4 + ylat/180*pi/2))

plt_Mer <- ggplot() +
  theme_bw() +
  scale_x_continuous(limits = range(xtick), breaks = xtick, labels = xlat, expand = c(0, 0)) +
  scale_y_continuous(limits = range(ytick), breaks = ytick, labels= ylat, expand = c(0, 0)) + 
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), )+
  coord_fixed(ratio = 1)

# central meridian --------------------------------------------------------
# see Figure_3_uCT.R -> 3 axes

# choose central lines with 3D plot, all 0s if aligned to medcol, eq+chiasm
clv <- 0
clh <- 0
clp <- 0
clq <- 0

# 2023
ind_axis <- vaxis_gen(0, ixy = lens_ixy)
xy <- data.frame(ucl_rot_Mo_right)[ind_axis,] # full right
cmer_right <- xy[order(xy$yM),]

ind_axis <- vaxis_gen(0, ixy = ind_xy)
xy <- data.frame(ucl_rot_Mo)[ind_axis,] #matched right
cmer <- xy[order(xy$yM),]

# plt <- plt + geom_path(data = cmer, aes(x=xM, y=yM), colour = 'grey50', lwd=1)
# xy <- data.frame(ucl_rot_Merc)[match(ind_axis, eyemap[,2]),]
xy <- data.frame(ucl_rot_Merc)[ind_axis,]
cmer_Merc <- xy[order(xy$y),]

ind_axis <- vaxis_gen(7, ixy = lens_ixy)
xy <- data.frame(ucl_rot_Mo)[ind_axis,]
mer10_right <- xy[order(xy$yM),]

ind_axis <- vaxis_gen(9, ixy = ind_xy)
xy <- data.frame(ucl_rot_Mo)[ind_axis,]
mer10 <- xy[order(xy$yM),]

# # 2021
# ind_axis <- vaxis_gen(-1) # -11 midline
# xy <- data.frame(ucl_rot_Mo)[match(ind_axis, eyemap[,2]),]
# cmer <- xy[order(xy$yM),]
# # plt <- plt + geom_path(data = cmer, aes(x=xM, y=yM), colour = 'grey50', lwd=1)
# xy <- data.frame(ucl_rot_Merc)[match(ind_axis, eyemap[,2]),]
# cmer_Merc <- xy[order(xy$y),]
# 
# ind_axis <- vaxis_gen(9)
# xy <- data.frame(ucl_rot_Mo)[match(ind_axis, eyemap[,2]),]
# mer10 <- xy[order(xy$yM),]

