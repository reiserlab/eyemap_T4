# some useful functions


# -- (10) -- 2D geographic projection  -------------------------------------------------------------------------------------------

cart2sph2tp <- function(xyz) {
  xyz <- as.matrix(xyz)
  colnames(xyz) <- c('x','y','z')
  xyz %<>% as_tibble() %>%  
    mutate(y = -y) %>% #change the view to inside-out
    mutate(theta = acos(z), ang = x/sin(acos(z))) %>%
    mutate(ang = if_else((1-abs(ang)) < .Machine$double.eps, 1, ang)) %>%
    mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(ang)) %>%
    select(x,y,z,theta, phi) %>%
    mutate(t = theta / pi * 180, p = phi/pi*180) %>%  
    mutate(p = if_else(p > 180, p - 360, p)) %>% #move to [-pi pi]
    as.data.frame()
  return(xyz)
}


# - polar to Mollweide, input thetaphi is N x 2 [0<=theta<=180, -180<=phi<=180] in degree, output [x,y]
# https://mathworld.wolfram.com/MollweideProjection.html
Mollweide <- function(thetaphi){
  N0 <- nrow(thetaphi)
  xy0 <- matrix(ncol = 2, nrow = N0)
  # deal with NA
  ii_NA <- is.na(thetaphi[,1]) | is.na(thetaphi[,2])
  thetaphi <- thetaphi[!ii_NA,]
  
  N <- nrow(thetaphi)
  lambda <- thetaphi[,2]/180*pi #longitude
  # phi <- (thetaphi[,1] - 90)/180*pi #latitude
  phi <- (90 - thetaphi[,1])/180*pi #latitude, NEW 2020-12-20
  xy <- matrix(ncol = 2, nrow = N)
  for (j in 1:N) {
    theta <- asin(2*phi[j]/pi) #initial guess
    if (abs(abs(theta) - pi/2) < 1e-3) {
      xy[j,] <- c(2*sqrt(2)/pi*lambda[j]*cos(theta), sqrt(2)*sin(theta))
    } else {
      dtheta <- 1
      while (dtheta > 1e-3) {
        theta_new <- theta - (2*theta + sin(2*theta) - pi*sin(phi[j])) / (2 + 2*cos(2*theta))
        dtheta <- abs(theta_new - theta)
        theta <- theta_new
      }
      xy[j,] <- c(2*sqrt(2)/pi*lambda[j]*cos(theta), sqrt(2)*sin(theta))
    }
  }
  xy0[!ii_NA,] <- xy
  return(xy0)
}

# - Mercator
#' Mercator projection with aninside-out flip
#'
#' @param xyz Nx3 matrix
#'
#' @return Nx2 [x, y] in Mercator coord
#' @export
#'
#' @examples https://mathworld.wolfram.com/MercatorProjection.html
cart2Mercator <- function(xyz){
  xyz[,2] <- -xyz[,2] # left-right flip
  thetaphi <- cart2sphZ(xyz)[,2:3,drop=F]
  azimelev <- data.frame(azim = thetaphi[,2],
                         elev = pi/2 - thetaphi[,1]) %>%
    mutate(azim = ifelse(azim > pi, azim - 2*pi, azim)) %>%
    transmute(x = azim, y = log(tan(pi/4 + elev/2)) )  # y is not linear
  return(azimelev)
}

# -- (12) -- polar <-> cartesian ----------------------------------------------------------------------------------

# theta is from z-axis
cart2sphZ <- function(xyz){
  if (is.null(dim(xyz)) | dim(xyz)[2] != 3) {
    stop("wrong dim or not matrix")
  }
  xyz <- as.matrix(xyz)
  # if (length(xyz) == 3) {
  #   xyz <- matrix(xyz, ncol = 3)
  # }
  r <- sqrt(rowSums(xyz^2))
  x <- xyz[,1] / r
  y <- xyz[,2] / r
  z <- xyz[,3] / r
  theta <- acos(z)
  # val_acos <- if_else(1- x/sin(theta) < .Machine$double.eps, 1, x/sin(theta))
  val_acos <- x/sin(theta)
  val_acos[1 - x/sin(theta) < .Machine$double.eps] <- 1
  val_acos[x/sin(theta) - (-1) < .Machine$double.eps] <- -1
  phi <- 2*pi*(y < 0) + (-1)^(y < 0)*acos(val_acos)
  
  return(cbind(r, theta, phi)) #theta [0,pi], phi [0, 2*pi]
}

sph2cartZ <- function(rtp){
  if (is.null(dim(rtp)) | dim(rtp)[2] != 3) {
    stop("wrong dim or not matrix")
  }
  rtp <- as.matrix(rtp)
  # if (length(rtp) == 3) {
  #   trp <- matrix(trp, ncol = 3)
  # }
  r <- rtp[,1]
  t <- rtp[,2]
  p <- rtp[,3]
  x <- r * cos(p) * sin(t)
  y <- r * sin(p) * sin(t)
  z <- r * cos(t)
  
  return(cbind(x,y,z))
}

sph2elaz <- function(rtp){
  rtp <- as.matrix(rtp)
  elev <- 90 - rtp[,2] / pi * 180 #theta to elev, elev=0 for eq, elev=90 for north pole
  azim <- rtp[,3]/pi*180
  azim <- if_else(azim >= 180, azim - 360, azim) # left side is +y in sph, same left side is positive for azim
  
  return(cbind(elev, azim))
}

elaz2sph <- function(elaz){
  elaz <- as.matrix(elaz)
  theta <- (90 - elaz[,1] ) / 180 * pi 
  azim <- if_else(elaz[,2] < 0, elaz[,2] + 360, elaz[,2])
  phi <- azim /180 *pi
  
  
  return(cbind(1, theta, phi))
}

# -- (eyemap 1) ---vertical axis -----------------------------------------------------------

# add var ixy to distinguish lens_ixy vs ind_xy, 2023
# load('data/lens_ixy.RData') # 2021
# 2021, ind going dorsal -> ventral, vaxis_gen(-1) goes thru center ucl_rot_sm[361,]

vaxis_gen <- function(x0, ixy = lens_ixy){
  dx <- 1
  y0 <- 0; dy <- 1
  # positive direction
  nb_ixy <- matrix(ncol = 3, nrow = 0)
  cc <- 0
  repeat{
    nb <- ixy
    # nb <-  nb[nb[,2] %in% (x0 + dx*cc), ,drop=F]
    # nb <-  nb[nb[,3] %in% (y0 + dy*cc), ,drop=F]
    nb <-  nb[(nb[,2] %in% (x0 + dx*cc)) & (nb[,3] %in% (y0 + dy*cc)), ,drop=F]
    if (length(nb) == 0) {
      break
    }
    nb_ixy <- rbind(nb_ixy, nb)
    cc <- cc + 1
  }
  ind_axis <- rev(nb_ixy[,1])
  # negative
  nb_ixy <- matrix(ncol = 3, nrow = 0)
  cc <- 1
  repeat{
    nb <- ixy
    # nb <-  nb[nb[,2] %in% (x0 - dx*cc), ,drop=F ]
    # nb <-  nb[nb[,3] %in% (y0 - dy*cc), ,drop=F ]
    nb <-  nb[(nb[,2] %in% (x0 - dx*cc)) & (nb[,3] %in% (y0 - dy*cc)), ,drop=F]
    if (length(nb) == 0) {
      break
    }
    nb_ixy <- rbind(nb_ixy, nb)
    cc <- cc + 1
  }
  ind_axis <- c(ind_axis, nb_ixy[,1])
  
  # # rid off omitted pt, 2021 only
  # ind_axis <- ind_axis[ind_axis %in% rownames(ucl_rot_sm)]
  
  return(ind_axis)
}

# -- (eyemap 2) -- horizontal axis -----------------------------------------------------------
# front -> back, ie. left->right on 2D
haxis_gen <- function(v0, ixy = lens_ixy){
  # # 2021
  # if (v0 %% 2 == 0) {
  #   ii <- match(361, vaxis_gen(-1)) # start from the central vertical axis
  #   ixy <- lens_ixy[match(vaxis_gen(-1)[ii - (v0/2)], lens_ixy[,1]), ]
  # }
  # if (v0 %% 2 != 0) {
  #   ii <- match(388, vaxis_gen(0)) # start from the central vertical axis
  #   ixy <- lens_ixy[match(vaxis_gen(0)[ii - (v0+1)/2], lens_ixy[,1]), ]
  # }
  
  # 2023
  if (v0 %% 2 == 0) {
    ii <- match(ixy[ixy[,2]==0 & ixy[,3]==0, 1], vaxis_gen(0, ixy)) # start from the central vertical axis
    ixy2 <- ixy[match(vaxis_gen(0, ixy)[ii - (v0/2)], ixy[,1]), ]
  }
  if (v0 %% 2 != 0) {
    ii <- match(ixy[ixy[,2]== -1 & ixy[,3]==0, 1], vaxis_gen(-1, ixy)) # start from the central vertical axis
    ixy2 <- ixy[match(vaxis_gen(-1, ixy)[ii - (v0+1)/2], ixy[,1]), ]
  }
  
  x0 <- ixy2[2]
  y0 <- ixy2[3]
  dx <- -1 
  dy <- 1
  # positive direction
  nb_ixy <- matrix(ncol = 3, nrow = 0)
  cc <- 0
  repeat{
    nb <- ixy
    nb <-  nb[nb[,2] %in% (x0 + dx*cc), ,drop=F]
    nb <-  nb[nb[,3] %in% (y0 + dy*cc), ,drop=F ]
    if (length(nb) == 0) {
      break
    }
    nb_ixy <- rbind(nb_ixy, nb)
    cc <- cc + 1
  }
  ind_axis <- rev(nb_ixy[,1])
  # negative
  nb_ixy <- matrix(ncol = 3, nrow = 0)
  cc <- 1
  repeat{
    nb <- ixy
    nb <-  nb[nb[,2] %in% (x0 - dx*cc), ,drop=F ]
    nb <-  nb[nb[,3] %in% (y0 - dy*cc), ,drop=F ]
    if (length(nb) == 0) {
      break
    }
    nb_ixy <- rbind(nb_ixy, nb)
    cc <- cc + 1
  }
  ind_axis <- c(ind_axis, nb_ixy[,1])
  
  # # rid off omitted pt, 2021 only
  # ind_axis <- ind_axis[ind_axis %in% eyemap[,2]]
  
  return(ind_axis)
}

# --- (2) --- child_node ----------------------------------------------------------------------------------

# To extract all nodes of a neuron 'xneu'<nat neuron> downstream of a node(s) 'xroot'<data.frame>
# cp function 'subset'

child_node <- function(xneu, xroot){
  df_head <- xroot
  i_node <- logical(length = dim(xneu$d)[1])
  repeat{
    ii <- match(xneu$d[,"Parent"], df_head[, "PointNo"])
    ii[is.na(ii)] <- 0
    ii <- as.logical(ii)
    df_head <- xneu$d[ii,]
    i_node <- i_node | ii
    if (sum(ii) == 0) {
      break
    }
  }
  df_node <- xneu$d[i_node,]
  return(df_node)
}

# -- (7) -- 3D cross product; 3D point-line dist ----------------------------------------------------------------------------------------

# a, b are two 3D vectors
cross3D <- function(a, b){
  a <- as.matrix(a)
  b <- as.matrix(b)
  prod <- matrix(c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1]), ncol = 1)
  return(prod)
}


# p is a point, a,b are points on a line. Slow in a loop, alt use sweep()
dist_pl3D <- function(p, a, b){
  p <- as.matrix(p)
  a <- as.matrix(a)
  b <- as.matrix(b)
  vc <- as.vector(cross3D(p-a,b-a))
  return(sqrt(sum(vc^2)) / sqrt(sum((b-a)^2)))
}

# -- (9) -- return a quaternion rotation matrix -----------------------------------------------------------------------

# cf. https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix

# vr is the axis of rotation, ang is angle in deg
# position after rotation:  pt1 <- pt0 %*% t(R_mat), where R_mat <- quaternion3D(...)
quaternion3D <- function(vr, ang){
  ang <- ang / 180 * pi
  vr <- vr / sqrt(sum(vr^2))
  
  qr <- cos(ang/2)
  qi <- vr[1]*sin(ang/2)
  qj <- vr[2]*sin(ang/2)
  qk <- vr[3]*sin(ang/2)
  
  R <- matrix(c(
    1-2*(qj^2+qk^2), 2*(qi*qj-qk*qr), 2*(qi*qk+qj*qr),
    2*(qi*qj+qk*qr), 1-2*(qi^2+qk^2), 2*(qj*qk-qi*qr),
    2*(qi*qk-qj*qr), 2*(qj*qk+qi*qr), 1-2*(qi^2+qj^2)),
    ncol = 3, byrow = T)
  
  return(R)
}

#  --(13)-- rigid/Euclidean transform neurons by pc  ---------------------------------------------------------------------------
xEucl_neu <- function(n, Rm, Tm){
  if (!is.neuronlist(n)) {
    n <- neuronlist(n)
  }
  n_xform <- n
  for (j in 1:length(n)) {
    xyz <- n[[j]]$d[, c("X","Y","Z")]
    xyz_xform <- sweep(as.matrix(xyz), 2, Tm, '-') %*% Rm
    n_xform[[j]]$d[, c("X","Y","Z")] <- xyz_xform
  }
  return(n_xform)
}

#  angle from cos ----------------------------------------------------------

# single angle
angcos <- function(v,u){
  ang <- acos(v %*% u / sqrt(sum(v^2)) / sqrt(sum(u^2))) /pi*180
  return(ang)
}
# two sets of vectors
angcos_vf <- function(vf1, vf2){
  if (nrow(vf1) != nrow(vf2)) {
    stop("wrong input")
  }
  vf1 <- as.matrix(vf1)
  vf2 <- as.matrix(vf2)
  ang <- matrix(ncol = 1, nrow = nrow(vf1))
  for (j in 1:nrow(vf1)) {
    v1 <- vf1[j,]
    v2 <- vf2[j,]
    if ((!is.na(v1[1])) & (!is.na(v2[1]))) {
      ang[j] <- angcos(v1,v2)
    }
  }
  return(ang)
}

# PD field fit, 1/ brute sampling -----------------------------------------------------------------------------------------------
# rotation and translation flow generating function ---------------------------------------------------------------

# NB, va is motion axis now (2020 06 25)

# va -- axis of rotation, pt0 -- positions of field on a unit sphere, return a field of optic flow
R_gen <- function(va, pt0){
  # pt0 <- sweep(pt0, 1, sqrt(rowSums(pt0^2)), '/')
  if (sum(sqrt(rowSums(pt0^2)) != 1) != 0) {
    warning("point positions not normalized")
    pt0 <- sweep(pt0, 1, sqrt(rowSums(pt0^2)), '/')
  }
  va <- -va  # from motion axis to vector field axis
  
  ## ## quaternion rep, less accurate
  # angle of rot, will be normalized later, so a small value
  ang <- 0.1/180*pi
  R_mat <- quaternion3D(va, ang)
  # position after rotation
  pt1 <- pt0 %*% t(R_mat)
  # flow vec
  v <- pt1 - sweep(pt0, 1, rowSums(pt1 * pt0), '*')
  v <- v / max(sqrt(rowSums(v^2)))# norm
  
  
  # ## ## cross prod, slow
  # v <- t(apply(pt0, 1, function(x) cross3D(va, x)))
  # v <- v / max(sqrt(rowSums(v^2)))# norm
  
  return(v)
}


t_gen <- function(va, pt0){
  pt0 <- sweep(pt0, 1, sqrt(rowSums(pt0^2)), '/')
  va <- -va  # from motion axis to vector field axis
  
  # ## ##
  # # translation, will be normalized later, so a small value
  # ds <- 0.1
  # # position after translation
  # pt1 <- sweep(pt0, 2, ds*va, '+')
  # # flow vec
  # # v3D <- pt1 - pt0
  # v <- pt1 - sweep(pt0, 1, rowSums(pt1 * pt0), '*')
  # v <- v / max(sqrt(rowSums(v^2)))# norm
  
  ## ##
  # v <- t(apply(pt0, 1, function(x) va - x * c(va %*% x)) ) #slow
  
  v <- sweep(-sweep(pt0, 1, pt0 %*% va, '*'), 2, va, '+')
  # v <- va - sweep(pt0, 1, pt0 %*% va, '*') # WRONG !
  v <- v / max(sqrt(rowSums(v^2)))# normalize
  
  return(v)
}

# sample axis angles ----------------------------------------------------------------------------------------------

# da is the angular differential in degree, output [t(degreee), p, x,y,z] on unit S2
axis_sampling <- function(da = 2) {
  # da <- 1 #radius = 2
  arcL <- da* 1/180*pi # arc length of dist between samples
  maxis <- matrix(ncol = 2) # motion axis
  for (theta in seq(0,180,by = da)) {
    if (theta == 0 | theta == 180) {
      maxis <- rbind(maxis, c(theta, 0))
    } else {
      phiN <- ceiling( 2*pi*sin(theta/180*pi) / arcL )
      if (phiN %% 2 == 0) {
        maxis <- rbind(maxis, c(theta, 0))
        for (phi_ind in 1:(phiN/2-1)) {
          maxis <- rbind(maxis, c(theta, arcL*phi_ind/sin(theta/180*pi)/pi*180))
          maxis <- rbind(maxis, c(theta, 360 - arcL*phi_ind/sin(theta/180*pi)/pi*180))
        }
        maxis <- rbind(maxis, c(theta, 180))
      } else {
        maxis <- rbind(maxis, c(theta, 0))
        for (phi_ind in 1:((phiN-1)/2)) {
          maxis <- rbind(maxis, c(theta, arcL*phi_ind/sin(theta/180*pi)/pi*180))
          maxis <- rbind(maxis, c(theta, 360 - arcL*phi_ind/sin(theta/180*pi)/pi*180))
        }
      }
    }
  }
  maxis <- maxis[-1,]
  
  maxis_xyz <- matrix(ncol = 3, nrow = dim(maxis)[1])
  for (j in 1:dim(maxis)[1]) {
    maxis_xyz[j,] <- c(sin(maxis[j,1]/180*pi)*cos(maxis[j,2]/180*pi),
                       sin(maxis[j,1]/180*pi)*sin(maxis[j,2]/180*pi),
                       cos(maxis[j,1]/180*pi))
  }
  
  return(cbind(maxis,maxis_xyz))
}

# --- (5) mkpoly -----------------------------------------------------------------------------------------------------

# make a single polygon from ashape$edges, correct indices order and add connecting rows
# input a matrix (ashape$edges), output a matrix (x,y) as a list of polygon vertices

# # help with randomizing inline points, eg.
# ash <- ashape(simdata_ct[[j]]+matrix(runif(dim(simdata_ct[[j]])[1]*2, 1e-9, 2e-9), ncol = 2), alpha = 5)
# simdata_poly[[j]] <- mkpoly(ash$edges)[[1]][,3:4]


# - v2
mkpoly <- function(xy) {
  # --- first step, separate into connected groups
  # xy <- xy_pj[[j]]$ashape[[k]]$edges
  xy <- as.data.frame(xy)[,1:6]
  xy_2 <- xy
  # if (dim(xy)[2] < 6) {
  #   xy <- t(xy)
  # }
  
  # remove singlet
  ind_edge <- c(as.matrix(xy[, 1:2]))
  ind_u <- unique(ind_edge)
  ind_s <- ind_u[sapply(ind_u, function(x) sum(ind_edge %in% x)) == 1]
  while (length(ind_s) > 0) {
    ind_m <- match(ind_s[1], xy[,1]) #find this index
    if (!is.na(ind_m)) {
      xy <- xy[-ind_m,]
    } else {
      ind_m <- match(ind_s[1], xy[,2]) #maybe in the second column
      xy <- xy[-ind_m,]
    }
    # search singlet again
    ind_edge <- c(as.matrix(xy[, 1:2]))
    ind_u <- unique(ind_edge)
    ind_s <- ind_u[sapply(ind_u, function(x) sum(ind_edge %in% x)) == 1]
  }
  
  xyset <- list() # vertices for a list of polygons
  if (dim(xy)[1] > 2) {
    # xy is of [x0 y0 x1 y1]
    N <- 1
    xyset[[N]] <- xy[1,]
    xy <- xy[-1, ]
    ii <- c()
    while (dim(xy)[1] >= 1) {
      ii[1] <- match(tail(xyset[[N]], 1)[2], xy[, 1])
      ii[2] <- match(tail(xyset[[N]], 1)[2], xy[, 2])
      if (!is.na(ii[1])) {
        xyset[[N]] <- rbind(xyset[[N]], xy[ii[1], ])
        xy <- xy[-ii[1], ]
      } else if (!is.na(ii[2])){
        xytmp <- xy[ii[2], c(2,1,5,6,3,4)]
        colnames(xytmp) <- colnames(xyset[[N]])
        xyset[[N]] <- rbind(xyset[[N]], xytmp)
        xy <- xy[-ii[2], ]
      } else {
        # xyset[[N]] <- xyset[[N]][,3:4] 
        N <- N + 1
        xyset[[N]] <- xy[1, ]
        xy <- xy[-1, ]
      }
    }
  }
  return(xyset)
}

# build hex grid ----------------------------------------------------------

#' Title
#' @details The point set need to be mostly regular and convex for this to work well.
#' Auto-7pts selection assumes leveled head z=[0,0,1].
#' Use lefteye=T for left eye.
#' @param pt N x 3 point set
#' @param j7 the first 1 (center) + 6 (neighbors) points. Supersedes j0
#' @param j0 the first 1 (center) point.
#' @param lefteye=FALSE assume right eye
#'
#' @return a list of 3 variables, ind_nb contains the indices of all hexagons; 
#' ind_xy contains the indices and xy coord for each point; checked indices.
#' @export
#'
#' @examples

# pt2hex <- function(pt, j7 = c()){
reghex <- function(pt, j7 = c(), j0 = c(), lefteye=FALSE){
  
  pt <- as.matrix(pt)
  
  # data var
  ind_nb <- matrix(ncol = 7, nrow = nrow(pt)) # ind of self and 6 nb in [p,v,q,-p,-v,-q]
  ind_xy <- matrix(ncol = 3, nrow = nrow(pt)) # [ind, x, y], use integer coord
  dv <- c(1, 1) # delta v
  dp <- c(1, 0)
  dq <- c(0, 1)
  dhex <- rbind(dp,dv,dq,-dp,-dv,-dq)
  
  ind_done <- c() # done checking
  ind_cur <- c() # need to check
  ind_new <- c() # save for next round
  ind_all <- seq(1, nrow(pt))
  
  # is j7 provided
  if (length(j7) == 0) {
    # to get the first hex
    if (length(j0) == 0) {
      j0 <- rowSums(sweep(pt,2,colMeans(pt))^2) %>% which.min()
    }
    zout <- pt[j0,] - colMeans(pt) # assume convex
    j7 <- rowSums(sweep(pt, 2, pt[j0,])^2) %>% order() %>% head(7)
    pt7 <- pt[j7,]
    
    # origin
    ind_nb[1,1] <- j7[1]
    ind_xy[1,] <- c(j7[1], 0, 0)
    # use vn = c(0,0,1) and zout to assign the other 6 pts
    vv <- sweep(pt7[-1,], 2, pt7[1,])
    # vv %*% vn
    iimax <- which.max(vv[,3])
    ind_nb[1,1+2] <- j7[iimax+1]
    ind_xy[1+2,] <- c(j7[iimax+1], dhex[2,])
    iip <- vv[,3] > 0
    iip[iimax] <- FALSE
    iip <- which(iip) + 1 #ind in j7
    if (length(iip) != 2) {
      print("iip error")
    }
    #iip[1] is +p for right eye
    if (zout %*% cross3D(vv[iimax,], vv[iip[1]-1,]) * (-1)^lefteye > 0) {
      ind_nb[1,1+1] <- j7[iip[1]]
      ind_xy[1+1,] <- c(j7[iip[1]], dhex[1,])
      ind_nb[1,1+3] <- j7[iip[2]]
      ind_xy[1+3,] <- c(j7[iip[2]], dhex[3,])
    } else { #iip[1] is +q
      ind_nb[1,1+1] <- j7[iip[2]]
      ind_xy[1+1,] <- c(j7[iip[2]], dhex[1,])
      ind_nb[1,1+3] <- j7[iip[1]]
      ind_xy[1+3,] <- c(j7[iip[1]], dhex[3,])
    }
    
    iimin <- which.min(vv[,3])
    ind_nb[1,1+5] <- j7[iimin+1]
    ind_xy[1+5,] <- c(j7[iimin+1], dhex[5,])
    iim <- vv[,3] < 0
    iim[iimin] <- FALSE
    iim <- which(iim) + 1
    if (length(iim) != 2) {
      print("iim error")
    }
    if (zout %*% cross3D(vv[iimin,], vv[iim[1]-1,]) * (-1)^lefteye > 0) {
      ind_nb[1,1+4] <- j7[iim[1]]
      ind_xy[1+4,] <- c(j7[iim[1]], dhex[4,])
      ind_nb[1,1+6] <- j7[iim[2]]
      ind_xy[1+6,] <- c(j7[iim[2]], dhex[6,])
    } else {
      ind_nb[1,1+4] <- j7[iim[2]]
      ind_xy[1+4,] <- c(j7[iim[2]], dhex[4,])
      ind_nb[1,1+6] <- j7[iim[1]]
      ind_xy[1+6,] <- c(j7[iim[1]], dhex[6,])
    }
    
  } else { #if j7 provided
    ind_nb[1,] <- j7
    ind_xy[1,] <- c(j7[1], 0, 0)
    ind_xy[2:7,] <- cbind(j7[-1], dhex)
  }
  
  # main loop
  ind_done <- j7[1]
  ind_cur <- ind_nb[1,-1]
  rnb <- 2 # row number for ind_nb and ind_xy
  rxy <- 8
  while (length(ind_cur) > 0) {
    for (j in ind_cur) {
      ixy <- ind_xy[match(j, ind_xy[,1]), ]
      ind_nb[rnb, 1] <- ixy[1] # home ind
      for (k in 1:3) { #check pairs of nbs in order of p,v,q
        xy <- ixy[2:3] + dhex[k,] #get expected xy
        ii1 <- which(ind_xy[,2] %in% xy[1] & ind_xy[,3] %in% xy[2]) #?exist nb
        xy <- ixy[2:3] - dhex[k,]
        ii2 <- which(ind_xy[,2] %in% xy[1] & ind_xy[,3] %in% xy[2])
        if (length(ii1) + length(ii2) == 2) { #already in ind_xy
          # ind_nb[rnb, 2*k+(0:1)] <- c(ind_xy[ii1,1], ind_xy[ii2,1]) # add to ind_nb, [+p,-p,+v,-v,+q,-q]
          ind_nb[rnb, k+ c(1,4)] <- c(ind_xy[ii1,1], ind_xy[ii2,1]) # add to ind_nb, in dhex order
        } 
        if (length(ii1) + length(ii2) == 1) { #one is missing
          ii <- ifelse(length(ii1) == 1, ii1, ii2)
          # 19 nbs for pca and local search
          ii19 <- sweep(pt, 2, pt[ixy[1],] )^2 %>% rowSums() %>% order() %>% head(1+8)
          pt19 <- pt[ii19, ]
          # PCA 
          pc19 <- prcomp(pt19)
          
          pt_pca <- sweep(pt, 2, pc19$center) %*% pc19$rotation
          
          vv <- pt_pca[ixy[1], 1:2] - pt_pca[ind_xy[ii,1], 1:2] 
          thr <- sqrt(sum(vv^2)) #dist
          xyz <- pt_pca[ixy[1], 1:2] + vv #expected position
          dd <- sweep(pt_pca[ii19, 1:2], 2, xyz)^2 %>% rowSums() %>% sqrt()
          ii3 <- ii19[which.min(dd)] #nearest
          if (min(dd) < thr*0.4 & !(ii3 %in% ind_xy[,1])) { # not too far off AND new 
            if (identical(ii, ii1)) {
              # ind_nb[rnb, 2*k+(0:1)] <- c(ind_xy[ii1,1], ii3)
              ind_nb[rnb, k+ c(1,4)] <- c(ind_xy[ii1,1], ii3)
              ind_xy[rxy,] <- c(ii3, ixy[2:3] - dhex[k,])
            } else {
              # ind_nb[rnb, 2*k+(0:1)] <- c(ii3, ind_xy[ii2,1])
              ind_nb[rnb, k+ c(1,4)] <- c(ii3, ind_xy[ii2,1])
              ind_xy[rxy,] <- c(ii3, ixy[2:3] + dhex[k,])
            }
            rxy <- rxy + 1
            ind_new <- c(ind_new, ii3)
          } else { #only register the existing nb
            if (identical(ii, ii1)) {
              # ind_nb[rnb, 2*k+0] <- ind_xy[ii1,1]
              ind_nb[rnb, k+1] <- ind_xy[ii1,1]
            } else {
              # ind_nb[rnb, 2*k+1] <- ind_xy[ii2,1]
              ind_nb[rnb, k+4] <- ind_xy[ii2,1]
            }
          }
        }
      }
      rnb <- rnb + 1
    }
    ind_done <- c(ind_done, ind_cur)
    ind_cur <- ind_new
    ind_new <- c()
  }
  
  return(list(ind_nb, ind_xy, ind_done))
}


# -- (6) -- rotate a point set ------------------------------------------------------------------------------------
# update 2022-01-28, (1)zv, xv is 1x3. (2) orig=c(0,0,0)
# rotate a point set (eg T4 dendrite nodes) assuming a plane whose normal as z-axis (eg. plane normal to pc3) and +x direction
# input: Nx3 point set - pts, 1X3 norm vec - zv, 1x3 +x dir - xv
# zv[1:3] is the origin, xv is corrected to be perpendicular to zv if not already
# output: dataframe [x,y,z] 

rot_zx0 <- function(pts, zv, xv, orig=c(0,0,0)){
  
  zv <- as.numeric(zv)
  xv <- as.numeric(xv)
  if (length(zv) != 3 | length(xv) != 3) {
    stop("zv and xv are 1x3")
  }
  
  zvu <- zv/sqrt(sum(zv^2))
  xvu <- xv/sqrt(sum(xv^2))
  if (zvu %*% xvu != 0) {
    warning("xv not perpendicular to zv, re-defind xv")
    xvu <- xvu - c(xvu %*% zvu) * zvu
    xvu <- xvu/sqrt(sum(xvu^2))
  }
  
  # find the rotation that makes zv3 z-axis and xv3 x-axis
  # first rotation about z to align zv3 to xz plane
  t1 <- 2*pi*(zvu[2]<0) + (-1)^(zvu[2]<0)*acos(zvu[1]/sqrt(sum(zvu[c(1,2)]^2))) #[0 2pi]
  t1 <- -t1
  R1 <- matrix(c(cos(t1), -sin(t1), 0, sin(t1), cos(t1), 0, 0, 0, 1), ncol = 3, byrow = T)
  zvu1 <- R1 %*% zvu # = R1 %*% matrix(zvu, nrow=3)
  
  # second, rotate about y to align zvu to z
  t2 <- 2*pi*(zvu1[1]<0) + (-1)^(zvu1[1]<0)*acos(zvu1[3]/sqrt(sum(zvu1[c(3,1)]^2)))
  t2 <- -t2
  R2 <- matrix(c(cos(t2), 0, sin(t2), 0, 1, 0, -sin(t2), 0, cos(t2)), ncol = 3, byrow = T)
  zvu2 <- R2 %*% zvu1
  xvu2 <- R2 %*% R1 %*% xvu
  
  # third rotation about z to align central point on x
  t3 <- 2*pi*(xvu2[2]<0) + (-1)^(xvu2[2]<0)*acos(xvu2[1]/sqrt(sum(xvu2[c(1,2)]^2)))
  t3 <- -t3
  R3 <- matrix(c(cos(t3), -sin(t3), 0, sin(t3), cos(t3), 0, 0, 0, 1), ncol = 3, byrow = T)
  xvu3 <- R3 %*% xvu2
  
  # apply to all points
  pts0 <- sweep(pts, 2, orig) # vec of pts wrt the orig
  
  pts_rot <- t(R3 %*% R2 %*% R1 %*% t(pts0))
  pts_rot <- as.data.frame(pts_rot)
  colnames(pts_rot) <- c("X","Y","Z")
  
  return(pts_rot)
}


# process Eyal's data -----------------------------------------------------
# return tb_v = [theta(base), phi(base), spk, spk, sub, sub, dark=0/bright=1]

Eyal_arena_2023 <- function(tb){
  # arena radius
  A <- (4+10/32 + 4+20/32)/2  
  # area center wrt to fly
  Dx <- (5+6/32) - A  
  Dy <- (4+20/32) - A
  Dz <- - (1+ 13/32)
  
  # LED angular size [inch]
  DL <- A * 1.25/180*pi
  
  # response amp, 
  resfac <- 12 # st. resfac * largest DSI ~ 10 pixel
  
  # - loop
  # [elev azim] for base, spk, subthr
  tb_v <- matrix(ncol = 2+2+2, nrow = nrow(tb)) 
  for (j in 1:nrow(tb)) {
    # angles
    t1 <- tb$meanThetaSpk[j] #mean theta spike
    t2 <- tb$meanThetaSub[j] #subthr
    tt <- c(0, t1, t2) #base, head(spk), head(sub)
    
    # amplitudes
    DSI1 <- tb$meanDSISpk[j] # spike
    DSI2 <- tb$meanDSISub[j] #subthr
    DD <- c(0, DSI1, DSI2)
    # DD <- c(0, 1, 1) # discard amplitude = DSI
    
    # azim position of LED in [LED count]
    nx <- tb$stimPosX[j] 
    # elev position of LED
    ny <- tb$stimPosY[j] 
    # arena rotation ang
    arot <- tb$arenaAng[j] /180*pi
    # equator below holder
    afly <- tb$headAng[j] /180*pi 
    
    if (!is.na(afly)) {
      for (k in 1:3) {
        # -- in arena 2D coord
        # 96 pixel * 1.25 = 120 deg per side
        # angle from midline
        ax <- (120 - nx*1.25) /180*pi - cos(tt[k])* DL/A * DD[k] * resfac
        # dist from bottom
        ay <- DL * ny + sin(tt[k]) * DL * DD[k] * resfac
        
        # -- go to lab coord via a translation
        # origin at arena center, xy-plane is level
        x <- A * cos(ax)
        y <- A * sin(ax)
        z <- ay - (4 + 24/32)/2
        
        xyz <- matrix(c(x,y,z), ncol=3)
        
        # -- arena rotation, 
        # arot is angle below ground, use -arot
        xyz <- matrix(c(cos(arot), 0, -sin(-arot),
                        0, 1, 0,
                        sin(-arot), 0, cos(arot)), ncol=3, byrow = T) %*% t(xyz)
        xyz <- t(xyz)
        
        # -- go to fly's coord
        xyz <- xyz + c(Dx, Dy, Dz)
        
        # -- fly rotation 
        # afly is angle below ground, use afly since it's ref that's rotating
        xyz <- matrix(c(cos(afly), 0, -sin(afly),
                        0, 1, 0, 
                        sin(afly), 0, cos(afly)), ncol = 3, byrow = T) %*% t(xyz)
        xyz <- t(xyz)
        
        # to spherical coord
        thetaphi <- cart2sphZ(xyz) /pi*180
        
        tb_v[j, (k-1)*2+(1:2)] <- round(c(90 - thetaphi[2], thetaphi[3]), 2)
      }
    }
  }
  tb_v <- cbind(tb_v, tb$edgeVal)
  
  return(tb_v)
}

# compute radius of curvature from a circle fit -----------------

roc_fit <- function(xyz){
  # flatten to 2D
  pca <- prcomp(xyz)
  pc3 <- pca$rotation[,3]
  com <- pca$center
  xyz_pc <- sweep(xyz,2,com) %*% pca$rotation
  xy_pc <- xyz_pc[, 1:2]
  # roc
  xyz_data <- xy_pc
  colnames(xyz_data) <- c('x','y')
  xyz_data2 <- xyz_data^2
  Y <- rowSums(xyz_data2)
  X <- cbind(2*xyz_data, rep(1,dim(xyz_data)[1]))
  cir <- lsfit(X,Y,intercept = FALSE)
  r_fit <- sqrt(cir[[1]][3]+sum(cir[[1]][1:2]^2)) #radius cir[[1]] == cir[["coefficients"]]
  
  return(r_fit)
}


# # (1) neu_dir -----------------------------------------------------------------------------------------------------
# 
# # To calculate the orientation of a subset of a neuron, several methods
# # Reroot the subset at the given tag/node, which is also the tail of the output vector 
# # input: xneu - given neuron, xroot - reroot with 'subset', node (int) or tag (str)
# # output: 1x6 vector
# 
# L <- 5000 #normalization factor ~ inter-ommatidium distance
# 
# # -- method 1
# # Calculate for each segment between branch points and end points, then take the mean
# neu_dir_seg <- function(xneu, xroot){
#   # check if node or tag
#   if (is.character(xroot)) {
#     ind_R = match(xneu$tags[[xroot]] , xneu$d$PointNo)
#     node_R = xneu$d[ind_R,]
#     xyz_R <- xyzmatrix(node_R)
#   } else if (is.numeric(xroot)) {
#     ind_R = xroot
#     node_R = xneu$d[ind_R,]
#     xyz_R <- xyzmatrix(node_R)
#   } else {
#     paste("wrong terminal selector")
#     return()
#   }
#   
#   # get the subset, reroot to xroot
#   xneu_g <- as.ngraph(xneu)
#   sub_nodes <- igraph::graph.dfs(xneu_g, root = ind_R, unreachable=FALSE, neimode='out')$order
#   sub_tree <- subset(xneu, sub_nodes)
#   # nopen3d()
#   # plot3d(sub_tree, add=TRUE, col='red', lwd=2)
#   
#   ap <- sub_tree$d # all nodes
#   ep <- sub_tree$d[sub_tree$EndPoints, ] #end nodes
#   bp <- sub_tree$d[sub_tree$BranchPoints, ] #branch nodes
#   
#   # get segments starting from ep
#   c_node <- ep #child node
#   p_node <- rbind(bp, node_R) #parent node
#   a_node <- ap # all nodes
#   seg_ends <- matrix(ncol = 8, nrow = dim(c_node)[1]) #list of segments
#   # colnames(seg_ends) <- c("PointNo","X","Y","Z","PointNo","X","Y","Z")
#   for (k in 1:dim(c_node)[1]) {
#     ii_node <- c_node[k,]
#     repeat{
#       ii <- match(ii_node$Parent, a_node$PointNo)
#       ii_node <- a_node[ii,]
#       isfound <- match(ii_node$PointNo, p_node$PointNo)
#       if (!is.na(isfound)) {
#         seg_ends[k,] <- unlist(c(c_node[k,c("PointNo","X","Y","Z")], p_node[isfound,c("PointNo","X","Y","Z")]))
#         break
#       }
#     }
#   }
#   seg_eb <- cbind(seg_ends, sqrt(rowSums((seg_ends[,6:8]-seg_ends[,2:4])^2))) #end to branch
#   
#   # get segments starting from bp
#   c_node <- bp #child node
#   p_node <- rbind(bp, node_R) #parent node
#   a_node <- ap # all nodes
#   seg_ends <- matrix(ncol = 8, nrow = dim(c_node)[1]) #list of segments
#   # colnames(seg_ends) <- c("PointNo","X","Y","Z","PointNo","X","Y","Z")
#   for (k in 1:dim(c_node)[1]) {
#     ii_node <- c_node[k,]
#     repeat{
#       ii <- match(ii_node$Parent, a_node$PointNo)
#       ii_node <- a_node[ii,]
#       isfound <- match(ii_node$PointNo, p_node$PointNo)
#       isfound
#       if (!is.na(isfound)) {
#         seg_ends[k,] <- unlist(c(c_node[k,c("PointNo","X","Y","Z")], p_node[isfound,c("PointNo","X","Y","Z")]))
#         break
#       }
#     }
#   }
#   seg_bb <- cbind(seg_ends, sqrt(rowSums((seg_ends[,6:8]-seg_ends[,2:4])^2))) #branch to branch
#   
#   # combine all seg and apply cutoff
#   v_seg <- rbind(seg_eb, seg_bb)
#   v_seg <- v_seg[order(v_seg[,9]),]
#   
#   # # cutoff length
#   # v_seg <- v_seg[v_seg[,9] > quantile(v_seg[,9], 0.3),]
#   
#   # # filter out seg with orientation opposite to seg-base
#   # v_seg_bool <- matrix(0, nrow = dim(v_seg)[1], ncol = 1)
#   # for (k in 1:dim(v_seg)[1]) {
#   #   v_s <- v_seg[k,2:4] - v_seg[k,6:8]
#   #   v_b <- xyz_D_avg - v_seg[k,6:8]
#   #   v_seg_bool[k] <- (sum(v_s*v_b) < 0)
#   # }
#   # v_seg <- v_seg[as.logical(v_seg_bool),]
#   
#   v_seg_avg <- colSums(v_seg[,2:4]-v_seg[,6:8])
#   v_seg_norm <- v_seg_avg/sqrt(sum(v_seg_avg^2)) * L
#   
#   
#   v_dir <- c(xyz_R, xyz_R + v_seg_norm)
#   return(v_dir)
# }
# 
# 
# # -- method 2
# # use Strahler order, eg. set so = c(2,3)
# # scaling -> multiply the vector by this value
# 
# neu_dir_SO <- function(xneu, xroot, so, sumType = c('avg', 'all', 'so_avg'), scaling = 1){
#   sumType <- match.arg(sumType)
#   # check if node or tag
#   if (is.character(xroot)) {
#     ind_R = match(xneu$tags[[xroot]] , xneu$d$PointNo)
#     node_R = xneu$d[ind_R,]
#     xyz_R <- colMeans(xyzmatrix(node_R))
#   } else if (is.numeric(xroot)) {
#     ind_R = xroot
#     node_R = xneu$d[ind_R,]
#     xyz_R <- colMeans(xyzmatrix(node_R))
#   } else {
#     paste("wrong terminal selector")
#     return()
#   }
#   
#   # TODO, multiple roots
#   # get the subset, reroot to xroot
#   xneu_g <- as.ngraph(xneu)
#   sub_nodes <- igraph::graph.dfs(xneu_g, root = ind_R, unreachable=FALSE, neimode='out')$order
#   sub_tree <- subset(xneu, sub_nodes)
#   sub_tree <- resample(sub_tree, stepsize = 100)
#   # nopen3d()
#   # plot3d(sub_tree, add=TRUE, col='red', lwd=2)
#   
#   # calculate Strahler order
#   taro <- sub_tree
#   taro_so <- strahler_order(taro)
#   
#   # get vector from segments of a given so
#   seg_v_so <- matrix(ncol = 3, nrow = length(so)) 
#   ii_row <- 1
#   for (j in so) {
#     seg_ls <- taro$SegList[taro_so$segments == j]
#     vec <- matrix(0, ncol = 3, nrow = 1)
#     for (k in 1:length(seg_ls)) {
#       node_ii <- seg_ls[[k]]
#       if (length(node_ii) > 1) {
#         vec <- vec + (taro$d[tail(node_ii,1), c("X","Y","Z")] - taro$d[head(node_ii,1), c("X","Y","Z")])
#       }
#     }
#     if (sumType == 'so_avg') { #order avg
#       seg_v_so[ii_row,] <- as.numeric(vec)/length(seg_ls) 
#     } else {
#       seg_v_so[ii_row,] <- as.numeric(vec)
#     }
#     ii_row <- ii_row + 1
#   }
#   
#   vso <- colSums(seg_v_so)
#   if (sumType == 'avg' | sumType == 'so_avg') {
#     # vso <- colMeans(seg_v_so)
#     vso_norm <- vso / sqrt(sum(vso^2))
#     v_dir <- c(xyz_R, xyz_R + vso_norm * scaling)
#     return(v_dir)
#   } 
#   if (sumType == 'all') {
#     v_dir <- c(xyz_R, xyz_R + vso)
#     return(v_dir)
#   }
# }
# 
# 
# # -- method 3, multiple branches
# # use Strahler order, eg. set so = c(2,3), 
# # scaling -> multiply the vector by this value
# neu_dir_SO_mul <- function(xneu, xroot, so, sumType = c('avg', 'all'), scaling = 1){
#   sumType <- match.arg(sumType)
#   # check if node or tag
#   if (is.character(xroot)) {
#     ind_R = match(xneu$tags[[xroot]] , xneu$d$PointNo)
#     node_R = xneu$d[ind_R,]
#     xyz_R <- colMeans(xyzmatrix(node_R))
#   } else if (is.numeric(xroot)) {
#     ind_R = xroot
#     node_R = xneu$d[ind_R,]
#     xyz_R <- colMeans(xyzmatrix(node_R))
#   } else {
#     paste("wrong terminal selector")
#     return()
#   }
#   
#   xneu_g <- as.ngraph(xneu)
#   
#   # - multiple branches
#   # -- get so for each branch
#   sub_tree_br <- list() #sub_tree(s) 
#   sub_tree_br_g <- list()
#   sub_tree_br_so <- list() # sub_tree(s) SO list
#   so_max_br <- c()
#   for (ii_br in 1:length(ind_R)) {
#     sub_points_br <- igraph::graph.dfs(xneu_g, root = ind_R[ii_br], unreachable=FALSE, neimode='out')$order
#     sub_points <- c(sub_points, sub_points_br) # collect all subtree nodes
#     sub_tree_br[[ii_br]] <- subset(tar, sub_points_br)
#     sub_tree_br_so[[ii_br]] <- strahler_order(sub_tree_br[[ii_br]]) # Strahler order
#     sub_tree_br_g[[ii_br]] <- as.ngraph(sub_tree_br[[ii_br]], weights = T)
#     so_max_br <- c(so_max_br, max(sub_tree_br_so[[ii_br]]$segments) )
#   }
# 
#   # -- define orientation based on strahler order 
#   so_L <- matrix(0, nrow = max(so), ncol = length(ind_R)) # cable length of SO branches
#   so_N <- matrix(0, nrow = max(so), ncol = length(ind_R)) #segment number
#   so_v_mat <- matrix(nrow = max(so), ncol = 3) #normalized dir vec
#   so_v_mat_raw <- matrix(nrow = max(so), ncol = 3) #in raw unit
#   
#   for (ii_so in so) {
#     seg_v <- matrix(ncol = 3, nrow = 0)
#     for (ii_br in 1:length(ii_root)) {
#       seg_ls <- sub_tree_br[[ii_br]]$SegList[sub_tree_br_so[[ii_br]]$segments == ii_so]
#       # add up so segments
#       if (length(seg_ls) > 0) {
#         for (ii_seg in 1:length(seg_ls)) {
#           node_ii <- seg_ls[[ii_seg]]
#           if (length(node_ii) > 1) {
#             seg_v <- rbind(seg_v, sub_tree_br[[ii_br]]$d[tail(node_ii,1), c("X","Y","Z")] - sub_tree_br[[ii_br]]$d[head(node_ii,1), c("X","Y","Z")])
#             so_L[ii_so, ii_br] <- so_L[ii_so, ii_br] + igraph::distances(sub_tree_br_g[[ii_br]], v = node_ii[1], to = tail(node_ii, 1) )
#           }
#         }
#       } 
#       so_N[ii_so, ii_br] <- length(seg_ls)
#     }
#     so_v_mat[ii_so, ] <- as.numeric(colSums(seg_v)) / sqrt(sum(colSums(seg_v)^2))
#     so_v_mat_raw[ii_so, ] <- as.numeric(colSums(seg_v))
#   }
# }
# 
# 
# # -- (11) -- arc length on unit sphere  ---------------------------------------------------------------------------
# # p1, p2 = [theta, phi] in radian, on unit sphere
# arcLength <- function(p1, p2) {
#   if (length(p1) != 2) {
#     stop("wrong input")
#   }
#   p1_xyz <- c(sin(p1[1])*cos(p1[2]), sin(p1[1])*sin(p1[2]), cos(p1[1]))
#   p2_xyz <- c(sin(p2[1])*cos(p2[2]), sin(p2[1])*sin(p2[2]), cos(p2[1]))
#   c2 <- sum((p1_xyz - p2_xyz)^2)
#   arc_ang <- acos((2-c2)/2)
#   return(arc_ang*1)
# }
# 
# 
# # -- (8) -- correct  CT1 kind -------------------------------------------------------------------------------------
# 
# # correct kink, break into 3 segments
# 
# kinkCT1 <- function(neu){
#   # p1 <- tar$d[tar$d$PointNo==15957800,]  # 120865 15957800     0 285640 305788 243000 -2 15957759
#   # p2 <- tar$d[tar$d$PointNo==17033465,] # 193597 17033465     0 287646 304732 244040 -2 17033464
#   # p3 <- tar$d[tar$d$PointNo==17081035,] # 195007 17081035     0 287767 304207 250320 -2 17081034
#   # p4 <- tar$d[tar$d$PointNo==17081048,] # 195081 17081048     0 283885 306765 251360 -2 17081047
#   
#   # load("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/R_common/p1p2p3p4_20200304.RData")
#   load("data/p1p2p3p4_20200304.RData")
#   
#   s1 <- xyzmatrix(p1)-xyzmatrix(p2) #note p1-p2 here
#   s2 <- xyzmatrix(p3)-xyzmatrix(p2)
#   s3 <- xyzmatrix(p4)-xyzmatrix(p3)
#   dz=40
#   k2x <- (s3[1]-s1[1]) / ((p3$Z-p2$Z)/dz)
#   k2y <- (s3[2]-s1[2]) / ((p3$Z-p2$Z)/dz)
#   k1x <- -s1[1]  / ((p2$Z-p1$Z)/dz)
#   k1y <- -s1[2]  / ((p2$Z-p1$Z)/dz)
#   k3x <- s3[1]  / ((p4$Z-p3$Z)/dz)
#   k3y <- s3[2]  / ((p4$Z-p3$Z)/dz)
#   
#   if (is.neuron(neu)) {
#     tar <- neu
#     
#     modd <- tar$d
#     modd2 <- tar$d
#     
#     for (zz in seq(p2$Z,p3$Z,dz)) {
#       ii <- modd[,5]==zz & modd[,3] < 320000
#       modd2[ii,3] <- modd[ii,3] + s1[1] + k2x*(zz-p2$Z)/dz
#       modd2[ii,4] <- modd[ii,4] + s1[2] + k2y*(zz-p2$Z)/dz
#     }
#     for (zz in seq(p1$Z,p2$Z-dz,dz)) {
#       ii <- modd[,5]==zz & modd[,3] < 320000
#       modd2[ii,3] <- modd[ii,3] - k1x*(zz-p1$Z)/dz
#       modd2[ii,4] <- modd[ii,4] - k1y*(zz-p1$Z)/dz
#     }
#     for (zz in seq(p3$Z+dz,p4$Z,dz)) {
#       ii <- modd[,5]==zz & modd[,3] < 320000
#       modd2[ii,3] <- modd[ii,3] + s3[1] - k3x*(zz-p3$Z)/dz
#       modd2[ii,4] <- modd[ii,4] + s3[2] - k3y*(zz-p3$Z)/dz
#     }
#     tar$d <- modd2
#     
#     return(tar)
#   } else if (is.neuronlist(neu)) {
#     tar_ls <- neu
#     for (j in 1:length(tar_ls)) {
#       tar <- tar_ls[[j]]
#       modd <- tar$d
#       modd2 <- tar$d
#       for (zz in seq(p2$Z,p3$Z,dz)) {
#         ii <- modd[,5]==zz & modd[,3] < 320000
#         modd2[ii,3] <- modd[ii,3] + s1[1] + k2x*(zz-p2$Z)/dz
#         modd2[ii,4] <- modd[ii,4] + s1[2] + k2y*(zz-p2$Z)/dz
#       }
#       for (zz in seq(p1$Z,p2$Z-dz,dz)) {
#         ii <- modd[,5]==zz & modd[,3] < 320000
#         modd2[ii,3] <- modd[ii,3] - k1x*(zz-p1$Z)/dz
#         modd2[ii,4] <- modd[ii,4] - k1y*(zz-p1$Z)/dz
#       }
#       for (zz in seq(p3$Z+dz,p4$Z,dz)) {
#         ii <- modd[,5]==zz & modd[,3] < 320000
#         modd2[ii,3] <- modd[ii,3] + s3[1] - k3x*(zz-p3$Z)/dz
#         modd2[ii,4] <- modd[ii,4] + s3[2] - k3y*(zz-p3$Z)/dz
#       }
#       tar$d <- modd2
#       tar_ls[[j]] <- tar
#     }
#     return(tar_ls)
#   } else {
#     print("invalid input")
#   }
#   
# }
# 
# # -- (9) -- return a quaternion rotation matrix -----------------------------------------------------------------------
# 
# # cf. https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
# 
# # vr is the axis of rotation [x,y,z], ang is angle in deg
# # position after rotation:  pt1 <- pt0 %*% t(R_mat), where R_mat <- quaternion3D(...)
# quaternion3D <- function(vr, ang){
#   ang <- ang / 180 * pi
#   vr <- vr / sqrt(sum(vr^2))
#   
#   qr <- cos(ang/2)
#   qi <- vr[1]*sin(ang/2)
#   qj <- vr[2]*sin(ang/2)
#   qk <- vr[3]*sin(ang/2)
#   
#   R <- matrix(c(
#     1-2*(qj^2+qk^2), 2*(qi*qj-qk*qr), 2*(qi*qk+qj*qr),
#     2*(qi*qj+qk*qr), 1-2*(qi^2+qk^2), 2*(qj*qk-qi*qr),
#     2*(qi*qk-qj*qr), 2*(qj*qk+qi*qr), 1-2*(qi^2+qj^2)),
#     ncol = 3, byrow = T)
#   
#   return(R)
# }
# 
# 
# 
# # -- (16) -- kernel regression ---------------------------------
# 
# # - 3d np, written for med_xyz -> ucl_rot_sm 
# # input, iv indep var (med), dv dep var (ucl), sp new sample, Nx3 [x,y,z]
# eyemap_np <- function(iv, dv, sp){
#   sp <- as.data.frame(sp)
#   colnames(sp) <- c('mc.x','mc.y','mc.z')
#   iv <- as.matrix(iv)
#   colnames(iv) <- c('x','y','z')
#   dv <- as.matrix(dv)
#   colnames(dv) <- c('x','y','z')
#   
#   kr <- matrix(nrow = nrow(sp), ncol = ncol(sp)) # kernel reg result
#   
#   for (k in 1:3) { #single variable np
#     npdata <- data.frame(mc = iv, ec = dv[,k])
#     # bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'generalized_nn', regtype= 'll')
#     # bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'generalized_nn', regtype= 'lc') 
#     # bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'fixed', regtype= 'll')
#     bw <- npregbw(formula= ec~mc.x+mc.y+mc.z, data= npdata, bwtype= 'adaptive_nn', regtype= 'll')
#     model_np <- npreg(bw)
#     for (j in 1:nrow(sp)) {
#       kr[j,k] <- predict(model_np, newdata = sp[j,])
#     }
#   }
#   colnames(kr) <- c('x','y','z')
#   
#   return(kr)
# }
# 
# # from Emil, read google sheet --------------------------------------------
# 
# # the dataframe that read_sheet("1hU_Wbb-uPclLBuE5s-BKmYN5snu7Io4jYKcfMdXgO2I", na = "NA",sheet = "left_lobe_v783")  returnsthe dataframe that read_sheet("1hU_Wbb-uPclLBuE5s-BKmYN5snu7Io4jYKcfMdXgO2I", na = "NA",sheet = "left_lobe_v783")  
# # returns x and y in the csv files are not directly indices though (y has negative values), more like Cartesian coordinates.
# 
# col_hex_grid <- function(google_sheet){
#   matrix = data.matrix(google_sheet)
#   
#   df = data.frame(col = as.numeric(),
#                   x = as.numeric(),
#                   y = as.numeric(),
#                   v_pos = as.numeric(),
#                   v_neg = as.numeric(),
#                   p_pos = as.numeric(),
#                   p_neg = as.numeric(),
#                   q_pos = as.numeric(),
#                   q_neg = as.numeric(),
#                   h_pos = as.numeric(),
#                   h_neg = as.numeric())
#   
#   cc = 0
#   
#   for (i in 3:(nrow(matrix)-2)) {
#     for (j in 3:(ncol(matrix)-2)) {
#       df[cc+j-1, "col"] = matrix[i,j]
#       df[cc+j-1, "x"] = j
#       df[cc+j-1, "y"] = -i
#       df[cc+j-1, "v_pos"] = matrix[i-2,j]
#       df[cc+j-1, "v_neg"] = matrix[i+2,j]
#       df[cc+j-1, "p_pos"] = matrix[i-1,j-1]
#       df[cc+j-1, "p_neg"] = matrix[i+1,j+1]
#       df[cc+j-1, "q_pos"] = matrix[i-1,j+1]
#       df[cc+j-1, "q_neg"] = matrix[i+1,j-1]
#       df[cc+j-1, "h_pos"] = matrix[i,j+2]
#       df[cc+j-1, "h_neg"] = matrix[i,j-2]
#     }
#     cc =  cc + length(2:(ncol(matrix)-1)) 
#   }
#   
#   df = df %>% filter(!is.na(col)) 
#   return(df)
# }
# 
# 
# # spatial profile of T4 dendrites -----------------------------------------
# 
# 
# PDOD_T4 <- function(LL){
#   com_xyz <- dir_type[[LL]][, c('comx','comy','comz')] %>% as.matrix()
#   v0 <- dir_type[[LL]][ ,c('rsx0','rsy0','rsz0')] %>% as.matrix()
#   v1 <- dir_type[[LL]][ ,c('rsxd','rsyd','rszd')] %>% as.matrix()
#   v3 <- dir_type[[LL]][ ,c('odx0','ody0','odz0')] %>% as.matrix()
#   v4 <- dir_type[[LL]][ ,c('odxd','odyd','odzd')] %>% as.matrix()
#   ellabc <- ell_type[[LL]][, c('ea','eb','ec')] %>% as.matrix()
#   
#   # nopen3d()
#   # points3d(med_xyz)
#   # for (j in 1:14) {
#   #   ii <- match(vaxis_gen(j), eyemap[,2])
#   #   points3d(med_xyz[ii,], col=j, size=10)
#   #   ii <- match(vaxis_gen(-j), eyemap[,2])
#   #   points3d(med_xyz[ii,], col=j, size=10)
#   # }
#   # for (j in 1:nrow(v0)) {
#   #   arrow3d(v0[j,], v1[j,], theta = pi/6, n = 8, col =pal_T4[LL], type = "rotation", lit=F)
#   # }
#   
#   # - size, normalized to 2-up/down, 4-sides, 6-all, 22-row dist
#   ellabc_uL6 <- matrix(ncol = 3, nrow = nrow(ellabc))
#   PD_uL6 <- matrix(ncol = 1, nrow = nrow(v0))
#   OD_uL6 <- matrix(ncol = 1, nrow = nrow(v0))
#   uL6 <- matrix(ncol = 1, nrow = nrow(v0))
#   for (j in 1:nrow(v0)) {
#     ii <- sweep(med_xyz, 2, com_xyz[j,])^2 %>% rowSums() %>% which.min()
#     ii7 <- nb_ind[ii,]
#     uL <- sweep(med_xyz[na.omit(ii7[-1]),,drop=F], 2, med_xyz[ii7[1],])^2 %>% rowSums() %>% sqrt() %>% mean()
#     uL6[j] <- uL
#     ellabc_uL6[j,] <- ellabc[j, ] / uL
#     PD_uL6[j] <- sqrt(sum((v1[j,]-v0[j,])^2)) / uL
#     OD_uL6[j] <- sqrt(sum((v4[j,]-v3[j,])^2)) / uL
#   }
#   
#   # -- size, normalized to side edges
#   ellabc_uL4 <- matrix(ncol = 3, nrow = nrow(ellabc))
#   PD_uL4 <- matrix(ncol = 1, nrow = nrow(v0))
#   OD_uL4 <- matrix(ncol = 1, nrow = nrow(v0))
#   uL4 <- matrix(ncol = 1, nrow = nrow(v0))
#   for (j in 1:nrow(v0)) {
#     ii <- sweep(med_xyz, 2, com_xyz[j,])^2 %>% rowSums() %>% which.min()
#     ii7 <- nb_ind[ii,]
#     # uL <- sweep(med_xyz[na.omit(ii7[2:5]),,drop=F], 2, med_xyz[ii7[1],])^2 %>% rowSums() %>% sqrt() %>% mean() #2021
#     uL <- sweep(med_xyz[na.omit(ii7[c(2,4,5,7)]),,drop=F], 2, med_xyz[ii7[1],])^2 %>% rowSums() %>% sqrt() %>% mean() #2023
#     uL4[j] <- uL
#     ellabc_uL4[j,] <- ellabc[j, ] / uL
#     PD_uL4[j,] <- sqrt(sum((v1[j,]-v0[j,])^2)) / uL
#     OD_uL4[j,] <- sqrt(sum((v4[j,]-v3[j,])^2)) / uL
#   }
#   
#   # up/down
#   ellabc_uL2 <- matrix(ncol = 3, nrow = nrow(ellabc))
#   PD_uL2 <- matrix(ncol = 1, nrow = nrow(v0))
#   OD_uL2 <- matrix(ncol = 1, nrow = nrow(v0))
#   uL2 <- matrix(ncol = 1, nrow = nrow(v0))
#   for (j in 1:nrow(v0)) {
#     ii <- sweep(med_xyz, 2, com_xyz[j,])^2 %>% rowSums() %>% which.min()
#     ii7 <- nb_ind[ii,]
#     # uL <- sweep(med_xyz[na.omit(ii7[6:7]),,drop=F], 2, med_xyz[ii7[1],])^2 %>% rowSums() %>% sqrt() %>% mean() # 2021
#     uL <- sweep(med_xyz[na.omit(ii7[c(3,6)]),,drop=F], 2, med_xyz[ii7[1],])^2 %>% rowSums() %>% sqrt() %>% mean() # 2023
#     uL2[j] <- uL
#     ellabc_uL2[j,] <- ellabc[j, ] / uL
#     # PD_uL2[j,] <- sqrt(sum((v1[j,]-v0[j,])^2)) / uL
#     # OD_uL2[j,] <- sqrt(sum((v4[j,]-v3[j,])^2)) / uL
#     PD_uL2[j,] <- sqrt(sum((v1[j,]-v0[j,])^2)) / ((2*uL2[j]+uL4[j]))
#     OD_uL2[j,] <- sqrt(sum((v4[j,]-v3[j,])^2)) / ((2*uL2[j]+uL4[j]))
#   }
#   
#   # -- size, normalized to side-side dist
#   ellabc_uL22 <- matrix(ncol = 3, nrow = nrow(ellabc))
#   PD_uL22 <- matrix(ncol = 1, nrow = nrow(v0))
#   OD_uL22 <- matrix(ncol = 1, nrow = nrow(v0))
#   uL22 <- matrix(ncol = 1, nrow = nrow(v0))
#   for (j in 1:nrow(v0)) {
#     ii <- sweep(med_xyz, 2, com_xyz[j,])^2 %>% rowSums() %>% which.min()
#     ii7 <- nb_ind[ii,]
#     if (!is.na(sum(med_xyz[ii7[2:5],]))) {
#       # uL <- sqrt(rowSums((med_xyz[ii7[c(2,3)],]- med_xyz[ii7[c(5,4)],])^2)) %>%mean()#2021
#       uL <- sqrt(rowSums((med_xyz[ii7[c(2,7)],]- med_xyz[ii7[c(4,5)],])^2)) %>%mean()#2023
#       uL22[j] <- uL
#       ellabc_uL22[j,] <- ellabc[j, ] / uL
#       # PD_uL22[j,] <- sqrt(sum((v1[j,]-v0[j,])^2)) / uL
#       # OD_uL22[j,] <- sqrt(sum((v4[j,]-v3[j,])^2)) / uL
#       PD_uL22[j,] <- sqrt(sum((v1[j,]-v0[j,])^2)) / (uL+uL4[j])
#       OD_uL22[j,] <- sqrt(sum((v4[j,]-v3[j,])^2)) / (uL+uL4[j])
#     }
#   }
#   
#   ###
#   if (LL == 2) {
#     PDOD_T4 <- cbind(PD_uL2, OD_uL2,
#                       PD_uL4, OD_uL4,
#                       PD_uL6, OD_uL6,
#                       PD_uL22, OD_uL22 )
#   }
#   if (LL ==4) {
#     PDOD_T4 <- cbind(PD_uL2, OD_uL2,
#                       PD_uL4, OD_uL4,
#                       PD_uL6, OD_uL6,
#                       PD_uL22, OD_uL22 )
#   }
#   
#   return(PDOD_T4)
# }
