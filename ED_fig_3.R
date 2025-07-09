# see Fig_2_EM.R for ED Fig.3B, ED Fig.3C, ED Fig.3D, 

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
## ## choose type
LL <- 1

ssn <- dplyr::bind_rows(seg_summ_type[[LL]])
df <- matrix(ncol = 3, nrow = 0)
for (ii_h in 1:4) {
  if (nrow(ssn[ssn$so == ii_h,]) > 0 ) {
    if (LL <= 2) {
      hh <- hist(ssn[ssn$so == ii_h,]$angH, breaks = seq(0,180,dtheta), plot = F)
      df <- rbind(df, cbind(hh$mids, hh$density/h_norm, ii_h))
    } else {
      hh <- hist(ssn[ssn$so == ii_h,]$angV, breaks = seq(0,180,dtheta), plot = F)
      df <- rbind(df, cbind(hh$mids, hh$density/h_norm, ii_h))
    }
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


