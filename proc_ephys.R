# process H2 recordings

library(tidyverse)
library(readxl)

source("eyemap_func.R")

# process 2021 grating data  ---------------------------------------------------------------------

tb <- read_excel("data_Eyal/H2_2021/movGrtTabwEZ_ThetaCorr.xlsx", sheet = 1, col_names = T ) %>% data.frame()
tb_v <- Eyal_arena_2023(tb)

save(tb_v, tb, file = 'data/H2_tuning.rda')

# process new data 2023 edge -----------------------------------------------------------

tb <- read_excel("data_Eyal/H2_2023/movEdgeTabwPY_ThetaCorr_fullResp.xlsx", sheet = 1, col_names = T ) %>% data.frame()
names(tb)[names(tb) == "headAngPitch"] <- "headAng"

tb_v <- Eyal_arena_2023(tb)

# tb_v [theta(base), phi(base), spk, spk, sub, sub, dark=0/bright=1]
save(tb_v, tb, file = 'data/H2_tuning_2023.rda')

# process 2023 grating data ----------------------------------------------------------------

tb <- read_excel("data_Eyal/H2_2023/movGrtTabwPY_ThetaCorr.xlsx", sheet = 1, col_names = T ) %>% data.frame()
names(tb)[names(tb) == "headAngPitch"] <- "headAng"

tb_v <- Eyal_arena_2023(tb)

save(tb_v, tb, file = 'data/H2_tuning_2023_grating.rda')



