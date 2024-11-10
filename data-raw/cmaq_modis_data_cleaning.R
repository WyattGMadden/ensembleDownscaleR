
###STAGE 1###

#cmaq 
rel_loc <- "~/Documents/school/emory/research/chang/nasa/code/onedrive_code/"
Y.input.cmaq <- read.csv(paste0(rel_loc, "Stage 1 Input Data/CTM/Y.csv"))
X.input.cmaq <- read.csv(paste0(rel_loc, "Stage 1 Input Data/CTM/X.csv"))
L.input.cmaq <- as.matrix(read.csv(paste0(rel_loc, "Stage 1 Input Data/CTM/L.csv")))
M.input.cmaq <- as.matrix(read.csv(paste0(rel_loc, "Stage 1 Input Data/CTM/M.csv")))
monitor.locs.cmaq <- as.matrix(read.csv(paste0(rel_loc, "Stage 1 Input Data/CTM/Monitor_XY.csv")))
cmaq.dateinfo <- read.csv(paste0(rel_loc, "Stage 3 Input Data/CTM_Date_Mon_ID.csv"))

cmaq_aqs_matched <- Y.input.cmaq
cmaq_aqs_matched$ctm <- X.input.cmaq$CTM

L_cmaq_names <- colnames(L.input.cmaq)[4:ncol(L.input.cmaq)]
cmaq_aqs_matched[, L_cmaq_names] <- L.input.cmaq[, L_cmaq_names]

M_cmaq_names <- colnames(M.input.cmaq)[4:ncol(M.input.cmaq)]
cmaq_aqs_matched[, M_cmaq_names] <- M.input.cmaq[, M_cmaq_names]

cmaq_aqs_matched$Date <- as.Date(cmaq.dateinfo$Date)

monitor.locs.cmaq <- as.data.frame(monitor.locs.cmaq)

cmaq_aqs_matched <- merge(cmaq_aqs_matched, 
                 monitor.locs.cmaq, 
                 by.x = "Space_ID", 
                 by.y = "ID",
                 all.x = T,
                 sort = F)

names(cmaq_aqs_matched) <- tolower(names(cmaq_aqs_matched))
names(cmaq_aqs_matched) <- gsub("\\.", "_", names(cmaq_aqs_matched))

cmaq_aqs_matched <- cmaq_aqs_matched[, c(c("time_id", "space_id", "spacetime_id"), 
                         names(cmaq_aqs_matched)[4:ncol(cmaq_aqs_matched)])]
cmaq_aqs_matched <- cmaq_aqs_matched[with(cmaq_aqs_matched, order(time_id, space_id, spacetime_id)), ]

#filter for only space.id's that are in every spacetime
n.space <- max(cmaq_aqs_matched$space_id)
n.spacetime <- max(cmaq_aqs_matched$spacetime_id)
all_ids <- paste0(rep(1:n.spacetime, each = n.space), "_", rep(1:n.space, times = n.spacetime))
all_ids <- expand.grid(space_id = 1:n.space, spacetime_id = 1:n.spacetime)
obs_ids <- unique(cmaq_aqs_matched[, c("space_id", "spacetime_id")])
obs_ids$dummy <- 1
all_ids <- merge(all_ids, obs_ids, by = c("space_id", "spacetime_id"), all.x = T)
to_remove <- all_ids[is.na(all_ids$dummy), c("space_id", "spacetime_id")]
cmaq_aqs_matched <- cmaq_aqs_matched[!(cmaq_aqs_matched$space_id %in% to_remove$space_id), ]
#reset space_id
cmaq_aqs_matched$space_id <- as.numeric(as.factor(cmaq_aqs_matched$space_id))



usethis::use_data(cmaq_aqs_matched, overwrite = TRUE)


#modis
Y.input.modis <- read.csv(paste0(rel_loc, "Stage 1 Input Data/MAIA/Y.csv"))
X.input.modis <- read.csv(paste0(rel_loc, "Stage 1 Input Data/MAIA/X.csv"))
L.input.modis <- as.matrix(read.csv(paste0(rel_loc, "Stage 1 Input Data/MAIA/L.csv")))
M.input.modis <- as.matrix(read.csv(paste0(rel_loc, "Stage 1 Input Data/MAIA/M.csv")))
monitor.locs.modis <- as.matrix(read.csv(paste0(rel_loc, "Stage 1 Input Data/MAIA/Monitor_XY.csv")))
modis.dateinfo <- read.csv(paste0(rel_loc, "Stage 3 Input Data/MAIA_Date_Mon_ID.csv"))

modis_aqs_matched <- Y.input.modis
modis_aqs_matched$aod <- X.input.modis$aod

L_modis_names <- colnames(L.input.modis)[4:ncol(L.input.modis)]
modis_aqs_matched[, L_modis_names] <- L.input.modis[, L_modis_names]

M_modis_names <- colnames(M.input.modis)[4:ncol(M.input.modis)]
modis_aqs_matched[, M_modis_names] <- M.input.modis[, M_modis_names]

modis_aqs_matched$Date <- as.Date(modis.dateinfo$Date)

monitor.locs.modis <- as.data.frame(monitor.locs.modis)

modis_aqs_matched <- merge(modis_aqs_matched, 
                 monitor.locs.cmaq, 
                 by.x = "Space_ID", 
                 by.y = "ID",
                 all.x = T,
                 sort = F)


names(modis_aqs_matched) <- tolower(names(modis_aqs_matched))
names(modis_aqs_matched) <- gsub("\\.", "_", names(modis_aqs_matched))

modis_aqs_matched <- modis_aqs_matched[, c(c("time_id", "space_id", "spacetime_id"), 
                         names(modis_aqs_matched)[4:ncol(modis_aqs_matched)])]
modis_aqs_matched <- modis_aqs_matched[with(modis_aqs_matched, order(time_id, space_id, spacetime_id)), ]

#filter for only space.id's that are in every spacetime
n.space <- max(modis_aqs_matched$space_id)
n.spacetime <- max(modis_aqs_matched$spacetime_id)
all_ids <- paste0(rep(1:n.spacetime, each = n.space), "_", rep(1:n.space, times = n.spacetime))
all_ids <- expand.grid(space_id = 1:n.space, spacetime_id = 1:n.spacetime)
obs_ids <- unique(modis_aqs_matched[, c("space_id", "spacetime_id")])
obs_ids$dummy <- 1
all_ids <- merge(all_ids, obs_ids, by = c("space_id", "spacetime_id"), all.x = T)
to_remove <- all_ids[is.na(all_ids$dummy), c("space_id", "spacetime_id")]
modis_aqs_matched <- modis_aqs_matched[!(modis_aqs_matched$space_id %in% to_remove$space_id), ]
#reset space_id
modis_aqs_matched$space_id <- as.numeric(as.factor(modis_aqs_matched$space_id))
sort(unique(modis_aqs_matched$space_id))

usethis::use_data(modis_aqs_matched, overwrite = TRUE)


###STAGE 2###

#cmaq
pred.locs.cmaq <- as.matrix(read.csv(paste0(rel_loc, "Stage 2 Input Data/CTM/Cell_XY.csv")))
L.pred.cmaq <- as.matrix(read.csv(paste0(rel_loc, "Stage 2 Input Data/CTM/L.csv")))
M.pred.cmaq <- as.matrix(read.csv(paste0(rel_loc, "Stage 2 Input Data/CTM/M.csv")))
X.pred.cmaq <- read.csv(paste0(rel_loc, "Stage 2 Input Data/CTM/X.csv"))
cmaq.dateinfo.pred <- read.csv(paste0(rel_loc, "Stage 4 Input Data/CTM_Pred_Date_Cell_ID.csv"))

cmaq_full <- X.pred.cmaq

cmaq_full[, M_cmaq_names] <- M.pred.cmaq[, M_cmaq_names]
cmaq_full[, L_cmaq_names] <- L.pred.cmaq[, L_cmaq_names]

cmaq_full$Date <- as.Date(cmaq.dateinfo.pred$Date)

pred.locs.cmaq <- as.data.frame(pred.locs.cmaq)

cmaq_full <- merge(cmaq_full, 
                 pred.locs.cmaq, 
                 by.x = "Space_ID", 
                 by.y = "Cell",
                 all.x = T,
                 sort = F)


names(cmaq_full) <- tolower(names(cmaq_full))
names(cmaq_full) <- gsub("\\.", "_", names(cmaq_full))

cmaq_full <- cmaq_full[, c(c("time_id", "space_id", "spacetime_id"), 
                         names(cmaq_full)[4:ncol(cmaq_full)])]
cmaq_full <- cmaq_full[with(cmaq_full, order(time_id, space_id, spacetime_id)), ]

#filter predictions for only June 2004 to save package memory. 
#this month has a lot of missing aod and no missing cmaq
cmaq_full <- cmaq_full[cmaq_full$date < as.Date("2004-07-01") & 
                       cmaq_full$date >= as.Date("2004-06-01"), ]


usethis::use_data(cmaq_full, overwrite = TRUE)



#modis
pred.locs.modis = as.matrix(read.csv(paste0(rel_loc, "Stage 2 Input Data/MAIA/Cell_XY.csv")))
L.pred.modis = as.matrix(read.csv(paste0(rel_loc, "Stage 2 Input Data/MAIA/L.csv")))
M.pred.modis = as.matrix(read.csv(paste0(rel_loc, "Stage 2 Input Data/MAIA/M.csv")))
X.pred.modis = read.csv(paste0(rel_loc, "Stage 2 Input Data/MAIA/X.csv"))
modis.dateinfo.pred <- read.csv(paste0(rel_loc, "Stage 4 Input Data/MAIA_Pred_Date_Cell_ID.csv"))

modis_full <- X.pred.modis

modis_full[, M_modis_names] <- M.pred.modis[, M_modis_names]
modis_full[, L_modis_names] <- L.pred.modis[, L_modis_names]

modis_full$Date <- as.Date(modis.dateinfo.pred$Date)

pred.locs.modis <- as.data.frame(pred.locs.modis)

modis_full <- merge(modis_full, 
                 pred.locs.modis, 
                 by.x = "Space_ID", 
                 by.y = "Cell",
                 all.x = T,
                 sort = F)


names(modis_full) <- tolower(names(modis_full))
names(modis_full) <- gsub("\\.", "_", names(modis_full))

modis_full <- modis_full[, c(c("time_id", "space_id", "spacetime_id"), 
                         names(modis_full)[4:ncol(modis_full)])]
modis_full <- modis_full[with(modis_full, order(time_id, space_id, spacetime_id)), ]

#filter predictions for only June 2004 to save package memory. 
#this month has a lot of missing aod and no missing cmaq
modis_full <- modis_full[modis_full$date < as.Date("2004-07-01") & 
                         modis_full$date >= as.Date("2004-06-01"), ]

usethis::use_data(modis_full, overwrite = TRUE)



