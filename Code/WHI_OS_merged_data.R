# Load packages ################################################################
library(dplyr)
library(fastDummies)
library(match2C)
library(gtsummary)
library(fastDummies)

# Set base path for data files.
# Not setting WD to avoid accidentally outputting
# files to the wrong place.
path_in <- 'C:/Users/chenk/Dropbox/kan_dylan/attributable_effect_case_study/datasets/data_needed'


# Path to output processed data set
path_out <- 'C:/Users/chenk/Dropbox/kan_dylan/attributable_effect_case_study/datasets/data_needed'

################################################################################
### Import data
################################################################################



dem_os_pub <- read.csv(paste0(path_in, "/dem_os_pub.csv"))
outc_adj_os_pub <- read.csv(paste0(path_in, "/outc_adj_os_pub.csv"))
f34_os_pub <- read.csv(paste0(path_in, "/f34_os_pub.csv"))
f43_os_pub <- read.csv(paste0(path_in, "/f43_os_pub.csv"))
f80_os_pub <- read.csv(paste0(path_in, "/f80_os_pub.csv"))
f80_os_pub <- f80_os_pub[f80_os_pub$F80DAYS < 0, ]
outc_bc_os_pub <- read.csv(paste0(path_in, "/outc_bc_os_pub.csv"))


################################################################################
### Merge data
################################################################################


d1 <- dem_os_pub[c('ID','AGE','ETHNIC','EDUC','INCOME')]
d2 <- outc_adj_os_pub[c('ID','BREAST','BREASTINV')]
d3 <- f34_os_pub[c('ID','ALCNOW','SMOKNOW')]
d4 <- f43_os_pub[c('ID','TOTH')]
d5 <- f80_os_pub[c('ID','BMIX','SYSTBP1','DIASBP1')]
d5 <- d5[complete.cases(d5),]
d6 <- outc_bc_os_pub[c('ID','STAGE','BEHAVIOR','ERASSAY','PRASSAY')]

data <- merge(d1,d2,by = 'ID')
data <- merge(data,d3,by = 'ID')
data <- merge(data,d4,by = 'ID')
data <- merge(data,d5,by = 'ID')


merged_data <- merge(data,d6,by = 'ID', all = TRUE)
merged_data$SENSITIVE <- NA
merged_data$SENSITIVE <- (merged_data$ERASSAY == 1 | merged_data$PRASSAY == 1)
merged_data <- merged_data[!is.na(merged_data$ALCNOW),] 

################################################################################
### Match data
### Matching based on the outcome (invasive breast cancer)
################################################################################




data = dummy_cols(merged_data, select_columns = c("ETHNIC","EDUC","INCOME"
                                 ,"SMOKNOW"),remove_selected_columns = TRUE)


data[is.na(data$ETHNIC_1),]$ETHNIC_1 = 0
data[is.na(data$ETHNIC_2),]$ETHNIC_2 = 0
data[is.na(data$ETHNIC_3),]$ETHNIC_3 = 0
data[is.na(data$ETHNIC_4),]$ETHNIC_4 = 0
data[is.na(data$ETHNIC_5),]$ETHNIC_5 = 0
data[is.na(data$ETHNIC_8),]$ETHNIC_8 = 0

data[is.na(data$EDUC_1),]$EDUC_1 = 0
data[is.na(data$EDUC_2),]$EDUC_2 = 0
data[is.na(data$EDUC_3),]$EDUC_3 = 0
data[is.na(data$EDUC_4),]$EDUC_4 = 0
data[is.na(data$EDUC_5),]$EDUC_5 = 0
data[is.na(data$EDUC_6),]$EDUC_6 = 0
data[is.na(data$EDUC_7),]$EDUC_7 = 0
data[is.na(data$EDUC_8),]$EDUC_8 = 0
data[is.na(data$EDUC_9),]$EDUC_9 = 0
data[is.na(data$EDUC_10),]$EDUC_10 = 0
data[is.na(data$EDUC_11),]$EDUC_11 = 0

data[is.na(data$INCOME_1),]$INCOME_1 = 0
data[is.na(data$INCOME_2),]$INCOME_2 = 0
data[is.na(data$INCOME_3),]$INCOME_3 = 0
data[is.na(data$INCOME_4),]$INCOME_4 = 0
data[is.na(data$INCOME_5),]$INCOME_5 = 0
data[is.na(data$INCOME_6),]$INCOME_6 = 0
data[is.na(data$INCOME_7),]$INCOME_7 = 0
data[is.na(data$INCOME_8),]$INCOME_8 = 0
data[is.na(data$INCOME_9),]$INCOME_9 = 0


data[is.na(data$SMOKNOW_0),]$SMOKNOW_0 = 0
data[is.na(data$SMOKNOW_1),]$SMOKNOW_1 = 0


Z = data$BREASTINV
X = as.matrix(data[c(2,5,6,7,8,9,15:20,22:32,34:42,44,45)])
dataset_in_dataframe = data
# Estimate pscore
propens = glm(Z ~ X, family = 'binomial')$fitted.value

dist_list_pscore_maha_soft = create_list_from_scratch(Z, X, exact = NULL, 
                                                      p = propens, 
                                                      caliper_low = 10, 
                                                      k = 10000,
                                                      method = 'robust maha', 
                                                      penalty = Inf) 

matching_output = match_2C_list(Z, dataset_in_dataframe, 
                                dist_list_pscore_maha_soft, 
                                dist_list_2 = NULL, 
                                controls = 1)


matched_data = matching_output$matched_data_in_order
matched_data = matched_data[!is.na(matched_data$matched_set),]

matched_data$SENSITIVE = (matched_data$ERASSAY == 1 | matched_data$PRASSAY == 1)

index_set_sensitive = unique(na.omit(matched_data[matched_data$SENSITIVE 
                                                  == TRUE,]$matched_set))
index_set_insensitive = unique(na.omit(matched_data[matched_data$SENSITIVE 
                                                    == FALSE,]$matched_set))

data_sensitive = matched_data[matched_data$matched_set 
                              == index_set_sensitive[1],]

k = 1
for (i in index_set_sensitive[2:length(index_set_sensitive)]) {
  d = matched_data[matched_data$matched_set == i,]
  data_sensitive = rbind(data_sensitive,  d[order(-d$BREASTINV),])
  k = k + 1
  print(k)
}


data_insensitive = matched_data[matched_data$matched_set == 
                                  index_set_insensitive[1],]
k = 1
for (i in index_set_insensitive[2:length(index_set_insensitive)]) {
  d = matched_data[matched_data$matched_set == i,]
  data_insensitive = rbind(data_insensitive, d[order(-d$BREASTINV),])
  k = k + 1
  print(k)
}




matched_data1 = data_sensitive
matched_data2 = data_insensitive
final_data_used = rbind(data_sensitive,data_insensitive)

check_balance(Z, matching_output, 
              cov_list = colnames(final_data_used)[c(2,6,7,8,9,15:20,
                                                     22:32,34:42,44,45)],
              plot_propens = FALSE)

matched_data1 = matched_data1 %>%
  left_join(f34_os_pub, by = c('ID','ALCNOW')) 


matched_data2 = matched_data2 %>%
  left_join(f34_os_pub, by = c('ID','ALCNOW')) 



matched_data1$ALCO7MORE = (matched_data1$ALCSWK >= 18)
matched_data2$ALCO7MORE = (matched_data2$ALCSWK >= 18)

matched_data1[is.na(matched_data1$ALCO7MORE),]$ALCO7MORE = FALSE
matched_data2[is.na(matched_data2$ALCO7MORE),]$ALCO7MORE = FALSE

matched_data1$Z = matched_data1$ALCO7MORE
matched_data1$R = matched_data1$BREASTINV

matched_data2$Z = matched_data2$ALCO7MORE
matched_data2$R = matched_data2$BREASTINV





# Data characteristic

data = data %>%
  left_join(f34_os_pub, by = c('ID','ALCNOW')) 

data$ALCO7MORE = (data$ALCSWK >= 18)

data[is.na(data$ALCO7MORE),]$ALCO7MORE = FALSE


final_data_used = final_data_used %>%
  left_join(f34_os_pub, by = c('ID','ALCNOW')) 

final_data_used$ALCO7MORE = (final_data_used$ALCSWK >= 18)

final_data_used[is.na(final_data_used$ALCO7MORE),]$ALCO7MORE = FALSE



data2 <- data %>% select(AGE,BREAST,BREASTINV,SENSITIVE,ALCNOW,ALCO7MORE,TOTH,
                         BMIX,SYSTBP1,DIASBP1,EDUC_1,
                         EDUC_2,EDUC_3,EDUC_4,EDUC_5,
                         EDUC_6,EDUC_7,EDUC_8,EDUC_9,EDUC_10,EDUC_11,EDUC_NA,
                         ETHNIC_1,ETHNIC_2,ETHNIC_3,ETHNIC_4,ETHNIC_5,ETHNIC_8,
                         ETHNIC_NA,
                         INCOME_1,INCOME_2,INCOME_3,INCOME_4,INCOME_5,INCOME_6,
                         INCOME_7,INCOME_8,INCOME_9,INCOME_NA,
                         SMOKNOW_0,SMOKNOW_1,SMOKNOW_NA)


data3 <- final_data_used %>% select(AGE,BREAST,BREASTINV,SENSITIVE,ALCNOW,
                                    ALCO7MORE,TOTH,BMIX,SYSTBP1,DIASBP1,EDUC_1,
                                    EDUC_2,EDUC_3,EDUC_4,EDUC_5,
                                    EDUC_6,EDUC_7,EDUC_8,EDUC_9,EDUC_10,
                                    EDUC_11,EDUC_NA,
                                    ETHNIC_1,ETHNIC_2,ETHNIC_3,ETHNIC_4,
                                    ETHNIC_5,ETHNIC_8,ETHNIC_NA,
                                    INCOME_1,INCOME_2,INCOME_3,INCOME_4,INCOME_5
                                    ,INCOME_6,
                                    INCOME_7,INCOME_8,INCOME_9,INCOME_NA,
                                    SMOKNOW_0,SMOKNOW_1,SMOKNOW_NA)


data2 %>% tbl_summary(by = BREASTINV)


data3 %>% tbl_summary(by = BREASTINV)


write.csv(matched_data1,'matched_data1.csv')
write.csv(matched_data2,'matched_data2.csv')


