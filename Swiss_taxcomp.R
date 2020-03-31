library(dplyr)
library(plm)
library(lfe)
library(stargazer)
library(stringr)
library(tidyr)
library(spdep)
library(splm)


#################################################################################
#       Strategic interaction in Swiss local income tax setting revisited       #
#################################################################################
#                          Fabian Muny and Lukas Peham                          #
#           Master Course in Spatial Economics - Winter Term 2019/20            #
#                  Vienna University of Economics and Business                  #
#################################################################################

###
# 1.) Generate samples
###

#load data
load("data.RData")
data.s <- data
data <- data.s[[1]]

#Create variables for municipality and canton-year fixed effects
data$mFE <- as.factor(data$gdenr)
data$cyFE <- as.factor(paste(data$canton, data$year, sep=""))

#Create datasets:
#Full sample
pdata <- pdata.frame(data, index = c("gdenr", "year"))
#Generate variable that gives proportion of neighbors in another canton
pdata$intens_10 <- (pdata$nb_z_10 / pdata$nb_10)*100
#Generate lag of the 10km-neighbor average cantonal tax rate of unmarried singles with income of CHF 1,000,000
pdata$W_10_z_i_s_1000_l <- lag(pdata$W_10_z_i_s_1000, k=1)
pdata$W_10_z_i_s_1000_l <- pdata$W_10_z_i_s_1000_l*(pdata$intens_10/100)
#Generate sample of border municipalities (delete observation with missing values for W_10_z_i_s_1000_l)
pdata_border <- pdata.frame(filter(pdata, !is.na(W_10_z_i_s_1000_l)), index = c("gdenr", "year"))

#balanced sample:
load("data.RData")
#p.df <- data[[1]] %>% filter(gdenr %in% shp$bfsnr)
p.df <- data[[1]]
#Formula just for model.matrix; equal to select all the variables of RHS of formula
fm <- canton ~ gdenr + year + p_i_s_1000  + W_10_p_i_s_1000 + population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + 
  employment_pc + pop_unemployment + left + movie + W_10_population + W_10_pop_foreign + W_10_pop_young + W_10_pop_old +
  W_10_pop_secondary + W_10_pop_tertiary + W_10_employment_pc + W_10_pop_unemployment + W_10_left + W_10_movie -1

p.df <- model.matrix(fm , p.df)
#complete.cases(p.df)
sum(!complete.cases(p.df)) #no NAs
p.df <- data.frame(p.df)

is.pbalanced(p.df) #our data is not balanced

#make p.df balanced
range_years <- p.df %>% group_by(gdenr) %>% count()
range(range_years$n)
dim(range_years[range_years$n!=30,])
dim(range_years[range_years$n==12,])
dim(range_years[range_years$n==27,])
dim(range_years[(range_years$n!=12 & range_years$n!=27 & range_years$n!=30),])
range_years[(range_years$n!=12 & range_years$n!=27 & range_years$n!=30),]
#for 395 munis the data is not over the full time period
#99: 12years, 289: 27y, 7: 28y

#delete the 12year munis 
#m12 <- range_years[range_years$n==12,]$gdenr
p.df <- p.df %>% filter(!(gdenr %in% range_years[range_years$n==12,]$gdenr))

munis_per_y <- p.df %>% group_by(year) %>% count()
#only 2015 of 2304 munis have data for 1983-1985; delete these years
p.df <- p.df %>% filter(!(year %in% 1983:1985))

#finally delete the 7 municipalities which have no data for 2011:2012
p.df <- p.df %>% filter(!(gdenr %in% range_years[range_years$n==28,]$gdenr))

#Create CANTON_YEAR dummy, which we need to include canton-year fixed effects
p.df <- left_join(p.df, data[[1]] %>% select(gdenr, year, canton))
p.df$cyFE <- as.factor(paste(p.df$canton, p.df$year, sep=""))
p.df$mFE <- as.factor(p.df$gdenr)

is.pbalanced(p.df) #CHECK :)

pdata_balanced <- pdata.frame(p.df, index = c("gdenr", "year"))


###
# 2.) Replication of weights matrix
###

#Balanced weights matrix
#a) Approach Fabian:
load("wdf.RData")
df_10 <- filter(w.df, dist_m<=10)
df_10$value <- case_when(df_10$dist_m<=0~0, df_10$dist_m>0~1)
df_10 <- select(df_10, gdenr_fr, gdenr_to, value)
df_10_wide <- spread(df_10, gdenr_to, value)
rownames(df_10_wide) <- df_10_wide$gdenr_fr
gdenr_fr <- df_10_wide$gdenr_fr
W_mat_10 <- df_10_wide[-1]
W_mat_10[is.na(W_mat_10)] <- 0
W_mat_10$rsum <- rowSums(W_mat_10)
W_mat_10$gdenr_fr <- gdenr_fr
W_mat_10 <- filter(W_mat_10, rsum > 0)
W_mat_10 <- filter(W_mat_10, gdenr_fr %in% unique(pdata_balanced$gdenr))
rownames(W_mat_10) <- W_mat_10$gdenr_fr
W_mat_10 <- select(W_mat_10, as.character(W_mat_10$gdenr_fr))
W_listw_10 <- mat2listw(as.matrix(W_mat_10), row.names = as.numeric(row.names(as.matrix(W_mat_10))), style="W")
attributes(W_listw_10)
nblist_10 <- W_listw_10[["neighbours"]]

#Check
length(unique(pdata_balanced$gdenr)) 
#-> islands that have been dropped from weights matrix are still in the balanced panel
length(nblist_10)

pdata_balanced <- pdata_balanced %>% filter(as.character(pdata_balanced$gdenr) %in% rownames(W_mat_10))
length(unique(pdata_balanced$gdenr)) 
identical(as.character(unique(pdata_balanced$gdenr)),rownames(W_mat_10))
#-> works!

######
# 3. Add spatial lags to data sets
######

#Function to compute spatial lags in panel data
panel_lag <- function(X_vars, listw, data, join_data, join_vars){
  lag_list <- list()
  for (j in 1:length(X_vars)){
    lag_x_all <- data.frame()
    for (i in 1986:2012){
      pdata_W <- data
      pdata_W_year <- filter(pdata_W, year==i)
      XX <- select(pdata_W_year, gdenr, year, x=(X_vars[j]))
      lag_x <- lag.listw(listw, XX$x, NAOK = TRUE)
      XX$lag_x <- lag_x
      lag_x_all <- rbind(lag_x_all, XX)
    }
    names(lag_x_all) <- c("gdenr", "year", X_vars[j], paste0("W_10_", X_vars[j]))
    lag_list[[X_vars[j]]]<-lag_x_all
  }
  df_lag <- data.frame(lag_list[[1]][1:2])
  for(j in 1:length(X_vars)){
    df_lag <- left_join(df_lag, lag_list[[j]], by = c("gdenr", "year"))
  }
  join <- select(join_data, join_vars)
  df_lag <- left_join(df_lag, join, by = c("gdenr", "year"))
  df_lag <- pdata.frame(df_lag, index = c("gdenr", "year"))
  return(df_lag)
}

X_names <- c("p_i_s_1000", "population", "pop_foreign", "pop_young", "pop_old", "pop_secondary", "pop_tertiary", "employment_pc",
             "pop_unemployment", "left", "movie")
idvars <- c("gdenr", "canton", "year", "cyFE", "mFE")
pdata_lag <- panel_lag(X_vars = X_names, listw = W_listw_10, data = pdata_balanced, join_data = pdata_balanced, join_vars = idvars)

###Generate 2nd and 3rd spatial lags
pdata_lag2 <- panel_lag(X_vars = paste0("W_10_",X_names), listw = W_listw_10, data = pdata_lag, join_data = pdata_lag, join_vars = idvars)
pdata_lag2 <- rename_at(pdata_lag2, vars(starts_with("W_10_W_10_")), funs(str_replace(., "W_10_W_10_", "W2_10_")))
pdata_lag3 <- panel_lag(X_vars = paste0("W2_10_",X_names), listw = W_listw_10, data = pdata_lag2, join_data = pdata_lag, join_vars = idvars)
pdata_lag3 <- rename_at(pdata_lag3, vars(starts_with("W_10_W2_10_")), funs(str_replace(., "W_10_W2_10_", "W3_10_")))
#add second and third spatial lag to pdata_lag
pdata_lag <- left_join(pdata_lag, select(pdata_lag2, gdenr, year, starts_with("W2_10")), by = c("gdenr", "year"))
pdata_lag <- left_join(pdata_lag, select(pdata_lag3, gdenr, year, starts_with("W3_10")), by = c("gdenr", "year"))
rm(X_names, idvars, pdata_lag2, pdata_lag3)

###
# 4.) Replication
###

#Specification OLS
model_OLS <- (p_i_s_1000 ~ W_10_p_i_s_1000 + population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + employment_pc 
              + pop_unemployment + left + movie | cyFE + mFE | 0 | mFE + year)

#a1) OLS full sample - no WX
OLS_all <- felm(model_OLS, data = pdata, cmethod ='cgm2')
summary(OLS_all)

#a2) OLS balanced sample - no WX
OLS_balanced <- felm(model_OLS, data = pdata_balanced, cmethod ='cgm2')
summary(OLS_balanced)

#a3) OLS border sample - no WX
OLS_border <- felm(model_OLS, data = pdata_border, cmethod ='cgm2')
summary(OLS_border)

stargazer(OLS_all, OLS_balanced, OLS_border, type = "text")

model_OLS_WX <- (p_i_s_1000 ~ W_10_p_i_s_1000 + population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + 
                   employment_pc + pop_unemployment + left + movie + W_10_population + W_10_pop_foreign + W_10_pop_young + 
                   W_10_pop_old + W_10_pop_secondary + W_10_pop_tertiary + W_10_employment_pc + W_10_pop_unemployment + 
                   W_10_left + W_10_movie | cyFE + mFE | 0 | mFE + year)

#a4) OLS full sample WX
OLS_all_WX <- felm(model_OLS_WX, data = pdata, cmethod ='cgm2')
summary(OLS_all_WX)

#a5) OLS balanced sample WX
OLS_balanced_WX <- felm(model_OLS_WX, data = pdata_balanced, cmethod ='cgm2')
summary(OLS_balanced_WX)

#a6) OLS border sample WX
OLS_border_WX <- felm(model_OLS_WX, data = pdata_border, cmethod ='cgm2')
summary(OLS_border_WX)

stargazer(OLS_all_WX, OLS_balanced_WX, OLS_border_WX, type = "text")

#a7) OLS balanced sample - no WX with replicated W
OLS_balanced_Wnew <- felm(model_OLS, data = pdata_lag, cmethod ='cgm2')
summary(OLS_balanced_Wnew)

#a8) OLS balanced sample WX with replicated W
OLS_balanced_WX_Wnew <- felm(model_OLS_WX, data = pdata_lag, cmethod ='cgm2')
summary(OLS_balanced_WX_Wnew)

stargazer(OLS_balanced, OLS_balanced_Wnew, OLS_balanced_WX, OLS_balanced_WX_Wnew, type = "text")


#Specification S2SLS
model_S2SLS <- (p_i_s_1000 ~ population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + employment_pc 
                + pop_unemployment + left + movie | 
                  cyFE + mFE | 
                  (W_10_p_i_s_1000 ~ W_10_population + W_10_pop_foreign + W_10_pop_young + W_10_pop_old + W_10_pop_secondary 
                   + W_10_pop_tertiary + W_10_employment_pc + W_10_pop_unemployment + W_10_left + W_10_movie) | 
                  mFE + year)

#b1) S2SLS full sample
S2SLS_all <- felm(model_S2SLS, data = pdata, cmethod ='cgm2')
summary(S2SLS_all)

#b2) S2SLS balanced sample
S2SLS_balanced <- felm(model_S2SLS, data = pdata_balanced, cmethod ='cgm2')
summary(S2SLS_balanced)

#b3) S2SLS border sample
S2SLS_border <- felm(model_S2SLS, data = pdata_border, cmethod ='cgm2')
summary(S2SLS_border)

stargazer(S2SLS_all, S2SLS_balanced, S2SLS_border, type = "text")

#Table 2 of the paper
stargazer(OLS_all, OLS_border, S2SLS_border, type = "text")

#b4 S2SLS balanced sample with replicated W
S2SLS_balanced_Wnew <- felm(model_S2SLS, data = pdata_lag, cmethod ='cgm2')
summary(S2SLS_balanced_Wnew)

#table for paper
stargazer(OLS_all, S2SLS_all, OLS_balanced, S2SLS_balanced, OLS_balanced_Wnew, S2SLS_balanced_Wnew, OLS_border, S2SLS_border, keep = c("W_10_p_i_s_1000", "W_10_p_i_s_1000(fit)"), keep.stat=c("n","adj.rsq"), star.cutoffs = c(0.05, 0.01, 0.001))

#b5 S2SLS balanced sample with replicated W - two lags
model_iv2 <- as.formula(p_i_s_1000 ~ population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + employment_pc + 
                          pop_unemployment + left + movie | cyFE + mFE | 
                          (W_10_p_i_s_1000 ~ W_10_population + W_10_pop_foreign + W_10_pop_young + W_10_pop_old + W_10_pop_secondary +
                             W_10_pop_tertiary + W_10_employment_pc + W_10_pop_unemployment + W_10_left + W_10_movie + W2_10_population +
                             W2_10_pop_foreign + W2_10_pop_young + W2_10_pop_old + W2_10_pop_secondary + W2_10_pop_tertiary + 
                             W2_10_employment_pc + W2_10_pop_unemployment + W2_10_left + W2_10_movie) | mFE + year)
S2SLS_balanced_Wnew_lag2 <- felm(model_iv2, data = pdata_lag, cmethod ='cgm2')
summary(S2SLS_balanced_Wnew_lag2)


#b5 S2SLS balanced sample with replicated W - three lags
model_iv3 <- as.formula(p_i_s_1000 ~ population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + employment_pc 
                        + pop_unemployment + left + movie | 
                          cyFE + mFE | (W_10_p_i_s_1000 ~ W_10_population + W_10_pop_foreign + W_10_pop_young + W_10_pop_old + 
                                          W_10_pop_secondary + W_10_pop_tertiary + W_10_employment_pc + W_10_pop_unemployment + 
                                          W_10_left + W_10_movie + W2_10_population + W2_10_pop_foreign + W2_10_pop_young + 
                                          W2_10_pop_old + W2_10_pop_secondary + W2_10_pop_tertiary + W2_10_employment_pc + 
                                          W2_10_pop_unemployment + W2_10_left + W2_10_movie + W3_10_population + W3_10_pop_foreign + 
                                          W3_10_pop_young + W3_10_pop_old + W3_10_pop_secondary + W3_10_pop_tertiary + W3_10_employment_pc + 
                                          W3_10_pop_unemployment + W3_10_left + W3_10_movie) | mFE + year)

S2SLS_balanced_Wnew_lag3 <- felm(model_iv3, data = pdata_lag, cmethod ='cgm2')
summary(S2SLS_balanced_Wnew_lag3)

stargazer(S2SLS_balanced_Wnew, S2SLS_balanced_Wnew_lag2, S2SLS_balanced_Wnew_lag3, type = "text")

#c) IVP border sample
model_IVP <- (p_i_s_1000 ~ population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + employment_pc 
              + pop_unemployment + left + movie + W_10_population + W_10_pop_foreign + W_10_pop_young + W_10_pop_old 
              + W_10_pop_secondary + W_10_pop_tertiary + W_10_employment_pc + W_10_pop_unemployment + W_10_left + W_10_movie | 
                cyFE + mFE | 
                (W_10_p_i_s_1000 ~ W_10_z_i_s_1000_l) | 
                mFE + year)

IVP_border <- felm(model_IVP, data = pdata_border, cmethod ='cgm2')
summary(IVP_border)

#f) IVP full sample without WX
model_IVP_noWX <- (p_i_s_1000 ~ population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + employment_pc 
                   + pop_unemployment + left + movie | 
                     cyFE + mFE | 
                     (W_10_p_i_s_1000 ~ W_10_z_i_s_1000_l) | 
                     mFE + year)

IVP_border_noWX <- felm(model_IVP_noWX, data = pdata_border, cmethod ='cgm2')
summary(IVP_border_noWX)

#table for paper
stargazer(IVP_border, IVP_border_noWX, S2SLS_balanced_Wnew, S2SLS_balanced_Wnew_lag2, S2SLS_balanced_Wnew_lag3, IVP_border, IVP_border_noWX,
          keep.stat=c("n","adj.rsq"), star.cutoffs = c(0.05, 0.01, 0.001))

###
#5.) Extension spml
###

#a) individual and cyFE via Within transformation - no WX
fm_a <- p_i_s_1000 ~ Within(population, cyFE) +  Within(pop_foreign, cyFE) +  Within(pop_young, cyFE) +  
  Within(pop_old, cyFE) +  Within(pop_secondary, cyFE) +  Within(pop_tertiary, cyFE) + 
  Within(employment_pc, cyFE) +  Within(pop_unemployment, cyFE) +  Within(left, cyFE) +  Within(movie, cyFE)

spml_a <- spml(fm_a, data = pdata_lag, index = c("gdenr", "year"),
               listw = W_listw_10, lag=TRUE, spatial.error="b",
               model="within", effect="individual")
summary(spml_a)

#b) individual and cyFE via panel dimensions - no WX
p.df.pdf.cy <- pdata.frame(data.frame(pdata_lag), index = c("gdenr", "cyFE"))
fm_b <- p_i_s_1000 ~ population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + 
  employment_pc + pop_unemployment + left + movie

spml_b <- spml(fm_b, data = p.df.pdf.cy, index = c("gdenr", "cyFE"),
               listw = W_listw_10, lag=TRUE, spatial.error="b",
               model="within", effect="twoway")
summary(spml_b)

#c) individual and cyFE via Within transformation - including WX
fm_c <- p_i_s_1000 ~ Within(population, cyFE) +  Within(pop_foreign, cyFE) +  Within(pop_young, cyFE) +  
  Within(pop_old, cyFE) +  Within(pop_secondary, cyFE) +  Within(pop_tertiary, cyFE) + 
  Within(employment_pc, cyFE) +  Within(pop_unemployment, cyFE) +  Within(left, cyFE) +  Within(movie, cyFE) +
  Within(W_10_population, cyFE) +  Within(W_10_pop_foreign, cyFE) +  Within(W_10_pop_young, cyFE) +  
  Within(W_10_pop_old, cyFE) +  Within(W_10_pop_secondary, cyFE) +  Within(W_10_pop_tertiary, cyFE) + 
  Within(W_10_employment_pc, cyFE) +  Within(W_10_pop_unemployment, cyFE) +  Within(W_10_left, cyFE) +  Within(W_10_movie, cyFE)

spml_c <- spml(fm_c, data = pdata_lag, index = c("gdenr", "year"),
               listw = W_listw_10, lag=TRUE, spatial.error="b",
               model="within", effect="individual")
summary(spml_c)

#d) individual and cyFE via panel dimensions - including WX
fm_d <- p_i_s_1000 ~ population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + 
  employment_pc + pop_unemployment + left + movie + W_10_population + W_10_pop_foreign + W_10_pop_young + 
  W_10_pop_old + W_10_pop_secondary + W_10_pop_tertiary + W_10_employment_pc + W_10_pop_unemployment + 
  W_10_left + W_10_movie

spml_d <- spml(fm_d, data = p.df.pdf.cy, index = c("gdenr", "cyFE"),
               listw = W_listw_10, lag=TRUE, spatial.error="b",
               model="within", effect="twoway")
summary(spml_d)

#update 20.02.20:
is.pbalanced(p.df.pdf.cy) #-> FALSE!!
#This might be the reason for the misleading result, panel is not balanced anymore if we use cyFE as year
#BUT: why does it still run through and there is no error message?

#Show that approach id = ("gdenr", "cyFE") works for the coefficient (not for Standard error) (plm supports unbalanced panels)
fm_a2 <- p_i_s_1000 ~ W_10_p_i_s_1000 + population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + 
  employment_pc + pop_unemployment + left + movie
plm_a <- plm(fm_a2, data = p.df.pdf.cy, index = c("gdenr", "cyFE"), model="within", effect="twoways")
summary(plm_a)
#-> exactly the same coefficient as line 185 #a2) OLS balanced sample - no WX
summary(spml_a)




#summary table
library(summarytools)
library(stringi)
library(xtable)
pdata_balanced <- pdata.frame(pdata_balanced, index = c("gdenr", "year"))
ols_all1 <- lm(p_i_s_1000 ~ W_10_p_i_s_1000 + population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + employment_pc 
               + pop_unemployment + left + movie, pdata)
ols_all2 <- lm(p_i_s_1000 ~ W_10_p_i_s_1000 + population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + employment_pc 
               + pop_unemployment + left + movie, pdata_balanced)
ols_all3 <- lm(p_i_s_1000 ~ W_10_p_i_s_1000 + population + pop_foreign + pop_young + pop_old + pop_secondary + pop_tertiary + employment_pc 
               + pop_unemployment + left + movie, pdata_border)

x1 <- data.frame(gdenr = (str_sub(names(ols_all1$fitted.values), end=-6)), year = (stri_sub(names(ols_all1$fitted.values), -4)))
x2 <- data.frame(gdenr = (str_sub(names(ols_all2$fitted.values), end=-6)), year = (stri_sub(names(ols_all2$fitted.values), -4)))
x3 <- data.frame(gdenr = (str_sub(names(ols_all3$fitted.values), end=-6)), year = (stri_sub(names(ols_all3$fitted.values), -4)))

xx1 <- left_join(x1, pdata, by = c("gdenr","year"))
xx2 <- left_join(x2, pdata, by = c("gdenr","year"))
xx3 <- left_join(x3, pdata, by = c("gdenr","year"))

df_select1 <- select(xx1, p_i_s_1000, population, pop_foreign, pop_young, pop_old,
                     pop_secondary, pop_tertiary, employment_pc, pop_unemployment, left, movie)
df_select2 <- select(xx2, p_i_s_1000, population, pop_foreign, pop_young, pop_old,
                     pop_secondary, pop_tertiary, employment_pc, pop_unemployment, left, movie)
df_select3 <- select(xx3, p_i_s_1000, population, pop_foreign, pop_young, pop_old,
                     pop_secondary, pop_tertiary, employment_pc, pop_unemployment, left, movie)

sum_tab1 <- descr(df_select1, stats = c("mean", "sd"), transpose = TRUE)
sum_tab2 <- descr(df_select2, stats = c("mean", "sd"), transpose = TRUE)
sum_tab3 <- descr(df_select3, stats = c("mean", "sd"), transpose = TRUE)

sum_tab1
sum_tab2
sum_tab3

#Combine tables
sum_tab <- cbind(sum_tab1, sum_tab2, sum_tab3)
sum_tab

#Adjust order of rows
ord <- c("p_i_s_1000", "population", "pop_foreign", "pop_young", "pop_old", "pop_secondary", 
         "pop_tertiary", "employment_pc", "pop_unemployment", "left", "movie")

sum_tab <- sum_tab[match(ord, rownames(sum_tab)), ]
rm(countries, ord)

#Print and export to LaTex
print(sum_tab)
xtable(sum_tab)

length(unique(xx1$gdenr))
length(unique(xx2$gdenr))
length(unique(xx3$gdenr))

length(unique(xx1$year))
length(unique(xx2$year))
length(unique(xx3$year))