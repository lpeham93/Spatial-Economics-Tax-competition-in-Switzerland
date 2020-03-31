# Generate weights matrix that Parchet used
#library(tidyr)
library(dplyr)
library(spdep)

#load data
load("data.RData")
df <- data[[1]]

#load weights data
load("wdf.RData")

#Distance should be zero only between the municipality itself
a <- w.df %>% filter(dist_m == 0)
ind <- a %>% transmute(check = ifelse(gdenr_fr == gdenr_to, F, T))
View(a[ind[,1],]) # these municipalities have neighbours with zero distance!!
gd_fail <- a[ind[,1],]$gdenr_fr
check_fails <- df %>% filter(gdenr %in% gd_fail, year %in% c(1983, 2012))

#for some reason that I don't understand, there are municipalities which have
#neighours 0 meters away
#Maybe unification / separation of municipalities over time is a reason

#still include this municipality relations further on
#and get rid of the rows where the gdenr_from == gdenr_to
#w.df <- w.df %>% filter(gdenr_fr != gdenr_to) #that way we keep the weird zero distance neighbours

#exclude all zero distance observations
w.df <- w.df %>% filter(dist_m > 0)

w5 <- w.df %>% filter(dist_m <=5)
w10 <- w.df %>% filter(dist_m <=10)

#make neighbourhood lists
nb <- as.list(unique(df$gdenr))
names(nb) <- unique(df$gdenr)

get.nb <- function(nr, dat){
  sub <- dat %>% filter(gdenr_fr == nr)
  if(nrow(sub) == 0){return(as.vector(0L))}
  as.integer(unique(sub$gdenr_to))
}

nb_5 <- lapply(nb, get.nb, dat=w5)
class(nb_5) <- "nb"
nb_10 <- lapply(nb, get.nb, dat=w10)
class(nb_10) <- "nb"
#save(nb_5, nb_10, file = "nblists_parchet.RData")

#I think the next question is: "Which municipalities do we actually use for the regression?"
#then we subset the nb lists for only the ones we actually use

#get islands; just as an information; there are a lot of them!!
cond <- lapply(nb_5, function(x) length(x)==1 & all(x>0))
cond <- unlist(cond)
islands_5 <- as.integer(names(nb_5)[cond])

cond <- lapply(nb_10, function(x) length(x)==1 & all(x>0))
cond <- unlist(cond)
islands_10 <- as.integer(names(nb_10)[cond])

#W_mat - Zero.policy Error!! Don't know whats going wrong!! Maybe too many islands
wlist_5 <- nb2listw(nb_5, style = "W", zero.policy = TRUE)
Wmat_5 <- listw2mat(wlist_5)

#Compare wLists
load("gdenr_fabi.Rdata")
gdenr_luki <- as.numeric(names(nb_5))
notinw.df <- setdiff(gdenr_luki, gdenr_fabi)

View(df %>% filter(gdenr %in% notinw.df))

check <- df %>% select(c(gdenr, nb_5))
check <- unique(check)

f <- function(x){
  if (length(x)==1 && x[1]==0L){return(NA)}
  length(x)
}
c.list <- sapply(nb_5, f) #islands appear as having 1 neighbour!!!
c.list <- data.frame(gdenr = as.numeric(names(c.list)), c.list)
check <- left_join(check, c.list)
ind <- check$nb_5 == check$c.list

nrow(check[is.na(ind),]) #95