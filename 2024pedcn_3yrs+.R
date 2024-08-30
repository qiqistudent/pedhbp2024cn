######################################################################################
### DATA PREPERATION
#Change variable names as follows: height as Ht, sbp as msbp, dbp as mdbp
#Trunc age and name as age_trun, and categorize sex as 1 for male and 2 for female
#Download 'ht.lms.csv' and 'bp.refcn.csv' to a basic directory [dir] where you want to use.
######################################################################################


#Import heigt LMS data 
ht.ref<-read.csv('ht.lms.csv')
d1<-merge(data,ht.ref,by=c("age_trun","sex"),all.x=T)

#Adjustment for extreme values
d1$sd2ps<-with(d1,ht_m*((1+2*ht_l*ht_s)^(1/ht_l)))
d1$sd3ps<-with(d1,ht_m*((1+3*ht_l*ht_s)^(1/ht_l)))
d1$sd3ng<-with(d1,ht_m*((1-3*ht_l*ht_s)^(1/ht_l)))
d1$sd2ng<-with(d1,ht_m*((1-2*ht_l*ht_s)^(1/ht_l)))

#Calculate the zscore for height
d1$zht<-with(d1,ifelse(((Ht/ht_m)^ht_l-1)/(ht_s*ht_l)>3,3+(Ht-sd3ps)/(sd3ps-sd2ps),
                     ifelse(((Ht/ht_m)^ht_l-1)/(ht_s*ht_l)< -3,-3-(sd3ng-Ht)/(sd2ng-sd3ng),
                            ((Ht/ht_m)^ht_l-1)/(ht_s*ht_l))))

#Convert to height percentile categories as p5,p10,p25,p50,p75,p90,p95

d1$htpct<-with(t,ifelse(zht< -1.4632,5,# -1.4632 is midpoint between Zp5 and Zp10
                        ifelse(zht< -0.97802,10,# -0.97802 is midpoint between Zp10 and Zp25
                               ifelse(zht< -0.33724,25,# -0.33724 is midpoint between Zp25 and Zp50
                                      ifelse(zht< 0.33724,50,# -0.33724 is midpoint between Zp25 and Zp50
                                             ifelse(zht< 0.97802,75,# 0.97802 is midpoint between Zp75 and Zp90
                                                    ifelse(zht< 1.4632,90,# 1.4632 is midpoint between Zp90 and Zp95
                                                           95)))))))
#Import heigt-specific bp data
bp.ref<-read.csv('bp.refcn.csv')
d2<-merge(d1,bp.refcn,by=c("age_trun","sex","htpct"),all.x=T)

#bp grade:1 for normotension;2 for high normal bp; 3 for stage 1/mild hypertension;4 for stage 2/moderate to severe hypertension
d2$hbp.ch<-with(d2,ifelse(msbp<sbp90 &mdbp<dbp90,1,ifelse((sbp90<=msbp & msbp<sbp95) | (dbp90<=mdbp & mdbp<dbp95),2,ifelse((sbp95<=msbp & msbp<sbp99+5) | (dbp95<=mdbp & mdbp<dbp99+5),3,4))))
d2$hbp.ad<-with(d2,ifelse(msbp<120 &mdbp<80,1,ifelse((120<=msbp & msbp<140) | (80<=mdbp & mdbp<90),2,ifelse((140<=msbp & msbp<160) | (90<=mdbp & mdbp<100),3,4))))

library(dplyr)
d.new<- d2%>% rowwise() %>% mutate(hbpcn.g4 = max(hbp.ch,hbp.ad))

#bp grade:1 for normotension;2 for high normal bp; 3 for stage 1/mild hypertension;4 for stage 2/moderate to severe hypertension
d.new$hbpcn.g4<-with(d.new,ifelse(age_trun<=16,hbpcn.g4,hbp.ad))


