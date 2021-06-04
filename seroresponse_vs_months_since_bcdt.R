# re-analysis of Achiron et al 2021 
# Humoral immune response to COVID-19 mRNA vaccine in patients with multiple sclerosis treated with high-efficacy disease-modifying therapies
# https://journals.sagepub.com/doi/10.1177/17562864211012835
# and
# Deepak et al 2021
# Glucocorticoids and B Cell Depleting Agents Substantially Impair Immunogenicity of mRNA Vaccines to SARS-CoV-2
# https://www.medrxiv.org/content/10.1101/2021.04.05.21254656v2
#
# convalescent sera titers corresponding to Deepak from 
# Turner et al 2021 SARS-CoV-2 infection induces long-lived bone marrow plasma cells in humans
# https://www.nature.com/articles/s41586-021-03647-4
#
# b-cell data from Baker et al 2021
# COVID-19 vaccine-readiness for ocrelizumab and other anti-CD20-depleting therapies in multiple sclerosis and other autoimmune diseases
# https://www.authorea.com/users/334776/articles/461938-covid-19-vaccine-readiness-for-ocrelizumab-and-other-anti-cd20-depleting-therapies-in-multiple-sclerosis-and-other-autoimmune-diseases
#
# COVID vaccine seroresponse probability as function of time and b-cell population size since last Ocrevus dose
#
# Mike Famulare 28 May 2021
# all opinions and analysis my own.
#

library(tidyverse)
library(mgcv)
library(betareg)
library(matrixStats)

set.seed(10)

# helper functions
logit <- function(p){log(p/(1-p))}
invlogit <- function(lo){exp(lo)/(1+exp(lo))}


# load digitized data on covid vaccine seroresponse
d <- read.table('BCDT_covid_vac_seroconversion_data.csv',sep=',',header=TRUE, stringsAsFactors = TRUE)
  
  # add numeric outcome variable
  d$outcome=as.numeric(d$serostatus)-1
  
  # edit label for later plot 
  levels(d$serostatus)<-paste('vaccine-',levels(d$serostatus),sep='')
  

# logistic regression
summary(s_model <- glm(outcome ~ months, data=d, family='binomial')) #intercept insignificantly different by source (p=0.5)

# fit spline model
summary(s_model2 <- gam(outcome ~ s(months), data=d, family='binomial'))

  # construct predictor with intervals
  pred_dat<-data.frame(months=seq(0,15,b=0.2))
  
  tmp<-predict(s_model,newdata=pred_dat,se=TRUE)
  
  pred_dat$fit <- invlogit(tmp$fit)
  pred_dat$fit_lower <- invlogit(tmp$fit-2*tmp$se.fit)
  pred_dat$fit_upper <- invlogit(tmp$fit+2*tmp$se.fit)
  
  pred_dat$fit_lower20 <- invlogit(tmp$fit-0.84*tmp$se.fit)
  pred_dat$fit_upper80 <- invlogit(tmp$fit+0.84*tmp$se.fit)
  
  
  # for the spline too
  tmp<-predict(s_model2,newdata=pred_dat,se=TRUE)
  
  pred_dat$fit_gam <- invlogit(tmp$fit)
  pred_dat$fit_gam_lower <- invlogit(tmp$fit-2*tmp$se.fit)
  pred_dat$fit_gam_upper <- invlogit(tmp$fit+2*tmp$se.fit)
  
ggplot() +
  geom_ribbon(data=pred_dat,aes(x=months,ymin=fit_lower,ymax=fit_upper),alpha=0.2) +
  geom_ribbon(data=pred_dat,aes(x=months,ymin=fit_lower20,ymax=fit_upper80),alpha=0.2) +
  geom_line(data=pred_dat,aes(x=months, y=fit))+
  geom_point(data=d,aes(x=months,y=outcome,color=source)) +
  scale_x_continuous(breaks=seq(3,15,by=3)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  ylab('seroconversion probability') + xlab('months between BCDT and first COVID vax') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  theme(panel.border=element_blank(), 
        axis.line = element_line(), 
        axis.ticks = element_line(colour='black')) 

ggsave('covid_vax_seroconversion_probability_vs_time_since_BCDT.png',units='in',width=6,height=3)


# with gam
ggplot() +
  geom_ribbon(data=pred_dat,aes(x=months,ymin=fit_lower,ymax=fit_upper),alpha=0.2) +
  geom_ribbon(data=pred_dat,aes(x=months,ymin=fit_lower20,ymax=fit_upper80),alpha=0.2) +
  geom_line(data=pred_dat,aes(x=months, y=fit))+
  geom_point(data=d,aes(x=months,y=outcome,color=source)) +
  geom_line(data=pred_dat,aes(x=months, y=fit_gam),color='darkblue')+
  scale_x_continuous(breaks=seq(3,15,by=3)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  ylab('seroconversion probability') + xlab('months between BCDT and first COVID vax') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  theme(panel.border=element_blank(), 
        axis.line = element_line(), 
        axis.ticks = element_line(colour='black')) 

ggsave('covid_vax_seroconversion_probability_vs_time_since_BCDT_with_gam.png',units='in',width=6,height=3)



# load data on B-cell restoration after Ocvrevus or rituximab
b <- read.table('Baker2020 - COVID-19 vaccine-readiness for ocrelizumab - Fig3.csv',sep=',',header=TRUE, stringsAsFactors = TRUE)
b2 <- read.table('Ellrichmann2019_Fig1a.csv',sep=',',header=TRUE, stringsAsFactors = TRUE) %>%
  filter(cell!='cut_axis_200') %>%
  filter(months>1)

# select reconstituting population
# I think it's preferable to include both cell subsets, which parallel each other,
# for a more complete picture of uncertainty.  But it's more common to measure just cd19+.
# Qualitative conclusions do not depend on this choice since they parallel each other.
# (of course, individual-level data on seroconversion and cd19+ b-cell reconstitution
# would be much better. Someone should do a study!)
b <- b %>% filter(cell!="Memory B-cells") %>% filter(months>1) # include both CD19 B cell and naive b-cell data
# b <- b %>% filter(cell=='CD19 B-cells') %>% filter(months>1) # just cd19 b cell data

# join tables

b <- b %>% select(months,fraction_of_baseline,cell) %>%
  rbind(b2 %>% select(months,fraction_of_baseline,cell))



# explore how B-cell reconstitution maps onto response
# reformat data and the previous model for plotting together
plot_dat <- b %>%
  select(months,fraction_of_baseline,cell) %>%
  rbind(pred_dat %>% mutate(fraction_of_baseline=fit,
                            cell='probability of seroconversion') %>%
          select(months,fraction_of_baseline,cell)
  ) %>%
  mutate(data=cell) %>%
  select(months,fraction_of_baseline,data)

# reformatting when using cd19+ and naive b-cell data
plot_dat$data=factor(plot_dat$data,levels=c('probability of seroconversion','Naive B-cells','CD19 B-cells',levels(b2$cell)))
levels(plot_dat$data)=c('seroconversion probability','mean naive B-cell count relative to baseline','mean cd19+ B-cell count relative to baseline',levels(b2$cell))

# reformating when using cd19+ data only
# plot_dat$data=factor(plot_dat$data,levels=c('probability of seroconversion','CD19 B-cells'))
# levels(plot_dat$data)=c('seroconversion probability','mean cd19+ B-cell count relative to baseline')


ggplot() +
  geom_ribbon(data=pred_dat,aes(x=months,ymin=fit_lower,ymax=fit_upper),alpha=0.2) +
  geom_ribbon(data=pred_dat,aes(x=months,ymin=fit_lower20,ymax=fit_upper80),alpha=0.2) +
  geom_line(data=plot_dat,aes(x=months,y=fraction_of_baseline,group=data,linetype=data)) +  
  scale_x_continuous(breaks=seq(3,21,by=3)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  ylab('fraction') + xlab('months between BCDT and first COVID vax') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  theme(panel.border=element_blank(), 
        axis.line = element_line(), 
        axis.ticks = element_line(colour='black')) 

ggsave('covid_vax_seroconversion_probability_and_b_cell_reconstitution.png',units='in',width=6,height=3)

# seroconversion time parallels mean b-cell repopulation time, as expected!


# smooth extrapolation of repopulation time data 
# fit beta regression for b-cell repopulation vs time since the data are noisy and missing early values
summary(b_model <- betareg('fraction_of_baseline ~ months',data=b))

# fit logit-linear regression as sanity check. I'm surprised the variance is as high as it is.
# summary(b_model_2 <- lm('logit(fraction_of_baseline) ~ months',data=b))


# smooth predictor
b_pred <- data.frame(months=seq(0,30,b=0.2))

b_pred$fit <- predict(b_model,newdata=b_pred,type='response')
b_pred$fit_lower <- predict(b_model,newdata=b_pred,type='quantile',at=0.025)
b_pred$fit_upper <- predict(b_model,newdata=b_pred,type='quantile',at=0.975)
b_pred$fit_lower20 <- predict(b_model,newdata=b_pred,type='quantile',at=0.2)
b_pred$fit_upper80 <- predict(b_model,newdata=b_pred,type='quantile',at=0.8)

# variance checks out!  the jaggedness of the naive b-cell data induces a large variance
# tmp <- predict(b_model_2,newdata=b_pred,type='response',se=TRUE)
# b_pred$fit_2 <- invlogit(tmp$fit)
# b_pred$fit_2_lower <- invlogit(tmp$fit-2*tmp$se.fit) 
# b_pred$fit_2_upper <- invlogit(tmp$fit+2*tmp$se.fit)

ggplot() +
  geom_ribbon(data=b_pred,aes(x=months,ymin=fit_lower,ymax=fit_upper),alpha=0.2) +
  geom_ribbon(data=b_pred,aes(x=months,ymin=fit_lower20,ymax=fit_upper80),alpha=0.2) +
  # geom_ribbon(data=b_pred,aes(x=months,ymin=fit_2_lower,ymax=fit_2_upper),alpha=0.2,fill='darkgreen') +
  geom_line(data=b,aes(x=months,y=fraction_of_baseline,group=cell,color=cell)) +  
  geom_line(data=b_pred,aes(x=months, y=fit))+
  # geom_line(data=b_pred,aes(x=months, y=fit_2),color='darkgreen')+
  scale_x_continuous(breaks=seq(3,30,by=3)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  ylab('fraction of baseline') + xlab('months since last BCDT') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  theme(panel.border=element_blank(), 
        axis.line = element_line(), 
        axis.ticks = element_line(colour='black')) 

ggsave('mean_b-cell_repopulation_vs_time_since_BCDT.png',units='in',width=6,height=3)



#####
# translate seroconversion probability vs time to probability vs b-cell
# the idea here is that there's a clear signal of COVID vax seroconversion vs time
# and a clear correlate of b-cell population size relative to baseline that parallels 
# seroconversion in time.
# Thus in lieu of having individual-level b-cell population size and seroconversion data,
# I hypothesize that time is a proxy for b-cell population size, and use a little 
# arithmetic to relate the seroconversion model to b-cell population.
# This can give an estimate of the individual-level correlate of seroconversion probability
# that can be used to optimize BCDT and vaccine dose timing.


# construct model
seroconversion_vs_b_cell <- function(newdata,s_model_coef,b_model_coef){
  invlogit(s_model_coef[1] + 
             s_model_coef[2]*( # months coefficient
             (logit(newdata$relative_to_baseline)-b_model_coef[1])/b_model_coef[2]) # this inverts b-cell population relative to mean vs months
          )
}

# prediction range (can't quite go to 1 because of beta regression assumption, but that's fine for this purpose)
pred_dat<-data.frame(relative_to_baseline= seq(0.001,0.999,by=0.001))

# MLE
pred_dat$fit <- seroconversion_vs_b_cell(pred_dat,coef(s_model),coef(b_model))

# confidence interval from parametric bootstrap
n_rep <- 1e3
mod_samples <- matrix(NA,nrow=nrow(pred_dat),ncol=n_rep)
colnames(mod_samples) <- paste('replicate_',1:n_rep,sep='')

for(k in 1:n_rep){
  
  # sample model parameters
  s_coef <- coef(s_model) + chol(vcov(s_model))%*%rnorm(2)
  b_coef <- coef(b_model)[1:2] + chol(vcov(b_model)[1:2,1:2])%*%rnorm(2)
  
  # model curve
  mod_samples[,k] = seroconversion_vs_b_cell(pred_dat,s_coef,b_coef)
}

# convex hull CI
pred_dat$fit_lower <- rowQuantiles(mod_samples,probs=0.025)
pred_dat$fit_upper <- rowQuantiles(mod_samples,probs=0.975)

pred_dat$fit_lower20 <- rowQuantiles(mod_samples,probs=0.2)
pred_dat$fit_upper80 <- rowQuantiles(mod_samples,probs=0.8)

# join in model trajectories
pred_dat_2 <- pred_dat %>% select(relative_to_baseline,fit) %>%
  cbind(mod_samples[,1:50]) %>% 
  gather(key='model',value='probability',-relative_to_baseline)

ggplot() +
  geom_ribbon(data=pred_dat,aes(x=relative_to_baseline,ymin=fit_lower,ymax=fit_upper),alpha=0.2) +
  geom_ribbon(data=pred_dat,aes(x=relative_to_baseline,ymin=fit_lower20,ymax=fit_upper80),alpha=0.2) +
  # geom_line(data=pred_dat_2 %>% filter(model!='fit'),aes(x=relative_to_baseline, y=probability,group=model),alpha=0.1)+
  geom_line(data=pred_dat,aes(x=relative_to_baseline, y=fit))+
  scale_x_continuous(breaks=seq(0,1,by=0.1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  ylab('seroconversion probability') + xlab('CD19+ B-cell population relative to baseline') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  theme(panel.border=element_blank(), 
        axis.line = element_line(), 
        axis.ticks = element_line(colour='black')) 

ggsave('covid_vax_seroconversion_probability_vs_CD19+_population_relative_to_baseline.png',units='in',width=6,height=3)




## look at titer levels for seroconverters and time since bcdt from deepak
# 
t_data <- d %>% filter(!is.na(titer)) %>% filter(serostatus=='vaccine-positive') # there is no significant difference in slope if the negatives are included, but it's clear that negatives after a few months are different than positives for reasons unknown

# convalescent plasma from Turner et al 2021 for comparison
c_data <- read.table('Turner2021_Fig1b.csv',sep=',',header = TRUE, stringsAsFactors = TRUE) 

# first timepoint is average 44 days since symptom onset, so that is closest comparator for post-vaccination titers in Deepak
c_data %>% group_by(timepoint) %>% summarize(mean_collection_time=mean(days_post_symptom_onset))

c_data <- c_data %>% filter(timepoint==1)


# log-linear
summary(t_model <- lm(log(titer) ~ months, data=t_data)) 

# construct predictor with intervals
pred_dat<-data.frame(months=seq(0,15,b=0.2))

tmp<-predict(t_model,newdata=pred_dat,se=TRUE)

pred_dat$fit <- exp(tmp$fit)
pred_dat$fit_lower <- exp(tmp$fit-2*tmp$se.fit)
pred_dat$fit_upper <- exp(tmp$fit+2*tmp$se.fit)
pred_dat$fit_lower20 <- exp(tmp$fit-0.84*tmp$se.fit)
pred_dat$fit_upper80 <- exp(tmp$fit+0.84*tmp$se.fit)

c_plot_data <- data.frame(x=17,
                          ymin=quantile(c_data$titer,0.025),
                          lower=quantile(c_data$titer,0.20),
                          middle=median(c_data$titer),
                          upper=quantile(c_data$titer,0.80),
                          ymax=quantile(c_data$titer,0.975)
                          )

ggplot() +
  geom_hline(aes(yintercept=40),linetype='dashed')+
  geom_ribbon(data=pred_dat,aes(x=months,ymin=fit_lower,ymax=fit_upper),alpha=0.2) +
  geom_ribbon(data=pred_dat,aes(x=months,ymin=fit_lower20,ymax=fit_upper80),alpha=0.2) +
  geom_line(data=pred_dat,aes(x=months, y=fit))+
  geom_point(data=d ,aes(x=months,y=titer,color=serostatus)) +
  geom_boxplot(data=c_plot_data, aes(color='convalescent',x=x,ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax),stat='identity')+
  scale_x_continuous(breaks=seq(0,15,by=3)) +
  scale_y_continuous(trans='log', breaks=10^(0:4),limits=c(1,4e4))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  ylab('Anti-S IgG Titer') + xlab('months between BCDT and first COVID vax') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  theme(panel.border=element_blank(), 
        axis.line = element_line(), 
        axis.ticks = element_line(colour='black')) 

ggsave('covid_vax_IgG_titer_vs_time_since_BCDT.png',units='in',width=6,height=3)



# compare slope on months between seroconversion and b-cell repopulation
slope_pred <- data.frame(model=factor(c('seroconversion','anti-S IgG titer','b-cell repopulation'),levels=c('seroconversion','anti-S IgG titer','b-cell repopulation')),
                         months_coef = c(coef(s_model)[2],coef(t_model)[2],coef(b_model)[2]),
                         months_se=c(vcov(s_model)[2,2],vcov(t_model)[2,2],vcov(b_model)[2,2])) 
ggplot(slope_pred) +
  geom_errorbar(aes(x=model,ymin=months_coef-2*months_se,ymax=months_coef+2*months_se),width=0.25)+
  geom_point(aes(x=model,y=months_coef)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  ylab('months coefficient') + xlab('model') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border=element_blank(), 
    axis.line = element_line(), 
    axis.ticks = element_line(colour='black')) 

ggsave('months_coef_seroresponse_and_b-cell_repopulation.png',units='in',width=4,height=3)

# b-cell repopulation slope is slower than seroconversion slope
# would seem there is nonlinearity to seroconversion vs mean b-cell pop, such that
# seroconversion recovers faster than b-cell population. This wouldn't surprise me,
# but I don't have the expertise to know for sure.


