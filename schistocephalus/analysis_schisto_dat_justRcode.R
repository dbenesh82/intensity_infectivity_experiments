## ---- message=FALSE, warning=FALSE---------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(boot)
library(MuMIn)

sdat <- read.csv(file = "Schisto_double_exposure.csv", header = TRUE, sep = ',')
head(sdat)

## ---- message=TRUE, warning=TRUE-----------------------------------------
sdat <- select(sdat, worm.fam = Clutch, trt = Treatment, tl = Total.Length, fw = Fish.weight, 
               age_diss = Age_at_diss1, fsex = Sex, fdead = Dead, intensity = Inf.)
sdat$dose <- 2 # each fish got two worms

## ------------------------------------------------------------------------
sapply(sdat, function(x) sum(is.na(x))) # missing values in each variable

## ------------------------------------------------------------------------
sd_avg <- filter(sdat, !is.na(intensity))%>% #only select those with intensity data (one fish with missing data)
  group_by(trt)%>%
  summarize(n = n(), tdose = sum(dose, na.rm=T), tint = sum(intensity, na.rm=T))%>%
  mutate(inf.rate = tint/tdose)%>%
  select(trt, n, inf.rate)
sd_avg 

## ------------------------------------------------------------------------
mdat <- filter(sdat, !is.na(fw), !is.na(tl), !is.na(intensity))

global.mod <- glm(cbind(intensity, dose - intensity) ~ trt + (tl * fw) + fsex + worm.fam, 
                  data = mdat, family = 'binomial', na.action = "na.fail") # global model

## ------------------------------------------------------------------------
model.set <- dredge(global.mod, fixed = 'trt') # only include models with treatment in set
importance(model.set)

## ---- message=FALSE, warning=FALSE---------------------------------------
ggplot(data = sdat, aes(y = intensity/dose, x = tl, color = trt)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = 'lm') 

## ---- message=FALSE, warning=FALSE---------------------------------------
ggplot(data = sdat, aes(y = intensity/dose, x = fw, color = trt)) + 
  geom_point() + 
  geom_smooth(method = 'lm') 

## ---- warning=FALSE------------------------------------------------------
mdat <- filter(sdat, !is.na(tl), !is.na(intensity))

global.mod <- glm(cbind(intensity, dose - intensity) ~ trt + tl + fsex + worm.fam, 
                  data = mdat, family = 'binomial', na.action = "na.fail")
model.set <- dredge(global.mod, fixed = 'trt') # only include models with treatment in set
importance(model.set)

## ------------------------------------------------------------------------
model.set

## ------------------------------------------------------------------------
# make an age factor variable
mdat <- mutate(mdat, age_fac = if_else(age_diss < 100, "dead",
                                       if_else(age_diss < 110, "most fish", "late dissection")))

# plot it
ggplot(data = mdat, aes(y = intensity/dose, x = tl, color = age_fac)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = 'lm')

## ------------------------------------------------------------------------
# just fit to fish of same age
modb <- glm(cbind(intensity, dose - intensity) ~ tl + trt, 
            data = filter(mdat, age_fac == "most fish"),
            family = 'binomial')
summary(modb)

## ------------------------------------------------------------------------
filter(mdat, age_fac == "most fish")%>%
  group_by(trt)%>%
  summarize(n = n(), mean_length = mean(tl, na.rm=T), std_dev_length = sd(tl, na.rm=T))

## ------------------------------------------------------------------------
mod1 <- glm(cbind(intensity, dose - intensity) ~ trt, data = sdat, family = 'binomial')
summary(mod1)

## ------------------------------------------------------------------------
mod1 <- glm(cbind(intensity, dose - intensity) ~ trt, data = sdat, family = 'quasibinomial')
summary(mod1)
anova(mod1, test = 'F')

## ------------------------------------------------------------------------
mod2 <- glm(cbind(intensity, dose - intensity) ~ fsex + trt, data = sdat, family = 'quasibinomial')
anova(mod1, mod2, test = "LRT")

## ------------------------------------------------------------------------
mod3 <- glm(cbind(intensity, dose - intensity) ~ fsex + worm.fam + trt, data = sdat, family = 'quasibinomial')
anova(mod2, mod3, test = "LRT")

## ---- message=FALSE, warning=FALSE---------------------------------------
# get the expected binomial distribution, given the overall mean infection rate
mean.inf.rate <- mean(sdat$intensity/sdat$dose, na.rm=T)
expected <- dbinom(0:2, size = 2, prob = mean.inf.rate)

# calculate the proportion of fish with different infection intensities for each treatment
observed <- data.frame( prop.table(table(sdat$intensity)) )
names(observed) <- c('intensity', 'observed')

# combine expected and observed distributions and plot
binom_exp <- cbind(observed, expected)%>%
  gather('dist', 'freq', expected:observed)

ggplot(binom_exp, aes(x = intensity, y = freq, fill = dist)) + 
  geom_bar(stat = 'identity', position = position_dodge()) +
  labs(y = 'frequency', fill = NULL) +
  theme_bw()

## ------------------------------------------------------------------------
# make one-dimensional contingency table
cont.table <- table(sdat$intensity)

# calculate chi-square test; simulation used to get p-value, given small sample sizes
chisq.test(cont.table, p = expected, simulate.p.value = TRUE, B = 10000)

## ---- message=FALSE, warning=FALSE---------------------------------------
# get the expected binomial distribution, given the mean in each treatments
expected <- c(dbinom(0:2, size = 2, prob = sd_avg$inf.rate[1]),
              dbinom(0:2, size = 2, prob = sd_avg$inf.rate[2]))
expected <- data.frame(trt = rep(sd_avg$trt, each = 3), 
                       intensity = factor(rep(0:2, 2)),
                       expected)

# calculate the proportion of fish with different infection intensities for each treatment
observed <- data.frame( prop.table(table(sdat$trt, sdat$intensity), 1) )
names(observed) <- c('trt', 'intensity', 'observed')

# combine expected and observed distributions and plot
binom_exp <- left_join(expected, observed)%>%
  gather('dist', 'freq', expected:observed)

ggplot(binom_exp, aes(x = intensity, y = freq, fill = dist)) + 
  geom_bar(stat = 'identity', position = position_dodge()) +
  labs(y = 'frequency', fill = NULL) +
  facet_wrap(~trt) + theme_bw()

## ------------------------------------------------------------------------
# contingency table
cont.table <- table(sdat$trt, sdat$intensity)

# test if the observed and expected distributions differ for the uncrowded treatment
chisq.test(x = cont.table[1,],
           p = dbinom(0:2, size = 2, prob = sd_avg$inf.rate[1]))
# test if the observed and expected distributions differ for the crowded treatment
chisq.test(x = cont.table[2,],
           p = dbinom(0:2, size = 2, prob = sd_avg$inf.rate[2]))

# using the option to get a p-value by Monte Carlo simulation yields comparable results

## ------------------------------------------------------------------------
chisq.test(cont.table)

## ------------------------------------------------------------------------
# set an aesthetic theme for plots
theme.o<-theme_update(
  axis.text.y = element_text(colour="black", size = 12),
  axis.text.x = element_text(colour="black", size = 12, face = "bold"),
  axis.title.y = element_text(colour="black", size = 15, angle = 90, face = "bold", lineheight=0.4),
  axis.title.x = element_blank(),
  axis.ticks = element_line(colour="black"),
  panel.border = element_rect(colour = "black",fill=NA),
  panel.grid.minor=element_blank(),
  panel.grid.major.x=element_blank(),
  panel.grid.major=element_line(color="gray",linetype = "dotted"),
  panel.background= element_rect(fill = NA))

## ------------------------------------------------------------------------
# create function to extract means and 95% conf. intervals from glm results
get_ci_log_reg <- function(mod) { #takes a fitted logistic regression
  mod_sum <- summary(mod)
  
  # extract parameters and calculate 95% cis
  param <- mod_sum$coefficients[,1]
  se <- mod_sum$coefficients[,2]
  cil <- param - 1.96*se
  ciu <- param + 1.96*se
  
  sd_avg <- data.frame(param, cil, ciu) # make data frame
  sd_avg <- mutate_all(sd_avg, inv.logit) # convert from logistic to proportion
  sd_avg$trt <- names(param) # add parameter names
  
  return(sd_avg)
}

## ---- message=FALSE, warning=FALSE---------------------------------------
# plot just the average infection rates; no correction for fish length
mod1.1 <- glm(cbind(intensity, dose - intensity) ~ trt - 1, data = sdat, family = 'quasibinomial')

sd_avg <- get_ci_log_reg(mod1.1)
sd_avg <- mutate(sd_avg, trt = factor(trt, labels = c("One/Copepod", "Two/Copepod")))

sdat <- mutate(sdat, trtf = factor(trt, labels = c("One/Copepod", "Two/Copepod"))) # need this var for overlaying data on means

ggplot(sd_avg, aes(x = trt, y = param)) +
  geom_dotplot(data = sdat, aes(y = intensity/dose, x = trtf),
               fill = 'red', alpha = 0.5,
               binaxis = 'y', stackdir = 'center') +
  geom_point(size = 10, shape = "-") +
  geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.25) +
  labs(y = "Infection rate") +
  scale_y_continuous(limits = c(0,1)) + 
  scale_x_discrete(expand = c(0.25, 0))

## ---- message=FALSE, warning=FALSE---------------------------------------
# fit model
mod5 <- glm(cbind(intensity, dose - intensity) ~ trt-1 + tl, data = sdat, family = 'quasibinomial')

# make new dataset for predictions to plot
newdat <- data.frame(tl = rep(seq(min(sdat$tl, na.rm = T), max(sdat$tl, na.rm = T), 0.1), 2))
newdat$trt <- factor(rep(c("1+1", "2+0"), each = length(newdat$tl)/2 ))

# make predictions
preddat <- as.data.frame(predict(mod5, newdata = newdat, se = T))
preddat <- cbind(newdat, preddat)
preddat <- mutate(preddat,
                  ciu = fit + 1.96 * se.fit,
                  cil = fit - 1.96 * se.fit,
                  trt = factor(trt, labels = c("One/Copepod", "Two/Copepod")))

sdat$trtf <- factor(sdat$trt, labels = c("One/Copepod", "Two/Copepod"))


ggplot(sdat, aes(x = tl, y = intensity/dose, color = trtf)) + 
  geom_point(size = 3, alpha = 0.5, 
             position = position_jitter(width = 0.1, height = 0)) +
  geom_line(data = preddat, aes(x = tl, y = inv.logit(fit), color = trt), size = 1.5) +
  geom_ribbon(data = preddat, alpha = 0.5,
              aes(x = tl, 
                  y = inv.logit(fit), 
                  ymin = inv.logit(cil),
                  ymax = inv.logit(ciu),
                  fill = trt, color = NULL)) +
  labs(y = "Infection rate", x = "\nFish length (mm)",
       color = NULL, fill = NULL) +
  theme(axis.title.x = element_text(size = 15, face = "bold", lineheight=0.4),
        axis.text.x = element_text(colour="black", size = 12, face = "plain"))


