## ---- message=FALSE, warning=FALSE---------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(boot)
library(MuMIn)

cdat <- read.csv(file = "cam_fish_inf.csv", header = TRUE, sep = ',')
head(cdat)

## ------------------------------------------------------------------------
cdat <- select(cdat, block = Block, trt = Treatment, tl = Length, fw = Weight,
               ffam = Fish_family, fsex = Fish_sex, 
               uneaten = Uneaten, dose = Dose, intensity = Intensity)%>%
  mutate(block = factor(block))
# note: there is one fish with an uneaten copepod. This copepod was uninfected, so the dose should have still been 6.

## ------------------------------------------------------------------------
sapply(cdat, function(x) sum(is.na(x))) # missing values in each variable

## ------------------------------------------------------------------------
cd_avg <- group_by(cdat, trt)%>%
  summarize(n = n(), tdose = sum(dose, na.rm=T), tint = sum(intensity, na.rm=T))%>%
  mutate(inf.rate = tint/tdose)%>%
  select(trt, n, inf.rate)
cd_avg 

## ------------------------------------------------------------------------
mdat <- filter(cdat, !is.na(fsex)) # remove rows with missing fish sex data; sex was only variable with missing data

global.mod <- glm(cbind(intensity, dose - intensity) ~ trt + block + (tl * fw) + fsex + ffam, 
                  data = mdat, family = 'binomial', na.action = "na.fail") # global model

## ------------------------------------------------------------------------
model.set <- dredge(global.mod, fixed = 'trt') # only include models with treatment in set
importance(model.set)

## ------------------------------------------------------------------------
group_by(cdat, fsex)%>%
  summarize(n = n(), tdose = sum(dose, na.rm=T), tint = sum(intensity, na.rm=T))%>%
  mutate(inf.rate = tint/tdose)%>%
  select(fsex, n, inf.rate)

## ------------------------------------------------------------------------
# use full data and re-fit global model
global.mod <- glm(cbind(intensity, dose - intensity) ~ trt + block + (tl * fw) + ffam, 
                  data = cdat, family = 'binomial', na.action = "na.fail")
model.set <- dredge(global.mod, fixed = 'trt') # only include models with treatment in set
importance(model.set)

## ------------------------------------------------------------------------
ggplot(data = cdat, aes(x = tl, y = intensity/dose, color = is.na(fsex))) + 
  geom_point() + 
  labs(color = "No sex?", y = "infection rate", x = "total length") +
  geom_smooth(se=F, method='lm') + # for subgroups
  geom_smooth(aes(color = NULL), se = FALSE, color = 'black', method = 'lm') # for full dataset

## ---- message=FALSE, warning=FALSE---------------------------------------
ggplot(cdat, aes(y = intensity/dose, x = block)) + geom_boxplot() + 
  geom_dotplot(fill = 'red', binaxis = 'y', stackdir = 'center')

## ------------------------------------------------------------------------
model.set[1:5,]

## ------------------------------------------------------------------------
mod0 <- glm(cbind(intensity, dose - intensity) ~ block + tl + trt, data = cdat, family = 'binomial')
summary(mod0)

## ---- message=FALSE, warning=FALSE---------------------------------------
# get expected frequencies, based on binomial, given mean infection rate
mean.inf.rate <- mean(cdat$intensity/cdat$dose, na.rm=T)
expected <- dbinom(0:6, size = 6, prob = mean.inf.rate)

# calculate the proportion of fish with different infection intensities for each treatment
observed <- data.frame( prop.table(table(cdat$intensity)) )
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
cont.table <- table(cdat$intensity)

# calculate chi-square test; simulation used to get p-value, given small sample sizes
chisq.test(cont.table, p = expected, simulate.p.value = TRUE, B = 10000)

## ---- message=FALSE, warning=FALSE---------------------------------------
# get the expected binomial distribution, given the mean in each treatments
expected <- c(dbinom(0:6, size = 6, prob = cd_avg$inf.rate[1]),
              dbinom(0:6, size = 6, prob = cd_avg$inf.rate[2]),
              dbinom(0:6, size = 6, prob = cd_avg$inf.rate[3]))
expected <- data.frame(trt = rep(cd_avg$trt, each = 7), 
                       intensity = factor(rep(0:6, 3)),
                       expected)

# calculate the proportion of fish with different infection intensities for each treatment
observed <- data.frame( prop.table(table(cdat$trt, cdat$intensity), 1) )
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
cont.table <- table(cdat$trt, cdat$intensity)

# test if the observed and expected distributions differ for the one per copepod treatment
chisq.test(x = cont.table[1,],
           p = dbinom(0:6, size = 6, prob = cd_avg$inf.rate[1]),
           simulate.p.value = TRUE, B = 10000)
# test if the observed and expected distributions differ for the two per copepod treatment
chisq.test(x = cont.table[2,],
           p = dbinom(0:6, size = 6, prob = cd_avg$inf.rate[2]),
           simulate.p.value = TRUE, B = 10000)
# test if the observed and expected distributions differ for the three per copepod treatment
chisq.test(x = cont.table[3,],
           p = dbinom(0:6, size = 6, prob = cd_avg$inf.rate[3]),
           simulate.p.value = TRUE, B = 10000)
# using the option to get a p-value by Monte Carlo simulation is needed given small sample sizes

## ------------------------------------------------------------------------
chisq.test(cont.table, simulate.p.value = TRUE, B = 10000)

## ------------------------------------------------------------------------
mod0 <- glm(cbind(intensity, dose - intensity) ~ block + tl + trt, data = cdat, family = 'quasibinomial')
anova(mod0, test = "F")

## ------------------------------------------------------------------------
mod1 <- glm(cbind(intensity, dose - intensity) ~ trt, data = cdat, family = 'quasibinomial')
anova(mod1, test = "F")

## ------------------------------------------------------------------------
mod2 <- glm(cbind(intensity, dose - intensity) ~ block + trt, data = cdat, family = 'quasibinomial')
anova(mod2, test = "F")

## ------------------------------------------------------------------------
mod3 <- glm(cbind(intensity, dose - intensity) ~ tl + ffam + trt, data = cdat, family = 'quasibinomial')
anova(mod3, test = "F")

## ------------------------------------------------------------------------
mod4 <- glm(cbind(intensity, dose - intensity) ~ block + tl + ffam + trt, data = cdat, family = 'quasibinomial')
anova(mod4, test = "F")

## ------------------------------------------------------------------------
anova(mod1, mod4, test = "F") 

## ------------------------------------------------------------------------
modp <- glm(cbind(intensity, dose - intensity) ~ trt - 1 + block, data = cdat, family = 'quasibinomial')

# make new dataset for predictions to plot
newdat <- data.frame(trt = rep(levels(cdat$trt), each = 2),
                     block = rep(c("1", "2"), 3))
newdat <- filter(newdat, trt != "1x6" | block != 2) # remove point where no data exist

# make prediction
preddat <- as.data.frame(predict(modp, newdata = newdat, se = T))
preddat <- cbind(newdat, preddat)
preddat <- mutate(preddat,
                  ciu = fit + 1.96 * se.fit,
                  cil = fit - 1.96 * se.fit,
                  trt = factor(trt, labels = c("One/Copepod", "Two/Copepod", "Three/Copepod")))


# need this to overlay data points over means
cdat$trtf <- factor(cdat$trt, labels = c("One/Copepod", "Two/Copepod", "Three/Copepod"))

## ---- message=FALSE, warning=FALSE---------------------------------------
#set a aesthetic theme for plots
theme.o<-theme_update(
  axis.text.y = element_text(colour="black", size = 12),
  axis.text.x = element_text(colour="black", size = 10, face = "bold"),
  axis.title.y = element_text(colour="black", size = 15, angle = 90, face = "bold", lineheight=0.4),
  axis.title.x = element_blank(),
  axis.ticks = element_line(colour="black"),
  panel.border = element_rect(colour = "black",fill=NA),
  panel.grid.minor=element_blank(),
  panel.grid.major.x=element_blank(),
  panel.grid.major=element_line(color="gray",linetype = "dotted"),
  panel.background= element_rect(fill = NA))


ggplot(preddat, aes(x = trt, y = inv.logit(fit))) + 
  geom_dotplot(data = cdat, aes(y = intensity/dose, x = trtf), 
               fill = 'red', alpha = 0.5,
               binaxis = 'y', stackdir = 'center') +
  geom_point(size = 10, shape = "-") +
  geom_errorbar(aes(ymin = inv.logit(cil), ymax = inv.logit(ciu)), width = 0.25) +
  labs(y = "Infection rate") +
  facet_grid(~block)

## ---- message=FALSE, warning=FALSE---------------------------------------
modp <- glm(cbind(intensity, dose - intensity) ~ trt - 1, data = cdat, family = 'quasibinomial')
# make new dataset for predictions to plot
newdat <- data.frame(trt = levels(cdat$trt))
# make prediction
preddat <- as.data.frame(predict(modp, newdata = newdat, se = T))
preddat <- cbind(newdat, preddat)
preddat <- mutate(preddat,
                  ciu = fit + 1.96 * se.fit,
                  cil = fit - 1.96 * se.fit,
                  trt = factor(trt, labels = c("One/Copepod", "Two/Copepod", "Three/Copepod")))

# make a plot
ggplot(preddat, aes(x = trt, y = inv.logit(fit))) + 
  geom_dotplot(data = cdat, aes(y = intensity/dose, x = trtf), 
               fill = 'red', alpha = 0.5,
               binaxis = 'y', stackdir = 'center') +
  geom_point(size = 10, shape = "-") +
  geom_errorbar(aes(ymin = inv.logit(cil), ymax = inv.logit(ciu)), width = 0.25) +
  labs(y = "Infection rate") 

