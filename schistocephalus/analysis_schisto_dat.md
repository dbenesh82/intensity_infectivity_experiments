Analysis of Schisto double infections
================

**Background**: these data are from an experiment in which stickleback fish were infected with a tapeworm. Fish were each given two tapeworm larvae via an intermediate host (a copepod). There were two treatments. Fish were either given two larvae in one copepod (the crowded treatment) or they were given two larvae in two copepods (one larvae per copepod, the uncrowded treatment). The goal of the experiment was to determine if crowding in the intermediate host affects a worm's chance of successfully infecting a fish.

**Analysis** Preliminary stuff. Set the wd, import libraries, and look at the data.

``` r
library(ggplot2)
library(dplyr)
library(tidyr)
library(boot)
library(MuMIn)

setwd("C:/Users/phosp/OneDrive/Documents/Benesh/Research/Intensity_infectivity/analysis/schistocephalus/")

sdat <- read.csv(file = "Schisto_double_exposure.csv", header = TRUE, sep = ',')
head(sdat)
```

    ##   Fish_no Clutch Treatment Total.Length Std.Length Fish.weight Inf. Sex
    ## 1       1      S       2+0           45         40       0.774    1   f
    ## 2       2      S       2+0           37         33       0.391    2   f
    ## 3       3      S       2+0           40         35       0.444    2   f
    ## 4       4      S       2+0           37         33       0.367    1   m
    ## 5       5      S       1+1           40         36       0.540    2   f
    ## 6       6      S       1+1           38         34       0.410    2   m
    ##   Worm1 Worm2 Dead Dissected Early.exposure Later.exposure Age_at_diss1
    ## 1 184.7    NA    0 10/5/2010      6/23/2010      6/28/2010          104
    ## 2  73.6  66.0    0 10/5/2010      6/23/2010      6/28/2010          104
    ## 3  62.2  64.8    0 10/5/2010      6/23/2010      6/28/2010          104
    ## 4 113.2    NA    0 10/5/2010      6/23/2010      6/28/2010          104
    ## 5  93.5  62.3    0 10/5/2010      6/23/2010      6/28/2010          104
    ## 6  73.7  58.5    0 10/5/2010      6/23/2010      6/28/2010          104
    ##   Age_at_diss2 Remark
    ## 1           99       
    ## 2           99       
    ## 3           99       
    ## 4           99       
    ## 5           99       
    ## 6           99

Clean up the data a bit.

``` r
sdat <- select(sdat, worm.fam = Clutch, trt = Treatment, tl = Total.Length, fw = Fish.weight, 
               age_diss = Age_at_diss1, fsex = Sex, fdead = Dead, intensity = Inf.)
sdat$dose <- 2 # each fish got two worms
```

Look at which variables have missing data.

``` r
sapply(sdat, function(x) sum(is.na(x))) # missing values in each variable
```

    ##  worm.fam       trt        tl        fw  age_diss      fsex     fdead 
    ##         0         0         1        10         1         0         1 
    ## intensity      dose 
    ##         1         0

One exposed fish was lost and thus has empty values for 'fdead' and 'intensity'. It can be excluded. Beware of using fish weight in models, as it is missing for 10 fish. These fish died before dissections, but months after exposure, so their worms would not have been too small to be overlooked.

Calculate infection rates for the two treatments.

``` r
sd_avg <- filter(sdat, !is.na(intensity))%>% #only select those with intensity data (one fish with missing data)
  group_by(trt)%>%
  summarize(n = n(), tdose = sum(dose, na.rm=T), tint = sum(intensity, na.rm=T))%>%
  mutate(inf.rate = tint/tdose)%>%
  select(trt, n, inf.rate)
sd_avg 
```

    ## # A tibble: 2 × 3
    ##      trt     n  inf.rate
    ##   <fctr> <int>     <dbl>
    ## 1    1+1    27 0.3703704
    ## 2    2+0    30 0.3500000

Does not look like there is a difference, 37% successfully infected in the uncrowded treatment while 35% successfully infected in the crowded treatment.

Let's start fitting models. The response variable (infection) is a proportion, so we will run a logistic regression. To get an idea of which variables are important, we'll fit many models. All the models will include treatment, because that's what we want to test. But first we need to filter rows with missing values, so that the sample size is constant, regardless of the variables in the model.

``` r
mdat <- filter(sdat, !is.na(fw), !is.na(tl), !is.na(intensity))

global.mod <- glm(cbind(intensity, dose - intensity) ~ trt + (tl * fw) + fsex + worm.fam, 
                  data = mdat, family = 'binomial', na.action = "na.fail") # global model
```

Some fish variables like length and weight might be important, but fish sex appears unimportant. Worm family also does not seem important.

``` r
model.set <- dredge(global.mod, fixed = 'trt') # only include models with treatment in set
```

    ## Fixed terms are "trt" and "(Intercept)"

``` r
importance(model.set)
```

    ##                      trt  tl   fw   fw:tl worm.fam fsex
    ## Importance:          1.00 1.00 1.00 0.47  0.25     0.23
    ## N containing models:   20   12   12    4    10       10

Here's the plot of fish length vs infectivity. Larger fish had lower infection rates. This is either because small fish are more susceptible or because uninfected fish grow longer.

``` r
ggplot(data = sdat, aes(y = intensity/dose, x = tl, color = trt)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = 'lm') 
```

![](analysis_schisto_dat_files/figure-markdown_github/unnamed-chunk-7-1.png)

Here's the same plot, but with fish weight.

``` r
ggplot(data = sdat, aes(y = intensity/dose, x = fw, color = trt)) + 
  geom_point() + 
  geom_smooth(method = 'lm') 
```

![](analysis_schisto_dat_files/figure-markdown_github/unnamed-chunk-8-1.png)

It increases in heavier fish. This is probably caused by the fact that the fish weight variable also includes worm weight. Let's exclude fish weight, given that this pattern is not convincing and that there was missing data for this variable. The model still includes fish length, which is highly correlated with weight.

``` r
mdat <- filter(sdat, !is.na(tl), !is.na(intensity))

global.mod <- glm(cbind(intensity, dose - intensity) ~ trt + tl + fsex + worm.fam, 
                  data = mdat, family = 'binomial', na.action = "na.fail")
model.set <- dredge(global.mod, fixed = 'trt') # only include models with treatment in set
```

    ## Fixed terms are "trt" and "(Intercept)"

``` r
importance(model.set)
```

    ##                      trt  tl   fsex worm.fam
    ## Importance:          1.00 1.00 0.26 0.24    
    ## N containing models:    8    4    4    4

Now total fish length is the only variable that seems to be important. When we look at the top models, we see that the first three are by far the best (cumulative weight ~ 0.93). They all include fish length. The top model only includes fish length.

``` r
model.set
```

    ## Global model call: glm(formula = cbind(intensity, dose - intensity) ~ trt + tl + 
    ##     fsex + worm.fam, family = "binomial", data = mdat, na.action = "na.fail")
    ## ---
    ## Model selection table 
    ##     (Int) fsx      tl trt wrm.fam df  logLik  AICc delta weight
    ## 3 10.1100     -0.2625   +          3 -53.337 113.1  0.00  0.559
    ## 4 10.7500   + -0.2754   +          4 -53.208 115.2  2.06  0.200
    ## 7 10.1500     -0.2631   +       +  4 -53.335 115.4  2.31  0.176
    ## 8 10.8300   + -0.2767   +       +  5 -53.202 117.6  4.45  0.060
    ## 1 -0.5306               +          2 -59.885 124.0 10.87  0.002
    ## 2 -0.6810   +           +          3 -59.583 125.6 12.49  0.001
    ## 5 -0.6179               +       +  3 -59.780 126.0 12.89  0.001
    ## 6 -0.7809   +           +       +  4 -59.460 127.7 14.56  0.000
    ## Models ranked by AICc(x)

Interpreting this fish length effect is complicated though. Small fish may be more susceptible to *Schistocephalus* infection. Some of my previous studies found negative, but non-significant associations between fish weight at exposure and infectivity (see [here](http://onlinelibrary.wiley.com/doi/10.1111/j.1600-0706.2011.19925.x/full) and [here](http://onlinelibrary.wiley.com/doi/10.1111/evo.12388/full)). However, in this experiment the fish were measured at dissection and not at exposure, so it may simply be that uninfected fish grew longer than infected fish, causing the negative relationship. Also, fish were not dissected at the same time. Most fish were dissected around 100 days after exposure, but some died earlier and a few were dissected after ~170 days to have large worms for breeding.

Thus, fish length is a mixture of fish size at exposure, fish growth rate, and fish growth duration. We can at least look at the effect of fish size at a given age at dissection. We can define three age groups: 1) fish dissected about 3 months post exposure (the majority), 2) fish that died before the planned dissection date, and 3) a handful of fish that were dissected after about 5.5 months. When we plot infection rate vs length for these three groups, there is still a negative relationship for the main group ('most fish').

``` r
# make an age factor variable
mdat <- mutate(mdat, age_fac = if_else(age_diss < 100, "dead",
                                       if_else(age_diss < 110, "most fish", "late dissection")))

# plot it
ggplot(data = mdat, aes(y = intensity/dose, x = tl, color = age_fac)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = 'lm')
```

![](analysis_schisto_dat_files/figure-markdown_github/unnamed-chunk-11-1.png)

And this negative relationship is significant when we fit a logistic regression with only the majority of fish dissected at about the same time. Treatment effects are still non-significant, though.

``` r
# just fit to fish of same age
modb <- glm(cbind(intensity, dose - intensity) ~ tl + trt, 
            data = filter(mdat, age_fac == "most fish"),
            family = 'binomial')
summary(modb)
```

    ## 
    ## Call:
    ## glm(formula = cbind(intensity, dose - intensity) ~ tl + trt, 
    ##     family = "binomial", data = filter(mdat, age_fac == "most fish"))
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.6993  -0.9800  -0.6064   0.4934   2.7588  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  11.2817     4.2908   2.629  0.00856 **
    ## tl           -0.2997     0.1064  -2.817  0.00484 **
    ## trt2+0       -0.4370     0.5142  -0.850  0.39534   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 69.857  on 43  degrees of freedom
    ## Residual deviance: 60.328  on 41  degrees of freedom
    ## AIC: 82.963
    ## 
    ## Number of Fisher Scoring iterations: 4

Controlling for fish length did not change the treatment differences, which is not surprising given that average fish length was very comparable in the two treatment groups.

``` r
filter(mdat, age_fac == "most fish")%>%
  group_by(trt)%>%
  summarize(n = n(), mean_length = mean(tl, na.rm=T), std_dev_length = sd(tl, na.rm=T))
```

    ## # A tibble: 2 × 4
    ##      trt     n mean_length std_dev_length
    ##   <fctr> <int>       <dbl>          <dbl>
    ## 1    1+1    20    40.80000       2.483631
    ## 2    2+0    24    40.41667       2.888006

So fish size variables are difficult to interpret and have no impact on our main interest, treatment effects. Thus, we do not consider further models including fish size.

Another approach to model-building is to think of appropriate models a priori and then compare them. Let's start with the simplest model, with only treatment.

*Model 1*: Just treatment

``` r
mod1 <- glm(cbind(intensity, dose - intensity) ~ trt, data = sdat, family = 'binomial')
summary(mod1)
```

    ## 
    ## Call:
    ## glm(formula = cbind(intensity, dose - intensity) ~ trt, family = "binomial", 
    ##     data = sdat)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.3603  -1.3127   0.3730   0.4343   2.0492  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept) -0.53063    0.28180  -1.883   0.0597 .
    ## trt2+0      -0.08841    0.39073  -0.226   0.8210  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 90.709  on 56  degrees of freedom
    ## Residual deviance: 90.658  on 55  degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## AIC: 123.77
    ## 
    ## Number of Fisher Scoring iterations: 4

There is no effect of crowding. Overdispersion seems mild. The ratio of residual deviance to residual degrees of freedom is about 1.6, so we would expect the results to be unaffected by using the argument `family = 'quasibinomial'`.

``` r
mod1 <- glm(cbind(intensity, dose - intensity) ~ trt, data = sdat, family = 'quasibinomial')
summary(mod1)
```

    ## 
    ## Call:
    ## glm(formula = cbind(intensity, dose - intensity) ~ trt, family = "quasibinomial", 
    ##     data = sdat)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.3603  -1.3127   0.3730   0.4343   2.0492  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept) -0.53063    0.31471  -1.686   0.0974 .
    ## trt2+0      -0.08841    0.43637  -0.203   0.8402  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasibinomial family taken to be 1.247232)
    ## 
    ##     Null deviance: 90.709  on 56  degrees of freedom
    ## Residual deviance: 90.658  on 55  degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

The results are similar. Below we examine the distribution of the data and whether the binomial distribution fits the data well. But first we'll fit a couple more models. To be conservative we will use `family = 'quasibinomial'`

*Model 2*: treatment plus fish variables (just sex, given concerns with fish length discussed above)

``` r
mod2 <- glm(cbind(intensity, dose - intensity) ~ fsex + trt, data = sdat, family = 'quasibinomial')
anova(mod1, mod2, test = "LRT")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: cbind(intensity, dose - intensity) ~ trt
    ## Model 2: cbind(intensity, dose - intensity) ~ fsex + trt
    ##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
    ## 1        55     90.658                     
    ## 2        54     90.054  1  0.60332   0.4912

This is not a significantly better model. Finally, we'll add one more term, worm family (there were two used in the experiment).

*Model 3*: treatment, fish sex, worm family

``` r
mod3 <- glm(cbind(intensity, dose - intensity) ~ fsex + worm.fam + trt, data = sdat, family = 'quasibinomial')
anova(mod2, mod3, test = "LRT")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: cbind(intensity, dose - intensity) ~ fsex + trt
    ## Model 2: cbind(intensity, dose - intensity) ~ fsex + worm.fam + trt
    ##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
    ## 1        54     90.054                     
    ## 2        53     89.808  1   0.2461    0.664

This model is not an improvement, and treatment is still not significant (not shown). Thus, we can conclude that the treatment effect (or lack thereof) is not sensitive to model specifications.

So I think two models can be justified. I prefer a simple model with just `trt`, but a larger model including several control variables (`trt + fsex + worm.fam`) is also reasonable. Neither suggests that crowding affected establishment probability.

We noted some overdispersion in the logistic regression, so let's look at how well the binomial approximates the observed data distribution.

``` r
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
```

![](analysis_schisto_dat_files/figure-markdown_github/unnamed-chunk-18-1.png)

For single infections (uncrowded treatment) it fits quite well. But for fish exposed to doubly-infected copepods, we see more null and double infections than expected. This suggests that when worms share a copepod host, their chance to infect is not independent. That is, both infect or both fail to infect more frequently than expected. We can do chi-square tests to see if these differences are significant. First, we'll confirm that the uncrowded treatment is well-approximated by the binomial.

``` r
# contingency table
cont.table <- table(sdat$trt, sdat$intensity)

# test if the observed and expected distributions differ for the uncrowded treatment
chisq.test(x = cont.table[1,],
           p = dbinom(0:2, size = 2, prob = sd_avg$inf.rate[1]))
```

    ## Warning in chisq.test(x = cont.table[1, ], p = dbinom(0:2, size = 2, prob =
    ## sd_avg$inf.rate[1])): Chi-squared approximation may be incorrect

    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  cont.table[1, ]
    ## X-squared = 0.33727, df = 2, p-value = 0.8448

``` r
# using the option to get a p-value by Monte Carlo simulation yields comparable results
```

The difference is non-significant. Now, we'll check in the crowded treatment. It is marginally significant, suggesting that null and double infections are overrepresented, relative to expectations.

``` r
# test if the observed and expected distributions differ for the crowded treatment
chisq.test(x = cont.table[2,],
           p = dbinom(0:2, size = 2, prob = sd_avg$inf.rate[2]))
```

    ## Warning in chisq.test(x = cont.table[2, ], p = dbinom(0:2, size = 2, prob =
    ## sd_avg$inf.rate[2])): Chi-squared approximation may be incorrect

    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  cont.table[2, ]
    ## X-squared = 7.1203, df = 2, p-value = 0.02843

``` r
# using the option to get a p-value by Monte Carlo simulation yields comparable but slightly more conservative results
```

Instead of comparing the distribution of each treatment to the expected binomial distribution, we could compare them to each other. That is are the data distributed differently in the two treatments, even though the means are not different?

``` r
chisq.test(cont.table)
```

    ## Warning in chisq.test(cont.table): Chi-squared approximation may be
    ## incorrect

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  cont.table
    ## X-squared = 5.1744, df = 2, p-value = 0.07523

This chi-square test is marginally non-significant. So, there is suggestive evidence that the data are distributed differently in the two treatments, despite their similar mean values.

Let's make a couple nice plots that visualize the data.

``` r
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
```

``` r
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
```

First make a plot comparing treatment groups.

``` r
# plot just the average infection rates; no correction for fish length
mod1.1 <- glm(cbind(intensity, dose - intensity) ~ trt - 1, data = sdat, family = 'quasibinomial')

sd_avg <- get_ci_log_reg(mod1.1)
sd_avg <- mutate(sd_avg, trt = factor(trt, labels = c("Uncrowded", "Crowded")))

sdat <- mutate(sdat, trtf = factor(trt, labels = c("Uncrowded", "Crowded"))) # need this var for overlaying data on means

ggplot(sd_avg, aes(x = trt, y = param)) +
  geom_dotplot(data = sdat, aes(y = intensity/dose, x = trtf),
               fill = 'red', alpha = 0.5,
               binaxis = 'y', stackdir = 'center') +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.25) +
  labs(y = "Infection rate") +
  scale_y_continuous(limits = c(0,1)) + 
  scale_x_discrete(expand = c(0.25, 0))
```

![](analysis_schisto_dat_files/figure-markdown_github/unnamed-chunk-24-1.png)

The means (black points) are very similar in the two groups, but the red points show how the distribution of the data might be different in the two treatments.

Given that fish length appeared important, if not easily interpretable, we can plot the effect of both treatment and fish length simultaneously.

``` r
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
                  trt = factor(trt, labels = c("Uncrowded", "Crowded")))

sdat$trtf <- factor(sdat$trt, labels = c("Uncrowded", "Crowded"))


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
```

![](analysis_schisto_dat_files/figure-markdown_github/unnamed-chunk-25-1.png)

**Conclusion**: There is little support for the hypothesis that crowding reduces the infection rate of these worms. However, it is possible that worms from crowded intermediate hosts are more likely to either both establish or neither establish. This needs further study.
