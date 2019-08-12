# Marine subsidies drive patterns of avian island biogeography: 

# Read in data, load packages:
set.seed(123)

library(dplyr)
library(lme4)
library(AICcmodavg)
library(MuMIn)
library(visreg)
library(DHARMa)

dat <- read.csv("birds-ibt-subsidies.csv", stringsAsFactors = TRUE)

# Check if habitat heterogeneity and isolation (at 90% of asymptote - see SI) contribute more than just area alone?
area.mod <- lmer(l.S.chao1 ~ l.area.std + 
                   year.std + 
                   (1|node),
                 data = dat, 
                 REML = FALSE) 

hh.mod <- lmer(l.S.chao1 ~ l.area.std + 
                 hab.het.std + 
                 year.std + 
                 (1|node),
               data = dat, 
               REML = FALSE) 

ibt.mod <- lmer(l.S.chao1 ~ l.area.std + 
                  l.isolation.90.std +
                  year.std + 
                  (1|node),
                data = dat, 
                REML = FALSE) 

full.mod <- lmer(l.S.chao1 ~ l.area.std + 
                   l.isolation.90.std +
                   hab.het.std + 
                  year.std + 
                   (1|node),
                 data = dat, 
                 REML = FALSE) 

null.mod <-lmer(l.S.chao1 ~ 1 + 
                  year.std + 
                  (1|node),
                data = dat, 
                REML = FALSE) 

table.rich <- aictab(c(area.mod, 
                  hh.mod, 
                  ibt.mod,
                  full.mod,
                  null.mod), 
                modnames = c("area",
                             "area + hh", 
                             "area + iso", 
                             "area + iso + hh",
                             "null")); table.rich

# Density: 
dens.area.mod <- glmer.nb(n.birds.50 ~ l.area.std + 
                            year.std + 
                            (1|node) + 
                            offset(log(buf50_area)),
                          data = dat,
                          na.action = "na.fail", 
                          control = glmerControl(optimizer = "bobyqa", 
                                                 optCtrl = list(maxfun = 2e6)))

dens.hh.mod <- glmer.nb(n.birds.50 ~ l.area.std + hab.het.std + 
                          year.std + 
                          (1|node) + 
                          offset(log(buf50_area)),
                        data = dat,
                        na.action = "na.fail", 
                        control = glmerControl(optimizer = "bobyqa", 
                                               optCtrl = list(maxfun = 2e6)))

dens.ibt.mod <- glmer.nb(n.birds.50 ~ l.area.std + l.isolation.90.std + 
                           year.std + 
                           (1|node) + 
                           offset(log(buf50_area)),
                         data = dat,
                         na.action = "na.fail", 
                         control = glmerControl(optimizer = "bobyqa", 
                                                optCtrl = list(maxfun = 2e6)))

dens.full.mod <- glmer.nb(n.birds.50 ~ l.area.std + 
                            l.isolation.90.std + 
                            hab.het.std +
                            year.std + 
                            (1|node) + 
                            offset(log(buf50_area)),
                          data = dat,
                          na.action = "na.fail", 
                          control = glmerControl(optimizer = "bobyqa", 
                                                 optCtrl = list(maxfun = 2e6)))

dens.null.mod <- glmer.nb(n.birds.50 ~ 1 + 
                            year.std +
                            (1|node) +
                            offset(log(buf50_area)),
                          data = dat,
                          na.action = "na.fail",
                          control = glmerControl(optimizer = "bobyqa",
                                                 optCtrl = list(maxfun = 2e6)))

table.dens <- aictab(c(dens.area.mod, 
                                dens.hh.mod, 
                                dens.ibt.mod,
                                dens.full.mod,
                                dens.null.mod), 
                              modnames = c("area",
                                           "area + hh", 
                                           "area + iso", 
                                           "area + iso + hh",
                                           "null")); table.dens

#### Richness ~ subsidies ####
rich.global <- lmer(l.S.chao1 ~ l.area.std + 
                      soild15N.std + 
                      sqrt.wrack.std + 
                      rocky.prop.std + 
                      l.area.std:soild15N.std + 
                      l.area.std:sqrt.wrack.std + 
                      sqrt.wrack.std:rocky.prop.std + 
                      soild15N.std:rocky.prop.std +
                      year.std + 
                      (1|node),
                    data = dat,
                    REML = FALSE,
                    na.action = "na.fail")

# Check VIFs: I dropped l.area.std:rocky.prop.std because Zuur, Ieno & elphick 2010 recommend dropping variables with VIF over 3. Once this is removed, none are over 3. 
car::vif(rich.global) 

# Analyze all possible combinations of remaining variables to calculate meaningful RVIs (Arnold et al. 2010, Burnham & Anderson 2001). 
# Keep year parameter in all possible combinations. 
rich.allsubsets <- dredge(rich.global, subset = ~ "year.std")
avg.rich <- model.avg(rich.allsubsets, fit = TRUE) 

summary(avg.rich) # Equal proportions between fixefs and between interactions.

# Model diagnostics on global model:
hist(summary(rich.global)$residuals) # residuals normally distributed
plot(rich.global) # not bad 
qqnorm(resid(rich.global))
qqline(resid(rich.global))

# Predictions with averaged richness model:
# Average island: 
avg.isl.rich <- expand.grid(l.area.std = mean(dat$l.area.std),
                        soild15N.std = mean(dat$soild15N.std),
                        rocky.prop.std = mean(dat$rocky.prop.std),
                        sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                        year.std = mean(dat$year.std))

pred.avg.isl.rich <- predict(avg.rich, avg.isl.rich, type = "response", se.fit = TRUE, re.form = NA, full = T)

10^pred.avg.isl.rich$fit # 8.23 species for an average island.
10^(pred.avg.isl.rich$se.fit *1.96) # +/- 95% CI: 1.13

# Predictions for island 1 SD bigger:
big.isl.rich <- expand.grid(l.area.std = mean(dat$l.area.std) + sd(dat$l.area.std), 
                        soild15N.std = mean(dat$soild15N.std),
                        rocky.prop.std = mean(dat$rocky.prop.std),
                        sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                        year.std = mean(dat$year.std))

pred.big.isl.rich <- predict(avg.rich, big.isl.rich, type = "response", se.fit = TRUE, re.form = NA, full = T)

10^pred.big.isl.rich$fit # 13.26938 species for an island bigger by 1 sd. 
10^(pred.big.isl.rich$se.fit*1.96) # +/- 1.1913 
10^pred.big.isl.rich$fit - 10^pred.avg.isl.rich$fit # 5.035039 more birds for one sd larger than average area. 

# Predictions for an island 1 SD smaller: 
small.isl.rich <- expand.grid(l.area.std = mean(dat$l.area.std) - sd(dat$l.area.std), 
                          soild15N.std = mean(dat$soild15N.std),
                          rocky.prop.std = mean(dat$rocky.prop.std),
                          sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                          year.std = mean(dat$year.std))

pred.small.isl.rich <- predict(avg.rich, small.isl.rich, type = "response", se.fit = TRUE, re.form = NA, full = T)

10^pred.small.isl.rich$fit # 5.10984 species for an island smaller than average by 1 sd. 
10^(pred.small.isl.rich$se.fit*1.96) # +/- 1.18477

10^pred.big.isl.rich$fit - 10^pred.small.isl.rich$fit # 8.15

# Calculate richness based on soil d15N: 
# Lowest d15N value: 1.3, highest = 14.7. Standardize these for comparison: 
(1.3 - mean(dat$soild15N)) / sd(dat$soild15N) # -1.844
(14.7 + mean(dat$soild15N)) / sd(dat$soild15N) # 7.192

low.d15n.rich <- expand.grid(l.area.std = mean(dat$l.area.std), 
                     soild15N.std = -1.844,
                     rocky.prop.std = mean(dat$rocky.prop.std),
                     sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                     year.std = mean(dat$year.std))

pred.low.d15n.rich <- predict(avg.rich, low.d15n.rich, type = "response", se.fit = TRUE, re.form = NA, full = T)

10^pred.low.d15n.rich$fit # 10.259

high.d15n.rich <- expand.grid(l.area.std = mean(dat$l.area.std), # This is 1.
                     soild15N.std = 7.192,
                     rocky.prop.std = mean(dat$rocky.prop.std),
                     sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                     year.std = mean(dat$year.std))

pred.high.d15n.rich <- predict(avg.rich, high.d15n.rich, type = "response", se.fit = TRUE, re.form = NA, full = T)

10^pred.low.d15n.rich$fit - 10^pred.high.d15n.rich$fit # 6.76 more species on low d15n islands than high d15n islands of equal (average) size.

#### Density ~ subsidies: ####
# Include same interactions as with richness model: 
density.global <- glmer.nb(n.birds.50 ~ l.area.std + 
                             soild15N.std + 
                             sqrt.wrack.std + 
                             rocky.prop.std +
                             l.area.std:soild15N.std + 
                             l.area.std:sqrt.wrack.std +
                             sqrt.wrack.std:rocky.prop.std + 
                             soild15N.std:rocky.prop.std +
                             year.std + 
                             (1|node) + 
                             offset(log(buf50_area)),
                           data = dat,
                           na.action = "na.fail", 
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))

# Check VIFs: none are over 2.
car::vif(density.global)

# Check global model: 
hist(summary(density.global)$residuals) # residuals normally distributed
plot(density.global) # Weird - check DHARma diagnostics (below) instead since glmm. 
qqnorm(resid(density.global))
qqline(resid(density.global)) # QQplot looks good.

# Test all subsets for meaningful RVIs: include the both year and the offset term in all models.
density.mod <- dredge(density.global, subset = ~ "year.std" & "offset(log(buf50_area))")
avg.density <- model.avg(density.mod, fit = TRUE) 

summary(avg.density)

# GLMM Model diagnostics: use DHARMa package.

simulationOutput.dens <- simulateResiduals(fittedModel = density.global, n = 250)

# Residual checks on these simulated residuals:
plotSimulatedResiduals(simulationOutput = simulationOutput.dens)
qqnorm(resid(density.global))
qqline(resid(density.global))

# Goodness of fit tests:
testUniformity(simulationOutput = simulationOutput.dens)

# test for zero inflation in data
testZeroInflation(simulationOutput.dens)

# # Check individual predictors:
plotResiduals(dat$l.area.std, simulationOutput.dens$scaledResiduals)
plotResiduals(dat$soild15N.std, simulationOutput.dens$scaledResiduals)
plotResiduals(dat$sqrt.wrack.std, simulationOutput.dens$scaledResiduals)
plotResiduals(dat$rocky.prop.std, simulationOutput.dens$scaledResiduals)
plotResiduals(dat$year.std, simulationOutput.dens$scaledResiduals)

ggplot(data.frame(l.area.std=dat$l.area.std,pearson=residuals(density.global, type="pearson")),
       aes(x=l.area.std,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(soild15N.std=dat$soild15N.std,pearson=residuals(density.global, type="pearson")),
       aes(x=soild15N.std,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(sqrt.wrack.std=dat$sqrt.wrack.std,pearson=residuals(density.global, type="pearson")),
       aes(x=sqrt.wrack.std,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(rocky.prop.std=dat$rocky.prop.std,pearson=residuals(density.global, type="pearson")),
       aes(x=rocky.prop.std,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(year.std=dat$year.std,pearson=residuals(density.global, type="pearson")),
       aes(x=year.std,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(year.std=dat$year.std, pearson=residuals(density.global, type="pearson")),
       aes(x=year.std,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(buf50_area=dat$buf50_area, pearson=residuals(density.global, type="pearson")),
       aes(x=buf50_area,y=pearson)) +
  geom_point() +
  theme_bw()

#### Density predictions: ####
# Density for an average island:
avg.isl.dens <- expand.grid(l.area.std = mean(dat$l.area.std),
                      soild15N.std = mean(dat$soild15N.std),
                      rocky.prop.std = mean(dat$rocky.prop.std),
                      sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                      year.std = mean(dat$year.std),
                      buf50_area = 1)

pred.avg.isl.dens <- predict(avg.density, avg.isl.dens, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred.avg.isl.dens$fit * 10000 # 20.70 individuals/hectare if everything else is at its mean. 
pred.avg.isl.dens$se.fit*1.96 * 10000 # 4.09

# Add one sd of area: 
big.isl.dens<- expand.grid(l.area.std = mean(dat$l.area.std) + sd(dat$l.area.std),
                       soild15N.std = mean(dat$soild15N.std),
                       rocky.prop.std = mean(dat$rocky.prop.std),
                       sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                       year.std = mean(dat$year.std),
                       buf50_area = 1)

pred.big.isl.dens <- predict(avg.density, big.isl.dens, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred.big.isl.dens$fit * 10000 # 17.17 
pred.small.isl.dens$fit * 10000 - pred.big.isl.dens$fit * 10000 # 3.5 individuals / sd area.
pred.big.isl.dens$se.fit*1.96*10000 

# Subtract one SD of area: 
small.isl.dens <- expand.grid(l.area.std = mean(dat$l.area.std) - sd(dat$l.area.std),
                            soild15N.std = mean(dat$soild15N.std),
                            rocky.prop.std = mean(dat$rocky.prop.std),
                            sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                            year.std = mean(dat$year.std),
                            buf50_area = 1)

pred.small.isl.dens <- predict(avg.density, small.isl.dens, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred.small.isl.dens$fit * 10000 # 24.96823
pred.small.isl.dens$se.fit*1.96*10000 # 5.75

# Add one sd of d15N: 
d15n.dens <- expand.grid(l.area.std = mean(dat$l.area.std),
                             soild15N.std = mean(dat$soild15N.std) + sd(dat$soild15N.std),
                             rocky.prop.std = mean(dat$rocky.prop.std),
                             sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                             year.std = mean(dat$year.std),
                             buf50_area = 1)

pred.d15n.dens <- predict(avg.density, low.d15n.dens, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred.d15n.dens$fit * 10000 
pred.d15n.dens$se.fit * 10000 * 1.96

