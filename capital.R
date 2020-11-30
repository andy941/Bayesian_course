# small/big
# low/high
# Missing data are indicated by -99.99 or -999.

library(dplyr)
library(ggplot2)
library(rjags)
library(tidyr)
library(gridExtra)

dat <- read.table("6_Portfolios_2x3.csv",
                  sep = ',',
                  header = F)
Size <- rep(c("SmallCap", "LargeCap"), each = 3)
Value <- rep(c("Growth", "Core", "Value"), 2)
colnames(dat) <- c("Date", paste0(Size, "_", Value))

return_calc <- function(x, capital = 10000) {
  
  final = c()
  
  for (i in 1:length(x)) {
    final[i] = capital + (capital/100) * x[i]
    capital = capital + (capital/100) * x[i] 
  }
  return(final)
}

return_analysis <- function(x,
                            capital = 10000,
                            analysis = "window",
                            window = 120,
                            n_bootstraps = 1000,
                            replace = T) {
  
  
  if (analysis == "window") {
    FUN <- function(x) {
      result = c()
      for (i in 1:(length(x) - window)) {
        data = x[i:(i + window)]
        result[i] = tail(return_calc(data, capital = capital), 1)
      }
      return(result)
    }
    final <- apply(x[,-1], FUN = FUN, MARGIN = 2) %>%
      as.data.frame()
    final$Date = paste0(x$Date[1:(length(x$Date)-window)], "-",x$Date[(window+1):length(x$Date)])
    
    return(final)
  } 
  
  
  if (analysis == "bootstrap") {
    final = data.frame()
    for (i in 1:n_bootstraps) {
      sample = sample(1:length(x[,1]), size = window, replace = replace)
      tmp = tail(apply(x[sample,-1],
                       FUN = return_calc,
                       capital = capital, 
                       MARGIN = 2),
                 1)
      final = rbind.data.frame(final, tmp)
    }
    return(final)
  }
  return(final)
}

# Get the dataframes for 10 years return widows
windows <- return_analysis(dat, analysis = "window", capital = 100)

windows_Long <- data.frame()
for (i in 1:6) {
  long = data.frame(Return = windows[,-7][,i])
  long$Date = windows$Date
  long$Size = Size[i]
  long$Value = Value[i]
  windows_Long = rbind.data.frame(windows_Long, long)
}

windows_Long$Group <- paste0(windows_Long$Size, "_", windows_Long$Value)
levels(windows_Long$Size)  <- c("LargeCap", "SmallCap")
levels(windows_Long$Value)  <- c("Growth", "Core", "Value")

windows_Long$Size <- as.factor(windows_Long$Size)
windows_Long$Value <- as.factor(windows_Long$Value)

# Drop Core from the model
windows_Long <- windows_Long %>%
  filter(Value != "Core")  %>%
  droplevels.data.frame()

# Explore data
theme <- theme(axis.text.x = element_blank(),
               axis.ticks.x =  element_blank(),
               panel.background = element_blank(),panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"))

plot(windows$SmallCap_Value - windows$LargeCap_Growth,
     type = 'l',
     ylab = "%SmallCap_Value - %LargeCap_Growth", 
     xlab = "Historical 10 year periods")
abline(a = 0, b = 0, col='red')

ggplot(windows_Long, aes(x=Date, y=Return, group = Group)) +
  geom_line(aes(colour = Group)) + theme

theme <- theme(
               panel.background = element_blank(),panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"))

ggplot(windows_Long, aes(x=Return, colour = Group)) +
  geom_density() + theme
  
# Bayesian model using rjags
mod_string = " model {
    for (i in 1:length(Return)) {
        Return[i] ~ dgamma(shape[i], rate[i])
        
        shape[i] = b0_s + bs[1]*Size[i] + bs[2]*Value[i]
        rate[i] = b0_r + br[1]*Size[i] + br[2]*Value[i]
    }
    
    for (m in 1:2) {
        br[m] ~ dunif(0,1e5)
    }
    
    for (n in 1:2) {
        bs[n] ~ dunif(0,1e5)
    }
    
    b0_s ~ dunif(0,1e5)
    b0_r ~ dunif(0,1e5)
    
}"

data_jags = list()
data_jags$Return = windows_Long$Return
data_jags$Size = as.numeric(windows_Long$Size) -1
data_jags$Value = as.numeric(windows_Long$Value) -1


params = c("b0_s", "bs", "b0_r", "br")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=4)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params, 
                       n.iter=6e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
par(mar=c(2,2,2,2))
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
#autocorr.plot(mod_sim)

effectiveSize(mod_sim)

# compute the DIC for the model
dic = dic.samples(mod, n.iter=1e3)

# Try to fit nother model with interaction terms
#################################################----------------------------
mod_string2 = " model {
    for (i in 1:length(Return)) {
        Return[i] ~ dgamma(shape[i], rate[i])
        
        shape[i] = b0_s + bs[1]*Size[i] + bs[2]*Value[i] + bs[3]*Size[i]*Value[i]
        rate[i] = b0_r + br[1]*Size[i] + br[2]*Value[i] + br[3]*Size[i]*Value[i]
    }
    
    for (m in 1:3) {
        br[m] ~ dunif(0,1e5)
    }
    
    for (n in 1:3) {
        bs[n] ~ dunif(0,1e5)
    }
    
    b0_s ~ dunif(0,1e5)
    b0_r ~ dunif(0,1e5)
    
}"


mod2 = jags.model(textConnection(mod_string2), data=data_jags, n.chains=4)
update(mod, 1e3)

mod_sim2 = coda.samples(model=mod2,
                       variable.names=params, 
                       n.iter=6e3)
mod_csim2 = as.mcmc(do.call(rbind, mod_sim2))

## convergence diagnostics
par(mar=c(2,2,2,2))
plot(mod_sim2)

gelman.diag(mod_sim2)
autocorr.diag(mod_sim2)
autocorr.plot(mod_sim2)

effectiveSize(mod_sim2)

# compute the DIC for the model
dic2 = dic.samples(mod2, n.iter=1e3)


# ANALYSIS
###________________________________________________________________________________
p <- colMeans(mod_csim2) %>% 
  t() %>%
  as.data.frame()
size = c(0,1)
value = c(0,1)
data = c()
means = c()
var = c()
mod_csim2 = as.data.frame(mod_csim2)
mean_estimates = c()

for (s in size) {
  for (v in value) {
    shape = p$b0_s + s*p$`bs[1]` + v*p$`bs[2]` + s*v*p$`bs[3]`
    rate = p$b0_r + s*p$`br[1]` + v*p$`br[2]` + s*v*p$`br[3]`
    draws = rgamma(1000000, shape, rate)
    data = cbind(data, draws)
    
    mean_estimates = append(mean_estimates, shape/rate)
    
    shape_mm = mod_csim2$b0_s + s*mod_csim2$`bs[1]` + v*mod_csim2$`bs[2]` + s*v*mod_csim2$`bs[3]`
    rate_mm = mod_csim2$b0_r + s*mod_csim2$`br[1]` + v*mod_csim2$`br[2]` + s*v*mod_csim2$`br[3]`
    means = cbind(means, shape_mm/rate_mm)
    var = cbind(var, shape_mm/(rate_mm**2))
  }
}

colnames(data) = levels(as.factor(windows_Long$Group))
data = as.data.frame(data)
colnames(means) = levels(as.factor(windows_Long$Group))
means = as.data.frame(means)
colnames(var) = levels(as.factor(windows_Long$Group))
var = as.data.frame(var)
names(mean_estimates) = levels(as.factor(windows_Long$Group))

# Residual analysis
residuals = windows[, c(4,6,1,3)] - mean_estimates
residuals = pivot_longer(residuals, 
                         cols = everything(), 
                         names_to = "Group",
                         values_to = "res")

ggplot(residuals, aes(Group, res, color=Group)) +
  geom_violin() +
  theme

data_long = pivot_longer(data, 
                         cols = everything(), 
                         names_to = "Group",
                         values_to = "Return")
theme <- theme(
               panel.background = element_blank(),panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"))

ggplot(data_long, aes(x=Return, colour = Group)) +
  geom_density() + theme

means_long = means %>% pivot_longer(cols = everything(), 
                       names_to = "Group",
                       values_to = "Means")

p1 <- ggplot(means_long, aes(x=Group, y=Means, color = Group)) +
  geom_violin(show.legend = F) + theme

var_long = var %>% pivot_longer(cols = everything(), 
                       names_to = "Group",
                       values_to = "var")

p2 <- ggplot(var_long, aes(x=Group, y=var, color = Group)) +
  geom_violin(show.legend = F) + theme

grid.arrange(p1,p2, nrow = 1)

plot(density(data$SmallCap_Value))
points(density(windows$SmallCap_Value), type = 'l', col = 'red')

prob_LG = mean(data$LargeCap_Growth > windows$LargeCap_Growth[1010])
prob_SV = mean(data$SmallCap_Value < windows$SmallCap_Value[1010])
prob = mean(data$SmallCap_Value > data$LargeCap_Growth)
  
