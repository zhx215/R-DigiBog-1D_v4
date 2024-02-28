### The R code for 1-D DigiBog model (in Fortran) presented in Morris et al. (2015), modified to include seasonality, surface pond, snowmelt, and variable density ###
### Author: Zhengyu Xia ###

## start ##
setwd("C:/Users/zhyxi/Dropbox/digibog")
rm(list=ls())
set.seed(888)

## parameter ##
t_extent <- 3000 # yr
lateral_extent <- 50000 # cm, radius

oxic_decay_base <- 0.042 # yr-1 (from Morris et al., 2015)
anoxic_decay_base <- 0.0002 # yr-1 (from Morris et al., 2015)
base_temp <- 6.29 # degree C (from Morris et al., 2015)
Q10_oxic <- 2.5 # (from Morris et al., 2015)
Q10_anoxic <- 2.5 # (from Morris et al., 2015)

porosity <- 0.3 # (from Morris et al., 2015)
k_param_a <- 31740 # cm yr-1, hydraulic conductivity model parameter (from Morris et al., 2015)
k_param_b <- 8 # hydraulic conductivity model parameter (from Morris et al., 2015)
drain_effi <- 1 # draining efficiency, circular = 2, ellipse = 1 
x_factor <- 0.5 # recalcitrance (from Morris et al., 2015)
ddf <- 0.15 # mm degree C-1 day-1, degree day factor to calculate snowmelt
pond_size <- 2 # cm

density_mode <- "variable"
initial_density <- 0.04 # g cm-3, it is the initial density in variable mode or the constant density in constant mode

## climate forcing data ##
source("R-DigiBog-1D_v4 climate input.R")

## function ##
prod <- function(z,t){ # water table depth cm, air temp degree C
  if ((t>0)&(z>=-6.3)&(z<=66.8)) {
    return(0.0001*(9.3+1.33*z-0.022*z^2)^2*(0.1575*t+0.009132)) # g cm-2 yr-1
  } else {
    return(0.0000001) # prevent no production
  }
}
dHdt <- function(p,et,melt,theta_d,K_sum,H,L,childs_factor){ #precip cm/yr, ET cm/yr, snowmelt cm/yr, porosity, hydraulic conductivity sum cm yr-1, water table height cm, radius cm, draining efficiency 
  return(((p+melt-et)-(childs_factor*K_sum*H)/(L^2))/theta_d)
}
aet_pet_ratio <- function(wtd,threshold_wtd,least_ratio){ # when AET/PET ratio decreases cm, the minimum ratio
  if (wtd<=threshold_wtd) {
    return(1)
  } else {
    a <- (1-least_ratio)/(threshold_wtd-66.8)
    b <- ((threshold_wtd-66.8)*1-(1-least_ratio)*threshold_wtd)/(threshold_wtd-66.8)
    return(a*wtd+b)
  }
}
derived_density <- function(layer_k){
  if (density_mode=="variable"){
    initial_k <- k_param_a*exp(k_param_b*1)/3153600000 # m s-1
    current_k <- layer_k/3153600000 #m s-1
    density_change <- (log10(current_k)-log10(initial_k))*-20 # kg m-3
    return(initial_density+density_change/1000) # g cm-3
  } else {
    return(initial_density)
  }
}

## run model ##
source("R-DigiBog-1D_v4 algorithm.R")

## plot results ##
source("R-DigiBog-1D_v4 plot.R")