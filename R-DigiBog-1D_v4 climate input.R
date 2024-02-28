## read climate station data ##
monthly_record <- read.csv("monthly_record.csv")

## set k_value for Thornthwaite PET calculation ##
k_value <- c(71,84,98,114,128,136,133,121,106,90,76,68)/100 # 50N
# k_value <- c(80,89,99,110,120,125,123,115,104,93,83,78)/100 # 40N

## set monthly days ##
day <- c(31,28,31,30,31,30,31,31,30,31,30,31)

## get 12-monthly climatology ##
temp_12mn <- precip_12mn <- vector()
temp_12mn_sd <- precip_12mn_sd <- vector()
for (i in 1:12){
  temp_12mn[i] <- mean(monthly_record[,3][which(monthly_record[,2]==i)])
  precip_12mn[i] <- mean(monthly_record[,4][which(monthly_record[,2]==i)])
  temp_12mn_sd[i] <- sd(monthly_record[,3][which(monthly_record[,2]==i)])
  precip_12mn_sd[i] <- sd(monthly_record[,4][which(monthly_record[,2]==i)])
}
# note: you can input your own 12-monthly temperature and precipitation data

## calculate monthly pet ##
monthly_pet <- function(x1){
  temp <- x1
  temp[temp<0] <- 0
  I <- (temp/5)^1.514
  J <- sum(I)
  c <- 0.000000675*J^3-0.0000771*J^2+0.01792*J+0.49239
  pet <- k_value*1.6*(10*temp/J)^c # cm
  return(pet*10) # mm
}

## produce daily temp ##
generate_daily_temp_pet <- function(x1,x2){ # 12-monthly temp, 12-monthly temp sd
  t_12 <- t_365 <- vector()
  for (i in 1:12){
    t_12[i] <- rnorm(1,x1[i],x2[i])
  }
  for (i in 1:12){
    if (i != 12){
      coef <- summary(lm(c(t_12[i],t_12[i+1])~c(1,day[i])))$coefficients
      t_365 <- c(t_365,seq(1,day[i],1)*coef[2]+coef[1])
    } else {
      coef <- summary(lm(c(t_12[i],t_12[1])~c(1,day[i])))$coefficients
      t_365 <- c(t_365,seq(1,day[i],1)*coef[2]+coef[1])
    }
  }
  pet_12 <- monthly_pet(t_12)
  pet_365 <- vector()
  for (i in 1:12){
    pet_365 <- c(pet_365,rep(pet_12[i]/day[i],day[i]))
  }
  return(c(t_365[351:365],t_365[1:350],pet_365)) # mm, daily tepmerature + daily PET
}

## produce daily precip ##
generate_daily_precip <- function(x1,x2){  # 12-monthly precip, 12-monthly precip sd
  p_12 <- p_365 <- vector()
  for (i in 1:12){
    p_12[i] <- rnorm(1,x1[i],x2[i])
  }
  p_12[p_12<0] <- 0
  p_12[p_12>=750] <- 749.9
  max_mm <- 50 # maximum possible daily precipitation events
  min_mm <- 1 # minimum possible daily precipitation events
  for (i in 1:12){
    if (as.integer(p_12[i]/max_mm)>2){
      min_event <- as.integer(p_12[i]/max_mm)
    } else {
      min_event <- 2
    }
    if (as.integer(p_12[i]/min_mm)<15){
      max_event <- as.integer(p_12[i]/min_mm)
    } else {
      max_event <- 15
    }
    if (p_12[i]<2) {min_event <- 1; max_event <- 1}
    p_generator <- rep(0,day[i])
    event <- sample(min_event:max_event,1)
    p_generator[sort(sample(1:day[i],event))] <- p_12[i]/event
    p_365 <- c(p_365,p_generator)
  }
  return(p_365)
}