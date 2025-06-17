################################################################################
# Successfull simulation with 100,200,300 sample size and (z_1,y_1), y_1 is    #
# time dependent covariate and simulation scenerio considers low missing values# 
# in three different missing data mechanism MCAR, MAR, MNAR. functions we have #
# used are jm2s, jmLOCF, jmMI, jmwMI.                                          #
#                                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
################################################################################
library(MASS)
library(dplyr)
library(boot)
library(mice)
library(nlme)
library(survival)
library(JM)
library(JMtwostage)
library(Matrix)
library(stringr)
library(lme4)
library(insight)
library(timereg)
library(rsimsum)
library(ggpubr)
library(latex2exp)



simcode<-function(id=1:100,
                  t0=rep(0,10),
                  z_1=rep(100,0),
                  sd_b=rep(1,3),
                  rchange=0.4
){
  nid<-length(id)
  beta_1<-c(1.1,-0.003,-0.025)
  sd_b_1<-sd_b[1]
  y_1<-beta_1[1]+beta_1[2]*z_1+beta_1[3]*t0+rnorm(1,0,sd_b_1)+rnorm(nid,0,1)
  beta_surv<-c(-5,0.01,-0.5)#(intercept,z_1,y_1)
  lambda<-beta_surv[1]+beta_surv[2]*z_1+beta_surv[3]*y_1
  tfail<-rexp(length(id),exp(lambda))
  tchange<-rexp(length(id),rchange)
  tfail<-round(tfail,3)
  tchange<-round(tchange,3)
  tevent<-pmin(tfail,tchange)
  fevent<-tfail==tevent
  if(isTRUE(all(fevent)))
    return(list(cbind('start'=t0,
                      'stop'=t0+tevent,
                      'event'=fevent,
                      'id'=id,
                      'z_1'=z_1,
                      'y_1'=y_1)))
  
  c(list(cbind('start'=t0,'stop'=t0+tevent,'event'=fevent,'id'=id,'z_1'=z_1,'y_1'=y_1)),
    simcode(id=id[!fevent],t0=(t0+tevent)[!fevent],z_1=z_1[!fevent],sd_b=sd_b,rchange=0.02))
  
}


#n=100;Miss="Low";Miss_type="MAR"
#n=200;Miss="Low";Miss_type="MAR"
#n=100;Miss="Low";Miss_type="MCAR"
#n=300;Miss="Low";Miss_type="MNAR"
simdata<-function(n=100,Miss="Low",Miss_type="MCAR"){
  n<-n
  simulated_data<-simcode(id=1:n,
                          t0=rep(0,n),
                          z_1=rnorm(n,40,10),
                          sd_b=rep(1,3),
                          rchange=0.02)
  
  simulated_data <- as.data.frame(do.call(rbind, simulated_data))
  simulated_data<-simulated_data[order(simulated_data$id),]
  dat<-simulated_data[simulated_data$start!=simulated_data$stop,]
  
  
  if(Miss!="None"){
    if(Miss=="Low"){
      if(Miss_type=="MCAR"){
        dat$Prob_1<-(dat$start!=0)*rep(0.2,nrow(dat))
      }
      if(Miss_type=="MAR"){
        #miss_beta<-list(c(-1.9,1.2,-0.01))
        miss_beta<-list(c(-0.6,-0.001,-0.01))
        dat$Prob_1<-(dat$start!=0)*inv.logit(miss_beta[[1]][1]+miss_beta[[1]][2]*dat$z_1+
                                               miss_beta[[1]][3]*dat$start)
      }
      if(Miss_type=="MNAR"){
        miss_beta<-list(c(-1.5,0.01,0.4))
        dat$Prob_1<-(dat$start!=0)*inv.logit(miss_beta[[1]][1]+miss_beta[[1]][2]*dat$z_1+
                                               miss_beta[[1]][3]*dat$y_1)
      }
    }
    
    
    if(Miss=="High"){
      if(Miss_type=="MCAR"){
        dat$Prob_1<-(dat$start!=0)*rep(0.45,nrow(dat))
      }
      if(Miss_type=="MAR"){
        #miss_beta<-list(c(-1.9,1.2,-0.01))
        miss_beta<-list(c(-0.08,-0.001,-0.001))
        dat$Prob_1<-(dat$start!=0)*inv.logit(miss_beta[[1]][1]+miss_beta[[1]][2]*dat$z_1+
                                               miss_beta[[1]][3]*dat$start)
      }
      if(Miss_type=="MNAR"){
        if(n==100){
          miss_beta<-list(c(-1.4,0.01,-1.7))
        }else if(n==200){
          miss_beta<-list(c(-1.1,0.01,-1.6)) 
        }else{
          miss_beta<-list(c(-1.3,0.01,-1.7))  
        }
        dat$Prob_1<-(dat$start!=0)*inv.logit(miss_beta[[1]][1]+miss_beta[[1]][2]*dat$z_1+
                                               miss_beta[[1]][3]*dat$y_1)
      }
    }
    
    
    flag<-0
    while(flag==0)
    {
      dat$mis_1<-rbinom(nrow(dat),1,dat$Prob_1)
      gm<-sum(dat$mis_1)/nrow(dat)
      #print(gm)
      if(Miss=="Low"){
        if(round(gm,2)>=0.15&round(gm,2)<=0.28) flag<-1
      }else if(Miss=="High"){
        #round(gm,2)
        if(round(gm,2)>=0.29&round(gm,1)<=0.37) flag<-1
      }
    }
    dat[dat$mis_1==1,"y_1"]<-NA
    
  }
  
  dat<-dat[names(dat)%in%c("start","stop","event","id","z_1","y_1")]
  dat<-dat[order(dat$id),]
  id_freq<-data.frame(table(dat$id))
  survival_time<-c()
  survival_status<-c()
  for(i in 1:nrow(id_freq)){
    data1<-tail(dat[dat$id==i,],1)
    survival_time[i]<-data1$stop
    survival_status[i]<-data1$event
  }
  
  dat1<-list()
  dat1$id<-dat$id
  dat1$Time<-dat$start
  dat1$z_1<-dat$z_1
  dat1$y_1<-dat$y_1
  
  
  dat1$survival_time<-rep(survival_time,times=id_freq$Freq)
  dat1$survival_status<-rep(survival_status,times=id_freq$Freq)
  dat1<-as.data.frame(dat1)
  return(dat1)
}

################################################################################
# This part will create list containing datatsets with different missing data  #
# setting like High or Low, MCAR or MAR or MNAR, sample size 100, 200, 300     #
# User can use different parameter as par their need                           #
# Note :Run the for loop one case at a time, it might break the session if you #
# run all at once.                                                             #
# If for some case the data generation takes time, then hit the stop button    #
# in the console and check how many dataset left, and rerun from that iteration#
# number                                                                       #
#                                                                              #
################################################################################
reps <- 1:100
data_set <- list()

for(i in 1:100){
  cat("Start simulation no: ",i,'\n')
   data_set[["n=100_high_mcar"]][[i]]<-simdata(n=100,Miss="High",Miss_type="MCAR")
  # data_set[["n=100_high_mar"]][[i]]<-simdata(n=100,Miss="High",Miss_type="MAR")
  # data_set[["n=100_high_mnar"]][[i]]<-simdata(n=100,Miss="High",Miss_type="MNAR")
  # data_set[["n=200_high_mcar"]][[i]]<-simdata(n=200,Miss="High",Miss_type="MCAR")
  # data_set[["n=200_high_mar"]][[i]]<-simdata(n=200,Miss="High",Miss_type="MAR")
  # data_set[["n=200_high_mnar"]][[i]]<-simdata(n=200,Miss="High",Miss_type="MNAR")
  # data_set[["n=300_high_mcar"]][[i]]<-simdata(n=300,Miss="High",Miss_type="MCAR")
  # data_set[["n=300_high_mar"]][[i]]<-simdata(n=300,Miss="High",Miss_type="MAR")
  # data_set[["n=300_high_mnar"]][[i]]<-simdata(n=300,Miss="High",Miss_type="MNAR")
  cat("End simulation no: ",i,'\n')
}


fit_models<-function(data,model){
  longdata<-data
  survdata<-data[!duplicated(data$id),]
  if(model=="locf"){
    fit<-jmLOCF(ldata=longdata,
                sdata=survdata,
                timeDep=c("y_1"),
                visitTime="Time",
                coxModel=Surv(survival_time,survival_status)~z_1+td(y_1),
                model="Cox",
                id="id")
    
  }else if(model=="2s"){
    fit<-jm2s(ldata=longdata,
              sdata=survdata,
              timeDep=c("y_1"),
              impModel=list(y_1~z_1+Time+(1|id)),
              visitTime="Time",
              coxModel=Surv(survival_time,survival_status)~z_1+td(y_1),
              model="Cox",
              id="id")
  }else if(model=="mi"){
    fit<-jmMI(ldata=longdata,
              sdata=survdata,
              timeDep=c("y_1"),
              impModel=list(y_1~z_1+Time+(1|id)),
              visitTime="Time",
              coxModel=Surv(survival_time,survival_status)~z_1+td(y_1),
              model="Cox",
              id="id")
  }else if(model=="wmi"){
    fit<-jmwMI(ldata=longdata,
               sdata=survdata,
               timeDep=c("y_1"),
               impModel=list(y_1~z_1+Time+(1|id)),
               visitTime="Time",
               ipwModel=list(~z_1+Time+(1|id)),
               coxModel=Surv(survival_time,survival_status)~z_1+td(y_1),
               model="Cox",
               id="id")
  }else if(model=="cc"){
    fit<-jmCC(ldata=longdata,
              sdata=survdata,
              timeDep=c("y_1"),
              visitTime="Time",
              coxModel=Surv(survival_time,survival_status)~z_1+td(y_1),
              model="Cox",
              id="id")
  }
  
  fit$res[c(1:2)]
}



results <- list()

#n=100, low missing data, MCAR
results[["n=100,cc,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mcar"]],
    FUN = fit_models,
    model = "cc"
  )
)

results[["n=100,locf,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mcar"]],
    FUN = fit_models,
    model = "locf"
  )
)

results[["n=100,twos,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mcar"]],
    FUN = fit_models,
    model = "2s"
  )
)

results[["n=100,mi,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mcar"]],
    FUN = fit_models,
    model = "mi"
  )
)

results[["n=100,wmi,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mcar"]],
    FUN = fit_models,
    model = "wmi"
  )
)


## n=100, low missing data, MAR
results[["n=100,cc,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mar"]],
    FUN = fit_models,
    model = "cc"
  )
)

results[["n=100,locf,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mar"]],
    FUN = fit_models,
    model = "locf"
  )
)

results[["n=100,twos,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mar"]],
    FUN = fit_models,
    model = "2s"
  )
)

results[["n=100,mi,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mar"]],
    FUN = fit_models,
    model = "mi"
  )
)

results[["n=100,wmi,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mar"]],
    FUN = fit_models,
    model = "wmi"
  )
)

## n=100, low missing data, MNAR
results[["n=100,cc,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mnar"]],
    FUN = fit_models,
    model = "cc"
  )
)

results[["n=100,locf,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mnar"]],
    FUN = fit_models,
    model = "locf"
  )
)

results[["n=100,twos,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mnar"]],
    FUN = fit_models,
    model = "2s"
  )
)

results[["n=100,mi,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mnar"]],
    FUN = fit_models,
    model = "mi"
  )
)

results[["n=100,wmi,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=100_high_mnar"]],
    FUN = fit_models,
    model = "wmi"
  )
)

################################################################################
## n=200, high missing data, MCAR

results[["n=200,cc,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mcar"]],
    FUN = fit_models,
    model = "cc"
  )
)

results[["n=200,locf,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mcar"]],
    FUN = fit_models,
    model = "locf"
  )
)

results[["n=200,twos,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mcar"]],
    FUN = fit_models,
    model = "2s"
  )
)

results[["n=200,mi,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mcar"]],
    FUN = fit_models,
    model = "mi"
  )
)

results[["n=200,wmi,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mcar"]],
    FUN = fit_models,
    model = "wmi"
  )
)

## n=200, high missing data, MAR

results[["n=200,cc,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mar"]],
    FUN = fit_models,
    model = "cc"
  )
)

results[["n=200,locf,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mar"]],
    FUN = fit_models,
    model = "locf"
  )
)

results[["n=200,twos,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mar"]],
    FUN = fit_models,
    model = "2s"
  )
)

results[["n=200,mi,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mar"]],
    FUN = fit_models,
    model = "mi"
  )
)

results[["n=200,wmi,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mar"]],
    FUN = fit_models,
    model = "wmi"
  )
)

## n=200, high missing data, MNAR

results[["n=200,cc,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mnar"]],
    FUN = fit_models,
    model = "cc"
  )
)

results[["n=200,locf,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mnar"]],
    FUN = fit_models,
    model = "locf"
  )
)

results[["n=200,twos,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mnar"]],
    FUN = fit_models,
    model = "2s"
  )
)

results[["n=200,mi,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mnar"]],
    FUN = fit_models,
    model = "mi"
  )
)

results[["n=200,wmi,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=200_high_mnar"]],
    FUN = fit_models,
    model = "wmi"
  )
)

################################################################################
# n=300, high missing data, MCAR

results[["n=300,cc,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mcar"]],
    FUN = fit_models,
    model = "cc"
  )
)

results[["n=300,locf,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mcar"]],
    FUN = fit_models,
    model = "locf"
  )
)

results[["n=300,twos,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mcar"]],
    FUN = fit_models,
    model = "2s"
  )
)

results[["n=300,mi,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mcar"]],
    FUN = fit_models,
    model = "mi"
  )
)

results[["n=300,wmi,high,MCAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mcar"]],
    FUN = fit_models,
    model = "wmi"
  )
)

# n=300, high missing data, MAR

results[["n=300,cc,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mar"]],
    FUN = fit_models,
    model = "cc"
  )
)

results[["n=300,locf,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mar"]],
    FUN = fit_models,
    model = "locf"
  )
)

results[["n=300,twos,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mar"]],
    FUN = fit_models,
    model = "2s"
  )
)

results[["n=300,mi,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mar"]],
    FUN = fit_models,
    model = "mi"
  )
)

results[["n=300,wmi,high,MAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mar"]],
    FUN = fit_models,
    model = "wmi"
  )
)

# n=300, high missing data, MNAR

results[["n=300,cc,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mnar"]],
    FUN = fit_models,
    model = "cc"
  )
)

results[["n=300,locf,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mnar"]],
    FUN = fit_models,
    model = "locf"
  )
)

results[["n=300,twos,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mnar"]],
    FUN = fit_models,
    model = "2s"
  )
)

results[["n=300,mi,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mnar"]],
    FUN = fit_models,
    model = "mi"
  )
)

results[["n=300,wmi,high,MNAR"]] <- do.call(
  rbind.data.frame,
  lapply(
    X = data_set[["n=300_high_mnar"]],
    FUN = fit_models,
    model = "wmi"
  )
)

################################################################################

relhaz <- do.call(rbind.data.frame, results)
relhaz$par <- rep(c("y_1","z_1"),times=9*5*length(reps))
relhaz$Model<-rep(c("cc","locf","twos","mi","wmi"),each=9*2*length(reps))
relhaz$SampleSize<-rep(c("n=100,mcar","n=100,mar","n=100,mnar",
                         "n=200,mcar","n=200,mar","n=200,mnar",
                         "n=300,mcar","n=300,mar","n=300,mnar"),
                       each=2*length(reps))
#beta_surv<-c(-5,0.01,-0.5)#(intercept,z_1,y_1)
ms<-multisimsum(data=relhaz,par="par",true=c(z_1=0.01,y_1=-0.5),
                estvarname = "logHR",se="SE",methodvar = "Model",x=TRUE,ref="twos",by="SampleSize")
sms<-summary(ms)


library(ggplot2)
autoplot(sms,par="z_1")
autoplot(sms,par="z_1",type="lolly")
autoplot(sms,par="z_1",type="zip")
autoplot(sms,par="z_1",type="est")
autoplot(sms,par="z_1",type="se")
autoplot(sms,par="z_1",type="est_ba")
autoplot(sms,par="z_1",type="se_ba")
autoplot(sms,par="z_1",type="est_density")
autoplot(sms,par="z_1",type="se_density")
autoplot(sms,par="z_1",type="se_hex")
autoplot(sms,par="z_1",type="est_ridge")
autoplot(sms,par="z_1",type="se_ridge")
autoplot(sms,par="z_1",type="heat")
autoplot(sms,par="z_1",type="nlp")


autoplot(sms,par="y_1")
autoplot(sms,par="y_1",type="lolly")
autoplot(sms,par="y_1",type="zip")
autoplot(sms,par="y_1",type="est")
autoplot(sms,par="y_1",type="se")
autoplot(sms,par="y_1",type="est_ba")
autoplot(sms,par="y_1",type="se_ba")
autoplot(sms,par="y_1",type="est_density")
autoplot(sms,par="y_1",type="se_density")
autoplot(sms,par="y_1",type="se_hex")
autoplot(sms,par="y_1",type="est_ridge")
autoplot(sms,par="y_1",type="se_ridge")
autoplot(sms,par="y_1",type="heat")
autoplot(sms,par="y_1",type="nlp")


plot_z1_cover_heat<-autoplot(sms, par = "z_1", type = "heat", stats = "cover")
ggsave("plot_z1_high/plot_z1_cover_heat.jpeg",plot=plot_z1_cover_heat,dpi=300,width=10,height=8,units=c('in'))
plot_z1_bias<-autoplot(sms, par = "z_1", type = "lolly", stats = "bias")
ggsave("plot_z1_high/plot_z1_bias.jpeg",plot=plot_z1_bias,dpi=300,width=10,height=8,units=c('in'))
plot_z1_mse<-autoplot(sms, par = "z_1", type = "lolly", stats = "mse")
ggsave("plot_z1_high/plot_z1_mse.jpeg",plot=plot_z1_mse,dpi=300,width=10,height=8,units=c('in'))
plot_z1_se<-autoplot(sms, par = "z_1", type="se_ridge")
ggsave("plot_z1_high/plot_z1_se.jpeg",plot=plot_z1_mse,dpi=300,width=10,height=8,units=c('in'))
plot_z1_empse<-autoplot(sms, par = "z_1", type = "lolly", stats = "empse")
ggsave("plot_z1_high/plot_z1_empse.jpeg",plot=plot_z1_empse,dpi=300,width=10,height=8,units=c('in'))



plot_y1_cover_heat<-autoplot(sms, par = "y_1", type = "heat", stats = "cover")
ggsave("plot_y1_high/plot_y1_cover_heat.jpeg",plot=plot_y1_cover_heat,dpi=300,width=10,height=8,units=c('in'))
plot_y1_bias<-autoplot(sms, par = "y_1", type = "lolly", stats = "bias")
ggsave("plot_y1_high/plot_y1_bias.jpeg",plot=plot_y1_bias,dpi=300,width=10,height=8,units=c('in'))
plot_y1_mse<-autoplot(sms, par = "y_1", type = "lolly", stats = "mse")
ggsave("plot_y1_high/plot_y1_mse.jpeg",plot=plot_y1_mse,dpi=300,width=10,height=8,units=c('in'))
plot_y1_se<-autoplot(sms, par = "y_1", type="se_ridge")
ggsave("plot_y1_high/plot_y1_se.jpeg",plot=plot_y1_mse,dpi=300,width=10,height=8,units=c('in'))
plot_y1_empse<-autoplot(sms, par = "y_1", type = "lolly", stats = "empse")
ggsave("plot_z1_high/plot_y1_empse.jpeg",plot=plot_y1_empse,dpi=300,width=10,height=8,units=c('in'))




com_cover_heat<-ggarrange(plot_z1_cover_heat,plot_y1_cover_heat, ncol=2, nrow=1, common.legend = TRUE, legend="right")
ggsave("plot_com_high/com_cover_heat.jpeg",plot=com_cover_heat,dpi=300,width=10,height=7,units=c('in'))

com_bias<-ggarrange(plot_z1_bias,plot_y1_bias, ncol=1, nrow=2, common.legend = TRUE, legend="right")
ggsave("plot_com_high/com_bias.jpeg",plot=com_bias,dpi=300,width=7,height=10,units=c('in'))

com_density<-ggarrange(plot_z1_ridge,plot_y1_ridge, ncol=2, nrow=1, common.legend = TRUE, legend="right")
ggsave("plot_com_high/com_density.jpeg",plot=com_density,dpi=300,width=9,height=12,units=c('in'))


Simulation_plot_data <- read.csv("C:/rproject/jm2stage_20_01_25/Simulation_plot_data_highMiss.csv")

mse_data_high<-Simulation_plot_data[Simulation_plot_data$stat=="empse",]

mse100_z1_high<-ggplot(mse_data_high[mse_data_high$SampleSize=="n=100"&mse_data_high$par=="z_1",],aes(y=est,x=Model))+
  geom_point(aes(shape=Model,color=Miss_type),size=2)+
  ylim(c(0,2))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   axis.line = element_line(colour = "black",linetype = 'solid',size = 0.7),
                   axis.text.x = element_text(face='bold',size=7),
                   axis.text.y = element_text(face='bold',size=7),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(face='bold',size=10)
                   #axis.ticks.y = element_blank()
  )+ylab(TeX("Emperical MSE of $\\gamma_1$"))

mse100_y1_high<-ggplot(mse_data_high[mse_data_high$SampleSize=="n=100"&mse_data_high$par=="y_1",],aes(y=est,x=Model))+
  geom_point(aes(shape=Model,color=Miss_type),size=2)+
  ylim(c(0,25))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   axis.line = element_line(colour = "black",linetype = 'solid',size = 0.7),
                   axis.text.x = element_text(face='bold',size=7),
                   axis.text.y = element_text(face='bold',size=7),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(face='bold',size=10)
                   #axis.ticks.y = element_blank()
  )+ylab(TeX("Emperical MSE of $\\gamma_2$"))



mse200_z1_high<-ggplot(mse_data_high[mse_data_high$SampleSize=="n=200"&mse_data_high$par=="z_1",],aes(y=est,x=Model))+
  geom_point(aes(shape=Model,color=Miss_type),size=2)+
  ylim(c(0,2))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   axis.line = element_line(colour = "black",linetype = 'solid',size = 0.7),
                   axis.text.x = element_text(face='bold',size=7),
                   axis.text.y = element_text(face='bold',size=7),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(face='bold',size=10)
                   #axis.ticks.y = element_blank()
  )+ylab(TeX("Emperical MSE of $\\gamma_1$"))

mse200_y1_high<-ggplot(mse_data_high[mse_data_high$SampleSize=="n=200"&mse_data_high$par=="y_1",],aes(y=est,x=Model))+
  geom_point(aes(shape=Model,color=Miss_type),size=2)+
  ylim(c(0,25))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   axis.line = element_line(colour = "black",linetype = 'solid',size = 0.7),
                   axis.text.x = element_text(face='bold',size=7),
                   axis.text.y = element_text(face='bold',size=7),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(face='bold',size=10)
                   #axis.ticks.y = element_blank()
  )+ylab(TeX("Emperical MSE of $\\gamma_2$"))

mse300_z1_high<-ggplot(mse_data_high[mse_data_high$SampleSize=="n=300"&mse_data_high$par=="z_1",],aes(y=est,x=Model))+
  geom_point(aes(shape=Model,color=Miss_type),size=2)+
  ylim(c(0,2))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   axis.line = element_line(colour = "black",linetype = 'solid',size = 0.7),
                   axis.text.x = element_text(face='bold',size=7),
                   axis.text.y = element_text(face='bold',size=7),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(face='bold',size=10)
                   #axis.ticks.y = element_blank()
  )+ylab(TeX("Emperical MSE of $\\gamma_1$"))

mse300_y1_high<-ggplot(mse_data_high[mse_data_high$SampleSize=="n=300"&mse_data_high$par=="y_1",],aes(y=est,x=Model))+
  geom_point(aes(shape=Model,color=Miss_type),size=2)+
  ylim(c(0,25))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   axis.line = element_line(colour = "black",linetype = 'solid',size = 0.7),
                   axis.text.x = element_text(face='bold',size=7),
                   axis.text.y = element_text(face='bold',size=7),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(face='bold',size=10)
                   #axis.ticks.y = element_blank()
  )+ylab(TeX("Emperical MSE of $\\gamma_2$"))

com_mse_high<-ggarrange(mse100_z1_high,
                   mse100_y1_high,
                   mse200_z1_high,
                   mse200_y1_high,
                   mse300_z1_high,
                   mse300_y1_high,
                   ncol=2, 
                   nrow=3, 
                   common.legend = TRUE,
                   legend="right")
ggsave("plot_com_high/com_mse_high.jpeg",plot=com_mse_high,dpi=300,width=9,height=12,units=c('in'))




