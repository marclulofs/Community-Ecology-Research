setwd("~/Theo/codes+data")

library(picante)
library(pspline)
library(scales)

color = c("azure3", "deepskyblue2", "blue4", "sienna1", "red3","yellowgreen","green4")
###################################
## Temperature data and function ##
###################################

Temperature <- read.table("PhanerozoicTemperature.txt", header = T)
res<-sm.spline(Temperature[,"Age"],Temperature[,"Temperature"],df=100)

##############################
##   Amphibian phylogenies  ##
##############################
load("FamilyAmphibiaTrees_results.Rdata"); load("finalAmphibia_results.Rdata")

pdf("Speciation RTT Amphibia.pdf", width=12, height=15)
	
	par(mfrow=c(3,3), mar=c(4,4,1,1))

for(i in 1:length(FamilyAmphibiaTrees_res))
{
  agei<-FamilyAmphibiaTrees_res[[i]]$Clade_age; print(agei)
  plot(c(seq(-agei,0,by=1),0),rep(2,length(c(seq(agei,0,by=-1),0))),xlab="",ylab="Speciation rate (events/lineage/Myr)",type="l",ylim=c(0, 0.5), col="white",lwd="1",las="1",cex.axis="1",cex.main="1",cex.lab="1",bty="n", main=(FamilyAmphibiaTrees_res[[i]][[1]]))
  abline(v=c(-145,-100.5,-66,-56,-33.9,-23.03,-5.33,-2.58),col="black",lty="dotted",lwd="1")
  abline(h=0,col="grey",lwd="0.5")
  
  legend("topleft",bty="n", c("CST", "TimeVar_EXPO", "TimeVar_LIN", "TempVar_EXPO", "TempVAR_LIN", "CarbVar_EXPO", "CarbVar_LIN") ,lwd="2",col=c("chartreuse3", "dodgerblue", "firebrick", "purple", "orange","yellow2","plum"),cex=0.7)
  legend("top",bty="n",c("Best model:",finalAmphibia[[i]][1,1],finalAmphibia[[i]][1,5]), cex = 0.7)
  
  Temperatures <- subset(Temperature, Age<agei)
  AtmCarbon<- subset(Carbon, Age<agei)
	for (j in 1:7)
	{
	model<-finalAmphibia[[i]][j,]
	model.name<-model[1]
	values<-c(as.numeric(model[-1]),1,6)
	names(values)<-names(model[-1])
	
	if(model.name =="BCST"){
		lines(c(seq(-agei,0,by=1),0),rep(abs(values["Lambda"]),length(c(seq(agei,0,by=-1),0))),type="l",main='',ylim=c(0,(max(c(abs(values["Lambda"])),abs(values["Lambda"])+0.05))), lty=1,col=alpha(color[1], 0.5),lwd="2",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
	}
	if(model.name =="BCSTDCST"){
		lines(c(seq(-agei,0,by=1),0),rep(abs(values["Lambda"]),length(c(seq(agei,0,by=-1),0))),type="l",main='',ylim=c(0,(max(abs(values["Lambda"])+0.05))), lty=2,col=alpha(color[1], 0.5),lwd="2",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
	}
	if(model.name =="BTimeVar_EXPO"){
		lines(c(seq(-agei,0,by=1),0),abs(values["Lambda"])*exp(values["Alpha"]*c(seq(agei,0,by=-1),0)),type="l",main='', lty=1, col=alpha(color[2], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
	}
	if(model.name =="BTimeVarDCST_EXPO"){
		lines(c(seq(-agei,0,by=1),0),abs(values["Lambda"])*exp(values["Alpha"]*c(seq(agei,0,by=-1),0)),type="l",main='', lty=2, col=alpha(color[2], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
	}
	if(model.name =="BCSTDTimeVar_EXPO"){
	  lines(c(seq(-agei,0,by=1),0),rep(abs(values["Lambda"]),length(c(seq(agei,0,by=-1),0))),type="l",main='',ylim=c(0,(max(c(abs(values["Lambda"])),abs(values["Lambda"])+0.05))), lty=3,col=alpha(color[2], 0.5),lwd="2",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
	}
	if(model.name =="BTimeVarDTimeVar_EXPO"){
		lines(c(seq(-agei,0,by=1),0),abs(values["Lambda"])*exp(values["Alpha"]*c(seq(agei,0,by=-1),0)),type="l",main='',lty=4, col=alpha(color[2], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
		}
	if(model.name =="BTimeVar_LIN"){
		lines(c(seq(-agei,0,by=1),0),abs(values["Lambda"])+(values["Alpha"]*c(seq(agei,0,by=-1),0)),type="l",main='', lty=1, col=alpha(color[3], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
		}
	if(model.name =="BTimeVarDCST_LIN"){
		lines(c(seq(-agei,0,by=1),0),abs(values["Lambda"])+(values["Alpha"]*c(seq(agei,0,by=-1),0)),type="l",main='', lty=2, col=alpha(color[3], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
	}
	if(model.name =="BCSTDTimeVar_LIN"){
	  lines(c(seq(-agei,0,by=1),0),rep(abs(values["Lambda"]),length(c(seq(agei,0,by=-1),0))),type="l",main='',ylim=c(0,(max(c(abs(values["Lambda"])),abs(values["Lambda"])+0.05))), lty=3,col=alpha(color[3], 0.5),lwd="2",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
	}
	if(model.name =="BTimeVarDTimeVar_LIN"){
		lines(c(seq(-agei,0,by=1),0),abs(values["Lambda"])+(values["Alpha"]*c(seq(agei,0,by=-1),0)),type="l",main='', lty=4, col=alpha(color[3], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
	}
	
	Temp_fun<-function(x){predict(res,x)}
	x=0:500
	
	if(model.name=="BTempVar_EXPO"){
	  f.lamb<-function(x){abs(values["Lambda"])*exp(values["Alpha"]*Temp_fun(x))}
	  lines(-Temperatures[,"Age"], f.lamb(Temperatures[,"Age"]),type="l",main='', lty=1, col=alpha(color[4], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n") 
	  }
	if(model.name=="BTempVarDCST_EXPO"){
	  f.lamb<-function(x){abs(values["Lambda"])*exp(values["Alpha"]*Temp_fun(x))}
	  lines(-Temperatures[,"Age"], f.lamb(Temperatures[,"Age"]),type="l",main='', lty=2, col=alpha(color[4], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")	
	  }
	if(model.name =="BCSTDTempVar_EXPO"){
	  lines(c(seq(-agei,0,by=1),0),rep(abs(values["Lambda"]),length(c(seq(agei,0,by=-1),0))),type="l",main='',ylim=c(0,(max(c(abs(values["Lambda"])),abs(values["Lambda"])+0.05))), lty=3,col=alpha(color[4], 0.5),lwd="2",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
	}
	if(model.name=="BTempVarDTempVar_EXPO"){
	  f.lamb<-function(x){abs(values["Lambda"])*exp(values["Alpha"]*Temp_fun(x))}
	  lines(-Temperatures[,"Age"], f.lamb(Temperatures[,"Age"]),type="l",main='', lty=4, col=alpha(color[4], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")	
	  }
	if(model.name=="BTempVar_LIN"){
	  f.lamb<-function(x){abs(values["Lambda"])+(values["Alpha"]*Temp_fun(x))}
	  lines(-Temperatures[,"Age"], f.lamb(Temperatures[,"Age"]),type="l",main='', lty=1, col=alpha(color[5], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")	
	  }
	if(model.name=="BTempVarDCST_LIN"){
	  f.lamb<-function(x){abs(values["Lambda"])+(values["Alpha"]*Temp_fun(x))}
	  lines(-Temperatures[,"Age"], f.lamb(Temperatures[,"Age"]),type="l",main='', lty=2, col=alpha(color[5], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")	
	  }
	if(model.name=="BCSTDTempVar_LIN"){
	  lines(c(seq(-agei,0,by=1),0),rep(abs(values["Lambda"]),length(c(seq(agei,0,by=-1),0))),type="l",main='',ylim=c(0,(max(c(abs(values["Lambda"])),abs(values["Lambda"])+0.05))), lty=3,col=alpha(color[5], 0.5),lwd="2",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")
	}
	if(model.name=="BTempVarDTempVar_LIN"){
	  f.lamb<-function(x){abs(values["Lambda"])+(values["Alpha"]*Temp_fun(x))}
	  lines(-Temperatures[,"Age"], f.lamb(Temperatures[,"Age"]),type="l",main='', lty=4, col=alpha(color[5], 0.5),lwd="1",las="1",cex.axis="0.7",cex.main="0.7",cex.lab="0.7",bty="n")	
	}
	
	}
}

dev.off()
