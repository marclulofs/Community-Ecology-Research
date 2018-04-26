setwd("~/Documents/Stage_M2/codes+data")

### Load all results
Tetrapod=c("Amphibia","Bird","CrocoTurtle","Mammal","Squamate")
for(i in 1:length(Tetrapod)){load(sprintf('final%s_results.Rdata',Tetrapod[i]))}

finalTetrapod=list(finalAmphibia,finalSquamate,finalCrocoTurtle,finalBird,finalMammal)

################################################################
####  What are the models selected across phylogenies ?  #######
################################################################

# Blank dataset
Best_models <- matrix(NA,7,6)
colnames(Best_models) <- c("Tetrapoda","Amphibia","Squamata","Testudines+Croc","Aves","Mammalia")
rownames(Best_models)<- c("BCST",#"BCSTDCST",
                    "BTimeVar_EXPO",#"BTimeVarDCST_EXPO","BCSTDTimeVar_EXPO","BTimeVarDTimeVar_EXPO",
                    "BTimeVar_LIN",#"BTimeVarDCST_LIN","BCSTDTimeVar_LIN","BTimeVarDTimeVar_LIN",
                    "BTempVar_EXPO",#"BTempVarDCST_EXPO","BCSTDTempVar_EXPO","BTempVarDTempVar_EXPO",
                    "BTempVar_LIN",#"BTempVarDCST_LIN","BCSTDTempVar_LIN","BTempVarDTempVar_LIN",
                    "BCarbVar_EXPO",#"BCarbVarDCST_EXPO","BCSTDCarbVar_EXPO","BCarbVarDCarbVar_EXPO",
                    "BCarbVar_LIN"#"BCarbVarDCST_LIN","BCSTDCarbVar_LIN","BCarbVarDCarbVar_LIN"
                    )
Best_models[,1:6]<-0

# Count how many times each model is selected for each group
for (i in 1:length(finalTetrapod)){
  for(j in 1:length(finalTetrapod[[i]])){
    model.name <- finalTetrapod[[i]][[j]][1]
    
    if(model.name =="BCST"){
      Best_models[1,i+1] <- as.numeric(Best_models[1,i+1])+1
    }
    if(model.name =="BTimeVar_EXPO"){
      Best_models[2,i+1] <- as.numeric(Best_models[2,i+1])+1
    }
    if(model.name =="BTimeVar_LIN"){
      Best_models[3,i+1] <- as.numeric(Best_models[3,i+1])+1
    }
    if(model.name=="BTempVar_EXPO"){
      Best_models[4,i+1] <- as.numeric(Best_models[4,i+1])+1
    }
    if(model.name=="BTempVar_LIN"){
      Best_models[5,i+1] <- as.numeric(Best_models[5,i+1])+1
    }
    if(model.name=="BCarbVar_EXPO"){
      Best_models[6,i+1] <- as.numeric(Best_models[6,i+1])+1
    }
    if(model.name=="BCarbVar_LIN"){
      Best_models[7,i+1] <- as.numeric(Best_models[7,i+1])+1
    }
  }
}

# Sum for all Tetrapods
for(j in 1:7){Best_models[j,1] <- sum(as.numeric(Best_models[j,2:6]))}

#write.table(Best_models,file="Best_models.txt",quote=FALSE,sep="\t",row.names=TRUE)
#Best_models <- read.table("Best_models.txt", header=T)

##  Conversion in percentage
bm_prop <- Best_models
for (i in 1:7){
  for(j in 1:6){
    total_phyl <- sum(as.numeric(Best_models[1:7,j]))
    prop<-as.numeric(Best_models[i,j])/total_phyl
    bm_prop[i,j]<-prop
  }
  print(sum(bm_prop[1:7,j]))
}
#write.table(bm_prop,file = "bm_prop.txt",quote=FALSE,sep="\t",row.names=TRUE)

################
### Figure 5 ###
################
#pdf(file = "Fig5 barplot.pdf")
barplot(bm_prop,legend.text = T,beside=F,border=NA,horiz = F,space=0.2,col=c("azure2", "deepskyblue2", "blue4", "sienna1", "red3","yellowgreen","green4"),xaxt="n", main = "")
#dev.off()

#for (i in 1:7){for(j in 1:6){bm_prop[i,j]<-0}}
