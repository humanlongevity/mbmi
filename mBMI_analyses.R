library(readxl)
library(reshape2)
library(survival)
library(powerSurvEpi)
library(lattice)
library(logistf)
library(alluvial)

##
##functions
##
regress_multivar_o<-function(pertinentdataframe,ncovar,noutcomes,choice){
  line<-0
  D<-pertinentdataframe
  if(ncovar==0){
    out <- data.frame(variable = NA,outcome=NA, p= NA,coef=NA,r2=NA)   
    if(choice=="linear"){
      for (i in (ncovar+1+noutcomes):ncol(D)) { # features
        for (y in 1:noutcomes){
          if(!is.na(mean(na.omit(D[,i])))){
            line<-line+1
            m <- summary(lm(D[,c(y,i)])) # run model without adjustment 
            out[line, 1] <- colnames(D)[i]           # print variable name
            out[line, 2] <- colnames(D)[y]           # print outcome name
            out[line, 3] <- m$coefficients[2,4]      # print p value
            out[line, 4] <- m$coefficients[2,1]      # print coefficient
            out[line, 5] <- m$r.squared              # print r2
          }
        }
      }
    }
    else {
      for (i in (ncovar+1+noutcomes):ncol(D)) { # features
        for (y in 1:noutcomes){
          if(!is.na(mean(na.omit(D[,i])))){
            line<-line+1
            m <- summary(glm(D[,c(y,i)],family="binomial")) # run model without adjustment 
            out[line, 1] <- colnames(D)[i]           # print variable name
            out[line, 2] <- colnames(D)[y]           # print outcome name
            out[line, 3] <- m$coefficients[2,4]      # print p value
            out[line, 4] <- m$coefficients[2,1]      # print coefficient
          }
        }
      }
    }
    
    
  }
  else{
    out <- data.frame(variable = NA,outcome=NA, p= NA,coef=NA,r2=NA,p_withcovar= NA,coef_withcovar=NA)   
    if(choice=="linear"){
      for (i in (ncovar+1+noutcomes):ncol(D)) { # features
        for (y in 1:noutcomes){
          if(!is.na(mean(na.omit(D[,i])))){
            line<-line+1
            m <- summary(lm(D[,c(y,i)])) # run model without adjustment 
            out[line, 1] <- colnames(D)[i]           # print variable name
            out[line, 2] <- colnames(D)[y]           # print outcome name
            out[line, 3] <- m$coefficients[2,4]      # print p value
            out[line, 4] <- m$coefficients[2,1]      # print coefficient
            out[line, 5] <- m$r.squared              # print r2
            m <- summary(lm(D[,c(y,i,(1+noutcomes):(ncovar+noutcomes))])) # run model with adjustment 
            out[line, 6] <- m$coefficients[2,4]      # print p value
            out[line, 7] <- m$coefficients[2,1]      # print coefficient
          }
        }
      }
    }
    else {
      for (i in (ncovar+1+noutcomes):ncol(D)) { # features
        for (y in 1:noutcomes){
          if(!is.na(mean(na.omit(D[,i])))){
            line<-line+1
            m <- summary(glm(D[,c(y,i)],family="binomial")) # run model without adjustment 
            out[line, 1] <- colnames(D)[i]           # print variable name
            out[line, 2] <- colnames(D)[y]           # print outcome name
            out[line, 3] <- m$coefficients[2,4]      # print p value
            out[line, 4] <- m$coefficients[2,1]      # print coefficient
            m <- summary(glm(D[,c(y,i,(1+noutcomes):(ncovar+noutcomes))],family="binomial")) # run model with adjustment 
            out[line, 6] <- m$coefficients[2,4]      # print p value
            out[line, 7] <- m$coefficients[2,1]      # print coefficient
          }
        }
      }
    }
  }
  return(out)
}



##
##analyses
##

##
##read in the supplementary data:
df<-as.data.frame(read_excel("table s2.xlsx"))
#set factor levels
df$Category<-factor(df$Category,levels=c("Normal weight, metabolically healthy","Overweight, metabolically overweight","Obese, metabolically obese","Outlier: Metabolic BMI << BMI","Outlier: Metabolic BMI >> BMI"))
df$PGCategory<-factor(df$PGCategory,levels=c("PG 0-10%","PG 10-25%","PG 25-50%","PG 50-75%","PG 75-90%","PG 90-100%","MC4R carrier"))
df$OvObCategory<-factor(df$OvObCategory,levels=c("Normal weight, metabolically healthy","Overweight, metabolically overweight","Obese, metabolically obese","Obese, metabolically healthy","Obese, metabolically overweight","Overweight, metabolically healthy","Normal weight, metabolically obese","Normal weight, metabolically overweight","Overweight, metabolically obese"))


##I could see demonstrating how the outlier categories and residuals are calculated here in the code
#or maybe just check them

##
##compare mBMI to BMI
##

##
##Plotting mBMI/BMI groups versus various phenotypes:
##

#restrict to unrelated
unrelated_eur_df<-df[!is.na(df$PC1),]

#remove two extreme TG outliers to aid graphing:
unrelated_eur_df1<-unrelated_eur_df
unrelated_eur_df1$TG[unrelated_eur_df1$TG>6]<-NA

#put everything on the same scale to aid graphing:
scaled_unrelated_eur_df<-unrelated_eur_df1
for (i in c(6:15,17:27,29:31,33:37)){
  scaled_unrelated_eur_df[,i]<-(unrelated_eur_df1[,i]-min(na.omit(unrelated_eur_df1[,i])))/(max(na.omit(unrelated_eur_df1[,i]))-min(na.omit(unrelated_eur_df1[,i])))
}

#plot Figure 2B:
temp3a1<-subset(scaled_unrelated_eur_df,select=c(Category,BMI,Age,IR,WH,SYSBP,DIABP,PG))
temp3a1<-temp3a1[!is.na(temp3a1$Category),]
temp4a<-melt(data=temp3a1,id=1)
p <- ggplot(temp4a, aes(x=variable, y=value, fill=Category)) + 
  geom_boxplot(outlier.size=-1,cex=1)
p +scale_fill_manual(values=c("cornflowerblue", "orange", "red","cyan","magenta"))+theme(text = element_text(size=30),legend.title=element_blank())

#plot Figure 2C:
temp3a1mri_tg<-subset(scaled_unrelated_eur_df,select=c(Category,TG,Chol,LDL,HDL,AG,PFAT,VAT,SAT))
temp3a1mri_tg<-temp3a1mri_tg[!is.na(temp3a1mri_tg$Category),]
temp4amri_tg<-melt(data=temp3a1mri_tg,id=1)
p <- ggplot(temp4amri_tg, aes(x=variable, y=value, fill=Category)) + 
  geom_boxplot(outlier.size=-1,cex=1)
p +scale_fill_manual(values=c("cornflowerblue","orange","red","cyan","magenta"))+theme(text = element_text(size=30),legend.title=element_blank())


#and split up by genetic category:

#plot Figure 5A
temp3a1g<-subset(scaled_unrelated_eur_df,select=c(PGCategory,BMI,Age,IR,WH,SYSBP,DIABP,PG))
temp3a1g<-temp3a1g[!is.na(temp3a1g$PGCategory),]
temp4ag<-melt(data=temp3a1g,id=1)
p <- ggplot(temp4ag, aes(x=variable, y=value, fill=PGCategory)) + 
  geom_boxplot(outlier.size=-1,cex=1)
p + scale_fill_manual(values=c("dark gray","purple","blue","olivedrab3","yellow","orange","red"))+theme(legend.key.size = unit(2,"line"),legend.text=element_text(size=30),text = element_text(size=30),legend.title=element_blank())

#plot Figure 5B
temp3a1mrig<-subset(scaled_unrelated_eur_df,select=c(PGCategory,TG,Chol,LDL,HDL,AG,PFAT,VAT,SAT))
temp3a1mrig<-temp3a1mrig[!is.na(temp3a1mrig$PGCategory),]
temp4amrig<-melt(data=temp3a1mrig,id=1)
p <- ggplot(temp4amrig_tg, aes(x=variable, y=value, fill=PGCategory)) + 
  geom_boxplot(outlier.size=-1,cex=1)
p + scale_fill_manual(values=c("dark gray","purple","blue","olivedrab3","yellow","orange","red"))+theme(legend.key.size = unit(2,"line"),legend.text=element_text(size=30),text = element_text(size=30),legend.title=element_blank())

#and split them up by mBMI/BMI category when it's broken into obesity and overweight predictions as opposed to all outliers

#plot Figure S3A
temp3a1ovob<-subset(scaled_unrelated_eur_df,select=c(OvObCategory,BMI,Age,IR,WH,SYSBP,DIABP,PG))
temp3a1ovob<-temp3a1ovob[!is.na(temp3a1ovob$OvObCategory),]
temp4aovob<-melt(data=temp3a1ovob,id=1)
p <- ggplot(temp4aovob, aes(x=variable, y=value, fill=OvObCategory)) + 
  geom_boxplot(outlier.size=-1,cex=1)
p + geom_jitter(shape=20,position=position_jitterdodge(0.20)) +scale_fill_manual(values=c("light gray","orange","red","pink1","pink2","pink3","yellow1","yellow2","yellow3"))+theme(text = element_text(size=30),legend.title=element_blank())

#plot Figure S3B
temp3a1ovobmri<-subset(scaled_unrelated_eur_df,select=c(OvObCategory,TG,Chol,LDL,HDL,AG,PFAT,VAT,SAT))
temp3a1ovobmri<-temp3a1ovobmri[!is.na(temp3a1ovobmri$OvObCategory),]
temp4aovobmri<-melt(data=temp3a1ovobmri,id=1)
p <- ggplot(temp4aovobmri, aes(x=variable, y=value, fill=OvObCategory)) + 
  geom_boxplot(outlier.size=-1,cex=1)
p + geom_jitter(shape=20,position=position_jitterdodge(0.20)) +scale_fill_manual(values=c("light gray","orange","red","pink1","pink2","pink3","yellow1","yellow2","yellow3"))+theme(text = element_text(size=30),legend.title=element_blank())


##
##calculate p-values for differences between outlier groups
##

temp5a1<-subset(scaled_unrelated_eur_df,select=c(Category,BMI,Age,IR,WH,SYSBP,DIABP,PG,original_PG,TG,Chol,LDL,HDL,AG,PFAT,VAT,SAT,HighBloodPressureMed_v1))
temp5a1$mbmi_below_bmi_vs_normal[temp5a1$Category=="Outlier: Metabolic BMI << BMI"]<-1
temp5a1$mbmi_below_bmi_vs_normal[temp5a1$Category=="Normal weight, metabolically healthy"]<-0
temp5a1$mbmi_above_bmi_vs_normal[temp5a1$Category=="Outlier: Metabolic BMI >> BMI"]<-1
temp5a1$mbmi_above_bmi_vs_normal[temp5a1$Category=="Normal weight, metabolically healthy"]<-0
temp5a1$mbmi_below_bmi_vs_ov[temp5a1$Category=="Outlier: Metabolic BMI << BMI"]<-1
temp5a1$mbmi_below_bmi_vs_ov[temp5a1$Category=="Overweight, metabolically overweight"]<-0
temp5a1$mbmi_above_bmi_vs_ov[temp5a1$Category=="Outlier: Metabolic BMI >> BMI"]<-1
temp5a1$mbmi_above_bmi_vs_ov[temp5a1$Category=="Overweight, metabolically overweight"]<-0
temp5a1$mbmi_below_bmi_vs_ob[temp5a1$Category=="Outlier: Metabolic BMI << BMI"]<-1
temp5a1$mbmi_below_bmi_vs_ob[temp5a1$Category=="Obese, metabolically obese"]<-0
temp5a1$mbmi_above_bmi_vs_ob[temp5a1$Category=="Outlier: Metabolic BMI >> BMI"]<-1
temp5a1$mbmi_above_bmi_vs_ob[temp5a1$Category=="Obese, metabolically obese"]<-0
temp5a1$mbmi_above_bmi_vs_mbmi_below[temp5a1$Category=="Outlier: Metabolic BMI >> BMI"]<-1
temp5a1$mbmi_above_bmi_vs_mbmi_below[temp5a1$Category=="Outlier: Metabolic BMI << BMI"]<-0
pert<-temp5a1[,c(19:25,2:18)]
temp5a1_out<-regress_multivar_o(pert,0,7,"logistic")
#temp5a1_out contains the p-values for all comparisons between groups

##
##check groups are as expected
##
df$check_group<-NA
df$check_group[df$BMI<24.95&df$BMI>18.45&df$Residuals<.5&df$Residuals>(-.5)]<-"Normal weight, metabolically healthy"
df$check_group[df$Residuals>.5]<-"Outlier: Metabolic BMI >> BMI"
df$check_group[df$Residuals<(-.5)]<-"Outlier: Metabolic BMI << BMI"
#cutoffs are chosen to match standard BMI guidelines but allow for reasonable rounding
df$check_group[df$BMI>24.95&df$Residuals<.5&df$Residuals>(-.5)]<-"Overweight, metabolically overweight"
df$check_group[df$BMI>30&df$Residuals<.5&df$Residuals>(-.5)]<-"Obese, metabolically obese"
table(df$check_group,df$Category)

df$check_ovob_group<-NA
df$check_ovob_group[df$mBMI<(-0.073)&df$BMI>18.45&df$BMI<24.95]<-"Normal weight, metabolically healthy"
df$check_ovob_group[df$mBMI>(-0.073)&df$mBMI<(0.314)&df$BMI>18.45&df$BMI<24.95]<-"Normal weight, metabolically overweight"
df$check_ovob_group[df$mBMI>(0.314)&df$BMI>18.45&df$BMI<24.95]<-"Normal weight, metabolically obese"
df$check_ovob_group[df$mBMI<(-0.073)&df$BMI>24.95&df$BMI<29.95]<-"Overweight, metabolically healthy"
df$check_ovob_group[df$mBMI>(-0.073)&df$mBMI<(0.314)&df$BMI>24.95&df$BMI<29.95]<-"Overweight, metabolically overweight"
df$check_ovob_group[df$mBMI>(0.314)&df$BMI>24.95&df$BMI<29.95]<-"Overweight, metabolically obese"
df$check_ovob_group[df$mBMI<(-0.073)&df$BMI>29.95]<-"Obese, metabolically healthy"
df$check_ovob_group[df$mBMI>(-0.073)&df$mBMI<(0.314)&df$BMI>29.95]<-"Obese, metabolically overweight"
df$check_ovob_group[df$mBMI>(0.314)&df$BMI>29.95]<-"Obese, metabolically obese"
table(df$check_ovob_group,df$OvObCategory)

##
##see if different groups are more likely to become obese
##
df$BMI_Category_v3<-NA
df$BMI_Category_v3[df$BMIv3<24.95&df$BMIv3>18.45]<-"Normal weight"
df$BMI_Category_v3[df$BMIv3>24.95&df$BMIv3<29.95]<-"Overweight"
df$BMI_Category_v3[df$BMIv3>29.95]<-"Obese"
df$BMI_Category_v3<-factor(df$BMI_Category_v3,levels=c("Normal weight","Overweight","Obese"))
#only looking at people who had all timepoints measured:
table(df$Category[!is.na(df$BMIv2)],df$BMI_Category_v3[!is.na(df$BMIv2)])
table(df$OvObCategory[!is.na(df$BMIv2)],df$BMI_Category_v3[!is.na(df$BMIv2)])

#and make a plot like Figure S5
melted<-melt(table(df$Category[!is.na(df$BMIv2)],df$BMI_Category_v3[!is.na(df$BMIv2)]))
p <- ggplot(melted, aes(x=Var1, y=value, fill=Var2)) + 
  geom_col(cex=1,position="dodge")
p +scale_fill_manual(values=c("cornflowerblue", "orange", "red"))+theme(text = element_text(size=30),legend.title=element_blank())
melted<-melt(table(df$OvObCategory[!is.na(df$BMIv2)],df$BMI_Category_v3[!is.na(df$BMIv2)]))
p <- ggplot(melted, aes(x=Var1, y=value, fill=Var2)) + 
  geom_col(cex=1,position="dodge")
p +scale_fill_manual(values=c("cornflowerblue", "orange", "red"))+theme(text = element_text(size=30),legend.title=element_blank())

#see if people who are normal weight but have an overweight or obese BMI at timepoint 1 are more likely to switch to overweight or obese at timepoint 3 than are those with normal weight and a healthy metabolome at v1:
df$ov_or_ob_at_v3[df$BMI_Category_v3=="Normal weight"]<-0
df$ov_or_ob_at_v3[df$BMI_Category_v3=="Overweight"]<-1
df$ov_or_ob_at_v3[df$BMI_Category_v3=="Obese"]<-1
df$normal_weight_ov_ob_mbmi_v1[df$OvObCategory=="Normal weight, metabolically healthy"]<-0
df$normal_weight_ov_ob_mbmi_v1[df$OvObCategory=="Normal weight, metabolically overweight"]<-1
df$normal_weight_ov_ob_mbmi_v1[df$OvObCategory=="Normal weight, metabolically obese"]<-1
table(df$normal_weight_ov_ob_mbmi_v1,df$ov_or_ob_at_v3)
fisher.test(df$normal_weight_ov_ob_mbmi_v1,df$ov_or_ob_at_v3)

#count proportion of obese metabolome individuals at v1 who still have obese metabolome at v3:
length(na.omit(df$mBMIv3[df$mBMI>0.314&df$mBMIv3>0.314]))/length(na.omit(df$mBMIv3[df$mBMI>0.314]))

#make Figure 4a:
df$BMI_Category_v1<-NA
df$BMI_Category_v1[df$BMI<24.95&df$BMI>18.45]<-"Normal weight"
df$BMI_Category_v1[df$BMI>24.95&df$BMI<29.95]<-"Overweight"
df$BMI_Category_v1[df$BMI>29.95]<-"Obese"
df$BMI_Category_v1<-factor(df$BMI_Category_v1,levels=c("Normal weight","Overweight","Obese"))

df$mBMI_Category_v1<-NA
df$mBMI_Category_v1[df$mBMI<(-0.073)]<-"Metabolically healthy"
df$mBMI_Category_v1[df$mBMI>(-0.073)&df$mBMI<0.314]<-"Metabolically overweight"
df$mBMI_Category_v1[df$mBMI>0.314]<-"Metabolically obese"
df$mBMI_Category_v1<-factor(df$mBMI_Category_v1,levels=c("Metabolically healthy","Metabolically overweight","Metabolically obese"))

df$mBMI_Category_v3<-NA
df$mBMI_Category_v3[df$mBMIv3<(-0.073)]<-"Metabolically healthy"
df$mBMI_Category_v3[df$mBMIv3>(-0.073)&df$mBMIv3<0.314]<-"Metabolically overweight"
df$mBMI_Category_v3[df$mBMIv3>0.314]<-"Metabolically obese"
df$mBMI_Category_v3<-factor(df$mBMI_Category_v3,levels=c("Metabolically healthy","Metabolically overweight","Metabolically obese"))

alluvial1<-subset(df[!is.na(df$BMI)&!is.na(df$BMIv2)&!is.na(df$BMIv3),],select=c(BMI_Category_v1,BMI_Category_v3,mBMI_Category_v1,mBMI_Category_v3))
alluvial1<-na.omit(alluvial1)

one<-data.frame(table(paste(alluvial1$mBMI_Category_v1,alluvial1$BMI_Category_v1,alluvial1$BMI_Category_v3,sep="_")))
onedf <- data.frame(do.call('rbind', strsplit(as.character(one$Var1),'_',fixed=TRUE)))

colnames(onedf)[1:3]<-c('Timepoint_1_mBMI','Timepoint_1_BMI','Timepoint_3_BMI')

foralluvialfreq<-data.frame(onedf,one$Freq)

colnames(foralluvialfreq)[4]<-"Freq"

foralluvialfreq$Timepoint_1_BMI<-factor(foralluvialfreq$Timepoint_1_BMI,levels=c("Normal weight","Overweight","Obese"))
foralluvialfreq$Timepoint_3_BMI<-factor(foralluvialfreq$Timepoint_3_BMI,levels=c("Normal weight","Overweight","Obese"))
foralluvialfreq$Timepoint_1_mBMI<-factor(foralluvialfreq$Timepoint_1_mBMI,levels=c("Metabolically healthy","Metabolically overweight","Metabolically obese"))

alluvial(foralluvialfreq[,c(2,3)], freq=foralluvialfreq$Freq,
         col = ifelse(foralluvialfreq$Timepoint_1_mBMI == "Metabolically obese", "red",ifelse(foralluvialfreq$Timepoint_1_mBMI == "Metabolically overweight", "orange", "cornflowerblue")),
         border = ifelse(foralluvialfreq$Timepoint_1_mBMI == "Metabolically obese", "red",ifelse(foralluvialfreq$Timepoint_1_mBMI == "Metabolically overweight", "orange", "cornflowerblue")),
         cex = 1.2
)

#and 4b:
one<-data.frame(table(paste(alluvial1$BMI_Category_v1,alluvial1$mBMI_Category_v1,alluvial1$mBMI_Category_v3,sep="_")))
onedf <- data.frame(do.call('rbind', strsplit(as.character(one$Var1),'_',fixed=TRUE)))

colnames(onedf)[1:3]<-c('Timepoint_1_BMI','Timepoint_1_mBMI','Timepoint_3_mBMI')

foralluvialfreq<-data.frame(onedf,one$Freq)

colnames(foralluvialfreq)[4]<-"Freq"

foralluvialfreq2$Timepoint_1_BMI<-factor(foralluvialfreq2$Timepoint_1_BMI,levels=c("Normal weight","Overweight","Obese"))
foralluvialfreq2$Timepoint_1_mBMI<-factor(foralluvialfreq2$Timepoint_1_mBMI,levels=c("Metabolically healthy","Metabolically overweight","Metabolically obese"))
foralluvialfreq2$Timepoint_3_mBMI<-factor(foralluvialfreq2$Timepoint_3_mBMI,levels=c("Metabolically healthy","Metabolically overweight","Metabolically obese"))

alluvial(foralluvialfreq2[,c(2,3)], freq=foralluvialfreq2$Count.of.count,
         col = ifelse(foralluvialfreq2$Timepoint_1_BMI == "Obese", "red",ifelse(foralluvialfreq2$Timepoint_1_BMI == "Overweight", "orange", "cornflowerblue")),
         border = ifelse(foralluvialfreq2$Timepoint_1_BMI == "Obese", "red",ifelse(foralluvialfreq2$Timepoint_1_BMI == "Overweight", "orange", "cornflowerblue")),
         cex = 0.7
)


##
##cardiovascular events
##

#age of twinsuk participants at each visit:
summary(df$Age[df$cohort=="TwinsUK"])
summary(df$Agev3[df$cohort=="TwinsUK"])

#get the counts for each type of event before v1:
table(df$Category,df$Stroke_before_v1)
table(df$Category,df$Infarction_before_v1)
table(df$Category,df$Angina_before_v1)
table(df$Category,df$Angioplasty_before_v1)
table(df$Category,df$Any_CV_or_stroke_before_v1)
table(df$Category,df$Any_CV_before_v1)
fisher.test(df$Category,df$Any_CV_before_v1)
table(df$Category,df$HighBloodPressureMed_v1)
fisher.test(df$Category[df$Category=="Outlier: Metabolic BMI << BMI"|df$Category=="Outlier: Metabolic BMI >> BMI"],df$HighBloodPressureMed_v1[df$Category=="Outlier: Metabolic BMI << BMI"|df$Category=="Outlier: Metabolic BMI >> BMI"])

#and after v1:
table(df$Category,df$Stroke_after_v1)
table(df$Category,df$Infarction_after_v1)
table(df$Category,df$Angina_after_v1)
table(df$Category,df$Angioplasty_after_v1)
table(df$Category,df$Any_CV_after_v1)
table(df$Category,df$Any_CV_or_stroke_after_v1)
table(df$OvObCategory,df$Any_CV_after_v1)
#event rate in people with healthy metabolomes:
length(na.omit(df$Any_CV_or_stroke_after_v1[df$Any_CV_or_stroke_after_v1==1&df$mBMI<(-.073)]))/length(na.omit(df$Any_CV_or_stroke_after_v1[df$mBMI<(-.073)]))
#and rate in people with unhealthy metabolomes:
 #who are normal weight/overweight
length(na.omit(df$Any_CV_or_stroke_after_v1[df$Any_CV_after_v1==1&df$mBMI>0.314&df$BMI<29.95]))/length(na.omit(df$Any_CV_after_v1[df$mBMI>0.314&df$BMI<29.95]))
 #and who are obese:
length(na.omit(df$Any_CV_or_stroke_after_v1[df$Any_CV_after_v1==1&df$mBMI>0.314&df$BMI>29.95]))/length(na.omit(df$Any_CV_after_v1[df$mBMI>0.314&df$BMI>29.95]))

##power check for survival analysis
 #for 5 groups
 #frequency of each group
add5<-c(0.38,0.29,0.12,0.10,0.10)
 #frequency of event in each group
add15<-c(0.010,0.036,0.030,0.000,0.037)
temp5<-c(1,1,1,1,1)
power.stratify(1750,13,gVec=add5,PVec=add15,lambda0Vec = temp5,HR=1.5)

 #and for 2 outlier groups
 #frequency of each group
add<-c(0.45,0.55)
 #frequency of event in each group
add1<-c(0.0140,0.0287)
temp<-c(1,1)
power.stratify(1750,13,gVec=add,PVec=add1,lambda0Vec = temp,HR=1.5)

#survival analysis
 #first all 5 categories, with age as a covariate
summary(coxph(Surv(df$Age_at_any_CV_or_censor-df$Age, df$Any_CV_after_v1,type='right') ~ df$Category+df$Age))
 #then just the two healthier groups vs. the two less healthy groups
df$twocat[df$Category=="Normal weight, metabolically healthy"|df$Category=="Outlier: Metabolic BMI << BMI"]<-"Healthier mBMI"
df$twocat[df$Category=="Obese, metabolically obese"|df$Category=="Overweight, metabolically overweight"|df$Category=="Outlier: Metabolic BMI >> BMI"]<-"Less healthy mBMI"
summary(coxph(Surv(df$Age_at_any_CV_or_censor-df$Age, df$Any_CV_after_v1,type='right') ~ df$twocat+df$Age))


#and make a graph like in Figure 4C:
melted<-melt(df[,c(3,39,40,42,43,46,51)],id=1)
melted<-na.omit(melted)
p <- ggplot(melted, aes(x=variable, y=value, fill=Category)) + 
  geom_bar(cex=1,stat="identity")
p +scale_fill_manual(values=c("cornflowerblue", "orange", "red","cyan","magenta"))+theme(text = element_text(size=30),legend.title=element_blank())

#and like in Figure S6
melted<-melt(df[,c(3,34,35,37,38,44)],id=1)
melted<-na.omit(melted)
p <- ggplot(melted, aes(x=variable, y=value, fill=Category)) + 
  geom_bar(cex=1,stat="identity")
p +scale_fill_manual(values=c("cornflowerblue", "orange", "red","cyan","magenta"))+theme(text = element_text(size=30),legend.title=element_blank())

#and make Figure 4D
km_agewithstartage_anycv_result <- survfit(Surv(df$Age,df$Age_at_any_CV_or_censor, df$Any_CV_after_v1) ~ df$Category, type="kaplan-meier", conf.type="log")
anycvlegend<-c("Normal weight, metabolically healthy","Overweight, metabolically overweight","Obese, metabolically obese","Outlier: Metabolic BMI << BMI","Outlier: Metabolic BMI >> BMI")
plot(km_agewithstartage_anycv_result,
     xlab="Age", ylab="Proportion without Any CV Event", lwd=4, col=c("cornflowerblue","red","cyan","magenta","orange"),ylim=c(0.7,1),xlim=c(40,90),cex.axis=1.5,cex.lab=1.5,cex=2)
legend(x="bottomleft", col=c("cornflowerblue","orange","red","cyan","magenta"), lwd=4, legend=anycvlegend,cex=1.5)

##
##compare twins to each other
##

#create a dataframe of twins
pairdf<-df[!is.na(df$familyID),]
#sort by family ID and BMI
pairdf<-pairdf[order(pairdf$familyID,pairdf$BMI),]
#separate out pairs of twins to merge them again
pairdf$twinnumber<-c(1,2)
pairdf1<-pairdf[pairdf$twinnumber==1,]
pairdf2<-pairdf[pairdf$twinnumber==2,]
pairdf12<-data.frame(pairdf1,pairdf2)

#now put twins into categories:
pairdf12$twin_weight_pattern<-NA
pairdf12$twin_weight_pattern[pairdf12$BMI.1>18.45&pairdf12$BMI>18.45]<-"Both normal weight"
pairdf12$twin_weight_pattern[pairdf12$BMI.1>24.95&pairdf12$BMI>24.95]<-"Both overweight"
pairdf12$twin_weight_pattern[pairdf12$BMI.1>29.95&pairdf12$BMI>29.95]<-"Both obese"
pairdf12$twin_weight_pattern[pairdf12$BMI.1>29.95&pairdf12$BMI>18.45&pairdf12$BMI<24.95]<-"Y normal weight, X obese"
pairdf12$twin_weight_pattern[pairdf12$BMI.1>24.95&pairdf12$BMI.1<29.95&pairdf12$BMI>18.45&pairdf12$BMI<24.95]<-"Y normal weight, X overweight"
pairdf12$twin_weight_pattern[pairdf12$BMI.1>29.95&pairdf12$BMI<29.95&pairdf12$BMI>24.95]<-"Y overweight, X obese"
#there are some NA instances here b/c the lower weight twin is below 18.5; remove them
pairdf12<-pairdf12[!is.na(pairdf12$twin_weight_pattern),]
#and make it a factor
pairdf12$twin_weight_pattern<-factor(pairdf12$twin_weight_pattern,levels=c("Both normal weight","Both overweight","Both obese","Y normal weight, X obese","Y normal weight, X overweight","Y overweight, X obese"))

#make separate MZ and DZ categories
justmz<-pairdf12[pairdf12$Zygosity=="MZ",]
justdz<-pairdf12[pairdf12$Zygosity=="DZ",]

#and remove the pairs where either twin is overweight or mBMI overweight
justmz_simple<-justmz[justmz$twin_weight_pattern!="Y overweight, X obese",]
justmz_simple<-justmz_simple[justmz_simple$twin_weight_pattern!="Y normal weight, X overweight",]

justdz_simple<-justdz[justdz$twin_weight_pattern!="Y overweight, X obese",]
justdz_simple<-justdz_simple[justdz_simple$twin_weight_pattern!="Y normal weight, X overweight",]

justmz_simpler<-justmz_simple[justmz_simple$twin_weight_pattern!="Both overweight",]
justdz_simpler<-justdz_simple[justdz_simple$twin_weight_pattern!="Both overweight",]

#and make figure s2:
col<-c("light gray","dim gray","black","spring green","deep sky blue","medium blue")
col2<-c("light gray","spring green","black")
xyplot(justmz_simpler$mBMI~justmz_simpler$mBMI.1, groups = justmz_simpler$twin_weight_pattern, xlab=list(label="X Twin Prediction",cex=2),ylab=list(label="Y Twin Prediction",cex=2),
       key = list(text = list(as.character(unique(justmz_simpler$twin_weight_pattern))), 
                  points = list(pch = 19, col = col2),cex=2,corner = c(.05, .93)), pch = 19, col = col,cex=2,scales=list(cex=2))

col<-c("light gray","dim gray","black","spring green","deep sky blue","medium blue")
col2<-c("light gray","black","spring green")
xyplot(justdz_simpler$mBMI~justdz_simpler$mBMI.1, groups = justdz_simpler$twin_weight_pattern, xlab=list(label="X Twin Prediction",cex=2),ylab=list(label="Y Twin Prediction",cex=2),
       key = list(text = list(as.character(unique(justdz_simpler$twin_weight_pattern))), 
                  points = list(pch = 19, col = col2),cex=2,corner = c(.05, .93)), pch = 19, col = col,cex=2,scales=list(cex=2.1),ylim=c(-1.3,2.1),xlim=c(-1.4,2.1))

#get the r2:
summary(lm(justmz_simpler$mBMI~justmz_simpler$mBMI.1))
summary(lm(justdz_simpler$mBMI~justdz_simpler$mBMI.1))

##
##genetics
##

#polygenic risk score vs. BMI
#here we use the original polygenic risk score as opposed to the percentile
summary(lm(unrelated_eur_df$normal_BMI[unrelated_eur_df$cohort=="TwinsUK"]~unrelated_eur_df$original_PG[unrelated_eur_df$cohort=="TwinsUK"]))
summary(lm(unrelated_eur_df$normal_BMI[unrelated_eur_df$cohort=="HN"]~unrelated_eur_df$original_PG[unrelated_eur_df$cohort=="HN"]))
summary(lm(unrelated_eur_df$normal_BMIv2~unrelated_eur_df$original_PG))
summary(lm(unrelated_eur_df$normal_BMIv3~unrelated_eur_df$original_PG))

#make Figure S7A
tempgraph<-subset(unrelated_eur_df,select=c(original_PG,BMI,BMIv2,BMIv3))
tempgraph<-tempgraph[unrelated_eur_df$cohort=="TwinsUK",]
p <- ggplot(data = tempgraph)
for (i in 2:4){
  cat=as.character(i)
  loop_input<-paste("geom_smooth(se=FALSE,color=cat,aes(tempgraph[,",i,"],tempgraph$original_PG))")
  p <- p + eval(parse(text=loop_input))
}
p+xlab("BMI")+ylab("Polygenic risk score")+theme(axis.text=element_text(size=16),axis.title=element_text(size=20))

#if you want to put the actual points in it too:
tempgraph<-subset(unrelated_eur_df,select=c(original_PG,BMI,BMIv2,BMIv3))
tempgraph<-tempgraph[unrelated_eur_df$cohort=="TwinsUK",]
p <- ggplot(data = tempgraph)
for (i in 2:4){
  cat=as.character(i)
  loop_input<-paste("geom_smooth(se=FALSE,color=cat,aes(tempgraph[,",i,"],tempgraph$original_PG))")
  p <- p + eval(parse(text=loop_input))
  loop_input<-paste("geom_point(color=cat,aes(tempgraph[,",i,"],tempgraph$original_PG))")
  p <- p + eval(parse(text=loop_input))
}
p

#compare polygenic risk score to each trait
polygen_vs_traits<-regress_multivar_o(unrelated_eur_df[,c(17,6:14,18:21,23:27,29:31,33:36)],0,1,"linear")

#compare mc4r carriers to others
mc4r_vs_traits<-regress_multivar_o(unrelated_eur_df[,c(32,6:14,18:21,23:27,29:31,33:36)],0,1,"logistic")

#compare mc4r twins
#look at twin pairs where at least one carries an MC4R variant
subset(pairdf12[pairdf12$MC4R==1|pairdf12$MC4R.1==1,],select=c(BMI,BMI.1,PC1,PC1.1,MC4R,MC4R.1))
#note that the second to the last one wasn't used in the analysis b/c the family is not European ancestry
#and the 4th to last one is not used b/c it was found in the twin w/o PC1, meaning the twin who wasn't randomly assigned to the unrelated whites group

#tabulate polygenic risk score vs. mc4r vs. obesity
unrelated_eur_df$obese[unrelated_eur_df$BMI<29.95]<-0
unrelated_eur_df$obese[unrelated_eur_df$BMI>29.95]<-1
unrelated_eur_df$normalweight[unrelated_eur_df$BMI>24.95]<-0
unrelated_eur_df$normalweight[unrelated_eur_df$BMI<24.95&unrelated_eur_df$BMI>18.45]<-1
unrelated_eur_df$overweight[unrelated_eur_df$BMI>29.95|(unrelated_eur_df$BMI<24.95&unrelated_eur_df$BMI>18.45)]<-0
unrelated_eur_df$overweight[unrelated_eur_df$BMI>24.95&unrelated_eur_df$BMI<29.95]<-1
unrelated_eur_df$firstPGquartile[unrelated_eur_df$PG>0.25]<-0
unrelated_eur_df$firstPGquartile[unrelated_eur_df$PG<0.25]<-1
#number of people in each category (obese/not obese vs. first quartile polygenic risk score/not first quartile)
table(unrelated_eur_df$obese,unrelated_eur_df$MC4R)
table(unrelated_eur_df$obese,unrelated_eur_df$firstPGquartile)
#percent of people in each category who are mc4r carriers
table(unrelated_eur_df$obese[unrelated_eur_df$MC4R==1],unrelated_eur_df$firstPGquartile[unrelated_eur_df$MC4R==1])/table(unrelated_eur_df$obese,unrelated_eur_df$firstPGquartile)

#now the proportion in normal weight individuals:
table(unrelated_eur_df$normalweight,unrelated_eur_df$MC4R)
table(unrelated_eur_df$normalweight,unrelated_eur_df$firstPGquartile)
table(unrelated_eur_df$normalweight[unrelated_eur_df$MC4R==1],unrelated_eur_df$firstPGquartile[unrelated_eur_df$MC4R==1])/table(unrelated_eur_df$normalweight,unrelated_eur_df$firstPGquartile)
table(unrelated_eur_df$normalweight[unrelated_eur_df$MC4R==1])/table(unrelated_eur_df$normalweight)

#make a plot like Figure S7B:
forbarplot<-c(sum(na.omit(unrelated_eur_df$MC4R[unrelated_eur_df$obese==1&unrelated_eur_df$firstPGquartile==1]))/length(na.omit(unrelated_eur_df$MC4R[unrelated_eur_df$obese==1&unrelated_eur_df$firstPGquartile==1])),sum(na.omit(unrelated_eur_df$MC4R[unrelated_eur_df$overweight==1&unrelated_eur_df$firstPGquartile==1]))/length(na.omit(unrelated_eur_df$MC4R[unrelated_eur_df$overweight==1&unrelated_eur_df$firstPGquartile==1])),sum(na.omit(unrelated_eur_df$MC4R[unrelated_eur_df$obese==1&unrelated_eur_df$firstPGquartile==0]))/length(na.omit(unrelated_eur_df$MC4R[unrelated_eur_df$obese==1&unrelated_eur_df$firstPGquartile==0])),sum(na.omit(unrelated_eur_df$MC4R[unrelated_eur_df$overweight==1&unrelated_eur_df$firstPGquartile==0]))/length(na.omit(unrelated_eur_df$MC4R[unrelated_eur_df$overweight==1&unrelated_eur_df$firstPGquartile==0])),sum(na.omit(unrelated_eur_df$MC4R[unrelated_eur_df$normalweight==1]))/length(na.omit(unrelated_eur_df$MC4R[unrelated_eur_df$normalweight==1])))
barplot(forbarplot,cex.axis=1.5)


