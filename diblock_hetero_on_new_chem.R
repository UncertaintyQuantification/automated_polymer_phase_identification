library(stringr)
library(ggplot2)
library(readxl)
library(Rcpp)
library(RcppEigen)
library(randomForest)
library(caret)
library(xgboost)
library(FastGaSP)
library(latex2exp)
library(tidyverse)
sourceCpp(file='src/fgasp_modify.cpp')  ###
source('functions/functions_fgasp_modify.R')
source('functions/functions.R')
#first four sheets are diblock, and the last two sheets are triblock
sheets_name=excel_sheets('data/MachineLearning_Data.xlsx')
#change the format of sheets_name
for(i in 1:length(sheets_name)){
  sheets_name[i]=strsplit(sheets_name[i]," ",fixed=T)[[1]][1]
}
#change '-' to '_' in the sheet name
data_name=gsub("-","_",sheets_name)
#read data
column_names=c('name','Morphology','f_F','Mn_D','Mn_F','Mn_all')
for(i in 1:4){
  temp=read_xlsx("data/MachineLearning_Data.xlsx",sheet = i)[,1:6]
  #change column names
  colnames(temp)=column_names
  #remove the empty rows
  temp=temp[which(!is.na(temp[,1])),]
  #treat "---" and "" the same as NA
  temp$Morphology=factor(temp$Morphology,exclude = c("---",""),
                         levels = c("DIS","BCC","Sigma","HEX","GYR","LAM","HEX/LAM","A15"))
  temp$chemistry = sheets_name[i]
  assign(data_name[i],temp)
}


#combine dataset from different sheets
all_data=as.data.frame(do.call(rbind,mget(data_name[1:4])))

sum(str_sub(all_data$name,-1,-1)=='C')

table(all_data[,2])

files=list.files("data/Diblocks")

header_rec=rep(NA,length(files))
p_rec = rep(NA,length(files))

#most files do not have header, some have header
#check which file has the header, and the dimension of each file
for(i in 1:length(files)){
  dat=read.table(paste("data/Diblocks/",files[i],sep = ""), header = FALSE)
  header_rec[i]=is.character(dat$V1)
  if(is.character(dat$V1)){
    p_rec[i]=dim(dat)[1]-1
  }else{
    p_rec[i]=dim(dat)[1]
  }
}


sum(header_rec) ## 7 files have header
table(p_rec)
### 261 files have 485 rows
### 137 files have 491 rows 
### 7 files have 500 rows (with 1 row for column name)

p1=485;n1=261
p2=491;n2=137
p3=500;n3=7

path="data/Diblocks/"
data_list=load_data(p1=p1,n1=n1,p2=p2,n2=n2,p3=p3,n3=n3,path=path,files=files,header_rec=header_rec,all_data=all_data,morph_levels=levels(all_data$Morphology))

###combine all three sets
q_all=data_list$p1_q[1,]
I_q=combine_data_as_same_dim(xout=q_all,sets1=data_list$p1_I_q,x2=data_list$p2_q[1,],sets2=data_list$p2_I_q,x3=data_list$p3_q[1,],sets3=data_list$p3_I_q)
error=combine_data_as_same_dim(xout=q_all,sets1=data_list$p1_error,x2=data_list$p2_q[1,],sets2=data_list$p2_error,x3=data_list$p3_q[1,],sets3=data_list$p3_error)
# log_var=matrix(NA,dim(I_q)[1],p1)
# for(i in 1:dim(I_q)[1]){
#   log_var[i,]=(error[i,]^2)/(I_q[i,]^2) #variance in log space
# }
log10_var=matrix(NA,dim(I_q)[1],p1)
for(i in 1:dim(I_q)[1]){
  log10_var[i,]=(error[i,]^2)/(I_q[i,]^2*log(10)^2) #variance in log space
}
label=factor(c(data_list$p1_label,data_list$p2_label,data_list$p3_label))
name=c(data_list$p1_name,data_list$p2_name,data_list$p3_name)
frac=c(data_list$p1_frac,data_list$p2_frac,data_list$p3_frac)
weight=c(data_list$p1_Mn,data_list$p2_Mn,data_list$p3_Mn)
weight_D=c(data_list$p1_Mn_D,data_list$p2_Mn_D,data_list$p3_Mn_D)
weight_F=c(data_list$p1_Mn_F,data_list$p2_Mn_F,data_list$p3_Mn_F)
chemistry=factor(c(data_list$p1_chemistry,data_list$p2_chemistry,data_list$p3_chemistry),levels = sheets_name[1:4])

sum(str_sub(name,-1,-1)=='C',na.rm = T)

#not consider na, HEX/LAM and A15 group
na_index=which(is.na(label))
HEX_LAM_index=which(label=="HEX/LAM")
A15_index=which(label=='A15')
rm_ind=c(na_index,HEX_LAM_index,A15_index)

I_q=I_q[-rm_ind,]
error=error[-rm_ind,]
#log_var=log_var[-rm_ind,]
log10_var=log10_var[-rm_ind,]
label=factor(label[-rm_ind])
name=name[-rm_ind]
frac=frac[-rm_ind]
weight=weight[-rm_ind]
weight_D=weight_D[-rm_ind]
weight_F=weight_F[-rm_ind]
chemistry=chemistry[-rm_ind]

table(chemistry,label)

####choose the group of new chemistry
new_chem='D-12F'#c('D-9F','D-12F')
new_chem_ind=which(chemistry==new_chem)


##1. smooth 
homo=F 
kernel_type='matern_5_2' #c('matern_5_2','exp')
truncate_range_param_est=c(0.01,0.2) ##this is only used for range parameters?
truncate_range=c(0.01,0.2) ##these two can be different, but I let them be the same
truncate_ind=which(q_all>truncate_range_param_est[1] & q_all<=truncate_range_param_est[2]) # truncate as in peak detection step
#truncate_ind=1:length(q_all)
#estimate parameters with training data
est_method="DIS_samples" #c("one_sample","all_samples","DIS_samples","manually_set")
if(est_method=="one_sample"){
  ##use one sample estimate parameters
  ind_optim=143#36
  param=get_param(input=q_all[truncate_ind],output=log10(I_q[ind_optim,truncate_ind]),kernel_type=kernel_type,homo=homo,error_var=log10_var[ind_optim,truncate_ind])
}else if(est_method=="all_samples"){
  ##use all training samples to estimate parameters
  outputs_train=log10(I_q[-new_chem_ind,])
  m_all=optim(c(log(1),log(.1)),get_log_lik_sum,input=q_all[truncate_ind],outputs=outputs_train[,truncate_ind],kernel_type=kernel_type,
              homo=homo,error_vars=log10_var[-new_chem_ind,truncate_ind],method="L-BFGS-B",control = list(fnscale=-1))
  param=m_all$par
}else if(est_method=="DIS_samples"){
  ##use all training DIS samples to estimate parameters
  outputs_train=log10(I_q[-new_chem_ind,])
  DIS_train_ind=which(label[-new_chem_ind]=='DIS')
  m_DIS=optim(c(log(1),log(.1)),get_log_lik_sum,input=q_all[truncate_ind],outputs=outputs_train[DIS_train_ind,truncate_ind],kernel_type=kernel_type,
              homo=homo,error_vars=log10_var[-new_chem_ind,truncate_ind][DIS_train_ind,],method="L-BFGS-B",control = list(fnscale=-1))
  param=m_DIS$par
}else if(est_method=="manually_set"){
  ##manually set parameters
  param=c(5,5)
}
param
#smooothing
smoothed_log_hetero=smooth_intensity_fix_param(q_all=q_all,log_I=log10(I_q),param=param,kernel_type=kernel_type,homo=homo,log_var=log10_var)


##2. peak detect
peak_height_cutoff=.02
peak12_not_close=T
peak_res_hetero=peak_detect(q_all=q_all,smoothed_log=smoothed_log_hetero,peak_height_cutoff=peak_height_cutoff,peak12_not_close=peak12_not_close,truncate_range=truncate_range)#,peak_slope_cutoff=peak_slope_cutoff)
dim(peak_res_hetero$peak)


remove_small_peak=T
small_peak_ratio_cut=3
if(remove_small_peak){
  peak_res_hetero=update_peak_with_small_peak_removed(peak_res=peak_res_hetero,small_peak_ratio_cut=small_peak_ratio_cut)
}


# par(mfrow=c(1,1))
# i=336#316#228
# label[i]
# plot(log10(I_q[i,]),pch=20)#,main=paste0('log(beta): ',param[1],'  log(nugget): ',param[2]))
# lines(smoothed_log_hetero[i,],col=3)
# abline(v=peak_res_hetero$peak[i,],col=2)
# lines(log(I_q[i,])-10*sqrt(log_var[i,]),lty=2,col=4)
# lines(log(I_q[i,])+10*sqrt(log_var[i,]),lty=2,col=4)




##get information of gradient
grad_data_hetero=get_grad_data(peak_res_hetero$numeri_grad_record,peak_res_hetero$numeri_grad2_record,all_peak=peak_res_hetero$peak)

# data_grad=data.frame(label=label,grad2=-grad_data_hetero$peak1_grad2_value*diff(q_all)[1]^2,width=grad_data_hetero$peak1_grad2_width_to_zero[,3])
# ggplot(data_grad, aes(x=width, fill=label)) +
#   geom_histogram( color="#e9ecef", alpha=0.6) +
#   scale_fill_manual(values=colors_plot) + facet_wrap(~label) + xlim(c(0,40))


###3.two-step prediction

#settings for the model
add_ratios=T
add_peak1_grad2=T
add_peak1_grad2_width=T
add_frac=F
add_weight=F
add_weight_D=F
add_weight_F=F


n_peak=3
dat=get_data(model_peak=peak_res_hetero$peak,model_label=label,n_peak=n_peak,
             add_ratios,add_peak1_grad2,add_peak1_grad2_width,add_frac,add_weight,add_weight_D,add_weight_F,
             peak1_grad2_value=grad_data_hetero$peak1_grad2_value,peak1_grad2_width_to_zero=grad_data_hetero$peak1_grad2_width_to_zero[,3],
             frac=frac,weight=weight,weight_D=weight_D,weight_F=weight_F)
#split training and testing
test_ind=which(chemistry==new_chem)
train_data=dat$model_data[-test_ind,]
train_label=dat$model_label[-test_ind]
test_data=dat$model_data[test_ind,]
test_label=dat$model_label[test_ind]
# two-step prediction model
## get the index of the samples that have peaks for training and testing
train_have_peak_ind=which(apply(peak_res_hetero$peak[-test_ind,],1,function(x) sum(!is.na(x)))>0)
test_have_peak_ind=which(apply(peak_res_hetero$peak[test_ind,],1,function(x) sum(!is.na(x)))>0)
##for samples without peaks detected, predict as DIS
pred_label=factor(rep(NA,length(test_ind)),levels=levels(dat$model_label))
pred_label[-test_have_peak_ind]='DIS'
pred_prob=matrix(NA,length(test_ind),length(levels(label)),dimnames = list(NULL,levels(label)))
pred_prob[-test_have_peak_ind,]=rep(c(0,0,0,0,0,1),each=length(test_ind)-length(test_have_peak_ind))

##for the remaining samples, build classification models
set.seed(1)
model_used='rforest' #c('bagging','rforest','xgboost')
model_hetero=classification_train_test(train_input=train_data[train_have_peak_ind,],train_output=train_label[train_have_peak_ind],
                                    test_input=test_data[test_have_peak_ind,],test_output=test_label[test_have_peak_ind],model_used=model_used)
pred_label[test_have_peak_ind]=model_hetero$class_pred_label
pred_prob[test_have_peak_ind,]=model_hetero$pred_prob
mean(pred_label==test_label)
sum(pred_label!=test_label)
wrong_ind=which(test_label!= pred_label)


par(mfrow=c(1,1), mgp=c(2,.5,0), mar=c(3,3,2,.5)+.1)
if(model_used != 'xgboost'){
  #model_hetero$model$importanceSD
  var.name=c(expression(p[1]),expression(p[2]/p[1]),expression(p[3]/p[1]),expression(p[3]/p[2]),expression(paste(y^"''",(p[1]))),expression(w[1]))
  import=model_hetero$model$importance[,7]*100
  import_sort=sort(import,index.return=T)
  if(model_used=='rforest') title_model='random forest'
  if(model_used=='bagging') title_model='bagging'
  barplot(import_sort$x, horiz=TRUE, las=1,col='#A9F5B0',cex.main=1,
          main=paste('Variable Importance of',title_model,'\nfor predicting',new_chem),names.arg=var.name[import_sort$ix])
  #varImpPlot(model_hetero$model)
}else{
  importance_matrix = xgb.importance(model = model_hetero$model)
  #xgb.plot.importance(importance_matrix)
  var.name=c(expression(p[1]),expression(p[2]/p[1]),expression(p[3]/p[1]),expression(p[3]/p[2]),expression(paste(y^"''",(p[1]))),expression(w[1]))
  var_name=importance_matrix$Feature
  for(i in 1:length(var_name)) var_name[which(var_name==colnames(train_data)[i])]=var.name[i]
  barplot(rev(importance_matrix$Gain*100), horiz=TRUE, las=1,col='#A9F5B0',cex.main=1,
          main=paste('Variable Importance of extreme gradient boosting \nfor predicting',new_chem),names.arg=rev(var_name))
}





wrong_ind=which(test_label!= pred_label)
length(wrong_ind)
pred_prob[wrong_ind,]
max_prob_incorrect=apply(matrix(pred_prob[wrong_ind,],length(wrong_ind)),1,max)
max_prob_correct=apply(pred_prob[-wrong_ind,],1,max)
sort(max_prob_incorrect)
sort(max_prob_correct)

sum(max_prob_correct>.7)/length(max_prob_correct)

sum(c(max_prob_incorrect,max_prob_correct)<.8)

max_pred_prob_dat=data.frame(max_prob=c(max_prob_incorrect,max_prob_correct),
                             Type=rep(c('Incorrect','Correct'),c(length(max_prob_incorrect),length(max_prob_correct))))

pdf(paste0('plots/pred_prob_group',which(levels(chemistry)==new_chem),'.pdf'),width=3,height=1.8)
#pdf(paste0('plots/xg_pred_prob_D_12F.pdf'),width=3,height=1.8)
ggplot(max_pred_prob_dat, aes(Type, max_prob))+
  #geom_boxplot(alpha=.3)+
  geom_violin(aes(color=Type),scale = "area")+#,alpha=.3
  scale_color_manual(values = c("#14775a", "#E69F00"))+
  geom_dotplot(aes(fill=Type),stroke=.3,binaxis = "y", stackdir ="center",binwidth=.008,dotsize=3.0,stackratio=.6)+
  scale_fill_manual(values = c("#8dd3c7","#f8e66a"))+
  ylab("Max Prob")+
  theme_linedraw()+
  theme(legend.position = "none",axis.title.x=element_blank(),
        legend.text=element_text(size=10),axis.text=element_text(size=12),axis.title=element_text(size=13),
        panel.grid.major = element_line(colour = "gray80",linetype='dotted'),panel.grid.minor = element_blank())
dev.off()


get_conf_mat(model_label=test_label,pred_label=pred_label)
plot_ind=which(test_label != pred_label)
pdf(paste0(new_chem,"_missclass.pdf"), width=5*1.5, height=1.5)
par(mfrow=c(1,5), mgp=c(2,1,0), mar=c(3,3,2,.5)+.1)
for(i in 1:length(plot_ind)){
  plot_i=plot_ind[i]
  plot(q_all, log10(I_q[test_ind,][plot_i,]), type='l', lwd=1.5,
       xlab='q', ylab=expression(paste(log[10], " [I(q)] (a.u.)")), cex.main=1,
       main=TeX(paste("True:",ifelse(test_label[plot_i]=='Sigma','$sigma$',as.character(test_label[plot_i])),
                      ', Pred:', ifelse(pred_label[plot_i]=='Sigma','$sigma$',as.character(pred_label[plot_i])) )),
       cex.main=1)
  #main=paste(plot_i,"- True:",test_label[plot_i],", Pred:", pred_label[plot_i]))
  lines(q_all,smoothed_log_hetero[test_ind,][plot_i,],col=5)
  abline(v=q_all[peak_res_hetero$peak[test_ind,][plot_i,]],col=7,lty=4,lwd=1.5)
}
dev.off()
#c(1:length(label))[test_ind][plot_ind][c(83,84,86)]
peak_res_hetero$peak_height[test_ind,][plot_ind,]





#plot all samples with max prob smaller than some threshold:
small_prob_ind=which(apply(pred_prob,1,max)<.6)
par(mfrow=c(3,4), mgp=c(2,1,0), mar=c(3,3,2,.5)+.1)
for(i in 1:length(small_prob_ind)){
  plot_i=small_prob_ind[i]
  plot(q_all,log10(I_q[test_ind,][plot_i,]),pch=20,cex=.5,xlab='q',ylab="log I(q)")
  title(paste("True:",test_label[plot_i],", Pred:", pred_label[plot_i]),col.main=ifelse(plot_i %in% wrong_ind,6,1),cex.main=1)
  #main=paste(plot_i,"- True:",test_label[plot_i],", Pred:", pred_label[plot_i]))
  lines(q_all,smoothed_log_hetero[test_ind,][plot_i,],col=4)
  abline(v=q_all[peak_res_hetero$peak[test_ind,][plot_i,]],col=2,lty=4)
}

#maximum prob of correctly predicted by morphology
correct_pred_prob_dat=data.frame(max_prob=max_prob_correct,phase=test_label[-wrong_ind])
ggplot(correct_pred_prob_dat, aes(phase, max_prob))+
  #geom_boxplot(alpha=.3)+
  geom_violin(aes(color=phase),scale = "area")+#,alpha=.3
  #scale_color_manual(values = c("#14775a", "#E69F00"))+
  geom_dotplot(aes(fill=phase),stroke=.3,binaxis = "y", stackdir ="center",binwidth=.008,dotsize=2.0,stackratio=.7)+
  #scale_fill_manual(values = c("#8dd3c7","#f8e66a"))+
  ylab("Max Prob")+
  theme_linedraw()+
  theme(legend.position = "none",axis.title.x=element_blank(),
        legend.text=element_text(size=10),axis.text=element_text(size=12),axis.title=element_text(size=13))


##plot all testing
par(mfrow=c(4,5), mgp=c(2,1,0), mar=c(3,3,2,.5)+.1)
for(i in 1:dim(I_q[test_ind,])[1]){
  plot_i=i
  plot(q_all,log10(I_q[test_ind,][plot_i,]),pch=20,cex=.5,xlab='q',ylab="I(q)",cex.main=1,
       main=paste("True:",label[test_ind][plot_i]),type='l')
  lines(q_all,smoothed_log_hetero[test_ind,][plot_i,],col=4)
  abline(v=q_all[peak_res_hetero$peak[test_ind,][plot_i,]],col=2,lty=4)
}


##plot all
par(mfrow=c(4,5), mgp=c(2,1,0), mar=c(3,3,2,.5)+.1)
for(i in 1:dim(I_q)[1]){
  plot_i=i
  plot(q_all,log10(I_q[plot_i,]),pch=20,cex=.5,xlab='q',ylab="I(q)",cex.main=1,
       main=paste("True:",label[plot_i]),type='l')
  lines(q_all,smoothed_log_hetero[plot_i,],col=4)
  abline(v=q_all[peak_res_hetero$peak[plot_i,]],col=2,lty=4)
}


##plot one class
class_select='Sigma'
phase_ind=which(label==class_select)
phase_ind=phase_ind[!phase_ind %in% test_ind]
pdf(paste0("all_",class_select,"_train_samples_",new_chem,".pdf"))
par(mfrow=c(5,6), mgp=c(2,1,0), mar=c(3,3,2,.5)+.1)
for(i in 1:length(phase_ind)){
  plot_i=phase_ind[i]
  plot(q_all,log10(I_q[plot_i,]),pch=20,cex=.5,xlab='q',ylab="I(q)",cex.main=1,
       main=paste0(plot_i,"-True: ",label[plot_i]),type='l')
  lines(q_all,smoothed_log_hetero[plot_i,],col=4)
  abline(v=q_all[peak_res_hetero$peak[plot_i,]],col=2,lty=4)
}
dev.off()

label_plot = c("BCC","HEX","LAM", "Sigma", "DIS")
#label_ind=c(8,19,,,)
n_plot = 2
set.seed(6)
for(i in 1:n_plot){
  pdf(paste0("true_sample",i,".pdf"), width=5*1.5, height=1.5)
  par(mfrow=c(1,5), mgp=c(2,1,0), mar=c(3,3,2,.5)+.1)
  for (label_here in label_plot){
    label_ind = which(train_label==label_here)
    plot_i=sample(label_ind,1)
    plot(q_all,log10(I_q[-test_ind,][plot_i,]),type='l',lwd=1.5,xlab='q',ylab=expression(paste(log[10]," [I(q)] (a.u.)")),cex.main=1,
         main=TeX(paste("True:",ifelse(train_label[plot_i]=='Sigma','$sigma$',as.character(train_label[plot_i])) )),cex.main=1)
    #main=paste(plot_i,"- True:",test_label[plot_i],", Pred:", pred_label[plot_i]))
    lines(q_all,smoothed_log_hetero[-test_ind,][plot_i,],col=5)
    abline(v=q_all[peak_res_hetero$peak[-test_ind,][plot_i,]],col=7,lty=4,lwd=1.5)
  }
  dev.off()
}




#combine dataset from different sheets
colors_plot=c('#f0ba40','#ed7118','#13528f','#117833','#008b85','#9b1129','#a04dba','#969696')
#colors_plot=c('#ed7118','#13528f','#117833','#9b1129','#a04dba','#969696')

data_plot=data.frame(f_F=frac,Mn=weight,Morphology=label,Chemistry=chemistry)#[-rm_ind,]
levels(data_plot$Morphology)
data_plot=data_plot[which(!is.na(label)),]

# ggplot(data_plot,aes(x=f_F,y=Mn))+geom_point(aes(color=Morphology,fill=Morphology,shape=Chemistry))+
#   scale_color_manual(values = colors_plot,na.translate = F)+
#   scale_fill_manual(values = alpha(colors_plot,.3),na.translate = F)+
#   scale_shape_manual(values=21:24,na.translate = F)+ #values = 1:length(unique(all_data$Morphology)),
#   #xlim(c(0,1))+
#   ylim(c(0,30))+theme_bw()+
#   scale_x_continuous(name=expression(f[F]), breaks=seq(0,1,.2),limits=c(0,1))+
#   facet_wrap(~Chemistry,nrow=1)+
#   #xlab(expression(f[F]))+
#   ylab(expression(M[n]))+
#   guides(colour = guide_legend(order = 1), fill = guide_legend(order = 1), 
#          shape = "none")+
#   theme(legend.text=element_text(size=10),axis.text=element_text(size=12),axis.title=element_text(size=15),strip.text.x = element_text(size = 10)) 
# 
# table(data_plot$chemistry,data_plot$Morphology)


data_plot$new_group = factor(ifelse(data_plot$Chemistry!=new_chem,"Train",paste0('Test: ',new_chem)),levels=c("Train",paste0('Test: ',new_chem)))
ggplot(data_plot,aes(x=f_F,y=Mn))+geom_point(aes(color=Morphology,fill=Morphology,shape=Chemistry))+
  scale_color_manual(values = colors_plot,na.translate = F)+
  scale_fill_manual(values = alpha(colors_plot,.3),na.translate = F)+
  scale_shape_manual(values=21:24,na.translate = F)+ #values = 1:length(unique(all_data$Morphology)),
  xlim(c(0,1))+ylim(c(0,30))+theme_bw()+
  facet_wrap(~new_group)+
  xlab(expression(f[F]))+ylab(expression(M[n]))+
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 1), 
         shape = guide_legend(order = 2))

ggplot(data_plot, aes(x=f_F,y=Mn, group = Morphology)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Morphology),
                  #bins = 2,
                  show_guide = FALSE,
                  breaks=.08)+
  geom_point(aes(color=Morphology,fill=Morphology,shape=Chemistry),size=.8)+
  scale_fill_manual(values = alpha(colors_plot[-c(1,5)],.5),na.translate = F)+
  scale_color_manual(values = colors_plot[-c(1,5)],na.translate = F)+
  scale_shape_manual(values=21:24,na.translate = F)+
  facet_wrap(~new_group)+
  xlim(c(0,1))+ylim(c(0,30))+theme_bw()+
  xlab(expression(f[F]))+ylab(expression(M[n]))+
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 1), 
         shape = F)


feature_data = dat$model_data
feature_data[which(is.na(feature_data))]=0
feature_data=as.data.frame(feature_data)
feature_data$phase = as.character(dat$model_label)
feature_data$Group = 'Train'
feature_data$Group[test_ind] = 'Correct test'
feature_data$Group[test_ind][wrong_ind] = 'Incorrect test'
feature_data$phase = factor(feature_data$phase, levels = levels(label))
feature_data$Group = factor(feature_data$Group, levels = c('Incorrect test','Correct test','Train'))


set.seed(2)
var.name=c(expression(p[1]),expression(p[2]/p[1]),expression(p[3]/p[1]),expression(p[3]/p[2]),expression(paste(y^"''",(p[1]))),expression(w[1]))
for(i in 1:6){
  pdf(paste0(new_chem,"_feature",i,".pdf"), width = 4, height = 2.5)
  y_axis = sym(names(feature_data)[i])
  print(feature_data %>% arrange(desc(Group)) %>%
          ggplot()+
          geom_jitter(aes(x=phase,y=!!y_axis,color=Group,shape = phase),stroke=1)+
          xlab('Phase')+ylab(var.name[i])+
          scale_y_continuous(breaks=c(0))+
          scale_x_discrete(labels=c("DIS","BCC",expression(sigma),"HEX","GYR","LAM"))+
          scale_shape_manual(values = 1:6,na.translate = F)+
          theme_linedraw()+
          theme(legend.position = if(i==1){c(0.82,0.8)}else{'none'},legend.text=element_text(size=10),
                axis.text=element_text(size=12),axis.title=element_text(size=13),
                legend.background=element_rect(fill = alpha("white", 0)),
                legend.key = element_rect(colour = "transparent", fill = "transparent"),
                legend.title=element_blank(),
                panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
          guides(shape = 'none')
        )
  dev.off()
}


