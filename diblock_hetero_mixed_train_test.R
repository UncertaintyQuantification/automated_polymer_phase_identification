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
#### For I_q and error: linear interpolate the second and third groups to have same dimensions as the first group
I_q=combine_data_as_same_dim(xout=q_all,sets1=data_list$p1_I_q,x2=data_list$p2_q[1,],sets2=data_list$p2_I_q,x3=data_list$p3_q[1,],sets3=data_list$p3_I_q)
error=combine_data_as_same_dim(xout=q_all,sets1=data_list$p1_error,x2=data_list$p2_q[1,],sets2=data_list$p2_error,x3=data_list$p3_q[1,],sets3=data_list$p3_error)
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


#not consider na, HEX/LAM and A15 group, removing those samples
na_index=which(is.na(label))
HEX_LAM_index=which(label=="HEX/LAM")
A15_index=which(label=='A15')
rm_ind=c(na_index,HEX_LAM_index,A15_index)

I_q=I_q[-rm_ind,]
error=error[-rm_ind,]
log10_var=log10_var[-rm_ind,]
label=factor(label[-rm_ind])
name=name[-rm_ind]
frac=frac[-rm_ind]
weight=weight[-rm_ind]
weight_D=weight_D[-rm_ind]
weight_F=weight_F[-rm_ind]
chemistry=chemistry[-rm_ind]

#0. n-fold cross validation
nfold=5 #number of folds to use
set.seed(1)
N=length(label)
folds_all=get_folds(nfold=nfold,size=N) #get the samples for each fold

table(label,folds_all)
correct_misidentified = T
if(correct_misidentified){
  label[c(138,139,165)]=c('DIS','DIS','HEX') #correcting labels
}

model_used='rforest' #c('bagging','rforest','xgboost')
## peak height cutoff 
peak_height_cutoff=.02
## decide whether to remove small peaks
remove_small_peak=T
small_peak_ratio_cut=3

#record predicted label, prob, and the accuracy on each fold
pred_label_all=factor(rep(NA,N),levels=levels(label)) 
pred_prob_all=matrix(NA,N,length(levels(label)),dimnames = list(NULL,levels(label)))
accuracy_fold=rep(NA,nfold)
smoothed_log_with_test=matrix(NA,N,length(q_all))
peak_with_test=matrix(NA,N,15)
for(i_fold in 1:nfold){
  print(i_fold)
  train_ind=which(folds_all!=i_fold)
  test_ind=which(folds_all==i_fold)
  
  #1. smooth 
  homo=F 
  kernel_type='matern_5_2' #c('matern_5_2','exp')
  truncate_range=c(.01,.2) #select the truncate range
  truncate_ind=which(q_all>truncate_range[1] & q_all<=truncate_range[2]) # truncate as in peak detection step
  
  ##estimate parameters with training data
  est_method="DIS_samples" #c("one_sample","all_samples","DIS_samples","manually_set")
  if(est_method=="one_sample"){
    ##use one sample estimate parameters
    ind_optim=143
    param=get_param(input=q_all[truncate_ind],output=log10(I_q[ind_optim,truncate_ind]),kernel_type=kernel_type,homo=homo,error_var=log10_var[ind_optim,truncate_ind])
  }else if(est_method=="all_samples"){
    ##use all training samples to estimate parameters
    outputs_train=log10(I_q[train_ind,])
    m_all=optim(c(log(1),log(.1)),get_log_lik_sum,input=q_all[truncate_ind],outputs=outputs_train[,truncate_ind],kernel_type=kernel_type,
                homo=homo,error_vars=log10_var[train_ind,truncate_ind],method="L-BFGS-B",control = list(fnscale=-1))
    param=m_all$par
  }else if(est_method=="DIS_samples"){
    ##use all training DIS samples to estimate parameters
    outputs_train=log10(I_q[train_ind,])
    DIS_train_ind=which(label[train_ind]=='DIS')
    m_DIS=optim(c(log(1),log(.1)),get_log_lik_sum,input=q_all[truncate_ind],outputs=outputs_train[DIS_train_ind,truncate_ind],kernel_type=kernel_type,
                homo=homo,error_vars=log10_var[train_ind,truncate_ind][DIS_train_ind,],method="L-BFGS-B",control = list(fnscale=-1))
    param=m_DIS$par
  }else if(est_method=="manually_set"){
    ##manually set parameters
    param=c(4,-0.5)
  }
  ##smooothing with estimated parameters
  smoothed_log_hetero=smooth_intensity_fix_param(q_all=q_all,log_I=log10(I_q),param=param,kernel_type=kernel_type,homo=homo,log_var=log10_var)
  smoothed_log_with_test[test_ind,]=smoothed_log_hetero[test_ind,]
  
  
  #2. peak detect & feature selection
  peak12_not_close=T #the first two peaks shouldn't be too close
  peak_res_hetero=peak_detect(q_all=q_all,smoothed_log=smoothed_log_hetero,peak_height_cutoff=peak_height_cutoff,peak12_not_close=peak12_not_close,truncate_range=truncate_range)#,peak_slope_cutoff=peak_slope_cutoff)
  #dim(peak_res_hetero$peak)
  
  ## decide whether to remove small peaks
  if(remove_small_peak){
    peak_res_hetero=update_peak_with_small_peak_removed(peak_res=peak_res_hetero,small_peak_ratio_cut=small_peak_ratio_cut)
  }
  peak_with_test[test_ind,1:dim(peak_res_hetero$peak)[2]]=peak_res_hetero$peak[test_ind,]
  
  ##get information of gradient
  grad_data_hetero=get_grad_data(peak_res_hetero$numeri_grad_record,peak_res_hetero$numeri_grad2_record,all_peak=peak_res_hetero$peak)
  
  #3. prediction model
  ##settings for the model: the features used in the classification models
  add_ratios=T
  add_peak1_grad2=T
  add_peak1_grad2_width=T
  add_frac=F
  add_weight=F
  add_weight_D=F
  add_weight_F=F
  add_chem=F
  
  n_peak=3
  ##obtain the data for modeling
  dat=get_data(model_peak=peak_res_hetero$peak,model_label=label,n_peak=n_peak,
               add_ratios,add_peak1_grad2,add_peak1_grad2_width,add_frac,add_weight,add_weight_D,add_weight_F,add_chem,
               peak1_grad2_value=grad_data_hetero$peak1_grad2_value,peak1_grad2_width_to_zero=grad_data_hetero$peak1_grad2_width_to_zero[,3],
               frac=frac,weight=weight,weight_D=weight_D,weight_F=weight_F,chemistry = chemistry)
  train_data=dat$model_data[train_ind,]
  train_label=dat$model_label[train_ind]
  test_data=dat$model_data[-train_ind,]
  test_label=dat$model_label[-train_ind]
  
  ### get the index of the samples that have peaks for training and testing
  train_have_peak_ind=which(apply(peak_res_hetero$peak[train_ind,],1,function(x) sum(!is.na(x)))>0)
  test_have_peak_ind=which(apply(peak_res_hetero$peak[-train_ind,],1,function(x) sum(!is.na(x)))>0)
  ### for samples without peaks detected, predict as DIS
  pred_label=factor(rep(NA,length(test_ind)),levels=levels(dat$model_label))
  pred_label[-test_have_peak_ind]='DIS'
  
  ###for the remaining samples, build classification models
  model_hetero=classification_train_test(train_input=train_data[train_have_peak_ind,],train_output=train_label[train_have_peak_ind],
                                      test_input=test_data[test_have_peak_ind,],test_output=test_label[test_have_peak_ind],model_used=model_used)
  pred_label[test_have_peak_ind]=model_hetero$class_pred_label
  pred_label_all[test_ind]=pred_label
  pred_prob_all[test_ind,][test_have_peak_ind,]=model_hetero$pred_prob
  pred_prob_all[test_ind,][-test_have_peak_ind,]=rep(c(0,0,0,0,0,1),each=length(test_ind)-length(test_have_peak_ind))
  
  accuracy_fold[i_fold]=mean(pred_label==test_label)
  
}

mean(accuracy_fold) #accuracy
sum(label != pred_label_all) #number of missclassified samples

#remove columns with all NAs
peak_with_test=peak_with_test[,apply(peak_with_test, 2, function(x) any(!is.na(x)))]


###4. analyze results
get_conf_mat(model_label=label,pred_label=pred_label_all)
wrong_ind=which(label != pred_label_all)
#wrong_ind=wrong_ind[order(label[wrong_ind])] #sort the index by the levels of label
length(wrong_ind)
par(mfrow=c(4,5), mgp=c(2,1,0), mar=c(3,3,1.2,.5)+.1)
for(i in 1:length(wrong_ind)){
  plot_i=wrong_ind[i]
  plot(q_all,log10(I_q[plot_i,]),pch=20,cex=.5,xlab='q',ylab="I(q)",cex.main=1,
       main=paste(plot_i,"- True:",label[plot_i],", Pred:", pred_label_all[plot_i]))
  lines(q_all,smoothed_log_with_test[plot_i,],col=4)
  abline(v=q_all[peak_with_test[plot_i,]],col=2,lty=4)
  max_prob=max(pred_prob_all[plot_i,])
  text(0.2, sum(par('usr')[3:4]*c(.2,.8)), paste0('prob: \n',max_prob),col=ifelse(max_prob>.8,6,4))
}#c(138,139,165) are the incorrectly labeled samples


#prediction prob of incorrectly predicted samples
pred_prob_all[wrong_ind,]
apply(pred_prob_all[wrong_ind,],1,max) #max prob for each sample
max_prob_incorrect=apply(pred_prob_all[wrong_ind,],1,max) 
max_prob_correct=apply(pred_prob_all[-wrong_ind,],1,max) 
sort(max_prob_incorrect)#sort the max prob
sort(max_prob_correct)#sort the max prob for correct prediction

incorrect_labeled_ind=which(wrong_ind %in% c(138,139,165))
group_incorrect=rep('Misclassified',length(wrong_ind))
group_incorrect[incorrect_labeled_ind]='Mislabeled'
group_correct=rep('Correctly classified',length(max_prob_correct))

max_pred_prob_dat=data.frame(max_prob=c(max_prob_incorrect,max_prob_correct),
                             Type=rep(c('Incorrect','Correct'),c(length(max_prob_incorrect),length(max_prob_correct))),
                             Group=c(group_incorrect,group_correct))

#pdf('plots/pred_prob_all_after.pdf',width=4.6,height=1.8)
ggplot(max_pred_prob_dat, aes(Type, max_prob))+
  geom_violin(aes(color=Type),scale = "area")+#,alpha=.3
  scale_color_manual(values = c("#14775a", "#E69F00"))+
  geom_dotplot(aes(fill=Group),stroke=.1,binaxis = "y", stackdir ="center",binwidth=.002,dotsize=7,stackratio=.5)+#
  scale_fill_manual(values = c("#8dd3c7","#f8e66a" , "#e41a1c"))+
  ylab("Max Prob")+
  theme_linedraw()+
  theme(axis.title.x=element_blank(),#legend.position = "none",
        legend.text=element_text(size=10),axis.text=element_text(size=12),axis.title=element_text(size=13),
        panel.grid.major = element_line(colour = "gray80",linetype='dotted'),panel.grid.minor = element_blank())+
  guides(fill = guide_legend(order = 1),colour = 'none')
#dev.off()




