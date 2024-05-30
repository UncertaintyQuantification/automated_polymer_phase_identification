library(stringr)
library(ggplot2)
library(readxl)
library(Rcpp)
library(RcppEigen)
library(randomForest)
library(caret)
library(xgboost)
library(FastGaSP)
#library(ggtext)
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

rm_not_used_ind=T
if(rm_not_used_ind){
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
}




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
  outputs_train=log10(I_q)
  m_all=optim(c(log(1),log(.1)),get_log_lik_sum,input=q_all[truncate_ind],outputs=outputs_train[,truncate_ind],kernel_type=kernel_type,
              homo=homo,error_vars=log10_var[,truncate_ind],method="L-BFGS-B",control = list(fnscale=-1))
  param=m_all$par
}else if(est_method=="DIS_samples"){
  ##use all training DIS samples to estimate parameters
  outputs_train=log10(I_q)
  DIS_train_ind=which(label=='DIS')
  m_DIS=optim(c(log(1),log(.1)),get_log_lik_sum,input=q_all[truncate_ind],outputs=outputs_train[DIS_train_ind,truncate_ind],kernel_type=kernel_type,
              homo=homo,error_vars=log10_var[DIS_train_ind,truncate_ind],method="L-BFGS-B",control = list(fnscale=-1))
  param=m_DIS$par
}else if(est_method=="manually_set"){
  ##manually set parameters
  param=c(4,-0.5)
}
param
#smooothing
smoothed_log_hetero=smooth_intensity_fix_param(q_all=q_all,log_I=log10(I_q),param=param,kernel_type=kernel_type,homo=homo,log_var=log10_var)

##2. peak detect
peak_height_cutoff=.02
#peak_height_cutoff=m_cutoff$par
#delta_q=mean(q_all[-1]-q_all[-length(q_all)]) ##they are the same
#grad2_cutoff=-.003/(delta_q^2)
#peak_slope_cutoff=0#0.008/delta_q
peak12_not_close=T
peak_res_hetero=peak_detect(q_all=q_all,smoothed_log=smoothed_log_hetero,peak_height_cutoff=peak_height_cutoff,peak12_not_close=peak12_not_close,truncate_range=truncate_range)#,peak_slope_cutoff=peak_slope_cutoff)
dim(peak_res_hetero$peak)

remove_small_peak=T
small_peak_ratio_cut=3
if(remove_small_peak){
  peak_res_hetero=update_peak_with_small_peak_removed(peak_res=peak_res_hetero,small_peak_ratio_cut=small_peak_ratio_cut)
}

##get information of gradient
grad_data_hetero=get_grad_data(peak_res_hetero$numeri_grad_record,peak_res_hetero$numeri_grad2_record,all_peak=peak_res_hetero$peak)

# #not consider na, HEX/LAM and A15 group
# na_index=which(is.na(label))
# HEX_LAM_index=which(label=="HEX/LAM")
# A15_index=which(label=='A15')
# rm_ind=c(na_index,HEX_LAM_index,A15_index)

################ plots ################
#colors_plot=c('#f0ba40','#ed7118','#13528f','#117833','#008b85','#9b1129','#a04dba','#969696')
colors_plot=c('#969696','#ed7118','#a04dba','#117833','#13528f','#9b1129')
#c("DIS","BCC","Sigma","HEX","GYR","LAM","HEX/LAM","A15")

#update correct label
label[c(138,139,165)]=c('DIS','DIS','HEX')
data_plot=data.frame(f_F=frac,Mn=weight,Morphology=label,Chemistry=chemistry)#[-rm_ind,]
levels(data_plot$Morphology)
data_plot=data_plot[which(!is.na(label)),]

###phase diagram
ggplot(data_plot,aes(x=f_F,y=Mn))+geom_point(aes(color=Morphology,fill=Morphology,shape=Chemistry))+
  scale_color_manual(values = colors_plot,labels=c("DIS","BCC",expression(sigma),"HEX","GYR","LAM"))+
  scale_fill_manual(values = alpha(colors_plot,.3),labels=c("DIS","BCC",expression(sigma),"HEX","GYR","LAM"))+
  scale_shape_manual(values=21:24,labels=c("DIS","BCC",expression(sigma),"HEX","GYR","LAM"))+ #values = 1:length(unique(all_data$Morphology)),
  #xlim(c(0,1))+
  ylim(c(0,30))+theme_bw()+
  scale_x_continuous(name=expression(f[F]), breaks=seq(0,1,.2),limits=c(0,1))+
  facet_wrap(~Chemistry,nrow=2)+
  #xlab(expression(f[F]))+
  ylab(expression(M[n]))+
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 1), 
         shape = "none")+
  theme(legend.text=element_text(size=10),axis.text=element_text(size=12),axis.title=element_text(size=15),strip.text.x = element_text(size = 10),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.text.align = 0) 

pdf('plots/phase_diagram.pdf',width = 5,height = 3.8)
ggplot(data_plot, aes(x=f_F,y=Mn, group = Morphology)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = Morphology),
                  #bins = 2,
                  show_guide = FALSE,
                  breaks=.08)+
  geom_point(aes(color=Morphology,fill=Morphology,shape=Chemistry),size=1.2)+
  scale_fill_manual(values = alpha(colors_plot,.4),labels=c("DIS","BCC",expression(sigma),"HEX","GYR","LAM"))+
  scale_color_manual(values = colors_plot,labels=c("DIS","BCC",expression(sigma),"HEX","GYR","LAM"))+
  scale_shape_manual(values=21:24,labels=c("DIS","BCC",expression(sigma),"HEX","GYR","LAM"))+
  scale_alpha_continuous(range=c(0.3))+
  facet_wrap(~Chemistry)+
  scale_x_continuous(name=expression(italic(f)[F]), breaks=seq(0,1,.2),limits=c(0,1))+
  #xlab(expression(f[F]))+xlim(c(0,1))+
  ylab(expression(italic(M)[n]~' (kDa)'))+ylim(c(0,30))+
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 1), 
         shape = F)+
  theme_bw()+
  theme(legend.text=element_text(size=10),axis.text=element_text(size=12),
        axis.title=element_text(size=15,face = 'italic'),strip.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.spacing.x = unit(.7, "lines"),
        legend.text.align = 0)
dev.off()



###sampled curves from each phase
set.seed(17)
par(mfrow=c(2,4),  mgp=c(2.2,1,0), mar=c(3.2,1.8,2,1)+.1)
for(i in 1:length(levels(label))){
  label_plot=levels(label)[i]
  ind_plot = which(label==label_plot)
  if(length(ind_plot)>0){
    if(length(ind_plot)<=10){
      #matplot(q_all,t((I_q[ind_plot,])),type="l",main=label_plot,lwd=.5,ylab="I(q)",xlab="q",log = "y") 
      matplot(q_all,t(log(I_q[ind_plot,])),type="l",lwd=.3,xlab=expression(paste('q (',ring(A),')')),ylab="",yaxt = "n",
              main=paste0("#",label_plot,": ",length(ind_plot)) ) 
      title(ylab=expression(paste(log[10]," [I(q)] (a.u.)")), mgp=c(.6,1,0),cex.lab=1.2) 
    }else{
      #matplot(q_all,t((I_q[ind_plot,][sample(length(ind_plot),10),])),type="l",ylab="I(q)",xlab="q",main=label_plot,lwd=.5,log = "y") 
      matplot(q_all,t(log(I_q[ind_plot,][sample(length(ind_plot),10),])),type="l",lwd=.3,xlab=expression(paste('q (',ring(A),')')),ylab="",yaxt = "n",
              main=paste0("#",label_plot,": ",length(ind_plot)) ) 
      title(ylab=expression(paste(log[10]," [I(q)] (a.u.)")), mgp=c(.6,1,0),cex.lab=1.2) 
    }
  }
}


#one example sample from each phase
#set.seed(12)
pdf('plots/one_sample_each_phase.pdf',width = 6,height = 3.5)
par(mfrow=c(2,3),mgp=c(2.2,1,0), mar=c(3.2,2,2,.5)+.1)
#c(15,NA,13,12,13,16)
for(i in 1:6){
  label_plot=levels(label)[i]
  ind_plot = which(label==label_plot)
  ind_plot_select=12#sample(length(ind_plot),1)
  
  logI_here=log10(I_q[ind_plot,][ind_plot_select,])
  smoothed_I_here=smoothed_log_hetero[ind_plot,][ind_plot_select,]
  
  peaks_x=q_all[peak_res_hetero$peak[ind_plot,][ind_plot_select,]]
  peaks_y=smoothed_I_here[peak_res_hetero$peak[ind_plot,][ind_plot_select,]]
  add_y=diff(range(smoothed_I_here))*.1
  
  if(label_plot == 'Sigma') label_plot = expression(sigma)
  plot(q_all,logI_here,xlab=expression(paste('q (',ring(A),')')),ylab='',main=label_plot,yaxt = "n",
       type='l',lwd=1.5,ylim=c(min(logI_here),max(logI_here)+add_y),
       cex.lab=1.4, cex.axis=1.2,cex.main=1.5)
  title(ylab=expression(paste(log[10]," [I(q)] (a.u.)")), mgp=c(.6,1,0),cex.lab=1.2) 
  # yticks <- seq(-2L,3L,1L)
  # for (i in 1:length(yticks)){
  #   axis(2, at=yticks[i], labels=bquote(10^.(yticks[i])),cex.axis=1.3)
  # }
  lines(q_all,smoothed_I_here,col=5,lwd=1.2)
  abline(v=peaks_x,col='orange',lty=2,lwd=1.5)
  
  points(peaks_x,peaks_y+add_y,pch=25,col=2,lwd=1.2,bg=2)
}
dev.off()

###tSNE plot
library(Rtsne)
set.seed(2)
tsne_intensity=Rtsne(log10(I_q))
tsne_intensity_data=data.frame(V1 =tsne_intensity$Y[,1],V2 =tsne_intensity$Y[,2],Morphology=label)
#colors_plot=c('#ed7118','#13528f','#117833','#9b1129','#a04dba','#969696')
ggplot(tsne_intensity_data,aes(V1,V2))+geom_point(aes(color=Morphology,shape=Morphology),na.rm=T,alpha=0.7,size=1.2,stroke=1)+
  #scale_colour_discrete(na.translate = F)+
  scale_color_manual(values = colors_plot,
                     labels=c("DIS","BCC",expression(sigma),"HEX","GYR","LAM"))+
  #scale_shape_discrete(na.translate = F)+
  scale_shape_manual(values = 1:length(unique(label)),
                     labels=c("DIS","BCC",expression(sigma),"HEX","GYR","LAM"))+
  theme_linedraw()+#ggtitle("t-SNE")
  xlab("t-SNE 1") + ylab("t-SNE 2")+
  theme(legend.text=element_text(size=10),axis.text=element_text(size=12),axis.title=element_text(size=13),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.text.align = 0) 


###first three locations plot
all_peak=peak_res_hetero$peak
data_peak=data.frame(Morphology=label,peak1=q_all[all_peak[,1]],
                     peak2_peak1=q_all[all_peak[,2]-all_peak[,1]], # location of second all_peak - location of the first all_peak
                     peak2=q_all[all_peak[,2]],
                     peak3_peak2=q_all[all_peak[,3]-all_peak[,2]], # location of third all_peak - location of the second all_peak
                     peak1_grad2=grad_data_hetero$peak1_grad2_value,
                     peak1_width=grad_data_hetero$peak1_grad2_width_to_zero[,3])
## peak1 v.s. peak2-peak1
ggplot(data_peak,aes(peak1,peak2_peak1))+geom_point(aes(color=Morphology,shape=Morphology),na.rm=T,alpha=0.7,size=1.2,stroke=1)+
  #scale_colour_discrete(na.translate = F)+
  scale_color_manual(values = colors_plot,na.translate = F)+
  #scale_shape_discrete(na.translate = F)+
  scale_shape_manual(values = 1:length(unique(data_peak$Morphology)),na.translate = F)+
  xlim(c(0.02,.11))+
  #ylim(c(0,.15))+
  ylab("peak2 - peak1")+theme_linedraw()+
  theme(legend.text=element_text(size=10),axis.text=element_text(size=12),axis.title=element_text(size=13),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
## peak2 v.s. peak3-peak2
ggplot(data_peak,aes(peak2,peak3_peak2))+geom_point(aes(color=Morphology,shape=Morphology),na.rm=T,alpha=0.7,size=1.2,stroke=1)+
  #scale_colour_discrete(na.translate = F)+
  scale_color_manual(values = colors_plot,na.translate = F)+
  #scale_shape_discrete(na.translate = F)+
  scale_shape_manual(values = 1:length(unique(data_peak$Morphology)),na.translate = F)+
  xlim(c(0.04,.18))+
  #ylim(c(0,.12))+
  ylab("peak3 - peak2")+theme_linedraw()+
  theme(legend.text=element_text(size=10),axis.text=element_text(size=12),axis.title=element_text(size=13),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 




