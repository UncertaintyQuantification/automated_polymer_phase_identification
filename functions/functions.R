load_data=function(p1,n1,p2,n2,p3,n3,path,files,header_rec,all_data,morph_levels){
  p1_q = matrix(NA,n1,p1)
  p1_I_q = matrix(NA,n1,p1)
  p1_error = matrix(NA,n1,p1)
  p1_label = factor(rep(NA,n1),levels = morph_levels)
  p1_name = rep(NA,n1)
  p1_frac = rep(NA,n1)
  p1_Mn = rep(NA,n1)
  p1_Mn_D=rep(NA,n1)
  p1_Mn_F=rep(NA,n1)
  p1_chemistry=rep(NA,n1)
  
  p2_q = matrix(NA,n2,p2)
  p2_I_q = matrix(NA,n2,p2)
  p2_error = matrix(NA,n2,p2)
  p2_label = factor(rep(NA,n2),levels = morph_levels)
  p2_name = rep(NA,n2)
  p2_frac = rep(NA,n2)
  p2_Mn = rep(NA,n2)
  p2_Mn_D=rep(NA,n2)
  p2_Mn_F=rep(NA,n2)
  p2_chemistry=rep(NA,n2)
  
  p3_q = matrix(NA,n3,p3)
  p3_I_q = matrix(NA,n3,p3)
  p3_error = matrix(NA,n3,p3)
  p3_label = factor(rep(NA,n3),levels = morph_levels)
  p3_name = rep(NA,n3)
  p3_frac = rep(NA,n3)
  p3_Mn = rep(NA,n3)
  p3_Mn_D=rep(NA,n3)
  p3_Mn_F=rep(NA,n3)
  p3_chemistry=rep(NA,n3)

  count1=count2=count3=0
  for(i in 1:length(files)){
    if(header_rec[i]){
      dat=read.table(paste(path,files[i],sep = ""), header = TRUE)
      if(dim(dat)[1]==p3){
        count3=count3+1
        p3_q[count3,]=as.numeric(str_sub(dat[,1],1,-2))#change string to numeric
        p3_I_q[count3,]=as.numeric(str_sub(dat[,2],1,-2))
        p3_error[count3,]=dat[,3]
        ind_label = which(all_data[,1]==str_sub(files[i],1,-5))
        if (length(ind_label)==1){
          p3_label[count3] = all_data[ind_label,2]
          p3_name[count3] = all_data[ind_label,1]
          p3_frac[count3] = all_data[ind_label,3]
          p3_Mn[count3]=all_data[ind_label,6]
          p3_Mn_D[count3]=all_data[ind_label,4]
          p3_Mn_F[count3]=all_data[ind_label,5]
          p3_chemistry[count3]=all_data[ind_label,7]
        }
      }else if(dim(dat)[1]==p1){
        count1=count1+1
        p1_q[count1,]=as.numeric(str_sub(dat[,1],1,-2))#change string to numeric
        p1_I_q[count1,]=as.numeric(str_sub(dat[,2],1,-2))
        p1_error[count1,]=dat[,3]
        ind_label = which(all_data[,1]==str_sub(files[i],1,-5))
        if (length(ind_label)==1){
          p1_label[count1] = all_data[ind_label,2]
          p1_name[count1] = all_data[ind_label,1]
          p1_frac[count1] = all_data[ind_label,3]
          p1_Mn[count1]=all_data[ind_label,6]
          p1_Mn_D[count1]=all_data[ind_label,4]
          p1_Mn_F[count1]=all_data[ind_label,5]
          p1_chemistry[count1]=all_data[ind_label,7]
        }
      }else if(dim(dat)[1]==p2){
        count2=count2+1
        p2_q[count2,]=as.numeric(str_sub(dat[,1],1,-2))#change string to numeric
        p2_I_q[count2,]=as.numeric(str_sub(dat[,2],1,-2))
        p2_error[count2,]=dat[,3]
        ind_label = which(all_data[,1]==str_sub(files[i],1,-5))
        if (length(ind_label)==1){
          p2_label[count2] = all_data[ind_label,2]
          p2_name[count2] = all_data[ind_label,1]
          p2_frac[count2] = all_data[ind_label,3]
          p2_Mn[count2]=all_data[ind_label,6]
          p2_Mn_D[count2]=all_data[ind_label,4]
          p2_Mn_F[count2]=all_data[ind_label,5]
          p2_chemistry[count2]=all_data[ind_label,7]
        }
      }
    }else{
      dat=read.table(paste(path,files[i],sep = ""), header = FALSE)
      if(dim(dat)[1]==p1){
        count1=count1+1
        p1_q[count1,]=dat[,1]
        p1_I_q[count1,]=dat[,2]
        p1_error[count1,]=dat[,3]
        ind_label = which(all_data[,1]==str_sub(files[i],1,-5))
        if (length(ind_label)==1){
          p1_label[count1] = all_data[ind_label,2]
          p1_name[count1] = all_data[ind_label,1]
          p1_frac[count1] = all_data[ind_label,3]
          p1_Mn[count1]=all_data[ind_label,6]
          p1_Mn_D[count1]=all_data[ind_label,4]
          p1_Mn_F[count1]=all_data[ind_label,5]
          p1_chemistry[count1]=all_data[ind_label,7]
        }
      }else if(dim(dat)[1]==p2){
        count2=count2+1
        p2_q[count2,]=dat[,1]
        p2_I_q[count2,]=dat[,2]
        p2_error[count2,]=dat[,3]
        ind_label = which(all_data[,1]==str_sub(files[i],1,-5))
        if (length(ind_label)==1){
          p2_label[count2] = all_data[ind_label,2]
          p2_name[count2] = all_data[ind_label,1]
          p2_frac[count2] = all_data[ind_label,3]
          p2_Mn[count2]=all_data[ind_label,6]
          p2_Mn_D[count2]=all_data[ind_label,4]
          p2_Mn_F[count2]=all_data[ind_label,5]
          p2_chemistry[count2]=all_data[ind_label,7]
        }
      }else if(dim(dat)[1]==p3){
        count3=count3+1
        p3_q[count3,]=dat[,1]
        p3_I_q[count3,]=dat[,2]
        p3_error[count3,]=dat[,3]
        ind_label = which(all_data[,1]==str_sub(files[i],1,-5))
        if (length(ind_label)==1){
          p3_label[count3] = all_data[ind_label,2]
          p3_name[count3] = all_data[ind_label,1]
          p3_frac[count3] = all_data[ind_label,3]
          p3_Mn[count3]=all_data[ind_label,6]
          p3_Mn_D[count3]=all_data[ind_label,4]
          p3_Mn_F[count3]=all_data[ind_label,5]
          p3_chemistry[count3]=all_data[ind_label,7]
        }
      }
    }
  }
  res=list(p1_q=p1_q,p1_I_q=p1_I_q,p1_error=p1_error,p1_label=p1_label,p1_name=p1_name,p1_frac=p1_frac,p1_Mn=p1_Mn,p1_Mn_D=p1_Mn_D,p1_Mn_F=p1_Mn_F,p1_chemistry=p1_chemistry,
           p2_q=p2_q,p2_I_q=p2_I_q,p2_error=p2_error,p2_label=p2_label,p2_name=p2_name,p2_frac=p2_frac,p2_Mn=p2_Mn,p2_Mn_D=p2_Mn_D,p2_Mn_F=p2_Mn_F,p2_chemistry=p2_chemistry,
           p3_q=p3_q,p3_I_q=p3_I_q,p3_error=p3_error,p3_label=p3_label,p3_name=p3_name,p3_frac=p3_frac,p3_Mn=p3_Mn,p3_Mn_D=p3_Mn_D,p3_Mn_F=p3_Mn_F,p3_chemistry=p3_chemistry)
  return(res)
}


combine_data_as_same_dim=function(xout,sets1,x2,sets2,x3=NULL,sets3=NULL){
  sets2_new=t(apply(sets2, 1, function(z) approx(x=x2,y=z,xout=xout)$y))
  comb_set=rbind(sets1,sets2_new)
  if(!is.null(sets3)){
    sets3_new=t(apply(sets3, 1, function(z) approx(x=x3,y=z,xout=xout)$y))
    comb_set=rbind(comb_set,sets3_new)
  }
  return(comb_set)
}

get_param=function(input,output,kernel_type,homo=T,error_var=NULL){
  fgasp.model_select=fgasp(input=input, output=output,kernel_type=kernel_type)
  if(homo){
    est_all=optim(c(log(1),log(.02)),log_lik,object=fgasp.model_select,method="L-BFGS-B",
                  control = list(fnscale=-1))  
  }else{
    est_all=optim(c(log(1),log(.02)),log_lik_hetero,
                  object=fgasp.model_select,error_var=error_var,
                  method="L-BFGS-B",control = list(fnscale=-1))  
  }
  return(est_all$par)
}

get_log_lik_sum=function(param,input,outputs,kernel_type,homo=T,error_vars=NULL){
  n = dim(outputs)[1]
  log_lik_sum=0
  for(i in 1:n){
    fgasp.model_i=fgasp(input=input, output=outputs[i,],kernel_type=kernel_type)
    if(homo){
      log_lik_i=log_lik(param,object=fgasp.model_i)
    }else{
      log_lik_i=log_lik_hetero(param,object=fgasp.model_i,error_var=error_vars[i,])
    }
    log_lik_sum=log_lik_sum+log_lik_i
  }
  return(log_lik_sum)
}


smooth_intensity_fix_param=function(q_all,log_I,param,kernel_type,homo=T,log_var=NULL){
  smoothed_log=matrix(NA,dim(log_I)[1],dim(log_I)[2])
  
  for(i in 1:dim(log_I)[1]){
    ##let me just use a fixed parameter
    fgasp.model=fgasp(input=q_all, output=log_I[i,],kernel_type=kernel_type)
    if(homo){
      pred.model=predict(param=param,object=fgasp.model,testing_input=q_all)
    }else{
      pred.model=predict.fgasp_hetero(param=param,object=fgasp.model,testing_input=q_all,error_var=log_var[i,])
    }
    smoothed_log[i,]=pred.model@mean
  }
  
  return(smoothed_log=smoothed_log)
}

peak_detect=function(q_all,smoothed_log,peak_height_cutoff,peak12_not_close=F,truncate_range){#,peak_slope_cutoff=0){#grad2_cutoff
  peak=matrix(NA,dim(smoothed_log)[1],30) 
  peak_height=matrix(NA,dim(smoothed_log)[1],30) 
  numeri_grad_record=matrix(NA,dim(smoothed_log)[1],dim(smoothed_log)[2]-1)
  numeri_grad2_record=matrix(NA,dim(smoothed_log)[1],dim(smoothed_log)[2]-2)
  peak1_domain=matrix(NA,dim(smoothed_log)[1],2)
  
  delta_q=mean(q_all[-1]-q_all[-length(q_all)]) ##they are the sam
  
  for(i in 1:dim(smoothed_log)[1]){
    diff_log_I=diff(smoothed_log[i,])#smoothed_log[i,2:dim(smoothed_log)[2]]-smoothed_log[i,1:(dim(smoothed_log)[2]-1)]
    numeri_grad_record[i,]=diff_log_I/delta_q
    diff_log_I2=diff(diff_log_I)#diff_log_I[-1]-diff_log_I[-length(diff_log_I)]
    numeri_grad2_record[i,]=diff_log_I2/(delta_q^2)
    
    sign_diff_log_I=sign(diff_log_I)
    ##change from negative to positive
    start_of_increase=which(diff(sign_diff_log_I)>0)+1
    if(sign_diff_log_I[1]>0){ #if directly increase, then first q is a start_of_increase
      start_of_increase=c(1,start_of_increase)
    }
    end_of_decrease=start_of_increase[-1]
    ##change from positive to negative
    potential_peak=which(diff(sign_diff_log_I)<0)+1
    if(sign_diff_log_I[length(sign_diff_log_I)]>0){ #if last one is increasing, then remove the last start_of_increase
      start_of_increase=start_of_increase[-length(start_of_increase)]
    }else{
      end_of_decrease=c(end_of_decrease,length(smoothed_log[i,]))
    }
    
    #remove the peak after q>.2, which are mostly noise, and remove the peak before q<.01
    keep=q_all[potential_peak]<=truncate_range[2] & q_all[potential_peak]>truncate_range[1]
    #keep=potential_peak<=400 & potential_peak>10
    potential_peak=potential_peak[keep]
    start_of_increase=start_of_increase[keep]
    end_of_decrease=end_of_decrease[keep] 
    
    # #remove the super flat peak
    # keep=abs(diff_log_I2[potential_peak-1])>10^(-3)
    # potential_peak=potential_peak[keep]
    # start_of_increase=start_of_increase[keep]
    # end_of_decrease=end_of_decrease[keep] 
    
    if(length(potential_peak)>0){
      # #make sure the first peak is the highest peak
      # while (which.max(smoothed_log[i,potential_peak])!=1) {
      #   potential_peak=potential_peak[-1]
      #   #keep the smallest log(I) as the start_of_increase
      #   if(smoothed_log[i,start_of_increase[1]]<smoothed_log[i,start_of_increase[2]]){
      #     start_of_increase=start_of_increase[-2]
      #   }else{
      #     start_of_increase=start_of_increase[-1]
      #   }
      # }
      
      height_increase=smoothed_log[i,(potential_peak)]-smoothed_log[i,start_of_increase]
      #height_decrease=smoothed_log[i,(potential_peak)]-smoothed_log[i,end_of_decrease]
      decrease_width=end_of_decrease-potential_peak
      #remove the potential peaks with width < 3
      smooth_out_int=which(decrease_width<3)
      if(length(smooth_out_int)>0){
        for(i_smooth in rev(smooth_out_int)){
          if(i_smooth<length(potential_peak)){
            start_of_increase=start_of_increase[-(i_smooth+1)]
            potential_peak=potential_peak[-i_smooth]
            end_of_decrease=end_of_decrease[-i_smooth]
          }
        }
        height_increase=smoothed_log[i,(potential_peak)]-smoothed_log[i,start_of_increase]
        #height_decrease=smoothed_log[i,(potential_peak)]-smoothed_log[i,end_of_decrease]
        decrease_width=end_of_decrease-potential_peak
      }
      
      
      #grad2=numeri_grad2_record[i,potential_peak-1]
      #peak_slope=height_increase/((potential_peak-start_of_increase)*delta_q)
      index_peak=which(height_increase>peak_height_cutoff)#  & peak_slope>peak_slope_cutoff)# & grad2<grad2_cutoff) 
      
      if(length(index_peak)>0){
        potential_peak=potential_peak[index_peak]
        start_of_increase=start_of_increase[index_peak]
        height_increase=height_increase[index_peak]
        
        #the first two peaks can't have close log(I)
        #if(peak12_not_close & length(potential_peak)>1 & abs(smoothed_log[i,potential_peak[1]]-smoothed_log[i,potential_peak[2]])<0.3){
        if(peak12_not_close & length(potential_peak)>1 & (potential_peak[2]-potential_peak[1])<10){
          #only keep higher peak
          if(smoothed_log[i,potential_peak[1]]<smoothed_log[i,potential_peak[2]]){
            potential_peak=potential_peak[-1]
          }else{
            potential_peak=potential_peak[-2]
          }
          #only keep lower start_of_increase
          if(smoothed_log[i,start_of_increase[1]]<smoothed_log[i,start_of_increase[2]]){
            start_of_increase=start_of_increase[-2]
          }else{
            start_of_increase=start_of_increase[-1]
          }
          height_increase=height_increase[-1]
          height_increase[1]=smoothed_log[i,potential_peak[1]]-smoothed_log[i,start_of_increase[1]]
        }
        
        if(length(potential_peak)>dim(peak)[2]){
          peak=cbind(peak,matrix(NA,dim(smoothed_log)[1],length(potential_peak)-dim(peak)[2]))
          peak_height=cbind(peak_height,matrix(NA,dim(smoothed_log)[1],length(potential_peak)-dim(peak)[2]))
        }
        peak[i,1:length(potential_peak)]=potential_peak
        peak_height[i,1:length(potential_peak)]=height_increase
        peak1_domain[i,]=c(start_of_increase[1],potential_peak[1])
      }
    }
  }
  na_col=which(apply(peak, 2, function(x) sum(!is.na(x)))==0)
  if(length(na_col)>0) {
    peak=as.matrix(peak[,-na_col])
    peak_height=as.matrix(peak_height[,-na_col])
  }

  return(list(peak=peak,peak_height=peak_height,peak1_domain=peak1_domain,numeri_grad_record=numeri_grad_record,numeri_grad2_record=numeri_grad2_record))
}


update_peak_with_small_peak_removed=function(peak_res,small_peak_ratio_cut){
  n_all=dim(peak_res$peak)[1]
  n_peak_max=dim(peak_res$peak)[2]
  for(i_sample in 1:n_all){
    peak_height_i=as.numeric(na.omit(peak_res$peak_height[i_sample,]))
    peak_i=as.numeric(na.omit(peak_res$peak[i_sample,]))
    if(length(peak_height_i)>1){
      for(j_peak in length(peak_height_i):2){
        if(peak_height_i[j_peak]/peak_height_i[j_peak-1]>small_peak_ratio_cut){
          peak_height_i=peak_height_i[-c(j_peak-1)]
          peak_i=peak_i[-c(j_peak-1)]
        }
      }
      peak_res$peak[i_sample,]=rep(NA,n_peak_max)
      peak_res$peak[i_sample,1:length(peak_i)]=peak_i
      peak_res$peak_height[i_sample,]=rep(NA,n_peak_max)
      peak_res$peak_height[i_sample,1:length(peak_i)]=peak_height_i
    }
  }
  return(peak_res)
}

rerun_KF_peak1=function(q_all,log_I,smoothed_log,peak1_domain,ind_optim){
  peak1_domain_i=peak1_domain[ind_optim,1]:peak1_domain[ind_optim,2]
  fgasp.model_rerun=fgasp(input=q_all[peak1_domain_i], output=log_I[ind_optim,peak1_domain_i])
  est_rerun=optim(c(log(1),log(.02)),log_lik,object=fgasp.model_rerun,method="L-BFGS-B",
                control = list(fnscale=-1))  
  #smoothed_sig_1peak=list()
  for(i in 1:dim(log_I)[1]){
    if(!is.na(peak1_domain[i,1])){
      peak1_domain_i=peak1_domain[i,1]:peak1_domain[i,2]
      fgasp.model=fgasp(input=q_all[peak1_domain_i], output=log_I[i,peak1_domain_i])
      pred.model=predict(param=est_rerun$par,object=fgasp.model,testing_input=q_all[peak1_domain_i])
      smoothed_log[i,peak1_domain_i]=pred.model@mean
      #smoothed_sig_1peak[[i]]=pred.model@mean
    }
  }
  return(smoothed_log)
}


peak_add_shoulder=function(peak,smoothed_log,peak1_domain,delta_q,update_grad=F){
  add_shoulder_rec=rep(F,dim(smoothed_log)[1])
  peak=cbind(peak,rep(NA,dim(smoothed_log)[1]))
  if(update_grad){
    numeri_grad_record=matrix(NA,dim(smoothed_log)[1],dim(smoothed_log)[2]-1)
    numeri_grad2_record=matrix(NA,dim(smoothed_log)[1],dim(smoothed_log)[2]-2)
  }
  #add_shoulder=rep(NA,dim(smoothed_log)[1])
  for(i in 1:dim(smoothed_log)[1]){
    if(!is.na(peak1_domain[i,1])){
      peak1_domain_i=peak1_domain[i,1]:peak1_domain[i,2]
      peak1_log=smoothed_log[i,peak1_domain_i]
      #peak1_grad=numeri_grad_record[i,peak1_domain_i[-length(peak1_domain_i)]]
      diff_log_I=diff(smoothed_log[i,])#smoothed_log[i,2:dim(smoothed_log)[2]]-smoothed_log[i,1:(dim(smoothed_log)[2]-1)]
      numeri_grad_i=diff_log_I/delta_q
      peak1_grad=numeri_grad_i[peak1_domain_i[-length(peak1_domain_i)]]
      sign_peak1_grad=sign(peak1_grad)
      
      if(update_grad){
        numeri_grad_record[i,]=numeri_grad_i
        diff_log_I2=diff(diff_log_I)
        numeri_grad2_record[i,]=diff_log_I2/(delta_q^2)
      }

      if(sum(sign_peak1_grad)!=length(sign_peak1_grad)){#if -1 exists
        potential_shoulder=which(diff(sign_peak1_grad)<0)+1
        potential_shoulder_log=smoothed_log[i,peak1_domain_i[potential_shoulder]]
        shoulder_keep=abs(potential_shoulder_log-smoothed_log[i,peak1_domain[i,2]])>.3 & abs(potential_shoulder_log-smoothed_log[i,peak1_domain[i,1]])>.3
        potential_shoulder=potential_shoulder[shoulder_keep]
        potential_shoulder_log=potential_shoulder_log[shoulder_keep]
        
        #the new smoothed log is more wiggly, we set a cutoff for the width of the negative gradient
        if(length(potential_shoulder)>0){
          negative_length_rec=rep(NA,length(potential_shoulder))
          for(i_shoulder in 1:length(potential_shoulder)){
            negative_length=1
            negative_end=potential_shoulder[i_shoulder]
            while(sign_peak1_grad[negative_end+1]==-1){
              negative_length=negative_length+1
              negative_end=negative_end+1
            }
            negative_length_rec[i_shoulder]=negative_length
          }
          potential_shoulder=potential_shoulder[which(negative_length_rec>2)] #cutoff = 2?
          potential_shoulder_log=potential_shoulder_log[which(negative_length_rec>2)]
        }
        
        if(length(potential_shoulder)>0){
          add_shoulder_rec[i]=T
          if(length(potential_shoulder)>1){ #use the smallest second deriv as the shoulder
            potential_shoulder=potential_shoulder[which.min(peak1_grad[potential_shoulder]-peak1_grad[potential_shoulder-1])]
          }
          peak[i,]=c(peak1_domain_i[potential_shoulder],peak[i,-dim(peak)[2]])
        }
      }
    }
  }
  if(sum(!is.na(peak[,dim(peak)[2]]))==0){
    peak=peak[,-dim(peak)[2]]
  }
  res=list(peak=peak,add_shoulder_rec=add_shoulder_rec)
  if(update_grad){
    res$numeri_grad_record=numeri_grad_record
    res$numeri_grad2_record=numeri_grad2_record
  }
  return(res)
}

plot_by_label=function(label,all_labels,plot_path,plot_this,...){
  for(cat in 1:length(all_labels)){
    cat_ind=which(label==all_labels[cat])
    if(all_labels[cat]=='HEX/LAM'){
      pdf(file = paste0(plot_path,'HEX_LAM',"_1.pdf"), width = 9, height = 6)
    }else{
      pdf(file = paste0(plot_path,all_labels[cat],"_1.pdf"), width = 9, height = 6)
    }
    par(mfrow=c(3,4), mgp=c(2,1,0), mar=c(3,3,1,1)+.1)
    
    for(i in 1:length(cat_ind)){
      plot_i=cat_ind[i]
      plot_this(plot_i=plot_i,...)
      if(i %% 12==0 & i!=length(cat_ind)){
        dev.off()
        pdf(file = paste0(plot_path,all_labels[cat],"_",i %/% 12+1,".pdf"), width = 9, height = 6)
        par(mfrow=c(3,4), mgp=c(2,1,0), mar=c(3,3,1,1)+.1)
      }
    }
    dev.off()
  }
  
  cat_ind=which(is.na(label))
  if(length(cat_ind)>0){
    pdf(file = paste0(plot_path,'NA',"_1.pdf"), width = 9, height = 6)
    par(mfrow=c(3,4), mgp=c(2,1,0), mar=c(3,3,1,1)+.1)
    for(i in 1:length(cat_ind)){
      plot_i=cat_ind[i]
      plot_this(plot_i=plot_i,...)
      if(i %% 12==0 & i!=length(cat_ind)){
        dev.off()
        pdf(file = paste0(plot_path,'NA',"_",i %/% 12+1,".pdf"), width = 9, height = 6)
        par(mfrow=c(3,4), mgp=c(2,1,0), mar=c(3,3,1,1)+.1)
      }
    }
    dev.off()
  }
  
}


get_grad_data=function(numeri_grad_record,numeri_grad2_record,all_peak){
  peak1_grad2_value=rep(NA,dim(numeri_grad2_record)[1])
  peak1_grad2_width_to_zero=matrix(NA,dim(numeri_grad2_record)[1],3,dimnames = list(NULL,c('left','right','width')))
  #peak1_up_down_width=matrix(NA,dim(numeri_grad2_record)[1],3,dimnames = list(NULL,c('left','right','width')))

  for(i in 1:dim(numeri_grad2_record)[1]){
    peak1_grad2_value[i]=numeri_grad2_record[i,all_peak[i,1]-1]
    #peak1_grad2_value[i]=min(numeri_grad2_record[i,all_peak[i,1]],numeri_grad2_record[i,all_peak[i,1]-1],numeri_grad2_record[i,all_peak[i,1]+1])
    
    ##get peak1_grad2_width_to_zero
    left_index=all_peak[i,1]-1
    right_index=all_peak[i,1]-1
    value_left=numeri_grad2_record[i,left_index]
    value_right=numeri_grad2_record[i,right_index]
    if(!is.na(value_left)){
      while(value_left<0 & left_index>1){
        left_index=left_index-1
        value_left=numeri_grad2_record[i,left_index]
      }
      while(value_right<0 & right_index<dim(numeri_grad2_record)[2]){
        right_index=right_index+1
        value_right=numeri_grad2_record[i,right_index]
      }
      peak1_grad2_width_to_zero[i,]=c(left_index,right_index,right_index-left_index)
    }
    
    # ##get peak1_up_down_width
    # sign_grad_i=sign(numeri_grad_record[i,])
    # left_index=all_peak[i,1]-1
    # right_index=all_peak[i,1]-1
    # value_left=sign_grad_i[left_index+1]-sign_grad_i[left_index]
    # value_right=sign_grad_i[right_index+1]-sign_grad_i[right_index]
    # if(!is.na(numeri_grad2_record[i,left_index])){
    #   while(value_left<=0 & left_index>1){
    #     left_index=left_index-1
    #     value_left=sign_grad_i[left_index+1]-sign_grad_i[left_index]
    #   }
    #   while(value_right<=0 & right_index<dim(numeri_grad_record)[2]){
    #     right_index=right_index+1
    #     value_right=sign_grad_i[right_index+1]-sign_grad_i[right_index]
    #   }
    #   peak1_up_down_width[i,]=c(left_index,right_index,right_index-left_index)
    # }
  }
  return(list(peak1_grad2_value=peak1_grad2_value,peak1_grad2_width_to_zero=peak1_grad2_width_to_zero#,
              #peak1_up_down_width=peak1_up_down_width
              ))
}

plot_label_overlay=function(plot_vec,label,y_lab){
  na_ind=is.na(label)
  plot_data=data.frame(index=1:length(plot_vec[!na_ind]),plot_vec=plot_vec[!na_ind],label=label[!na_ind])
  ggplot(plot_data,aes(index,plot_vec))+geom_point(aes(color=label,shape=label),na.rm=T)+
    ylab(y_lab)+
    scale_colour_discrete(na.translate = F)+
    scale_shape_manual(values = 1:length(unique(plot_data$label)),na.translate = F)+
    theme_bw()
}



get_data=function(model_peak,model_label,n_peak,
                  add_ratios=T,add_peak1_grad2=F,add_peak1_grad2_width=F,#add_peak1_up_down_width=F,
                  add_frac=F,add_weight=F,add_weight_D=F,add_weight_F=F,add_chem=F,
                  peak1_grad2_value=NA,peak1_grad2_width_to_zero=NA,#peak1_up_down_width=NA,
                  frac=NA,weight=NA,weight_D=NA,weight_F=NA,chemistry=NA){
  #label_sub=label[-rm_ind]
  model_label=factor(model_label) #remove the level HEX/LAM in the factor
  #peak_sub=as.matrix(all_peak[-rm_ind,])
  
  model_data=NULL
  if(add_ratios){
    peak_ratio=model_peak[,2:n_peak]/model_peak[,1]
    colnames(peak_ratio)=paste("r",2:n_peak,"1",sep="")
    model_data=cbind(model_data,p1=model_peak[,1],peak_ratio)
    if(n_peak>2){
      for(i in 3:n_peak){
        new_ratio=as.matrix(model_peak[,i:n_peak]/model_peak[,(i-1)])
        colnames(new_ratio)=paste("r",i:n_peak,(i-1),sep="")
        model_data=cbind(model_data,new_ratio)
      }
    }
  }
  if(add_peak1_grad2) model_data=cbind(model_data,peak1_grad2=peak1_grad2_value)
  if(add_peak1_grad2_width) model_data=cbind(model_data,peak1_grad2_width=peak1_grad2_width_to_zero)
  #if(add_peak1_up_down_width)  model_data=cbind(model_data,peak1_up_down_width=peak1_up_down_width[-rm_ind])
  if(add_frac) model_data=cbind(model_data,frac=frac)
  if(add_weight) model_data=cbind(model_data,weight=weight)
  if(add_weight_D) model_data=cbind(model_data,weight_D=weight_D)
  if(add_weight_F) model_data=cbind(model_data,weight_F=weight_F)
  if(add_chem) model_data=cbind(model_data,chemistry=chemistry)
  
  return(list(model_data=model_data,model_label=model_label))
}

get_folds=function(nfold,size){
  folds = cut(1:size, breaks=nfold, labels=FALSE) %>% sample()
  return(folds)
}

classification_CV=function(class_data,class_label,folds,model_used=c('bagging','rforest','xgboost')){
  nfold=max(folds)
  N = length(class_label)
  pred_prob=matrix(NA,N,length(levels(class_label)),dimnames = list(NULL,levels(class_label)))
  class_pred_label=factor(rep(NA,N),levels=levels(class_label))
  
 
  if(model_used!='xgboost') class_data[which(is.na(class_data))]=0
  if(model_used=='xgboost') xgb_params = list(#https://xgboost.readthedocs.io/en/latest/parameter.html
    booster = "gbtree",
    eta = 0.1,
    max_depth = 5,
    gamma = 0,
    subsample = 0.75,
    colsample_bytree = 1,
    objective = "multi:softprob",
    eval_metric = "mlogloss",
    num_class = length(levels(class_label))
  )

  for(i in 1:nfold){
    train_ind=which(folds!=i)
    test_ind=which(folds==i)
    
    if(model_used=='bagging'){
      m_data = randomForest(x=class_data[train_ind,],y=class_label[train_ind],mtry=ncol(class_data), importance=T)
    }else if(model_used=='rforest'){
      m_data = randomForest(x=class_data[train_ind,],y=class_label[train_ind],mtry=round(sqrt(ncol(class_data))), importance=T)
    }else if(model_used=='xgboost'){
      label_train <- as.integer(class_label[train_ind]) - 1
      label_test <- as.integer(class_label[test_ind]) - 1
      xgboost_train = xgb.DMatrix(data=class_data[train_ind,], label=label_train)
      #m_data= xgboost(data = xgboost_train,max.depth=3,nrounds=50,verbose=0)                   
      m_data = xgb.train(params = xgb_params,data = xgboost_train,nrounds = 500,verbose = 0)
    }
    if(model_used=='xgboost'){
      xgboost_test = xgb.DMatrix(data=class_data[test_ind,])#, label=class_label[test_ind])
      m_pred = predict(m_data, xgboost_test,reshape=T)
      colnames(m_pred) = levels(class_label)
      pred_prob[test_ind,]=m_pred
      class_pred_label[test_ind]=factor(apply(m_pred,1,function(x) colnames(m_pred)[which.max(x)]),levels=levels(class_label))
    }else{
      m_pred = predict(m_data, newdata = class_data[test_ind,], type = "prob")
      pred_prob[test_ind,]=m_pred
      class_pred_label[test_ind]=factor(apply(m_pred,1,function(x) colnames(m_pred)[which.max(x)]),levels=levels(class_label))
    }
  }
  
  return(list(class_pred_label=class_pred_label,N=N,pred_prob=pred_prob))
}

get_CV_accuracy=function(model_label,pred_label,folds){
  nfold=max(folds)
  N = length(model_label)
  accuracy_overall=1-sum(pred_label!=model_label)/N
  accuracy_fold=rep(NA,nfold)
  for(i in 1:nfold){
    test_ind=which(folds==i)
    accuracy_fold[i]=1-sum(pred_label[test_ind]!=model_label[test_ind])/length(test_ind)
  }
  accuracy_average=mean(accuracy_fold)
  
  #create confusion matrix for the best performance group (the third column)
  conf=confusionMatrix(data=factor(pred_label,levels = levels(model_label)),reference=model_label)
  conf_mat=cbind(t(conf$table), conf$byClass[,"Sensitivity"])
  dimnames(conf_mat)=list(true=levels(model_label),pred=c(levels(model_label),"Accuracy"))
  
  return(list(accuracy_overall=accuracy_overall,accuracy_average=accuracy_average,accuracy_fold=accuracy_fold,conf_mat=conf_mat))
}

get_conf_mat=function(model_label,pred_label){
  if(length(unique(model_label))!=length(levels(model_label)) & length(unique(pred_label))!=length(levels(pred_label))){
    levels_all=sort(unique(c(model_label,pred_label)))
    model_label=factor(model_label,levels=levels_all)
    pred_label=factor(pred_label,levels = levels_all)
  }
  conf=confusionMatrix(data=factor(pred_label,levels = levels(model_label)),reference=model_label)
  conf_mat=cbind(t(conf$table), conf$byClass[,"Sensitivity"])
  dimnames(conf_mat)=list(true=levels(model_label),pred=c(levels(model_label),"Accuracy"))
  return(conf_mat)
}

classification_train_test=function(train_input,train_output,test_input,test_output,model_used=c('bagging','rforest','xgboost')){
  n_test = length(test_output)
  #pred_prob=matrix(NA,n_test,length(levels(train_output)),dimnames = list(NULL,levels(train_output)))
  class_pred_label=factor(rep(NA,n_test),levels=levels(train_output))
  
  if(model_used!='xgboost') {
    train_input[which(is.na(train_input))]=0
    test_input[which(is.na(test_input))]=0
  }
  if(model_used=='xgboost') xgb_params = list(#https://xgboost.readthedocs.io/en/latest/parameter.html
    booster = "gbtree",
    eta = 0.1,
    max_depth = 5,
    gamma = 0,
    subsample = 0.75,
    colsample_bytree = 1,
    objective = "multi:softprob",
    eval_metric = "mlogloss",
    num_class = length(levels(train_output))
  )
  
  if(model_used=='bagging'){
    m_data = randomForest(x=train_input,y=train_output,mtry=ncol(train_input), importance=T)
  }else if(model_used=='rforest'){
    m_data = randomForest(x=train_input,y=train_output,mtry=round(sqrt(ncol(train_input))), importance=T)
  }else if(model_used=='xgboost'){
    label_train <- as.integer(train_output) - 1
    label_test <- as.integer(test_output) - 1
    xgboost_train = xgb.DMatrix(data=train_input, label=label_train)
    #m_data= xgboost(data = xgboost_train,max.depth=3,nrounds=50,verbose=0)                   
    m_data = xgb.train(params = xgb_params,data = xgboost_train,nrounds = 500,verbose = 0)
  }
  if(model_used=='xgboost'){
    xgboost_test = xgb.DMatrix(data=test_input)#, label=test_output)
    m_pred = predict(m_data, xgboost_test,reshape=T)
    colnames(m_pred) = levels(train_output)
    #pred_prob=m_pred
    class_pred_label=factor(apply(m_pred,1,function(x) colnames(m_pred)[which.max(x)]),levels=levels(train_output))
  }else{
    m_pred = predict(m_data, newdata = test_input, type = "prob")
    #pred_prob=m_pred
    class_pred_label=factor(apply(m_pred,1,function(x) colnames(m_pred)[which.max(x)]),levels=levels(train_output))
  }

  return(list(class_pred_label=class_pred_label,pred_prob=m_pred,model=m_data))
}

CV_for_cutoff=function(peak_height_cutoff,nfold,train_ind,model_used,smooth_log_train,label_train,frac_train=NULL,weight_train=NULL,weight_D_train=NULL,weight_F_train=NULL){
  peak12_not_close=T
  peak_res_hetero=peak_detect(q_all=q_all,smoothed_log=smooth_log_train,peak_height_cutoff=peak_height_cutoff,peak12_not_close=peak12_not_close)#,peak_slope_cutoff=peak_slope_cutoff)
  ##get information of gradient
  grad_data_hetero=get_grad_data(peak_res_hetero$numeri_grad_record,peak_res_hetero$numeri_grad2_record,all_peak=peak_res_hetero$peak)
  #get data
  n_peak=3
  #settings for the data
  add_ratios=T
  add_peak1_grad2=T
  add_peak1_grad2_width=T
  add_frac=F
  add_weight=F
  add_weight_D=F
  add_weight_F=F
  dat=get_data(model_peak=peak_res_hetero$peak,model_label=label_train,n_peak=n_peak,
               add_ratios,add_peak1_grad2,add_peak1_grad2_width,add_frac,add_weight,add_weight_D,add_weight_F,
               peak1_grad2_value=grad_data_hetero$peak1_grad2_value,peak1_grad2_width_to_zero=grad_data_hetero$peak1_grad2_width_to_zero[,3],
               frac=frac_train,weight=weight_train,weight_D=weight_D_train,weight_F=weight_F_train)
  #cross-validation
  folds_all=get_folds(nfold=nfold,size=length(dat$model_label))
  #two step prediction
  have_peak_ind=which(apply(peak_res_hetero$peak,1,function(x) sum(!is.na(x)))>0)
  ##for samples without peaks detected, predict as DIS
  pred_label=factor(rep(NA,length(dat$model_label)),levels=levels(dat$model_label))
  pred_label[-have_peak_ind]='DIS'
  #random_forest
  rf_ratio_hetero=classification_CV(class_data=dat$model_data[have_peak_ind,],class_label=dat$model_label[have_peak_ind],folds=folds_all[have_peak_ind],model_used=model_used)
  pred_label[have_peak_ind]=rf_ratio_hetero$class_pred_label
  accuracy_fold=rep(NA,nfold)
  for(i in 1:nfold){
    test_ind=which(folds_all==i)
    accuracy_fold[i]=1-sum(pred_label[test_ind]!=dat$model_label[test_ind])/length(test_ind)
  }
  accuracy_average=mean(accuracy_fold)
  
  print(c(peak_height_cutoff,accuracy_average))
  return(accuracy_average)
}


get_theta_hat=function(param,input,output,kernel_type,homo=T,error_var=NULL){
  beta=exp(param[1])
  tau=exp(param[2])
  R0=abs(outer(input,(input),'-'))
  if(kernel_type=='exp'){
    R = exp(-beta*R0)
  }else if(kernel_type=='matern_5_2'){
    R = matern_5_2_funct(R0, beta)
  }
  if(homo){
    R_tilde=R+tau*diag(length(input))
  }else{
    R_tilde=R+tau*diag(error_var)
  }
  L_tilde=t(chol(R_tilde))
  X=matrix(rep(1,length(input)))
  R_inv_y=as.matrix(backsolve(t(L_tilde), forwardsolve(L_tilde,output)))
  R_inv_X=as.matrix(backsolve(t(L_tilde), forwardsolve(L_tilde,X)))
  theta_hat=t(X)%*%R_inv_y/(t(X)%*%R_inv_X)
  
  return(theta_hat)
}

get_theta_hat_with_KF=function(param,input,output,kernel_type,homo=T,error_var=NULL){
  beta=exp(param[1])
  gamma=1/beta
  VV=exp(param[2])
  
  delta_x=input[2:p1]-input[1:(p1-1)]
  
  if(kernel_type=="matern_5_2"){
    lambda=sqrt(5.0)/gamma;
    W0=Construct_W0_matern_5_2(1.0,lambda);   
    GG=Construct_G_matern_5_2(delta_x,lambda);  
    W=Construct_W_matern_5_2(1.0,delta_x,lambda,W0);
    
  }else if(kernel_type=="exp"){
    lambda=1.0/gamma;
    W0=Construct_W0_exp(1.0,lambda);  
    GG=Construct_G_exp(delta_x,lambda);  
    W=Construct_W_exp(1.0,delta_x,lambda,W0);
  }
  if(homo){
    Q_K=Get_Q_K(GG,W,W0,VV)
  }else{
    Q_K=Get_Q_K_hetero(GG,W,W0,VV,error_var=error_var)
  }
  
  Q=Q_K[[1]]
  K=Q_K[[2]]
  
  X=matrix(rep(1,length(input)))
  R_inv_y=R_inv_times_z_pre_constructed(z=output, sqrt_Q=sqrt(Q), K=K, GG=GG)
  R_inv_X=R_inv_times_z_pre_constructed(z=X, sqrt_Q=sqrt(Q), K=K, GG=GG)
  theta_hat=t(X)%*%R_inv_y/(t(X)%*%R_inv_X)
  
  return(theta_hat)
}

get_W0_GG_W = function(input,gamma,kernel_type){
  delta_x=input[2:p1]-input[1:(p1-1)]
  
  if(kernel_type=="matern_5_2"){
    lambda=sqrt(5.0)/gamma;
    W0=Construct_W0_matern_5_2(1.0,lambda);   
    GG=Construct_G_matern_5_2(delta_x,lambda);  
    W=Construct_W_matern_5_2(1.0,delta_x,lambda,W0);
    
  }else if(kernel_type=="exp"){
    lambda=1.0/gamma;
    W0=Construct_W0_exp(1.0,lambda);  
    GG=Construct_G_exp(delta_x,lambda);  
    W=Construct_W_exp(1.0,delta_x,lambda,W0);
  }
  
  return(list(W0=W0,GG=GG,W=W))
}

smooth_intensity_with_mean_fix_param=function(q_all,log_I,param,theta_hat,kernel_type,homo=T,log_var=NULL){
  beta=exp(param[1])
  tau=exp(param[2])
  
  smoothed_log=matrix(NA,dim(log_I)[1],dim(log_I)[2])
  
  R0=abs(outer(input,(input),'-'))
  if(kernel_type=='exp'){
    R = exp(-beta*R0)
  }else if(kernel_type=='matern_5_2'){
    R = matern_5_2_funct(R0, beta)
  }
  
  X=matrix(rep(1,length(q_all)))
  X_theta=X%*%theta_hat
  
  if(homo){
    R_tilde=R+tau*diag(length(input))
    L_tilde=t(chol(R_tilde))
  }
  
  for(i in 1:dim(log_I)[1]){
    if(!homo){
      R_tilde=R+tau*diag(log_var[i,])
      L_tilde=t(chol(R_tilde))
    }
    
    R_inv_res=as.matrix(backsolve(t(L_tilde), forwardsolve(L_tilde,log_I[i,]-X_theta)))
    smoothed_log[i,]=X_theta+R %*% R_inv_res
  }

  return(smoothed_log)

}

smooth_intensity_with_mean_fix_param_with_KF=function(q_all,log_I,param,theta_hat,kernel_type,homo=T,log_var=NULL){
  gamma=1/exp(param[1])
  VV=exp(param[2])
  
  smoothed_log=matrix(NA,dim(log_I)[1],dim(log_I)[2])
  
  W0_GG_W=get_W0_GG_W(input=q_all,gamma=gamma,kernel_type=kernel_type)
  
  W0=W0_GG_W[[1]]
  GG=W0_GG_W[[2]]
  W=W0_GG_W[[3]]
  
  X=matrix(rep(1,length(q_all)))
  X_theta=X%*%theta_hat
  
  if(homo){
    Q_K=Get_Q_K(GG=GG,W=W,C0=W0,VV=VV)
  }
  
  for(i in 1:dim(log_I)[1]){
    if(!homo){
      Q_K=Get_Q_K_hetero(GG=GG,W=W,C0=W0,VV=VV,error_var=log_var[i,])
    }
    
    sqrt_Q=sqrt(Q_K[[1]])
    K=Q_K[[2]]
    
    R_inv_res=R_inv_times_z_pre_constructed(z=log_I[i,]-X_theta, sqrt_Q=sqrt_Q, K=K, GG=GG)
    if(homo){
      smoothed_log[i,]=X_theta+R_times_z_pre_constructed(z=R_inv_res, sqrt_Q=sqrt_Q, K=K, GG=GG)-VV*R_inv_res
    }else{
      smoothed_log[i,]=X_theta+R_times_z_pre_constructed(z=R_inv_res, sqrt_Q=sqrt_Q, K=K, GG=GG)-VV*log_var[i,]*R_inv_res
    }

  }
  
  return(smoothed_log)
}

