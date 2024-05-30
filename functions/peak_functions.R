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
    
    #remove the peak outside the truncation range, which are mostly noise
    keep=q_all[potential_peak]<=truncate_range[2] & q_all[potential_peak]>truncate_range[1]
    #keep=potential_peak<=400 & potential_peak>10
    potential_peak=potential_peak[keep]
    start_of_increase=start_of_increase[keep]
    end_of_decrease=end_of_decrease[keep] 
    
    #remove the super flat peak
    #keep=apply(cbind(start_of_increase,end_of_decrease),1, function(x) max(abs(diff_log_I2[x[1]:x[2]-1]))>10^(-3))
    keep=abs(diff_log_I2[potential_peak-1])>10^(-3)
    potential_peak=potential_peak[keep]
    start_of_increase=start_of_increase[keep]
    end_of_decrease=end_of_decrease[keep] 
    
    #redefine start_of_increase and end_of_decrease based on smoothed diff_log_I2
    if(length(potential_peak)>0){
      for(t in length(potential_peak):1){
        if(potential_peak[t]-start_of_increase[t]<=1 | end_of_decrease[t]-potential_peak[t] <=1){
          potential_peak=potential_peak[-t]
          start_of_increase=start_of_increase[-t]
          end_of_decrease=end_of_decrease[-t] 
        }else{
          left_side_ind=which(abs(filter(diff_log_I2[start_of_increase[t]:potential_peak[t]-1],rep(1 / 3, 3),sides=2))>10^(-3)|
                                abs(filter(diff_log_I[start_of_increase[t]:potential_peak[t]-1],rep(1 / 3, 3),sides=2))>10^(-3))
          right_side_ind=which(abs(filter(diff_log_I2[potential_peak[t]:end_of_decrease[t]-1],rep(1 / 3, 3),sides=2))>10^(-3) |
                                 abs(filter(diff_log_I[potential_peak[t]:end_of_decrease[t]-1],rep(1 / 3, 3),sides=2))>10^(-3))
          if(length(left_side_ind)==0 | length(right_side_ind)==0){
            potential_peak=potential_peak[-t]
            start_of_increase=start_of_increase[-t]
            end_of_decrease=end_of_decrease[-t] 
          }else{
            start_of_increase[t]=c(start_of_increase[t]:potential_peak[t])[min(left_side_ind)]
            end_of_decrease[t]=c(potential_peak[t]:end_of_decrease[t])[max(right_side_ind)]
          }
        }
      }
    }
    
    
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
      height_decrease=smoothed_log[i,(potential_peak)]-smoothed_log[i,end_of_decrease]
      height_peak=pmax(height_increase,height_decrease)
      # decrease_width=end_of_decrease-potential_peak
      # #remove the potential peaks with width < 3
      # smooth_out_int=which(decrease_width<3)
      # if(length(smooth_out_int)>0){
      #   for(i_smooth in rev(smooth_out_int)){
      #     if(i_smooth<length(potential_peak)){
      #       start_of_increase=start_of_increase[-(i_smooth+1)]
      #       potential_peak=potential_peak[-i_smooth]
      #       end_of_decrease=end_of_decrease[-i_smooth]
      #     }
      #   }
      #   height_increase=smoothed_log[i,(potential_peak)]-smoothed_log[i,start_of_increase]
      #   #height_decrease=smoothed_log[i,(potential_peak)]-smoothed_log[i,end_of_decrease]
      #   decrease_width=end_of_decrease-potential_peak
      # }
      
      
      #grad2=numeri_grad2_record[i,potential_peak-1]
      #peak_slope=height_increase/((potential_peak-start_of_increase)*delta_q)
      index_peak=which(height_peak>peak_height_cutoff)#  & peak_slope>peak_slope_cutoff)# & grad2<grad2_cutoff) 
      
      if(length(index_peak)>0){
        potential_peak=potential_peak[index_peak]
        start_of_increase=start_of_increase[index_peak]
        end_of_decrease=end_of_decrease[index_peak]
        #height_increase=height_increase[index_peak]
        
        #the first two peaks can't have close log(I)
        if(peak12_not_close & length(potential_peak)>1){
          for(t in length(potential_peak):2){
            if((potential_peak[t]-potential_peak[t-1])<10){# &diff(smoothed_log[i,c(potential_peak[t-1],potential_peak[t])])<diff(range(smoothed_log[i,]))/20
              #only keep higher peak
              if(smoothed_log[i,potential_peak[t-1]]<smoothed_log[i,potential_peak[t]]){
                potential_peak=potential_peak[-(t-1)]
              }else{
                potential_peak=potential_peak[-t]
              }
              #only keep lower start_of_increase
              if(smoothed_log[i,start_of_increase[t-1]]<smoothed_log[i,start_of_increase[t]]){
                start_of_increase=start_of_increase[-t]
              }else{
                start_of_increase=start_of_increase[-(t-1)]
              }
              #only keep lower end_of_decrease
              if(smoothed_log[i,end_of_decrease[t-1]]<smoothed_log[i,end_of_decrease[t]]){
                end_of_decrease=end_of_decrease[-t]
              }else{
                end_of_decrease=end_of_decrease[-(t-1)]
              }
            }
          }
        }
        height_increase=smoothed_log[i,(potential_peak)]-smoothed_log[i,start_of_increase]
        height_decrease=smoothed_log[i,(potential_peak)]-smoothed_log[i,end_of_decrease]
        height_peak=pmax(height_increase,height_decrease)
        
        if(length(potential_peak)>dim(peak)[2]){
          peak=cbind(peak,matrix(NA,dim(smoothed_log)[1],length(potential_peak)-dim(peak)[2]))
          peak_height=cbind(peak_height,matrix(NA,dim(smoothed_log)[1],length(potential_peak)-dim(peak)[2]))
        }
        peak[i,1:length(potential_peak)]=potential_peak
        peak_height[i,1:length(potential_peak)]=height_peak
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
      if(length(peak_height_i)>2){
        for(j_peak in (length(peak_height_i)-1):2){
          if(peak_height_i[j_peak+1]/peak_height_i[j_peak]>small_peak_ratio_cut &
             (peak_i[j_peak+1]-peak_i[j_peak])/(peak_i[j_peak+1]-peak_i[j_peak-1])<1/2.2){
            peak_height_i=peak_height_i[-c(j_peak)]
            peak_i=peak_i[-c(j_peak)]
          }else if(peak_height_i[j_peak-1]/peak_height_i[j_peak]>small_peak_ratio_cut &
                   (peak_i[j_peak]-peak_i[j_peak-1])/(peak_i[j_peak+1]-peak_i[j_peak-1])<1/2.2){
            peak_height_i=peak_height_i[-c(j_peak)]
            peak_i=peak_i[-c(j_peak)]
          }
        }
      }

      j_peak=1
      if(peak_height_i[j_peak+1]/peak_height_i[j_peak]>small_peak_ratio_cut){
        peak_height_i=peak_height_i[-c(j_peak)]
        peak_i=peak_i[-c(j_peak)]
      }
      peak_res$peak[i_sample,]=rep(NA,n_peak_max)
      peak_res$peak[i_sample,1:length(peak_i)]=peak_i
      peak_res$peak_height[i_sample,]=rep(NA,n_peak_max)
      peak_res$peak_height[i_sample,1:length(peak_i)]=peak_height_i
    }
  }
  return(peak_res)
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
    value_left=mean(numeri_grad2_record[i,left_index+(-1:1)])
    value_right=mean(numeri_grad2_record[i,right_index+(-1:1)])
    if(!is.na(value_left)){
      while(value_left<(10^(-3))/delta_q^2 & left_index>2){
        left_index=left_index-1
        value_left=mean(numeri_grad2_record[i,left_index+(-1:1)])
      }
      while(value_right<(10^(-3))/delta_q^2 & right_index<dim(numeri_grad2_record)[2]-1){
        right_index=right_index+1
        value_right=mean(numeri_grad2_record[i,right_index+(-1:1)])
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


