##likelihood function of GaSP using fast computation
##the Filtering algorithm is implemented using Rcpp and RcppEigen
log_lik_hetero<-function(param,object,error_var){
  param=as.vector(param)
  if(object@have_noise==T && length(param)!=2){
       stop("Please specify both the log inverse range parameter 
            and log nugget parameter if there is noise. \n")
  }
  
  log_det_S2_hetero=Get_log_det_S2_hetero(param,object@have_noise,object@delta_x,
                                   object@output, object@kernel_type,error_var);
  ##log likelihood
  -log_det_S2_hetero[[1]]/2-(object@num_obs)/2*log(log_det_S2_hetero[[2]])
}


##prediction 
predict.fgasp_hetero<-function(param,object,testing_input,var_data=TRUE,sigma_2=NULL,error_var){
  param=as.vector(param)
  if(object@have_noise==T && length(param)!=2){
    stop("Please specify both the log inverse range parameter 
          and log nugget parameter. \n")
  }
  ## if sigma_2 is not given, I will estimate it
  if(is.null(sigma_2)){
    sigma_2=Get_log_det_S2_hetero(param,object@have_noise,object@delta_x,
                           object@output,object@kernel_type,error_var)[[2]]/object@num_obs
  }
  
  predictobj<- new("predictobj.fgasp")
  
  testing_input_sorted=sort(testing_input)
  predictobj@num_testing=length(testing_input_sorted)
  predictobj@testing_input=testing_input_sorted
  predictobj@param=param
  predictobj@var_data=var_data
  ## we don't support the repeated testing input in this version
  if(length(which(testing_input_sorted[2:predictobj@num_testing]-
                  testing_input_sorted[1:(predictobj@num_testing-1)]==0))>0){
    stop("Please delete the repeated testing inputs. \n")
  }
  ##combine the input and testing input
  input_all_ori=c(object@input,testing_input_sorted)
  input_all_sort=sort(input_all_ori,index.return=TRUE)
  
  input_all=input_all_sort[[1]] ##all sorted input with and without the observations
  delta_x_all=input_all[2:length(input_all_ori)]-input_all[1:(length(input_all_ori)-1)]
  
  ##create a vector of inputs where the one mean there is an observation
  ##and zero means no observation
  index_all=rep(0, length(input_all))
  index_all[1:object@num_obs]=1
  index_all=index_all[input_all_sort[[2]]]
  index_all=as.integer(index_all)
  
  
  
  testing_loc=which(index_all==0)
  delta_x_all_here=delta_x_all
  
  if(length(which(delta_x_all==0))>0){
    #input_all=input_all[-which(delta_x_all==0)]
    index_all=index_all[-(which(delta_x_all==0)+1)]
    delta_x_all_here=delta_x_all[-which(delta_x_all==0)]
  }
  
  
  
  KF_smoother_result=Kalman_smoother_hetero(param, object@have_noise,index_all, delta_x_all_here,
                                            object@output,sigma_2,object@kernel_type,error_var)
  if(length(which(delta_x_all==0))==0){
    predictobj@mean=KF_smoother_result[[1]][testing_loc]
    if(var_data==T | object@have_noise==F){
      predictobj@var=KF_smoother_result[[2]][testing_loc]
    }else{
      predictobj@var=KF_smoother_result[[2]][testing_loc]-sigma_2*exp(param[2])*error_var
    }
  }else{##enlarge it to the orginal long sequence that may contain knots. 
    res=rep(NA,length(input_all_ori))
    res[-(which(delta_x_all==0)+1)]=KF_smoother_result[[1]]
    res[which(delta_x_all==0)+1]= res[which(delta_x_all==0)]
    predictobj@mean=res[testing_loc]
    
    res[-(which(delta_x_all==0)+1)]=KF_smoother_result[[2]]
    res[which(delta_x_all==0)+1]= res[which(delta_x_all==0)]
    if(var_data==T | object@have_noise==F){
      predictobj@var=res[testing_loc]
    }else{
      predictobj@var=res[testing_loc]-sigma_2*exp(param[2])*error_var
    }
    
  }
  
  return(predictobj)
  
}



