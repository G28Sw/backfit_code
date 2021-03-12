require(tibble)
require(dplyr)
get_lambdas<-function(y,x,f1,f2){
  #This function estimates lambda based on the method of moments
  
  resids = lm(y~x)$residuals
  
  df=tibble(client=f1,item = f2,y_residual = resids)
  
  Ua0 = df %>% group_by(client) %>% summarise(Ua = (var(y_residual)*(length(y_residual)-1)))
  Ua0[is.na(Ua0)] <- 0
  
  Ub0 = df %>% group_by(item) %>% summarise(Ub = (var(y_residual)*(length(y_residual)-1)))
  Ub0[is.na(Ub0)] <- 0
  
  Ue0 = (var(df$y_residual)*(length(df$y_residual)-1)*length(df$y_residual))#sum((as.numeric(dist(df$rating)))^2)*0.5
  Ue0[is.na(Ue0)] <- 0
  
  Ua = sum(Ua0[,2])
  Ub = sum(Ub0[,2])
  Ue = Ue0
  
  Va0 = df %>% group_by(client) %>% summarise(Va = (length(y_residual))^2)
  Vb0 = df %>% group_by(item) %>% summarise(Vb = (length(y_residual))^2)
  
  Va = sum(Va0[,2])
  Vb = sum(Vb0[,2])
  
  N = nrow(df)
  R = length(unique(df$client))
  C = length(unique(df$item))
  
  U = c(Ua, Ub, Ue)
  
  M = matrix(rep(0,9),nrow=3)
  M[1,c(2,3)] = N-R
  M[2,c(1,3)] = N-C
  M[3,1] = N^2 - Va
  M[3,2] = N^2 - Vb
  M[3,3] = N^2-N
  
  s = solve(M, U)
  sigmasq_a =  s[1]
  sigmasq_b =  s[2]
  sigmasq_e =  s[3]
  
  lambda_a = sigmasq_e/sigmasq_a
  lambda_b = sigmasq_e/sigmasq_b
  if (lambda_a<0){
    print("lambda_a is negative")
    lambda_a = 0
  }
  if (lambda_b<0){
    print("lambda_b is negative")
    lambda_b = 0
  }
  
  return(c(lambda_a,lambda_b,sigmasq_e))
}


level_means<-function(x,f,denom){
  #This function uses tibble to do tapply on multiple columns of a matrix
  
  d=tibble(f,as_tibble(x))
  pwsum=d %>%
    group_by(f) %>%
    summarize_all(sum,na.rm=TRUE)
  num = data.matrix(pwsum[,-1])
  centering_weight = (1/denom)/sum(1/denom)
  estimate = scale(num, center= colSums(centering_weight*num),scale=F)/denom
  return(estimate)
  
}

backfitting_random_effect <-function(y,x,f1,f2,lambda_a,lambda_b,epsilon,max_iter,trace.it){
  #this function obtains blup estimates of random effects

  
  
  n_a = table(f1)
  a_hat_denom = as.vector(n_a+lambda_a)
  
  n_b = table(f2)
  b_hat_denom = as.vector(n_b+lambda_b)
  
  a_haty=rep(0,length(n_a))
  b_haty=rep(0,length(n_b))
  
  iter = 0
  relative_change = epsilon+2
  
  while((relative_change > epsilon)&(iter<=max_iter)){ 
    iter = iter+1
    deltaf=0
    
    
    residay = y - b_haty[f2]
    new_a_haty = as.numeric(level_means(residay,f1,a_hat_denom))
    
    
  
    residby = y - a_haty[f1]
    new_b_haty = as.numeric(level_means(residby,f2,b_hat_denom))

    

    deltaf=mean((new_a_haty[f1]-a_haty[f1])^2)+mean((new_b_haty[f2]-b_haty[f2])^2)
    dd = mean((a_haty[f1]+b_haty[f2])^2)
    relative_change = deltaf/dd
    if(trace.it){cat("relative change in blup estimate iterations : ",relative_change,"\n")}
    a_haty = new_a_haty
    b_haty = new_b_haty
    
  }
  
  
  l = list("a_hat"= new_a_haty, "b_hat" = new_b_haty)
  
  return(l)
  
}

backfitting_covariance_centering <-function(y,x,f1,f2,lambda_a,lambda_b,epsilon,max_iter,trace.it){
  #this function smooths out x by alternating between SA and SB using weighted centering
  
  preda=predb=x*0
  
  
  n_a = table(f1)
  a_hat_denom = as.vector(n_a+lambda_a)
  
  n_b = table(f2)
  b_hat_denom = as.vector(n_b+lambda_b)
  
  iter = 0
  relative_change = epsilon+2
  
  while((relative_change > epsilon)&(iter<=max_iter)){ 
    iter = iter+1
    deltaf=0
    
    preda_old = preda
    predb_old = predb
    
    resida = x - predb
    a_hat = level_means(resida,f1,a_hat_denom)
    preda = a_hat[f1,]

    
    residb = x - preda
    b_hat = level_means(residb,f2,b_hat_denom)
    predb = b_hat[f2,]
    
    
    
    deltaf=mean((preda_old-preda)^2)+mean((predb_old-predb)^2)
    dd = mean((preda_old+predb_old)^2)
    relative_change = deltaf/dd
    if(trace.it){cat("relative change in betahat and covariance betahat iterations is : ",relative_change,"\n")}
    
  }
  fitted_val = preda+predb
  
  
  l = list("Iteration" = iter, 
           "fitted_val" = fitted_val)  

  return(l)
  
}
backfitting_betahat_covariance <-function(y,x,f1,f2,epsilon,max_iter,random_ef=FALSE,trace.it=FALSE,var_lm_under_gls_ind=FALSE){
  #this function finds covariance of betahat_gls under gls using sandwich method. If user wants to obtain blup estimates
  #user neeeds to give random_ef = TRUE

  start_time = Sys.time()
  lambdas = get_lambdas(y,x,f1,f2) 
  lambda_a = lambdas[1]
  lambda_b = lambdas[2]
  sigma_sq_e = lambdas[3]
  out_cov_centering = backfitting_covariance_centering(y,x,f1,f2,lambda_a,lambda_b,epsilon,max_iter,trace.it)
  resid_sx = x - out_cov_centering$fitted_val#out_cov$fitted_val
  backfit_niter = out_cov_centering$Iteration
  xwx_inv = t(x)%*%resid_sx
  betahat = solve(xwx_inv)%*%t(resid_sx)%*%y
  

  cov_gls_through_backfitting = sigma_sq_e*solve(xwx_inv)
  #naive cov matrix with problem at (1,1)
  
  resid_sx_f1 = tibble(f1,as_tibble(resid_sx))
  ZA_t_resid_sx = resid_sx_f1 %>% group_by(f1) %>% summarize_all(sum,na.rm=TRUE)
  ZA_t_resid_sx = data.matrix(ZA_t_resid_sx[,-1])
  
  resid_sx_f2 = tibble(f2,as_tibble(resid_sx))
  ZB_t_resid_sx = resid_sx_f2 %>% group_by(f2) %>% summarize_all(sum,na.rm=TRUE)
  ZB_t_resid_sx = data.matrix(ZB_t_resid_sx[,-1])
  
  temp_cov = sigma_sq_e*t(resid_sx)%*%resid_sx + (sigma_sq_e/lambda_a) * t(ZA_t_resid_sx) %*% ZA_t_resid_sx +(sigma_sq_e/lambda_b) * t(ZB_t_resid_sx) %*% ZB_t_resid_sx 
  cov_gls_backfitting_sandwiched = solve(xwx_inv) %*% temp_cov %*% solve(xwx_inv)
  
  if(var_lm_under_gls_ind){V_lm_under_gls = var_betahat_lm_under_gls(x,f1,f2,lambdas)}
  end_time = Sys.time()
  
  time_to_compute = difftime(end_time,start_time,unit="secs")
  l = list("betahat"=betahat,"covariance_matrix"=cov_gls_backfitting_sandwiched,
           "Computational_Time" = time_to_compute, "Iteration" = backfit_niter, 
           "lambdas" = lambdas)  

  if(random_ef){
    rf_resu = backfitting_random_effect(y-x%*%betahat,x,f1,f2,lambda_a,lambda_b,epsilon=epsilon,
      max_iter=max_iter,trace.it=trace.it)
    l[["a_hat"]] = rf_resu$a_hat
    l[["b_hat"]] = rf_resu$b_hat
  }
  if(var_lm_under_gls_ind){
    l[["var_lm_under_gls"]] = V_lm_under_gls 
  }
  
return(l)
}

var_betahat_lm_under_gls <-function(x,f1,f2,lambdas){
  #this function finds variance of betahat ols under gls
  xt_x_inv = solve(t(x)%*%(x))
  
  
  df1 = tibble(f1,as_tibble(x))
  pwsum1 = df1 %>%group_by(f1) %>%summarize_all(sum)
  za_t_x=data.matrix(pwsum1[,-1])
  
  df2 = tibble(f2,as_tibble(x))
  pwsum2 = df2 %>%group_by(f2) %>%summarize_all(sum)
  zb_t_x=data.matrix(pwsum2[,-1])
  
  
  cov_z_x = (lambdas[3]/lambdas[1])*t(za_t_x)%*%za_t_x + (lambdas[3]/lambdas[2])*t(zb_t_x)%*%zb_t_x + lambdas[3]*t(x)%*%x
  
  cov_sandwiched = xt_x_inv%*%cov_z_x%*%xt_x_inv
  return(cov_sandwiched)
}

