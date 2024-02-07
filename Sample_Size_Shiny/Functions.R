#********************************************************************************************
#********************************************************************************************
#***** List of Gaussian correlation matrices, each one based on a univariate parameter ******
#********************************************************************************************
#********************************************************************************************
corr.matrix.list.gauss = function(theta, x){
  d = length(theta)
  D = as.matrix(dist(x))
  
  helper = function(i){
    exp(-D^2/theta[i])
  }
  
  return(lapply(1:d, helper))
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#***** List of Power Exponential  correlation matrices, each one based on a univariate ******
#***** parameter                                                                       ******
#********************************************************************************************
#********************************************************************************************
corr.matrix.list.pow.exp = function(theta, p, x){
  d = length(theta)
  D = as.matrix(dist(x))
  
  helper = function(i){
    exp(-D^p[i]/theta[i])
  }
  
  return(lapply(1:d, helper))
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#*****                       The Matern correlation function                           ******
#********************************************************************************************
#********************************************************************************************
Matern = function(nu, h, theta){
  ifelse(h == 0, 1, 
         (2*sqrt(nu)*h/theta)^nu*besselK(2*sqrt(nu)*h/theta, nu)/(gamma(nu)*2^(nu - 1)))
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#***** List of Gaussian correlation matrices, each one based on a univariate parameter ******
#********************************************************************************************
corr.matrix.list.matern = function(theta, nu, x){
  d = length(theta)
  D = as.matrix(dist(x))
  
  helper = function(i){
    Matern(nu, D, theta[i])
  }
  
  return(lapply(1:d, helper))
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#***** List of correlation matrices, each one based on a univariate parameter and the   *****
#***** correlation family chosen by the user                                            *****
#********************************************************************************************
#********************************************************************************************
corr.matrix.list = function(x, theta, p, nu, corr.family="gauss"){
  switch(corr.family,
         "Gaussian" = corr.matrix.list.gauss(theta, x),
         "Power Exponential" = corr.matrix.list.pow.exp(theta, p, x),
         "Matern" = corr.matrix.list.matern(theta, nu, x)
  )
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#*****                  Fredholm Integral Equation numerical solution                   *****
#********************************************************************************************
#********************************************************************************************
eigenvals = function(theta, p, nu, corr.family, n.grid=50){
  if((prod(!is.na(theta))==1) & 
     is.numeric(theta) & 
     min(theta)>0 & min(p)>=1 & 
     max(p)<=2 & 
     nu>0){
    d = length(theta)
    x = seq(.5/n.grid, 1-.5/n.grid, length=n.grid) #grid 
    R.list = corr.matrix.list(x, theta, p, nu, corr.family)
    
    decomp = eigen(R.list[[1]], symmetric=1)
    Lambda = decomp$values/n.grid
    if(d>1){
      for(i in 2:d){
        decomp = eigen(R.list[[i]], symmetric=1)
        temp.lambda = decomp$values/n.grid
        temp.lambda = temp.lambda[temp.lambda > 5e-10]
        Lambda = c(outer(Lambda, temp.lambda))
        Lambda = Lambda[Lambda > 5e-10]
      }
    }
    
    Lambda = Lambda[Lambda > 5e-10]
    Lambda = sort(Lambda, decreasing=T)
    
    return(Lambda)
  } else{
    return(0)
  }
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#*****                              Evaluating output arguments                        ******
#********************************************************************************************
#********************************************************************************************
unexplained_curve = function(theta, p, nu, corr.family, threshold, n.budget){
  d = length(theta)
  lambda = eigenvals(theta, p, nu, corr.family, 50)
  if(max(lambda)==0){
    return(list(lambda_curve=NA, unexplained=NA, 
                n.required=NA, delta=NA))
  } else{
    sum.lambda = cumsum(lambda)
    lambda_curve = sqrt(1-sum.lambda)
    unexplained = sqrt(1-sum.lambda[n.budget])
    temp = which(lambda_curve <= threshold)
    n.required = ifelse(length(temp)>0, min(temp), length(lambda))
    delta = max(((length(lambda)/20)%/%5)*5, 1)
    
    return(list(lambda_curve=lambda_curve, unexplained=unexplained, 
                n.required=n.required, delta=delta))
  }
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#*****                             Invalid input error message                         ******
#********************************************************************************************
#********************************************************************************************
error_message = function(goalInput, n.budget, threshold, theta, p, nu){
  if(sum(is.na(c(theta, n.budget, threshold, p, nu))>0)){
    message = "Make sure all inputs are numeric"
  } else if(min(theta)<=0){
    message = "All correlation parameters must be strictly positive!"
  } else if((goalInput) & (threshold <=0)){
    message = "The Root Average Unexplained Variation must be strictly positive!"
  } else if(!(goalInput) & (n.budget <= 0)){
    message = "Sample size must be strictly positive!"
  } else if(max(p)>2 | min(p)<1){
    message = "All powers must be between 1 and 2!"
  } else if(nu<0){
    message = "Smothness parameter must be strictly positive!"} else{
      message = ""
    }
  
  return(message)
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#*****                                     Output Graph                                ******
#********************************************************************************************
#********************************************************************************************
plot_graph = function(error, curve_object, goalInput, 
                      threshold, n.budget, p, nu){ 
  if(error==""){
    delta = curve_object$delta
    lambda_curve = curve_object$lambda_curve
    n.required = curve_object$n.required
    unexplained = curve_object$unexplained
    
    up.lim.y = ifelse(goalInput==1,
                      max(1-lambda_curve[1], threshold+.05),
                      max(1-lambda_curve[1], unexplained + .1))
    up.lim.x = ifelse(goalInput==1, length(lambda_curve),
                      min(length(lambda_curve), n.budget*5))
    
    if(goalInput==1){
      lab.locs =  c(-Inf, 1, 
                    seq(delta, 
                        length(lambda_curve) + delta, 
                        by=delta))
    } else{
      lab.locs = c(-Inf, n.budget, 
                   round(seq(1, 
                             up.lim.x, 
                             length=min(ceiling(up.lim.x/10), 100))))
    }
    
    df = data.frame(x = 1:length(lambda_curve), y = lambda_curve)
    
    if(goalInput==1){
      df2 = data.frame(x = c(1, n.required, n.required),
                       y = c(threshold, threshold, 0))
      
      x_point = n.required
      y_point = threshold
    } else{
      df2 = data.frame(x = c(1, n.budget, n.budget),
                       y = c(unexplained, unexplained, 0))
      
      x_point = n.budget
      y_point = unexplained
    }
    
    ggplot(df, aes(x, y)) + 
      geom_line(df2, mapping = aes(x, y), 
                col = 'navy blue', 
                linetype = 'dashed',
                linewidth = 1) + 
      geom_line(linewidth = 1.5, col = 'red') + 
      scale_x_continuous(expand = c(0, 0), 
                         lim = c(1, up.lim.x)) + 
      scale_y_continuous(expand = c(0, 0), 
                         lim = c(0, up.lim.y)) + 
      xlab('Sample Size') + 
      ylab('Root Average Unexplained Variation') + 
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_line(colour = "grey92"),
            axis.title.x = element_text(size=18, face="bold"),
            axis.title.y = element_text(size=18, face="bold"),
            axis.text.x = element_text(hjust=1, size=14),
            axis.text.y = element_text(hjust=1, size=14),
            plot.title = element_text(size=18, hjust = .5,
                                      face="bold.italic")) + 
      geom_point(aes(x = x_point, y= y_point), 
                 shape=21, size=5, fill="white", 
                 colour = "dodgerblue4", stroke = 2)
  }
  else{
    plot(c(0,1), c(0,1), xlab='', ylab='', axes=0, col=0)
    text(x=.5, y=1, error, cex=1.5, col=2)
  }
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#*****                                     Output Text                                 ******
#********************************************************************************************
#********************************************************************************************
text_to_print = function(error, curve_object, goalInput, n.budget, threshold){
  n.required = curve_object$n.required
  unexplained = curve_object$unexplained
  if(error==''){
    message = ifelse(goalInput==1, 
                     paste("A Root Average Unexplained Variation of", threshold, 
                           "requires a sample size of at least", n.required),
                     paste("The Root Average Unexplained Variation for a sample size of",
                           as.numeric(n.budget), "is at least", round(unexplained,3))
    )
  } else{
    message = ''
  }
  return(message)
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#*****                             Invalid input error message                         ******
#********************************************************************************************
#********************************************************************************************
error_message_path = function(theta, p, nu, n_path){
  if(sum(is.na(c(n_path, p, theta, nu)))>0){
    message = "Make sure all inputs are numeric"
  } else if(theta<=0){
    message = "Correlation length must be strictly positive!"
  } else if(n_path <= 0){
    message = "Number of realizations must be strictly positive!"
  } else if(p>2 | p<1){
    message = "Power must be between 1 and 2!"
  } else if(nu<0){
    message = "Smoothness parameter must be strictly positive!"} else{
      message = ""
    }
  
  return(message)
}
#********************************************************************************************
#********************************************************************************************





#********************************************************************************************
#********************************************************************************************
#*****                                     Output Graph                                ******
#********************************************************************************************
#********************************************************************************************
sample_path_graph = function(error, theta, p, nu, corr_family, n.paths){
  if(error==""){
    n.grid = 200
    x = seq(.5/n.grid, 1-.5/n.grid, length=n.grid) #grid 
    C = corr.matrix.list(x, theta, p, nu, corr_family)[[1]]
    y = matrix(MASS::mvrnorm(n.paths, rep(0, n.grid), C), ncol=n.paths, byrow=T)
    y = c(y)
    df = data.frame(path = rep(1:n.paths, each = n.grid), 
                    x = rep(x, n.paths),
                    z = y)
    df$path = as.factor(df$path)
    
    p = ggplot(df, aes(x, z, col = path)) + 
      geom_line(linewidth = 1, show.legend = FALSE) + 
      xlab('x') + 
      ylab('Z(x)') + 
      scale_x_continuous(expand = expansion(c(0,0))) + 
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_line(colour = "grey92"),
            axis.title.x = element_text(size=18, face="bold"),
            axis.title.y = element_text(size=18, face="bold"),
            axis.text.x = element_text(hjust=1, size=14),
            axis.text.y = element_text(hjust=1, size=14),
            plot.title = element_text(size=18, hjust = .5,
                                      face="bold.italic")) + 
      ylim(c(-3.5, 3.5)) 
    
    print(p)
  } else{
    plot(c(0,1), c(0,1), xlab='', ylab='', axes=0, col=0)
    text(x=.5, y=.5, error, cex=1.5, col=2)
  }
}
#********************************************************************************************
#********************************************************************************************



#********************************************************************************************
#********************************************************************************************
#****         Pulling the critical sample size out of the RAUV curve object             *****
#********************************************************************************************
#********************************************************************************************
min_samp_size = function(theta, threshold, corr_family, nu){
  d = length(theta)
  unexplained_curve(theta, rep(2,d), nu, corr_family, 
                    threshold, 10*d)$n.required
}
#********************************************************************************************
#********************************************************************************************



#********************************************************************************************
#********************************************************************************************
#****          a technical function used to print a progress bar on the screen           ****
#********************************************************************************************
#********************************************************************************************
apply_pb = function(X, MARGIN, FUN, ...){
  env = environment()
  pb_Total = sum(dim(X)[MARGIN])
  counter = 0
  pb = txtProgressBar(min = 0, max = pb_Total, style = 3)
  
  wrapper = function(...){
    curVal = get("counter", envir = env)
    assign("counter", curVal +1 ,envir= env)
    setTxtProgressBar(get("pb", envir= env), curVal +1)
    FUN(...)
  }
  res = apply(X, MARGIN, wrapper, ...)
  close(pb)
  res
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#*****                             Invalid input error message                         ******
#********************************************************************************************
#********************************************************************************************
error_message_robust = function(threshold, theta_low, theta_high, 
                                nu, quantile, MC_size){
  theta = c(theta_low, theta_high)
  rest = c(threshold, quantile, MC_size)
  
  if(sum(is.na(c(theta, rest))>0)){
    message = "Make sure all inputs are numeric"
  } else if(min(theta)<=0){
    message = "All correlation parameters must be strictly positive!"
  } else if(threshold <=0){
    message = "The Root Average Unexplained Variation must be strictly positive!"
  } else if(min(theta_high-theta_low)<0){
    message = "High correlation parameter values must be greater than low correlation parameter values"
  } else if(nu<0){
    message = "Smoothness parameter must be strictly positive!"
  } else if((quantile<0)|(quantile>1)){
    message = "Quantile must be between 0 and 1!"
  } else if(MC_size < 1000){
    message = "Monte Carlo sample size must be at least 1000!"
  } else{
    message = ""
  }
  
  return(message)
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#****         Calculating  the robust sample size and plotting the histogram            ****
#********************************************************************************************
#********************************************************************************************
samp_size_robust = function(lower, upper, threshold, 
                            MC_size, quant, corr_family, nu){
  unif.samp = function(lim){
    runif(MC_size)*diff(lim)+lim[1]
  }
  
  error = error_message_robust(threshold, lower, upper, nu, quant, MC_size)
  
  if(error==""){
    lims = cbind(lower, upper)
    samp = apply(lims, 1, unif.samp)	
    
    n.samp = apply_pb(samp, 1, min_samp_size, threshold=threshold, 
                      corr_family=corr_family, nu=nu)
    n.crit = ceiling(quantile(n.samp, quant))
    df = data.frame(SS = n.samp)
    breaks = pretty(range(n.samp), n = nclass.FD(n.samp), min.n = 1)
    bwidth = breaks[2] - breaks[1]
    
    ggplot(df, aes(x = SS)) + 
      geom_histogram(alpha = 0.7, 
                     position="identity", 
                     aes(y = ..density..), 
                     color="black", fill="royal blue",
                     binwidth = bwidth) + 
      geom_vline(xintercept = n.crit, color="black", 
                 linetype="dashed", linewidth=1) +
      xlab("Sample Size") + 
      ylab("Density") + 
      ggtitle(substitute(paste("Sample size distribution, ", n>=n0), 
                         list(n0=n.crit))) + 
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_line(colour = "grey92"),
            axis.title.x = element_text(size=18, face="bold"),
            axis.title.y = element_text(size=18, face="bold"),
            axis.text.x = element_text(hjust=1, size=14),
            axis.text.y = element_text(hjust=1, size=14),
            plot.title = element_text(size=18, hjust = .5,
                                      face="bold.italic")) + 
  scale_y_continuous(expand = c(0, 0))
  } else {
    plot(c(0,1), c(0,1), xlab='', ylab='', axes=0, col=0)
    text(x=.5, y=.5, error, cex=1.5, col=2)
  }
}
#********************************************************************************************
#********************************************************************************************