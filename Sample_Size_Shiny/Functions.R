#********************************************************************************************
#********************************************************************************************
#***** List of Gaussian correlation matrices, each one based on a univariate parameter ******
#********************************************************************************************
#********************************************************************************************
corr.matrix.list.gauss <- function(theta, x)
{
	d <- length(theta)
	D <- as.matrix(dist(x))
	R.list <- list()

	for(i in 1:d)
	{
		R.list[[i]] <- exp(-D^2/theta[i])
	}
	
	return(R.list)
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#***** List of Power Exponential  correlation matrices, each one based on a univariate ******
#***** parameter                                                                       ******
#********************************************************************************************
#********************************************************************************************
corr.matrix.list.pow.exp <- function(theta, p, x)
{
	d <- length(theta)
	D <- as.matrix(dist(x))
	R.list <- list()

	for(i in 1:d)
	{
		R.list[[i]] <- exp(-D^p[i]/theta[i])
	}
	
	return(R.list)
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#*****                       The Matern correlation function                           ******
#********************************************************************************************
#********************************************************************************************
Matern <- function(nu, h, theta)
{
  ifelse(h == 0, 1, (2*sqrt(nu)*h/theta)^nu*besselK(2*sqrt(nu)*h/theta, nu)/(gamma(nu)*2^(nu - 1)))
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#***** List of Gaussian correlation matrices, each one based on a univariate parameter ******
#********************************************************************************************
corr.matrix.list.matern <- function(theta, nu, x)
{
	d <- length(theta)
	D <- as.matrix(dist(x))
	R.list <- list()

	for(i in 1:d)
	{
		R.list[[i]] <- Matern(nu, D, theta[i])
	}
	
	return(R.list)
}
#********************************************************************************************
#********************************************************************************************




#********************************************************************************************
#********************************************************************************************
#***** List of correlation matrices, each one based on a univariate parameter and the   *****
#***** correlation family chosen by the user                                            *****
#********************************************************************************************
#********************************************************************************************
corr.matrix.list <- function(x, theta, p, nu, corr.family="gauss")
{
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
eigenvals <- function(theta, p, nu, corr.family, n.grid=50)
{
  if((prod(!is.na(theta))==1) & is.numeric(theta) & min(theta)>0 & min(p)>=1 & max(p)<=2 & nu>0)
  {
	  d <- length(theta)
	  x <- seq(.5/n.grid, 1-.5/n.grid, length=n.grid) #grid 
	  R.list <- corr.matrix.list(x, theta, p, nu, corr.family)

	  decomp <- eigen(R.list[[1]], symmetric=1)
	  Lambda <- decomp$values/n.grid
	  if(d>1)
	  {
		  for(i in 2:d)
		  {
			  decomp <- eigen(R.list[[i]], symmetric=1)
			  temp.lambda <- decomp$values/n.grid
			  temp.lambda <- temp.lambda[temp.lambda > 5e-10]
			  Lambda <- c(outer(Lambda, temp.lambda))
			  Lambda <- Lambda[Lambda > 5e-10]
		  }
	  }
	
	  Lambda <- Lambda[Lambda > 5e-10]
	  Lambda <- sort(Lambda, decreasing=T)

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
unexplained_curve <- function(theta, p, nu, corr.family, threshold, n.budget)
{
  d <- length(theta)
  lambda <- eigenvals(theta, p, nu, corr.family, 50)
  if(max(lambda)==0){
    return(list(lambda_curve=NA, unexplained=NA, 
                n.required=NA, delta=NA))
  } else{
    sum.lambda <- cumsum(lambda)
    lambda_curve <- sqrt(1-sum.lambda)
    unexplained <- sqrt(1-sum.lambda[n.budget])
    temp <- which(lambda_curve <= threshold)
    n.required <- ifelse(length(temp)>0, min(temp), length(lambda))
     delta <- max(((length(lambda)/20)%/%5)*5, 1)
  
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
error_message <- function(goalInput, n.budget, threshold, theta, p, nu)
{
    if(sum(is.na(c(theta, n.budget, threshold, p, nu))>0)){
    message <- "Make sure all inputs are numeric"
    } else if(min(theta)<=0){
    message <- "All correlation parameters must be strictly positive!"
    } else if((goalInput) & (threshold <=0)){
    message <- "The Root Average Unexplained Variation must be strictly positive!"
    } else if(!(goalInput) & (n.budget <= 0)){
    message <- "Sample size must be strictly positive!"
    } else if(max(p)>2 | min(p)<1){
    message <- "All powers must be between 1 and 2!"
    } else if(nu<0){
    message <- "Smothness parameter must be strictly positive!"} else{
    message <- ""
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
plot_graph <- function(error, curve_object, goalInput, threshold, n.budget, p, nu)
{ 
  if(error==""){
    delta <- curve_object$delta
    lambda_curve <- curve_object$lambda_curve
    n.required <- curve_object$n.required
    unexplained <- curve_object$unexplained
    
    up.lim.y <- ifelse(goalInput==1,
                     max(1-lambda_curve[1], threshold+.05),
                     max(1-lambda_curve[1], unexplained + .1))
    up.lim.x <- ifelse(goalInput==1, length(lambda_curve),
                       min(length(lambda_curve), n.budget*5))
    #up.lim.y <- max(1-lambda_curve[1], threshold+.05)
    #up.lim.x <- length(lambda_curve)
    if(goalInput==1){
      lab.locs <-  c(-Inf, 1, seq(delta, length(lambda_curve)+delta, by=delta))
    } else{
      lab.locs <- c(-Inf, n.budget, round(seq(1, up.lim.x, length=min(ceiling(up.lim.x/10), 100))))
    }
  
    par(mar=c(5,5,1,1))
    plot(lambda_curve, type='l', lwd=2, col=4,
         ylab='Root Average Unexplained Variation', xlab='Sample Size',
         cex.lab=1.4, axes=0, ylim=c(0, up.lim.y), xlim = c(1,up.lim.x))
    axis(1, pos=0, at = lab.locs)
    axis(2, pos=1, at=seq(0,1,by=.05))
  
    col.n <- ifelse(goalInput==0, 2, 0)
    col.var <- ifelse(goalInput==0, 0, 3)
  
    lines(c(-1,n.required), c(threshold,threshold), lty=2, lwd=2, col=col.var)
    lines(c(n.required,n.required), c(0,threshold), lty=2, lwd=2, col=col.var)
  
    lines(c(n.budget, n.budget), c(0,unexplained), lty=2, lwd=2, col=col.n)
    lines(c(0, n.budget), c(unexplained, unexplained), lty=2, lwd=2, col=col.n)
    lines(lambda_curve, lwd=2, col=4)
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
text_to_print <- function(error, curve_object, goalInput, n.budget, threshold)
{
  n.required <- curve_object$n.required
  unexplained <- curve_object$unexplained
  if(error==''){
    message <- ifelse(goalInput==1, 
           paste("A Root Average Unexplained Variation of", threshold, 
                  "requires a sample size of at least", n.required),
           paste("The Root Average Unexplained Variation for a sample size of",
                 as.numeric(n.budget), "is at least", round(unexplained,3))
         )
    } else{
      message <- ''
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
error_message_path <- function(theta, p, nu, n_path)
{
    if(sum(is.na(c(n_path, p, theta, nu)))>0){
    message <- "Make sure all inputs are numeric"
    } else if(theta<=0){
    message <- "Correlation length must be strictly positive!"
    } else if(n_path <= 0){
    message <- "Number of realizations must be strictly positive!"
    } else if(p>2 | p<1){
    message <- "Power must be between 1 and 2!"
    } else if(nu<0){
    message <- "Smoothness parameter must be strictly positive!"} else{
    message <- ""
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
sample_path_graph <- function(error, theta, p, nu, corr_family, n.paths)
{
	if(error=="")
	{
		n.grid <- 200
		x <- seq(.5/n.grid, 1-.5/n.grid, length=n.grid) #grid 
		C <- corr.matrix.list(x, theta, p, nu, corr_family)[[1]]
		y <- matrix(mvrnorm(n.paths, rep(0, n.grid), C), ncol=n.paths, byrow=T)
		lims <- c(min(y), max(y))
		lims <- c(min(lims[1],-3), max(lims[2],3))
	
		par(mar=c(5,5,4,1))
		plot(x, y[,1], type='l', lwd=2, xlab='x', ylab='Z(x)', cex.lab=1.5, 
	    	 xlim=c(0,1), ylim=lims, axes=0)
	   	axis(1, pos=lims[1]-.05)
	   	axis(2, pos=0, at = seq(round(lims[1])-1, round(lims[2])+1, by=1))
		for(i in 1:n.paths)
		{
			lines(x, y[,i], col=i, lwd=2)
		}
	}

	else
	{
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
min_samp_size <- function(theta, threshold, corr_family, nu)
{
  d <- length(theta)
  unexplained_curve(theta, rep(2,d), nu, corr_family, threshold, 10*d)$n.required
}
#********************************************************************************************
#********************************************************************************************



#********************************************************************************************
#********************************************************************************************
#****          a technical function used to print a progress bar on the screen           ****
#********************************************************************************************
#********************************************************************************************
apply_pb <- function(X, MARGIN, FUN, ...)
{
  env <- environment()
  pb_Total <- sum(dim(X)[MARGIN])
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
  
  wrapper <- function(...)
  {
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir= env)
    setTxtProgressBar(get("pb", envir= env), curVal +1)
    FUN(...)
  }
  res <- apply(X, MARGIN, wrapper, ...)
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
error_message_robust <- function(threshold, theta_low, theta_high, nu, quantile, MC_size)
{
  theta <- c(theta_low, theta_high)
  rest <- c(threshold, quantile, MC_size)
  
  if(sum(is.na(c(theta, rest))>0)){
    message <- "Make sure all inputs are numeric"
  } else if(min(theta)<=0){
    message <- "All correlation parameters must be strictly positive!"
  } else if(threshold <=0){
    message <- "The Root Average Unexplained Variation must be strictly positive!"
  } else if(min(theta_high-theta_low)<0){
    message <- "High correlation parameter values must be greater than low correlation parameter values"
  } else if(nu<0){
    message <- "Smoothness parameter must be strictly positive!"
  } else if((quantile<0)|(quantile>1)){
    message <- "Quantile must be between 0 and 1!"
  } else if(MC_size < 1000){
    message <- "Monte Carlo sample size must be at least 1000!"
  } else{
    message <- ""
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
samp_size_robust <- function(lower, upper, threshold, MC_size, quant, corr_family, nu)
{
  unif.samp <- function(lim)
  {
    runif(MC_size)*diff(lim)+lim[1]
  }
  
  error <- error_message_robust(threshold, lower, upper, nu, quant, MC_size)
  
  if(error=="")
  {
     lims <- cbind(lower, upper)
     samp <- apply(lims, 1, unif.samp)	
  
     n.samp <- apply_pb(samp, 1, min_samp_size, threshold=threshold, 
                     corr_family=corr_family, nu=nu)
     n.crit <- ceiling(quantile(n.samp, quant))
  
     par(mar=c(5,5,3,2))
     hist(n.samp, main = " ", col="light gray", xlab="n", cex.lab=1.6)
     abline(v=n.crit, lwd=4)
     title(main=substitute(paste("Sample size distribution, ", n>=n0), 
                           list(n0=n.crit)), cex.main=1.6)
  }else
  {
    plot(c(0,1), c(0,1), xlab='', ylab='', axes=0, col=0)
    text(x=.5, y=.5, error, cex=1.5, col=2)
  }
}
#********************************************************************************************
#********************************************************************************************