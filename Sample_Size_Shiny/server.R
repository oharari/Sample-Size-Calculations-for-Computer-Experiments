library(shiny)
library(MASS)
library(ggplot2)


source("Functions.R")


shinyServer(
  function(input, output){    
      
    output$action_thetas <- renderUI({
    	num_thetas <- as.integer(input$d)
    	theta.vec <- rep(1, num_thetas)
    	lapply(1:(num_thetas+1+(input$corr_family=="Matern")), function(i) {
      		if(i>(1+(input$corr_family=="Matern"))){
        		list(numericInput(paste0("theta", i-1-(input$corr_family=="Matern")), 
        		                  label=HTML(paste0("Correlation length &theta;",i-1-(input$corr_family=="Matern"),":")), 
                     value = 1, min=.1, step=.1)
                 )
      		} else if(i==1){
  
        	  	switch(input$goal,
              		   "var_to_n" = numericInput("var_target", label="Target RAUV", 
                                    value = 0.05, min = .005, max = .995, step = .005),
            			"n_to_var" = textInput("n_budget", label="Sampling budget", value=10*input$d)
        			   )
      		    }
        	  	else if(input$corr_family=="Matern")
        	  	{
        	  	  numericInput("nu", label=HTML(paste0("&nu; (smoothness parameter)")), value = 2.5, min=.5, max=100, step=.5)
        	  	}
            
    	}) #end of lapply    
      
  	})#end of renderUI
    
    
    output$action_p <- renderUI({
      num_p <- as.integer(input$d)
      p.vec <- rep(1, num_p)
      
      if(input$corr_family=="Power Exponential")
      {
      lapply(1:num_p, function(i){
        list(numericInput(paste0("p", i), label=paste0("Roughness parameter p", i, ":"), 
                          value = 1, min=1, max=2, step=.1)
        )
                }
             
        )
      }
    })
  
	goalInput <- eventReactive(input$go,{
    	goalInput <- ifelse(input$goal=="var_to_n", 1, 0)
  	})
  
	
 
	 threshold <- eventReactive(input$go,{
    	threshold <- ifelse(input$goal=="var_to_n", as.numeric(input$var_target), 0.05)
    })
  
  
  
  	n.budget <- eventReactive(input$go,{
    	n.budget <- ifelse(input$goal=="var_to_n", 10*input$d, as.integer(input$n_budget))
  	})
  
  	
  
  thetas <- eventReactive(input$go, {
    	d <- as.numeric(input$d)  
    	thetas <- sapply(1:d, function(i) {
      		as.numeric(input[[paste0("theta", i)]])
    	})
  	})
  	
  	
  	
  p_s <- eventReactive(input$go, {
  	  d <- as.numeric(input$d)  
  	  p_s <- sapply(1:d, function(i) {
  	    ifelse(input$corr_family=="Power Exponential", 
  	           as.numeric(input[[paste0("p", i)]]), 2)
  	  })
  	})
  	
  	
  	
  nu_par <- eventReactive(input$go, {
  	             ifelse(input$corr_family=="Matern", input$nu, 2.5)
  	                        })
  
  
  myerror <- eventReactive(input$go, {
		myerror <- error_message(goalInput(), n.budget(), threshold(), thetas(), p_s(), nu_par())
  	})
  	
  	
  my_corr_family <- eventReactive(input$go, {
  	  my_corr_family <- input$corr_family
  	})
    
  
  lambdas <- eventReactive(input$go, {
    	lambdas <- eigenvals(thetas(), p_s(), nu_par(), my_corr_family(), n.grid=50)
    })

  
  my_curve <- eventReactive(input$go, {    
    	my_curve <- unexplained_curve(thetas(), p_s(), nu_par(), my_corr_family(), threshold(), n.budget())
  	})
  
  
	my_text <- eventReactive(input$go, {
  		my_text <- text_to_print(myerror(), my_curve(), goalInput(), n.budget(), threshold())
	})
 
 
 	output$graph <- renderPlot({
   		plot_graph(myerror(), my_curve(), goalInput(), threshold(), n.budget(), p_s(), nu_par())
 	})
  
  
  	output$mytext <- renderText({
    	my_text()
  	})
  	
  	
  	
    output$p_or_nu <- renderUI({
      if(input$corr_family2=="Matern"){
        numericInput("samp_nu", HTML("Smothness parameter (&nu;)"), 
                     min=.5, value=2.5, step=.5)
      }
      else if(input$corr_family2=="Power Exponential"){
        numericInput("samp_p", "Roughness parameter (p)", 
                     min=1, value=1, max=2, step=.1)
      }
    })
    
    
    corr_family_path <- reactive({
      corr_family_path <- input$corr_family2
    })
    
    theta_path <- reactive({
      theta_path <- as.numeric(input$samp_theta)
    })
    
    p_path <- reactive({
      p_path <- ifelse(input$corr_family2=="Power Exponential", as.numeric(input$samp_p), 1)
    })
    
    nu_path <- reactive({
      nu_path <- ifelse(input$corr_family2=="Matern", as.numeric(input$samp_nu), 2.5)
    })
    
    n_path <- reactive({
      n_path <- as.integer(input$samp_n)
    })
    
    error_path <- reactive({
      error_path <- error_message_path(theta_path(), p_path(), nu_path(), n_path())
    })
    
    output$samp_path <- renderPlot({
      sample_path_graph(error_path(), theta_path(), p_path(), nu_path(), corr_family_path(), n_path())
    })
    
   output$action_nu <- renderUI({
      if(input$corr_family3=="Matern")
      {
        numericInput("nu2", label=HTML(paste0("&nu; (smoothness parameter)")), value = 2.5, min=.5, max=100, step=.5)
      }
   })
   
   output$action_thetas_low <- renderUI({
       num_thetas <- as.integer(input$d2)
       lapply(1:(num_thetas), function(i){
           list(numericInput(paste0("theta_low", i), 
                             label=HTML(paste0("Correlation length &theta;",i," (lower):")), 
                             value = .75, min=.1, step=.1))}
       )
       
   })
   
   output$action_thetas_high <- renderUI({
     num_thetas <- as.integer(input$d2)
     lapply(1:(num_thetas), function(i){
       list(numericInput(paste0("theta_high", i), 
                         label=HTML(paste0("Correlation length &theta;",i," (upper):")), 
                         value = 1.25, min=.1, step=.1))}
     )
     
   })
   
   
   thetas_low <- eventReactive(input$go2, {
     d <- as.numeric(input$d2)  
     thetas_low <- sapply(1:d, function(i) {
       as.numeric(input[[paste0("theta_low", i)]])
     })
   }) 
   
   
   thetas_high <- eventReactive(input$go2, {
     d <- as.numeric(input$d2)  
     thetas_high <- sapply(1:d, function(i) {
       as.numeric(input[[paste0("theta_high", i)]])
     })
   }) 
    
   
   nu_par2 <- eventReactive(input$go2, {
     ifelse(input$corr_family3=="Matern", input$nu2, 2.5)
   })
    
   
   my_corr_family2 <- eventReactive(input$go2, {
     my_corr_family2 <- input$corr_family3
   })
   
   
   my_MC_size <- eventReactive(input$go2, {
     my_MC_size <- input$MC_size
   })
   
   
   my_quantile <- eventReactive(input$go2,{
     my_quantile <- input$quantile
   })
   
   
   threshold2 <- eventReactive(input$go2,{
     threshold2 <- input$var_target2
   })
   
   
   output$histogram <- renderPlot({
     withProgress(message = 'Calculation in progress', value=NULL,
       samp_size_robust(thetas_low(), thetas_high(), threshold2(), my_MC_size(), my_quantile(), my_corr_family2(), nu_par2()))
   })
   
  }
)