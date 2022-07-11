simulation = function(population, R) {
  simulation_implementation = function(population, R, n) {
    ############################################
    ## Functions to help with the calculation
    ############################################
    cal_relative_bias = function(population_val, sample_val) {
      return((sample_val - population_val)/population_val)
    }
    is_in_CI = function(sample_mean, true_mean, sample_variance) {
      lower_limit_CI = sample_mean - 1.96 * sqrt(sample_variance)
      upper_limit_CI = sample_mean + 1.96 * sqrt(sample_variance)
      return(true_mean > lower_limit_CI & true_mean < upper_limit_CI)
    }
    #########################################
    # Variables related to the population
    #########################################
    # Population size
    N = length(population)
    # Population mean
    population_mean = mean(population)
    # Population variance of mean
    population_var = (1 - n/N) * var(population) / n
    ###############################################
    # Initialisation of variables to store simulation results
    #####################################################
    mean_point = numeric(R)
    relative_bias_point = numeric(R)
    variance_point = numeric(R)
    relative_bias_var = numeric(R)
    is_in_CI_indicator = numeric(R)
    
    # Main
    # Simulation: Simple Random Sampling Without Replacement (SRSWOR)
    for (r in 1:R) {
      # Generate indices for each sample (SRSWOR)
      index_of_samples = sample(c(1:N), size = n, replace = FALSE)
      
      # Calculate sample quantities
      sample = population[index_of_samples]
      sample_mean = mean(sample)
      sample_var = (1 - n/N) * var(sample) / n
      
      # Save sample variance for return
      mean_point[r] = sample_mean
      variance_point[r] = sample_var
      
      # Calculate and save simulated quantities for return
      relative_bias_point[r] = cal_relative_bias(population_val = population_mean,
                                                 sample_val = sample_mean)
      relative_bias_var[r] = cal_relative_bias(population_val = population_var,
                                               sample_val = sample_var)
      is_in_CI_indicator[r] = is_in_CI(sample_mean = sample_mean, 
                                       true_mean = population_mean, #population_mean, 
                                       sample_variance = sample_var)
    }
    ###########################
    # Return of the function
    ###########################
    rt_list = list(
      sample_size = n,
      population_mean = population_mean,
      mean_point_estimator = mean(mean_point),
      population_var = population_var,
      variance_point_estimator = mean(variance_point),
      relative_bias_point_estimator = mean(relative_bias_point) * 100,
      relative_bias_var_estimator = mean(relative_bias_var) * 100,
      coverage_prob = mean(is_in_CI_indicator)
    )
    return(rt_list)
  }
  
  output_df = function(output_1_n_10, output_1_n_25, 
                       output_1_n_50, output_1_n_100, output_1_n_500) {
    
    output_list = list(output_1_n_10, output_1_n_25,
                       output_1_n_50, output_1_n_100, 
                       output_1_n_500)
    Rel.Bias.Point.Est = c()
    Var.Point.Est = c()
    Rel.Bias.Var.Est = c()
    Coverage.Prob = c()
    for (index in 1:5){
      Rel.Bias.Point.Est = append(Rel.Bias.Point.Est, 
                                  output_list[[index]]$relative_bias_point_estimator)
      Var.Point.Est = append(Var.Point.Est, 
                             output_list[[index]]$variance_point_estimator)
      Rel.Bias.Var.Est = append(Rel.Bias.Var.Est, 
                                output_list[[index]]$relative_bias_var_estimator)
      Coverage.Prob = append(Coverage.Prob, 
                             output_list[[index]]$coverage_prob)
    }
    
    output.data <- data.frame(
      sample.size = c("n=10","n=25","n=50","n=100","n=500"),
      Rel.Bias.Point.Est = round(Rel.Bias.Point.Est, 3),
      Var.Point.Est      = round(Var.Point.Est, 3),
      Rel.Bias.Var.Est   = round(Rel.Bias.Var.Est, 3),
      Coverage.Prob      = round(Coverage.Prob, 3),
      stringsAsFactors = FALSE
    )
    return(output.data)
  }
  
  output_1_n_10 = simulation_implementation(population, R = R, n=10)
  output_1_n_25 = simulation_implementation(population, R = R, n=25)
  output_1_n_50 = simulation_implementation(population, R = R, n=50)
  output_1_n_100 = simulation_implementation(population, R = R, n=100)
  output_1_n_500 = simulation_implementation(population, R = R, n=500)
  
  output_dataframe = output_df(output_1_n_10, output_1_n_25, output_1_n_50, 
                               output_1_n_100,output_1_n_500)
  return(output_dataframe)
}
