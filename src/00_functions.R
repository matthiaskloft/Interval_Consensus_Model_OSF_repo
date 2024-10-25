################################################################################
# Helper Functions #
################################################################################

packages <- c("tidyverse", "extraDistr", "showtext")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(packages, update = F, character.only = T)

## Loading Google fonts (https://fonts.google.com/)
try(font_add_google("Inter", "google_font"))

#------------------------------------------------------------------------------>
# Transformation Functions ###
#------------------------------------------------------------------------------>

### Function to transform Simplex to bivariate normal --------------------------
simplex_to_bvn <-
  function(simplex, 
           type = "ilr") # c("ilr", "sb")
  {
    if (type == "ilr") {
      
    # ILR --------------------------------------------------------------------->  
      
      if (!is.data.frame(simplex) &
          !is.matrix(simplex)) {
        
        # vector
        
        Y <- rep(NA, 2)
        Y[1] = sqrt(1 / 2) * log(simplex[1] / simplex[3])
        Y[2] = sqrt(2 / 3) * log(simplex[2] / sqrt(simplex[1] * simplex[3]))
        
        return(Y)
        
        
      } else{
        
        # dataframe
        
        Y <- apply(
          X = simplex,
          MARGIN = 1,
          FUN = function(X) {
            Y <- rep(NA, 2)
            Y[1] = sqrt(1 / 2) * log(X[1] / X[3])
            Y[2] = sqrt(2 / 3) * log(X[2] / sqrt(X[1] * X[3]))
        
            return(Y)
          },
          simplify = FALSE
        )
        Y <- do.call(what = "rbind", args = Y)
        return(data.frame(x_bvn_1 = Y[, 1], x_bvn_2 = Y[, 2]))
      }
      
    } else{

    # SB ---------------------------------------------------------------------->
      
      if (!is.data.frame(simplex) &
          !is.matrix(simplex)) {
        
        # vector
        
        Y <- rep(NA, 2)
        Y[1] = log(simplex[1] / simplex[3])
        Y[2] = log(simplex[2] / (simplex[1] + simplex[3]))
        
        return(Y)
        
      } else{
        
        # dataframe
        
        Y <- apply(
          X = simplex,
          MARGIN = 1,
          FUN = function(X) {
            Y <- rep(NA, 2)
            Y[1] = log(X[1] / X[3])
            Y[2] = log(X[2] / (X[1] + X[3]))
            
            return(Y)
          },
          simplify = FALSE
        )
        Y <- do.call(what = "rbind", args = Y)
        return(data.frame(x_bvn_1 = Y[, 1], x_bvn_2 = Y[, 2]))
      }
    }
  }


### Convert bivariate normal to simplex --------------------------------------->
bvn_to_simplex <- function(bvn, type = "ilr") {
  
  if(type == "ilr"){
    ### ILR ------------------------------------------------------------------->
    
    if (!is.data.frame(bvn) &
        !is.matrix(bvn)) {
      
      # vector

      Y <- rep(NA, 3)
      Y[1] <- exp(sqrt(2) * bvn[1])
      Y[2] <- exp(sqrt(3/2) *  bvn[2] + bvn[1] / sqrt(2))
      Y[3] <- 1
      Y <- Y/sum(Y)
      
      names(Y) <- c("x_1", "x_2", "x_3")
      
      return(Y)
      
    } else {
      
      # dataframe
      
      Y <- apply(
        X = bvn,
        MARGIN = 1,
        FUN = function(X) {
          Y <- rep(NA, 3)
          Y[1] <- exp(sqrt(2) * X[1])
          Y[2] <- exp((sqrt(3 / 2) * X[2]) + (X[1] / sqrt(2)))
          Y[3] <- 1
          Y <- Y / sum(Y)
          
          names(Y) <- c("x_1", "x_2", "x_3")
          
          return(Y)
        },
        simplify = FALSE
      )
      Y <- do.call(what = "rbind", args = Y)
      
      return(data.frame(
        x_1 = Y[, 1],
        x_2 = Y[, 2],
        x_3 = Y[, 3]
      ))
    }
    
  } else{
    
    ### SB -------------------------------------------------------------------->
    
    if (!is.data.frame(bvn) &
        !is.matrix(bvn)) {
    
      # vector
      
      Y <- rep(NA, 3)
      Y[1] <- exp(bvn[1]) / ((exp(bvn[1]) + 1) * (exp(bvn[2]) + 1))
      Y[2] <- exp(bvn[2]) /  (exp(bvn[2]) + 1)
      Y[3] <- 1 /           ((exp(bvn[1]) + 1) * (exp(bvn[2]) + 1))
      names(Y) <- c("x_1", "x_2", "x_3")
      
      return(Y)
      
    } else {
      
      # dataframe
      
      Y <- apply(
        X = bvn,
        MARGIN = 1,
        FUN = function(X) {
          Y <- rep(NA, 3)
          Y[1] <- exp(X[1]) / ((exp(X[1]) + 1) * (exp(X[2]) + 1))
          Y[2] <- exp(X[2]) /  (exp(X[2]) + 1)
          Y[3] <- 1 /         ((exp(X[1]) + 1) * (exp(X[2]) + 1))
          
          names(Y) <- c("x_1", "x_2", "x_3")
          
          return(Y)
        },
        simplify = FALSE
      )
      Y <- do.call(what = "rbind", args = Y)
      
      return(data.frame(
        x_1 = Y[, 1],
        x_2 = Y[, 2],
        x_3 = Y[, 3]
      ))
    } # end if else vector or dataframe
  } # end if else ilr or sb
} # end function


#------------------------------------------------------------------------------>
# Simulation functions
#------------------------------------------------------------------------------>

# Generate interval data for one respondent (for illustration plots)
generate_itm_one_respondent <- function(
                                        Tr_loc,
                                        Tr_wid,
                                        lambda_loc,
                                        lambda_wid,
                                        E_loc,
                                        E_wid,
                                        a_loc,
                                        b_loc,
                                        b_wid,
                                        omega = 0) {
  idx <- 1:n_items
  
  Tr_splx <- cbind(Tr_loc, Tr_wid) %>% bvn_to_simplex()
  Tr_L <- Tr_splx[, 1]
  Tr_U <- 1 - Tr_splx[, 3]

  ### Model
  # precision
  sigma_loc <- exp(- log(E_loc) - log(lambda_loc))
  sigma_wid <- exp(-log(E_wid) - log(lambda_wid))
  
  # generate unbounded interval response data
  error <- extraDistr::rbvnorm(
    n = n_items,
    mean1 =  0,
    mean2 = 0,
    sd1 = sigma_loc,
    sd2 = sigma_wid,
    cor = omega
  )
  error_loc <- error[,1]
  error_wid <- error[,2]
  
  A_loc <- Tr_loc + error_loc
  A_wid <- Tr_wid + error_wid
  
  Y_loc <- A_loc * a_loc + b_loc
  Y_wid <- A_wid + b_wid
  
  # convert to bounded interval response data
  X <- bvn_to_simplex(cbind(Y_loc,Y_wid))
  
  # data frame of responses
  responses <- data.frame(
    idx = idx,
    Tr_loc = Tr_loc,
    Tr_wid = Tr_wid,
    error_loc = error_loc,
    error_wid = error_wid,
    A_loc = A_loc,
    A_wid = A_wid,
    Y_loc,
    Y_wid,
    x_splx_1 = X[, 1],
    x_splx_2 = X[, 2],
    x_splx_3 = X[, 3],
    sum = rowSums(X),
    x_L = X[, 1],
    x_U = 1 - X[, 3],
    Tr_L = Tr_L,
    Tr_U = Tr_U
  )
  return(responses)
}


# Generate interval data for Simulation Study --------------------------------->

generate_itm_data_sim_study <- 
  function(
    n_respondents,
    n_items,
    link = "ilr", # c("ilr", "sb")
    mu_Tr_loc = NULL,
    mu_Tr_wid = NULL,
    mu_E_loc = NULL,
    mu_E_wid = NULL,
    sigma_Tr_loc = NULL,
    sigma_Tr_wid = NULL,
    sigma_lambda_E_loc = NULL,
    sigma_lambda_E_wid = NULL,
    sigma_a_loc = NULL,
    sigma_b_loc = NULL,
    sigma_b_wid = NULL,
    omega_beta = NULL
    ){

    ### Hyper Priors ###
    
    # compute a benchmark for the mean and SD of the parameters
    if (link == "ilr") {
      
      mean_benchmark <- simplex_to_bvn(c(.4, .2, .4), type = link)
    } else{
      mean_benchmark <- simplex_to_bvn(c(.425, .15, .425), type = link)
    }
    sd_benchmark_loc <- simplex_to_bvn(c(.98, .01, .01), type = link)
    sd_benchmark_wid <- simplex_to_bvn(c(.495, .01, .495), type = link)
    
    # mean for Tr_loc
    mu_Tr_loc <- ifelse(is.null(mu_Tr_loc), mean_benchmark[1], mu_Tr_loc)
    # mean for Tr_wid
    mu_Tr_wid <- ifelse(is.null(mu_Tr_wid), mean_benchmark[2], mu_Tr_wid)
    # SD forTr_loc
    sigma_Tr_loc <- ifelse(is.null(sigma_Tr_loc), 
                           sd_benchmark_loc[1] / 4, 
                           sigma_Tr_loc)
    # SD Tr_wid
    sigma_Tr_wid <- ifelse(is.null(sigma_Tr_wid),
                           abs(sd_benchmark_wid[2] - mean_benchmark[2]) / 4, 
                           sigma_Tr_wid)
    
    # mean fro E_loc
    mu_E_loc <- ifelse(is.null(mu_E_loc),log(sigma_Tr_loc) , mu_E_loc)
    # mean for E_wid
    mu_E_wid <- ifelse(is.null(mu_E_wid), log(sigma_Tr_wid), mu_E_wid)
    
    # SDs for other parameters
    sigma_lambda_E_loc <- ifelse(is.null(sigma_lambda_E_loc), 
                                 .3, 
                                 sigma_lambda_E_loc)
    sigma_lambda_E_wid <- ifelse(is.null(sigma_lambda_E_wid),
                                 .3, 
                                 sigma_lambda_E_wid)
    sigma_a_loc <- ifelse(is.null(sigma_a_loc),
                          .3, 
                          sigma_a_loc)
    sigma_b_loc <- ifelse(is.null(sigma_b_loc),
                          sigma_Tr_loc / 3, 
                          sigma_b_loc)
    sigma_b_wid <- ifelse(is.null(sigma_b_wid), 
                          sigma_Tr_wid / 3, 
                          sigma_b_wid)
    if (is.null(omega_beta)) {
      omega_beta <- 3
    } else{
      omega_beta <- omega_beta
    }
    
    ### Indices ###
    n <- n_respondents * n_items
    ii <- rep(1:n_respondents, each = n_items)
    jj <- rep(1:n_items, times = n_respondents)
    
    ### True Parameters ###
    
    # true intervals
    Tr_loc <- rnorm(n_items, mu_Tr_loc, sigma_Tr_loc)
    Tr_wid <- rnorm(n_items, mu_Tr_wid, sigma_Tr_wid)
    if (link == "ilr") {
      Tr_splx <- cbind(Tr_loc, Tr_wid) %>% bvn_to_simplex(type = "ilr")
    } else {
      Tr_splx <- cbind(Tr_loc, Tr_wid) %>% bvn_to_simplex(type = "sb")
    }
    Tr_L <- Tr_splx[, 1]
    Tr_U <- 1 - Tr_splx[, 3]
    
    # discernibility
    lambda_loc <- 1 / exp(rnorm(n_items, 0, sigma_lambda_E_loc))
    lambda_wid <- 1 / exp(rnorm(n_items, 0, sigma_lambda_E_wid))
    # respondent proficiency
    E_loc <- 1 / exp(rnorm(n_respondents, mu_E_loc, sigma_lambda_E_loc))
    E_wid <- 1 / exp(rnorm(n_respondents, mu_E_wid, sigma_lambda_E_wid))
    # respondent scaling bias
    a_loc <- exp(rnorm(n_respondents, 0, sigma_a_loc))
    # respondent shifting bias
    b_loc <- rnorm(n_respondents, 0, sigma_b_loc)
    b_wid <- rnorm(n_respondents, 0, sigma_b_wid)
    
    ### Model ###
    
    # expected interval
    mu_loc <- Tr_loc[jj] * a_loc[ii] + b_loc[ii]
    mu_wid <- Tr_wid[jj]             + b_wid[ii]
    # precision
    sigma_loc <- a_loc[ii] / E_loc[ii] / lambda_loc[jj]
    sigma_wid <- 1 / E_wid[ii] / lambda_wid[jj]
    # residual correlation
    omega <- rbeta(n = n_items,
                   shape1 = omega_beta, 
                   shape2 = omega_beta) * 2 - 1

    # generate unbounded interval response data
    Y <- extraDistr::rbvnorm(n = n,
                             mean1 =  mu_loc, 
                             mean2 = mu_wid,
                             sd1 = sigma_loc,
                             sd2 = sigma_wid, 
                             cor = omega[jj])
    Y_loc <- Y[,1]
    Y_wid <- Y[,2]
    
    # convert to bounded interval response data
    if (link == "ilr") {
      X <- bvn_to_simplex(Y, type = "ilr")
    } else {
      X <- bvn_to_simplex(Y, type = "sb")
    }
    
    ### Save Objects ###
    
    # data frame of responses
    responses <- data.frame(
      ii = ii,
      jj = jj,
      Y_loc,
      Y_wid,
      x_splx_1 = X[, 1],
      x_splx_2 = X[, 2],
      x_splx_3 = X[, 3],
      sum = rowSums(X),
      x_L = X[, 1],
      x_U = 1 - X[, 3]
    )
    
    # list of true parameters
    parameters <- list(
      Tr_loc = Tr_loc,
      Tr_wid = Tr_wid,
      Tr_splx = Tr_splx,
      Tr_L = Tr_L,
      Tr_U = Tr_U,
      lambda_loc = lambda_loc,
      lambda_wid = lambda_wid,
      E_loc = E_loc,
      E_wid = E_wid,
      a_loc = a_loc,
      b_loc = b_loc,
      b_wid = b_wid,
      omega_beta = omega_beta,
      omega = omega,
      mu_Tr_loc = mu_Tr_loc,
      mu_Tr_wid = mu_Tr_wid,
      sigma_Tr_loc = sigma_Tr_loc,
      sigma_Tr_wid = sigma_Tr_wid,
      sigma_lambda_E_loc = sigma_lambda_E_loc,
      sigma_lambda_E_wid = sigma_lambda_E_wid,
      sigma_a_loc = sigma_a_loc,
      sigma_b_loc = sigma_b_loc,
      sigma_b_wid = sigma_b_wid,
      link = link
    )
    sim_data <- list(responses = responses, parameters = parameters)
    return(sim_data)
  }

#------------------------------------------------------------------------------>
# Plotting Functions
#------------------------------------------------------------------------------>

# helper: gather interval values for cumulative density plot ------------------>
gather_values <- function(lower,upper,item_id, step_size) {
  
  df_out <- map_dfr(1:length(lower), function(.x) {
    # gather values between bounds by step size
    samples <- seq(lower[.x], upper[.x], by = step_size) %>% as.double()
    item_id <- rep(item_id[.x], length(samples))
    
    return(data.frame(
      samples = samples, 
      item_id = item_id))
  })
  
  return(df_out)
}


#------------------------------------------------------------------------------>

# Plot Intervals on the unbounded 2D scale

plot_example_2Dscatter <- function(data) {
  plot <-
    data %>%
    ggplot() +
    geom_point(aes(x = Y_loc, y = Y_wid), size = 1.7, alpha = .2, shape = 16) +
    geom_vline(xintercept = 0, linetype = 2, alpha = .3) +
    geom_hline(yintercept = 0, linetype = 2, alpha = .3) +
    scale_x_continuous(limits = c(-5, 5), breaks = seq(-3, 3, 3)) +
    scale_y_continuous(limits = c(-6, 4), breaks = seq(-3, 3, 3)) +
    labs(x = "Unbounded Location", y = "Unbounded Width") +
    theme_itm() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "#6d6d6e", size = .3),
          axis.ticks = element_line(colour = "#6d6d6e", size = .3),
          plot.margin = margin(.1, .3, .1, .1, "cm")
          )
  
  return(plot)
}


#------------------------------------------------------------------------------>

# Plot Intervals on the bounded response scale

plot_example_intervals <- function(data) {
  plot <-
    data %>%
    ggplot() +
    geom_errorbar(
      aes(
        y = (idx),
        xmin = Tr_L,
        xmax = Tr_U
      ),
      width = 0,
      linewidth = 4,
      color = "grey80"
    ) +
    geom_errorbar(aes(
      y = (idx),
      xmin = x_L,
      xmax = x_U
    ),
    width = .6,
    linewidth = .7) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, .25),
      labels = c("0", ".25", ".50", ".75", "1"),
      expand = expansion(0, 0)
    ) +
    scale_y_continuous(expand = expansion(0, .2)) +
    labs(x = "Bounded Response Scale", y = "Item Index") +
    theme_itm() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = "#6d6d6e", size = .3),
      axis.ticks.x = element_line(colour = "#6d6d6e", size = .3),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
      plot.margin = margin(.1, .3, .1, .1, "cm")
    )
  return(plot)
}


#------------------------------------------------------------------------------>

# Plot Cumulative Density of Intervals & Consensus Interval

plot_intvls_aggregated <- function(lower,
                                   upper,
                                   item_id,
                                   item_name = NA,
                                   lower_mean,
                                   upper_mean,
                                   lower_mean_logit,
                                   upper_mean_logit,
                                   truth = NA,
                                   multiple = FALSE,
                                   min,
                                   max,
                                   facet_wrap = FALSE,
                                   design = NULL,
                                   step_size,
                                   show_quantiles = TRUE,
                                   ncol = 2) {
  # gather values between bounds
  gathered <-
    gather_values(lower, upper, item_id, step_size = step_size) %>%
    group_by(item_id) %>%
    # compute the maximum density
    mutate(
      max_density = max(table(round(
        samples, -floor(log10(step_size))
      )), na.rm = TRUE) %>%
        as.double(),
      median = median(samples, na.rm = TRUE),
      q_05 = quantile(samples, probs = .05),
      q_95 = quantile(samples, probs = .95)
    ) %>%
    ungroup() %>%
    full_join(
      data.frame(
        lower_mean = lower_mean,
        upper_mean = upper_mean,
        lower_mean_logit = lower_mean_logit,
        upper_mean_logit = upper_mean_logit,
        truth = as.numeric(truth),
        item_id = item_id,
        item_name = item_name
      ) %>%
        distinct()
    ) %>%
    mutate(item_id = factor(item_id),
           item_name = factor(item_name))
  
  plot <-
    try(    ggplot_aggregated_single(
      data = gathered,
      binwidth = step_size,
      min = min,
      max = max,
      facet_wrap = facet_wrap,
      design = design,
      show_quantiles = show_quantiles,
      ncol = ncol
    ))
  return(plot)
}

# helper: gather interval values -----------------------------------------------
gather_values <- function(lower,upper,item_id, step_size) {
  # identify columns for lower and upper response value
  
  
  df_out <- map_dfr(1:length(lower), function(.x) {
    # gather values between bounds by step size
    samples <- seq(lower[.x], upper[.x], by = step_size) %>% as.double()
    item_id <- rep(item_id[.x], length(samples))
    
    return(data.frame(samples = samples, item_id = item_id))
  })
  
  return(df_out %>% arrange(item_id))
}

# helper: ggplot for aggregated responses --------------------------------------
ggplot_aggregated_single <-
  function(data,
           min = min,
           max = max,
           binwidth,
           facet_wrap = FALSE,
           design = NULL,
           show_quantiles = TRUE,
           ncol = 2) {
    # avoid warnings
    #samples <- NULL
    
    # plot
    plot <-
      ggplot2::ggplot(data, ggplot2::aes(x = samples)) +
      
      ggplot2::geom_area(
        stat = "bin",
        color = "black", fill = "gray95", binwidth = binwidth, linewidth = .7
      ) +
      ggplot2::geom_vline(aes(xintercept = truth), color = ggokabeito::palette_okabe_ito(order = 1), line_width = 1) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(
          xmin = lower_mean_logit,
          xmax = upper_mean_logit,
          y = max_density * 1.15,
          height = max_density * .05
        ),
        col = "grey70",
        height = 0,
        linetype = 1,
        linewidth = 8
      ) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(
          xmin = lower_mean,
          xmax = upper_mean,
          y = max_density * 1.15,
          height = max_density * .1
        ),
        col = "black",
        linetype = 1,
        linewidth = 1.5
      ) +
      ggplot2::scale_x_continuous(
        limits = c(min, max),
        labels = c(c("0", ".25", ".50", ".75", "1")),
        breaks = seq(min, max, length.out = 5),
        expand = ggplot2::expansion()
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .05))) +
      ggplot2::labs(x = "Response Value",
                    y = "Cumulative Pointwise Frequency") +
      theme_itm() +
      theme(
        plot.margin = margin(.2, .5, .2, .2, "cm"),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "#6d6d6e", size = .3),
        axis.ticks = element_line(colour = "#6d6d6e", size = .3),
        axis.text = element_text(size = 10),
        strip.text = element_text(margin = margin(0, 0, .3, 0, "cm")),
        ) +
          
      ggplot2::theme(
        #strip.text = ggplot2::element_text(size = 14),
        #strip.background = ggplot2::element_rect(fill = "gray95"), 
        #strip.clip = "on"
      )
    # add facet wrap
    if (facet_wrap == TRUE) {
      if (is.null(design)) {
        plot <- plot + ggplot2::facet_wrap(~ item_name, scales = "free", ncol = ncol)
      } else {
        plot <- plot +
          facet_manual(~ item_name, design = design, scales = "free")
      }
    }
    
    if (show_quantiles == TRUE) {
      # add quantiles to the plot
      plot <- plot +
        ggplot2::geom_vline(
          ggplot2::aes(xintercept = median),
          col = "black",
          linetype = 1,
          linewidth = .6
        ) +
        ggplot2::geom_vline(
          ggplot2::aes(xintercept = q_05),
          col = "black",
          linetype = 2,
          linewidth = .7
        ) +
        ggplot2::geom_vline(
          ggplot2::aes(xintercept = q_95),
          col = "black",
          linetype = 2,
          linewidth = .7
        )
    }
    return(plot)
  }



# Simulation Study Helpers ------------------------------------------------
# Import location and width estimates
prep_locwid <- function(path,
                        # also obtain simplemeans?
                        simplemeans = FALSE){
  # list all rds files
  rds_files <- list.files(path, pattern = ".rds", full.names = TRUE)
  
  
  # Function to compare point estimates
  point_comparison <- function(res = res,
                               method = "model",
                               estimate_name,
                               truth_name,
                               summary,
                               pm) {
    tmp_summary <- sapply(res$results, function(x) {
      if(method == "model") {  
        x$fit_summary |>
          dplyr::filter(!grepl("_beta", variable)) |> 
          dplyr::filter(!grepl("_splx", variable)) |> 
          dplyr::filter(stringr::str_detect(variable,
                                            # ensure that it is at the beginning of the string
                                            # to avoid issues with a_loc and lambda_loc
                                            paste0("^", estimate_name))) |>
          dplyr::select(all_of({{summary}}))
      } else if(method == "simple") {
        x$simple_means |>
          dplyr::select(all_of({{estimate_name}}))
      }})
    # compare against true parameters
    ret <- list()
    for(i in 1:length(tmp_summary)){
      ret[[i]] <- pm(tmp_summary[[i]], res$results[[i]]$true_parameters[[truth_name]])
    }
    return(ret)
  }
  
  # Function to compute absolute bias
  fn_abs_bias <- function(est_param,
                          true_param,
                          average = TRUE) {
    if (isTRUE(average)) {
      mean(abs(est_param - true_param), na.rm = TRUE)
    } else {
      abs(est_param - true_param)
    }
  }
  
  l_out <- list()
  
  # loop over files to read them in, extracts estimates, computes bias 
  for(i in 1:length(rds_files)){
    # read in rds file
    res <- readRDS(rds_files[i])
    
    # compute inidividual biases
    loc_bias <- point_comparison(res = res,
                                 method = "model",
                                 estimate_name = "Tr_loc",
                                 truth_name = "Tr_loc",
                                 summary = "mean",
                                 pm = fn_abs_bias)
    wid_bias <- point_comparison(res = res,
                                 method = "model",
                                 estimate_name = "Tr_wid",
                                 truth_name = "Tr_wid",
                                 summary = "mean",
                                 pm = fn_abs_bias)
    
    if(isTRUE(simplemeans)){
      smloc_bias <- point_comparison(res = res,
                                     method = "simple",
                                     estimate_name = "simplemean_loc",
                                     truth_name = "Tr_loc",
                                     summary = "mean",
                                     pm = fn_abs_bias)
      smwid_bias <- point_comparison(res = res,
                                     method = "simple",
                                     estimate_name = "simplemean_wid",
                                     truth_name = "Tr_wid",
                                     summary = "mean",
                                     pm = fn_abs_bias)
    }
    
    
    if(isTRUE(simplemeans)){
      l_out[[i]] <- data.frame(loc_bias = unlist(loc_bias),
                               wid_bias = unlist(wid_bias),
                               smloc_bias = unlist(smloc_bias),
                               smwid_bias = unlist(smwid_bias),
                               iteration = rep(i, length(loc_bias)))
    }
    else{
      # combine into a data frame
      l_out[[i]] <- data.frame(loc_bias = unlist(loc_bias),
                               wid_bias = unlist(wid_bias),
                               iteration = rep(i, length(loc_bias)))
    }
    
    
    
    
  }
  # combine all data frames
  out <- do.call(rbind, l_out)
  return(out)
  
}


# Create common theme for all plots ---------------------------------------

theme_itm <- function(hide_axis_text_y = FALSE,
                      base_size = 12) {
  showtext_auto()
  # theme
  #ggplot2::theme_minimal(base_family = "sans") +
  #try(ggplot2::theme_minimal(base_family = "google_font")) +
  theme <-
    ggplot2::theme_minimal(base_family = "google_font", base_size = base_size) +
    ggplot2::theme(
      # remove minor grid
      panel.grid.minor = ggplot2::element_blank(),
      # Title and Axis Texts
      plot.title = ggplot2::element_text(
        face = "plain",
        size = ggplot2::rel(1.25),
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1.3), hjust = 0.5),
      axis.text.x = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.1)),
      axis.text.y = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.1)),
      axis.title.x = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.25)),
      axis.title.y = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.25),
        vjust = .3
      ),
      axis.line = element_line(colour = "#6d6d6e"),
      
      # Faceting
      strip.text = ggplot2::element_text(
        face = "plain",
        size = ggplot2::rel(1.1),
        hjust = 0.5
      ),
      strip.text.x.top = ggplot2::element_text(
        face = "plain",
        size = ggplot2::rel(1.2),
        hjust = 0.5
      ),
      # strip.text.y = element_blank(),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      # Grid
      panel.grid = ggplot2::element_line(colour = "#F3F4F5"),
      # Legend
      legend.title = ggplot2::element_text(face = "plain"),
      legend.position = "top",
      legend.justification = 1,
      # Panel/Facets
      panel.spacing.x = ggplot2::unit(1.6, "lines"),
      panel.spacing.y = ggplot2::unit(1.6, "lines"),
      # Remove vertical grid lines
      panel.grid.major.x = ggplot2::element_blank()
      
    )
  # hide y axis text
  if (isTRUE(hide_axis_text_y)) {
    theme <- theme +
      theme(axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
  } # end if
  
  return(theme)
}

