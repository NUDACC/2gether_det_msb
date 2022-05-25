
###------------------------------###
###------------------------------###
### Minimal Sufficient Balance
### Adaptive Randomization
###------------------------------###
###------------------------------###


###--------------------###
### Libraries
###--------------------###
#library(tidyverse)

###--------------------###
### Functions to compute
### balance
###--------------------###



#' @name get_votes
#' @param center string indicating variable name for center
#' @param covariates list of strings indicating variables on which to conduct balance
#' @param treamtent string indicating treatment assignment variable
#' @param min_n_adapt integer indicating # of already randomized individuals required to run adaptive algorithm
#' @param prob_vote probability of assigning to treatment if it would improve balance
#' @return treatment assignment probability
#' @note 1 = treatment
get_votes <- function(data, new_data, # data
                      covariates, treatment, # variable names
                      single = F,
                      min_n_adapt = 10,
                      prob_vote = 0.7, # probability split for majority arm
                      show_votes = F){ # output options
  
  
  if(nrow(data) < min_n_adapt){
    
    prob = .5
    
    if(show_votes){
      
      return(
        list(
          prob = prob,
          votes = NULL,
          majority = NULL
        )
      )
      
    } else {
      
      return(prob)
      
    }
  } else {
    overall_vote <- list()
    # center_vote <- list()
    mean_vote <- list()
    
    dtrt <- data %>%
      summarize(pct_trt = mean(get(treatment)))
    
    vote <- ifelse(dtrt > 0.5, "Arm 0",
                   ifelse(dtrt < 0.5, "Arm 1", "Neutral"))
    
    overall_vote[["overall"]] <- vote
    
    # new_center <- new_data %>% select(get(center)) %>% pull()
    
    # Summarize assignment within center of interest
    # dc <- data %>%
    #   filter(.data[[center]] == new_center) %>%
    #   summarize(x = sum(get(treatment)), n = n())
    
    # Check if this is a brand new center
    # or if patients have already been
    # randomized in this center.
    # if(dc$n == 0){
    #   
    #   vote <- "Neutral"
    #   center_vote[[paste("center")]] <- vote
    #   
    #   
    # } else {
    #   
    #   # Get p-value
    #   bt <- binom.test(dc$x, dc$n)
    #   
    #   # Compute and return vote for arm
    #   v0 <- (bt$p.value < 0.3 & dc$x/dc$n < 0.5)
    #   v1 <- (bt$p.value < 0.3 & dc$x/dc$n > 0.5)
    #   vote <- ifelse(v0, "Arm 0",
    #                  ifelse(v1, "Arm 1", "Neutral"))
    #   center_vote[["center"]] <- vote
    # }
    # 
    for(j in covariates){
      
      if(nrow(data) == 1){
        
        vote <- "Neutral"
        
      } else if(is.numeric(data[[j]])){
        
        if(nrow(data) <= 4 | length(unique(data[,paste(treatment)])) == 1){
          
          vote <- "Neutral"
          
        } else {
          
          # Compute balance using linear model
          mod <- lm(data[[j]] ~ data[[treatment]])
          res <- summary(mod)$coefficients
          
          # Formulate output as tibble
          out <- tibble(
            var = j, # variable name
            diff = res[2, "Estimate"], # Difference
            se = res[2, "Std. Error"], # SE of difference
            stat = res[2, "t value"], # test statistic
            p = res[2, "Pr(>|t|)"] # p-value
          )
          
          x_sum <- data %>%
            group_by(get(treatment)) %>%
            summarize(x = mean(get(j), na.rm = T)) %>%
            rename(treatment = `get(treatment)`)
          
          x1 <- x_sum %>%
            filter(treatment == 1) %>%
            pull(x)
          
          x0 <- x_sum %>%
            filter(treatment == 0) %>%
            pull(x)
          
          # Get record IDs value
          new_col_name <- new_data %>% pull(get(j))
          
          if(single){
              
            # Compute whether new single patient gets voted to arm 0 or arm 1
            v1 <- ((out$p < 0.3) & (out$diff > 0) & (new_col_name < x0)) |
              ((out$p < 0.3) & (out$diff < 0) & (new_col_name > x0))
            v0 <- ((out$p < 0.3) & (out$diff > 0) & (new_col_name < x1)) |
              ((out$p < 0.3) & (out$diff < 0) & (new_col_name > x1))
          
          } else {
            
            # Get partner's value
            new_col_name_p <- new_data %>% pull(get(paste0(j, "_p")))
            
            v1 <- ((out$p < 0.3) & (out$diff > 0) & ((new_col_name + new_col_name_p)/2 < x0)) |
              ((out$p < 0.3) & (out$diff < 0) & ((new_col_name + new_col_name_p)/2 > x0))
            v0 <- ((out$p < 0.3) & (out$diff > 0) & ((new_col_name + new_col_name_p)/2 < x1)) |
              ((out$p < 0.3) & (out$diff < 0) & ((new_col_name + new_col_name_p)/2 > x1))
              
          }
          
          # Return vote
          vote <- ifelse(v1, "Arm 1", ifelse(v0, "Arm 0", "Neutral"))
          
        }
        
      } else if(is.factor(data[[j]])){
        
        if((length(unique(data[[j]])) == 1) |
           (length(unique(data[[treatment]])) == 1)){
          
          vote <- "Neutral"
          
        } else {
          
          mod <- chisq.test(data[[j]], data[[treatment]])
          
          # Format output as table
          out <- tibble(
            var = j,
            observed = list(mod$observed),
            expected = list(mod$expected),
            diff = list(mod$observed - mod$expected),
            stat = mod$statistic,
            p = mod$p.value)
          
          new_col_name <- new_data %>% pull(get(j)) %>% as.character()
          
          # If this is not the first person with a given category, 
          # they get votes as normal.
          if(new_col_name %in% rownames(out$expected[[1]])){
            v0 <- (out$p < 0.3 & 
                     out$expected[[1]][new_col_name, 1] > out$observed[[1]][new_col_name, 1])
            v1 <- (out$p < 0.3 & 
                     out$expected[[1]][new_col_name, 2] > out$observed[[1]][new_col_name, 2])
          } else {
            # If this is the first person with a given category, they vote neutral
            
            v0 <- FALSE
            v1 <- FALSE
            
          }
          
          
          
          # For partnered participants, compute votes for each person. 
          # Then, reconcile if they vote in the same direction.
          if(!single){
            
            new_col_name_p <- new_data %>% 
              pull(get(paste0(j, "_p"))) %>%
              as.character()
            
            # person 1 votes
            v01 <- v0
            v11 <- v1
            
            # person 2 votes
            v02 <- (out$p < 0.3 & 
                     out$expected[[1]][new_col_name_p, 1] > out$observed[[1]][new_col_name_p, 1])
            v12 <- (out$p < 0.3 & 
                     out$expected[[1]][new_col_name_p, 2] > out$observed[[1]][new_col_name_p, 2])
            
            
            if(v11 & v12){ # both vote for treatment?
              v1 <- TRUE; v0 <- FALSE
            } else if(v01 & v02){ # both vote for control?
              v0 <- TRUE; v1 <- FALSE
            } else { # If partner votes conflict, vote neutral
              v1 <- FALSE; v0 <- FALSE
            }
           
            
            
          }
          
          vote <- ifelse(v1, "Arm 1", ifelse(v0, "Arm 0", "Neutral"))
          
        }
        
      }
      
      mean_vote[[paste(j)]] <- vote
      
    }
    
    total_vote <- t(as.data.frame(c(overall_vote, mean_vote)))
    #return(t(as.data.frame(total_vote)))
    
    vote_tab <- as.data.frame(table(total_vote))
    
    vt <- vote_tab %>% filter(total_vote != "Neutral")
    
    if(nrow(vt) == 0){
      
      majority = NULL
      
    } else {
      
      majority <- as.character(vt$total_vote[vote_tab$Freq == max(vt$Freq)])
      majority <- majority[!is.na(majority)]
      
    }
    
    
    if(is.null(majority) | length(majority) > 1){
      prob <- 0.5
    } else if(majority == "Neutral"){
      prob <- 0.5
    } else if(majority == "Arm 1"){
      prob <- prob_vote
    } else {
      prob <- 1 - prob_vote
    }
    
    if(show_votes){
      
      return(
        list(
          prob = prob,
          votes = total_vote,
          majority = majority
        )
      )
      
    } else {
      
      return(prob)
      
    }
  }
  
}

# prob <- get_votes(center, covariates, treatment, data, new_data)

#' @name generate_assignment
#' @param prob numeric [0,1]: Probability of assignment to arm 1 (treatment)
#' @return treatment assigment as integer
#' @note 1 = treatment
generate_assignment <- function(prob){
  return(
    rbinom(1, 1, prob)
  )
}
