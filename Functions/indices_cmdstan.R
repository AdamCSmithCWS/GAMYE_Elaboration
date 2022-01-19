## generic annual indices function for
## cmdstanr fit object that includes some
## parameters that estimate annual indices of abundance

index_function <- function(fit = stanfit,
                           parameter = "n",
                           strat = NULL,#"Stratum_Factored",
                           year = "Year",
                           quant = 0.95,
                           weights_df = NULL,
                           area = NULL,#"Area",
                           summary_regions = NULL # optional column to provide summary regions, if not supplied summarises across all strata
                            ){
  
  lu <- ((1-(quant))/2)
  uu <- 1-((1-(quant))/2)
  
  if(is.null(strat)){
    dims <- year

  }else{
    if(!is.null(weights_df)){
      if(is.null(area)){
      stop("Weights were supplied but no stratum dimension")
      
      return(NULL)
      }
      if(!any(grepl(pattern = strat,x = names(weights_df)))){
        stop("weights_df must include a column with a name that matches strat")
        return(NULL)
      }
    }

    dims <- c(strat,year)
  }
  
  smpls <- posterior_samples(fit = fit,
                             parm = parameter,
                             dims = dims)
    
  if(!is.null(weights_df)){
    # weights_df <- data.frame(Stratum_Factored = 1:nstrata_fit,
    # Area = abs(rnorm(nstrata_fit,1000,500)),
    # province = rep_len(c("Ontario","Quebec","Saskatchewan"),
    #                nstrata_fit))
    nstrata_fit <- length(unique(smpls[[strat]]))
    nstrata_w <- nrow(weights_df)
    if(nstrata_fit != nstrata_w){
      stop("Lengths of strata and weights are different")
      return(NULL)
    }
    if(!is.null(summary_regions)){
    
    weights_df <- weights_df %>% 
      rename_with(.,~gsub(pattern = area,
                          replacement = "a",
                  x = .x,
                  fixed = TRUE)) %>% 
    rename_with(., ~gsub(pattern = summary_regions,
                         replacement = "rrr",.x,
                         fixed = TRUE))
    tmp <- weights_df %>% 
      group_by(rrr) %>% 
    summarise(suma = sum(a),.groups = "keep")
    
    weights_df <- weights_df %>% 
      left_join(tmp,by = "rrr") %>% 
      mutate(w = a/suma)
    
    
    }else{
      weights_df <- weights_df %>% 
        rename_with(.,~gsub(pattern = area,
                            replacement = "a",
                            x = .x,
                            fixed = TRUE)) %>% 
        mutate(w = a/sum(a))  
    }

    smpls <- smpls %>% 
      left_join(.,weights_df,
                by = strat)
    
    if(!is.null(summary_regions)){
      inds <- smpls %>% 
        mutate(.value = .value*w) %>% 
        rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                             fixed = TRUE)) %>% 
        group_by(yyy,rrr,.draw) %>% 
        summarise(.vsum = sum(.value)) %>% 
        group_by(yyy,rrr) %>% 
        summarise(mean = mean(.vsum),
                  median = median(.vsum),
                  lci = quantile(.vsum,lu),
                  uci = quantile(.vsum,uu),
                  .groups = "keep") %>% 
        rename_with(., ~gsub(replacement = summary_regions,pattern = "rrr",.x,
                              fixed = TRUE)) %>% 
        rename_with(., ~gsub(replacement = year,pattern = "yyy",.x,
                             fixed = TRUE))
      
    }else{
    inds <- smpls %>% 
      mutate(.value = .value*w) %>% 
      #rename_with(., ~gsub(pattern = strat,replacement = "s",.x,
                          # fixed = TRUE)) %>% 
      rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                           fixed = TRUE)) %>% 
      group_by(yyy,.draw) %>% 
      summarise(.vsum = sum(.value)) %>% 
      group_by(yyy) %>% 
      summarise(mean = mean(.vsum),
                median = median(.vsum),
                lci = quantile(.vsum,lu),
                uci = quantile(.vsum,uu),
                .groups = "keep") %>% 
      # rename_with(., ~gsub(replacement = strat,pattern = "s",.x,
      #                      fixed = TRUE)) %>% 
      rename_with(., ~gsub(replacement = year,pattern = "yyy",.x,
                           fixed = TRUE))
    
}#end summary regions
  }else{
  inds <- smpls %>% 
    rename_with(., ~gsub(pattern = strat,replacement = "sss",.x,
                         fixed = TRUE)) %>% 
    rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                         fixed = TRUE)) %>% 
    group_by(sss,yyy) %>% 
    summarise(mean = mean(.value),
              median = median(.value),
              lci = quantile(.value,lu),
              uci = quantile(.value,uu),
              .groups = "keep") %>% 
    rename_with(., ~gsub(replacement = strat,pattern = "sss",.x,
                         fixed = TRUE)) %>% 
    rename_with(., ~gsub(replacement = year,pattern = "yyy",.x,
                         fixed = TRUE))
}
  
  return(list = c("inds","smpls"))
  }



# helper functions --------------------------------------------------------



### function to extract the dimension values from an bayesian model fit
### works within the gather_samples function
dim_ext <- function(dim = 1,
                    var = "",
                    cl = "Parameter",
                    dat = NULL){
  ##3 function to extract the indicator values from cmdstanr output
  require(stringr)
  
  pat = paste0("(?<=",var,"\\[")
  
  if(dim > 1){
    for(j in 1:(dim-1)){
      
      pat2 = paste0(pat,")[:digit:]+")
      cl2 = str_extract(unlist(dat[,cl]),pattern = pat2)
      
      d = max(nchar(cl2))
      
      pat = paste0(pat,"[:digit:]{1,",d,"}[:punct:]")
    }
  }
  
  
  pat = paste0(pat,")[:digit:]+")
  dds = as.integer(str_extract(unlist(dat[,cl]),pattern = pat))
  return(dds)
  
}


### function to generate the same tidy output as gather-draws in tidybayes package
## dims should be a character vector defining the dimensions of the parameter
## e.g., parm = "nsmooth", dims = c("stratum","year"),
## function works with cmdstanr output by default and rstan fits
## with is_rstan == TRUE
posterior_samples <- function(fit = cmdstanfit,
                              parm = "nsmooth",
                              dims = NULL,
                              is_rstan = FALSE,
                              is_mcmc = FALSE){
  require(posterior)
  require(tidyverse)
  
  if(length(dims) > 0){
    parm_ex <- paste0(parm,"[")
  }else{
    parm_ex <- parm
  }
  if(class(fit)[1] == "stanfit"){is_rstan <- TRUE}
  
  if(class(fit)[1] == "mcmc"){is_mcmc <- TRUE}
  
  if(is_rstan | is_mcmc){
    samples <- as_draws_df(as.array(fit)) %>% 
      dplyr::select(starts_with(parm_ex,ignore.case = FALSE),
                    .chain,
                    .iteration,
                    .draw)#,pars = c(parm)))
  }else{
    
    samples <- as_draws_df(fit$draws(variables = c(parm)))
    
  }
  
  
  
  plong <- suppressWarnings(samples %>% pivot_longer(
    cols = starts_with(parm_ex,ignore.case = FALSE),
    names_to = c(".variable"),
    values_to = ".value",
    values_drop_na = TRUE
  )) 
  
  for(dn in 1:length(dims)){
    dd = dims[dn]
    plong[,dd] = dim_ext(dim = dn,
                         var = parm,
                         cl = ".variable",
                         dat = plong)
    
  }
  
  plong <- plong %>% mutate(.variable = parm)
  return(plong)
  
}





posterior_sums <- function(samples = n_samples,
                           quantiles = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),
                           ci = 0.95,
                           dims = NULL){
  
  cis = c((1-ci)/2,1-((1-ci)/2))
  
  if(!is.null(dims)){
    for(i in 1:length(dims)){
      samples <- samples %>% 
        rename_with(~gsub(dims[i],paste0("d",i),.x,fixed = TRUE))
    }
    if(length(dims) == 1){
      sums = samples %>% 
        group_by(d1) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value,cis[1]),
                  uci = quantile(.value,cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))
      
      if(!is.null(quantiles)){
        qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
        for(i in 1:length(quantiles)){
          qq = quantiles[i]
          qn = qs[i]
          
          sumt = samples %>% 
            group_by(d1) %>% 
            summarise(tt = as.numeric(quantile(.value,qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
            rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) 
          
          sums = left_join(sums,sumt,by = dims)
        }
        
      }
      
    }
    if(length(dims) == 2){
      sums = samples %>% 
        group_by(d1,d2) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value,cis[1]),
                  uci = quantile(.value,cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))
      
      
      if(!is.null(quantiles)){
        qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
        for(i in 1:length(quantiles)){
          qq = quantiles[i]
          qn = qs[i]
          
          sumt = samples %>% 
            group_by(d1,d2) %>% 
            summarise(tt = as.numeric(quantile(.value,qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
            rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))
          
          
          sums = left_join(sums,sumt,by = dims)
        }
        
      }
      
      
    }
    
    if(length(dims) == 3){
      sums = samples %>% 
        group_by(d1,d2,d3) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value,cis[1]),
                  uci = quantile(.value,cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d2",dims[2],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d3",dims[3],.x,fixed = TRUE))
      
      
      if(!is.null(quantiles)){
        qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
        for(i in 1:length(quantiles)){
          qq = quantiles[i]
          qn = qs[i]
          
          sumt = samples %>% 
            group_by(d1,d2,d3) %>% 
            summarise(tt = as.numeric(quantile(.value,qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
            rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d3",dims[3],.x,fixed = TRUE))
          
          
          sums = left_join(sums,sumt,by = dims)
        }
        
      }
      
      
    }
    
    
    if(length(dims) == 4){
      sums = samples %>% 
        group_by(d1,d2,d3,d4) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value,cis[1]),
                  uci = quantile(.value,cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d2",dims[2],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d3",dims[3],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d4",dims[4],.x,fixed = TRUE))
      
      
      if(!is.null(quantiles)){
        qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
        for(i in 1:length(quantiles)){
          qq = quantiles[i]
          qn = qs[i]
          
          sumt = samples %>% 
            group_by(d1,d2,d3,d4) %>% 
            summarise(tt = as.numeric(quantile(.value,qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
            rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d3",dims[3],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d4",dims[4],.x,fixed = TRUE))
          
          
          sums = left_join(sums,sumt,by = dims)
        }
        
      }
      
      
    }
    
    if(length(dims) > 4){stop("ERROR this function cannot handle more than 4 dimensions, but it could be easily modified")}
    
  }else{
    
    
    sums = samples %>% 
      summarise(mean = mean(.value),
                median = median(.value),
                sd = sd(.value),
                lci = quantile(.value,cis[1]),
                uci = quantile(.value,cis[2]),
                .groups = "keep")
    
    if(!is.null(quantiles)){
      qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
      for(i in 1:length(quantiles)){
        qq = quantiles[i]
        qn = qs[i]
        
        sumt = samples %>% 
          summarise(tt = as.numeric(quantile(.value,qq)),
                    .groups = "keep") %>% 
          rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE))
        sums = bind_cols(sums,sumt)
      }
      
    }
    
  }
  return(sums)
  
}









