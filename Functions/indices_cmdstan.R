## generic annual indices function for
## cmdstanr fit object that includes some
## parameters that estimate annual indices of abundance


# weights_df <- data.frame(Stratum_Factored = 1:nstrata_fit,
# Area = abs(rnorm(nstrata_fit,1000,500)),
# province = rep_len(c("Ontario","Quebec","Saskatchewan"),
#                nstrata_fit))
# ttinds <- index_function(fit = stanfit,
#                       parameter = "n",
#                       strat = "Stratum_Factored",
#                       year = "Year",
#                       first_dim = "s",
#                       quant = 0.95,
#                       weights_df = weights_df,
#                       area = "Area",
#                       summary_regions = "province",
#                      year_1 = 1966)
# 



index_function <- function(fit = stanfit,
                           parameter = "n",
                           strat = NULL,#"Stratum_Factored",
                           year = "Year",
                           first_dim = "s",
                           quant = 0.95,
                           weights_df = NULL, # ifsupplied a composite summary estimate is generated
                           area = NULL,#"Area",
                           summary_regions = NULL, # optional column to provide summary regions, if not supplied summarises across all strata if weights_df is upplied
                           year_1 = 1970,
                           to_summarise = FALSE # if TRUE provides summary indices across regions, otherwise jsut strata estimates
                            ){
  
  lu <- ((1-(quant))/2)
  uu <- 1-((1-(quant))/2)
  year_1 <- year_1-1
  if(is.null(strat)){
    dims <- year
    
  }else{
    if(!is.null(weights_df)){
      if(is.null(strat)){
      stop("Weights were supplied but no stratum dimension")
      
      return(NULL)
      }
      if(!any(grepl(pattern = strat,x = names(weights_df)))){
        stop("weights_df must include a column with a name that matches strat")
        return(NULL)
      }
      if(!to_summarise){
        to_summarise <- TRUE
        warning("weights were supplied by to_summarise is FALSE, changing to_summaries to TRUE")
      }
    }
    if(to_summarise){
     if(is.null(weights_df) | is.null(strat)){
       stop("to_summarise is TRUE but no weights were supplied or strat is NULL")
       return(NULL)
     }
    }
  if(tolower(substr(first_dim,1,1)) == "s"){
    dims <- c(strat,year)
    }
    if(tolower(substr(first_dim,1,1)) == "y"){
    dims <- c(year,strat)
    }
  }
  
  smpls <- posterior_samples(fit = fit,
                             parm = parameter,
                             dims = dims)
    
  if(!is.null(weights_df)){
    nstrata_fit <- length(unique(smpls[[strat]]))
    nstrata_w <- nrow(weights_df)
    if(nstrata_fit != nstrata_w){
      stop("Lengths of strata and weights are different")
      return(NULL)
    }

   
    if(!is.null(summary_regions)){
    
    weights_df <- weights_df %>% 
      rename_with(.,~gsub(pattern = area,
                          replacement = "stratum_area",
                  x = .x,
                  fixed = TRUE)) %>% 
    rename_with(., ~gsub(pattern = summary_regions,
                         replacement = "summary_region",.x,
                         fixed = TRUE))
    tmp <- weights_df %>% 
      group_by(summary_region) %>% 
    summarise(regional_area_sum = sum(stratum_area),.groups = "keep")
    
    weights_df <- weights_df %>% 
      left_join(.,tmp,by = "summary_region") %>% 
      mutate(stratum_weight = stratum_area/regional_area_sum)
    
    
    }else{
      weights_df <- weights_df %>% 
        rename_with(.,~gsub(pattern = area,
                            replacement = "stratum_area",
                            x = .x,
                            fixed = TRUE)) %>% 
        mutate(stratum_weight = stratum_area/sum(stratum_area),
               summary_region = "Survey_wide") 
      summary_regions <- "Survey_wide"
      
      }

    smpls <- smpls %>% 
      left_join(.,weights_df,
                by = strat)
    
   # if(summary_regions != "Survey_wide"){
      inds <- smpls %>% 
        mutate(.value = .value*stratum_weight) %>% 
        rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                             fixed = TRUE)) %>% 
        group_by(yyy,summary_region,.draw) %>% 
        summarise(.vsum = sum(.value),
                  .groups = "keep") %>% 
        group_by(yyy,summary_region) %>% 
        summarise(mean = mean(.vsum),
                  median = median(.vsum),
                  lci = quantile(.vsum,lu),
                  uci = quantile(.vsum,uu),
                  .groups = "keep") %>% 
        mutate(true_year = yyy+year_1) %>% 
        rename_with(., ~gsub(replacement = summary_regions,pattern = "summary_region",.x,
                              fixed = TRUE)) %>% 
        rename_with(., ~gsub(replacement = year,pattern = "yyy",.x,
                             fixed = TRUE))
      
      smpls <- smpls %>% 
        mutate(.value = .value*stratum_weight) %>% 
        rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                             fixed = TRUE)) %>% 
        group_by(yyy,summary_region,.draw) %>% 
        summarise(.value = sum(.value),
                  .groups = "keep") %>%
        rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                             fixed = TRUE)) %>% 
        rename_with(., ~gsub(replacement = summary_regions,
                             pattern = "summary_region",.x,
                             fixed = TRUE)) %>% 
        mutate(true_year = yyy+year_1) %>% 
        rename_with(., ~gsub(pattern = "yyy",replacement = year,.x,
                             fixed = TRUE))
      
      

      
      
      weights_df <- weights_df %>% 
        rename_with(., ~gsub(replacement = area,
                             pattern = "stratum_area",.x,
                             fixed = TRUE))%>% 
        rename_with(., ~gsub(replacement = summary_regions,
                             pattern = "summary_region",.x,
                             fixed = TRUE))
      
    # }else{
    # inds <- smpls %>% 
    #   mutate(.value = .value*stratum_weight) %>% 
    #   #rename_with(., ~gsub(pattern = strat,replacement = "s",.x,
    #                       # fixed = TRUE)) %>% 
    #   rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
    #                        fixed = TRUE)) %>% 
    #   group_by(yyy,.draw) %>% 
    #   summarise(.vsum = sum(.value),
    #             .groups = "keep") %>% 
    #   group_by(yyy) %>% 
    #   summarise(mean = mean(.vsum),
    #             median = median(.vsum),
    #             lci = quantile(.vsum,lu),
    #             uci = quantile(.vsum,uu),
    #             .groups = "keep") %>%
    #   mutate(true_year = yyy+year_1,
    #          summary_region = "Survey_wide") %>%  
    #   # rename_with(., ~gsub(replacement = strat,pattern = "s",.x,
    #   #                      fixed = TRUE)) %>% 
    #   rename_with(., ~gsub(replacement = year,pattern = "yyy",.x,
    #                        fixed = TRUE))
    # 
    # weights_df <- weights_df %>% 
    #   rename_with(., ~gsub(replacement = area,
    #                        pattern = "stratum_area",.x,
    #                        fixed = TRUE)) %>% 
    #   mutate(summary_region = "Survey_wide")
    # summary_regions <- "Survey_wide"
#}#end summary regions
    
  }else{ ## everything up to here requires a weights_df data frame
    if(!is.null(strat)){
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
    mutate(true_year = yyy+year_1) %>% 
    rename_with(., ~gsub(replacement = strat,pattern = "sss",.x,
                         fixed = TRUE)) %>% 
    rename_with(., ~gsub(replacement = year,pattern = "yyy",.x,
                         fixed = TRUE))
  summary_regions <- strat
  
  smpls <- smpls %>% 
    rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                         fixed = TRUE)) %>% 
    rename_with(., ~gsub(replacement = strat,pattern = "sss",.x,
                         fixed = TRUE)) %>% 
    mutate(true_year = yyy+year_1) %>% 
    rename_with(., ~gsub(pattern = "yyy",replacement = year,.x,
                         fixed = TRUE))
  
    }else{
      summary_regions <- "Complete_pooling"

      inds <- smpls %>%
        rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                             fixed = TRUE)) %>%
        group_by(yyy) %>%
        summarise(mean = mean(.value),
                  median = median(.value),
                  lci = quantile(.value,lu),
                  uci = quantile(.value,uu),
                  .groups = "keep") %>%
        mutate(true_year = yyy+year_1,
               Survey_wide = "Survey_wide") %>%
        rename_with(., ~gsub(replacement = year,pattern = "yyy",.x,
                             fixed = TRUE))


      smpls <- smpls %>%
        rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                             fixed = TRUE)) %>%
        mutate(true_year = yyy+year_1,
               Survey_wide = "Survey_wide") %>%
        rename_with(., ~gsub(pattern = "yyy",replacement = year,.x,
                             fixed = TRUE))

}
}
  
inds <- inds %>% 
  mutate(Region_type = summary_regions)

  return(list(indices = inds,
              samples = smpls,
              parameter = parameter,
              strat = strat,
              year = year,
              dims = dims,
              quant = quant,
              weights_df = weights_df,
              area = area,#"Area",
              summary_regions = summary_regions,
              to_summarise = to_summarise,
              year_1 = year_1+1# optional column to provide summary regions, if not supplied summarises across all strata
  ))
}



# Trends function ---------------------------------------------------------
#tr_tmp <- trends_function(ind_list = ttinds)
trends_function <- function(ind_list = ind_list,
                            start_year = NULL,
                            end_year = NULL,
                            quant = 0.95){

  indices = ind_list$indices
  samples = ind_list$samples
  parameter = ind_list$parameter
  strat = ind_list$strat
  year = ind_list$year
  dims = ind_list$dims
  weights_df = ind_list$weights_df
  area = ind_list$area
  summary_regions = ind_list$summary_regions
  to_summarise = ind_list$to_summarise
  
  if(is.null(end_year)){
    end_year <- max(samples$true_year)
  }
  if(is.null(start_year)){
    start_year <- min(samples$true_year)
  }
  
  nyrs <- end_year-start_year
  lu <- ((1-(quant))/2)
  uu <- 1-((1-(quant))/2)
  
  if(!is.null(weights_df) & to_summarise){
    
    
    indt <- samples %>% 
      filter(true_year %in% c(start_year,end_year)) %>% 
      ungroup() %>% 
      select(-matches(match = year,ignore.case = FALSE)) %>% 
      pivot_wider(names_from = true_year,
                  values_from = .value,
                  names_prefix = "Y") %>% 
      rename_with(., ~gsub(replacement = "start",
                           pattern = paste0("Y",start_year),.x,
                           fixed = TRUE))%>% 
      rename_with(., ~gsub(replacement = "end",
                           pattern = paste0("Y",end_year),.x,
                           fixed = TRUE))%>% 
      rename_with(., ~gsub(replacement = "stratum_trend",
                           pattern = summary_regions,.x,
                           fixed = TRUE))
    
    
    
    tt <- indt %>% 
      group_by(.draw,stratum_trend) %>% 
      summarise(end = sum(end),
                start = sum(start),
                t = texp(end/start,ny = nyrs),
                ch = chng(end/start),
                .groups = "keep") %>% 
      group_by(stratum_trend) %>% 
      summarise(trend = mean(t),
                lci = quantile(t,lu,names = FALSE),
                uci = quantile(t,uu,names = FALSE),
                percent_change = median(ch),
                p_ch_lci = quantile(ch,lu,names = FALSE),
                p_ch_uci = quantile(ch,uu,names = FALSE),
                prob_decline = prob_dec(ch,0),
                prob_decline_GT30 = prob_dec(ch,-30),
                prob_decline_GT50 = prob_dec(ch,-50),
                prob_decline_GT70 = prob_dec(ch,-70))%>% 
      rename_with(., ~gsub(replacement = summary_regions,
                           pattern = "stratum_trend",.x,
                           fixed = TRUE))
    
    
    
    
  }else{ #else is.null weights_df
    
    if(!is.null(strat)){
    indt <- samples %>% 
      filter(true_year %in% c(start_year,end_year)) %>% 
      #ungroup() %>% 
      select(-matches(match = year,ignore.case = FALSE)) %>% 
      pivot_wider(names_from = true_year,
                  values_from = .value,
                  names_prefix = "Y") %>% 
      rename_with(., ~gsub(replacement = "start",
                           pattern = paste0("Y",start_year),.x,
                           fixed = TRUE))%>% 
      rename_with(., ~gsub(replacement = "end",
                           pattern = paste0("Y",end_year),.x,
                           fixed = TRUE))%>% 
      rename_with(., ~gsub(replacement = "stratum_trend",
                           pattern = strat,.x,
                           fixed = TRUE))
    
    
    
    tt <- indt %>% 
      group_by(.draw,stratum_trend) %>% 
      summarise(t = texp(end/start,ny = nyrs),
                ch = chng(end/start),
                .groups = "keep") %>% 
      group_by(stratum_trend) %>% 
      summarise(trend = mean(t),
                lci = quantile(t,lu,names = FALSE),
                uci = quantile(t,uu,names = FALSE),
                percent_change = median(ch),
                p_ch_lci = quantile(ch,lu,names = FALSE),
                p_ch_uci = quantile(ch,uu,names = FALSE),
                prob_decline = prob_dec(ch,0),
                prob_decline_GT30 = prob_dec(ch,-30),
                prob_decline_GT50 = prob_dec(ch,-50),
                prob_decline_GT70 = prob_dec(ch,-70))%>% 
      rename_with(., ~gsub(replacement = strat,
                           pattern = "stratum_trend",.x,
                           fixed = TRUE))
    
    }else{
      
      indt <- samples %>% 
        filter(true_year %in% c(start_year,end_year)) %>% 
        #ungroup() %>% 
        select(-matches(match = year,ignore.case = FALSE)) %>% 
        pivot_wider(names_from = true_year,
                    values_from = .value,
                    names_prefix = "Y") %>% 
        rename_with(., ~gsub(replacement = "start",
                             pattern = paste0("Y",start_year),.x,
                             fixed = TRUE))%>% 
        rename_with(., ~gsub(replacement = "end",
                             pattern = paste0("Y",end_year),.x,
                             fixed = TRUE))
      
      
      
      tt <- indt %>% 
        group_by(.draw) %>% 
        summarise(t = texp(end/start,ny = nyrs),
                  ch = chng(end/start),
                  .groups = "keep") %>% 
        summarise(trend = mean(t),
                  lci = quantile(t,lu,names = FALSE),
                  uci = quantile(t,uu,names = FALSE),
                  percent_change = median(ch),
                  p_ch_lci = quantile(ch,lu,names = FALSE),
                  p_ch_uci = quantile(ch,uu,names = FALSE),
                  prob_decline = prob_dec(ch,0),
                  prob_decline_GT30 = prob_dec(ch,-30),
                  prob_decline_GT50 = prob_dec(ch,-50),
                  prob_decline_GT70 = prob_dec(ch,-70))
    }
    
    
  }
  
  tt <- tt %>% 
    mutate(Region_type = summary_regions)
  return(tt)
}


# helper functions --------------------------------------------------------

p_lt <- function(x,th){
  length(which(x < th))/length(x)
}


p_neg <- function(x){
  length(which(x < 0))/length(x)
}

texp <- function(x,ny = 2019-1974){
  (x^(1/ny)-1)*100
}




chng <- function(x){
  (x-1)*100
}

prob_dec <- function(ch,thresh){
  
  length(which(ch < thresh))/length(ch)
}

# 
# slope_trend <- function(x,y){
#   x = log(x)
#   n = length(y)
#   sx = sum(x)
#   sy = sum(y)
#   ssy = sum(y^2)
#   sxx = sum(x*y)
#   b = (n*sxx - sy*sx)/(n*ssy - sy^2)
#   return(b)
# }
# 
# 
# sltexp <- function(x){((exp(x)-1)*100)} 


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









