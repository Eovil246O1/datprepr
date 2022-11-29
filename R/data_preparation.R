#' Creates boolean vars for NaN observations - 1 if Nan, 0 - otherwise
#'
#' @param .dt - input data.table to which the transformation to be applied
#' @param range - range of columns of dt to which the transformation to be applied
#' @param threshold - a NaN share of input population, sufficient for creation of a boolean var
#' @param verbose - indicates if detailed messages being shown in console
#'
#' @return dt - output data.table with boolean vars for NaNs created
#' @export
#'
#' @examples
create_nan_vars = function(.dt, range = 1:ncol(.dt), threshold = .05, verbose = T){
  if(!is.data.table(.dt)) error("dt parameter is expected to be a data.table")
  dt = copy(.dt)
  col_num = ncol(dt)
  for (i in range) {
    colname = names(dt)[i]
    if (nrow(dt[is.na(eval(as.name(colname)))|is.nan(eval(as.name(colname)))|is.null(eval(as.name(colname)))|is.infinite(eval(as.name(colname)))])
        < threshold*nrow(dt)) {
      if(verbose) message("Column ", colname, ": skipped NaN variable creation as the number of NaN observations is less than ", threshold)
      next
    }
    dt[is.na(eval(as.name(colname)))|is.nan(eval(as.name(colname)))|is.null(eval(as.name(colname)))|is.infinite(eval(as.name(colname))),
          paste0(colname, '_NaN') := 1]
    dt[is.na(eval(as.name(paste0(colname, '_NaN')))), paste0(colname, '_NaN') := 0]
    if(verbose) message(paste0("Column ", colname, ": variable ", colname, '_NaN', " created..."))
  }
  return(dt)
}

#' Removes outliers from the chosen characteristics of input data table
#' using a formula Lower quantile - H, Upper quantile + H, where H = coeff * (interquantile distance)
#'
#' @param .dt - input data.table to which the transformation to be applied
#' @param range - range of columns of dt to which the transformation to be applied
#' @param probs - a vector of 2 values of lower an upper quantiles
#' @param coeff - coefficient for interquartile distance multiplication
#' @param replace - indicates an approach for outlier replacement, can be one of 'mean', 'na', 'maxmin'
#' @param create_variable - indicates a threshold for a ratio of outlier to full population, starting from which an additional variable will be created to store outlier values
#' @param verbose - indicates if detailed messages being shown in console
#'
#' @return dt - output data.table with boolean vars for NaNs created
#' @export
#'
#' @examples
remove_outliers = function(.dt, range = 1:ncol(.dt), probs = c(.25, .75), coeff = 1.5, replace = 'mean', create_variable = .05, verbose = T) {
  if(!replace %in% c('mean', 'na', 'maxmin')) stop("Parameter 'replace' should be one of 'mean', 'na', 'maxmin'")
  if(!is.data.table(.dt)) stop("dt parameter is expected to be a data.table")
  dt = copy(.dt)
  is_column_numeric = sapply(dt, is.numeric)
  
  for (i in range) {
    dt[, res := NA_integer_]
    colname = names(dt)[i]
    if(!is_column_numeric[i]) {
      if(verbose) message(paste0("Column ", colname, " is not numeric - variable will be skipped"))
      next
    }
    qnt = quantile(dt[,eval(as.name(colname))], probs, na.rm = T)
    H = coeff * (qnt[2] - qnt[1])
    
    dt[eval(as.name(colname)) < (qnt[1] - H), res := -1]
    dt[eval(as.name(colname)) > (qnt[2] + H), res := 1]
    dt[, res := dplyr::coalesce(res, 0)]
    
    if(create_variable < mean(abs(dt$res))) {
      dt[, paste0(colname, '_outliers') := abs(res)]
      if(verbose) message(paste0("Column ", paste0(colname, '_outliers'), " created"))
    }
    
    if(replace == 'maxmin') {
      dt[, eval(colname) := ifelse(res == 1,
                                   max(eval(as.name(colname))),
                                   ifelse(res == -1,
                                          min(eval(as.name(colname))),
                                          eval(as.name(colname))))]
    } else if (replace == 'na') {
      dt[res != 0, eval(colname) := NA]
    } else if (replace == 'mean') {
      dt[, eval(colname) := ifelse(res != 0, mean(eval(as.name(colname))), eval(as.name(colname)))]
    }
  }
  dt[, res := NULL]
  return(dt)
}

#' Transforms text variables of input data.table within the given range to a set of dummy variables
#'
#' @param .dt - input data.table to which the transformation to be applied
#' @param range - range of columns of dt to which the transformation to be applied
#' @param delete_origin - indicates if origin text varibles to be deleted
#' 
#' @return dt - output data.table with dummy variables for text
#' @export
#'
#' @examples

dummy_transform = function(.dt, range = 1:ncol(.dt), delete_origin = T) {
  if(!is.data.table(.dt)) stop("dt parameter is expected to be a data.table")
  dt = copy(.dt)
  nums = sapply(dt[, ..range], is.numeric)
  not_nums = dt[, ..range][, !nums, with = F]
  dummy_transform = dummyVars(~. , data = not_nums)
  dummy_data = predict(dummy_transform, not_nums)

  if (delete_origin) {
    dt = data.table(dt[, !range, with = F], dt[, ..range][, nums, with = F], dummy_data)
  } else {
    dt = data.table(dt, dummy_data)
  }
  return(dt)
}

#' Makes data replacement (for NA, NAN, NULL, INF elements).
#' For numeric parameters for replacement is used a median, for characters - mode
#'
#' @param overall_sample_num - input object to which the transformation to be applied
#' 
#' @return overall_sample_mutated - output object with dummy variables for text
#' @export
#'
#' @examples

replace_nan = function(overall_sample_num) {
  overall_sample_mutated = overall_sample_num %>%
    mutate_if(is.numeric,
              .funs = ~
                ifelse(is.na(.)|is.nan(.)|is.null(.)|is.infinite(.),
                       median(., na.rm = TRUE),
                       .)) %>%
    mutate_if(is.character,
              .funs = ~
                ifelse(is.na(.)|is.nan(.)|is.null(.)|is.infinite(.),
                       mode_value(.),
                       .))
  return(overall_sample_mutated)
}

#' Get the ‘mode’ (a kind of value) of an R object.
#'
#' @param x - input object
#' @param na_rm - indicates if NA, NaN, Inf values should be removed. Default is TRUE
#' 
#' @return res - output value of a statistic mode
#' @export
#'
#' @examples

mode_value = function(x, na_rm = T) {
  if (na_rm) {
    ux = unique(x[!is.na(x)&!is.nan(x)&!is.null(x)&!is.infinite(x)])
  } else {
    ux = unique(x)
  }
  res = ux[which.max(tabulate(match(x, ux)))]
  return(res)
}

#' Scale input dataset.
#'
#' @param input - input object
#' @param range - range of columns to which the transformation to be applied
#' 
#' @return res - output object
#' @export
#'
#' @examples

scale_custom = function(input, range = 1:ncol(input)){
  range_num = sapply(input[, range], is.numeric)
  res = as.data.frame(apply(input[, range][, range_num],  2, function (x) (x - mean(x)) / sd(x)))
  res = cbind(res, input[, range][!range_num])
  return(res)
}

#' Determine number of clusters
#' http://www.bagualu.net/wordpress/wp-content/uploads/2015/10/A_Handbook_of_Statistical_Analyses_Using_R__Second_Edition.pdf
#' page 329
#'
#' @param input - input object
#' @param init_claster_num - initial number of clusters
#' 
#' @return count_claster - evaluated optimal count of clusters
#' @export
#'
#' @examples

cluster_determine = function(input, init_claster_num = 15) {
  wss = rep(0, init_claster_num)
  for (i in 1:init_claster_num) wss[i] = sum(kmeans(input[,-c(1,2, ncol(input))], centers = i)$withinss)
  plot(1:init_claster_num, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

  #choose the optimal cluster number (when margin weighted squares effect higher than margin growth of cluster count)
  for (i in 1:length(wss)) {
    if (-diff(wss, lag = 1)[i]/max(wss) > 1/length(wss)) {
      count_claster = i + 1
      print(c(count_claster, wss, diff(wss, lag = 1)[i]/max(wss)))
    }
    else {break}
  }
  return(count_claster)
}
