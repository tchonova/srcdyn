
### Prioritise and rank profiles ####

#' Prioritise and rank profiles
#'
#' This function takes a dataframe with time profiles containing profile IDs (numeric or character),
#' dates (date format) and intensity (numeric) and returns profiles prioritised based on the following
#' criteria: maximum intensity, cumulative intensity, frequency of detection.
#'
#' @param x a tibble in long format.
#' @param profile_id column name of the profile IDs (default is profile_id).
#' @param dates column name of the dates (default is dates).
#' @param intensity column name of the intensity (default is intensity).
#' @param max_int selects profiles with at least one detection higher than x (default is 7).
#' @param cum_int selects x percent of profiles with highest cumulative intensity (default is 5).
#' @param freq_det selects profiles detected in more than x percent of the samples (default is 70).
#'
#' @export
#'
#' @return List with the following elements
#'
#' * min_one_criteria: a tibble containing profile IDs, dates and intensities of
#' prioritised profiles that fulfill at least one criteria.
#'
#' * all_criteria: a tibble containing profile IDs, dates and intensities of
#' prioritised profiles that fulfill all criteria.
#'
#' * rank: a tibble containing profile IDs, a logical column for each criteria
#' indicating if it is fulfilled for the given profile and priority rank
#' for profiles that fulfill all criteria.
#'
#' @details Priority rank: profiles that meet all three criteria are selected.
#' In a first step, profiles are ranked by importance for each of the three
#' criteria (e.g., highest maximum contribution = highest importance). The
#' summed-up ranks are then used to assign an overall priority rank to each
#' profile.
#'
#' @references Chonova, T., Honti, M., Loos, M., Ruppe, S., Langlois, I.,
#' Griesshaber, D.,Fenner, K., Singer, H. Unveiling Industrial Contamination in
#' the Rhine: Insights from Data Mining of High-Frequency Measurements (in prep)
#'
#'
prioritise_profiles <- function(x, profile_id = profile_id, dates = dates, intensity = intensity,
                                max_int = 7, freq_det = 70, cum_int = 5) {
  library(dplyr); library(forcats)

  x <- x %>% select(dates = {{dates}}, profile_id = {{profile_id}}, intensity = {{intensity}})

  sampled_days <- x %>% distinct(dates) %>% nrow()

  priority_profiles <- list()

  priority_profiles$rank <- x %>%
    filter(!is.na(intensity)) %>%
    group_by(profile_id) %>%
    summarise(int_max = max(round(log10(intensity), 2)),
              det_freq = n(),
              int_cum = sum(intensity)) %>%
    ungroup() %>%
    mutate(det_freq_perc = round(det_freq*100/sampled_days));

  priority_profiles$rank <- priority_profiles$rank %>%
    arrange(desc(int_cum)) %>%
    mutate(int_cum_perc = (row_number()/nrow(priority_profiles$rank))*100) %>%
    mutate() %>%
    select(-int_cum, -det_freq)

  priority_rank <- priority_profiles$rank %>%
    filter(int_max >= max_int, det_freq_perc >= freq_det, int_cum_perc <= cum_int) %>%
    mutate(Rank1 = as.numeric(fct_rev(as.factor(int_cum_perc))),
           Rank2 = as.numeric(fct_rev(as.factor(det_freq_perc))),
           Rank3 = as.numeric(fct_rev(as.factor(int_max)))) %>%
    rowwise() %>%
    mutate(Rank_sum = sum(Rank1, Rank2, Rank3, na.rm = T)) %>%
    ungroup() %>%
    arrange(Rank_sum) %>%
    mutate(ranks = as.numeric(as.factor(Rank_sum))) %>%
    mutate(ranks = row_number()) %>%
    select(profile_id, ranks)

  priority_profiles$all_criteria <- suppressMessages({
    priority_profiles$rank %>%
      filter(int_max >= max_int, det_freq_perc >= freq_det, int_cum_perc <= cum_int) %>%
      distinct(profile_id) %>%
      left_join(x) %>% select(profile_id, dates, intensity)
  })

  priority_profiles$min_one_criteria <- suppressMessages({
    priority_profiles$rank %>%
      filter(int_max >= max_int | det_freq_perc >= freq_det | int_cum_perc <= cum_int) %>%
      distinct(profile_id) %>%
      left_join(x) %>% select(profile_id, dates, intensity)
  })

  suppressMessages({
    priority_profiles$rank <- priority_profiles$rank %>%
      filter(int_max >= max_int | det_freq_perc >= freq_det | int_cum_perc <= cum_int) %>%
      left_join(priority_rank) %>%
      mutate(maximum_intensity = ifelse(int_max >= max_int, TRUE, FALSE),
             detection_frequency = ifelse(det_freq_perc >= freq_det, TRUE, FALSE),
             cumulative_intensity = ifelse(int_cum_perc <= cum_int, TRUE, FALSE)) %>%
      select(profile_id, maximum_intensity, detection_frequency,
             cumulative_intensity, ranks) %>%
      arrange(ranks)
  })

  return(priority_profiles)
}


### Intensity spread ####

#' Compute intensity spread
#'
#' This function takes a dataframe with profile IDs (numeric or character) and
#' intensity (numeric) and computes intensity spread of time profiles as described
#' in Anliker et al. 2020.
#'
#' @param x a tibble in long format.
#' @param profile_id column name of the profile IDs (default is profile_id).
#' @param intensity column name of the intensity (default is intensity).
#'
#' @export
#'
#' @return a tibble with profile IDs and their corresponding
#' intensity spread.
#'
#' @references Anliker, S., Loos, M., Comte, R., Ruff, M., Fenner, K.,
#' Singer, H., 2020. Assessing Emissions from Pharmaceutical Manufacturing
#' Based on Temporal High-Resolution Mass Spectrometry Data. Environ. Sci.
#' Technol. 54, 4110–4120. https://doi.org/10.1021/acs.est.9b07085
#'
#'
compute_int_spr <- function(x, profile_id = profile_id, intensity = intensity) {

  x <- x %>% select(profile_id = {{profile_id}}, intensity = {{intensity}})

  library(dplyr);
  int_spread <- x %>%
    group_by(profile_id) %>%
    summarise(min = min(intensity),
              max_intensity = max(intensity),
              q05 = quantile(intensity, probs = 0.05, na.rm = TRUE),
              q95 = quantile(intensity, probs = 0.95, na.rm = TRUE),
              intensity_spread = q95/q05) %>%
    ungroup() %>%
    select(profile_id, intensity_spread)
  return(int_spread)
}


### Detect breaks ####

#' Detect breaks in time series
#'
#' This function takes a dataframe with profile IDs (numeric or character),
#' dates (date format) and intensity (numeric) and detects activity breaks in
#' time profiles.
#'
#' @param x a tibble in long format.
#' @param profile_id column name of the profile IDs (default is profile_id).
#' @param dates column name of the dates (default is dates).
#' @param intensity column name of the intensity (default is intensity).
#'
#' @export
#'
#' @return A two-columns tibble with profile IDs and
#' the corresponding longest break detected.
#'
#' @references Chonova, T., Honti, M., Loos, M., Ruppe, S., Langlois, I., Griesshaber, D.,
#' Fenner, K., Singer, H. Unveiling Industrial Contamination in the Rhine:
#' Insights from Data Mining of High-Frequency Measurements (in prep)
#'
#'
detect_breaks <- function(x, profile_id = profile_id, dates = dates,
                          intensity = intensity) {
  library(dplyr);

  x <- x %>% select(dates = {{dates}}, profile_id = {{profile_id}}, intensity = {{intensity}})


  count_length <- function(x){
    r <- rle(x)
    res <- unlist(sapply(r$lengths, function(x) rep(x, x)))
  }

  detect_br <- x %>%
    filter(!is.na(intensity)) %>%
    arrange(profile_id, as.Date(dates)) %>%
    group_by(profile_id) %>%
    mutate(limit_break = min(intensity)+
             diff(quantile(intensity, probs=c(0.001,0.999),
                           na.rm = TRUE))/30, # range
           below_limit_break = ifelse(intensity < limit_break, 1, NA),
           break_length = count_length(below_limit_break)) %>%
    distinct(profile_id, break_length) %>%
    slice_max(break_length) %>%
    ungroup()

  return(detect_br)
}


### Plot time profile ####

#' Plot time profile
#'
#' This function takes a dataframe with profile IDs (numeric or character),
#' dates (date format) and intensity (numeric) and plots a selected time profile.
#'
#' @param x a tibble in long format.
#' @param profile_id column name of the profile IDs (default is profile_id).
#' @param dates column name of the dates (default is dates).
#' @param intensity column name of the intensity (default is intensity).
#' @param sel_prof the ID of the selected profile.
#'
#' @export
#'
#' @return A line plot of the selected time profile
#'
#'
plot_profile <- function(x, profile_id = profile_id, dates = dates, intensity = intensity, sel_prof) {

  library(dplyr); library(tidyr); library(ggplot2);

  x <- x %>% select(dates = {{dates}}, profile_id = {{profile_id}}, intensity = {{intensity}})

  fake_prof <- seq(min(x$dates), max(x$dates), by="days") %>%
    enframe(value = "dates") %>%
    select(dates) %>%
    mutate(profile_id = "fake")

  x %>%
    filter(profile_id == sel_prof) %>%
    full_join(fake_prof) %>%
    expand(dates, profile_id) %>% left_join(x) %>%
    filter(profile_id != "fake") %>%
    ggplot(aes(x = as.Date(dates), y = intensity)) +
    geom_line() + geom_point() +
    labs(title = sel_prof, x = "Date", y = "intensity")+
    theme_bw()
}


### Calculate distances ####

#' Calculate distance measures
#'
#' This function takes a dataframe with profile IDs (numeric or character),
#' dates (date format) and intensity (numeric) and calculates distance measure
#' on the scaled data.
#'
#'
#' @param x a tibble in long format.
#' @param x a tibble in long format.
#' @param profile_id column name of the profile IDs (default is profile_id).
#' @param dates column name of the dates (default is dates).
#' @param intensity column name of the intensity (default is intensity).
#' @param method the distance measure to be used. This must be one of "euclidean",
#' "maximum", "manhattan", "canberra", "binary" or "minkowski" (default is "euclidean").
#'
#' @export
#'
#' @return Distance measure calculated on the scaled data.
#'
calculate_dist <- function(x, profile_id = profile_id, dates = dates,
                         intensity = intensity, method = "euclidean") {

  library(dplyr); library(tidyr);

  x <- x %>% select(dates = {{dates}}, profile_id = {{profile_id}}, intensity = {{intensity}})

  x1 <- x %>%
    arrange(as.Date(dates)) %>%
    pivot_wider(id_cols = dates, names_from = profile_id,
                values_from = intensity, values_fill = 0)
  x1 <- as.data.frame(x1)
  rownames(x1) <- x1[,1 , drop = TRUE] ; x1 <- x1[,-1]
  x1 <- x1 %>% as.matrix()
  dt_sc <- scale(x1, center = FALSE, scale = apply(x1, 2, sd, na.rm = TRUE))
  dist_mat <- dist(t(dt_sc), method = method)

  return(dist_mat)

}


### Show similar profiles ####

#' Show profiles with similar time patterns
#'
#' This function takes a dataframe with profile IDs (numeric or character),
#' dates (date format) and intensity (numeric) as well as its distance measure,
#' selects profiles with similar time pattern and returns line plots.
#'
#' @param x a tibble in long format.
#' @param dist_mat a distance measure (output from find_friends() or dist() function).
#' @param profile_id column name of the profile IDs (default is profile_id).
#' @param dates column name of the dates (default is dates).
#' @param intensity column name of the intensity (default is intensity).
#' @param sel_prof the ID of the selected profile.
#' @param nr the number of similar profiles to plot (default is 3).
#'
#' @export
#'
#' @return List with the following elements
#'
#' List with the following elements
#'
#' * plt: line plots of the selected time profile and profiles with time patterns similar to it.
#'
#' * prof: a table with IDs of similar profiles and their respective distances to the selected profile.
#'
#'
show_friends <- function(x, dist_mat, profile_id = profile_id, dates = dates,
                         intensity = intensity, sel_prof, nr = 3) {

  library(dplyr); library(tibble); library(ggplot2)

  x <- x %>% select(dates = {{dates}}, profile_id = {{profile_id}}, intensity = {{intensity}})

  prof <- as.matrix(dist_mat)
  prof <- enframe(prof[sel_prof,]) %>% slice_min(value, n = nr+1) %>% select(prof_id = name, distance = value);

  fake_prof <- seq(min(x$dates), max(x$dates), by="days") %>%
    enframe(value = "dates") %>%
    select(dates) %>%
    mutate(profile_id = "fake")

   p1 <- x %>%
    filter(profile_id %in% prof$prof_id) %>%
    full_join(fake_prof) %>%
    expand(dates, profile_id) %>% left_join(x) %>%
    filter(profile_id != "fake") %>%
    ggplot(aes(x = as.Date(dates), y = intensity)) +
    geom_line(size = 0.5) + geom_point(size = 0.5) +
    facet_wrap(vars(factor(profile_id, levels = prof$prof_id)), scales = "free", ncol = 1) +
    labs(x = "Date", y = "") +
    theme_bw()

  theme_cust <- theme(legend.position="bottom", legend.direction = "vertical",
                      legend.key.height = unit(0.3, "cm"), # reduce place between legend elements
                      legend.box.margin = margin(t = -20, r = 0, b = 0, l = 0), # plot legend closer to graph
                      legend.background = element_rect(fill = "transparent", color = NA), # remove legend box background
                      legend.title = element_blank(), legend.text = element_text(size = 8))

  p2 <- x %>%
    filter(profile_id %in% prof$prof_id) %>%
    full_join(fake_prof) %>%
    expand(dates, profile_id) %>% left_join(x) %>%
    filter(profile_id != "fake") %>%
    ggplot(aes(x = as.Date(dates), y = intensity, color = profile_id), si) +
    geom_line(size = 0.5) + geom_point(size =0.5) +
    labs(title = sel_prof, x = "Date", y = "intensity") +
    theme_bw() + theme_cust


  p3 <- x %>%
    filter(profile_id %in% prof$prof_id) %>%
    full_join(fake_prof) %>%
    expand(dates, profile_id) %>% left_join(x) %>%
    filter(profile_id != "fake") %>%
    ggplot(aes(x = as.Date(dates), y = log10(intensity), color = profile_id), si) +
    geom_line(size = 0.5) + geom_point(size =0.5) +
    labs(title = sel_prof, x = "Date", y = "log10(intensity)") +
    theme_bw() + theme_cust

  plt <- ggpubr::ggarrange(p1, ggpubr::ggarrange(p2, p3, ncol = 1, nrow = 2, labels = c("B", "C")),
                           ncol = 2, nrow = 1, labels = c("A"))

  res <- list(plt = plt, prof = prof)

  return(res)

}
