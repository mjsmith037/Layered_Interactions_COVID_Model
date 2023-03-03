library(tidyverse)

# source("common_functions.R")

#### Summarizing Functions ####
# Total number infected
total_infected <- . %>%
  filter(time == max(time), status != "Susceptible") %>%
  pull(n) %>%
  sum()

# Total number of deaths
total_died <- . %>%
  filter(time == max(time), status == "Dead") %>%
  pull(n) %>%
  sum()

# Average number infectious
average_infectious <- . %>%
  filter(status == "Asymptomatically Infectious" | status == "Symptomatically Infectious") %>%
  group_by(time) %>%
  summarise(n = sum(n), .groups="drop") %>%
  pull(n) %>%
  mean()

# Time to peak
time_to_peak <- . %>%
  filter(status == "Asymptomatically Infectious" | status == "Symptomatically Infectious") %>%
  group_by(time) %>%
  summarise(n = sum(n), .groups="drop") %>%
  slice_max(n, with_ties=FALSE) %>%
  pull(time)

# Time to peak death (proxy for hospitalization need)
time_to_peak_death <- . %>%
  filter(status == "Dead") %>%
  mutate(change = n - lag(n)) %>% 
  slice_max(change, with_ties=FALSE) %>%
  pull(time) %>%
  {ifelse(length(.) == 0, Inf, .)}
# Maximum number infectious
peak_infectious <- . %>%
  filter(status == "Asymptomatically Infectious" | status == "Symptomatically Infectious") %>%
  group_by(time) %>%
  summarise(n = sum(n), .groups="drop") %>%
  pull(n) %>%
  max()


#### Get all summaries ####
summarize_epidemic <- . %>%
  clean_output() %>%
  summarise(!!!list("Total Infected"     = total_infected(.),
                    "Average Infectious" = average_infectious(.),
                    "Total Died"         = total_died(.),
                    "Time to Peak"       = time_to_peak(.),
                    "Time to Peak Death" = time_to_peak_death(.),
                    "Maximum Infectious" = peak_infectious(.)))


# total number of vulnerable people infected
vulnerable_infected <- function(output, network) {
  clean_output(output, grouping_vars=c("status", "vulnerable"),
               vulnerable=network %N>% as_tibble() %>% pull(vulnerable)) %>%
    filter(time == max(time), vulnerable == "Vulnerable", status != "Susceptible") %>%
    pull(n) %>%
    sum()
}

# time to first vulnerable person infected
time_to_vulnerable <- function(output, network) {
  clean_output(output, grouping_vars=c("status", "vulnerable"),
               vulnerable=network %N>% as_tibble() %>% pull(vulnerable)) %>%
    filter(vulnerable == "Vulnerable", status == "Exposed", n > 0) %>%
    pull(time) %>%
    {ifelse(length(.) == 0, Inf, min(.))}
}

# number of vulnerable households infected
vuln_households_inf <- function(output, network) {
  tibble(status = output %>%
           set_colnames(1:ncol(.)) %>%
           as_tibble() %>%
           mutate(across(everything(), . %>% factor(levels=0:5, labels=c("S", "E", "Ia", "Is", "R", "D")))) %>%
           pull(),
         house_id = network %N>%
           as_tibble() %>%
           pull(household_id),
         house_vuln = network %N>%
           as_tibble() %>%
           group_by(household_id) %>%
           mutate(any(vulnerable)) %>%
           pull()) %>%
    group_by(house_id) %>% 
    summarise(house_vuln = unique(house_vuln),
              house_inf = ifelse(all(status == "S"), "safe", "infected") %>%
                factor(levels=c("safe", "infected")), .groups="drop") %>%
    filter(house_vuln) %>%
    tabyl(house_inf) %>%
    filter(house_inf == "infected") %>%
    pull(percent)
}
