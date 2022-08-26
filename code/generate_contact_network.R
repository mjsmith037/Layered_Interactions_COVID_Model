library(magrittr)
# library(tidycensus) #this is better for the future
library(tidygraph)
library(tidyverse)
library(assertthat)


ensure_rates_are_functions <- function(rate) {
  if (is.numeric(rate) && length(rate) == 1) {
    rate <- eval(substitute(function(n) {rep(rr, n)}, list(rr=rate)))
  }
  assert_that(is.function(rate))
  return(rate)
}

intelligently_assign_ages <- function(node_df, locale="United States") {
  # TODO: more likely to be children in larger households (?)
  age_distribution <- read_csv("../data/age_demographics.csv", skip=1, col_types="dddddddddddddc") %>%
    filter(Locale == locale) %>%
    select(-Locale) %>%
    pivot_longer(everything(), names_to="age", values_to="proportion")
  ## ensure no households of just young children
  node_df %>%
    mutate(age_class = sample(age_distribution$age, n(), replace=TRUE,
                              prob=age_distribution$proportion)) %>%
    mutate(is_child = age_class %in% c("Under 5 years", "5 to 9 years", "10 to 14 years")) %>%
    group_nest(household_id) %>%
    mutate(size = map_int(data, nrow),
           children = map_int(data, ~sum(.$is_child))) %>%
    mutate(too_many_children = size == children,
           extra_adults = size > children + 2) %>%
    unnest(data) %>%
    mutate(age_class = replace(
      age_class,
      c(node_id[too_many_children], node_id[extra_adults & !is_child]),
      age_class[c(node_id[extra_adults & !is_child], node_id[too_many_children])]
    )) %>%
    ## shouldn't be necessary, but just to make sure nodes align upon returning
    arrange(node_id) %>%
    pull(age_class)
}

# generate the household units (currently fully connected within households)
# note: this allows for further fine-tuning of household structure
generate_household_unit <- function(household_size) {
  play_erdos_renyi(household_size, 1, directed=FALSE)
}

# add coworker links
add_employment_layer <- function(network, coworker_transmission_rate, risky_behavior_level=2) {
  # remove links previously assigned to this route of transmission
  network <- network %E>% filter(type != "coworker")

  ## mostly arbitrary
  ## unverified source: https://www.naics.com/business-lists/counts-by-company-size/
  ## Additional caveats: omitting case of <5 employees
  ##                     omitting cases of >100 employees (for now, bc small populations and bc does not reflect typical contacts)
  ##                     stochastic asignment means actual sizes could differ substantially from specified size (from lower)
  workplace_sizes <- tibble(lower=c(5, 10, 20, 50),
                            upper=c(9, 19, 49, 99),
                            business_count=c(1803143, 779088, 471372, 169239)) %>%
    mutate(employee_count=lower * business_count,
           proportion=employee_count / sum(employee_count))
  workplace_assignment <- network %N>% as_tibble() %>%
    group_by(household_id) %>% mutate(household_vuln = any(vulnerable)) %>% ungroup() %>%
    # extract working age adults (US retirement age is 67)
    filter(age_class %in% c("20 to 24 years",
                            "25 to 34 years",
                            "35 to 44 years",
                            "45 to 54 years",
                            "55 to 59 years",
                            "60 to 64 years"),
           case_when(risky_behavior_level == 0 ~ !household_vuln,
                     risky_behavior_level == 1 ~ !vulnerable,
                     risky_behavior_level == 2 ~ TRUE)) %>%
    # account for unemployment
    # MN unemployment rate (May 2020): 9.4% (https://mn.gov/deed/data/current-econ-highlights/state-national-employment.jsp)
    filter(rbernoulli(n(), 0.9)) %>%
    # assing to workplace of a particular size
    mutate(workplace_size = workplace_sizes %>%
    {sample(.$lower, n(), replace=TRUE, prob=.$proportion)}) %>%
    # asign individuals to workplaces
    group_by(workplace_size) %>%
    mutate(workplace_id = str_c(workplace_size,
                                sample(1:(n() %/% unique(workplace_size) + 1), n(), replace=TRUE),
                                sep="_")) %>%
    ungroup()
  coworker_transmission <- network %N>% as_tibble() %>%
    # add workplace id to node info
    left_join(workplace_assignment %>% select(node_id, workplace_id), by = "node_id") %>%
    # add links (fully connected within workspaces)
    tbl_graph(group_by(workplace_assignment, workplace_id) %>%
                group_modify(~tidyr::crossing(from=.x$node_id, to=.x$node_id) %>% filter(from < to)) %>%
                ungroup() %>%
                transmute(from, to, weight=coworker_transmission_rate(n()), type="coworker"))
  return(coworker_transmission)
}

# add school layer
add_school_layer <- function(network, classmate_transmission_rate, approx_classroom_size, risky_behavior_level=2) {
  # remove links previously assigned to this route of transmission
  network <- network %E>% filter(type != "classmate")

  # until other layers are added as well, we can just randomly assign groups
  # TODO: consider having stronger and weaker ties for classrooms, schools respectively
  # TODO: include "under 5 years" as daycare?
  school_assignment <- network %N>% as_tibble() %>%
    group_by(household_id) %>% mutate(household_vuln = any(vulnerable)) %>% ungroup() %>%
    filter(age_class %in% c("5 to 9 years", "10 to 14 years", "15 to 19 years"),
           case_when(risky_behavior_level == 0 ~ !household_vuln,
                     risky_behavior_level == 1 ~ !vulnerable,
                     risky_behavior_level == 2 ~ TRUE)) %>%
    group_by(age_class) %>%
    mutate(school_count = n() %/% approx_classroom_size + 1) %>%
    group_by(age_class, household_id) %>%
    mutate(school_id = str_c(age_class, sample(1:unique(school_count), 1), sep="_")) %>%
    ungroup()
  classmate_transmission <- network %N>% as_tibble() %>%
    # add school id to node info
    left_join(school_assignment %>% select(node_id, school_id), by = "node_id") %>%
    # add links (fully connected within classrooms)
    tbl_graph(group_by(school_assignment, school_id) %>%
                group_modify(~tidyr::crossing(from=.x$node_id, to=.x$node_id) %>% filter(from < to)) %>%
                ungroup() %>%
                transmute(from, to, weight=classmate_transmission_rate(n()), type="classmate"))
  return(classmate_transmission)
}

# formation of socially monogomous household units
add_household_merging_layer <- function(network, between_household_transmission_rate, number_households_merged) {
  # remove links previously assigned to this route of transmission
  network <- network %E>% filter(type != "between_household")

  # TODO: merge by physical proximity?
  between_household_transmission <- network %N>%
    as_tibble() %>%
    nest(tmp=-household_id) %>%
    mutate(household_merger_id = sample(ceiling(1:max(household_id) / number_households_merged))) %>%
    unnest(cols=tmp) %>%
    tbl_graph(group_by(., household_merger_id) %>%
                group_modify(~tidyr::crossing(from=.x$node_id, to=.x$node_id) %>% filter(from < to)) %>%
                ungroup() %>%
                transmute(from, to, weight=between_household_transmission_rate(n()), type="between_household"))
  return(between_household_transmission)
}

add_bar_layer <- function(network, bar_transmission_rate, approx_bar_size, percent_bar_goers,
                          risky_behavior_level=1) {
  # risky behavior level:
  #   0: no one from households with vulnerable members goes to bar
  #   1: no vulnerable people go to bars
  #   2: vulnerable do not act differently

  # remove links previously assigned to this route of transmission
  network <- network %E>% filter(type != "bar")

  bar_assignment <- network %N>% as_tibble() %>%
    group_by(household_id) %>% mutate(household_vuln = any(vulnerable)) %>% ungroup() %>%
    # no minors in the bars
    filter(age_class %>%
             is_in(c("Under 5 years", "5 to 9 years", "10 to 14 years", "15 to 19 years")) %>%
             not(),
           # implement vulnerable behavior re. bar-going
           case_when(risky_behavior_level == 0 ~ !household_vuln,
                     risky_behavior_level == 1 ~ !vulnerable,
                     risky_behavior_level == 2 ~ TRUE)) %>%
    sample_frac(percent_bar_goers) %>%
    mutate(bar_id = sample(1:(n() %/% approx_bar_size + 1), n(), replace=TRUE)) %>%
    ungroup()
  bar_transmission <- network %N>% as_tibble() %>%
    # add school id to node info
    left_join(bar_assignment %>% select(node_id, bar_id), by = "node_id") %>%
    # add links (fully connected within classrooms)
    tbl_graph(group_by(bar_assignment, bar_id) %>%
                group_modify(~tidyr::crossing(from=.x$node_id, to=.x$node_id) %>% filter(from < to)) %>%
                ungroup() %>%
                transmute(from, to, weight=bar_transmission_rate(n()), type="bar"))
  return(bar_transmission)
}

generate_contact_network <- function(number_of_households=100,
                                     # core layer: household transmission
                                     within_household_transmission_rate=function(n) rep(0.4, n),
                                     # socially monogomous household units
                                     between_household_transmission_rate=NULL, #function(n) rep(0.3, n),
                                     number_households_merged=2,
                                     # school transmission
                                     classmate_transmission_rate=NULL, #function(n) rep(0.2, n),
                                     approx_classroom_size=20,
                                     # work layer
                                     coworker_transmission_rate=NULL, #function(n) rep(0.1, n),
                                     # bar layer
                                     bar_transmission_rate=NULL, #function(n) rep(0.3, n),
                                     approx_bar_size=20,
                                     percent_bar_goers=0.1,
                                     # additional params
                                     risky_behavior_level=2,
                                     locale="United States",
                                     remove_unused_layers=TRUE) {
  # get a distribution of household sizes
  household_size_distribution <-
    read_csv("../data/household_size_demographics.csv",
             skip=3, col_types="cddddddddl",
             col_names=c("Locale", "1", "2", "3", "4", "5", "6",
                         "7", # note actual heading here is 7 **or more**
                         "average_household_size", "extrapolated")) %>%
    filter(Locale == locale) %>%
    select(matches("\\d")) %>%
    summarise(across(matches("\\d"),
                     # note that this rounding can introduce slight deviations from
                     # the prescribed number of households
                     ~divide_by(., sum(c_across())) %>%
                       multiply_by(number_of_households) %>%
                       round())) %>%
    as.vector() %>%
    rep(as.integer(names(.)), times=.)

  # generate the household units and combine
  household_network <- lapply(household_size_distribution, generate_household_unit) %>%
    do.call(bind_graphs, .) %>%
    mutate(household_id = rep(1:length(household_size_distribution), household_size_distribution),
           node_id = 1:n()) %>%
    activate(edges) %>% mutate(weight = within_household_transmission_rate(n()), type = "household") %>%
    # since age is likely important, first give all individuals age (demographic data)
    activate(nodes) %>% mutate(age_class = intelligently_assign_ages(.N(), locale=locale))

  # assign vulnerability based on co-morbidities
  # MN obesity rate: 13% for children, 30% for adults (https://data.web.health.state.mn.us/web/mndata/obesity_basic)
  # or just use hospitalizations (or icu admission or case-fatality) directly
  # case fatality rates by age (http://www.ecie.com.ar/images/paginas/COVID-19/4MMWR-Severe_Outcomes_Among_Patients_with_Coronavirus_Disease_2019_COVID-19-United_States_February_12-March_16_2020.pdf)
  case_fatality_distribution <- read_csv("../data/case_fatality.csv", skip=1, col_types="cccc") %>%
    mutate_at(vars(-age), ~str_split_fixed(., "-", n=2) %>%
                apply(1, . %>% as.numeric() %>% divide_by(100) %>% mean(na.rm=TRUE)))

  household_network %<>%
    activate(nodes) %>%
    mutate(vulnerable = rbernoulli(
      n(), case_when(age_class %in% c("Under 5 years", "5 to 9 years", "10 to 14 years", "15 to 19 years") ~
                       case_fatality_distribution %>% filter(age == "0-19") %>% pull(hospitalization),
                     age_class %in% c("20 to 24 years", "25 to 34 years", "35 to 44 years") ~
                       case_fatality_distribution %>% filter(age == "20-44") %>% pull(hospitalization),
                     age_class == "45 to 54 years" ~
                       case_fatality_distribution %>% filter(age == "45-54") %>% pull(hospitalization),
                     age_class %in% c("55 to 59 years", "60 to 64 years") ~
                       case_fatality_distribution %>% filter(age == "55-64") %>% pull(hospitalization),
                     age_class == "65 to 74 years" ~
                       case_fatality_distribution %>% filter(age == "65-74") %>% pull(hospitalization),
                     age_class == "75 to 84 years" ~
                       case_fatality_distribution %>% filter(age == "75-84") %>% pull(hospitalization),
                     age_class == "85 years and over" ~
                       case_fatality_distribution %>% filter(age == "85+") %>% pull(hospitalization))))

  # there are several ways we can model an underlying baseline risk of infection:
  # 1) we can add a fully connected network layer (with relatively small weights)
  # background_transmission <- play_erdos_renyi(sum(household_size_distribution), 1, directed=FALSE) %>%
  #   mutate(node_id = 1:n()) %>%
  #   activate(edges) %>% mutate(weight=background_transmission_rate(n()), type = "background")
  # 2) we can make the network explicitly spatial and use a distance kernel

  # add employment layer
  # possible implementation: employment rate * job_distribution -> number of employees in each job
  # select individuals in a given age range to fill jobs randomly
  if (!is.null(coworker_transmission_rate)) {
    coworker_transmission <- add_employment_layer(household_network,
                                                  ensure_rates_are_functions(coworker_transmission_rate),
                                                  risky_behavior_level)
  }

  # add school layer
  if (!is.null(classmate_transmission_rate)) {
    classmate_transmission <- add_school_layer(household_network,
                                               ensure_rates_are_functions(classmate_transmission_rate),
                                               approx_classroom_size,
                                               risky_behavior_level)
  }

  if (!is.null(between_household_transmission_rate) && (number_households_merged > 1)) {
    between_household_transmission <- add_household_merging_layer(household_network,
                                                                  ensure_rates_are_functions(between_household_transmission_rate),
                                                                  number_households_merged,
                                                                  risky_behavior_level)
  }

  if (!is.null(bar_transmission_rate)) {
    bar_transmission <- add_bar_layer(household_network,
                                      ensure_rates_are_functions(bar_transmission_rate),
                                      approx_bar_size,
                                      percent_bar_goers,
                                      risky_behavior_level)
  }

  #### join the layers ####
  layers <- ls() %>% .[str_detect(., "_transmission$")]
  full_contact_network <- household_network
  tmp <- lapply(layers, function(i) {
    full_contact_network <<-
      graph_join(full_contact_network, get(i),
                 by=intersect(colnames(get(i) %N>% as_tibble()),
                              colnames(full_contact_network %N>% as_tibble())))
  })
  full_contact_network %<>% to_undirected()

  return(full_contact_network)
}
