source("generate_contact_network.R")
source("quantify_disease_outcomes.R")
library(janitor)
library(Rcpp)
sourceCpp("simulate_disease_on_network.cpp")

# intended to be run from the command line, but replace cArgs[x] with default
# values to run interactively
cArgs <- commandArgs(TRUE)

locale <- cArgs[1]                 # "Texas"
background <- as.numeric(cArgs[2]) # 0.001
nsims <- as.integer(cArgs[3])      # 10000
pop_size <- as.integer(cArgs[4])   # 3000
pseudo <- as.logical(cArgs[5])     # FALSE
rnaught <- as.logical(cArgs[6])    # FALSE
set.seed(as.integer(cArgs[7]))     # 0

# a function to clean up the C++ simulation output
clean_output <- function(output, grouping_vars="status", vulnerable=NULL) {
  # convert the disease simulation output to a more manageable format
  timeseries <- output %>%
    set_colnames(1:ncol(.)) %>%
    as_tibble() %>%
    mutate(across(everything(), . %>% factor(levels=0:5, labels=c("S", "E", "Ia", "Is", "R", "D"))),
           vulnerable = vulnerable) %>%
    # wide -> long format
    pivot_longer(names_to="time", values_to="status", matches("\\d+")) %>%
    # clean up types
    mutate(time = as.integer(time),
           status = factor(status)) %>%
    # summarise number of each class at each timestep
    count(!!!syms(grouping_vars), time) %>%
    # add 0's for missing classes at each timestep
    complete(!!!syms(grouping_vars), nesting(time), fill=list(n=0)) %>%
    # clean up vulnerable flag for plotting (if necessary)
    mutate(across(any_of("vulnerable"), ~case_when(. ~ "Vulnerable", TRUE ~ "Not Vulnerable"))) %>%
    # clean up statuses
    mutate(status = factor(status, levels=c("S", "E", "Ia", "Is", "R", "D"),
                           labels=c("Susceptible", "Exposed", "Asymptomatically Infectious",
                                    "Symptomatically Infectious", "Recovered", "Dead")))
  # return the cleaned output
  return(timeseries)
}

# workhorse simulation function adjusts parameters and calls C++ simulation code
run_simulations <- function(rseed=0,
                            number_of_households=1000,
                            maxtime=1000,
                            within_household_transmission_rate=0.13,
                            background_transmission_rate=0.001,
                            classmate_transmission_rate=0.01,
                            coworker_transmission_rate=0.01,
                            risky_behavior_level=2,
                            locale="United States",
                            pseudo_risk=FALSE,
                            rnaught=FALSE) {
  #### Run the Simulation ####
  set.seed(rseed)
  ## first generate a contact network
  full_contact_network <-
    generate_contact_network(number_of_households=number_of_households,
                             within_household_transmission_rate=function(n) rep(within_household_transmission_rate, n),
                             classmate_transmission_rate=function(n) rep(classmate_transmission_rate, n),
                             coworker_transmission_rate=function(n) rep(coworker_transmission_rate, n),
                             risky_behavior_level=risky_behavior_level,
                             locale=locale)

  # pick a patient 0
  initial_infected_node <- full_contact_network %N>%
    filter(centrality_degree() > 0) %>%
    as_tibble() %>%
    pull(node_id) %>%
    sample(1)
  initial_condition <- rep(0, full_contact_network %>% with_graph(graph_order()))
  initial_condition[initial_infected_node] <- 2

  # optionally down-sample links to better compare to risk-tolerance regimes
  if (pseudo_risk) {
    downsample <- function(data, n) {
      if (is.na(n)) {
        return(tibble(from=integer(), to=integer(), weight=numeric()))
      } else {
        return(slice_sample(data, n=n))
      }
    }
    set.seed(rseed)
    subsampled_edges <-
      generate_contact_network(number_of_households=number_of_households,
                               within_household_transmission_rate=function(n) rep(within_household_transmission_rate, n),
                               classmate_transmission_rate=function(n) rep(classmate_transmission_rate, n),
                               coworker_transmission_rate=function(n) rep(coworker_transmission_rate, n),
                               risky_behavior_level=2,
                               locale=locale) %E>%
      as_tibble() %>%
      nest(data=!type) %>%
      left_join(full_contact_network %E>% as_tibble() %>% tabyl(type), by="type") %>%
      rowwise() %>%
      mutate(n = min(n, nrow(data))) %>%
      ungroup() %>%
      mutate(data = map2(data, n, downsample)) %>%
      select(!c(n, percent)) %>%
      unnest(data)
    full_contact_network <- tbl_graph(nodes=full_contact_network %N>% as_tibble(), edges=subsampled_edges)
  }

  # run the simulation
  simulation_output <- run_disease_simulation_on_network(
    maxtime,
    full_contact_network %E>% select(from, to) %>% as_tibble() %>% as.matrix() %>% subtract(1),
    initial_condition,
    full_contact_network %N>% as_tibble() %>% pull(vulnerable),
    full_contact_network %E>% as_tibble() %>% pull(weight),
    background_transmission_rate, ifelse(rnaught, 0, 1/5.5), 0.35, 1/4.5, 1/27000, 1/1000
  )

  # clean up the output and summarize the epidemic using several measures
  summarize_epidemic(simulation_output) %>%
    mutate(`Vulnerable Infected` = vulnerable_infected(simulation_output, full_contact_network),
           `Time to Vulnerable Infection` = time_to_vulnerable(simulation_output, full_contact_network),
           `Percent Vulnerable Households Infected` = vuln_households_inf(simulation_output, full_contact_network)) %>%
    mutate(pop_size = full_contact_network %>% with_graph(graph_order()),
           vuln_pop_size = full_contact_network %N>% as_tibble() %>% pull(vulnerable) %>% sum())
}

# generate a list of parameter combinations to simulate over, note transmission
# rates are sampled uniformly in log-space
param_list <- tibble(classmate_transmission_rate  = 10^-runif(nsims, 0.5, 3),
                     coworker_transmission_rate   = 10^-runif(nsims, 0.5, 3)) %>%
  mutate(background_transmission_rate = background,
         locale=locale,
         # convert the (approximate) population size to the number of households necessary
         number_of_households=case_when(locale == "United States" ~ round(pop_size / 2.43),
                                        locale == "Florida"       ~ round(pop_size / 2.48),
                                        locale == "Texas"         ~ round(pop_size / 2.84))) %>%
  # perform for all three risk-tolerance regimes
  {bind_rows(mutate(., risky_behavior_level=0),
             mutate(., risky_behavior_level=1),
             mutate(., risky_behavior_level=2))} %>%
  mutate(pseudo_risk = pseudo,
         rnaught     = rnaught,
         rseed       = 1:n()) %>%
  rowwise() %>%
  group_split()

# allow function to be called with list of parameters
run_simulations <- lift(run_simulations)
# and loop through the parameter list; note the mc.cores parameter for parallelization
results <- parallel::mclapply(param_list, mc.cores=1, function(x) {
  run_simulations(x) %>% bind_cols(x)
}) %>% bind_rows()

write_csv(results, file=str_glue("../results/{background}_{ifelse(pseudo, 'pseudorisk', 'realrisk')}_{ifelse(rnaught, 'R0', 'full')}_{locale}.csv"))

