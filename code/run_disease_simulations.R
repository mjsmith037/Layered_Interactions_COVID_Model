source("common_functions.R")
source("common_parameters.R")
source("generate_contact_network.R")
source("quantify_disease_outcomes.R")

library(Rcpp)
sourceCpp("simulate_disease_on_network.cpp")

library(janitor)

cArgs <- commandArgs(TRUE)

locale <- cArgs[1]                 # "Texas"
background <- as.numeric(cArgs[2]) # 0.1
nsims <- as.integer(cArgs[3])      # 10000
set.seed(as.integer(cArgs[4]))     # 0
pseudo <- as.logical(cArgs[5])     # FALSE

run_simulations <- function(rseed=0,
                            number_of_households=1000,
                            maxtime=1000,
                            within_household_transmission_rate=0.13,
                            background_transmission_rate=0.001,
                            classmate_transmission_rate=0.01,
                            coworker_transmission_rate=0.01,
                            risky_behavior_level=2,
                            locale="United States",
                            pseudo_risk=FALSE) {
  #### Run the Simulation ####
  set.seed(rseed)
  full_contact_network <-
    generate_contact_network(number_of_households=number_of_households,
                             within_household_transmission_rate=function(n) rep(within_household_transmission_rate, n),
                             classmate_transmission_rate=function(n) rep(classmate_transmission_rate, n),
                             coworker_transmission_rate=function(n) rep(coworker_transmission_rate, n),
                             risky_behavior_level=risky_behavior_level,
                             locale=locale)

  initial_infected_node <- full_contact_network %N>%
    filter(centrality_degree() > 0) %>%
    as_tibble() %>%
    pull(node_id) %>%
    sample(1)
  initial_condition <- rep(0, full_contact_network %>% with_graph(graph_order()))
  initial_condition[initial_infected_node] <- 2

  # optionally down-sample links to better compare to risk-tolerance regimes
  if (pseudo_risk) {
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
      mutate(data = map2(data, n, ~slice_sample(.x, n=.y))) %>%
      select(!c(n, percent)) %>% 
      unnest(data)
    full_contact_network <- tbl_graph(nodes=full_contact_network %N>% as_tibble(), edges=subsampled_edges)
  }

  simulation_output <- run_disease_simulation_on_network(
    maxtime,
    full_contact_network %E>% select(from, to) %>% as_tibble() %>% as.matrix() %>% subtract(1),
    initial_condition,
    full_contact_network %N>% as_tibble() %>% pull(vulnerable),
    full_contact_network %E>% as_tibble() %>% pull(weight),
    background_transmission_rate, 1/5.5, 0.35, 1/4.5, 1/27000, 1/1000
  )

  # plot_network(full_contact_network)
  # plot_area(simulation_output, with_graph(full_contact_network, graph_order()))
  # plot_line(simulation_output, with_graph(full_contact_network, graph_order()))

  summarize_epidemic(simulation_output) %>%
    mutate(`Vulnerable Infected` = vulnerable_infected(simulation_output, full_contact_network),
           `Time to Vulnerable Infection` = time_to_vulnerable(simulation_output, full_contact_network)) %>%
    mutate(pop_size = full_contact_network %>% with_graph(graph_order()),
           vuln_pop_size = full_contact_network %N>% as_tibble() %>% pull(vulnerable) %>% sum())
}

param_list <- tibble(classmate_transmission_rate  = 10^-runif(nsims, 0.5, 3),
                     coworker_transmission_rate   = 10^-runif(nsims, 0.5, 3)) %>%
  mutate(background_transmission_rate = background,
         locale=locale,
         number_of_households=case_when(locale == "United States" ~ 1235,
                                        locale == "Florida"       ~ 1211,
                                        locale == "Texas"         ~ 1056,
                                        locale == "US_1960"       ~ 921)) %>%
  {bind_rows(mutate(., risky_behavior_level=0),
             mutate(., risky_behavior_level=1),
             mutate(., risky_behavior_level=2))} %>%
  mutate(pseudo_risk = pseudo,
         rseed       = 1:n()) %>%
  rowwise() %>%
  group_split()

run_simulations <- lift(run_simulations)
results <- parallel::mclapply(param_list, mc.cores=128, function(x) {
  run_simulations(x) %>% bind_cols(x)
}) %>% bind_rows()

write_csv(results, file=str_glue("../results/{background}_{ifelse(pseudo, 'PseudoRisk', 'Background_Rate')}_Results_{locale}.csv"))
# write_csv(results, file="../results/PopulationStructure_by_PseudoRiskBehavior_Results.csv")

################################################################################

# risk_results <- read_csv("../results/Differential_Behavior_Results.csv", col_types="dddddddddddd")
# risk_results %>% filter(`Maximum Infectious` > 5) %>%
#   mutate(risky_behavior_level = factor(risky_behavior_level, levels=0:2,
#                                        labels=c("Vulnerable Households Avoid Work/School",
#                                                 "Vulnerable Individuals Avoid Work/School",
#                                                 "No Difference in Behavior"))) %>%
#   pivot_longer(!c(risky_behavior_level, `Total Infected`, `Average Infectious`, `Total Died`, `Time to Peak`, `Maximum Infectious`, `Vulnerable Infected`, `Time to Vulnerable Infection`), values_to="metric") %>%
#   pivot_longer(!c(name, metric, risky_behavior_level), names_to="variable") %>%
#   ggplot() +
#   aes(x=metric, y=value, colour=factor(risky_behavior_level)) +
#   geom_point() +
#   geom_smooth(se=FALSE) +
#   # geom_boxplot(alpha=0.75, width=0.33) +
#   facet_grid(variable~name, scales="free") +
#   scale_colour_manual(values=my_cols[c(1,3,5)] %>% setNames(NULL)) +
#   # coord_trans(y="log1p") +
#   scale_y_log10() +
#   scale_x_log10() +
#   theme_bw() +
#   theme(legend.position="bottom",
#         legend.title=element_blank(),
#         axis.title.x=element_blank())
#
# risk_results %>% filter(`Maximum Infectious` > 5) %>%
#   mutate(risky_behavior_level = factor(risky_behavior_level, levels=0:2,
#                                        labels=c("Vulnerable Households Avoid Work/School",
#                                                 "Vulnerable Individuals Avoid Work/School",
#                                                 "No Difference in Behavior"))) %>%
#   pivot_longer(!c(background_transmission_rate, classmate_transmission_rate, coworker_transmission_rate, risky_behavior_level), names_to="variable") %>%
#   group_by(variable) %>%
#   group_modify(~aov(value~background_transmission_rate + classmate_transmission_rate + coworker_transmission_rate + risky_behavior_level,
#                     data=.) %>% broom::tidy()) %>%
#   transmute(variable, term, sig=ifelse(p.value < 0.05, "*", " ")) %>%
#   pivot_wider(names_from=term, values_from=sig)
#
# risk_results %>% filter(`Maximum Infectious` > 5) %>%
#   mutate(risky_behavior_level = factor(risky_behavior_level, levels=0:2,
#                                        labels=c("Vulnerable Households Avoid Work/School",
#                                                 "Vulnerable Individuals Avoid Work/School",
#                                                 "No Difference in Behavior"))) %>%
#   pivot_longer(!c(background_transmission_rate, classmate_transmission_rate, coworker_transmission_rate, risky_behavior_level),
#                names_to="variable") %>%
#   ggplot() +
#   aes(x=classmate_transmission_rate, y=coworker_transmission_rate, colour=log(value)) +
#   geom_point() +
#   facet_grid(risky_behavior_level~variable, scales="free") +
#   scale_x_log10() +
#   scale_y_log10() +
#   theme_bw() +
#   theme(legend.position="bottom")

#risk_results <- read_csv("../results/Differential_Behavior_Results.csv", col_types="dddddddddddd")
#risk_results %>%
#  mutate(risky_behavior_level = factor(risky_behavior_level, levels=0:2,
#                                       labels=c("Vulnerable Households Avoid Work/School",
#                                                "Vulnerable Individuals Avoid Work/School",
#                                                "No Difference in Behavior"))) %>%
#  pivot_longer(!c(risky_behavior_level, `Total Infected`, `Average Infectious`, `Total Died`, `Time to Peak`, `Maximum Infectious`, `Vulnerable Infected`, `Time to Vulnerable Infection`), values_to="metric") %>%
#  pivot_longer(!c(name, metric, risky_behavior_level), names_to="variable") %>%
#  ggplot() +
#  aes(x=metric, y=value, colour=factor(risky_behavior_level)) +
#  # ggbeeswarm::geom_quasirandom(groupOnX=TRUE, dodge.width=0.67, varwidth=TRUE) +
#  geom_point() +
#  # geom_boxplot(alpha=0.75, width=0.33) +
#  facet_grid(variable~name, scales="free") +
#  scale_colour_manual(values=my_cols[c(1,3,5)] %>% setNames(NULL)) +
#  coord_trans(y="log1p") +
#  theme_bw() +
#  theme(legend.position="bottom",
#        legend.title=element_blank(),
#        axis.title.x=element_blank())

# read_csv("../results/MDH_Bar_Risk_Results.csv") %>%
#   # risk_results %>%
#   mutate(risky_behavior_level = factor(risky_behavior_level, levels=0:2,
#                                        labels=c("Vulnerable Households Avoid Bars",
#                                                 "Vulnerable Individuals Avoid Bars",
#                                                 "No Difference in Behavior"))) %>%
#   filter(`Vulnerable Infected` > 0) %>%
#   pivot_longer(c(-percent_bar_goers, -approx_bar_size, -risky_behavior_level, -bar_transmission_rate),
#                values_to="metric") %>%
#   pivot_longer(c(-name, -metric, -risky_behavior_level), names_to="variable") %>%
#   group_by(variable) %>% mutate(value = as.factor(value) %>% as.integer()) %>% ungroup() %>%
#   filter(name == "Time to Vulnerable Infection" |
#            name == "Vulnerable Infected") %>%
#   ggplot() +
#   aes(x=value, y=metric, colour=factor(risky_behavior_level)) +
#   ggbeeswarm::geom_quasirandom(groupOnX=TRUE, dodge.width=0.67, varwidth=TRUE) +
#   geom_smooth(method="lm", formula=y~x, se=FALSE) +
#   # geom_boxplot(aes(group=str_c(risky_behavior_level, value)), alpha=0.75, width=0.33) +
#   facet_grid(name~., scales="free") +
#   scale_colour_manual(values=my_cols[c(1,3,5)] %>% setNames(NULL)) +
#   xlab("- - - - - - - - - - - - - - - Increasing Bar Use/Transmission - - - - - - - - - - - - ->") +
#   theme_bw() +
#   theme(legend.position="bottom",
#         legend.title=element_blank(),
#         axis.text.x=element_blank(),
#         axis.title.y=element_blank())
# ggsave("~/Desktop/Bar_Simulations_Summary.pdf", width=7, height=9)
