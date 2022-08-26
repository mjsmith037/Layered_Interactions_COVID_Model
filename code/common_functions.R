library(magrittr)
library(tidyverse)
library(ggraph)
library(patchwork)

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

get_background <- function(descriptor) {
  case_when(descriptor == "Very High" ~ 100,
            descriptor == "High"      ~ 10,
            descriptor == "Medium"    ~ 1,
            descriptor == "Low"       ~ 0.1,
            descriptor == "Very Low"  ~ 0.01)
}

plot_network <- function(contact_network) {
  subsetted <- FALSE
  g_order <- with_graph(contact_network, graph_order())
  if (g_order > 500) {
    plot_graph <- contact_network %N>% filter(node_id %in% sample(node_id, 500))
    subsetted <- TRUE
  } else {
    plot_graph <- contact_network
  }

  plot_graph <- plot_graph %E>% mutate(type = recode(type,
                                                     household = "Housemates",
                                                     between_household = "Extended Housemates",
                                                     classmate = "Classmates",
                                                     coworker = "Coworkers",
                                                     bar = "Bar"))

  net <- ggraph(plot_graph, layout=plot_graph %>% igraph::layout_nicely(weight=. %E>% pull(weight))) +
    geom_edge_fan(aes(colour=type),
                  edge_alpha=0.67, edge_width=0.5,
                  end_cap=circle(radius=5, unit="pt"), start_cap=circle(radius=5, unit="pt")) +
    geom_node_point(size=3) +
    scale_edge_colour_manual(values=my_cols) +
    guides(edge_colour=guide_legend(nrow=2, override.aes=list(edge_width=2, edge_alpha=1))) +
    scale_colour_manual(values=my_cols) +
    theme_graph(base_family=MAIN_FONT) +
    theme(
      text=element_text(size=16, family=MAIN_FONT),
      plot.title=element_text(face="bold", size=18),
      plot.subtitle=element_text(size=12),
      plot.caption=element_text(size=12),
      legend.title=element_blank(),
      legend.text=element_text(size=18),
      legend.key.size=unit(18, "pt"),
      legend.position="bottom",
      legend.margin=margin(l=18, r=18)
    ) +
    ggtitle("Social Connectivity", "Note: background transmission links are not shown for clarity")
  if (subsetted) {
    net <- net +
      labs(caption=str_glue("* NOTE: only displaying 500 randomly selected individuals for clarity. Full graph contains {g_order} individuals."))
  }
  return(net)
}

plot_area <- function(timeseries, net_size) {
  area <- ggplot(clean_output(timeseries)) +
    aes(x=time, y=n/net_size*100, fill=status, colour=status) + geom_area() +
    scale_colour_manual(values=my_cols, aesthetics=c("fill", "colour")) +
    coord_cartesian(expand=FALSE) +
    scale_y_continuous(name="Percentage of Population",
                       labels=. %>% str_c("%")) +
    guides(fill=guide_legend(nrow=2)) +
    # sec.axis = sec_axis(~.*g_order/100, name="Number of individuals")) +
    theme_bw() +
    theme(
      text=element_text(size=16, family=MAIN_FONT),
      plot.title=element_text(face="bold", size=18),
      plot.subtitle=element_text(size=12),
      plot.caption=element_text(size=12),
      legend.title=element_blank(),
      legend.text=element_text(size=18),
      legend.key.size=unit(18, "pt"),
      legend.position="bottom",
      legend.margin=margin(l=18, r=18)
    ) +
    ggtitle("Change in proportion of population in each class through time")

  return(area)
}

plot_line <- function(timeseries, net_size) {
  line <- ggplot(timeseries %>%
                   # mutate(vulnerable = sim_output[[1]] %N>% pull(vulnerable)) %>%
                   clean_output() %>% #grouping_vars=c("status", "vulnerable")) %>%
                   filter(status != "Susceptible", status != "Recovered")) +
    aes(x=time, y=n, colour=status) +
    geom_line(size=1.5) +
    scale_colour_manual(values=my_cols) +
    guides(colour=guide_legend(nrow=2)) +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(name="Number of individuals",
                       limits=c(0, NA),
                       sec.axis = sec_axis(~./net_size*100,
                                           labels=. %>% str_c("%"),
                                           name="Percentage of Population")) +
    # facet_wrap(~vulnerable, scale="free_y") +
    theme_bw() +
    theme(
      text=element_text(size=16, family=MAIN_FONT),
      plot.title=element_text(face="bold", size=18),
      plot.subtitle=element_text(size=12),
      plot.caption=element_text(size=12),
      legend.title=element_blank(),
      legend.text=element_text(size=18),
      legend.key.size=unit(18, "pt"),
      legend.position="bottom",
      legend.margin=margin(l=18, r=18)
    ) +
    ggtitle("Proportion of population in disease classes through time")

  return(line)
}
