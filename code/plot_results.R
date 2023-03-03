library(magrittr)
library(tidyverse)
library(patchwork)
library(broom)
library(scales)

background <- as.numeric(commandArgs(TRUE)[1]) # 0.001
dir.create(file.path("../figures", background))

theme_set(theme_bw())
axis_breaks <- c(0.001, 0.01, 0.1)

my_cols <- c(
  # locale colors
  "Florida" = "#0b6884",
  "Texas"   = "#bc4b51",
  # risk-tolerance regime colors (for correlations matrix plot)
  "III"     = "#bc4b51",
  "II"      = "#7b678e",
  "I"       = "#0b6884",
  "Overall" = "#000000"
)

pos_lm <- possibly(~lm(...) %>% tidy() %>% filter(term != "(Intercept)"), otherwise=tibble())

get_results <- function(background) {
  bind_rows(read_csv(str_glue("../results/{background}_realrisk_full_Florida.csv"),
                     col_types="ddddddddddddddfddffd"),
            read_csv(str_glue("../results/{background}_realrisk_full_Texas.csv"),
                     col_types="ddddddddddddddfddffd")) %>%
    # remove because constant
    select(!any_of(c("background_transmission_rate", "pseudo_risk", "rnaught", "rseed", "number_of_households"))) %>%
    mutate(risky_behavior_level = factor(risky_behavior_level, levels=0:2,
                                         labels=c("Vulnerable Households\nAvoid Work/School",
                                                  "Vulnerable Individuals\nAvoid Work/School",
                                                  "No Difference\nin Behaviour")) %>% fct_rev()) %>%
    relocate(pop_size, vuln_pop_size, .after=everything()) %>%
    rename("Time to\nVulnerable Infection"="Time to Vulnerable Infection",
           "Proportion Vulnerable\nHouseholds Infected"="Percent Vulnerable Households Infected",
           "Peak Infectious"="Maximum Infectious")
}

make_locale_comp_fig <- function(background) {
  get_results(background)  %>%
    filter(`Total Infected` / pop_size > 0.05,
           risky_behavior_level == "No Difference\nin Behaviour") %>%
    mutate(`Epidemic Rate`                   = `Peak Infectious` / `Time to Peak`,
           `Maximum\nProportion Infectious`  = `Peak Infectious` / pop_size,
           `Average\nProportion Infectious`  = `Average Infectious` / pop_size,
           `Proportion Died`                 = `Total Died` / pop_size,
           `Proportion Infected`             = `Total Infected` / pop_size,
           `Proportion\nVulnerable Infected` = `Vulnerable Infected` / vuln_pop_size,
           across(matches("Time to"), . %>% add(1) %>% log()), .keep="unused") %>%
    rename_with(~str_glue("Log {.}"), matches("Time to")) %>%
    pivot_longer(!c(classmate_transmission_rate, coworker_transmission_rate,
                    risky_behavior_level, locale), names_to="variable") %>%
    arrange(variable) %>%
    mutate(variable = factor(variable) %>% fct_inorder()) %>%
    ggplot() +
    aes(y=value, x=locale, colour=locale, fill=locale) +
    geom_violin(show.legend=FALSE) +
    geom_boxplot(fill="white", colour="black", width=0.2, alpha=0.67, outlier.colour=NA, show.legend=FALSE) +
    facet_wrap(~variable, scales="free_y", nrow=2) +
    scale_colour_manual(values=my_cols, aesthetics=c("colour", "fill")) +
    theme(axis.title=element_blank())
  ggsave(str_glue("../figures/{background}/locale_comparison.pdf"), width=10, height=5)
}

make_risk_comp_fig <- function(background) {
  risk_behavior_data <- get_results(background) %>%
    mutate(`Epidemic Rate` = `Peak Infectious` / `Time to Peak`) %>%
    select(!c(locale, pop_size, vuln_pop_size)) %>%
    mutate(across(!c(`Epidemic Rate`, classmate_transmission_rate,
                     coworker_transmission_rate, risky_behavior_level),
                  . %>% add(1) %>% log())) %>%
    rename_with(~str_glue("Log {.}"),
                !c(`Epidemic Rate`, classmate_transmission_rate,
                   coworker_transmission_rate, risky_behavior_level)) %>%
    pivot_longer(!c(classmate_transmission_rate, coworker_transmission_rate,
                    risky_behavior_level), names_to="variable") %>%
    arrange(variable) %>%
    mutate(variable = factor(variable) %>% fct_inorder()) %>%
    group_by(variable) %>%
    mutate(value=rescale(value),
           classmate_transmission_rate = log(classmate_transmission_rate),
           coworker_transmission_rate = log(coworker_transmission_rate)) %>%
    ungroup()
  risk_plot <- function(dat) {
    ggplot(dat) +
      aes(x=classmate_transmission_rate, y=coworker_transmission_rate, z=value, colour=after_stat(value)) +
      stat_summary_hex(fun=mean, bins=30) +
      facet_grid(risky_behavior_level~variable) +
      coord_equal() +
      scale_x_continuous(name="Classmate Transmission Rate",
                         breaks=log(axis_breaks), labels=axis_breaks, expand=expansion(add=c(-0.13, 0.013))) +
      scale_y_continuous(name="Co-worker Transmission Rate",
                         breaks=log(axis_breaks), labels=axis_breaks, expand=expansion(add=c(-0.13, 0.013))) +
      scale_fill_viridis_c(na.value="transparent", aesthetics=c("colour", "fill"),
                           breaks=c(0.05, 0.35, 0.65, 0.95), labels=c("Low", "", "", "High"),
                           guide=guide_colourbar(title="Rescaled\noutcome\nmeasure")) +
      theme(legend.position="none")
  }
  risk_plot(risk_behavior_data %>% filter(as.integer(variable) <= 5)) +
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank()) +
    risk_plot(risk_behavior_data %>% filter(as.integer(variable) > 5)) +
    plot_layout(ncol=1) &
    theme(plot.margin=margin())

  ggsave(str_glue("../figures/{background}/risk_comparison.pdf"), width=10.5, height=12)
}

make_inttype_comp_fig <- function(background) {
  get_results(background) %>%
    filter(`Total Infected` / pop_size > 0.05) %>%
    mutate(`Epidemic Rate` = `Peak Infectious` / `Time to Peak`) %>%
    select(!c(pop_size, vuln_pop_size)) %>%
    mutate(across(!c(`Epidemic Rate`, locale, classmate_transmission_rate,
                     coworker_transmission_rate, risky_behavior_level),
                  . %>% add(1) %>% log())) %>%
    rename_with(~str_glue("Log {.}"),
                !c(`Epidemic Rate`, locale, classmate_transmission_rate,
                   coworker_transmission_rate, risky_behavior_level)) %>%
    mutate(`Log Time to Peak Death` = ifelse(is.infinite(`Log Time to Peak Death`), NA, `Log Time to Peak Death`),
           `Log Time to\nVulnerable Infection` = ifelse(is.infinite(`Log Time to\nVulnerable Infection`), NA, `Log Time to\nVulnerable Infection`)) %>%
    pivot_longer(!c(classmate_transmission_rate, coworker_transmission_rate,
                    risky_behavior_level, locale), names_to="variable") %>%
    arrange(variable) %>%
    mutate(variable = factor(variable) %>% fct_inorder()) %>%
    group_by(locale, risky_behavior_level, variable) %>%
    group_modify(~pos_lm(value~log(classmate_transmission_rate) + log(coworker_transmission_rate), data=.)) %>%
    mutate(term = factor(term,
                         levels=c("log(classmate_transmission_rate)", "log(coworker_transmission_rate)"),
                         labels=c("classmate", "co-worker"))) %>%
    na.omit() %>%
    ggplot() +
    aes(x=term, y=estimate, ymax=estimate + std.error, ymin=estimate - std.error,
        colour=locale, shape=risky_behavior_level) +
    geom_line(aes(x=as.integer(as.factor(term))), position=position_dodge2(width=0.15, reverse=TRUE)) +
    geom_linerange(linewidth=1, position=position_dodge2(width=0.15, reverse=TRUE)) +
    geom_point(size=3, position=position_dodge2(width=0.15, reverse=TRUE)) +
    facet_wrap(~variable, scales="free_y", nrow=2) +
    ylab(expression(atop("Linear model slope estimates",
                         paste(( y %~% m[1]*log(beta[classmate]) + m[2]*log(beta[co-worker]) + b))))) +
    scale_colour_manual(values=my_cols) +
    scale_x_discrete(expand=expansion(add=0.5), labels=parse_format()) +
    theme(legend.title=element_blank(),
          legend.position="bottom",
          axis.title.x=element_blank())
  ggsave(str_glue("../figures/{background}/inttype_comparison.pdf"), width=10, height=5)
}

make_pseudorisk_fig <- function(background) {
  risk_behavior_data <- bind_rows(read_csv(str_glue("../results/{background}_pseudorisk_full_Florida.csv"),
                                           col_types="ddddddddddddddfddffd"),
                                  read_csv(str_glue("../results/{background}_pseudorisk_full_Texas.csv"),
                                           col_types="ddddddddddddddfddffd")) %>%
    select(!c(background_transmission_rate, pseudo_risk, rnaught, pop_size, locale,
              vuln_pop_size, number_of_households, rseed)) %>%
    mutate(risky_behavior_level = factor(risky_behavior_level, levels=0:2,
                                         labels=c("Vulnerable Households\nAvoid Work/School",
                                                  "Vulnerable Individuals\nAvoid Work/School",
                                                  "No Difference\nin Behaviour")) %>% fct_rev()) %>%
    rename("Time to\nVulnerable Infection"="Time to Vulnerable Infection",
           "Proportion Vulnerable\nHouseholds Infected"="Percent Vulnerable Households Infected",
           "Peak Infectious"="Maximum Infectious") %>%
    mutate(`Epidemic Rate` = `Peak Infectious` / `Time to Peak`) %>%
    mutate(across(!c(`Epidemic Rate`, classmate_transmission_rate,
                     coworker_transmission_rate, risky_behavior_level),
                  . %>% add(1) %>% log())) %>%
    rename_with(~str_glue("Log {.}"),
                !c(`Epidemic Rate`, classmate_transmission_rate,
                   coworker_transmission_rate, risky_behavior_level)) %>%
    pivot_longer(!c(classmate_transmission_rate, coworker_transmission_rate,
                    risky_behavior_level), names_to="variable") %>%
    arrange(variable) %>%
    mutate(variable = factor(variable) %>% fct_inorder()) %>%
    group_by(variable) %>%
    mutate(value=rescale(value),
           classmate_transmission_rate = log(classmate_transmission_rate),
           coworker_transmission_rate = log(coworker_transmission_rate)) %>%
    ungroup() %>%
    mutate(risky_behavior_level=factor(risky_behavior_level,
                                       levels=c("No Difference\nin Behaviour",
                                                "Vulnerable Individuals\nAvoid Work/School",
                                                "Vulnerable Households\nAvoid Work/School"),
                                       labels=c("No links removed",
                                                "Links removed randomly\n(as many as when Vulnerable\nindividuals avoid work/school)",
                                                "Links removed randomly\n(as many as when Vulnerable\nhouseholds avoid work/school)")))
  risk_plot <- function(dat) {
    ggplot(dat) +
      aes(x=classmate_transmission_rate, y=coworker_transmission_rate, z=value, colour=after_stat(value)) +
      stat_summary_hex(fun=mean, bins=30) +
      facet_grid(risky_behavior_level~variable) +
      coord_equal() +
      scale_x_continuous(name="Classmate Transmission Rate",
                         breaks=log(axis_breaks), labels=axis_breaks, expand=expansion(add=c(-0.13, 0.013))) +
      scale_y_continuous(name="Co-worker Transmission Rate",
                         breaks=log(axis_breaks), labels=axis_breaks, expand=expansion(add=c(-0.13, 0.013))) +
      scale_fill_viridis_c(na.value="transparent", aesthetics=c("colour", "fill"),
                           breaks=c(0.05, 0.35, 0.65, 0.95), labels=c("Low", "", "", "High"),
                           guide=guide_colourbar(title="Rescaled\noutcome\nmeasure")) +
      theme(legend.position="none")
  }
  risk_plot(risk_behavior_data %>% filter(as.integer(variable) <= 5)) +
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank()) +
    risk_plot(risk_behavior_data %>% filter(as.integer(variable) > 5)) +
    plot_layout(ncol=1) &
    theme(plot.margin=margin())
  ggsave(str_glue("../figures/{background}/pseudorisk_comparison.pdf"), width=10.5, height=12)
}

make_outcome_corplot <- function(background) {
  my_dens <- function(data, mapping, ...) {ggplot(data) + mapping + geom_density(..., alpha=0.4)}
  my_points <- function(data, mapping, ...) {ggplot(data) + mapping + geom_point(..., size=0.1)}
  mycorrelations <- function(data,mapping,...){
    data2 = data
    data2$x = as.numeric(data[,as_label(mapping$x)])
    data2$y = as.numeric(data[,as_label(mapping$y)])
    data2$group = data[,as_label(mapping$colour)]

    correlation_df = data2 %>%
      bind_rows(data2 %>% mutate(group="Overall")) %>%
      group_by(group) %>%
      filter(sum(!is.na(x),na.rm=T)>1) %>%
      filter(sum(!is.na(y),na.rm=T)>1) %>%
      summarize(estimate = round(as.numeric(cor.test(x,y,method="pearson")$estimate),2),
                pvalue = cor.test(x,y,method="pearson")$p.value,
                pvalue_star = as.character(symnum(pvalue, corr = FALSE, na = FALSE,
                                                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                                  symbols = c("***", "**", "*", "'", " "))))%>%
      group_by() %>%
      mutate(group = factor(group, levels=c(as.character(unique(sort(data[,as_label(mapping$colour)]))), "Overall")))

    ggplot(data=correlation_df) +
      aes(x=1, y=group, color=group) +
      geom_text(aes(label=str_glue("{group}: {estimate}{pvalue_star}")), size=4)
  }
  corplot <- get_results(background) %>%
    filter(`Total Infected` / pop_size > 0.05) %>%
    mutate(`Epidemic Rate` = `Peak Infectious` / `Time to Peak`,
           `Time to Peak Death` = ifelse(is.infinite(`Time to Peak Death`), NA, `Time to Peak Death`),
           `Time to\nVulnerable Infection` = ifelse(is.infinite(`Time to\nVulnerable Infection`),
                                                    NA, `Time to\nVulnerable Infection`)) %>%
    select(!c(locale, pop_size, vuln_pop_size, classmate_transmission_rate, coworker_transmission_rate)) %>%
    mutate(across(!c(`Epidemic Rate`, risky_behavior_level), . %>% add(1) %>% log())) %>%
    rename_with(~str_glue("Log {.}"), !c(`Epidemic Rate`, risky_behavior_level)) %>%
    mutate(risky_behavior_level = factor(risky_behavior_level, labels=c("I", "II", "III"))) %>%
    relocate(risky_behavior_level, .after=everything()) %>%
    # na.omit() %>%
    # we do not load the GGally library because it annoyingly loads plyr, which
    # does not play well with dplyr
    GGally::ggpairs(mapping    = aes(fill=risky_behavior_level, colour=risky_behavior_level),
            diag       = list(continuous=my_dens),
            lower      = list(continuous=my_points),
            upper      = list(continuous=mycorrelations),
            axisLabels = "show",
            columns=1:10) +
    scale_colour_manual(values=my_cols, aesthetics=c("colour", "fill"))
  ggsave(corplot, filename=str_glue("../figures/{background}/outcome_correlations.png"), width=17, height=17)
}

make_R0_fig <- function(background) {
  bind_rows(read_csv(str_glue("../results/{background}_realrisk_R0_Florida.csv"),
                     col_types="ddddddddddddddfddffd"),
            read_csv(str_glue("../results/{background}_realrisk_R0_Texas.csv"),
                     col_types="ddddddddddddddfddffd")) %>%
    select(`Total Infected`, classmate_transmission_rate, coworker_transmission_rate, pop_size, locale, risky_behavior_level) %>%
    mutate(risky_behavior_level = factor(risky_behavior_level, levels=0:2,
                                         labels=c("Vulnerable Households\nAvoid Work/School",
                                                  "Vulnerable Individuals\nAvoid Work/School",
                                                  "No Difference in Behavior")) %>% fct_rev()) %>%
    mutate(value=`Total Infected` - 1, # %>% add(1) %>% log() %>% rescale(),
           classmate_transmission_rate = log(classmate_transmission_rate),
           coworker_transmission_rate = log(coworker_transmission_rate)) %>%
    ungroup() %>%
    mutate(value = ifelse(value > 20, 20, value)) %>%
    ggplot() +
    aes(x=classmate_transmission_rate, y=coworker_transmission_rate, z=value, colour=after_stat(value)) +
    stat_summary_hex(fun=mean, bins=30) +
    facet_grid(risky_behavior_level~.) +
    coord_equal() +
    scale_x_continuous(name="Classmate Transmission Rate",
                       breaks=log(axis_breaks), labels=axis_breaks, expand=expansion(add=c(-0.13, 0.013))) +
    scale_y_continuous(name="Co-worker Transmission Rate",
                       breaks=log(axis_breaks), labels=axis_breaks, expand=expansion(add=c(-0.13, 0.013))) +
    scale_fill_viridis_c(na.value="transparent", aesthetics=c("colour", "fill"),
                         option="B", name="Number of\nsecondary\ninfections",
                         breaks=c(0, 5, 10, 15, 20), labels=c("0", "5", "10", "15", expression(phantom()>=20)))
  ggsave(str_glue("../figures/{background}/R0_heatmap.pdf"), width=3.8, height=6)
}

make_locale_comp_fig(background)
make_risk_comp_fig(background)
make_inttype_comp_fig(background)
make_outcome_corplot(background) # note this figure uses the ggpairs function from GGally

if (file.exists(str_glue("../results/{background}_pseudorisk_full_Florida.csv")) &&
    file.exists(str_glue("../results/{background}_pseudorisk_full_Texas.csv"))) {
  make_pseudorisk_fig(background)}

if (file.exists(str_glue("../results/{background}_realrisk_R0_Florida.csv")) &&
    file.exists(str_glue("../results/{background}_realrisk_R0_Texas.csv"))) {
  make_R0_fig(background)}
