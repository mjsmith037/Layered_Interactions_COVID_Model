library(magrittr)
library(tidyverse)

MAIN_FONT <- "Open Sans"

# 3d5a6c-568ea3-c84630-e63946-3e885b-433a3f
# b8b3e9-ffb997-885a5a-90b494-718f94
# Colors for plotting
my_cols <- c(
  "S"  = "#0b6884", "Susceptible"                 = "#0b6884",
  "E"  = "#50b99a", "Exposed"                     = "#50b99a",
  "Ia" = "#ff5e5b", "Asymptomatically Infectious" = "#ff5e5b",
  "Is" = "#ffc847", "Symptomatically Infectious"  = "#ffc847",
  "R"  = "#7b678e", "Recovered"                   = "#7b678e",
  "D"  = "#bc4b51", "Dead"                        = "#bc4b51",
  "Vulnerable"     = "#ff5e5b",
  "Not Vulnerable" = "#0b6884",
  "Housemates"          = "#0b6884", "household"         = "#0b6884",
  "Extended Housemates" = "#bc4b51", "between_household" = "#bc4b51",
  "Classmates"          = "#50b99a", "classmate"         = "#50b99a",
  "Coworkers"           = "#7b678e", "coworker"          = "#7b678e",
  "Bar"                 = "#ff5e5b", "bar"               = "#ff5e5b",
  "Interaction"         = "#0b6884",
  "intervention1" = "#b8b3e9",
  "intervention2" = "#ffb997",
  "fullnet"       = "#885a5a",
  "stayathome"    = "#90b494",
  "intervention3" = "#718f94"
)

# my_alphas <- c("household"         = 1,
#                "between_household" = 0.5,
#                "classmate"         = 0.25,
#                "coworker"          = 0.1,
#                "bar"               = 1)

MAXTIME <- 1000
WITHIN_HOUSEHOLD_TRANSMISSION_RATE <- 0.33
