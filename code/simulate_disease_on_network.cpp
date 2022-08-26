#include <iostream>
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerMatrix run_disease_simulation_on_network(int maxtime, // maximum time to run simulation for (will stop early if appropriate)
                                                IntegerMatrix contact_network, // contact network adjacency matrix
                                                IntegerVector initial_condition, // initial disease status for each node
                                                LogicalVector vulnerable, // vector indicating which nodes are vulnerable
                                                NumericVector beta, // interaction transmission rates
                                                double eta,   // background transmission rate
                                                double sigma, // disease progression rate
                                                double rho,   // proportion asymptomatic
                                                double gamma, // recovery rate
                                                double mu,    // baseline mortality rate
                                                double nu) {  // disease-induced mortality rate

  int npop = vulnerable.size();
  int E = contact_network.nrow();

  int i, j, t;
  double tmp;

  // keep track of number of exposed/infectious individuals
  int current_number_exposed = 0;
  int current_number_infectious = 0;

  // initialize timeseries matrix
  IntegerMatrix timeseries(npop, maxtime);
  for (i=0; i<npop; i++) {
    timeseries(i,0) = initial_condition[i];
    if (initial_condition[i] == 1) {
      ++current_number_exposed;
    } else if (initial_condition[i] == 2 || initial_condition[i] == 3) {
      ++current_number_infectious;
    }
  }

  // progress through the discrete time epidemic simulation
  for (t=1; t<maxtime; t++) {

    // initialize the next timestep with the previous one
    for (i=0; i<npop; i++) {
      timeseries(i,t) = timeseries(i,t-1);
    }

    // if no exposed or infectious individuals, then we are at a steady state
    if (current_number_exposed == 0 && current_number_infectious == 0) {
      break;
    }

    // Loop through all individuals and update status if necessary
    for (i=0; i<npop; i++) {

      // status: 0=Susceptible, 1=Exposed, 2=Infectious (asymptomatic),
      // 3=Infectious (symptomatic), 4=Recovered, 5=Dead

      // If there are no infectious individuals, there can be no transmission
      if (current_number_infectious != 0 && timeseries(i,t-1) == 0) {
        // loop through all edges
        for (j=0; j<E; j++) {
          // if an edge connects our focal individual
          if ((contact_network(j,0) == i) &&
              // and the connected node is infectious
              (timeseries(contact_network(j,1),t-1) == 2 ||
              timeseries(contact_network(j,1),t-1) == 3) &&
              // then there is a potential for infection
              (runif(1)[0] < beta[j])) {
            ++timeseries(i,t);
            ++current_number_exposed;
            // only need to be infected once to matter
            break;
            // edge direction doesn't matter
          } else if ((contact_network(j,1) == i) &&
            // and the connected node is infectious
            (timeseries(contact_network(j,0),t-1) == 2 ||
            timeseries(contact_network(j,0),t-1) == 3) &&
            // then there is a potential for infection
            (runif(1)[0] < beta[j])) {
            ++timeseries(i,t);
            ++current_number_exposed;
            // only need to be infected once to matter
            break;
          }
        }
        // if no connection has lead to infection, also check for background infection (scaled to population size)
        if (timeseries(i,t) == 0 &&
            runif(1)[0] < (1 - pow(1 - eta / (double) npop, current_number_infectious))) {
          ++timeseries(i,t);
          ++current_number_exposed;
        }
      }

      // since these are all exclusive if statements, we can re-use a single random draw
      tmp = runif(1)[0];
      if (timeseries(i,t-1) == 1 && (tmp < sigma)) {
        if (runif(1)[0] < rho) {
          // E -> Ia
          ++timeseries(i,t);
        } else {
          // E -> Is
          ++ ++timeseries(i,t);
        }
        --current_number_exposed;
        // keep track of number of infectious individuals
        ++current_number_infectious;
      }

      if (timeseries(i,t-1) == 2) {
        if (vulnerable[i]) {
          // Ia -> R (vulnerable)
          if (tmp < (1/3 * gamma)) {
            ++ ++timeseries(i,t);
            --current_number_infectious;
          }
          // Ia -> R (not vulnerable)
        } else if (tmp < gamma) {
          ++ ++timeseries(i,t);
          --current_number_infectious;
        }
      }

      if (timeseries(i,t-1) == 3) {
        if (vulnerable[i]) {
          if (tmp < (1/3 * gamma + mu + nu)) {
            if (tmp < (mu + nu)) {
              // Is -> D (vulnerable)
              ++ ++timeseries(i,t);
            } else {
              // Is -> R (vulnerable)
              ++timeseries(i,t);
            }
            --current_number_infectious;
          }
        } else {
          if (tmp < (gamma + mu)) {
            if (tmp < mu) {
              // Is -> D (not vulnerable)
              ++ ++timeseries(i,t);
            } else {
              // Is -> R (not vulnerable)
              ++timeseries(i,t);
            }
            --current_number_infectious;
          }
        }
      }
    }
  }
  // no need to return timepoints after steady state
  IntegerMatrix::Sub sub_timeseries = timeseries(_,Range(0,t-1));
  return(sub_timeseries);
}
