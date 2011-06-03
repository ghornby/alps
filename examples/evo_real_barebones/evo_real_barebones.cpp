/*******************************************************

  File: evo_real_barebones.cpp

  This file serves as a minimalist barebones example of how to setup
  and run an ALPS population.

 --

  This is part of a software package which implements the
  Age-Layered Population Structure (ALPS) algorithm.

  Copyright (C) [2009] [Gregory S. Hornby]
********************************************************/

#include "alps.h"
#include "alps_gen.h"
#include "alps_layer.h"
#include "alps_individ_real.h"
using namespace alps;

#include <cmath>
#include <cstdio>
#include <string>
using namespace std;

/*
  Functions
  =========
*/
double f102(double x, double y)
{
  // F102 or Rana Function.
  // f102(x,y) = xsin(sqrt(|y+1-x|))cos(sqrt|x+y+1|)
  //             + (y+1)cos(sqrt|y+1-x|)sin(sqrt|x+y+1|)
  // Note this implementation is NOT shifted or rotated.
  // An improved test function would do both.

  double sum1 = y+1.0-x;
  sum1 = sqrt(abs(sum1));

  double sum2 = x+y+1.0;
  sum2 = sqrt(abs(sum2));

  double fitness = x*sin(sum1)*cos(sum2);
  fitness += (y+1.0)*cos(sum1)*sin(sum2);
  //  fitness += 511.707735; // To put the minimum at 0.
  return fitness;
}

double f102_func(const vector<double> params)
{
  double num = 1.0;
  double fitness = f102(params.back(), params[0]);

  for (vector<double>::size_type i=1; i<params.size(); i++) {
    num += 1.0;
    double val = f102(params[i-1], params[i]);
    fitness += val;
  }

  fitness /= num;

  return fitness;  
}

bool evaluate_indreal(vector<double>& fitness, const vector<double>& genes)
{
  fitness[0] = f102_func(genes);

  return true;
}

bool evaluate_individ(vector<double>& fitness, Individual* individ)
{
  if (individ == NULL) {
    cerr << "evo_real :: evaluate_individ() - error, null individ sent.\n";
    while (1) ;
    return false;
  }

  bool res = evaluate_indreal(fitness, ((Individ_Real*)individ)->get_genes());

  return res;
}

void setup_pop_gen(Individual* individ_config, AlpsGen* pop)
{
  (void) individ_config;
  pop->set_save_best(true);
  AlpsLayer layer_def;
  int type = 1;
  int age_gap = 5;
  int age_scheme = ALPS_AGING_FIBONACCI1;
  int Number_Layers = -1;
  if (type == 1) {
    // Configuration for a regular EA/GA:
    Number_Layers = 1;
    layer_def.set_select_type(ALPS_SELECT_TOURN);
    layer_def.set_size(400);
    layer_def.set_elitism(5);
    layer_def.set_tourn_size(5);
    pop->set_recomb_prob(0.5);
    pop->set_rec_rand2_prob(1.0); // 1.0
    pop->set_print_results_rate(400); // 400

  } else if (type == 2) {
    Number_Layers = 10;
    age_gap = 3; //4
    age_scheme = ALPS_AGING_EXP3;
    layer_def.set_select_type(ALPS_SELECT_TOURN);
    //    layer_def.set_select_type(ALPS_SELECT_DC);
    layer_def.set_size(40);
    layer_def.set_elitism(5);
    layer_def.set_tourn_size(5);
    layer_def.set_prob_select_prev(0.2);
    pop->set_recomb_prob(0.5);
    pop->set_rec_rand2_prob(1.0);
    pop->set_print_results_rate(1);

  } else if (type == 3) {
    Number_Layers = 10;
    age_gap = 3; //4
    age_scheme = ALPS_AGING_EXP3;
    //    layer_def.set_select_type(ALPS_SELECT_TOURN);
    layer_def.set_select_type(ALPS_SELECT_DC);
    layer_def.set_size(40);
    layer_def.set_elitism(0);
    layer_def.set_tourn_size(5);
    layer_def.set_prob_select_prev(0.25);
    pop->set_recomb_prob(0.5);
    pop->set_rec_rand2_prob(1.0);
    pop->set_print_results_rate(400);

  } else {
    cerr << "evo_real_barebones :: setup_pop_gen() - error, invalid EA type:"
	 << type << "\n";
    return;
  }

  pop->config_layers_same(age_scheme, age_gap,
			  Number_Layers, layer_def);
  pop->print_layers();
  pop->set_num_runs(1);
  pop->set_max_gen(2525);
}

void *ea_engine(void *arg1)
{
  (void) arg1; // getting rid of unused parameter warning.
  cout << "EA engine started.\n";

  int Number_Genes = 20;
  
  Individ_Real *individ_config = new Individ_Real(Number_Genes);
  individ_config->set_init_minmax(-512.0, 511.0);
  individ_config->set_minmax(-512.0, 511.0);

  vector<double> fitness;
  fitness.resize(1);

  // Configure a generational ALPS population:
  Alps *Population = new AlpsGen("real", individ_config);
  setup_pop_gen(individ_config, (AlpsGen*)Population);

  // Population->set_print_debug(true);
  //Population->set_minimize();
  Population->set_maximize();
  Population->write_header(cout);

  while (!Population->is_finished()) {
    int index;
    Individual* individ;
    int res = Population->get_next_individ(index, individ);
    if (res == -1) {
      continue; // Get another index / Try again.

    } else if (res == -2) {
      // Evolution is over.
      break;
    }

    vector<double> fitness;
    fitness.resize(1);
    int result = evaluate_individ(fitness, individ);
    if (result == false) {
      // Error evaluating this individual.
      Population->evaluate_error(index, individ);
    } else {
      // Evaluated successfully.
      Population->insert_evaluated(fitness, index, individ, 0);
    }
  }

  printf("EA engine ended.\n");

  return 0;
}

int main(int argc, char **argv) {
    
  (void) argc;
  (void) argv;
  ea_engine(0);

  return 0;
}


