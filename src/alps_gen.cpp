/*******************************************************

  File: alps_gen.cpp
  This is part of a software package which implements the
  Age-Layered Population Structure (ALPS) algorithm.

  Copyright (C) [2009] [Gregory S. Hornby]

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

********************************************************/

#include <algorithm>
#include <stdexcept>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
using namespace std;


#include "alps/alps_random_mt.h"
using namespace alps_random;

#include "alps/alps_utils.h"
using namespace alps_utils;

#include "alps/alps.h"
#include "alps/alps_gen.h"
#include "alps/alps_individ.h"
#include "alps/alps_layer.h"
#include "alps/alps_history.h"

namespace alps {


const int MAX_TRIES = 1;
const int ALPS_GEN_VERSION = 1;


// Evaluation States
const int ES_NOT_USED = -5;
const int ES_NOT_CREATED = -4;
const int ES_EVALUATING = 1;
const int ES_EVALUATED = 2;


/***** Initialization functions *****/


AlpsGen :: AlpsGen(const char *n, Individual* ind_config) :
  Alps(n, ind_config)
{
  update_correlation_func = 0;
  set_name(n);
  
  seed_individ_ = 0;

  init_variables();
}


AlpsGen :: AlpsGen(Individual* ind_config) :
  Alps("evolve", ind_config)
{
  update_correlation_func = 0;
  set_name("evolve");

  seed_individ_ = 0;

  init_variables();
}


AlpsGen :: ~AlpsGen()
{
  for (vector<Individual*>::size_type i=0; i<pop_parents_.size(); i++) {
    delete pop_parents_[i];
  }

  for (vector<Individual*>::size_type i=0; i<pop_childs_.size(); i++) {
    delete pop_childs_[i];
  }
}



// Probably want to log:
// a. layer_def stats
// b. structural organization stats.


void AlpsGen :: set_pop_size(unsigned int size)
{
  //  cout << "alps_gen :: set_pop_size(" << size << ")\n";

  /*
  if (size < 2) {
    cout << "alps_gen -- error, size < 2.\n";
    while (1) ;
  }
  */
  max_layer_evolving_ = num_layers_-1; // 0;

  evaluation_list_.resize(size);
  tries_list_.resize(size);
  create_type_list_.resize(size);
  parent1_list_.resize(size);
  parent2_list_.resize(size);
  tourn_list_.resize(size);
  evaluation_state_.resize(size);

  vector<Individual*>::size_type old_size = pop_parents_.size();

  for (vector<Individual*>::size_type i=size; i<pop_parents_.size(); i++) {
    delete pop_parents_[i];
    delete pop_childs_[i];
  }

  pop_parents_.resize(size);
  pop_childs_.resize(size);
  pop_order_.resize(size);

  for (vector<Individual*>::size_type i=old_size; i<size; i++) {
    pop_parents_[i] = sample_individ_->new_instance();
    pop_childs_[i] = sample_individ_->new_instance();
  }
  num_individs_ = size;
}


void AlpsGen :: init_variables()
{
  save_best_ = true;
  save_log_ = false;

  reval_duplicates_ = false;
  moveup_dead_parents_ = true;
  moveup_old_elite_ = true;

  num_individs_ = 0;
  max_generations_ = 100;
  prob_recomb_ = 0.5; // 0.5
  run_sample_ = false;
  num_runs_ = 1;
  run_number_ = 0;
  print_results_rate_ = 100;

  num_layers_ = 0;
  max_layer_evolving_ = 0;

  AlpsLayer def;
  def.set_elitism(0);
  def.set_select_type(ALPS_SELECT_TOURN);
  def.set_max_age(-1);
  def.set_start_index(0);
  def.set_size(1);

  vector<int> prev_layers;
  def.set_prev_layers(0, prev_layers);
  config_layer(def);

  save_log_ = false;
}


void AlpsGen :: initialize()
{
  //  cout << "alps_gen :: initialize()\n";
  num_evaluations_ = 0;
  num_evaluations_init_ = 0;
  num_evaluations_evolve_ = 0;
  
  // Set state to start evaluating first individual.
  tries_ = 0;
  evaluated_ = false;
  create_type_ = CREATE_RANDOM;
  current_individual_ = 0;

  for (int i=0; i<num_individs_; i++) {
    pop_childs_[i]->initialize();
    pop_parents_[i]->initialize();
  }

  generation_ = 0;
  if (seed_individ_) {
    generation_ = 1;
    //    setup_next_individual(); // <= delete this function.

  } else if (!run_sample_) {
    pop_childs_[0]->make_random();
  }

  for (int i=0; i<num_individs_; i++) {
    create_type_list_[i] = CREATE_NOT;
  }
  setup_next_generation_layers();  
}



int AlpsGen :: get_size()
{
  return num_individs_;
}

int AlpsGen :: get_num_evals()
{
  return num_evaluations_;
}


void AlpsGen :: set_moveup_dead_parents(bool b)
{
  moveup_dead_parents_ = b;
}

void AlpsGen :: set_moveup_old_elite(bool b)
{
  moveup_old_elite_ = b;
}


void AlpsGen :: set_tourn_size(unsigned int size)
{
  for (vector<AlpsLayer>::size_type i=0; i<layer_def_.size(); i++) {
    layer_def_[i].set_tourn_size(size);
  }
}

void AlpsGen :: set_max_gen(int n)
{
  max_generations_ = n;
}

void AlpsGen :: set_num_runs(int n)
{
  num_runs_ = n;
}

void AlpsGen :: set_run_number(int r)
{
  run_number_ = r;
}


void AlpsGen :: set_recomb_prob(double prob)
{
  prob_recomb_ = prob;
}

void AlpsGen :: set_rec_rand2_prob(double prob)
{
  // This variable is the probability that the second individual
  // selected for recombination is selected at random instead of
  // using the same method for choosing the first individual
  // (eg the first is likely chosen via tournament or roulette wheel).
  prob_rec_rand2_ = prob;
}


//bool AlpsGen :: configure(char *fname, int verbose, int (*translator_func)(char*))
bool AlpsGen :: configure(char *, int, int (*)(char*))
{
  cerr << "AlpsGen :: configure() - no longer supported..\n";
  while (1) ;

  return true;
}


void AlpsGen :: print_pop()
{
  int p = cout.precision();
  int w = cout.width();

  for (int i=0; i<num_individs_; i++) {
    cout << " " << setw(3) << i << " : "
	 << "P: " << setw(7) << pop_parents_[i]->get_fitness()
	 << " (" << pop_parents_[i]->get_age()
	 << ") - ";

    if (create_type_list_[i] == CREATE_RANDOM) {
      cout << " Random      ";
    } else if (create_type_list_[i] == CREATE_DUPLICATE) {
      cout << " Duplicate   ";
    } else if (create_type_list_[i] == CREATE_MUTATION) {
      cout.width(3);
      cout << " Mut(" << setw(3) << parent1_list_[i] << ")    ";
    } else if (create_type_list_[i] == CREATE_RECOMBINATION) {
      cout.width(3);
      cout << " Rec(" << setw(3) << parent1_list_[i] << ":"
	   << setw(3) << parent2_list_[i] << ")";
    } else if (create_type_list_[i] == CREATE_NOT) {
      cout << " Not Created ";
    } else if (create_type_list_[i] == CREATE_WAITING) {
      cout << " Waiting     ";
    } else {
      cout << " Unknown type";
      while (1) ;
    }

    cout << " - C: " << setw(7) << pop_childs_[i]->get_fitness()
	 << " (" << pop_childs_[i]->get_age()
	 << ")\n";
  }

  cout.precision(p);
  cout.width(w);
}



void AlpsGen :: print_parents()
{
  cout << "alps :: print_parents()\n";
  for (int i=0; i<num_individs_; i++) {
    cout << " " << setw(3) << i << " : ";

    if (create_type_list_[i] == CREATE_RANDOM) {
      cout << "Rand";

    } else if (create_type_list_[i] == CREATE_DUPLICATE) {
      cout << "Dup - ";
      int p1 = parent1_list_[i];
      cout << setw(3) << p1 << ":" << setw(7) << pop_parents_[p1]->get_fitness() << ":"
	   << setw(2) << pop_parents_[p1]->get_age();

    } else if (create_type_list_[i] == CREATE_MUTATION) {
      cout << "Mut - ";
      int p1 = parent1_list_[i];
      cout << setw(3) << p1 << ":" << setw(7) << pop_parents_[p1]->get_fitness() << ":"
	   << setw(2) << pop_parents_[p1]->get_age();

    } else if (create_type_list_[i] == CREATE_RECOMBINATION) {
      cout << "Rec - ";
      int p1 = parent1_list_[i];
      cout << setw(3) << p1 << ":" << setw(7) << pop_parents_[p1]->get_fitness() << ":"
	   << setw(2) << pop_parents_[p1]->get_age() << " + ";
      int p2 = parent2_list_[i];
      cout << setw(3) << p2 << ":" << setw(7) << pop_parents_[p2]->get_fitness() << ":"
	   << setw(2) << pop_parents_[p2]->get_age();

    } else if (create_type_list_[i] == CREATE_NOT) {
      cout << "Not";

    } else if (create_type_list_[i] == CREATE_WAITING) {
      cout << "Wtg";

    } else {
      cout << "Unknown\n";
      while (1) ;
    }

    cout << "\n";
  }
}


void AlpsGen :: print_parent_pop()
{
  cout << "alps :: print_parent_pop()\n";
  for (int i=0; i<num_individs_; i++) {
    cout << " " << setw(3) << i << " : "
	 << setw(7) << pop_parents_[i]->get_fitness()
	 << " (" << pop_parents_[i]->get_age()
	 << ") - ";

    if (create_type_list_[i] == CREATE_RANDOM) {
      cout << "Rand";

    } else if (create_type_list_[i] == CREATE_DUPLICATE) {
      int index1 = parent1_list_[i];
      cout << "Dup[" << setw(3) << index1 << "] : ";

    } else if (create_type_list_[i] == CREATE_MUTATION) {
      int index1 = parent1_list_[i];
      cout << "Mut[" << setw(3) << index1 << "] : ";

    } else if (create_type_list_[i] == CREATE_RECOMBINATION) {
      int index1 = parent1_list_[i];
      int index2 = parent2_list_[i];
      cout << "Rec[" << setw(3) << index1 << ":"
	   << setw(3) << index2 << "] : ";

    } else if (create_type_list_[i] == CREATE_NOT) {
      cout << "Not";

    } else if (create_type_list_[i] == CREATE_WAITING) {
      cout << "Wtg";

    } else {
      cout << "Unknown\n";
      while (1) ;
    }

    cout << "\n";
  }
}



void AlpsGen :: print_evaluation_state()
{
  cout << "AlpsGen :: print_evaluation_state()\n";
  for (unsigned int i=0; i< layer_def_[max_layer_evolving_].get_start_index()
	 + layer_def_[max_layer_evolving_].get_size(); i++) {
    cout << i << " evaluation state: " << evaluation_state_[i]
	 << "(create-type:" << create_type_list_[i] << ").\n";
  }
}



AlpsGen::const_iterator AlpsGen :: begin()
{
  return pop_parents_.begin();
}

AlpsGen::const_iterator AlpsGen :: end()
{
  return pop_parents_.end();
}



void AlpsGen:: print_write_gen_stats()
{
  if (print_gen_stats_ && (generation_ > 0)) {
    if ((num_evaluations_ % print_results_rate_) == 0)  {
      if (save_log_) {
	write_log_layers();
      }

      // Print some stats.
      print_generation_stats();
    }
  }
}

// Returns the index of the next individual to evaluate.
// Returns -1 if next individual is not ready yet.
// Returns -2 if all evaluations are done.
int AlpsGen :: get_next_individ(int& index, Individual*& individ)
{
  //  cout << "alps :: get_next_index()\n";

  if (print_debug_info_) {
    cout << "alps_gen :: get_next_individ : pgs=" << print_gen_stats_
	 << ", gen=" << generation_ << ", prr=" << print_results_rate_
	 << ", #Eval=" << num_evaluations_ << "\n";
  }

  if (is_finished_) {
    if (num_being_evaluated_ > 0) {
      // Waiting for an outstanding evaluation.
      return -1;
    }
    if (!setup_next_generation()) {
      print_write_gen_stats();
      return -2;
    }
  }

  if (run_number_ >= num_runs_) {
    if (save_log_) {
      log_file_.close();
      save_log_ = false;
    }
    print_write_gen_stats();
    return -2; // No more evaluations.
  }

  if (next_to_evaluate_ >= num_to_evaluate_) {
    if (num_being_evaluated_ > 0) {
      // Waiting for an outstanding evaluation.
      return -1;
    }
    if (!setup_next_generation()) {
      print_write_gen_stats();
      return -2;
    }

    print_write_gen_stats();
  }

  /*
  cout << "AlpsGen :: get_next_index()  Next: "
       << next_to_evaluate_ << ":" << evaluation_list_[next_to_evaluate_]
       << "\n";
  */

  /*
  if (save_log_ && (num_evaluations_ == 0)) {
    write_header(log_file_);
  }
  */

  int next_index = evaluation_list_[next_to_evaluate_];
  next_to_evaluate_++;
  num_being_evaluated_++;
  make_individual(next_index);

  if (create_type_list_[next_index] == CREATE_DUPLICATE) {
    if (!reval_duplicates_) {
      evaluation_state_[next_index] = ES_EVALUATED;
      num_being_evaluated_--;
      return get_next_individ(index, individ);
    }
  }

  evaluation_state_[next_index] = ES_EVALUATING;

  /*
  cout << "alps :: get_next_index() - index:" << next_index
       <<", evaluation_state=" << evaluation_state_[next_index] << "\n";
  */

  index = next_index;
  individ = pop_childs_[index];
  return 0;
}


Individual *AlpsGen :: get_individual(unsigned int index)
{
  if (index >= pop_childs_.size()) {
    cout << "pop :: get_individual() - invalid index: " << index << "\n";
    while (1) ;
    return 0;
  }

  return pop_childs_[index];
}


Individual *AlpsGen :: get_parent(int index)
{
  if ((index < 0) || (index >= num_individs_)) {
    printf("alps :: get_parent() - invalid index: %d.\n", index);
    while (1) ;
    return 0;
  }

  return pop_parents_[index];
}



Individual *AlpsGen :: get_current_individual()
{
  if (!run_sample_) {
    return pop_childs_[current_individual_];
  }

  //  printf("get_current() returning sample_individ_\n");
  return sample_individ_;
}



// Returns true if no more evaluations to be done.
bool AlpsGen :: is_finished()
{
  if (Alps::is_finished()) {
    return true;
  }

  if (num_runs_ <= 1) {
    if (generation_ > max_generations_) {
      set_finished(true);
      return true;
    }
  } else if (run_number_ >= num_runs_) {
    set_finished(true);
    return true;
  }
  return false;
}


Individual* AlpsGen :: get_current(int index)
{
  if ((index < 0) || (index > num_individs_)) {
    printf("alps :: get_current() - index: %d, is invalid.\n", index);
    while (1) ;
    return 0;
  }

  if (evaluation_state_[index] == ES_EVALUATED) {
    if (print_debug_info_) {
      cout << " get_current("<<index<<") - childs.\n";
    }
    return pop_childs_[index];
  }

  if (print_debug_info_) {
    cout << " get_current("<<index<<") - parent.\n";
  }
  return pop_parents_[index];
}



double AlpsGen :: get_fitness(int index, int layer)
{
  if ((index < 0) || (index > num_individs_)) {
    cout << "alps :: get_fitness() - index: " << index
	 << ", is invalid.\n";
    while (1) ;
    return 0.0;
  }

  if (layer_def_[layer].get_select_type() != ALPS_SELECT_DC) {
    return get_current(index)->get_fitness();
  }

  if ((pop_parents_[index]->get_age() > layer_def_[layer].get_max_age()) &&
      (layer_def_[layer].get_max_age() >= 0)) {
    return pop_childs_[index]->get_fitness();
  }

  if (pop_childs_[index]->compare_fitness(is_maximizing_,
					 pop_parents_[index]) == -1) {
    return pop_parents_[index]->get_fitness();
  }

  return pop_childs_[index]->get_fitness();
}




void AlpsGen :: print_generation_stats()
{
  if (pop_parents_.size() == 0) {
    return;
  }

  int best_index = 0;
  int best_layer = 0;
  double best_overall = pop_parents_[0]->get_fitness();
  vector<double> layer_avg(num_layers_);
  vector<double> layer_best(num_layers_);
  vector<int> layer_best_index(num_layers_);

  /*
  cout << "alps_gen :: print_generation_stats() - #evals:"
       << num_evaluations_ << ":" << num_evaluations_evolve_ << "\n";
  */

  for (int i=0; i<num_layers_; i++) {
    layer_avg[i] = 0.0;
    layer_best[i] = 0.0;
  }

  for (int i=0; i<=max_layer_evolving_; i++) {
    int start_index = layer_def_[i].get_start_index();
    int stop_index = layer_def_[i].get_start_index() + layer_def_[i].get_size();

    double best_fitness = pop_parents_[start_index]->get_fitness();
    double avg_fitness = best_fitness;

    if ((is_maximizing_ && (best_fitness > best_overall)) ||
	(!is_maximizing_ && (best_fitness < best_overall))) {
      best_overall = best_fitness;
      best_index = start_index;
      best_layer = i;
    }
    layer_best_index[i] = start_index;

    for (int j=start_index+1; j<stop_index; j++) {
      double fitness = pop_parents_[j]->get_fitness();
      avg_fitness += fitness;

      if ((is_maximizing_ && (fitness > best_fitness)) ||
	  (!is_maximizing_ && (fitness < best_fitness))) {
	best_fitness = fitness;
	layer_best_index[i] = j;

	if ((is_maximizing_ && (best_fitness > best_overall)) ||
	    (!is_maximizing_ && (best_fitness < best_overall))) {
	  best_overall = best_fitness;
	  best_index = j;
	  best_layer = i;
	}
      }
    }

    layer_avg[i] = avg_fitness;
    layer_best[i] = best_fitness;
  }

  double avg = layer_avg[0];
  int num_individs_ = layer_def_[0].get_size();
  for (int i=1; i<=max_layer_evolving_; i++) {
    avg += layer_avg[i];
    num_individs_ += layer_def_[i].get_size();
  }

  int best_age = pop_parents_[best_index]->get_age();


  if (print_debug_info_) {
    cout << "\n";
  }


  int prec = cout.precision();
  cout.precision(2);
  cout << fixed << num_evaluations_ << "  "
       << best_overall << " "
       << setw(3) << best_age << "  "
       << num_layers_;
  cout.precision(1);


  if (save_log_) {
    log_file_.precision(9);
    log_file_ << num_evaluations_ << " "
	     << best_overall << " " << best_age << " " 
	     << avg / double(num_individs_) << " " << num_layers_;
    log_file_.precision(5);
  }

  for (int i=0; i<=max_layer_evolving_; i++) {
    int index = layer_best_index[i];
    int layer_best_age = pop_parents_[index]->get_age();

    cout << ",  " << layer_best[i] << " " << setw(2) << layer_best_age;
    //    cout << "  " << layer_best[i] << " " << layer_best_age;
    if (save_log_) {
      log_file_ << " " << layer_best[i] << " " << layer_best_age
	       << " " << layer_avg[i]/double(layer_def_[i].get_size());
    }
  }

  for (int i=max_layer_evolving_+1; i<num_layers_; i++) {
    cout << "  0.0 0";
    if (save_log_) {
      log_file_ << " 0.0 0 0.0";
    }
  }
  if (save_log_) {
    log_file_ << endl;
  }

  /*
  pop_childs_[max_index]->print_history();
  pop_childs_[max_index]->print_history_ops();
  */

  cout << " - ";
  //  individual_print_op_history();

  /*
  printf("  : <%2d:%2d:%2d> ", generation_,
	 current_individual_, max_index);
  */
  cout << " ";
  print_time();
  cout << "\n";
  cout.precision(prec);

  if (print_debug_info_) {
    cout << "\n";
  }
}


void AlpsGen :: print_state()
{
  printf("<%2d:%2d>  evals:%d:%d:%d",
	 generation_, current_individual_,
	 num_evaluations_, num_evaluations_init_,
	 num_evaluations_evolve_);

  printf("  -  ");
  individual_print_op_history();
  printf("   ");
  print_time();
  printf("\n");  
}



int AlpsGen :: get_layer(int index)
{
  int count = 0;
  for (int i=0; i<num_layers_; i++) {
    count += layer_def_[i].get_size();
    if (index < count) {
      return i;
    }
  }

  cout << "alps_gen :: get_layer(1) - error, index exceeds maximum layer: "
       << index << ".\n";

  while (1) ;
  return num_layers_;
}



void AlpsGen :: insert_evaluated(vector<double>& fitness, int index,
				 Individual* individ, int print_individ_history)
{
  if ((index < 0) || (index >= num_to_evaluate_)) {
    printf("alps_gen :: insert_evaluated(4) - invalid index: %d (nte:%d.\n", index,
	   num_to_evaluate_);
    while (1) ;
    return;
  }

  if (evaluation_state_[index] != ES_EVALUATING) {
    printf("alps_gen :: insert_evaluated(4) not evaluating index: %d. State=%d.\n",
	   index, evaluation_state_[index]);
    while (1) ;
    return;
  }

  num_being_evaluated_--;


  if (generation_ > 0) {
    num_evaluations_evolve_++;
  } else {
    num_evaluations_init_++;
  }
  num_evaluations_++;

  // Set current individual's fitness.
  tries_list_[index]++;
  evaluated_ = true;
  pop_childs_[index]->set_fitness(is_maximizing_, fitness);
  if (print_individ_history) {
    if ((fitness[0] > 0.0) || ((generation_ > 0) &&
			       (create_type_list_[index] != CREATE_RANDOM))) {
      printf("%2d.%d: ", index, tries_list_[index]);
      pop_childs_[index]->print_history_full();
      printf("\n");
    }
  }


  if (create_type_list_[index] == CREATE_RANDOM) {
    // A randomly generated individual.
    if (pop_childs_[index]->keep(is_maximizing_)) {
      evaluation_state_[index] = ES_EVALUATED;
      if (print_debug_info_) {
	cout << "alps_gen :: insert_evaluated(4) - keep index:" << index << ".\n";
      }

    } else {
      // Do this individual again.
      next_to_evaluate_--;
      evaluation_list_[next_to_evaluate_] = index;
      evaluation_state_[index] = ES_NOT_CREATED;

      if (print_debug_info_) {
	cout << " alps_gen :: insert_evaluated(4) - doing index " << index
	     << " again: fitness=" << individ->get_fitness() << ".\n";
      }
    }

  } else {
    evaluation_state_[index] = ES_EVALUATED;
  }


  /*
  printf("Gen:%d, log:%d; NEE:%d, PRR:%d; NTE:%d, NBE:%d\n",
	 generation_, update_log, num_evaluations_evolve_, print_results_rate_,
	 next_to_evaluate_, num_being_evaluated_);
  if (((generation_ > 0) && (update_log == true) &&
       ((num_evaluations_evolve_ & print_results_rate_) == 0)) ||
      ((generation_ == 0) && (next_to_evaluate_ + num_being_evaluated_ == 0))) {
  */

  /*
  if (((generation_ > 0) && (update_log == true) &&
       ((num_evaluations_evolve_ & print_results_rate_) == 0)) ||
      ((generation_ == 0) && (next_to_evaluate_ + num_being_evaluated_ == 0))) {
  */


  /*
  printf("Gen:%d, log:%d; NEE:%d, PRR:%d; NTE:%d, NBE:%d\n",
	 generation_, update_log, num_evaluations_evolve_, print_results_rate_,
	 next_to_evaluate_, num_being_evaluated_);
  */



  /*
  if (save_log_) {
    if ((num_evaluations_ > 0) &&
	((num_evaluations_ % print_results_rate_) == 0)) {
      cout << "saving log.\n";
      if (!log_file_.is_open()) {
	cout << "Log file not open.\n";
      }
      write_log_layers();
    }
  }
  */
    

  if (save_log_ && (num_evaluations_ == 1)) {
    write_header(log_file_);
  }

}

void AlpsGen :: insert_evaluated(int index, Individual* individ,
				 int print_individ_history)
{
  insert_evaluated(individ->get_fitness_vec(), index, individ,
		   print_individ_history);
}



// Set fitness for networked version.
void AlpsGen :: evaluate_error(int index, Individual* individ)
{
  if (print_debug_info_) {
    cout << "alps_gen() :: evaluate_error(" << index << ") = "
	 << individ->get_fitness() << "\n";
  }


  if ((index < 0) || (index >= num_individs_)) {
    cout << "alps :: evaluate_error() - invalid index: " << index << "\n";
    return;
  }

  if (evaluation_state_[index] != ES_EVALUATING) {
    cout << "alps :: evaluate_error() not evaluating index: " << index << "\n";
    return;
  }

  // Do this individual again.
  num_being_evaluated_--;
  next_to_evaluate_--;  
  evaluation_list_[next_to_evaluate_] = index;
  evaluation_state_[index] = ES_NOT_CREATED;

  if (print_debug_info_) {
    cout << " :: Doing individual again: gen: " << generation_
	 << ", tries_: " << tries_ << "\n";
  }
}




void AlpsGen :: print_final_stats()
{

}


bool AlpsGen :: read_sample_individ(const char *fname)
{
  if (run_sample_) {
    printf("pop :: read_sample_individ() - already loaded individual.\n");
    return false;
  }

  if (sample_individ_ == 0) {
    sample_individ_ = new Individual();
  }

  return sample_individ_->read(fname);
}


bool AlpsGen :: read_seed_individ(const char *fname)
{
  if (sample_individ_ == 0) {
    seed_individ_ = new Individual();
  }
  
  return seed_individ_->read(fname);
}


void AlpsGen :: write_best_individ(std::ofstream& ostr)
{
  if (max_individ_ != 0) {
    printf("writing best:  ");
    max_individ_->write(ostr);
    max_individ_->print_history();

  } else {
    //  Pop_Vec[Max_New_Index]->write(ostr);
    printf("writing best: %d<%d>  ", max_index,
	   current_individual_);
    if (current_individual_ < max_index) {
      pop_parents_[max_index]->write(ostr);
      pop_parents_[max_index]->print_history();

    } else {
      pop_childs_[max_index]->write(ostr);
      pop_parents_[max_index]->print_history();
    }
  }
  printf("\n");
}


void AlpsGen :: write_best_individ()
{
  char fname[60];

  make_filename("_best.ind", fname);


  ofstream outfile(fname);
  write_best_individ(outfile);
  outfile.close();
}


// Writes best individual in the population after it
// has been ordered (thus best is in index pop_order_[0]).
void AlpsGen :: write_best_individ_gen()
{
  int best_index = 0;
  int start_index = 0;
  int stop_index = layer_def_[max_layer_evolving_].get_start_index()
    + layer_def_[max_layer_evolving_].get_size();

  for (int i=start_index; i<stop_index; i++) {
    if (compare_parents(i, best_index) == 1) {
      best_index = i;
    }
  }


  char fname[60];
  make_filename("_best.ind", fname);
  ofstream ofile(fname);

  pop_parents_[best_index]->write(ofile);
  write_header(ofile);

  ofile.close();
}





// Parses a line from the configuration file,
// returns false if a problem, otherwise true.
//bool AlpsGen :: parse_config(const char *config_line, int verbose,
//			     int (*translator_func)(char*))
bool AlpsGen :: parse_config(const char *config_line, int verbose,
			     int (*)(char*))
{
#ifdef OLD
  int i, i_val, loc;
  double f_val;


  printf("pop :: parse_config() - trying: %s", config_line);
  
  for (i=0; i<alps::Max_Line_Length; i++) {
    if (config_line[i] == '\n')
      return true; // Empty line, no error in parsing.
    if (config_line[i] != ' ')
      break;
  }

  printf("pop :: parse_config() - parsing: %s", config_line);

  // AlpsGen Population size.
  if (str_eq("pop size", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("AlpsGen :: Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 1) {
      printf("AlpsGen :: Error, number of individuals must be > 0.\n");
      return false;
    }

    num_individs_ = i_val;
    if (verbose) {
      printf("AlpsGen :: setting number of individuals to: %d.\n",
	     num_individs_);
    }
    return true;
  }

  // Number of generations
  if (str_eq("num generation", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("AlpsGen :: Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 1) {
      printf("AlpsGen :: Error, number of generations must be 0 < %d.\n",
	     i_val);
      return false;
    }

    max_generations_ = i_val;
    if (verbose) {
      printf("AlpsGen :: setting number of generations to: %d.\n", max_generations_);
    }
    return true;
  }

  // Number of runs
  if (str_eq("num runs", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("AlpsGen :: Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 1) {
      printf("AlpsGen :: Error, number of runs must be 0 < %d.\n",
	     i_val);
      return false;
    }

    num_runs_ = i_val;
    if (verbose) {
      printf("AlpsGen :: setting number of runs to: %d.\n", num_runs_);
    }
    return true;
  }

  // Starting run number
  if (str_eq("starting run", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("AlpsGen :: Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 1) {
      printf("AlpsGen :: Error, starting number, %d,  must be >0.\n",
	     i_val);
      return false;
    }

    run_number_ = i_val;
    if (verbose) {
      printf("AlpsGen :: setting run number to: %d.\n", run_number_);
    }
    return true;
  }

  // Recombine rate.
  if (str_eq("recomb rate", config_line)) {
    if (!get_double(config_line, &f_val)) {
      printf("AlpsGen :: Error parsing: %s\n", config_line);
      return false;
    }

    if ((f_val < 0.0f) || (f_val > 1.0f)) {
      printf("AlpsGen :: Error, recombination rate must be 0.0 < %.2f < 1.0.\n",
	     f_val);
      return false;
    }

    prob_recomb_ = f_val;
    if (verbose) {
      printf("AlpsGen :: setting recombination rate to: %f.\n",
	     prob_recomb_);
    }
    return true;
  }

  // Random number seed.
  if (str_eq("random seed", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("AlpsGen :: Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 0) {
      printf("AlpsGen :: Error, random seed, %d < 0.\n",
	     i_val);
      return false;
    }

    if (verbose) {
      printf("AlpsGen :: setting random seed %d.\n", i_val);
    }
    seed_random(i_val);
    srand(i_val);
    while (i_val-- > 0) {
      f_val = random_double();
      i = random_int(100);
    }
    return true;
  }


  // Seed individual.
  if (str_eq("seed ind", config_line)) {
    loc = find_field(config_line);
    i = loc;
    while (config_line[i] != '\n') {
      i++;
      if (i > Max_Line_Length) {
	printf("AlpsGen :: parsing - seed ind, max line length exceeded.\n");
	return false;
      }
    }


    /*
    config_line[i] = 0;
    if (verbose) {
      printf("Seed individual:  %s\n", &(config_line[loc]));
    }
    */

    if (!read_seed_individ(&(config_line[loc]))) {
      return false;
    }

    if (verbose) {
      seed_individ_->write("sample.ind");
      seed_individ_->print_history_full();
      printf("\n");
    }

    Individual *test_ind;
    test_ind = new Individual;
    test_ind->duplicate(seed_individ_);
    test_ind->write("duplicate.ind");

    return true;
  }


  // Run individual.
  if (str_eq("run", config_line)) {
    loc = find_field(config_line);
    i = loc;
    while (config_line[i] != '\n') {
      i++;
      if (i > Max_Line_Length) {
	printf("AlpsGen :: parsing - run, max line length exceeded.\n");
	return false;
      }
    }

    /*
    config_line[i] = 0;
    if (verbose) {
      printf("Run individual: %s\n", &(config_line[loc]));
    }
    */

    if (!read_sample_individ(&(config_line[loc]))) {
      printf("AlpsGen :: error reading sample individual.\n");
      return false;
    }

    sample_individ_->print_history_full();
    sample_individ_->write("run.ind");
    run_sample_ = true;

    return true;
  }


  // Restart population
  if (str_eq("restart pop", config_line)) {
    loc = find_field(config_line);
    i = loc;
    while (config_line[i] != '\n') {
      i++;
      if (i > Max_Line_Length) {
	printf("AlpsGen :: parsing - restart pop, max line length exceeded.\n");
	return false;
      }
    }

    /*
    config_line[i] = 0;
    if (verbose) {
      printf("Restarting population: %d:%d  %s\n",
	     loc, i, &(config_line[loc]));
    }
    */

    if (!read(&(config_line[loc]))) {
      return false;
    }

    write("sample.pop");
    if (verbose) {
      print_state();
      printf("\n");
    }

    return true;
  }
#endif


  if (Alps::parse_config(config_line, verbose)) {
    return true;
  }

  printf("AlpsGen :: Error with: %s\n", config_line);
  return false; // Not a recognized line type.
}



/********************************************************/


// Returns:
//   1 if index1 better fitness than index2
//   0 if equal
//   -1 if index1 strictly worse than index2
int AlpsGen :: compare_childs(int index1, int index2)
{
  Individual *ind1, *ind2;

  ind1 = pop_childs_[index1];
  ind2 = pop_childs_[index2];

  /*
  printf("alps :: compare_childs() : %d(%d) vs %d(%d) = %d.\n",
	 index1, (int)ind1->get_fitness(),
	 index2, (int)ind2->get_fitness(),
	 ind2->compare_fitness(ind1));
  */

  return ind1->compare_fitness(is_maximizing_, ind2);
}


AlpsGen *current_pop;
int compare_childs_func(const void *elt1, const void *elt2)
{
  int index1, index2;

  index1 = *(int*)elt1;
  index2 = *(int*)elt2;
  return current_pop->compare_childs(index1, index2);
}



void AlpsGen :: sort_layer_childs(int layer)
{
  int start_index = layer_def_[layer].get_start_index();
  int stop_index = layer_def_[layer].get_start_index() 
    + layer_def_[layer].get_size();


  /*
  vector<Individual*> inds(layer_def_[layer].get_size());
  for (unsigned int i=0; i<layer_def_[layer].get_size(); i++) {
    inds[i] = pop_childs_[i+start_index];
    if (inds[i] == 0) {
      cout << "Error, null individ! : " << i << endl;
      while (1) ;
    }
  }

  vector<Individual*>::iterator lstart = inds.begin();
  vector<Individual*>::iterator lstop = inds.end();

  if (is_maximizing_) {
    // !!!! Regular sort() has a bug in it.
    stable_sort(lstart, lstop, individ_compare_max);
  } else {
    // !!!!! Error in here when all optimization flags turned on.
    stable_sort(lstart, lstop, individ_compare_min);
  }

  for (unsigned int i=0; i<layer_def_[layer].get_size(); i++) {
    pop_childs_[i+start_index] = inds[i];
  }
  */

  vector<Individual*>::iterator lstart = pop_childs_.begin() 
    + start_index;
  vector<Individual*>::iterator lstop = pop_childs_.begin()
    + stop_index;

  // For some reason regular sort() crashes but stable_sort() works.
  if (is_maximizing_) {
    stable_sort(lstart, lstop, individ_compare_max);
  } else {
    stable_sort(lstart, lstop, individ_compare_min);
  }

  /*
  cout << "alps :: sort_layer_childs() layer " << layer << ", from "
       << start_index << " to " << stop_index << ".\n";
  */  
}




// Returns:
//   1 if index1 better fitness than index2
//   0 if equal
//   -1 if index1 strictly worse than index2
int AlpsGen :: compare_parents(int index1, int index2)
{
  Individual *ind1, *ind2;

  ind1 = pop_parents_[index1];
  ind2 = pop_parents_[index2];

  /*
  printf("alps :: compare_parents() : %d(%d) vs %d(%d) = %d.\n",
	 index1, (int)ind1->get_fitness(),
	 index2, (int)ind2->get_fitness(),
	 ind2->compare_fitness(ind1));
  */

  return ind1->compare_fitness(is_maximizing_, ind2);
  //  return ind2->compare_fitness(ind1);
}




void AlpsGen :: sort_parents(vector<Individual*>& inds, bool sort_is_maximizing_)
{
  vector<Individual*>::iterator lstart = inds.begin();
  vector<Individual*>::iterator lstop = inds.end();

  if (sort_is_maximizing_) {
    // !!!! Regular "sort()" has a bug in it.
    stable_sort(lstart, lstop, individ_compare_max);
  } else {
    stable_sort(lstart, lstop, individ_compare_min);
  }

  return;
}



// start_index will have worst fitness.
// stop_index-1 will have best fitness.
void AlpsGen :: sort_layer_parents(int layer)
{
  int start_index = layer_def_[layer].get_start_index();
  int stop_index = layer_def_[layer].get_start_index() + layer_def_[layer].get_size();

  vector<Individual*>::iterator lstart = pop_parents_.begin() 
    + start_index;
  vector<Individual*>::iterator lstop = pop_parents_.begin()
    + stop_index;

  if (is_maximizing_) {
    // !!!! Regular "sort()" has a bug in it.
    stable_sort(lstart, lstop, individ_compare_max);
  } else {
    stable_sort(lstart, lstop, individ_compare_min);
  }

  /*
  cout << "alps :: sort_layer_parents(" << layer << "), max:"
       << is_maximizing_ << ", layer " << layer
       << ", from " << start_index << " to " << stop_index << ".\n";
  */

  return;
}



// Sort the population, pop_parents_, from highest, pop_order_[0],
// to lowest, pop_order_[num_individs_-1].
void AlpsGen :: order_population()
{
  int i;

  // Order each layer.
  int max_layer = calc_max_layer(generation_);

  cout << "alps :: order_population() - max_layer: " << max_layer;

  for (i=0; i<=max_layer; i++) {
    sort_layer_parents(i);
  }

  cout << "alps :: layers sorted:\n";
  for (i=0; i<=max_layer; i++) {
    for (unsigned int j=layer_def_[i].get_start_index(); j<layer_def_[i].get_start_index()+layer_def_[i].get_size(); j++) {
      cout << j << " : " << create_type_list_[j];
      if (create_type_list_[j] != CREATE_NOT) {
	cout << "  " << pop_parents_[j]->get_fitness() << "\n";
      } else {
	cout << "  -\n";
      }
    }
  }
}


bool AlpsGen :: make_individual(int index)
{
  if (generation_ > max_generations_)
    return false;

  if (create_type_list_[index] == CREATE_RECOMBINATION) {
    pop_childs_[index]->duplicate(pop_parents_[parent1_list_[index]]);
    //    pop_childs_[index]->recombine(pop_parents_[parent2_list_[index]]);
    pop_childs_[index]->recombine(tourn_list_[index]);

    if (print_debug_info_) {
      cout << "alps :: make_individual()  "
	   << parent1_list_[index] << " " << parent2_list_[index]
	   << " -> " << index << "\n";
    }

  } else if (create_type_list_[index] == CREATE_MUTATION) {
    pop_childs_[index]->duplicate(pop_parents_[parent1_list_[index]]);
    pop_childs_[index]->mutate(tourn_list_[index]);
    if (print_debug_info_) {
      cout << "alps :: make_individual()  "
	   << parent1_list_[index] << " -> " << index << "\n";
    }

  } else if (create_type_list_[index] == CREATE_DUPLICATE) {
    pop_childs_[index]->duplicate(pop_parents_[parent1_list_[index]]);

  } else if (create_type_list_[index] == CREATE_RANDOM) {
    pop_childs_[index]->make_random();
    if (print_debug_info_) {
      cout << "alps :: make_individual(" << index 
	   << ") - create_random.\n";
    }

  } else {
    cout << "alps :: make_individual(" << index
	 << ") - invalid make type: " << create_type_list_[index] << "\n";
    while (1) ;
  }

  return true;
}



int AlpsGen :: random_parent()
{
  cout << "random_parent() - not implemented.\n";
  while (1) ;
  return 0;

#ifdef OLD
  int bot, mid, top;
  double val;

  bot = 0;
  top = num_individs_-1;
  mid = top/2;

  val = random_double() * Selection_Prob[num_individs_-1];

  //  printf("val= %f  ", val);


  while (bot != top) {
    mid = (top-bot)/2 + bot;
    if (Selection_Prob[mid] < val) {
      if (bot == mid)
	return top;
      bot = mid;
	continue;

    } else if (mid == bot) {
      return mid;

    } else if (Selection_Prob[mid-1] < val) {
      return mid;

    } else {
      top = mid;
      continue;
    }
  }

  return bot; // not used.
#endif
}


//int AlpsGen :: random_parent2(int index1)
int AlpsGen :: random_parent2(int)
{
  cout << "random_parent() - not implemented.\n";
  while (1) ;
  return 0;

#ifdef OLD
  int bot, mid, top, index=0;
  double val;

  bot = 0;
  top = num_individs_-2;
  mid = top/2;

  val = random_double() * (Selection_Prob[num_individs_-2]);
  while (bot != top) {
    mid = (top-bot)/2 + bot;
    if (Selection_Prob[mid] < val) {
      if (bot == mid) {
	index = top;
	break;
      }

      bot = mid;
      continue;

    } else if (mid == bot) {
      index = mid;
      break;

    } else if (Selection_Prob[mid-1] < val) {
      index =  mid;
      break;

    } else {
      top = mid;
      continue;
    }
  }
  
  if (top == bot)
    index = bot;

  if (index >= index1)
    index++;

  return index;
#endif
}



int AlpsGen :: calc_max_layer(int gen)
{
  if (num_layers_ == 1) {
    return 0;
  }

  if (gen >= layer_def_[num_layers_-2].get_max_age()) {
    return num_layers_-1;
  }

  for (int i=num_layers_-2; i>=0; i--) {
    if (gen >= layer_def_[i].get_max_age()) {
      for (int j=i+2; j<num_layers_; j++) {
	if (layer_def_[j].get_max_age() > layer_def_[i+1].get_max_age()) {
	  return j-1;
	}
      }
      return i+1;
    }
  }

  for (int j=1; j<num_layers_; j++) {
    if (layer_def_[j].get_max_age() > layer_def_[0].get_max_age()) {
      return j-1;
    }
  }

  return 0;
}


// calc_max_layer(1) is prob better.
int AlpsGen :: calc_max_layer()
{
  if (num_layers_ == 1) {
    return 0;
  }

  for (int i=num_layers_-1; i >= 0; i--) {
    if (layer_def_[i].get_max_age() != -1) {
      if (i == num_layers_-1) {
	printf("alps :: calc_max_layer() - error, last layer does not have an age of -1: %d.  generation:%d.\n",
	       layer_def_[i].get_max_age(), generation_);
	while (1) ;
      }
      if (layer_def_[i].get_max_age() < generation_) {
	return num_layers_-1;
      }
    }
  }


  if (generation_ > layer_def_[num_layers_-2].get_max_age()) {
    return num_layers_-1;
  }

  for (int i=0; i<num_layers_; i++) {
    if (layer_def_[i].get_max_age() > generation_) {
      for (int j=i+1; j<num_layers_; j++) {
	if (layer_def_[j].get_max_age() > layer_def_[i].get_max_age()) {
	  return j-1;
	}
      }
      return i;
    }
  }

  for (int i=num_layers_-1; i >= 0; i--) {
    if (layer_def_[i].get_max_age() != -1) {
      if (i == num_layers_-1) {
	printf("alps :: calc_max_layer(1) - error, last layer does not have an age of -1: %d.  generation:%d.\n",
	       layer_def_[i].get_max_age(), generation_);
	while (1) ;
      }
      if (layer_def_[i].get_max_age() < generation_) {
	return num_layers_-1;
      }
    }
  }

  printf("alps :: calc_max_layer(1) - error, reached end of function without calculating max_layer. Gen=%d.\n",
	 generation_);
  while (1) ;
}


void AlpsGen :: setup_layer_tourn(int layer)
{
  int start_index = layer_def_[layer].get_start_index();
  int stop_index = layer_def_[layer].get_start_index() + layer_def_[layer].get_size();
  int max_age = layer_def_[layer].get_max_age();


  if (print_debug_info_) {
    cout << "alps :: setup_layer_tourn(): layer " << layer
	 << " (" << start_index << "-" << stop_index << "), elitism:"
	 << layer_def_[layer].get_elitism() << "\n";
  }

  sort_layer_parents(layer);

  int num_selectable_in_layer = 0;
  vector<int> parent_list1(layer_def_[layer].get_size());
  for (int i=start_index; i<stop_index; i++) {
    if (create_type_list_[i] != CREATE_NOT) {
      create_type_list_[i] = CREATE_WAITING;
      if ((pop_parents_[i]->get_age() < max_age) ||
	  (max_age == -1)) {
	parent_list1[num_selectable_in_layer++] = i;
      }
    }
  }


  // Clear the tourn_list_. (actually a vector).
  for (int i=start_index; i<stop_index; i++) {
    tourn_list_[i].resize(0);
  }


  // Do Elitism.
  if (layer_def_[layer].get_elitism() > 0) {
    if (moveup_old_elite_) {
      // Here, if an "elite" individual is too old for its layer
      // it will be moved up to the next layer in the post process phase.
      for (int i = stop_index - layer_def_[layer].get_elitism();
	   i<stop_index; i++) {
	if ((pop_parents_[i]->get_age() < max_age) ||
	    (max_age < 0)) {
	  create_type_list_[i] = CREATE_DUPLICATE;
	  parent1_list_[i] = i;
	} else {
	  create_type_list_[i] = CREATE_WAITING;
	}
      }

      for (int i = stop_index - layer_def_[layer].get_elitism();
	   i<stop_index; i++) {
	if ((pop_parents_[i]->get_age() >= max_age) &&
	    (max_age >= 0) && (create_type_list_[i] == CREATE_DUPLICATE)) {
	  cout << "error, too old for elitism.\n";
	  while (1) ;
	}
      }

    } else {
      // The top elite individuals are kept but will not be moved up.
      unsigned int num_elite = 0;
      int last_elite = stop_index-1;

      vector<int> elite_index; // Individuals picked, to be moved.

      for (int i=stop_index-1; i>= start_index; i--) {
	if (create_type_list_[i] == CREATE_NOT) {
	  // No more individuals exist in this layer.
	  break;
	}
    
	if ((pop_parents_[i]->get_age() > max_age) &&
	    (max_age >= 0)) {
	  // Too old to keep in this layer.
	  continue;
	}

	num_elite++;
	create_type_list_[i] = CREATE_DUPLICATE;
	parent1_list_[i] = i;
	if (i < stop_index - (int)layer_def_[layer].get_elitism()) {
	  // An individ to move up.
	  elite_index.push_back(i);
	}
	last_elite = i;
	//    cout << "  " << i << " : " << pop_parents_[i]->get_fitness() << "\n";

	if (num_elite == layer_def_[layer].get_elitism()) {
	  break;
	}
      }

      int ctr = stop_index - 1;
      for (vector<int>::size_type i=0; i<elite_index.size(); i++) {
	while (create_type_list_[ctr] != CREATE_DUPLICATE) {
	  ctr--;
	}

	Individual* tmp = pop_parents_[ctr];
	pop_parents_[ctr] = pop_parents_[elite_index[i]];
	pop_parents_[elite_index[i]] = tmp;

	create_type_list_[elite_index[i]] = create_type_list_[ctr];
	create_type_list_[ctr] = CREATE_DUPLICATE;
	create_type_list_[ctr] = ctr;
      }
    }

    if (print_debug_info_) {
      cout << " L" << layer << " elite[" << max_age << "]: ";
      int start = layer_def_[layer].get_start_index();
      int stop = layer_def_[layer].get_start_index() + layer_def_[layer].get_size() - 1;

      for (int i=stop; i>=start; i--) {
	if (create_type_list_[i] == CREATE_DUPLICATE) {
	  cout << "  " << i << ":" << pop_parents_[i]->get_fitness()
	       << "(" << pop_parents_[i]->get_age() << ")";
	}
      }
      cout << "\n";
    }
  }

  //  cout << " -> num_elite=" << num_elite << "\n";


  int tourn_size = layer_def_[layer].get_tourn_size();
  for (int i=start_index; i<stop_index; i++) {
    if (create_type_list_[i] == CREATE_DUPLICATE) {
      continue;
    }

    int j=0;
    int k=0;
    vector<int> tournament_list(tourn_size);

    if ((pop_parents_[i]->get_age() <= max_age) ||
	(max_age == -1)) {
      tournament_list[0] = i;
      j = 1;
    }

    for (; j<tourn_size; j++) {
      double rval = random_double();
      if ((num_selectable_in_layer > j)
	  && (rval >= layer_def_[layer].get_prob_select_prev())) {
	// Select a potential parent from this layer.
	do {
	  if (num_selectable_in_layer > 1) {
	    tournament_list[j] = parent_list1[random_int(num_selectable_in_layer)];
	  } else {
	    tournament_list[j] = parent_list1[0];
	  }

	  if ((pop_parents_[tournament_list[j]]->get_age() > max_age) &&
	      (max_age >= 0)) {
	    printf("alps :: setup_layer_tourn() - error, individuals age %d, is greater than max_age: %d.\n",
		   pop_parents_[tournament_list[j]]->get_age(),
		   max_age);
	    while (1) ;
	  }


	  k = 0;
	  for (; k<j; k++) {
	    if (tournament_list[j] == tournament_list[k]) {
	      break;
	    }
	  }
	} while (k < j);

      } else if (layer_def_[layer].get_num_prev_layers() == 0) {
	int index = layer_def_[layer].get_start_index();
	if (layer_def_[layer].get_size() > 1) {
	  index += random_int(layer_def_[layer].get_size());
	}
	tournament_list[j] = index;

      } else {
	// Select a potential parent from a "previous" layer.
	int prev_layer = -1;
	if (layer_def_[layer].get_num_prev_layers() == 1) {
	  prev_layer = layer_def_[layer].get_prev_layer(0);
	} else {
	  prev_layer = layer_def_[layer].get_prev_layer(random_int(layer_def_[layer].get_num_prev_layers()));
	}

	do {
	  /*
	  do {
	    tournament_list[j] = random_int(layer_def_[prev_layer].get_size()) + 
	      layer_def_[prev_layer].get_start_index();
	  } while ((pop_parents_[tournament_list[j]]->get_age() > max_age) &&
		   (max_age >= 0));
	  */
	  tournament_list[j] = random_int(layer_def_[prev_layer].get_size()) + 
	    layer_def_[prev_layer].get_start_index();

	  k = 0;
	  for (; k<j; k++) {
	    if (tournament_list[k] == tournament_list[j]) {
	      break;
	    }
	  }
	} while (k < j);
      }
    }

    int parent1=-1, parent2=-1;

    if (random_double() < prob_recomb_) {
      create_type_list_[i] = CREATE_RECOMBINATION;
      
      // Tournament selection:
      // Have tournament of size tourn_size.
      // Take the two best individuals for recombination.
      if (compare_parents(tournament_list[0], tournament_list[1]) != -1) {
	parent1 = tournament_list[0];
	parent2 = tournament_list[1];

      } else {
	parent1 = tournament_list[1];
	parent2 = tournament_list[0];
      }

      for (int j=2; j<tourn_size; j++) {
	if (compare_parents(tournament_list[j], parent1) == 1) {
	  parent2 = parent1;
	  parent1 = tournament_list[j];

	} else if (compare_parents(tournament_list[j], parent2) == 1) {
	  parent2 = tournament_list[j];
	}
      }
      parent1_list_[i] = parent1;
      parent2_list_[i] = parent2;
      if (random_double() < prob_rec_rand2_) {
	if (tourn_size > 2) {
	  parent2_list_[i] = tournament_list[random_int(tourn_size-1)];
	  if (parent2_list_[i] == parent1_list_[i]) {
	    parent2_list_[i] = tournament_list[tourn_size-1];
	  }
	}
      }

      if (compare_parents(parent1, parent2) == -1) {
	cout << " => Error, parents not sorted.\n";
	while (1) ;
      }

    } else {
      create_type_list_[i] = CREATE_MUTATION;

      parent1 = tournament_list[0];
      for (int j=1; j<tourn_size; j++) {
	if (compare_parents(tournament_list[j], parent1) == 1) {
	  parent1 = tournament_list[j];
	}
      }
      parent1_list_[i] = parent1;
    }

    // Save the individuals in the tournament to pass to the
    // variation operator.
    vector<int> t_indices;
    t_indices.resize(0);
    for (int k=0; k<tourn_size; k++) {
      if (tournament_list[k] != parent1_list_[i]) {
	t_indices.push_back(tournament_list[k]);
	tourn_list_[i].push_back(pop_parents_[tournament_list[k]]);
      }
    }
    sort_parents(tourn_list_[i], !is_maximizing_);

    /*
    cout << i << " TL: " << parent1_list_[i] << ":"
	 << pop_parents_[parent1_list_[i]]->get_fitness();
    if (create_type_list_[i] == CREATE_RECOMBINATION) {
      cout << " P2:" << parent2_list_[i] << ":"
	   << pop_parents_[parent2_list_[i]]->get_fitness();
      if (compare_parents(parent1_list_[i], parent2_list_[i]) == -1) {
	cout << endl;
	cout << " !!! - - Error, parent2 is more fit.\n";
	while (1) ;
      }

    } else {
      cout << " mut               ";
    }

    for (unsigned int k=0; k<tourn_list_[i].size(); k++) {
      cout << "  " << t_indices[k] << ":"
	   << tourn_list_[i][k]->get_fitness();
    }
    cout << endl;
    */
  }
}



void AlpsGen :: setup_layer_roulette(int layer)
{
  int start_index = layer_def_[layer].get_start_index();
  int stop_index = layer_def_[layer].get_start_index() + layer_def_[layer].get_size();

  for (unsigned int i=0; i<layer_def_[layer].get_elitism(); i++) {
    create_type_list_[i+start_index] = CREATE_DUPLICATE;
    parent1_list_[i+start_index] = pop_order_[i+start_index];
  }

  //// The following code is only starter code and needs to be re-written.

  int parent1 = -1, parent2 = -1;
  for (int i=layer_def_[layer].get_start_index()+layer_def_[layer].get_elitism();
       i<stop_index; i++) {
    if (random_double() < prob_recomb_) {
      create_type_list_[i] = CREATE_RECOMBINATION;
      parent1 = random_parent();
      parent2 = random_parent2(parent1);

      parent1_list_[i] = pop_order_[parent1];
      parent2_list_[i] = pop_order_[parent2];
      /*
	printf("Parents: %2d(%2d)  %2d(%2d) -> %2d\n",
	parent1_list_[i], parent1,
	parent2_list_[i], parent2, i);
      */

    } else {
      create_type_list_[i] = CREATE_MUTATION;
      parent1 = random_parent();
      parent1_list_[i] = pop_order_[parent1];
      /*
	printf("Parent:  %2d(%2d) -> %2d\n",
	parent1_list_[i], parent1, i);
      */
    }
  }
}



void AlpsGen :: setup_layer_dc(int layer)
{
  int start_index = layer_def_[layer].get_start_index();
  int stop_index = layer_def_[layer].get_start_index() + layer_def_[layer].get_size();
  int max_age = layer_def_[layer].get_max_age();

  int prev_layer = -1;
  if (layer_def_[layer].get_num_prev_layers() == 1) {
    prev_layer = layer_def_[layer].get_prev_layer(0);

  } else if (layer_def_[layer].get_num_prev_layers() > 1) {
    prev_layer = layer_def_[layer].get_prev_layer(random_int(layer_def_[layer].get_num_prev_layers()));
  }

  int prevl_start = -1;
  int prevl_stop = -1;
  if (prev_layer != -1) {
    prevl_start = (int)layer_def_[prev_layer].get_start_index();
    prevl_stop = (int)layer_def_[prev_layer].get_stop_index();
  }


  // Clear the tourn_list_. (actually a vector).
  for (int i=start_index; i<stop_index; i++) {
    tourn_list_[i].resize(0);
  }


  vector<int> potential_parents;
  for (int i=start_index; i<stop_index; i++) {
    if (create_type_list_[i] != CREATE_NOT) {
      create_type_list_[i] = CREATE_WAITING;
      if ((pop_parents_[i]->get_age() <= max_age) ||
	  (max_age == -1)) {
	potential_parents.push_back(i);
      }
    }
  }




  int tourn_size = layer_def_[layer].get_tourn_size();
  if (tourn_size < 1) {
    tourn_size = 1;
  }

		       
  for (int i=start_index; i<stop_index; i++) {
    // Decide on the variation operator.
    if (random_double() < prob_recomb_) {
      create_type_list_[i] = CREATE_RECOMBINATION;
    } else {
      create_type_list_[i] = CREATE_MUTATION;
    }


    // Set Parent 1.
    if ((pop_parents_[i]->get_age() <= max_age) ||
	(max_age == -1)) {
      // Use the current individual.
      parent1_list_[i] = i;

    } else {
      // Current individ is too old, pick one at random.
      if (potential_parents.size() > 0) {
	parent1_list_[i] = potential_parents[0];
	if (potential_parents.size() > 1) {
	  parent1_list_[i] = potential_parents[random_int(potential_parents.size())];
	}

      } else if (prevl_start >= prevl_stop) {
	parent1_list_[i] = i;
	
      } else {
	parent1_list_[i] = prevl_start;
	int range = prevl_stop - prevl_start;
	if (range > 1) {
	  parent1_list_[i] = random_int(range)+prevl_start;
	}
      }
    }

    // Now, do tournament for parent2.
    int k=0;
    vector<int> tournament_list(tourn_size);

    for (int j=0; j<tourn_size; j++) {
      if ((((int)potential_parents.size() > j+1)
	   && (random_double() >= layer_def_[layer].get_prob_select_prev()))
	  || (layer_def_[layer].get_num_prev_layers() == 0)) {
	// Select a potential parent from this layer.
	do {
	  tournament_list[j] = potential_parents[0];
	  if (potential_parents.size() == 1) {
	    break;
	  }
	  
	  tournament_list[j] = potential_parents[random_int(potential_parents.size())];	
	  if (parent1_list_[i] == tournament_list[j]) {
	    continue;
	  }

	  k = 0;
	  for (; k<j; k++) {
	    if (tournament_list[j] == tournament_list[k]) {
	      break;
	    }
	  }
	} while (k < j);

      } else {
	// Select a potential parent from a "previous" layer.
	int prev_layer = -1;
	if (layer_def_[layer].get_num_prev_layers() == 1) {
	  prev_layer = layer_def_[layer].get_prev_layer(0);
	} else {
	  prev_layer = layer_def_[layer].get_prev_layer(random_int(layer_def_[layer].get_num_prev_layers()));
	}

	int num_tries = 2 *tourn_size * layer_def_[prev_layer].get_size();
	do {
	  do {
	    tournament_list[j] = random_int(layer_def_[prev_layer].get_size()) + 
	      layer_def_[prev_layer].get_start_index();
	  } while ((pop_parents_[tournament_list[j]]->get_age() > max_age) &&
		   (max_age >= 0));

	  k = 0;
	  for (; k<j; k++) {
	    if (tournament_list[k] == tournament_list[j]) {
	      break;
	    }
	  }
	} while ((k < j) && (num_tries-- > 0));
      }
    }


    // Save the individuals in the tournament to pass to the
    // variation operator.
    vector<int> t_indices;
    t_indices.resize(0);
    for (int k=0; k<tourn_size; k++) {
      if (tournament_list[k] != parent1_list_[i]) {
	t_indices.push_back(tournament_list[k]);
	tourn_list_[i].push_back(pop_parents_[tournament_list[k]]);
      }
    }
    sort_parents(tourn_list_[i], !is_maximizing_);

    if (print_debug_info_) {
      cout << "tourn_list L" << layer << " (" << tournament_list.size() << "): ";
      for (vector<int>::size_type k=0; k<tournament_list.size(); k++) {
	cout << " " << tournament_list[k];
      }
      cout << endl;
    }
  }
}



void AlpsGen :: setup_layer_random(int layer)
{
  int start_index = layer_def_[layer].get_start_index();
  int stop_index = layer_def_[layer].get_stop_index();

  for (int i=start_index; i<stop_index; i++) {
    create_type_list_[i] = CREATE_RANDOM;
  }
}




int AlpsGen :: setup_next_generation_layers()
{
  // Setup the layers.
  //  int max_layer = calc_max_layer(generation_);
  //  int max_layer2 = calc_max_layer();

  /*
  static int count = 0;
  printf("alps :: setup_next_generation_layers(): gen:%2d, maxlayer:%2d (%d).\n",
	 generation_, max_layer_evolving_, count);
  count++;
  if (count > 3) {
    while (1) ;
  }
  */

  /*
  static int count = 0;
  count++;
  cout << "alps :: setup_next_generation_layers(): gen:"
       << generation_ << ", maxlayer:" << max_layer_evolving_
       << " (" << count << ").\n";
  print_pop();
  */

  for (int i=0; i<=max_layer_evolving_; i++) {
    if (print_debug_info_) {
      printf("alps :: sngl(): setup layer:%2d, num-prev:%2d, max-age:%2d, mod:%2d.\n",
	     i, layer_def_[i].get_num_prev_layers(), layer_def_[i].get_max_age(),
	     (generation_ % layer_def_[i].get_max_age() == 0));	     
    }

    if (((layer_def_[i].get_num_prev_layers() == 0) &&
	 (layer_def_[i].get_max_age() > 0) &&
	 (generation_ % layer_def_[i].get_max_age() == 0)) ||
	(generation_ == 0)) {
      //      cout << "alps_gen :: sngl() - layer " << i << " setup random." << endl;
      setup_layer_random(i);

    } else if (layer_def_[i].get_select_type() == ALPS_SELECT_TOURN) {
      setup_layer_tourn(i);

    } else if (layer_def_[i].get_select_type() == ALPS_SELECT_ROULETTE) {
      setup_layer_roulette(i);

    } else if (layer_def_[i].get_select_type() == ALPS_SELECT_DC) {
      setup_layer_dc(i);

    } else {
      printf("alps :: setup_next_generation_layers() - error, layer %d, invalid selection type: %d.\n",
	     i, layer_def_[i].get_select_type());
      while (1) ;
    }
  }

  if (print_debug_info_) {
    print_parents();
  }

  // initialize variables.
  next_to_evaluate_ = 0;
  num_to_evaluate_ = 0;
  for (int i=0; i<=max_layer_evolving_; i++) {
    num_to_evaluate_ += layer_def_[i].get_size();
  }

  int i=0;
  for (; i<num_to_evaluate_; i++) {
    evaluation_list_[i] = i;
    tries_list_[i] = 0;
    evaluation_state_[i] = ES_NOT_CREATED;
  }
  for (; i<num_individs_; i++) {
    evaluation_list_[i] = -1;
    tries_list_[i] = 0;
    evaluation_state_[i] = ES_NOT_USED;
  }
  num_being_evaluated_ = 0;


  vector<int> increase_age(num_individs_);
  for (int i=0; i<num_individs_; i++) {
    increase_age[i] = false;
  }

  for (unsigned int i=0; i<layer_def_[max_layer_evolving_].get_start_index()+
	 layer_def_[max_layer_evolving_].get_size(); i++) {
    if (create_type_list_[i] == CREATE_RECOMBINATION) {
      increase_age[parent1_list_[i]] = true;
      increase_age[parent2_list_[i]] = true;
    } else if (create_type_list_[i] == CREATE_MUTATION) {
      increase_age[parent1_list_[i]] = true;
    }
  }

  for (unsigned int i=0; i<layer_def_[max_layer_evolving_].get_start_index()+
	 layer_def_[max_layer_evolving_].get_size(); i++) {
    /*
    if (increase_age[i]) {
      pop_parents_[i]->increase_age();
    }
    */
    pop_parents_[i]->increase_age();
  }

  /*
  printf("#2Eval:%d, Next2Eval:%d, #BingEval:%d\n",
	 num_to_evaluate_, next_to_evaluate_, num_being_evaluated_);
  */
  //  print_pop();
  //  print_parents();
  //  print_evaluation_state();

  /*
  if (save_log_) {
    ofstream individ_log(Individ_log_file_Name.c_str(), ios_base::app);
    if (!individ_log) {
      cerr << "alps :: setup_next_generation() - error opening individual history log.\n";
    }

    for (int i=0; i<=max_layer_evolving_; i++) {
      int index = layer_def_[i].get_start_index() + layer_def_[i].get_size() - 1;
      if ((pop_parents_[index]->get_num_evaluations() == 1) &&
	  (pop_parents_[index]->get_fitness() > 0.0)) {
	pop_parents_[index]->write_log(individ_log);
      }
    }

    individ_log.close();
  }
  */


  if (update_correlation_func != 0) {
    for (int i=0; i<=max_layer_evolving_; i++) {
      int index = layer_def_[i].get_start_index() + layer_def_[i].get_size() - 1;
      if (pop_parents_[index]->get_num_evaluations() == 1) {
	update_correlation_func(pop_parents_[index]);
      }
      //      update_correlation_func(index);
    }
  }

  return true;
}



// Checks to see if individual at index in pop_parents_ can be
//  moved into Layer layer in pop_childs.
// Assume layer is sorted, in ascending order.
// Inserts individual in correct spot, leaving the layer sorted.
bool AlpsGen :: move_up(int index, int layer)
{
  if (layer >= num_layers_) {
    /*
    printf("alps :: move_up() - index:%d, error, exceeded max layer: %d > %d.\n",
	   index, layer, num_layers_);
    */
    return false;
  }

  int start_index = layer_def_[layer].get_start_index();
  int stop_index = layer_def_[layer].get_start_index() + layer_def_[layer].get_size();

  /*
  printf("alps :: move_up() - index:%d (fit:%d) to layer:%d (%d - %d).\n",
	 index, (int)pop_parents_[index]->get_fitness(), layer,
	 start_index, stop_index);
  */

  // Does an insertion sort to insert the individual at index into this layer.
  int i=start_index;
  int can_replace = -1;
  for (; i<stop_index; i++) {
    /*
    if ((create_type_list_[i] != CREATE_NOT) &&
	(compare_ch
    */					       

    if (create_type_list_[i] == CREATE_NOT) {
      can_replace = i;
      start_index = i;
      //      printf(" : index %d is not created.\n", i);

    } else if (pop_parents_[index]->compare_fitness(is_maximizing_,
						   pop_childs_[i]) >= 0) {
      can_replace = i;
      /*
      printf(" : index %d (fit:%d) can replace since its less fit: compare(%d,%d)=%d.\n", i,
	     (int)pop_childs_[i]->get_fitness(),
	     index, i, compare_childs(index, i));
      */

    } else {
      // All individuals past here are "higher".
      /*
      printf(" : index %d (fit:%d) is more fit.\n", i,
	     (int)pop_childs_[i]->get_fitness());
      */
      break;
    }
  }

  if (can_replace == -1) {
    //    printf(" - can't replace.  exiting.\n");
    return false;
  }


  /*
  printf(" - can replace = %d (%d).\n", can_replace,
	 (int)pop_childs_[can_replace]->get_fitness());
  */

  if (create_type_list_[start_index] != CREATE_NOT) {
    // Try to move this individual up.
    move_up(start_index, layer+1);
  }

  Individual* tmp_child = pop_childs_[start_index];
  Individual* tmp_parent = pop_parents_[start_index];
  create_type_list_[start_index] = CREATE_WAITING;

  for (; start_index < can_replace; start_index++) {
    pop_childs_[start_index] = pop_childs_[start_index+1];
    pop_parents_[start_index] = pop_parents_[start_index+1];
  }

  /*
  printf(" - inserting %d:%d at index: %d.\n", index,
	 (int)pop_parents_[index]->get_fitness(), can_replace);
  */
  pop_childs_[can_replace] = tmp_child;
  pop_childs_[can_replace]->duplicate(pop_parents_[index]);
  pop_parents_[can_replace] = tmp_parent;
  pop_parents_[can_replace]->duplicate(pop_parents_[index]);
  create_type_list_[can_replace] = CREATE_WAITING;

  return true;
}


// Try to move up as many of individs from [l1_start, l1_end)
//  to the next layer.
void AlpsGen :: try_move_up_c2c(int l1_start, int l1_end,
				int layer)
{
  int l1_index = l1_end-1;

  int l2_start = layer_def_[layer].get_start_index();
  int l2_end = l2_start + layer_def_[layer].get_size();
  int l2_index = l2_start;


  int num_move = 0;
  do {
    if (pop_childs_[l1_index]->compare_fitness(is_maximizing_,
					      pop_childs_[l2_index]) == -1) {
      // Can't replace, done.
      break;
    }
    l1_index--;
    l2_index++;
    num_move++;
  } while ((l1_index >= l1_start) && (l2_index < l2_end));

  /*  
  cout << "try_move_up_c2c(" << layer << ") - " << l1_index << "("
       << pop_childs_[l1_index]->get_fitness() << ") - " << l1_end-1
       << "(" << pop_childs_[l1_end-1]->get_fitness() << "), to L"
       << layer << "\n";
  */

  if (num_move == 0) {
    //    cout << " - nothing worse to replace.\n";
    return;
  }

  /*
  cout << " - replace: " << l2_start << "("
       << pop_childs_[l2_start]->get_fitness() << ") - " << l2_index-1
       << "(" << pop_childs_[l2_index-1]->get_fitness() << ")\n";
  */

  int size = (l1_end-1) - l1_index;
  if (size != num_move) {
    cout << "alps_gen :: try_move_up() - size:" << size << " != num:"
	 << num_move << " WARNING!! " << endl;
    while (1) ;
  }


  int lparent = layer_def_[layer].get_parent_layer();
  if ((l2_start < l2_index) && (lparent != -1)) {
    try_move_up_c2c(l2_start, l2_index, lparent);
  }

  int l1_offset = l1_end - num_move;
  for (int i=0; i<num_move; i++) {
    Individual* tmp = pop_childs_[l2_start+i];
    pop_childs_[l2_start+i] = pop_childs_[l1_offset + i];
    pop_childs_[l1_offset + i] = tmp;

    if (l1_offset + i >= l1_end) {
      cout << "alps_gen :: try_move_up() - error, l1 index too high!: "
	   << l1_offset + i << " >= " << l1_end << endl;
      while (1) ;
    }
  }

  sort_layer_childs(layer);
}


// Copy individuals from Pop_Parent[l1_start, l1_end)
// to pop_childs_[], starting from the beginning.
// But also try to move up those individs in pop_childs_[].
void AlpsGen :: move_up_p2c(int l1_start, int l1_end, int layer)
{
  cout << "move_up_p2c(int, int, int) - not implemented.\n";
  while (1) ;

  int l2_start = layer_def_[layer].get_start_index();
  int l2_end = l2_start + (l1_end - l1_start);
  if (l2_end > l2_start + (int)layer_def_[layer].get_size()) {
    l2_end = l2_start + layer_def_[layer].get_size();
    l1_start = l1_end - (l2_end - l2_start);
  }

  /*
  cout << "p2c - move_up_p2c() - " << l1_start << " - " << l1_end
       << ", to L" << layer << " " << l2_start << " - " << l2_end << "\n";
  */


  int lparent = layer_def_[layer].get_parent_layer();
  if (lparent != -1) {
    try_move_up_c2c(l2_start, l2_end, lparent);
  }

  int size = l1_end - l1_start;
  for (int i=0; i<size; i++) {
    Individual* tmp = pop_childs_[l2_start+i];
    pop_childs_[l2_start+i] = pop_parents_[l1_start+i];
    pop_parents_[l1_start+i] = tmp;
  }
}


// Copy individuals from Pop_Parent[l1_start, l1_end)
// to pop_childs_[], starting from the beginning.
// But also try to move up those individs in pop_childs_[].
void AlpsGen :: move_up_p2c(vector<int>& toinsert, int layer)
{
  unsigned int max = toinsert.size();
  if (max == 0) {
    return;
  }

  if (max > layer_def_[layer].get_size()) {
    max = layer_def_[layer].get_size();
  }

  
  unsigned int l2_start = layer_def_[layer].get_start_index();
  unsigned int i=0;
  for (; i<max; i++) {
    if (pop_parents_[toinsert[i]]->compare_fitness(is_maximizing_,
						  pop_childs_[l2_start+i]) == -1) {
      // Doesn't replace.
      break;
    }
  }

  int lparent = layer_def_[layer].get_parent_layer();
  if (lparent != -1) {
    try_move_up_c2c(l2_start, l2_start+i, lparent);
  }

  /*  // !!!!! Debug:
  for (unsigned int j = layer_def_[layer].get_start_index();
       j < layer_def_[layer].get_start_index() + layer_def_[layer].get_size();
       j++) {
    if (pop_childs_[j] == 0) {
      cerr << "alps_gen :: move_up_p2c() - error, child " << j << " is null.\n";
      while (1) ;
    }
  }
  */

  for (unsigned int j=0; j<i; j++) {
    Individual* tmp = pop_childs_[l2_start+j];
    pop_childs_[l2_start+j] = pop_parents_[toinsert[j]];
    pop_parents_[toinsert[j]] = tmp;
  }

  /*  // !!!!! Debug:
  for (unsigned int j=0; j<toinsert.size(); j++) {
    if (pop_parents_[toinsert[j]] == 0) {
      cerr << "alps_gen :: move_up_p2c() - error, parent " << toinsert[j] << " is null.\n";
      while (1) ;
    }
  }
  */

  sort_layer_childs(layer);
}





// Assumes that the parent layer is sorted.
void AlpsGen :: post_process_layer_random(int layer)
{
  int lparent = layer_def_[layer].get_parent_layer();
  if (lparent == -1) {
    // No parent layer, can't do anything.
    return;
  }

  if ((layer == 0) && (generation_ == 0)) {
    // No parents to replace.
    return;
  }

  int l1_start = layer_def_[layer].get_start_index();
  int l1_end = layer_def_[layer].get_start_index() + layer_def_[layer].get_size();
  int l1_index = l1_end - 1;

  int l2_start = layer_def_[lparent].get_start_index();
  int l2_end = l2_start + layer_def_[lparent].get_size();
  int l2_index = l2_start;
    
  do {
    if (pop_parents_[l1_index]->compare_fitness(is_maximizing_,
					       pop_childs_[l2_index]) == -1) {
      // Can't replace, done.
      l1_index++;
      break;
    }
    l1_index--;
    l2_index++;
  } while ((l1_index >= l1_start) && (l2_index < l2_end));

  if (l1_index < l1_end) {
    // Move up individuals from [l1_index, l1_end)
    //  - include l1_index, but stop before l1_end.
    move_up_p2c(l1_index, l1_end, lparent);
  }

  /*
  printf("alps :: post_process_layer_random() - layer %d created randomly.\n",
	 layer);
  */

  return;
}


// Assumes that the parent layer is sorted.
void AlpsGen :: post_process_layer_tourn(int layer)
{
  if (print_debug_info_) {
    cout << "alps_gen :: post_process_layer_tourn(" << layer << "), Gen:"
	 << generation_ << ", MOE:" << moveup_old_elite_
	 << ", MDP:" << moveup_dead_parents_ << "\n";
  }

  if (!moveup_dead_parents_ && !moveup_old_elite_) {
    return;
  }

  int lparent = layer_def_[layer].get_parent_layer();
  if (lparent == -1) {
    // No parent layer, can't do anything.
    //    cout << "  - no parents, can't do anything.\n";
    return;
  }

  // tomove is a vector with elements sorted from best (index 0)
  // to worst (index last)


  vector<int> tomove;
  if (moveup_dead_parents_) {
    tomove.reserve(layer_def_[layer].get_size());
  } else if (moveup_old_elite_) {
    tomove.reserve(layer_def_[layer].get_elitism());
  }


  if (moveup_old_elite_) {
    int max_age = layer_def_[layer].get_max_age();

    int index = layer_def_[layer].get_start_index() + layer_def_[layer].get_size();
    for (unsigned int i=0; i<layer_def_[layer].get_elitism(); i++) {
      index--;
      if (pop_parents_[index]->get_age() > max_age) {
	if (create_type_list_[index] == CREATE_DUPLICATE) {
	  cout << " -- error, max_age:" << max_age << " elite:" << index
	       << ":" << pop_parents_[index]->get_age()
	       << " is too old to be duplicated.\n";
	  while (1) ;
	}
	    
	tomove.push_back(index);
      }
    }

    if (print_debug_info_) {
      cout << "  elite[" << tomove.size() << "]:";
      for (vector<int>::size_type i=0; i<tomove.size(); i++) {
	cout << " " << tomove[i];
      }
      cout << "\n";
    }
  }

  if (moveup_dead_parents_) {
    for (int i=layer_def_[layer].get_start_index() + layer_def_[layer].get_size()-1;
	 i >= (int)layer_def_[layer].get_start_index(); i--) {
      tomove.push_back(i);
    }

    if (print_debug_info_) {
      cout << "  deadP[" << tomove.size() << "]: "
	   << layer_def_[layer].get_start_index() << " - ";
      if (tomove.size() > 0) {
	cout << tomove[tomove.size()-1];
      }
      cout << "\n";
    }
  }

  if (tomove.size() > 0) {
    move_up_p2c(tomove, lparent);
  }

  /* // !!!!! Debug:
  for (unsigned int i = layer_def_[layer].get_start_index();
       i < layer_def_[layer].get_start_index() + layer_def_[layer].get_size();
       i++) {
    if (pop_parents_[i] == 0) {
      cerr << "alps_gen :: post_process_tourn() - error, parent " << i << " is null.\n";
      while (1) ;
    }
  }  
  */

  return;



  if (moveup_old_elite_) {
    int max_age = layer_def_[layer].get_max_age();

    int index = layer_def_[layer].get_start_index() + layer_def_[layer].get_size();
    for (unsigned int i=0; i<layer_def_[layer].get_elitism(); i++) {
      index--;
      if (pop_parents_[index]->get_age() > max_age) {
	tomove.push_back(index);
      }
    }

    if (print_debug_info_) {
      cout << "  elite[" << tomove.size() << "]:";
      for (vector<int>::size_type i=0; i<tomove.size(); i++) {
	cout << " " << tomove[i];
      }
      cout << "\n";
    }
  }



  // !!!!! BUG: Not trying to move up elite individuals if they're
  //   too old.

  int l1_start = layer_def_[layer].get_start_index();
  int l1_end = layer_def_[layer].get_start_index() + layer_def_[layer].get_size()
    - layer_def_[layer].get_elitism();
  int l1_index = l1_end - 1;

  int l2_start = layer_def_[lparent].get_start_index();
  int l2_end = l2_start + layer_def_[lparent].get_size();
  int l2_index = l2_start;
    
  do {
    if (pop_parents_[l1_index]->compare_fitness(is_maximizing_,
					       pop_childs_[l2_index]) == -1) {
      // Can't replace, done.
      l1_index++;
      break;
    }
    l1_index--;
    l2_index++;
  } while ((l1_index >= l1_start) && (l2_index < l2_end));

  /*
  cout << "  - move up L" << layer << ": "
       << l1_index << "(" << pop_parents_[l1_index]->get_fitness() << ") - "
       << l1_end-1 << "(" << pop_parents_[l1_end-1]->get_fitness() << "), L"
       << lparent << ": "
       << l2_start << "(" << pop_childs_[l2_start]->get_fitness() << ") - "
       << l2_index-1 << "(" << pop_childs_[l2_index-1]->get_fitness() << ")\n";
  */

  if (l1_index < l1_end) {
    // Move up individuals from [l1_index, l1_end)
    //  - include l1_index, but stop before l1_end.
    move_up_p2c(l1_index, l1_end, lparent);
  }

  return;
}


//void AlpsGen :: post_process_layer_roulette(int layer)
void AlpsGen :: post_process_layer_roulette(int)
{

}


void AlpsGen :: post_process_layer_dc(int layer)
{
  int start_index = layer_def_[layer].get_start_index();
  int stop_index = layer_def_[layer].get_start_index() + layer_def_[layer].get_size();
  int max_age = layer_def_[layer].get_max_age();
  
  for (int i=start_index; i<stop_index; i++) {
    if ((pop_parents_[i]->get_age() >= max_age) &&
	(max_age >= 0)) {
      move_up(i, layer+1);

    } else if (pop_childs_[i]->compare_fitness(is_maximizing_,
					      pop_parents_[i]) == -1) {
      // Parent is better than child, keep parent.
      /*
      printf(" - restoring parent at index %d (%3.0f:%d > %3.0f:%d).\n",
	     i, pop_parents_[i]->get_fitness(), pop_parents_[i]->get_age(),
	     pop_childs_[i]->get_fitness(), pop_childs_[i]->get_age());
      */
      pop_childs_[i]->duplicate(pop_parents_[i]);

    } else {
      // Child is better, move parent up.
      move_up(i, layer+1);
    }
  }
}


void AlpsGen :: post_process_layers()
{
  // !!!!! Moving up:
  //  - how did the original one-at-a-time version do it?
  //  - how about a flag to turn on/off moving up "dead" parents?
  //  - maybe only move up elite individuals when too old for their layer?


  // The parents are in: pop_parents_.
  // The new individuals are in pop_childs_

  /*
  cout << "alps_gen :: post_process_layers() - Gen: " << generation_
       << ",  MaxLayerEvolving: " << max_layer_evolving_ << "\n";
  */

  //  int max_layer = calc_max_layer(generation_);

  //  for (int i=max_layer_evolving_; i>=0; i--) {
  for (int i=0; i<=max_layer_evolving_; i++) {
    sort_layer_childs(i);
  }


  if (generation_ > 0) {
    for (int i=max_layer_evolving_; i>=0; i--) {
      if (layer_def_[i].get_select_type() == ALPS_SELECT_TOURN) {
	post_process_layer_tourn(i);

      } else if (layer_def_[i].get_select_type() == ALPS_SELECT_ROULETTE) {
	post_process_layer_roulette(i);
	
      } else if (layer_def_[i].get_select_type() == ALPS_SELECT_DC) {
	post_process_layer_dc(i);

      } else {
	cout << "AlpsGen :: post_process_layers() - error, layer " << i
	     << ", invalid selection type: " << layer_def_[i].get_select_type() << endl;
	while (1) ;
      }
    }
  }

  //  cout << "alps_gen :: doing post process" << endl;

  Individual *tmp;
  for (int i=0; i<num_individs_; i++) {
    tmp = pop_childs_[i];
    pop_childs_[i] = pop_parents_[i];
    pop_parents_[i] = tmp;
  }
  // Now the new pop is in pop_parents_ again.

  /*
  cout << "alps :: post_process() - final population: (MaxLayEv:"
       << max_layer_evolving_ << ")\n";
  print_pop();
  */

  return;
}


int AlpsGen :: setup_next_generation()
{
  if (print_debug_info_) {
    cout << "alps_gen :: setup_next_generation() - before post process.\n";
    print_pop();
  }

  post_process_layers();

  if (print_debug_info_) {
    cout << "alps_gen :: setup_next_generation() - after post process.\n";
    print_parent_pop();
  }


  /*
  if (!(generation_%1000) && (generation_ > 1)) {
    write_backup();
  }
  */

  generation_++;

  if (num_runs_ <= 1) {
    if (is_finished_ || ((generation_ > max_generations_) && (!run_sample_))) {
      // write final stats or something.
      return false; // Done!
    }

  } else {
    // May need to reset.
    if (is_finished_ ||
	((generation_ > max_generations_) &&
	 //	((num_evaluations_evolve_ > num_individs_*max_generations_) &&
	 (!run_sample_))) {
      // Write stats (if any).
      run_number_++;
      if (run_number_ >= num_runs_) {
	// Finished!
	return false;
      }
      initialize();
    }
  }

  if ((generation_ > 0) && (save_best_)) {
    //    order_population();
    if ((num_layers_ > 1) || (layer_def_[0].get_size() > 1)) {
      write_best_individ_gen();
    } else if ((num_evaluations_evolve_ % print_results_rate_) == 0) {
      write_best_individ_gen();
    }      
  }
  return setup_next_generation_layers();
}




istream& AlpsGen :: read(istream& istr)
{
  int version = 0;
  string s;
  istr >> s >> version;
  if ((s != "IndV") || (version != ALPS_GEN_VERSION)) {
    throw AlpsError("alps :: read() - error, bad version.");
  }  


  return istr;
}


bool AlpsGen :: read(const char *fname)
{
  ifstream file(fname);
  try {
    read(file);
    //  } catch (AlpsError e) {
    //    return false;
  } catch (...) {
    return false;
  }

  file.close();

  return true;
}



void AlpsGen :: write_log_layers()
{
  if (print_debug_info_) {
    cout << "alps_gen :: write_log_layers()\n";
  }

  int best_index = 0;
  int best_layer = 0;
  double best_overall = get_fitness(0, 0);
  vector<double> layer_avg(num_layers_);
  vector<double> layer_best(num_layers_);
  vector<int> layer_best_index(num_layers_);

  for (int i=0; i<num_layers_; i++) {
    layer_avg[i] = 0.0;
    layer_best[i] = 0.0;
  }

  for (int i=0; i<=max_layer_evolving_; i++) {
    int start_index = layer_def_[i].get_start_index();
    int stop_index = layer_def_[i].get_start_index() + layer_def_[i].get_size();

    double best_fitness = get_fitness(start_index, i);
    double avg_fitness = best_fitness;

    if ((is_maximizing_ && (best_fitness > best_overall)) ||
	(!is_maximizing_ && (best_fitness < best_overall))) {
      best_overall = best_fitness;
      best_index = start_index;
      best_layer = i;
    }
    layer_best_index[i] = start_index;

    for (int j=start_index+1; j<stop_index; j++) {
      double fitness = get_fitness(j, i);
      avg_fitness += fitness;

      if ((is_maximizing_ && (fitness > best_fitness)) ||
	  (!is_maximizing_ && (fitness < best_fitness))) {
	best_fitness = fitness;
	layer_best_index[i] = j;

	if ((is_maximizing_ && (best_fitness > best_overall)) ||
	    (!is_maximizing_ && (best_fitness < best_overall))) {
	  best_overall = best_fitness;
	  best_index = j;
	  best_layer = i;
	}
      }
    }

    layer_avg[i] = avg_fitness;
    layer_best[i] = best_fitness;
  }

  double avg = layer_avg[0];
  int num_individs_ = layer_def_[0].get_size();
  for (int i=1; i<=max_layer_evolving_; i++) {
    avg += layer_avg[i];
    num_individs_ += layer_def_[i].get_size();
  }


  int best_age = pop_childs_[best_index]->get_age();

  if (layer_def_[best_layer].get_select_type() == ALPS_SELECT_DC) {
    if (((pop_parents_[best_index]->get_age() < layer_def_[best_layer].get_max_age()) ||
	 (layer_def_[best_layer].get_max_age() >= 0)) &&
	(pop_childs_[best_index]->compare_fitness(is_maximizing_,
						 pop_parents_[best_index]) == -1)) {
      best_age = pop_parents_[best_index]->get_age();
    }
  }


  log_file_ << num_evaluations_evolve_ << " "
	   << scientific << best_overall << " " << best_age << " " 
	   << avg / double(num_individs_) << " " << num_layers_;

  for (int i=0; i<=max_layer_evolving_; i++) {
    int index = layer_best_index[i];
    int layer_best_age = pop_childs_[index]->get_age();
    if (layer_def_[i].get_select_type() == ALPS_SELECT_DC) {
      if (((pop_parents_[index]->get_age() < layer_def_[i].get_max_age()) ||
	   (layer_def_[i].get_max_age() >= 0)) &&
	  (pop_childs_[index]->compare_fitness(is_maximizing_,
					      pop_parents_[index]) == -1)) {
	layer_best_age = pop_parents_[index]->get_age();
      }
    }

    log_file_ << " " << layer_best[i] << " " << layer_best_age
	     << " " << layer_avg[i]/double(layer_def_[i].get_size());
  }

  for (int i=max_layer_evolving_+1; i<num_layers_; i++) {
    log_file_ << " 0.0 0.0 0";
  }

  //  log_file_ << endl;
  log_file_ << "\n";
}

ostream& AlpsGen :: write_header(ostream& ostr)
{
  ostr << "# AlpsGenV " << ALPS_GEN_VERSION << endl;
  ostr << "# " << num_individs_ << " "
       << max_generations_ << " "
       << num_runs_ << " "
       << prob_recomb_ << " "
       << num_layers_ << " "
       << get_random_seed()
       << endl;


  for (vector<int>::size_type i=0; i<layer_def_.size(); i++) {
    ostr << "# L" << i << ": ";
    layer_def_[i].write_header(ostr);
  }
  

  Alps::write_header(ostr);

  return ostr;
}

ostream& AlpsGen :: write(ostream& ostr)
{
  ostr << "AlpsGenV " << ALPS_GEN_VERSION << endl;

  Alps::write(ostr);

  ostr << num_individs_ << " "
       << max_generations_ << " "
       << prob_recomb_ << " "
       << num_layers_ << " "
       << max_layer_evolving_ << endl;



  for (vector<AlpsLayer>::iterator iter = layer_def_.begin();
       iter != layer_def_.end(); ++iter) {
    ostr << *iter;
  }

  ostr << current_individual_ << " "
       << tries_ << " "
       << evaluated_ << " "
       << create_type_ << " "
       << generation_ << endl;

  ostr << next_to_evaluate_ << " "
       << num_to_evaluate_ << " "
       << num_being_evaluated_ << endl;

  for (vector<int>::iterator iter = evaluation_list_.begin();
       iter != evaluation_list_.end(); ++iter) {
    ostr << *iter << " ";
  }
  ostr << endl;

  for (vector<int>::iterator iter = tries_list_.begin();
       iter != tries_list_.end(); ++iter) {
    ostr << *iter << " ";
  }
  ostr << endl;

  for (vector<int>::iterator iter = create_type_list_.begin();
       iter != create_type_list_.end(); ++iter) {
    ostr << *iter << " ";
  }
  ostr << endl;

  for (vector<int>::iterator iter = parent1_list_.begin();
       iter != parent1_list_.end(); ++iter) {
    ostr << *iter << " ";
  }
  ostr << endl;

  for (vector<int>::iterator iter = parent2_list_.begin();
       iter != parent2_list_.end(); ++iter) {
    ostr << *iter << " ";
  }
  ostr << endl;

  for (vector<int>::iterator iter = evaluation_state_.begin();
       iter != evaluation_state_.end(); ++iter) {
    ostr << *iter << " ";
  }
  ostr << endl;

  for (vector<int>::size_type i=0; i<pop_parents_.size(); i++) {
    pop_parents_[i]->write(ostr);
  }
  ostr << endl;

  return ostr;
}

bool AlpsGen :: write(const char *filename)
{
  ofstream ofile(filename);
  write(ofile);
  ofile.close();
  return true;
}


bool AlpsGen :: write_backup()
{
  char fname[60];

  make_filename("_bk.pop", fname);
  return write(fname);
}


} // namespace alps

/********************************************************/
