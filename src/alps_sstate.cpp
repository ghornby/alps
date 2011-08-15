/*******************************************************

  File: alps_sstate.cpp
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


#include "alps/alps_random_mt.h"
using namespace alps_random;

#include "alps/alps_utils.h"
using namespace alps_utils;

#include "alps/alps.h"
#include "alps/alps_sstate.h"
#include "alps/alps_layer.h"
#include "alps/alps_individ.h"
#include "alps/alps_history.h"

using namespace std;


namespace alps {

const int MAX_TRIES = 1;
const int ALPS_SSTATE_VERSION = 1;


// Evaluation States
const int ES_NOT_USED = -5;
const int ES_NOT_CREATED = -4;
const int ES_EVALUATING = 1;
const int ES_EVALUATED = 2;



// *************************************************** //
// ************ Initialization functions **************//

AlpsSState :: AlpsSState(const char *n, Individual* ind_config) :
  Alps(n, ind_config)
{
  set_name(n);
  
  seed_individ_ = 0;
  sample_individ_ = ind_config->new_instance();
  sample_individ_->duplicate_settings(ind_config);

  init_variables();
}


AlpsSState :: AlpsSState(Individual* ind_config) :
  Alps("evolve", ind_config)
{
  set_name("evolve");

  seed_individ_ = 0;
  sample_individ_ = ind_config->new_instance();
  sample_individ_->duplicate_settings(ind_config);

  init_variables();
}


AlpsSState :: ~AlpsSState()
{
}


void AlpsSState :: set_pop_size(unsigned int size)
{
  cout << "alps_sstate :: set_pop_size(" << size << ")\n";

  vector<Individual*>::size_type old_size = pop_parents_.size();

  for (vector<Individual*>::size_type i=size; i<pop_parents_.size(); i++) {
    delete pop_parents_[i];
  }

  pop_parents_.resize(size);

  for (vector<Individual*>::size_type i=old_size; i<size; i++) {
    pop_parents_[i] = sample_individ_->new_instance();
  }
  num_individs_ = size;
  next_index_ = num_individs_-1;
}



void AlpsSState :: init_variables()
{
  do_init_all_ = true;
  do_init_lyr_ = false;
  init_next_ = 0;

  select_sequential_ = false;
  select_sequential_ = true;
  next_index_ = num_individs_-1;

  max_evals_ = 100000;
  num_individs_ = 0;
  prob_recomb_ = 0.5;
  prob_rec_rand2_ = 0.5;
  num_layers_ = 0;
}


void AlpsSState :: initialize()
{
  num_evaluations_ = 0;
  num_evaluations_init_ = 0;
  num_evaluations_evolve_ = 0;
  
  for (int i=0; i<num_individs_; i++) {
    pop_parents_[i]->initialize();
  }

  elite_individs_.resize(layer_def_.size());
  for (unsigned int i=0; i<elite_individs_.size(); i++) {
    //    elite_individs_[i].resize(layer_def_[i].get_elitism());
    elite_individs_[i].resize(0);
  }
}



int AlpsSState :: get_size()
{
  return num_individs_;
}

void AlpsSState :: set_max_evals(int n)
{
  max_evals_ = n;
}

void AlpsSState :: set_max_gen(int n)
{
  max_evals_ = n * num_individs_;
}

void AlpsSState :: set_recomb_prob(double prob)
{
  prob_recomb_ = prob;
}

void AlpsSState :: set_rec_rand2_prob(double prob)
{
  prob_rec_rand2_ = prob;
}



// ********************************************** //
// *************** Print Functions ************** //

void AlpsSState :: print_generation_stats()
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

  for (int i=0; i<num_layers_; i++) {
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
  for (int i=1; i<num_layers_; i++) {
    avg += layer_avg[i];
    num_individs_ += layer_def_[i].get_size();
  }

  int best_age = pop_parents_[best_index]->get_age(num_evaluations_,
						  num_individs_);


  if (print_debug_info_) {
    cout << "\n";
  }

  int prec = cout.precision();
  cout.precision(2);
  //  cout << num_evaluations_ << " "
  cout << fixed << num_evaluations_ << "  "
       << scientific << setw(4) << best_overall << " "
       << setw(3) << best_age << "  "
       << fixed << num_layers_;
  cout.precision(1);

  if (save_log_) {
    log_file_.precision(9);
    log_file_ << num_evaluations_ << " "
	     << scientific << best_overall << fixed << " " << best_age << " " 
	     << setprecision(3) << avg / double(num_individs_) << " " << num_layers_;
    log_file_.precision(5);
  }

  for (int i=0; i<num_layers_; i++) {
    int index = layer_best_index[i];
    int layer_best_age = pop_parents_[index]->get_age(num_evaluations_, num_individs_);

    cout << ",  " << layer_best[i] << " " << setw(2) << layer_best_age;
    if (save_log_) {
      log_file_ << " " << layer_best[i] << " " << layer_best_age
	       << " " << layer_avg[i]/double(layer_def_[i].get_size());
    }
  }

  if (save_log_) {
    log_file_ << endl;
  }

  cout << " - ";
  //  individual_print_op_history();

  /*
  printf("  : <%2d:%2d:%2d> ", Generation,
	 Current_Individual, Max_Index);
  */
  cout << " ";
  print_time();
  cout << "\n";
  cout.precision(prec);

  if (print_debug_info_) {
    cout << "\n";
  }
}


void AlpsSState :: print_state(void)
{
  printf("<%2d>  evals:%d:%d:%d",
	 num_evaluations_/num_individs_,
	 num_evaluations_, num_evaluations_init_,
	 num_evaluations_evolve_);

  printf("  -  ");
  individual_print_op_history();
  printf("   ");
  print_time();
  printf("\n");  
}



// ************************************************************** //
// *************** Selection/Replacement Functions ************** //

// Returns:
//   1 if index1 better fitness than index2
//   0 if equal
//   -1 if index1 strictly worse than index2
int AlpsSState :: compare_parents(int index1, int index2)
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


int AlpsSState :: compare_elite(int layer, int index1, int index2)
{
  return compare_parents(elite_individs_[layer][index1],
			 elite_individs_[layer][index2]);
}


// start_index will have worst fitness.
// stop_index-1 will have best fitness.
void AlpsSState :: sort_layer_parents(int layer)
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



// Returns true if no more evaluations to be done.
bool AlpsSState :: is_finished()
{
  if (Alps::is_finished()) {
    return true;
  }

  if (num_evaluations_ > max_evals_) {
    set_finished(true);
    return true;
  }
  return false;
}



int AlpsSState :: get_age_move(int index)
{
  return pop_parents_[index]->get_age_move();
}


int AlpsSState :: current_age_move()
{
  return num_evaluations_/num_individs_;
}


void AlpsSState :: set_age_move(int index)
{
  return pop_parents_[index]->set_age_move(current_age_move());
}




// Similar to try_move_up() - but checks k(=?) individs and picks worst
void AlpsSState :: tryk_move_up(int index)
{
  if (print_debug_info_) {
    cout << "  tryk_move_up(" << index << ")\n";
  }

  int l = get_layer(index);  
  int lparent = layer_def_[l].get_parent_layer();
  if (lparent < 0) {
    return;
  }

  set_age_move(index);
  int max_age = layer_def_[lparent].get_max_age();

  int lstart = layer_def_[lparent].get_start_index();
  int lstop = layer_def_[lparent].get_size() + lstart;
  
  int istart = random_int(layer_def_[lparent].get_size()) + lstart;
  int i = istart;

  vector<int> repl_age;
  vector<int> repl_fitness;

  int worst_fit = -1;
  int worst_fit2 = -1;
  int worst_age = -1;
  
  int ctr = 20; // This is k.
  if ((int)layer_def_[lparent].get_size() < ctr) {
    ctr = layer_def_[lparent].get_size();
  }
  //  ctr = layer_def_[lparent].get_size(); // Uncomment to test entire layer.
  while (ctr-- > 0) {
    i++;
    if (i == lstop) {
      i = lstart;
    }

    if ((pop_parents_[i]->get_age(num_evaluations_, num_individs_) > max_age)
	&& (max_age > 0)) {
      // Replace first individ that is too old for this layer.
      worst_age = i;
      repl_age.push_back(i);
      break;

      repl_age.push_back(i);
      if (worst_age == -1) {
	worst_age = i;
      } else {
	if (compare_parents(i, worst_age) == -1) {
	  worst_age = i;
	}
      }

      continue;
//      break;
    }

    if (compare_parents(i, index) != 1) {
      // Replaceable:
      if (get_age_move(i)+1 < current_age_move()) {
	// Try to replace this first:
	if (worst_fit == -1) {
	  worst_fit = i;
	} else {
	  if (compare_parents(i, worst_fit) == -1) {
	    worst_fit = i;
	  }
	}

      } else {
	if (worst_fit2 == -1) {
	  worst_fit2 = i;
	} else {
	  if (compare_parents(i, worst_fit2) == -1) {
	    worst_fit2 = i;
	  }
	}
	repl_fitness.push_back(i);
      }

      continue;
    }
  }

  int r = -1;
  if (repl_age.size() > 0) {
    if (repl_age.size() == 1) {
      r = repl_age[0];
    } else {
      r = repl_age[random_int(repl_age.size())];
    }
    r = worst_age;

  } else if (worst_fit != -1) {
    r = worst_fit;

  } else if (repl_fitness.size() > 0) {
    if (repl_fitness.size() == 1) {
      r = repl_fitness[0];
    } else {
      r = repl_fitness[random_int(repl_fitness.size())];
    }
    r = worst_fit2;
  }

  if (r != -1) {
    tryk_move_up(r);
    Individual* tmp = pop_parents_[r];
    pop_parents_[r] = pop_parents_[index];
    pop_parents_[index] = tmp;
    update_elite(r); // i is where individ at index is moved to, update_elite on it.
  }
}


/*
void AlpsSState :: insert_evaluated(vector<double>& fitness,
				    int index, Individual* individ,
				    int print_individ_history)
*/
void AlpsSState :: insert_evaluated(vector<double>& fitness,
				    int index, Individual* individ,
				    int)
{
  if (individ == 0) {
    cout << "alps_ss :: insert_evaluated() - error, individ is null.\n";
    while (1) ;
  }

  if (print_debug_info_) {
    cout << "Insert_Evaluated(" << index << ") - fitness=" << fitness[0] << "\n";
  }

  individ->set_fitness(is_maximizing_, fitness);
  num_evaluations_++;

  if (!do_init_all_) {
    //    try_move_up(index);
    tryk_move_up(index);
  }
  pop_childs_.push_back(pop_parents_[index]);
  pop_parents_[index] = individ;

  update_elite(index);

  if (save_log_ && (num_evaluations_ == 1)) {
    write_header(log_file_);
  }

  if (print_gen_stats_ && (num_evaluations_ > 0)) {
    if ((num_evaluations_ % print_results_rate_) == 0)  {
      /*
      if (save_log_) {
      // Inside of "print_generation_stats()"
	write_log_layers();
      }
      */

      // Print some stats.
      print_generation_stats();
    }
  }
}


void AlpsSState :: insert_evaluated(int index, Individual* individ,
				    int print_individ_history)
{
  insert_evaluated(individ->get_fitness_vec(), index, individ,
		   print_individ_history);
}


//void AlpsSState :: evaluate_error(int index, Individual* individ)
void AlpsSState :: evaluate_error(int, Individual* individ)
{  
  num_evaluations_++;
  pop_childs_.push_back(individ);
}


Individual *AlpsSState :: get_individual(unsigned int index)
{
  if (index >= pop_parents_.size()) {
    printf("pop :: get_individual() - invalid index: %d.\n", index);
    while (1) ;
    return NULL;
  }

  return pop_parents_[index];
}



// Uses tournament selection to pick parent 2.
int AlpsSState :: stochastic_parent2(int l, unsigned int p1)
{
  // !!!!! This could be made to run faster by maintaining a
  //       list of valid parents for each layer.
  //       If a young individual becomes too old, it is removed
  //       If an old individual is replaced with a young child, it is added.

  int max_age = layer_def_[l].get_max_age();

  //  if (layer_def_[l].get_select_type() == ALPS_SELECT_TOURN) {

  vector<int> tourn;

  bool avail_this_layer = true;
  bool avail_prev_layer = true;
  if ((layer_def_[l].get_prob_select_prev() == 0.0) ||
      (l == 0)) {
    avail_prev_layer = false;
  }

  while ((tourn.size() < layer_def_[l].get_tourn_size()) &&
	 (avail_this_layer || avail_prev_layer)) {
    if (random_double() < layer_def_[l].get_prob_select_prev()) {
      // Select from previous layer:
      int l_prev = -1;
      if (layer_def_[l].get_num_prev_layers() == 1) {
	l_prev = layer_def_[l].get_prev_layer(0);
      } else {
	l_prev = layer_def_[l].get_prev_layer((unsigned int)random_int(layer_def_[l].get_num_prev_layers()));
      }

      if (get_layer(p1) == l_prev) {
	unsigned int i = random_int(layer_def_[l_prev].get_size()-1
				    + layer_def_[l_prev].get_start_index());
	if (i==p1) {
	  i = layer_def_[l_prev].get_size()-1
	    + layer_def_[l_prev].get_start_index();
	}
	tourn.push_back(i);

      } else {
	int i = random_int(layer_def_[l_prev].get_size()
			   + layer_def_[l_prev].get_start_index());
	tourn.push_back(i);
      }

    } else {
      // Select from this layer:
      unsigned int i = random_int(layer_def_[l].get_size())
	+ layer_def_[l].get_start_index();
      if (i == p1) {
	i++;
	if (i >= layer_def_[l].get_stop_index()) {
	  i = layer_def_[l].get_start_index();
	}
      }

      unsigned int istart = i;
      while ((pop_parents_[i]->get_age(num_evaluations_, num_individs_)
	      > max_age) && (max_age > 0)) {
	i++;
	if (i == p1) {
	  i++;
	}
	if (i >= layer_def_[l].get_stop_index()) {
	  i = layer_def_[l].get_start_index();
	}

	if (i == istart) {
	  // wrapped around, no more valid individs to select!
	  avail_this_layer = false;
	  break;
	}
      }
      
      if ((pop_parents_[i]->get_age(num_evaluations_, num_individs_)
	   <= max_age) || (max_age < 0)) {
	tourn.push_back(i);
      } else if (layer_def_[l].get_num_prev_layers() == 0) {
	break;
      }
    }
  }

  if (tourn.size() == 0) {
    return -1;
  }

  int p2 = tourn[0];
  for (unsigned int i=1; i<tourn.size(); i++) {
    if (compare_parents(p2, tourn[i]) != 1) {
      p2 = tourn[i];
    }
  }
  return p2;
}


// Completely random parent2.
int AlpsSState :: random_parent2(int l, unsigned int p1)
{
  //  if (layer_def_[l].get_select_type() == ALPS_SELECT_TOURN) {

  if (random_double() >= layer_def_[l].get_prob_select_prev()) {
    // Select from this layer:
    int max_age = layer_def_[l].get_max_age();
    unsigned int i = random_int(layer_def_[l].get_size()-1)
      + layer_def_[l].get_start_index();
    if (i == p1) {
      i = layer_def_[l].get_last_index();
    }

    unsigned int istart = i;
    while ((pop_parents_[i]->get_age(num_evaluations_, num_individs_)
	    > max_age) && (max_age > 0)) {
      i++;
      if (i == p1) {
	i++;
      }

      if (i >= layer_def_[l].get_stop_index()) {
	i = layer_def_[l].get_start_index();
      }

      if (i == istart) {
	// wrapped around, no more valid individs to select!
	break;
      }
    }
      
    if ((pop_parents_[i]->get_age(num_evaluations_, num_individs_)
	 <= max_age) || (max_age < 0)) {
      return i;
    }

    if (layer_def_[l].get_num_prev_layers() == 0) {
      // No previous layers, so select an individual at random.
      unsigned int i = random_int(layer_def_[l].get_size()-1)
	+ layer_def_[l].get_start_index();
      if (i == p1) {
	i = layer_def_[l].get_last_index();
      }
      return i;
    }

    // There is a previous layer, so randomly select from that:
  }


  // Select from previous layer:
  int l_prev = -1;
  if (layer_def_[l].get_num_prev_layers() == 1) {
    l_prev = layer_def_[l].get_prev_layer(0);
  } else {
    l_prev = layer_def_[l].get_prev_layer((unsigned int)random_int(layer_def_[l].get_num_prev_layers()));
  }

  if (get_layer(p1) == l_prev) {
    unsigned int i = random_int(layer_def_[l_prev].get_size()-1
				+ layer_def_[l_prev].get_start_index());
    if (i==p1) {
      i = layer_def_[l_prev].get_size()-1
	+ layer_def_[l_prev].get_start_index();
    }
    return i;
    
  } else {
    int i = random_int(layer_def_[l_prev].get_size()
		       + layer_def_[l_prev].get_start_index());
    return i;
  }
}



int AlpsSState :: stochastic_parent1(int index)
{
  // !!!!! This could be made to run faster by maintaining a
  //       list of valid parents for each layer.
  //       If a young individual becomes too old, it is removed
  //       If an old individual is replaced with a young child, it is added.

  int l = get_layer(index);
  int max_age = layer_def_[l].get_max_age();


  if (layer_def_[l].get_select_type() == ALPS_SELECT_DC) {
    if ((pop_parents_[index]->get_age(num_evaluations_, num_individs_)
	 <= max_age) || (max_age < 0)) {
      return index;
    }
    // Individ at this index is too old, so do tournament selection...
  }      


  vector<int> tourn;
  if ((max_age < 0) || (get_ind_age(index) <= max_age)) {
    if (print_debug_info_) {
      cout << "     stochastic_parent1() - push index:" << index << ", max_age:"
	   << max_age << ", age: "
	   << pop_parents_[index]->get_age(num_evaluations_, num_individs_) << "\n";
    }
    tourn.push_back(index);
  }

  bool avail_this_layer = true;
  bool avail_prev_layer = true;
  double prob_select_prev = layer_def_[l].get_prob_select_prev();
  if ((prob_select_prev == 0.0) || (l == 0)) {
    avail_prev_layer = false;
  }
  
  while ((tourn.size() < layer_def_[l].get_tourn_size()) &&
	 (avail_this_layer || avail_prev_layer)) {
    if (random_double() < prob_select_prev) {
      if (l == 0) {
	cout << " alps_sstate :: stochastic_parent1() -- error, no previous to layer 0.\n";
	while (1) ;
      }

      // Select from previous layer:
      int l_prev = -1;
      if (layer_def_[l].get_num_prev_layers() == 1) {
	l_prev = layer_def_[l].get_prev_layer(0);
      } else {
	l_prev = layer_def_[l].get_prev_layer((unsigned int)random_int(layer_def_[l].get_num_prev_layers()));
      }

      int i = random_int(layer_def_[l_prev].get_size()
			 + layer_def_[l_prev].get_start_index());
      tourn.push_back(i);

    } else {
      // Select from this layer:
      unsigned int i = random_int(layer_def_[l].get_size())
	+ layer_def_[l].get_start_index();

      unsigned int istart = i;
      while ((pop_parents_[i]->get_age(num_evaluations_, num_individs_)
	      > max_age) && (max_age > 0)) {
	i++;
	if (i == layer_def_[l].get_start_index() + layer_def_[l].get_size()) {
	  i = layer_def_[l].get_start_index();
	}

	if (i == istart) {
	  // wrapped around, no more valid individs to select!
	  avail_this_layer = false;
	  break;
	}
      }

      if ((pop_parents_[i]->get_age(num_evaluations_, num_individs_)
	   <= max_age) || (max_age < 0)) {
	tourn.push_back(i);
      } else if (layer_def_[l].get_num_prev_layers() > 0) {
	prob_select_prev = 1.0;
      } else {
	break;
      }
    }	  
  }

  if (tourn.size() == 0) {
    return -1;
  }
  int p1 = tourn[0];
  for (unsigned int i=1; i<tourn.size(); i++) {
    if (compare_parents(p1, tourn[i]) != 1) {
      p1 = tourn[i];
    }
  }
  return p1;
}


void AlpsSState :: start_init_layer0()
{
  do_init_lyr_ = true;
  init_next_ = 0;
  elite_individs_[0].resize(0);
}


Individual* AlpsSState :: make_mutate(int index)
{
  int p1 = stochastic_parent1(index);
  if (p1 == -1) {
    // No valid parents in layer 0.
    if (print_debug_info_) {
      cout << "alps_sstate :: make_mutate(" << index << ") - switch to init.\n";
    }
    return 0;
  }

  if (print_debug_info_) {
    cout << "alps_sstate :: make_mutate(" << index << ") - p1:" << p1
	 << "(" << pop_parents_[p1]->get_age(num_evaluations_, num_individs_) << ")\n";
  }

  Individual* individ = get_new_individ();
  individ->duplicate(pop_parents_[p1]);
  individ->mutate();

  return individ;
}


Individual* AlpsSState :: make_recombine(int index)
{
  int p1 = stochastic_parent1(index);
  if (p1 == -1) {
    if (print_debug_info_) {
      // No valid parents in layer 0.
      cout << "alps_sstate :: make_recombine(" << index << ") - switch to init.\n";
    }
    return 0;
  }

  int p2 = -1;
  double choice = random_double();
  if (choice < prob_rec_rand2_) {
    p2 = random_parent2(get_layer(index), p1);
  } else {
    p2 = stochastic_parent2(get_layer(index), p1);
  }
  if (p2 == -1) {
    // Only 1 valid parent in layer 0 so mutate it.
    Individual* individ = get_new_individ();
    individ->duplicate(pop_parents_[p1]);
    individ->mutate();
    return individ;
  }

  if (print_debug_info_) {
    cout << "alps_sstate :: make_recomb(" << index << ") - p1:"
	 << p1 << "(" << pop_parents_[p1]->get_age(num_evaluations_, num_individs_) << "), p2:"
	 << p2 << "(" << pop_parents_[p2]->get_age(num_evaluations_, num_individs_) << ")\n";
  }

  Individual* individ = get_new_individ();
  individ->duplicate(pop_parents_[p1]);

  if (choice < 1.5) {
    individ->recombine(pop_parents_[p2]);
  } else {
    individ->recombine_rand2(pop_parents_[p2]);
  }

  return individ;
}


Individual* AlpsSState :: get_new_individ()
{
  // Rather free'ing and then new'ing childs all the time,
  // the discarded individs are stored in pop_childs_.
  // It should never grow larger than the number of threads
  // using this object.

  if (pop_childs_.size() > 0) {
    Individual* ind = pop_childs_.back();
    pop_childs_.pop_back();
    if (ind == 0) {
      cout << "alps_ss :: get_new_individ() - error, individ is null.\n";
      while (1) ;
    }
    return ind;
  }

  return sample_individ_->new_instance();
}



void AlpsSState :: cout_parent(int index)
{
  cout << index << "::"
       << pop_parents_[index]->get_fitness()
       << ":" << get_ind_age(index);
}



// Correctly places index in lyr's elite.
// Assumes the vector is already filled.
void AlpsSState :: sort_elite(int lyr, int index)
{
  bool sort_up = true; // If true, check index+1, index+2, ...

  if (index == 0) {
    //    sort_up = true;
  } else if (index+1 == (int)elite_individs_[lyr].size()) {
    sort_up = false;
  } else {
    if (compare_elite(lyr, index-1, index) != 1) {
      //      sort_up = true;
    } else {
      sort_up = false;
    }
  }

  if (sort_up) {
    for (vector<int>::size_type i = index+1; 
	 i < elite_individs_[lyr].size(); i++) {
      if (compare_elite(lyr, i-1, i) != -1) {
	// Done sorting.
	break;
      }
      int tmp = elite_individs_[lyr][i-1];
      elite_individs_[lyr][i-1] = elite_individs_[lyr][i];
      elite_individs_[lyr][i] = tmp;
    }

  } else {
    for (vector<int>::size_type i = index; 
	 i > 0; i--) {
      if (compare_elite(lyr, i-1, i) != -1) {
	// Done sorting.
	break;
      }
      int tmp = elite_individs_[lyr][i-1];
      elite_individs_[lyr][i-1] = elite_individs_[lyr][i];
      elite_individs_[lyr][i] = tmp;
    }
  }
}


// Performs a single pass of trying to move the bottom individ up.
void AlpsSState :: sort_elite(int lyr)
{
  int max_age = layer_def_[lyr].get_max_age();

  if (max_age == -1) {
    // Most fit is first index, least fit is last index.
    for (vector<int>::size_type i = elite_individs_[lyr].size()-1;
	 i > 0; i--) {
      if (compare_elite(lyr, i-1, i) != -1) {
	// Done sorting.
	break;
      }
      int tmp = elite_individs_[lyr][i-1];
      elite_individs_[lyr][i-1] = elite_individs_[lyr][i];
      elite_individs_[lyr][i] = tmp;
    }

  } else {
    // Most fit is first index, least fit is last index.
    for (vector<int>::size_type i = elite_individs_[lyr].size()-1;
	 i > 0; i--) {
      if ((pop_parents_[i-1]->get_age(num_evaluations_, num_individs_) <= max_age) &&
	  (compare_elite(lyr, i-1, i) != -1)) {
	// Done sorting.
	if (print_debug_info_) {
	  cout << "  - done elite_insert:   ";
	  cout_parent(elite_individs_[lyr][i-1]);
	  cout << " is better then ";
	  cout_parent(elite_individs_[lyr][i]);
	  cout << "\n";

	    /*
	       << i-1 << ":" << pop_parents_[elite_individs_[lyr][i-1]]->get_fitness()
	       << " is better than: "
	       << i << ":" << pop_parents_[elite_individs_[lyr][i]]->get_fitness() << "\n";
	    */
	}
	break;
      }
      if (print_debug_info_) {
	cout << "  - swapping:            ";
	cout_parent(elite_individs_[lyr][i-1]);
	cout << " with ";
	cout_parent(elite_individs_[lyr][i]);
	cout << "\n";
	/*	  
	     << i-1 << ":" << pop_parents_[elite_individs_[lyr][i-1]]->get_fitness()
	     << " with: "
	     << i << ":" << pop_parents_[elite_individs_[lyr][i]]->get_fitness() << "\n";
	*/
      }
      int tmp = elite_individs_[lyr][i-1];
      elite_individs_[lyr][i-1] = elite_individs_[lyr][i];
      elite_individs_[lyr][i] = tmp;
    }
  }

  if (max_age > 0) {
    while (elite_individs_[lyr].size() > 0) {
      if (pop_parents_[elite_individs_[lyr].back()]->get_age(num_evaluations_, num_individs_) <= max_age) {
	break;
      }
      elite_individs_[lyr].pop_back(); // !!!! Had some errors with pop_back()
    }
  }
}


void AlpsSState :: print_elite(int lyr)
{
  if (!print_debug_info_) {
    return;
  }

  if (elite_individs_[lyr].size() > layer_def_[lyr].get_size()) {
    cout << "ERROR - invalid elitism size for layer: " << lyr << ".\n";
    while (1) ;
  }

  cout << "L" << lyr << ": elite(" << elite_individs_[lyr].size() << "):";
  for (vector<int>::size_type i=0; i<elite_individs_[lyr].size(); i++) {
    cout << " ";
    cout_parent(elite_individs_[lyr][i]);
  }
  cout << "\n";
}



void AlpsSState :: update_elite(int index)
{
  int l = get_layer(index);

  if (layer_def_[l].get_elitism() == 0) {
    return;
  }
  
  if (print_debug_info_) {
    cout << "Update_Elite() - (";
    cout_parent(index);
    cout << ")\n";
    /*
    << index << "("
	 << pop_parents_[index]->get_fitness() << ":"
	 << get_ind_age(index) << ")\n";
    */

    cout << "  Start:  ";
    print_elite(l);
  }

  // First, check to see if already on the list.
  // If already on list, then resort list and we're done.
  int max_age = layer_def_[l].get_max_age();
  for (vector<int>::size_type i=0; i<elite_individs_[l].size(); i++) {
    if ((get_ind_age(elite_individs_[l][i]) > max_age) &&
	(max_age > 0)) {
      // Remove individuals that are too old.
      if (i+1 < elite_individs_[l].size()) {
	elite_individs_[l][i] = elite_individs_[l][elite_individs_[l].size()-1];
      }
      elite_individs_[l].pop_back();
      continue;
    }	

    if (elite_individs_[l][i] == index) {
      sort_elite(l, i);
      if (print_debug_info_) {
	cout << "  Update: ";
	print_elite(l);
      }
      return;
    }
  }

  // Not already on list.
  if ((get_ind_age(index) > max_age) && (max_age > 0)) {
    return;
  }

  // If space on the list, add and sort.
  if (elite_individs_[l].size() < layer_def_[l].get_elitism()) {
    elite_individs_[l].push_back(index);
    sort_elite(l);
    if (print_debug_info_) {
      cout << "  - Add:  ";
      print_elite(l);
    }
    return;
  }

  elite_individs_[l].push_back(index);
  sort_elite(l);
  elite_individs_[l].resize(layer_def_[l].get_elitism());

  if (print_debug_info_) {
    cout << "  Insert: ";
    print_elite(l);
  }
}


bool AlpsSState :: is_elite(int index)
{
  int l = get_layer(index);

  if (print_debug_info_) {
    cout << "Is_Elite() - (";
    cout_parent(index);
    cout << ")\n";
  }

  update_elite(index);

  for (vector<int>::size_type i=0; i<elite_individs_[l].size(); i++) {
    if (elite_individs_[l][i] == index) {
      if ((get_ind_age(index) > layer_def_[l].get_max_age())
	  && (layer_def_[l].get_max_age() > 0)) {
	// Too old: remove from elite.
	if (i+1 != elite_individs_[l].size()) {
	  elite_individs_[l][i] = elite_individs_[l][elite_individs_[l].size()-1];
	}
	elite_individs_[l].pop_back();

	if (print_debug_info_) {
	  cout << "IsElite(): index: " << index << "::"
	       << pop_parents_[index]->get_fitness() << ":" << get_ind_age(index)
	       << " is too old (max:"
	       << layer_def_[l].get_max_age() << "), removed from elite: ";
	  print_elite(l);
	}

	return false;
      }

      if (print_debug_info_) {
	cout << "IsElite(): index: " << index << "::"
	     << pop_parents_[index]->get_fitness() << ":" << get_ind_age(index)
	     << " is elite: ";
	print_elite(l);
      }
      return true;
    }
  }

  if (print_debug_info_) {
    cout << "IsElite(): index: " << index << "::"
	 << pop_parents_[index]->get_fitness() << ":" << get_ind_age(index)
	 << " not elite: ";
    print_elite(l);
  }
  return false;
}


int AlpsSState :: get_ind_age(int index)
{
  return pop_parents_[index]->get_age(num_evaluations_, num_individs_);
}


Individual* AlpsSState :: make_random_individ()
{
  Individual* individ = get_new_individ();
  individ->make_random();
  //    individ->set_age(num_evaluations_);
  individ->set_age(-num_evaluations_);
  return individ;
}



// Returns the index of the next individual to evaluate.
// Returns -1 if next individual is not ready yet.
// Returns -2 if all evaluations are done.
int AlpsSState :: get_next_individ(int& index, Individual*& individ)
{
  if (print_debug_info_) {
    cout << "ALPS_SS :: get_next()\n";
  }


  if (do_init_all_) {
    // Making the original pop all random individs.
    index = init_next_;
    individ = make_random_individ();

    if (print_debug_info_) {
      cout << "  - do_init on:" << index
	   << "(" << individ->get_age(num_evaluations_, num_individs_) << ")\n";
    }

    init_next_++;
    if (init_next_ >= num_individs_) {
      init_next_ = 0;
      do_init_all_ = false;
      //	cout << "alps_sstate :: get_next() - done init_all\n";
    }

    return index;
  }

  if (select_sequential_) {
    // Go through the indices from num_individs_-1 down to 0.
    if (do_init_lyr_) {
      index = 0;

    } else {
      do {
	index = next_index_;
	next_index_--;
	if (next_index_ < 0) {
	  next_index_ = num_individs_-1;
	}
      } while (is_elite(index));
    }

  } else {
    do {
      index = random_int(num_individs_);
      if ((get_layer(index) == 0) && (do_init_lyr_)) {
	// If time to re-initialize layer 0, no elitism on it.
	break;
      }
    } while (is_elite(index));
  }

  if (print_debug_info_) {
    cout << "  - picked index:" << index
	 << "(" << individ->get_age(num_evaluations_, num_individs_) << ")\n";
  }

  if ((get_layer(index) == 0) && (do_init_lyr_)) {
    // Putting new random individs in bottom layer (sequentially).
    index = init_next_;

    individ = get_new_individ();
    individ->make_random();
    //    individ->set_age(num_evaluations_);
    individ->set_age(-num_evaluations_);

    if (print_debug_info_) {
      cout << "  - init layer on:" << index
	   << "(" << individ->get_age(num_evaluations_, num_individs_) << ")\n";
    }

    init_next_++;
    if (init_next_ >= (int)layer_def_[0].get_size()) {
      init_next_ = 0;
      do_init_lyr_ = false;
      next_index_ = pop_parents_.size() - 1;
      //	cout << "alps_sstate :: get_next() - done init_layer\n";
    }

    return index;
  }
  
  if (random_double() < prob_recomb_) {
    if (print_debug_info_) {
      cout << "alps_ss :: get_next() - mutate:" << index << "\n";
    }
    individ = make_mutate(index);
  } else {
    if (print_debug_info_) {
      cout << "alps_ss :: get_next() - recomb:" << index << "\n";
    }
    individ = make_recombine(index);
  }


  if (individ == 0) {
    start_init_layer0();
    index = init_next_++;
    if (init_next_ >= (int)layer_def_[0].get_size()) {
      init_next_ = 0;
      do_init_lyr_ = false;
    }

    individ = make_random_individ();

    if (print_debug_info_) {
      cout << "  - FORCED INIT:" << index
	   << "(" << individ->get_age(num_evaluations_, num_individs_) << ")\n";
    }
  }

  return index;
}










// ******************************************** //
// *************** I/O Functions ************** //


bool AlpsSState :: read_sample_individ(const char *fname)
{
  if (sample_individ_ == 0) {
    sample_individ_ = new Individual();
  }

  return sample_individ_->read(fname);
}


bool AlpsSState :: read_seed_individ(const char *fname)
{
  if (sample_individ_ == 0) {
    seed_individ_ = new Individual();
  }
  
  return seed_individ_->read(fname);
}


istream& AlpsSState :: read(istream& istr)
{
  int version = 0;
  string s;
  istr >> s >> version;
  if ((s != "IndV") || (version != ALPS_SSTATE_VERSION)) {
    throw AlpsError("alps_sstate :: read() - error, bad version.");
  }  


  return istr;
}


bool AlpsSState :: read(const char *fname)
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



void AlpsSState :: write_log_layers()
{
  if (print_debug_info_) {
    cout << "alps_sstate :: write_log_layers()\n";
  }

  int best_index = 0;
  int best_layer = 0;
  double best_overall = pop_parents_[0]->get_fitness();
  vector<double> layer_avg(num_layers_);
  vector<double> layer_best(num_layers_);
  vector<int> layer_best_index(num_layers_);

  for (int i=0; i<num_layers_; i++) {
    layer_avg[i] = 0.0;
    layer_best[i] = 0.0;
  }

  for (int i=0; i<(int)layer_def_.size(); i++) {
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
  for (unsigned int i=1; i<layer_def_.size(); i++) {
    avg += layer_avg[i];
    num_individs_ += layer_def_[i].get_size();
  }


  int best_age = pop_parents_[best_index]->get_age(num_evaluations_, num_individs_);

  log_file_ << num_evaluations_ << " "
	   << scientific << best_overall << " " << best_age << " " 
	   << avg / double(num_individs_) << " " << num_layers_;

  for (unsigned int i=0; i<layer_def_.size(); i++) {
    int index = layer_best_index[i];
    int layer_best_age = pop_parents_[index]->get_age(num_evaluations_, num_individs_);
    log_file_ << " " << fixed << layer_best[i] << " " << layer_best_age
	     << " " << layer_avg[i]/double(layer_def_[i].get_size());
  }

  //  log_file_ << endl;
  log_file_ << "\n";
}


ostream& AlpsSState :: write_header(ostream& ostr)
{
  ostr << "# AlpsSStateV " << ALPS_SSTATE_VERSION << endl;
  ostr << "# " << num_individs_ << " "
       << max_evals_ << " "
       << prob_recomb_ << " "
       << num_layers_ << " "
       << get_random_seed()
       << endl;


  /* 
  for (vector<AlpsLayer>::iterator iter = layer_def_.begin();
       iter != layer_def_.end(); ++iter) {
    ostr << "# " << *iter;
  }
  */
  for (vector<int>::size_type i=0; i<layer_def_.size(); i++) {
    ostr << "# L" << i << ": ";
    layer_def_[i].write_header(ostr);
  }
  

  Alps::write_header(ostr);

  return ostr;
}

ostream& AlpsSState :: write(ostream& ostr)
{
  ostr << "AlpsSStateV " << ALPS_SSTATE_VERSION << endl;

  Alps::write(ostr);

  ostr << num_individs_ << " "
       << max_evals_ << " "
       << prob_recomb_ << " "
       << num_layers_ << " "
       << get_random_seed() << endl;



  for (vector<AlpsLayer>::iterator iter = layer_def_.begin();
       iter != layer_def_.end(); ++iter) {
    ostr << *iter;
  }

  for (vector<int>::size_type i=0; i<pop_parents_.size(); i++) {
    pop_parents_[i]->write(ostr);
  }
  ostr << endl;

  return ostr;
}

bool AlpsSState :: write(const char *filename)
{
  ofstream ofile(filename);
  write(ofile);
  ofile.close();
  return true;
}


bool AlpsSState :: write_backup()
{
  char fname[60];

  make_filename("_bk.pop", fname);
  return write(fname);
}


} // namespace alps

/********************************************************/
