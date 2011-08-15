/*******************************************************

  File: alps.cpp
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

#include "alps/alps_utils.h"
using namespace alps_utils;

#include "alps/alps.h"
#include "alps/alps_layer.h"
#include "alps/alps_individ.h"

#include <stdexcept>
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>

using namespace std;


namespace alps {

ostream& operator<<(ostream& os, Alps& a)
{
  return a.write(os);
}

istream& operator>>(istream& is, Alps& a)
{
  return a.read(is);
}




const int MAX_TRIES = 1;
const int ALPS_VERSION = 1;


// Evaluation States
const int ES_NOT_USED = -1;
const int ES_NOT_CREATED = 0;
const int ES_EVALUATING = 1;
const int ES_EVALUATED = 2;



/***** Initialization functions *****/

Alps :: Alps(const char *n, Individual* ind_config)
  : alps_name_(n), is_maximizing_(true), is_finished_(false),
    print_gen_stats_(true), max_individ_(0), sample_individ_(0),
    seed_individ_(0), layer_def_(0)
{
  set_name(n);
  init_variables();

  sample_individ_ = ind_config->new_instance();
  sample_individ_->duplicate_settings(ind_config);
}

Alps :: Alps(Individual* ind_config)
  : is_maximizing_(true), is_finished_(false), print_gen_stats_(true),
    max_individ_(0), sample_individ_(0), seed_individ_(0)
{
  Alps("evolve", ind_config);
}


Alps :: ~Alps()
{
  if (log_file_.is_open()) {
    log_file_.close();
  }

  if (max_individ_) {
    delete max_individ_;
  }

  if (sample_individ_) {
    delete sample_individ_;
  }

  if (seed_individ_) {
    delete seed_individ_;
  }
}


// Probably want to log:
// a. Layer stats
// b. structural organization stats.


void Alps :: init_variables(void)
{
  save_log_ = false;

  run_sample_ = false;
  max_individ_ = 0;
  num_runs_ = 1;
  run_number_ = 0;
  print_results_rate_ = 100;

  save_log_ = false;

  seed_individ_ = NULL;
  sample_individ_ = NULL;
  //  sample_individ_ = new Individual();

  print_debug_info_ = false;
}


void Alps :: initialize(void)
{
  num_evaluations_ = 0;
  num_evaluations_init_ = 0;
  num_evaluations_evolve_ = 0;
}


void Alps :: set_name(const char *n)
{
  alps_name_.assign(n);
}


void Alps :: make_filename(const char *ending, char *fname)
{
  int i;

  for (i=0; i<40; i++) {
    fname[i] = alps_name_[i];
    if (fname[i] == 0)
      break;
  }

  if (num_runs_ > 1) {
    sprintf(&fname[i], ".%d%s", run_number_, ending);
  } else {
    strcpy(&fname[i], ending);
  }
}


int Alps :: set_logfile(const char *n)
{
  if (n == 0) {
    return false;
  }

  if (n[0] == 0) {
    return false;
  }

  log_file_name_.assign(n);
  log_file_name_ += ".log";

  save_log_ = true;

  // Open log file.
  log_file_.open(log_file_name_.c_str());
  cout << "Opening log: " << log_file_name_ << "\n";

  return true;
}


  //int Alps :: configure(const char *fname, int verbose)
int Alps :: configure(const char *, int)
{
  printf("pop :: configure() - no longer supported..\n");
  while (1) ;

  return true;
}




// ************************************************************ //
// *************** Layer Configuration Functions ************** //

void Alps :: print_layers()
{
  if (layer_def_.size() == 0) {
    cout << "Alps :: print_layers() - no layers allocated.\n";
    return;
  }

  cout << "Alps :: print_layers() - " << num_layers_ << " layers." << endl;

  int ctr = 0;
  for (vector<AlpsLayer>::iterator iter = layer_def_.begin();
       iter != layer_def_.end(); ++iter) {
    cout << "L" << ctr++ << ": ";
    iter->print_info();
  }
}


void Alps :: set_pop_size(unsigned int size)
{
  cout << "alps :: set_pop_size(" << size << ")\n";
}


int Alps :: get_layer(int index)
{
  int count = 0;
  for (unsigned int i=0; i<layer_def_.size(); i++) {
    count += layer_def_[i].get_size();
    if (index < count) {
      return i;
    }
  }

  cout << "alps :: get_layer(1) - error, index exceeds maximum layer: "
       << index << ".\n";

  while (1) ;
  return num_layers_;
}



int Alps :: calc_max_age_in_layer(int aging_scheme, int age_gap, int layer)
{
  if (age_gap <= 0) {
    age_gap = 1;
  }

  if (aging_scheme == ALPS_AGING_LINEAR) {
    return age_gap * (layer + 1);

  } else if (aging_scheme == ALPS_AGING_LINEAR1) {
    return age_gap * (layer + 1) - layer;

  } else if (aging_scheme == ALPS_AGING_POLY1) {
    //   0:1, 1:2, 2:4, 3:9, 4:16, 5:25 ....
    if (layer == 0) {
      return age_gap;
    } else if (layer == 1) {
      return age_gap + age_gap;
    }

    return layer * layer * age_gap;

  } else if (aging_scheme == ALPS_AGING_POLY2) {
    if (layer == 0) {
      if (age_gap > 1) {
	return age_gap + 2;
      }
      return 1;

    } else if (layer == 1) {
      return age_gap + age_gap + 1;

    } else if (layer == 2) {
      return 4 * age_gap;

    } else if (layer == 3) {
      return 8 * age_gap - 1;

    } else if (layer == 4) {
      return 16 * age_gap - 3;
    }

    return layer * layer * age_gap;

  } else if (aging_scheme == ALPS_AGING_EXP) {
    //   0:1, 1:2, 2:4, 3:8, 4:16, 5:32 ....
    if (layer == 0) {
      return age_gap;
    } else if (layer == 1) {
      return age_gap + age_gap;
    }

    int k = 4;
    for (int i=2; i<layer; i++) {
      k *= 2;
    }

    k *= age_gap;
    if (k < 0) {
      printf("alps :: calc_max_layer_age(3) - error doing exponential on layer %d: %d %d.\n",
	     k, age_gap, layer);
      while (1) ;
    }
    return k;

  } else if (aging_scheme == ALPS_AGING_EXP2) {
    //   0:1, 1:2, 2:4, 3:8, 4:16, 5:32 ....
    if (layer == 0) {
      return age_gap;
    } else if (layer == 1) {
      return age_gap + age_gap - 1;
    }

    int k = 1 << layer;
    k *= age_gap;
    int age_subt = 1 << layer;
    age_subt -= 1;
    k -= age_subt;
    if (k < 0) {
      printf("alps :: calc_max_layer_age(3) - error doing exponential on layer %d: %d %d.\n",
	     k, age_gap, layer);
      while (1) ;
    }
    return k;

  } else if (aging_scheme == ALPS_AGING_EXP3) {
    if (layer == 0) {
      if (age_gap > 1) {
	return age_gap;
      }
      return 1;

    } else if (layer == 1) {
      if (age_gap > 1) {
	return age_gap + age_gap - 1;
      }
    }

    int max = age_gap + age_gap - 1;
    for (int i=1; i<layer; i++) {
      max = max + max -1;
    }
    return max;

  } else if (aging_scheme == ALPS_AGING_FIBONACCI) {
    if (layer == 0) {
      return age_gap;
    } else if (layer == 1) {
      return age_gap + age_gap;
    }

    int num1 = age_gap;
    int num2 = age_gap;
    for (int i=1; i<=layer; i++) {
      int num3 = num1 + num2;
      num1 = num2;
      num2 = num3;
    }
    return num2;

  } else if (aging_scheme == ALPS_AGING_FIBONACCI1) {
    if (layer == 0) {
      return age_gap;
    }


    int num1 = age_gap;
    int num2 = age_gap;
    while (num2 <= 2) {
      int num3 = num2;
      num2 += num1;
      num1 = num3;
    }
    if (layer == 1) {
      return num1 + num2 - 1;
    }

    for (int i=1; i<=layer; i++) {
      int num3 = num2;
      num2 += num1 - 1;
      num1 = num3;
    }
    return num2;
  }

  printf("alps :: calc_max_layer_age(3) - error, invalid aging scheme: %d.\n",
	 aging_scheme);
  while (1) ;
}


void Alps :: set_layer_elitism(unsigned int l, unsigned int num)
{
  if (l >= layer_def_.size()) {
    throw invalid_argument("alps_gen :: set_layer_elitism() - error, invalid layer.");
  }

  layer_def_[l].set_elitism(num);
}


void Alps :: config_layer(AlpsLayer& ldef)
{
  layer_def_.resize(1);

  layer_def_[0] = ldef;
  layer_def_[0].set_max_age(-1);
  layer_def_[0].set_start_index(0);
  layer_def_[0].set_parent_layer(-1);

  vector<int> v;
  layer_def_[0].set_prev_layers(0, v);

  num_layers_ = 1;
  int num = layer_def_[0].get_size();

  set_pop_size(num);
  initialize();
}



void Alps :: config_layer(int size, int select_type, int elitism)
{
  layer_def_.resize(1);

  num_layers_ = 1;

  layer_def_[0].set_elitism(elitism);
  layer_def_[0].set_select_type(select_type);
  layer_def_[0].set_max_age(-1);
  layer_def_[0].set_start_index(0);
  layer_def_[0].set_size(size);
  layer_def_[0].set_parent_layer(-1);

  vector<int> v;
  layer_def_[0].set_prev_layers(0, v);
  layer_def_[0].set_prob_select_prev(0.0);

  int num = layer_def_[0].get_size();

  set_pop_size(num);
  initialize();
}


//void Alps :: config_layers_unique(int age_scheme, int age_gap,
//				  int num_layers_, vector<AlpsLayer>& def)
void Alps :: config_layers_unique(int, int, int nlayers,
				  vector<AlpsLayer>& def)
{
  if (num_layers_ < 1) {
    throw AlpsError("not enough layers.");
    return;
  }

  num_layers_ = nlayers;
  layer_def_.resize(num_layers_);

  for (int i=0; i<num_layers_; i++) {
    layer_def_[i] = def[i];
  }

  int num = 0;
  for (int i=0; i<num_layers_; i++) {
    num += layer_def_[i].get_size();
  }
  
  set_pop_size(num);
  initialize();
}


void Alps :: config_layers_same(int aging_scheme, int age_gap,
				      int nlayers, AlpsLayer& def)
{
  num_layers_ = nlayers;
  if (num_layers_ < 1) {
    cout << "alps_gen :: config_layers_same() - error, num_layers_ < 1 : " <<
      num_layers_ << endl;
    throw AlpsError("not enough layers.");
    return;
  }

  /*
  cout << "alps :: config_layers_same(" << aging_scheme << ", "
       << age_gap << ", " << nlayers << ", layer_def_)\n";
  */
  
  layer_def_.resize(num_layers_);

  for (int i=0; i<num_layers_; i++) {
    layer_def_[i] = def;

    layer_def_[i].set_max_age(calc_max_age_in_layer(aging_scheme, age_gap, i));
    layer_def_[i].set_parent_layer(i+1);
    if (layer_def_[i].get_parent_layer() == num_layers_) {
      layer_def_[i].set_parent_layer(-1);
    }

    vector<int> prev;

    if (i > 0) {
      prev.resize(1);
      prev[0] = i-1;      
      layer_def_[i].set_prev_layers(1, prev);
    } else {
      prev.resize(0);
      layer_def_[0].set_prev_layers(0, prev);
    }
  }

  layer_def_[0].set_start_index(0);
  layer_def_[0].set_prob_select_prev(0.0);
  for (int i=1; i<num_layers_; i++) {
    layer_def_[i].set_start_index(layer_def_[i-1].get_start_index()
			      + layer_def_[i-1].get_size());
  }

  int num = 0;
  for (int i=0; i<num_layers_; i++) {
    num += layer_def_[i].get_size();
  }

  for (int i=num_layers_-1; i>=0; i--) {
    if (layer_def_[i].get_parent_layer() == -1) {
      layer_def_[i].set_max_age(-1);
    }
  }

  set_pop_size(num);
  initialize();
}


/*
void Alps :: config_layers_tree(int aging_scheme, int age_gap,
				int num_single_layers, int num_tree_layers,
				int branch_factor, AlpsLayer& def)
*/
void Alps :: config_layers_tree(int, int,
				int num_single_layers, int num_tree_layers,
				int branch_factor, AlpsLayer&)
{
  num_layers_ = 0;

  vector<int> width(num_tree_layers);
  if (num_tree_layers > 0) {
    width[num_tree_layers-1] = branch_factor;
    num_layers_ += branch_factor;
  }

  for (int i=num_tree_layers-2; i>=0; i--) {
    width[i] = width[i+1]*branch_factor;
    num_layers_ += width[i];
  }

  num_layers_ += num_single_layers;

  layer_def_.resize(num_layers_);

  vector<int> prev_layers;
  for (int i=0; i<num_layers_; i++) {
    layer_def_[i].set_prev_layers(0, prev_layers);
  }


  layer_def_[0].set_start_index(0);
  for (int i=1; i<num_layers_; i++) {
    layer_def_[i].set_start_index(layer_def_[i-1].get_start_index() + layer_def_[i-1].get_size());

    if (layer_def_[i].get_num_prev_layers() <= 0) {
      layer_def_[i].set_prob_select_prev(0.0);
    }
  }

  int num = 0;
  for (int i=0; i<num_layers_; i++) {
    num += layer_def_[i].get_size();
  }

  for (int i=num_layers_-1; i>=0; i--) {
    if (layer_def_[i].get_parent_layer() == -1) {
      layer_def_[i].set_max_age(-1);
    }
  }

  set_pop_size(num);
  initialize();
}



// ************************************************************** //
// ************************************************************** //


void Alps :: set_num_runs(int n)
{
  num_runs_ = n;
}

void Alps :: set_run_number(int r)
{
  run_number_ = r;
}




void Alps :: set_print_debug(bool f)
{
  print_debug_info_ = f;
}

void Alps :: set_print_gen_stats(bool f)
{
  print_gen_stats_ = f;
}


void Alps :: set_print_results_rate(int r)
{
  print_results_rate_ = r;
}

void Alps :: set_run_sample(bool b)
{
  run_sample_ = b;
}

void Alps :: set_save_best(bool b)
{
  save_best_ = b;
}


int Alps :: get_num_evals()
{
  return num_evaluations_;
}

int Alps :: get_num_evals_evolve()
{
  return num_evaluations_evolve_;
}

int Alps :: get_num_evals_init()
{
  return num_evaluations_init_;
}

int Alps :: get_size()
{
  return 0;
}


int Alps :: get_num_layers()
{
  return layer_def_.size();
}


void Alps :: set_finished(bool b)
{
  is_finished_ = b;
}

bool Alps :: is_finished()
{
  return is_finished_;
}


// Returns the next individual to evaluate and its index.
// Returns -1 if next individual is not ready yet.
// Returns -2 if all evaluations are done.
//int Alps :: get_next_individ(int& index, Individual*& individ)
int Alps :: get_next_individ(int&, Individual*&)
{
  return -3;
}

/*
void Alps :: insert_evaluated(std::vector<double>& fitness, int index,
			      Individual* individ, int print_individ_history)
*/
void Alps :: insert_evaluated(std::vector<double>&, int,
			      Individual*, int)
{
}

/*
void Alps :: insert_evaluated(int index, Individual* individ,
			      int print_individ_history)
*/
void Alps :: insert_evaluated(int, Individual*, int)
{
}

//void Alps :: evaluate_error(int index, Individual* individ)
void Alps :: evaluate_error(int, Individual*)
{
}



//Individual *Alps :: get_individual(unsigned int index)
Individual *Alps :: get_individual(unsigned int)
{
  return 0;
}


Individual *Alps :: get_current_individual()
{
  return 0;
}



extern void print_time(void);

void Alps :: print_generation_stats(void)
{
}


void Alps :: print_generation_stats2(void)
{
}


void Alps :: print_state()
{
}



//void Alps :: add_to_stats(int index)
void Alps :: add_to_stats(int)
{
}


void Alps :: update_stats()
{
}


void Alps :: print_stats()
{
  printf("Evaluations %d :: #init: %d  #evolve: %d\n",
	 num_evaluations_, num_evaluations_init_,
	 num_evaluations_evolve_);
}


void Alps :: write_stats(void)
{
}


void Alps :: print_final_stats(void)
{

}



// Note: probably need to give this function an existing
//  individual to read into.
// For more ellaborate individual types, this individual
//  would already be properly configured with a translator function.
bool Alps :: read_sample_individ(const char *fname)
{
  if (run_sample_) {
    printf("pop :: read_sample_individ() - already loaded individual.\n");
    return false;
  }

  if (sample_individ_ == NULL) {
    sample_individ_ = new Individual();
  }

  return sample_individ_->read(fname);
}


bool Alps :: read_seed_individ(const char *fname)
{
  if (sample_individ_ == NULL) {
    seed_individ_ = new Individual();
  }

  return seed_individ_->read(fname);
}


//void Alps :: write_best_individ(const char* fname)
void Alps :: write_best_individ(const char*)
{
}


void Alps :: write_best_individ()
{
}


// Writes best individual in the population after it
// has been ordered (thus best is in index Order[0]).
void Alps :: write_best_individ_gen()
{
}



// Parses a line from the configuration file,
// returns false if a problem, otherwise true.
//bool Alps :: parse_config(const char *config_line, int verbose)
bool Alps :: parse_config(const char *, int)
{
  return false; // Not a recognized line type.
}



istream& Alps :: read(istream& istr)
{
  int version = 0;
  string s;
  istr >> s >> version;
  if ((s != "IndV") || (version != ALPS_VERSION)) {
    throw AlpsError("alps :: read() - error, bad version.");
  }  

  istr >> alps_name_;
  if (alps_name_ == "-") {
    alps_name_.clear();
  }

  istr >> log_file_name_;
  if (log_file_name_ == "-") {
    log_file_name_.clear();
  }

  istr >> num_runs_
       >> run_number_
       >> is_maximizing_
       >> run_sample_
       >> print_results_rate_
       >> save_best_
       >> save_log_;

  istr >> num_evaluations_
       >> num_evaluations_init_
       >> num_evaluations_evolve_;


  return istr;
}


bool Alps :: read(const char *fname)
{
  ifstream file(fname);
  try {
    read(file);
    //  } catch (AlpsError e) {
  } catch (...) {
    return false;
  }

  file.close();

  return true;
}


ostream& Alps :: write_header(ostream& ostr)
{
  ostr << "# AlpsV " << ALPS_VERSION << endl;

  ostr << "# ";
  if (alps_name_[0] != 0) {
    ostr << alps_name_ << endl;
  } else {
    ostr << "-" << endl;
  }

  ostr << "# ";
  if (log_file_name_[0] != 0) {
    ostr << log_file_name_ << endl;
  } else {
    ostr << "-" << endl;
  }

  ostr << "# ";
  ostr << num_runs_ << " "
       << run_number_ << " "
       << is_maximizing_ << " "    
       << run_sample_ << " "
       << print_results_rate_ << " "
       << save_best_ << " "
       << save_log_ << endl;

  ostr << "# ";
  ostr << num_evaluations_ << " "
       << num_evaluations_init_ << " "
       << num_evaluations_evolve_ << endl;

  return ostr;
}


ostream& Alps :: write(ostream& ostr)
{
  ostr << "AlpsV " << ALPS_VERSION << endl;

  if (alps_name_[0] != 0) {
    ostr << alps_name_ << endl;
  } else {
    ostr << "-" << endl;
  }

  if (log_file_name_[0] != 0) {
    ostr << log_file_name_ << endl;
  } else {
    ostr << "-" << endl;
  }

  ostr << num_runs_ << " "
       << run_number_ << " "
       << is_maximizing_ << " "    
       << run_sample_ << " "
       << print_results_rate_ << " "
       << save_best_ << " "
       << save_log_ << endl;

  ostr << num_evaluations_ << " "
       << num_evaluations_init_ << " "
       << num_evaluations_evolve_ << endl;

  return ostr;
}


bool Alps :: write(const char *fname)
{
  ofstream file(fname);
  write(fname);
  file.close();

  return true;
}



bool Alps :: write_backup(void)
{
  char fname[60];

  make_filename("_bk.pop", fname);
  return write(fname);
}


}
