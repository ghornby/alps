/*******************************************************
  Copyright Gregory S. Hornby.
  All Rights Reserved.
********************************************************/

#ifndef ALPS_SSTATE_HEADER_FILE
#define ALPS_SSTATE_HEADER_FILE

#include <iostream>
#include <vector>


#include "alps.h"


namespace alps {

class Individual;


class AlpsSState : public Alps
{
public:
  AlpsSState(const char *n, Individual* ind_config);
  AlpsSState(Individual* ind_config);
  ~AlpsSState();

  void init_variables(void);
  void initialize(void);

  int get_size();
  void set_max_evals(int n);
  void set_max_gen(int n);
  void set_recomb_prob(double prob);
  void set_rec_rand2_prob(double prob);


  void print_generation_stats();
  void print_state(void);


  int compare_parents(int index1, int index2);
  int compare_elite(int layer, int index1, int index2);
  bool is_finished();

  void insert_evaluated(std::vector<double>& fitness, int index,
			 Individual* individ, int print_individ_history);
  void insert_evaluated(int index, Individual* individ,
			int print_individ_history);

  void evaluate_error(int index, Individual* individ);



  Individual *get_individual(unsigned int index);
  int get_next_individ(int& index, Individual*& individ);


  bool read_sample_individ(const char *fname);
  bool read_seed_individ(const char *fname);

  std::istream& read(std::istream& istr);
  bool read(const char *filename);
  std::ostream& write_header(std::ostream& ostr);
  std::ostream& write(std::ostream& ostr);
  bool write(const char *filename);
  bool write_backup(void);


private:
  /* Variables */
  bool do_init_all_;
  bool do_init_lyr_;
  bool select_sequential_;

  int init_next_;
  int next_index_;

  int max_evals_;
  int num_individs_;

  double prob_recomb_;
  double prob_rec_rand2_;


  std::vector<bool> recalc_elitism_;

  std::vector<std::vector<int> > elite_individs_;
 
  Individual *sample_individ_; // Sample for configuring things, like number of genes/inputs.
  Individual *seed_individ_; // Seed for an evolutionary run (super-seeds Sample_Individ).

  std::vector<Individual*> pop_parents_;
  std::vector<Individual*> pop_childs_;

  
  /* Functions */
  void set_pop_size(unsigned int size);
  void sort_layer_parents(int layer);
  void write_log_layers(void);

  int get_age_move(int index);
  int current_age_move();
  void set_age_move(int index);
  void tryk_move_up(int index);

  int stochastic_parent1(int index);
  int stochastic_parent2(int l, unsigned int p1);
  int random_parent2(int l, unsigned int p1);
  void start_init_layer0();
  Individual* make_mutate(int index);
  Individual* make_recombine(int index);
  Individual* make_random_individ();

  void cout_parent(int index);
  void sort_elite(int lyr, int index);
  void sort_elite(int lyr);
  void print_elite(int lyr);
  void update_elite(int index);
  bool is_elite(int index);
  int get_ind_age(int index);
  Individual* get_new_individ();

};

}

#endif
