/*******************************************************
  Copyright Gregory S. Hornby.

********************************************************/

#ifndef ALPS_GEN_HEADER_FILE
#define ALPS_GEN_HEADER_FILE

#include "alps.h"

#include <iostream>
#include <vector>


namespace alps {

class Individual;
class AlpsLayer;


class AlpsGen : public Alps
{
public:
  //  typedef Individual* iterator;
  typedef std::vector<Individual*>::iterator iterator;
  //  typedef const Individual* const_iterator;
  typedef std::vector<Individual*>::const_iterator const_iterator;
  typedef size_t size_type;
  typedef Individual value_type;
  typedef std::ptrdiff_t difference_type;
  typedef Individual& reference;
  typedef const Individual& const_reference;


  AlpsGen(const char *n, Individual* ind_config);
  AlpsGen(Individual* ind_config);
  ~AlpsGen();

  void init_variables(void);
  void initialize(void);

  bool parse_config(const char *config_line, int verbose, int (*translator_func)(char*));
  bool configure(char *fname, int verbose, int (*translator_func)(char*));

  int get_size(void);
  int get_num_evals(void);

  const_iterator begin();
  const_iterator end();


  void set_moveup_dead_parents(bool b);
  void set_moveup_old_elite(bool b);

  void set_tourn_size(unsigned int size);
  void set_max_gen(int n);
  void set_num_runs(int n);
  void set_run_number(int n);
  void set_recomb_prob(double prob);
  void set_rec_rand2_prob(double prob);

  void print_pop();
  void print_parents();
  void print_parent_pop();
  void print_evaluation_state(void);

  void print_write_gen_stats();
  int get_next_individ(int& index, Individual*& individ);
  Individual *get_individual(unsigned int index);
  Individual *get_current_individual(void); // Data structure specific.
  Individual *get_parent(int index);

  bool is_finished();

  void print_generation_stats(void);
  void print_state(void);

  Individual* get_current(int index);
  int get_layer(int index);
  void insert_evaluated(std::vector<double>& fitness, int index,
			Individual* individ, int print_individ_history);
  void insert_evaluated(int index, Individual* individ,
			int print_individ_history);
  double get_fitness(int index, int layer);
  void evaluate_error(int index, Individual* individ);

  int compare_childs(int index1, int index2);
  int compare_parents(int index1, int index2);


  void print_final_stats(void);

  bool read_sample_individ(const char *fname);
  bool read_seed_individ(const char *fname);
  void write_best_individ(std::ofstream& outfile);
  void write_best_individ(void);
  void write_best_individ_gen(void);

  std::istream& read(std::istream& istr);
  bool read(const char *filename);
  std::ostream& write_header(std::ostream& ostr);
  std::ostream& write(std::ostream& ostr);
  bool write(const char *filename);
  bool write_backup(void);


private:
  // Static settings.
  int max_generations_;
  int num_individs_;
  double prob_recomb_;
  double prob_rec_rand2_;

  bool reval_duplicates_; // Re-evaluate duplicates.
  bool moveup_dead_parents_;
  bool moveup_old_elite_;
  int max_layer_evolving_;

  int (*update_correlation_func)(Individual*);


  std::vector<int> pop_order_;  // Not really used


  // Parent and Child populations.
  int max_index;
  std::vector<Individual*> pop_parents_;
  std::vector<Individual*> pop_childs_;


  // Current state of evolution.
  int current_individual_;
  int tries_;
  int evaluated_;
  int create_type_;
  int generation_;

  // State for network pop.
  int next_to_evaluate_;
  int num_to_evaluate_;
  int num_being_evaluated_;
  std::vector<int> evaluation_list_;
  std::vector<int> tries_list_;
  std::vector<int> create_type_list_;
  std::vector<int> parent1_list_;
  std::vector<int> parent2_list_;
  std::vector<std::vector<Individual*> > tourn_list_;
  std::vector<int> evaluation_state_;


  void set_pop_size(unsigned int size);

  void sort_layer_childs(int layer);
  void sort_parents(std::vector<Individual*>& inds, bool sort_maximize);
  void sort_layer_parents(int layer);
  void order_population(void);
  
  bool make_individual(int index);
  int random_parent(void);
  int random_parent2(int index1);
  int calc_max_layer(int gen);
  int calc_max_layer(void);

  void setup_layer_tourn(int layer);
  void setup_layer_roulette(int layer);
  void setup_layer_dc(int layer);
  void setup_layer_random(int layer);

  bool move_up(int index, int layer);
  void try_move_up_c2c(int l1_start, int l1_end, int layer);
  void move_up_p2c(int l1_start, int l1_end, int layer);

  void move_up_p2c(std::vector<int>& toinsert, int layer);

  void post_process_layer_random(int layer);
  void post_process_layer_tourn(int layer);
  void post_process_layer_roulette(int layer);
  void post_process_layer_dc(int layer);
  void post_process_layers(void);

  int setup_next_generation_layers(void);
  int setup_next_generation_alps(void);
  int setup_next_generation_helper(void);
  int setup_next_generation(void);

  void write_log_layers(void);
};

}

#endif
