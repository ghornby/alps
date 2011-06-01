/*******************************************************
  Copyright Gregory S. Hornby.
  All Rights Reserved.
********************************************************/

#ifndef POP_HEADER_FILE
#define POP_HEADER_FILE


#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace alps {

class Individual;
class AlpsLayer;


const int ALPS_PARADIGM_CANONICAL = 0;
const int ALPS_PARADIGM_DETERMINISTIC_CROWDING = 1;

const int ALPS_SELECT_ROULETTE = 1;
const int ALPS_SELECT_TOURN = 2;
const int ALPS_SELECT_DC = 3;

const int ALPS_AGING_LINEAR = 1;
const int ALPS_AGING_LINEAR1 = 2;
const int ALPS_AGING_POLY1 = 3;
const int ALPS_AGING_POLY2 = 4;
const int ALPS_AGING_EXP = 5;
const int ALPS_AGING_EXP2 = 6;
const int ALPS_AGING_EXP3 = 7;
const int ALPS_AGING_FIBONACCI = 8;
const int ALPS_AGING_FIBONACCI1 = 9;

const int ALPS_PARADIGM_ALPS = 2;




/*
Want to be able to select with:
- roulette wheel
- tournament
- deterministic crowding
*/



class Alps
{
  friend std::ostream& operator<<(std::ostream&, Alps&);
  friend std::istream& operator>>(std::istream&, Alps&);

public:
  Alps(void);
  Alps(const char *n, Individual* ind_config);
  Alps(Individual* ind_config);

  virtual ~Alps(void);

  virtual void init_variables(void);
  virtual void initialize(void);
  void set_name(const char *n);
  void make_filename(const char *ending, char *fname);
  int set_logfile(const char *n);

  void set_print_debug(bool f);
  void set_print_gen_stats(bool f);
  void set_print_results_rate(int r);
  void set_run_sample(bool b);
  void set_save_best(bool b);

  virtual void set_num_runs(int n);
  virtual void set_run_number(int n);


  bool get_maximize() const { return is_maximizing_; };
  void set_maximize() { is_maximizing_ = true; };
  void set_minimize() { is_maximizing_ = false; };

  virtual bool parse_config(const char *config_line, int verbose);
  virtual int configure(const char *fname, int verbose);

  int get_num_evals();
  int get_num_evals_evolve();
  int get_num_evals_init();

  virtual int get_size();
  int get_num_layers();

  // *** Layer functions: ***
  void print_layers();
  int get_layer(int index);
  int calc_max_age_in_layer(int paradigm, int age_gap, int layer);
  void set_layer_elitism(unsigned int layer, unsigned int elitism);

  void config_layer(AlpsLayer& def);
  void config_layer(int size, int select_type, int elitism);
  void config_layers_unique(int aging_scheme, int age_gap,
			    int num_layers, std::vector<AlpsLayer>& def);
  void config_layers_same(int aging_scheme, int age_gap, int num_layers,
			  AlpsLayer& def);
  void config_layers_tree(int aging_scheme, int age_gap, int num_single_layers,
			  int num_tree_layers, int branch_factor, AlpsLayer& def);
  // *** * ***




  virtual bool is_finished();
  void set_finished(bool b);

  virtual int get_next_individ(int& index, Individual*& individ);
  virtual void insert_evaluated(std::vector<double>& fitness, int index,
				Individual* individ, int print_individ_history);
  virtual void insert_evaluated(int index, Individual* individ,
				int print_individ_history);
  virtual void evaluate_error(int index, Individual* individ);

  virtual Individual *get_individual(unsigned int index);
  virtual Individual *get_current_individual(void); // Data structure specific.

  virtual void print_generation_stats(void);
  virtual void print_generation_stats2(void);
  virtual void print_state(void);
  virtual void add_to_stats(int index);
  virtual void update_stats(void);
  virtual void print_stats(void);
  virtual void write_stats(void);


  virtual void print_final_stats(void);

  bool read_sample_individ(const char *fname);
  bool read_seed_individ(const char *fname);

  virtual void write_best_individ(const char* fname);
  virtual void write_best_individ(void);
  virtual void write_best_individ_gen(void);

  virtual std::ostream& write_header(std::ostream& ostr);
  virtual std::istream& read(std::istream& file);
  virtual bool read(const char *fname);
  virtual std::ostream& write(std::ostream& file);
  virtual bool write(const char *fname);

  virtual bool write_backup(void);  

protected:
  std::string alps_name_;

  int num_runs_;
  int run_number_;
  bool is_maximizing_;

  bool is_finished_;

  int run_sample_;
  int print_results_rate_;
  bool print_gen_stats_;
  bool save_best_;
  bool save_log_;
  std::string log_file_name_;
  std::ofstream log_file_;
  std::string individ_log_file_name_;

  bool print_debug_info_;

  // Generated Statistics
  int num_evaluations_;
  int num_evaluations_evolve_;
  int num_evaluations_init_;

  Individual *max_individ_;
  Individual *sample_individ_; // Sample for configuring things, like number of genes/inputs.
  Individual *seed_individ_; // Seed for an evolutionary run (super-seeds Sample_Individ).


  // Layer Variables:
  int num_layers_;
  std::vector<AlpsLayer> layer_def_;
  virtual void set_pop_size(unsigned int size);
};

}

#endif
