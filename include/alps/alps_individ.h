/************************************************************
  Copyright Gregory S. Hornby.
  All Rights Reserved.
************************************************************/
#ifndef ALPS_INDIVIDUAL_HEADER_FILE
#define ALPS_INDIVIDUAL_HEADER_FILE

#include <iostream>
#include <fstream>

#include "alps_history.h"

const double Worst_Fitness = 1.5E+15;
const int NUM_HISTORY_OPS = 15;
const int NUM_FITNESS = 20;
const int HIST_NUM_COMPLEXITY = 27;


namespace alps {

extern int gNum_Fitness;
extern int gNum_Operators;
extern double gRec_Delta;
extern double gMutate_Size;



class Individual
{
public:
  Individual();
  Individual(Individual *sample);

  virtual Individual* new_instance();
  virtual ~Individual();


  /*** Interface Methods ***/
  virtual void initialize();
  virtual bool valid();


  double fitness_hist_diff(int index);
  double fitness_hist_diff();
  virtual int better(bool maximize);
  virtual bool keep(bool maximize);
  virtual int compare_fitness(bool maximizing, Individual *individ2);
  virtual bool same(Individual *individ2);
  bool created_by_variation();

  double get_fitness();
  void get_fitness(std::vector<double>& fitness_vec);
  std::vector<double>& get_fitness_vec();
  void set_fitness(double fitness);
  void set_fitness(bool maximize, std::vector<double>& fitness_vec);
  void set_fitness(std::vector<double>& fitness_vec);

  int get_num_evaluations();
  void incr_num_evaluations();

  int get_age() const;
  int get_age(int evals, int num_individs) const;
  void set_age(int age);
  void increase_age();

  int get_age_move() const;
  void set_age_move(int age);

  void set_take_age_older();
  void set_take_age_average();
  void set_take_age_younger();

  void set_creation(int type);
  int get_creation();
  
  void print_history();

  virtual void make_random();
  virtual void duplicate_settings(Individual *individ2);
  virtual void duplicate(Individual *ind);

  virtual bool mutate();
  virtual bool mutate(std::vector<Individual*>& ind2s);

  virtual bool recombine(Individual *parent2);
  virtual bool recombine_rand2(Individual *parent2);
  virtual bool recombine(std::vector<Individual*>& ind2s);

  virtual int phenotype_distance(Individual* ind2);
  virtual bool make_phenotype(void* arg);
  virtual bool make_phenotype();

  void write_log(char *string);
  std::ostream& write_log(std::ostream& log_stream);

  virtual std::istream& read(std::istream& istr);
  virtual bool read(const char *fname);
  virtual std::ostream& write(std::ostream& ostr);
  virtual bool write(const char *fname);

  void print_history_ops();
  void print_results();
  void print_history_full();

private:
  int evaluations_; // Count of number of times evaluated.
  int age_move_; // For steady-state ALPS: try_move_up()
  int assign_age_type_;

  IndHistory history_self_;
  IndHistory history_parent1_;
  IndHistory history_parent2_;
};

bool individ_compare_max(void const *elt1, void const *elt2);
bool individ_compare_min(const void *elt1, const void *elt2);

void individual_print_op_history();
void individual_print_op_history(char *string);
void individual_write_op_history(char *string);
void individual_read_op_history(const char *string);

}



#endif
