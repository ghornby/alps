/************************************************************
  Copyright Gregory S. Hornby.
  All Rights Reserved.
************************************************************/
#ifndef ALPS_HISTORY_HEADER_FILE
#define ALPS_HISTORY_HEADER_FILE

#include <iostream>
#include <fstream>
#include <vector>


namespace alps {

  const int CREATE_NOT = -1;

  const int CREATE_MUTATE = 1;
  const int CREATE_RECOMBINE = 2;
  const int CREATE_INIT = 3;

  const int CREATE_WAITING = 4;
  const int CREATE_RANDOM = 5;
  const int CREATE_MUTATION = 6;
  const int CREATE_RECOMBINATION = 7;
  const int CREATE_DUPLICATE = 8;



class IndHistory
{
  friend std::ostream& operator<<(std::ostream&, IndHistory&);
  friend std::istream& operator>>(std::istream&, IndHistory&);

 public:
  IndHistory();
  ~IndHistory();

  void initialize();


  void backup_fitness();
  void avg_fitness(std::vector<double>& fvec);
  void set_fitness(std::vector<double>& fvec);

  double get_fitness() const;
  double get_fitness(int index) const;
  std::vector<double>& get_fitness_vec();

  double get_fitness_prev() const;

  double get_complexity(int index) const;


  int get_how_created() const;
  void set_how_created(int h);

  int get_age() const;
  int get_age(int evals, int denom) const;
  void set_age(int a);
  void incr_age();
  void take_older(const IndHistory& h2);
  void take_younger(const IndHistory& h2);
  void take_average(const IndHistory& h2);

  void incr_num_mutate();
  void incr_num_recombine();

  std::vector<int>* get_op_count();


  void print_created_orig(std::ostream& ofs);
  void print_created(std::ostream& ofs);
  void print_ops(std::ostream& ofs);
  void print_fitnesses(std::ostream& ofs);
  void print(std::ostream& ofs);
  void print_f1(std::ostream& ofs);
  std::istream& read(std::istream& istr);
  std::ostream& write(std::ostream& ostr);

 private:
  int age_;
  int how_created_;
  int last_op_;
  int num_mutate_;
  int num_recombine_;

  std::vector<int> op_count_;
  std::vector<double> complexity_;
  std::vector<double> fitness_;
  std::vector<double> fitness_prev_;
};

}

#endif
