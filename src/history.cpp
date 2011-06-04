/*******************************************************

  File: history.cpp
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

#include "alps/history.h"
#include "alps/utils.h"
using namespace alps_utils;

#include <iostream>
#include <string>
#include <fstream>
using namespace std;


namespace alps {

const int Alps_Init_Age = 1;


ostream& operator<<(ostream& os, IndHistory& hist)
{
  return hist.write(os);
}


istream& operator>>(istream& is, IndHistory& hist)
{
  //  return hist.read(is);

  is >> hist.age_;
  is >> hist.num_mutate_;
  is >> hist.num_recombine_;
  is >> hist.how_created_;
  is >> hist.last_op_;

  int fit_size = 0;
  is >> fit_size;

  int fit_prev_size = 0;
  is >> fit_prev_size;

  int op_count_size = 0;
  is >> op_count_size;

  int complex_size = 0;
  is >> complex_size;

  hist.fitness_.clear();
  for (int i=0; i<fit_size; i++) {
    double val;
    is >> val;
    hist.fitness_.push_back(val);
  }

  int ctr=0;
  while ((is.get() != '\n') && (ctr < 10)) {
    ctr++;
  }

  hist.fitness_prev_.clear();
  for (int i=0; i<fit_prev_size; i++) {
    double val;
    is >> val;
    hist.fitness_prev_.push_back(val);
  }
  ctr=0;
  while ((is.get() != '\n') && (ctr < 10)) {
    ctr++;
  }

  hist.op_count_.clear();
  for (int i=0; i<op_count_size; i++) {
    int val;
    is >> val;
    hist.op_count_.push_back(val);
  }
  ctr=0;
  while ((is.get() != '\n') && (ctr < 10)) {
    ctr++;
  }

  hist.complexity_.clear();
  for (int i=0; i<complex_size; i++) {
    double val;
    is >> val;
    hist.complexity_.push_back(val);
  }
  ctr=0;
  while ((is.get() != '\n') && (ctr < 10)) {
    ctr++;
  }

  return is;
}



IndHistory :: IndHistory()
  : age_(Alps_Init_Age), how_created_(CREATE_INIT),
    num_mutate_(0), num_recombine_(0)
{
  op_count_.resize(1,0);
  fitness_.resize(1,0.0);
  fitness_prev_.resize(1,0.0);
  complexity_.resize(1,0.0);
}

IndHistory :: ~IndHistory()
{

}


void IndHistory :: initialize()
{
  age_ = Alps_Init_Age;
  how_created_ = CREATE_INIT;
  num_mutate_ = 0;
  num_recombine_ = 0;

  /*
  op_count_.resize(0);
  fitness_.resize(0);
  fitness_prev_.resize(0);
  complexity_.resize(0);  
  */


  op_count_.clear();
  op_count_.resize(1, 0);

  fitness_.clear();
  fitness_.resize(1, 0.0);

  fitness_prev_.clear();
  fitness_prev_.resize(1, 0.0);

  complexity_.clear();
  complexity_.resize(1, 0.0);
}



void IndHistory :: backup_fitness()
{
  fitness_prev_ = fitness_;
}

void IndHistory :: avg_fitness(vector<double>& fvec)
{
  if (fitness_.size() != fvec.size()) {
    throw AlpsError("history :: avg_fitness() - error, fitness vectors not same size.");
  }

  for (vector<double>::size_type i = 0; i<fitness_.size(); i++) {
    fitness_[i] = (fitness_[i] + fvec[i]) / 2.0;
  }
}

void IndHistory :: set_fitness(vector<double>& fvec)
{
  fitness_ = fvec;
}



double IndHistory :: get_fitness() const
{
  return fitness_[0];
}

double IndHistory :: get_fitness(int index) const
{
  return fitness_[index];
}

vector<double>& IndHistory :: get_fitness_vec()
{
  return fitness_;
}


double IndHistory :: get_fitness_prev() const
{
  return fitness_prev_[0];
}


double IndHistory :: get_complexity(int index) const
{
  return complexity_[index];
}


int IndHistory :: get_how_created() const
{
  return how_created_;
}

void IndHistory :: set_how_created(int h)
{
  how_created_ = h;
}


int IndHistory :: get_age() const
{
  return age_;
}

int IndHistory :: get_age(int evals, int denom) const
{
  //  return (evals - age_)/denom + Alps_Init_Age;
  return (age_+evals)/denom + Alps_Init_Age;
}

void IndHistory :: set_age(int a)
{
  age_ = a;
}

void IndHistory :: incr_age()
{
  age_++;
}

void IndHistory :: take_older(const IndHistory& h2)
{
  if (h2.age_ > age_) {
    age_ = h2.age_;
  }
}

void IndHistory :: take_younger(const IndHistory& h2)
{
  if (h2.age_ < age_) {
    age_ = h2.age_;
  }
}


void IndHistory :: take_average(const IndHistory& h2)
{
  age_ = (age_ + h2.age_) / 2;
}


void IndHistory :: incr_num_mutate()
{
  num_mutate_++;
}

void IndHistory :: incr_num_recombine()
{
  num_recombine_++;
}



vector<int>* IndHistory :: get_op_count()
{
  return &op_count_;
}


void IndHistory :: print_created_orig(ostream& ofs)
{
  if (how_created_ == CREATE_RANDOM) {
    ofs << "[*:";
  } else if (how_created_ == CREATE_INIT) {
    ofs << "[-:";
  } else {
    if (last_op_ != -1) {
      ofs << "[" << last_op_ << ":";
    } else {
      if (how_created_ == CREATE_MUTATE)
	ofs << "[m:";
      else if (how_created_ == CREATE_RECOMBINE)
	ofs << "[r:";
      else 
	ofs << "[?:";
    }
  }

  ofs << num_mutate_ << ":" << num_recombine_ << "]";
}


void IndHistory :: print_created(ostream& ofs)
{
  if (how_created_ == CREATE_RANDOM) {
    ofs << "[*(" << age_ << "):";
  } else if (how_created_ == CREATE_INIT) {
    ofs << "[-(" << age_ << "):";
  } else {
    if (how_created_ == CREATE_MUTATE)
      ofs << "[m(" << age_ << "):";
    else if (how_created_ == CREATE_RECOMBINE)
      ofs << "[r(" << age_ << "):";
    else 
      ofs << "[?(" << age_ << "):";
  }

  ofs << num_mutate_ << ":" << num_mutate_ << "]";
}


void IndHistory :: print_ops(ostream& ofs)
{
  if (op_count_.size() == 0)
    return;

  ofs << "[";
  for (vector<int>::size_type i=0; i<op_count_.size()-1; i++) {
    ofs << op_count_[i] << " ";
  }

  ofs << op_count_[op_count_.size()-1] << "]";
}


void IndHistory :: print_fitnesses(ostream& ofs)
{
  if (fitness_.size() == 1) {
    ofs << fitness_[0] << " ";
  } else {
    for (vector<double>::const_iterator iter = fitness_.begin();
	 iter != fitness_.end(); ++iter) {
      ofs << (*iter) << " ";
    }
  }
}

void IndHistory :: print(ostream& ofs)
{
  print_fitnesses(ofs);
  print_created(ofs);
}

void IndHistory :: print_f1(ostream& ofs)
{
  ofs << fitness_[0] << " ";
  print_created(ofs);
}


istream& IndHistory :: read(istream& is)
{
  is >> age_;
  is >> num_mutate_;
  is >> num_recombine_;
  is >> how_created_;
  is >> last_op_;

  int fit_size = 0;
  is >> fit_size;

  int fit_prev_size = 0;
  is >> fit_prev_size;

  int op_count_size = 0;
  is >> op_count_size;

  int complex_size = 0;
  is >> complex_size;


  fitness_.clear();
  for (int i=0; i<fit_size; i++) {
    double val;
    is >> val;
    fitness_.push_back(val);
  }
  
  fitness_prev_.clear();
  for (int i=0; i<fit_prev_size; i++) {
    double val;
    is >> val;
    fitness_prev_.push_back(val);
  }

  op_count_.clear();
  for (int i=0; i<op_count_size; i++) {
    int val;
    is >> val;
    op_count_.push_back(val);
  }

  complexity_.clear();
  for (int i=0; i<complex_size; i++) {
    int val;
    is >> val;
    complexity_.push_back(val);
  }

  return is;
}


ostream& IndHistory :: write(ostream& ostr)
{
  ostr << age_ << " ";
  ostr << num_mutate_ << " ";
  ostr << num_recombine_ << " ";
  ostr << how_created_ << " ";
  ostr << last_op_ << " ";
  ostr << fitness_.size() << " ";
  ostr << fitness_prev_.size() << " ";
  ostr << op_count_.size() << " ";
  ostr << complexity_.size();
  ostr << endl;

  {
    vector<double>::const_iterator iter = fitness_.begin();
    ostr << (*iter);
    for (++iter; iter != fitness_.end(); ++iter) {
      ostr << " " << (*iter);
    }
    ostr << endl;
  }

  {
    vector<double>::const_iterator iter = fitness_prev_.begin();
    ostr << (*iter);
    for (++iter; iter != fitness_prev_.end(); ++iter) {
      ostr << " " << (*iter);
    }
    ostr << endl;
  }

  if (!op_count_.empty()) {
    vector<int>::const_iterator iter = op_count_.begin();
    ostr << (*iter);
    for (++iter; iter != op_count_.end(); ++iter) {
      ostr << " " << (*iter);
    }
    ostr << endl;
  }

  if (!complexity_.empty()) {
    vector<double>::const_iterator iter = complexity_.begin();
    ostr << (*iter);
    for (++iter; iter != complexity_.end(); ++iter) {
      ostr << " " << (*iter);
    }
    ostr << endl;
  }

  return ostr;
}



}
