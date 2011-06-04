/*******************************************************

  File: individ_real.cpp
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


#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;


#include "alps/random_mt.h"
using namespace alps_random;

#include "alps/utils.h"
using namespace alps_utils;



#include "alps/individ_real.h"



namespace alps {

#define OP_MUTATE_ALL 0
#define OP_MUTATE_SOME 1
#define OP_RECOMB_ALL 2
#define OP_RECOMB_SOME 3
#define OP_RECOMB_SWAP 4



#define INDIVID_REAL_VERSION 1

/***********************************************************/

// Individual class functions.

Individ_Real :: Individ_Real(int num)
  : num_genes_(num)
{
  gene_vec_.resize(num_genes_);
  max_val_.resize(num_genes_);
  min_val_.resize(num_genes_);
  init_max_.resize(num_genes_);
  init_min_.resize(num_genes_);
  set_minmax(0.0, 1.0);
  set_init_minmax(0.0, 1.0);
}


Individ_Real :: Individ_Real(Individual *sample)
{
  Individ_Real *s = (Individ_Real*)sample;

  num_genes_ = s->num_genes_;
  gene_vec_.resize(num_genes_);

  max_val_ = s->max_val_;
  min_val_ = s->min_val_;
  init_max_ = s->init_max_;
  init_min_ = s->init_min_;

  //  initialize();

  if (sample != 0) {
    num_genes_ = s->num_genes_;
  }
}


Individ_Real :: ~Individ_Real()
{
}


Individ_Real* Individ_Real :: new_instance(void)
{
  //  printf("individ_real :: new_instance() - %d\n", num_genes_);
  Individ_Real* n = new Individ_Real(num_genes_);
  n->duplicate_settings(this);
  return n;
}


// Initialize variables.
void Individ_Real :: initialize(void)
{
  Individual::initialize();
}


int Individ_Real :: same(Individ_Real *individ2)
{
  int result = true;

  if (num_genes_ != individ2->num_genes_) {
    printf("individual :: same() - number of genes different %d != %d.\n",
	   num_genes_, individ2->num_genes_);
    result = false;
  }

  if (!result) {
    while (1) ;
  }

  return result;
}


void Individ_Real :: set_minmax(vector<double>& min_val_s,
				vector<double>& max_val_s)
{
  min_val_ = min_val_s;
  max_val_ = max_val_s;
}


void Individ_Real :: set_minmax(double min, double max)
{
  for (vector<double>::iterator iter = min_val_.begin();
       iter != min_val_.end(); iter++) {
    *iter = min;
  }

  for (vector<double>::iterator iter = max_val_.begin();
       iter != max_val_.end(); iter++) {
    *iter = max;
  }
}


void Individ_Real :: set_init_minmax(double min, double max)
{
  for (vector<double>::iterator iter = init_min_.begin();
       iter != init_min_.end(); iter++) {
    *iter = min;
  }

  for (vector<double>::iterator iter = init_max_.begin();
       iter != init_max_.end(); iter++) {
    *iter = max;
  }
}


void Individ_Real :: set_init_minmax(vector<double>& mins,
				     vector<double>& maxs)
{
  init_min_ = mins;
  init_max_ = maxs;
}



void Individ_Real :: zero_genes(void)
{
  for (int i=0; i<num_genes_; i++) {
    gene_vec_[i] = 0;
  }
}  



int Individ_Real :: get_num_genes(void)
{
  return num_genes_;
}


bool Individ_Real :: set_num_genes(int num)
{
  if (num <= 0) {
    return false;
  }

  num_genes_ = num;
  gene_vec_.resize(num);
  max_val_.resize(num_genes_);
  min_val_.resize(num_genes_);
  init_max_.resize(num_genes_);
  init_min_.resize(num_genes_);
  set_init_minmax(0.0, 1.0);
  set_minmax(0.0, 1.0);

  return true;
}

const vector<double>& Individ_Real :: get_genes()
{
  return gene_vec_;
}



void Individ_Real :: make_random(void)
{
  Individual::make_random();

  for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
    gene_vec_[i] = random_double()*(init_max_[i]-init_min_[i]) + init_min_[i];
  }
}


void Individ_Real :: duplicate_settings(Individ_Real *ind)
{
  Individual::duplicate_settings(ind);
  num_genes_ = ind->num_genes_;
  gene_vec_.resize(num_genes_);
  min_val_ = ind->min_val_;
  max_val_ = ind->max_val_;
  init_min_ = ind->init_min_;
  init_max_ = ind->init_max_;  
}
void Individ_Real :: duplicate_settings(Individual *individ2)
{
  duplicate_settings((Individ_Real*)individ2);
}



void Individ_Real :: duplicate(Individ_Real *ind)
{
  Individual::duplicate(ind);
  gene_vec_ = ind->gene_vec_;
  min_val_ = ind->min_val_;
  max_val_ = ind->max_val_;
  init_min_ = ind->init_min_;
  init_max_ = ind->init_max_;
}
void Individ_Real :: duplicate(Individual *ind)
{
  duplicate((Individ_Real*)ind);
}



bool Individ_Real :: mutate(double scale)
{
  int mutated = false;

  Individual::mutate();


  // Decide how many genes to mutate.
  const int k = 5;
  int max = (k <= (int)gene_vec_.size() ? k : (int)gene_vec_.size());
  int num = 1;
  if (max > 2) {
    num = random_int(max) + 1;
  }

  // Perform the mutations.
  for (int i=0; i<num; i++) {
    int index = 0;
    if (num_genes_ > 1) {
      index = random_int(num_genes_);
    }      
    double dist = max_val_[index] - min_val_[index];
    gene_vec_[index] += random_norm()*random_norm()*scale*dist;
    if (gene_vec_[index] < min_val_[index]) {
      gene_vec_[index] = min_val_[index];
    } else if (gene_vec_[index] > max_val_[index]) {
      gene_vec_[index] = max_val_[index];
    }
    mutated = true;
  }

  return mutated;
}

bool Individ_Real :: mutate()
{
  return mutate(0.01);
}



void Individ_Real :: get_unique(vector<Individual*>& inds2,
				vector<Individ_Real*>& indreals)
{
  indreals.resize(0);

  for (vector<Individual*>::size_type i=0; i<inds2.size(); i++) {
    Individ_Real* ind2 = (Individ_Real*)inds2[i];

    vector<double>::size_type j=0;
    for (; j<gene_vec_.size(); j++) {
      if (gene_vec_[j] != ind2->gene_vec_[j]) {
	break;
      }
    }
    if (j < gene_vec_.size()) {
      indreals.push_back(ind2);
    }
  }
}


bool Individ_Real :: mutate(vector<Individual*>& inds2)
{
  if (inds2.size() == 0) {
    return mutate();
  }

  /*
  int num_gene_same = 0;
  for (vector<Individual*>::size_type i=0; i<inds2.size(); i++) {
    vector<double>::size_type j=0;
    for (; j<gene_vec_.size(); j++) {
      Individ_Real* ind2 = (Individ_Real*)inds2[i];
      if (gene_vec_[j] != ind2->gene_vec_[j]) {
	break;
      }
    }
    if (j == gene_vec_.size()) {
      num_gene_same++;
    }
  }
  if (num_gene_same > 0) {
    cout << "Mutation: #gene same: " << num_gene_same << endl;
  }
  */


  vector<Individ_Real*> p2s;
  get_unique(inds2, p2s);
  if (p2s.size() == 0) {
    //    cout << "Mutation :: error, all in tournament are the same.\n";
    return mutate();
  }


  Individ_Real* ind2 = (Individ_Real*)p2s[0];
  if (p2s.size() > 1) {
    ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
  }

  //  double msize = 0.75 * random_norm();
  bool mutated = false;
  while (!mutated) {
    for (int i=0; i<num_genes_; i++) {
      if (random_double() < 0.4) {
	continue;
      }
      double range = gene_vec_[i] - ind2->gene_vec_[i];
      double dist = random_norm();
      //    double dist = random_double()*2.0 - 1.0;
      //    double dist = random_norm()*random_norm()*5.0;
      gene_vec_[i] += range*dist;
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }
    mutated = true;
  }

  Individual::mutate();

  return true;
}




// Clips correctly.
void Individ_Real :: recombine_linear_clip(double dist, vector<double>& parent1,
					   vector<double>& parent2)
{
  if (dist < 0.0) {
    // Extrapolating past P1, so need to make sure we don't go out of bounds.
    for (vector<double>::size_type i=0; i<parent1.size(); i++) {
      double val = gene_vec_[i] + (parent2[i] - parent1[i]) * dist;
      if (val > max_val_[i]) {
	dist = (max_val_[i] - gene_vec_[i]) / (parent2[i] - parent1[i]);
      } else if (val < min_val_[i]) {
	dist = (min_val_[i] - gene_vec_[i]) / (parent2[i] - parent1[i]);
      }
    }

  } else if (dist > 1.0) {
    // Extrapolating past P2, so need to make sure we don't go out of bounds.
    for (vector<double>::size_type i=0; i<parent1.size(); i++) {
      double val = gene_vec_[i] + (parent2[i] - parent1[i]) * dist;
      if (val > max_val_[i]) {
	dist = (max_val_[i] - gene_vec_[i]) / (parent2[i] - parent1[i]);
      } else if (val < min_val_[i]) {
	dist = (min_val_[i] - gene_vec_[i]) / (parent2[i] - parent1[i]);
      }
    }
  }

  /*
  for (vector<double>::size_type i=0; i<parent1.size(); i++) {
    //    parent1[i] += (parent2[i] - parent1[i]) * dist;
    gene_vec_[i] += (parent2[i] - parent1[i]) * dist;
    if (gene_vec_[i] < min_val_[i] - 0.000000001) {
      cout << "gene_vec_[" << i << "] = " << gene_vec_[i]
	   << " less than min: " << min_val_[i] << endl;
      while (1) ;
    } else if (gene_vec_[i] > max_val_[i] + 0.0000000001) {
      cout << "gene_vec_[" << i << "] = " << gene_vec_[i]
	   << " greater than max: " << max_val_[i] << endl;
      while (1) ;
    }
  }
  */

  for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
    gene_vec_[i] += (parent2[i] - parent1[i]) * dist;

    // Needed to keep gene values inside min-max limits since the above
    // equations can suffer from precision errors.
    if (gene_vec_[i] < min_val_[i]) {
      gene_vec_[i] = min_val_[i];
    } else if (gene_vec_[i] > max_val_[i]) {
      gene_vec_[i] = max_val_[i];
    }
  }
}


void Individ_Real :: recombine_linear(double dist, vector<double>& parent1,
				      vector<double>& parent2)
{
  for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
    gene_vec_[i] += (parent2[i] - parent1[i]) * dist;
    if (gene_vec_[i] < min_val_[i]) {
      gene_vec_[i] = min_val_[i];
    } else if (gene_vec_[i] > max_val_[i]) {
      gene_vec_[i] = max_val_[i];
    }
  }
}


bool Individ_Real :: recombine_rand2(Individ_Real *parent2)
{
  if (num_genes_ < 1) {
    return false;
  }

  Individual::recombine(parent2);

  Op_used(OP_RECOMB_ALL);

  vector<double> tmp(gene_vec_);  
  double dist = 2.0 * random_double() - 1.0;
  recombine_linear(dist, gene_vec_, parent2->gene_vec_);
  for (int i=0; i<num_genes_; i++) {
    if (random_double() < 0.25) {
      gene_vec_[i] = tmp[i];
    }
  }

  return true;
}

bool Individ_Real :: recombine_rand2(Individual *parent2)
{
  //  return recombine((Individ_Real*)parent2);
  return recombine_rand2((Individ_Real*)parent2);
}



bool Individ_Real :: recombine(Individ_Real *parent2)
{
  if (num_genes_ < 1) {
    return false;
  }

  Individual::recombine(parent2);

  Op_used(OP_RECOMB_ALL);    

  Individ_Real* ind2 = parent2;
  double msize = 0.75 * random_norm();
  for (int i=0; i<num_genes_; i++) {
    double dist = gene_vec_[i] - ind2->gene_vec_[i];
    gene_vec_[i] += msize*dist;
    if (gene_vec_[i] < min_val_[i]) {
      gene_vec_[i] = min_val_[i];
    } else if (gene_vec_[i] > max_val_[i]) {
      gene_vec_[i] = max_val_[i];
    }
  }

  return true;
}

bool Individ_Real :: recombine(Individual *parent2)
{
  return recombine((Individ_Real*)parent2);
}



bool Individ_Real :: recombine(vector<Individual*>& inds2)
{
  if (num_genes_ < 1) {
    return false;
  }

  if (inds2.size() == 0) {
    return false;
  }

  Op_used(OP_RECOMB_ALL);    

  vector<Individ_Real*> p2s;
  get_unique(inds2, p2s);
  if (p2s.size() == 0) {
    //    cout << "Recombine :: error, all in tournament are the same.\n";
    return mutate();
  }

  int recomb_type = 1;
  /*
  if (random_double() < 0.5) {
    recomb_type = 6;
  }
  */

  if (recomb_type == 1) {
    // Recomb 1
    Individ_Real* ind2 = (Individ_Real*)p2s[0];
    if (p2s.size() > 1) {
      int index = random_int(p2s.size());
      ind2 = (Individ_Real*)p2s[index];
    }

    /*
    //    double msize = 0.25 * random_norm();
    double msize = random_double() * 2.0 - 1.0;
    do {
      msize = 0.75 * random_norm();
    } while (msize > 1.0);
    */
      
    double msize = random_norm();
    recombine_linear(msize, gene_vec_, ind2->gene_vec_);
    /*
    for (int i=0; i<num_genes_; i++) {
      double dist = ind1->gene_vec_[i] - ind2->gene_vec_[i];
      gene_vec_[i] += msize*dist;
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }
    */

  } else if (recomb_type == 2) {
    // Recomb 2
    vector<double> mut_range(gene_vec_.size());
    {
      Individ_Real* ind2 = (Individ_Real*)p2s[0];
      if (p2s.size() > 1) {
	ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
      }

      for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
	mut_range[i] = gene_vec_[i] - ind2->gene_vec_[i];
      }
    }


    Individ_Real* ind2 = (Individ_Real*)p2s[0];
    if (p2s.size() > 1) {
      ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
    }

    //    double msize = random_double() * 2.0 - 1.0;
    //    double msize = random_double() * 1.8 - 0.9;
    double msize = 0.75 * random_norm();
    recombine_linear(msize, gene_vec_, ind2->gene_vec_);

    for (int i=0; i<num_genes_; i++) {
      gene_vec_[i] += 0.01 * random_norm() * random_norm() * mut_range[i];
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }


  } else if (recomb_type == 3) {
    // Recomb 3
    Individ_Real* ind2 = (Individ_Real*)p2s[0];
    if (p2s.size() > 1) {
      ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
    }

    //    double msize = random_double() * 2.0 - 1.0;
    //    double msize = random_double() * 1.8 - 0.9;
    double msize = 1.0 * random_norm();
    recombine_linear(msize, gene_vec_, ind2->gene_vec_);
    /*
    for (int i=0; i<num_genes_; i++) {
      double dist = gene_vec_[i] - ind2->gene_vec_[i];
      gene_vec_[i] += msize*dist;
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }
    */

  } else if (recomb_type == 4) {
    // Recomb 4
    // First: randomly select two DIFFERENT parents:
    Individ_Real* ind1 = this;
    Individ_Real* ind2 = (Individ_Real*)p2s[0];
    if (p2s.size() > 1) {
      int index = random_int(p2s.size()+1);
      if (index != (int)p2s.size()) {
	ind1 = (Individ_Real*)p2s[index];
      } else {
	// ind1 = this. // Set above.
      }

      index = random_int(p2s.size());
      ind2 = (Individ_Real*)p2s[index];
      if (ind2 == ind1) {
	ind1 = this;
      }
    }

    //    double msize = 0.25 * random_norm();
    /*
    double msize = random_double() * 2.0 - 1.0;
    for (int i=0; i<num_genes_; i++) {
      double dist = ind1->gene_vec_[i] - ind2->gene_vec_[i];
      gene_vec_[i] += msize*dist;
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }
    */
    double msize = random_double() * 2.0 - 1.0;
    recombine_linear(msize, ind1->gene_vec_, ind2->gene_vec_);


  } else if (recomb_type == 5) {
    // Recomb 5
    vector<double> mut_range(gene_vec_.size());
    {
      Individ_Real* ind2 = (Individ_Real*)p2s[0];
      if (p2s.size() > 1) {
	ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
      }

      for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
	mut_range[i] = gene_vec_[i] - ind2->gene_vec_[i];
      }
    }


    vector<double> gradient(gene_vec_.size());
    for (vector<double>::size_type j=0; j<gradient.size(); j++) {
      gradient[j] = gene_vec_[j];
    }
    int ctr = p2s.size();
    double num_pos = 1.0;
    for (vector<Individual*>::size_type i=0; i<p2s.size(); i++) {
      ctr -= 2;
      if (ctr == 0) {
	continue;
      } else if (ctr > 1) {
	num_pos += 1.0;
      }

      double mult = 1.0;
      if (ctr < 0.0) {
	mult = -1.0;
      }
      Individ_Real* ind2 = (Individ_Real*)p2s[i];
      for (vector<double>::size_type j=0; j<gradient.size(); j++) {
	gradient[j] += ind2->gene_vec_[j] * mult;
      }
    }


    //    double div = p2s.size() * (p2s.size() + 1.0) / 2.0;
    //    double div = num_pos;
    /*
    for (vector<double>::size_type j=0; j<gradient.size(); j++) {
      gradient[j] /= div;
    }
    */

    Individ_Real* ind1 = this;
    /*
    int index = random_int(inds2.size()+1);
    if (index < inds2.size()) {
      ind1 = (Individ_Real*)inds2[index];
    }
    */

    double dist = 0.75 * random_norm();
    //    double dist = 2.0*random_double() - 1.0;
    for (int i=0; i<num_genes_; i++) {
      gene_vec_[i] = ind1->gene_vec_[i] + gradient[i]*dist;
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }

    for (int i=0; i<num_genes_; i++) {
      gene_vec_[i] += 0.01 * random_norm() * random_norm() * mut_range[i];
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }

  } else if (recomb_type == 6) {
    // Recomb 6
    vector<double> mut_range(gene_vec_.size());
    {
      Individ_Real* ind2 = (Individ_Real*)p2s[0];
      if (p2s.size() > 1) {
	ind2 = (Individ_Real*)p2s[random_int(p2s.size())];
      }

      for (vector<double>::size_type i=0; i<gene_vec_.size(); i++) {
	mut_range[i] = gene_vec_[i] - ind2->gene_vec_[i];
      }
    }


    double mult = p2s.size();
    vector<double> gradient(gene_vec_.size());
    for (vector<double>::size_type j=0; j<gradient.size(); j++) {
      gradient[j] = gene_vec_[j] * mult;
    }
    for (vector<Individual*>::size_type i=0; i<p2s.size(); i++) {
      mult -= 2.0;
      if (mult == 0.0) {
	continue;
      }

      Individ_Real* ind2 = (Individ_Real*)p2s[i];
      for (vector<double>::size_type j=0; j<gradient.size(); j++) {
	gradient[j] += ind2->gene_vec_[j] * (double)mult;
      }
    }
    double div = p2s.size() * (p2s.size() + 1.0) / 2.0;
    for (vector<double>::size_type j=0; j<gradient.size(); j++) {
      gradient[j] /= div;
      //      gradient[j] = gradient[j] * 2.0 / div;
    }


    /*
    for (vector<double>::size_type i=0; i<gene_vec_.size();
	 i++) {
      cout << gradient[i] << ":" << gene_vec_[i] - p2s[0]->gene_vec_[i]
	   << ", ";
    }
    cout << endl;
    */

    /*
    Individ_Real* ind1 = this;
    int index = random_int(inds2.size()+1);
    if (index < inds2.size()) {
      ind1 = (Individ_Real*)inds2[index];
    }
    */

    double dist = 0.75 * random_norm();
    //    double dist = random_norm();
    //    double dist = random_double()*2.0 - 1.0;

    for (int i=0; i<num_genes_; i++) {
      gene_vec_[i] += dist * gradient[i];
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }

    /*
    for (int i=0; i<num_genes_; i++) {
      gene_vec_[i] += 0.1 * random_norm() * random_norm() * mut_range[i];
      if (gene_vec_[i] < min_val_[i]) {
	gene_vec_[i] = min_val_[i];
      } else if (gene_vec_[i] > max_val_[i]) {
	gene_vec_[i] = max_val_[i];
      }
    }
    */

    // Decide how many genes to mutate.
    const int k = 5;
    int max = (k <= (int)gene_vec_.size() ? k : (int)gene_vec_.size());
    int num = 1;
    if (max > 2) {
      num = random_int(max-1) + 1;
    }    
    for (int i=0; i<num; i++) {
      int index = random_int(num_genes_);
      gene_vec_[index] += 0.02 * random_norm() * mut_range[index];
      if (gene_vec_[index] < min_val_[index]) {
	gene_vec_[index] = min_val_[index];
      } else if (gene_vec_[index] > max_val_[index]) {
	gene_vec_[index] = max_val_[index];
      }
    }
  }



  return true;
}



int Individ_Real :: phenotype_distance(Individ_Real* ind2)
{
  int num = min(num_genes_, ind2->num_genes_);

  int dist = max(num_genes_, ind2->num_genes_) - num;
  for (int i=0; i<num; i++) {
    if (gene_vec_[i] != ind2->gene_vec_[i]) {
      num++;
    }
  }

  return dist;
}
int Individ_Real :: phenotype_distance(Individual* ind2)
{
  return phenotype_distance((Individ_Real*)ind2);
}


  //bool Individ_Real :: make_phenotype(void* arg)
bool Individ_Real :: make_phenotype(void*)
{
  return true;
}

bool Individ_Real :: make_phenotype(void)
{
  return true;
}


istream& Individ_Real :: read(std::istream& file)
{
  int version = 0;
  string s;
  file >> s >> version;
  if ((s != "IndRealV") || (version != INDIVID_REAL_VERSION)) {
    //    throw AlpsError("IndReal :: read() - error, bad version.");
    //    throw AlpsError("IndReal :: read() - error, bad version.");
    cout << "IndReal :: read() - error, bad version.\n";
  }

  Individual::read(file);

  file >> num_genes_;
  
  cout << "Read: " << num_genes_ << endl;
  gene_vec_.clear();
  for (int i=0; i<num_genes_; i++) {
    double val;
    file >> val;
    gene_vec_.push_back(val);
  }

  return file;
}


bool Individ_Real :: read(const char *fname)
{
  cout << "Individ_Real :: read() " << fname << endl;
  ifstream file(fname);
  try {
    read(file);
    //  } catch (AlpsError e) {
  } catch(...) {
    return false;
  }

  file.close();
  return true;  
}


ostream& Individ_Real :: write(std::ostream& ostr)
{
  ostr << "IndRealV " << INDIVID_REAL_VERSION << endl;

  Individual::write(ostr);

  if (num_genes_ <= 0) {
    ostr << "IndReal::write()  Error - invalid size: " << num_genes_ << endl;
    return ostr;
  }

  // Genes:
  ostr << num_genes_ << endl;
  int p = ostr.precision();
  ostr.precision(10);
  for (vector<double>::const_iterator iter = gene_vec_.begin();
       iter != gene_vec_.end(); ++iter) {
    ostr << *iter << " ";
  }
  ostr.precision(p);
  ostr << endl;

  return ostr;
}


bool Individ_Real :: write(const char *fname)
{
  ofstream file(fname);
  write(file);
  return true;
}


} // namespace alps

/********************************************************/
