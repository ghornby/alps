/*******************************************************

  File: alps_layer.cpp
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

#include "alps/alps_layer.h"
#include "alps/alps_gen.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <string>
using namespace std;

namespace alps {

ostream& operator<<(ostream& os, AlpsLayer& l)
{
  return l.write(os);
}

istream& operator>>(istream& is, AlpsLayer& l)
{
  return l.read(is);
}




typedef vector<AlpsLayer>::size_type layer_sz;

AlpsLayer :: AlpsLayer() :
  elitism_(2),
  max_age_(-1),
  parent_layer_(-1),
  num_prev_layers_(-1),
  select_type_(ALPS_SELECT_TOURN),
  size_(0), // -1
  start_index_(0),
  last_index_(-1),
  stop_index_(0),
  tourn_size_(3),
  prob_select_prev_(0.0)
{
}


AlpsLayer :: AlpsLayer(int s) :
  elitism_(2),
  max_age_(-1),
  parent_layer_(-1),
  num_prev_layers_(-1),
  select_type_(ALPS_SELECT_TOURN),
  size_(s),
  start_index_(0),
  last_index_(-1),
  stop_index_(0),
  tourn_size_(3),
  prob_select_prev_(0.0)
{
}


AlpsLayer :: ~AlpsLayer()
{
}


void AlpsLayer :: set_elitism(unsigned int v)
{
  if (v > size_) {
    elitism_ = size_;
    return;
  }

  elitism_ = v;
}


void AlpsLayer :: set_parent_layer(int p)
{
  parent_layer_ = p;
}

void AlpsLayer :: set_max_age(int a)
{
  max_age_ = a;
}

void AlpsLayer :: set_prev_layers(int n, vector<int>& l)
{
  if (n < 0) {
    return;
  }

  num_prev_layers_ = n;
  prev_layers_.resize(0);

  for (vector<int>::const_iterator iter = l.begin();
       iter != l.end(); ++iter) {
    prev_layers_.push_back(*iter);
  }
}


void AlpsLayer :: set_select_type(int t)
{
  select_type_ = t;
}


void AlpsLayer :: set_size(int s)
{
  if (s < 0)
    return;

  size_ = s;

  stop_index_ = start_index_ + size_;
  last_index_ = stop_index_ - 1;
}


void AlpsLayer :: set_start_index(int i)
{
  if (i < 0)
    return;

  start_index_ = i;

  stop_index_ = start_index_ + size_;
  last_index_ = stop_index_ - 1;
}



void AlpsLayer :: set_prob_select_prev(double p)
{
  if (p < 0.0)
    return;
  if (p > 1.0)
    return;

  prob_select_prev_ = p;
}

void AlpsLayer :: set_tourn_size(unsigned int t)
{
  if (t > 0) {
    tourn_size_ = t;
    return;
  }
}


unsigned int AlpsLayer :: get_elitism() const
{
  return elitism_;
}

int AlpsLayer :: get_max_age() const
{
  return max_age_;
}

int AlpsLayer :: get_num_prev_layers() const
{
  return num_prev_layers_;
}

int AlpsLayer :: get_parent_layer() const
{
  return parent_layer_;
}

int AlpsLayer :: get_prev_layer(unsigned int l) const
{
  if (l >= prev_layers_.size()) {
    throw AlpsError("AlpsLayer :: get_prev_layer() - invalid number of layers.");
  }

  return prev_layers_[l];
}

double AlpsLayer :: get_prob_select_prev() const
{
  return prob_select_prev_;
}

unsigned int AlpsLayer :: get_size() const
{
  return size_;
}

int AlpsLayer :: get_select_type() const
{
  return select_type_;
}

unsigned int AlpsLayer :: get_start_index() const
{
  return start_index_;
}

int AlpsLayer :: get_last_index() const
{
  return last_index_;
}

unsigned int AlpsLayer :: get_stop_index() const
{
  return stop_index_;
}


unsigned int AlpsLayer :: get_tourn_size() const
{
  return tourn_size_;
}


void AlpsLayer :: print() const
{
  cout << start_index_ << " - " << size_ << ", ma:" << max_age_
       << " el:" << elitism_ << ", #prev:" << num_prev_layers_
       << "(" << prob_select_prev_ << ")" << endl;
}


void AlpsLayer :: print_info() const
{
  cout << "select:" << select_type_;
  cout << "  elitism:" << elitism_;
  cout << "  age:" << max_age_;
  cout << "  start:" << start_index_;
  cout << "  size:" << size_;
  cout << "  prob:" << prob_select_prev_;
  cout << "  parent:" << parent_layer_;
  cout << "  Prev(" << num_prev_layers_ << "):";

  for (vector<int>::const_iterator iter = prev_layers_.begin();
       iter != prev_layers_.end(); ++iter) {
    cout << " " << (*iter);
  }
  cout << endl;
}


istream& AlpsLayer :: read(istream& is)
{
  is >> start_index_ >> size_ >> max_age_ >> select_type_
     >> elitism_ >> parent_layer_ >> prob_select_prev_
     >> num_prev_layers_;

  stop_index_ = start_index_ + size_;
  last_index_ = stop_index_ - 1;


  prev_layers_.clear();
  for (int i=0; i<num_prev_layers_; i++) {
    int val;
    is >> val;
    prev_layers_.push_back(val);
  }

  return is;
}


ostream& AlpsLayer :: write(ostream& ostr)
{
  ostr << start_index_ << " " << size_ << " "
       << max_age_ << " " << select_type_ << " "
       << elitism_ << " " << parent_layer_ << " "
       << prob_select_prev_ << " " << num_prev_layers_ << endl;

  if (prev_layers_.size() > 0) {
    for (vector<int>::const_iterator iter = prev_layers_.begin();
	 iter != prev_layers_.end(); ++iter) {
      ostr << (*iter) << " " ;
    }
    ostr << endl;
  }

  return ostr;
}


ostream& AlpsLayer :: write_header(ostream& ostr)
{
  ostr << "# " << start_index_ << " " << size_ << " "
       << max_age_ << " " << select_type_ << " "
       << elitism_ << " " << parent_layer_ << " "
       << prob_select_prev_ << " " << num_prev_layers_ << " : ";

  if (prev_layers_.size() > 0) {
    for (vector<int>::const_iterator iter = prev_layers_.begin();
	 iter != prev_layers_.end(); ++iter) {
      ostr << (*iter) << " " ;
    }
  }
  ostr << endl;

  return ostr;
}

}
