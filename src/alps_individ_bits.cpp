/*******************************************************

  File: individ_bits.cpp
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


#include "alps_utils.h"
using namespace alps_utils;

#include "alps_random_mt.h"
using namespace alps_random;

#include "alps_individ_bits.h"


namespace alps {


const int INDIVID_BITS_VERSION = 1;

double Prob_Inversion = 0.0;
double Prob_HTS = 0.0;
double Prob_Invert_HTS = 0.0;

/***********************************************************/

// Individual class functions.

Individ_Bits :: Individ_Bits(int num)
  : num_genes_(num), is_shuffled_(false)
{
  gene_vec_.resize(num);
}

Individ_Bits :: Individ_Bits(Individual *sample)
{
  if (sample == 0) {
    return;
  }

  Individ_Bits *s = (Individ_Bits*)sample;
  num_genes_ = s->num_genes_;
  gene_vec_.resize(num_genes_);

  is_shuffled_ = s->is_shuffled_;
  gene_loc_ = s->gene_loc_;
}


Individ_Bits :: ~Individ_Bits()
{
}


Individ_Bits* Individ_Bits :: new_instance()
{
  Individ_Bits* n = new Individ_Bits(gene_vec_.size());
  n->duplicate_settings(this);
  return n;
}


// Initialize variables.
void Individ_Bits :: initialize(void)
{
  Individual::initialize();
}


int Individ_Bits :: same(Individ_Bits *individ2)
{
  int result = true;

  if (gene_vec_.size() != individ2->gene_vec_.size()) {
    cout << "individual :: same() - number of genes different "
	 << gene_vec_.size() << " != " << individ2->gene_vec_.size() << endl;
    result = false;
  }

  if (!result) {
    while (1) ;
  }

  return result;
}


void Individ_Bits :: zero_genes(void)
{
  for (unsigned int i=0; i<gene_vec_.size(); i++) {
    gene_vec_[i] = 0;
  }
}  


int Individ_Bits :: get_num_genes()
{
  return gene_vec_.size();
}



bool Individ_Bits :: set_num_genes(int num)
{
  if (gene_vec_.size() <= 0) {
    return false;
  }

  num_genes_ = num;
  gene_vec_.resize(num);

  return true;
}

const vector<char>& Individ_Bits :: get_genes()
{
  return gene_vec_;
}



void Individ_Bits :: shuffle()
{
  // It is probably better to handle shuffling outside of the class.

  is_shuffled_ = true;

  gene_loc_.resize(gene_vec_.size());
  for (vector<char>::size_type i=0; i<gene_vec_.size(); i++) {
    gene_loc_[i] = i;
  }

}


void Individ_Bits :: make_random()
{
  Individual::make_random();

  for (vector<char>::size_type i=0; i<gene_vec_.size(); i++) {
    gene_vec_[i] = random_int(2);
  }
}


void Individ_Bits :: duplicate_settings(Individ_Bits *individ2)
{
  Individual::duplicate_settings(individ2);
  num_genes_ = individ2->num_genes_;
  gene_vec_.resize(individ2->gene_vec_.size());

  is_shuffled_ = individ2->is_shuffled_;
  gene_loc_ = individ2->gene_loc_;
}
void Individ_Bits :: duplicate_settings(Individual *individ2)
{
  duplicate_settings((Individ_Bits*)individ2);
}



void Individ_Bits :: duplicate(Individ_Bits *individ2)
{
  Individual::duplicate(individ2);
  gene_vec_ = individ2->gene_vec_;

  is_shuffled_ = individ2->is_shuffled_;
  gene_loc_ = individ2->gene_loc_;  
}
void Individ_Bits :: duplicate(Individual *ind)
{
  duplicate((Individ_Bits*)ind);
}


bool Individ_Bits :: mutate_virtual_tree(void)
{
  // Virtual Sub-Tree Mutation.
  // - works by picking a random subtree and flipping the bits.
  int tree_point = random_int(num_hts_points_);

  int level_pts = gene_vec_.size();
  int level = 1;
  int level_size = 1;
  while (tree_point >= level_pts) {
    tree_point -= level_pts;
    level_pts = level_pts >> 1;
    level++;
    level_size = level_size << 1;
  }
  
  int start_pt = level_size*tree_point;
  int stop_pt = start_pt + level_size;
  if (stop_pt > (int)gene_vec_.size()) {
    stop_pt = gene_vec_.size();
  }

  for (int i=start_pt; i<stop_pt; i++) {
    if (gene_vec_[i] == 1) {
      gene_vec_[i] = 0;
      } else {
      gene_vec_[i] = 1;
    }
  }

  return true;
}


bool Individ_Bits :: mutate()
{
  int mutated = false;

  //  write("mutating.ind");

  //  printf("ibits :: mutate()\n");

  if (gene_vec_.size() == 0) {
    return false;
  }

  Individual::mutate();


  {
    // Flips one randomly selected bit:
    //    int num = 1;
    int num = random_int(5) + 1;
    for (int i=0; i<num; i++) {
      int index = random_int(gene_vec_.size());
      if (gene_vec_[index] == 1) {
	gene_vec_[index] = 0;
      } else {
	gene_vec_[index] = 1;
      }
      mutated = true;
    }
  }

  /*
  {
    // Probability-based system that tests to flip every bit:
    bool is_flip = false;
    double prob_mut = 2.0 / double(gene_vec_.size());
    for (vector<char>::size_type i=0; i<gene_vec_.size(); i++) {
      if (random_double() < prob_mut) {
	is_flip = true;
	gene_vec_[i] = (gene_vec_[i]==0 ? 1 : 0);
      }
    }

    if (!is_flip) {
      int i = random_int(gene_vec_.size());
      gene_vec_[i] = (gene_vec_[i]==0 ? 1 : 0);
    }
  }
  */

  /*
  if (random_double() >= Prob_HTS) {
    int num = random_int(4) + 1;
    for (int i=0; i<num; i++) {
      int index = random_int(gene_vec_.size());
      if (gene_vec_[index] == 1) {
	gene_vec_[index] = 0;
      } else {
	gene_vec_[index] = 1;
      }
      mutated = true;
    }

  } else {
    mutated = mutate_virtual_tree();
  }
  */

  return mutated;
}



bool Individ_Bits :: recombine_tree(Individ_Bits *parent2)
{
  int recomb_point = random_int(num_hts_points_);
  //  int recomb_pt2 = recomb_point;

  int level_pts = gene_vec_.size();
  int level = 1;
  int level_size = 1;
  //  printf(" level_pts:%d, size:%d  recomb_point:%d.\n", level_pts, level_size, recomb_point);
  while (recomb_point >= level_pts) {
    recomb_point -= level_pts;
    level_pts = level_pts >> 1;
    level++;
    level_size = level_size << 1;
    //    printf(" level_pts:%d, size:%d  recomb_point:%d.\n", level_pts, level_size, recomb_point);
  }
  //  printf(" level_pts:%d, size:%d  recomb_point:%d.\n", level_pts, level_size, recomb_point);
  
  int start_pt = level_size*recomb_point;
  int stop_pt = start_pt + level_size;
  if (stop_pt > (int)gene_vec_.size()) {
    stop_pt = gene_vec_.size();
  }

  /*
  printf("recombine_tree() - #Genes:%d, #levels:%d, #pts:%d;  recomb-pt:%d(%d), level:%d, level-size:%d, start:%d, stop:%d.\n",
	 gene_vec_.size(), num_levels, num_points, recomb_pt2, recomb_point,
	 level, level_size, start_pt, stop_pt);
  */

  for (int i=start_pt; i<stop_pt; i++) {
    gene_vec_[i] = parent2->gene_vec_[i];
  }

  /*
  for (int i=start_pt; i<stop_pt; i++) {
    printf("%d ", (int)gene_vec_[i]);
  }
  printf("\n");

  for (int i=start_pt; i<stop_pt; i++) {
    printf("%d ", (int)parent2->gene_vec_[i]);
  }
  printf("\n\n");
  */  

  return true;
}



bool Individ_Bits :: recombine(Individ_Bits *parent2)
{
  if (gene_vec_.size() < 2) {
    return false;
  }

  Individual::recombine(parent2);

  //  recombine_tree(parent2);

  /*
  {
    bool is_change = false;
    double prob_rec = 0.5;
    for (vector<char>::size_type i=0; i<gene_vec_.size(); i++) {
      if (random_double() < prob_rec) {
	if (gene_vec_[i] != parent2->gene_vec_[i]) {
	  is_change = true;
	  gene_vec_[i] = parent2->gene_vec_[i];
	}
      }
    }

    if (!is_change) {
      int i = random_int(gene_vec_.size());
      gene_vec_[i] = (gene_vec_[i]==0 ? 1 : 0);
    }
  }
  */


  {
    unsigned int loc1 = random_int(gene_vec_.size());
    unsigned int loc2 = random_int(gene_vec_.size()-1);

    if (loc2 == loc1) {
      loc2 = gene_vec_.size() - 1;
    }

    unsigned int i=loc1;
    while (i != loc2) {
      gene_vec_[i] = parent2->gene_vec_[i];
      i++;
      if (i >= gene_vec_.size()) {
	i = 0;
      }
    }
  }
  


  /*
  if (random_double() >= Prob_HTS) {
    if (gene_vec_.size() == 2) {
      if (random_int(2)) {
	gene_vec_[0] = parent2->gene_vec_[0];
      }
      if (random_int(2)) {
	gene_vec_[1] = parent2->gene_vec_[1];
      }

    } else {
      unsigned int loc1 = random_int(gene_vec_.size());
      unsigned int loc2 = random_int(gene_vec_.size()-1);

      if (loc2 == loc1) {
	loc2 = gene_vec_.size() - 1;
      }

      unsigned int i=loc1;
      while (i != loc2) {
	gene_vec_[i] = parent2->gene_vec_[i];
	i++;
	if (i >= gene_vec_.size()) {
	  i = 0;
	}
      }
    }

  } else {
    recombine_tree(parent2);
  }
  */

  return true;
}

bool Individ_Bits :: recombine(Individual *parent2)
{
  return recombine((Individ_Bits*)parent2);
}


bool Individ_Bits :: recombine(vector<Individual*>& p2s)
{
  if (p2s.size() == 0) {
    return false;
  }

  return recombine(p2s[0]);
}


int Individ_Bits :: phenotype_distance(Individ_Bits* ind2)
{
  int num_genes = min(gene_vec_.size(), ind2->gene_vec_.size());

  int dist = max(gene_vec_.size(), ind2->gene_vec_.size()) - num_genes;
  for (int i=0; i<num_genes; i++) {
    if (gene_vec_[i] != ind2->gene_vec_[i]) {
      num_genes++;
    }
  }

  return dist;
}
int Individ_Bits :: phenotype_distance(Individual* ind2)
{
  return phenotype_distance((Individ_Bits*)ind2);
}


bool Individ_Bits :: make_phenotype(void*)
{
  return true;
}

bool Individ_Bits :: make_phenotype(void)
{
  return true;
}


std::istream& Individ_Bits :: read(std::istream& istr)
{
  int version = 0;
  string s;
  istr >> s >> version;
  if ((s != "IndBitsV") || (version != INDIVID_BITS_VERSION)) {
    //    throw AlpsError("IndBits :: read() - error, bad version.");
    cout << "IndBits :: read() - error, bad version.\n";
  }

  Individual::read(istr);

  istr >> num_genes_;
  
  cout << "Read: " << num_genes_ << endl;
  gene_vec_.clear();
  for (int i=0; i<num_genes_; i++) {
    int val;
    istr >> val;
    gene_vec_.push_back(val);
  }

  return istr;


#ifdef OLD
  int i, j, num, version;

  // if (Created) {
  //    uncreate();
  // }


  file_read(file, gFile_Buffer, gFILE_BUFFER_LEN);
  num = sscanf(gFile_Buffer, "Individ_Bits Version: %d\n", &version);
  if (num != 1) {
    printf("Individ_Bits: error in reading version\n");
    return false;
  } else if (version != INDIVID_BITS_VERSION) {
    printf("Individ_Bits: Error reading unsupported version #%d\n", version);
    return false;
  }    

  Evaluated = 2;
  
  file_read(file, gFile_Buffer, gFILE_BUFFER_LEN);  
  num = sscanf(gFile_Buffer, "%d\n", &num_genes_);
  if (num != 1) {
    printf("Individ_Bits: error in reading num_genes_\n");
    printf("<%s>\n", gFile_Buffer);
    return false;
  }

  if (num_genes_ < 1) {
    printf("Individ_Bits: error in reading num_genes_: %d < 1.\n",
	   num_genes_);
    printf("<%s>\n", gFile_Buffer);
    return false;
  }    

  create();


  if (!file_read(file, gFile_Buffer, gFILE_BUFFER_LEN)) {
    printf("Individ_Bits :: read() - error, empty line.");
    return false;
  }

  j = 0;
  for (i=0; i<num_genes_; i++) {
    gene_vec_[i] = read_integer(&(gFile_Buffer[j]));
    while (gFile_Buffer[j] != ' ')
      j++;
    j++;
  }

  if (!file_read(file, gFile_Buffer, gFILE_BUFFER_LEN)) {
    // end of data file.
    return true;
  }

#endif

  return istr;
}


bool Individ_Bits :: read(const char *fname)
{
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


ostream& Individ_Bits :: write(std::ostream& ostr)
{
  ostr << "IndBitsV " << INDIVID_BITS_VERSION << endl;

  Individual::write(ostr);

  if (gene_vec_.size() <= 0) {
    ostr << "IndReal::write()  Error - invalid size: " << gene_vec_.size() << endl;
    return ostr;
  }

  // Genes:
  ostr << gene_vec_.size() << endl;
  for (vector<char>::const_iterator iter = gene_vec_.begin();
       iter != gene_vec_.end(); ++iter) {
    ostr << (int)*iter << " ";
  }
  ostr << endl;

  return ostr;


#ifdef OLD
  int i, loc;

  sprintf(gFile_Buffer, "Individ_Bits Version: %d\n", INDIVID_BITS_VERSION);
  file_write(file, gFile_Buffer);

  sprintf(gFile_Buffer, "%d\n", gene_vec_.size());
  file_write(file, gFile_Buffer);

  if (gene_vec_.size() <= 0) {
    printf("Ind::write()  Error - size is %d.\n", gene_vec_.size());
    return;
  }

  loc = 0;
  for (i=0; i<gene_vec_.size(); i++) {
    //    sprintf(gFile_Buffer+loc, "%d ", gene_vec_[i]);
    //    while (gFile_Buffer[++loc] != ' ')
    //      ;
    sprintf(gFile_Buffer+loc, "%d", gene_vec_[i]);
    loc++;
  }
  gFile_Buffer[loc-1] = '\n';
  gFile_Buffer[loc] = '\n';
  gFile_Buffer[loc+1] = 0;
  file_write(file, gFile_Buffer);

  loc = 0;
  for (i=0; i<gene_vec_.size(); i++) {
    sprintf(gFile_Buffer+loc, "%d", gene_vec_[location_[i]]);
    //    while (gFile_Buffer[++loc] != ' ')
    //      ;
    loc++;
  }
  gFile_Buffer[--loc] = '\n';
  file_write(file, gFile_Buffer);


  if (Bit_Value) {
    sprintf(gFile_Buffer, "DeRandomized: ");
    loc = 14;
    for (i=0; i<gene_vec_.size(); i++) {
      sprintf(gFile_Buffer+loc, "%d", gene_vec_[i] ^ Bit_Value[i]);
      //    while (gFile_Buffer[++loc] != ' ')
      //      ;
      loc++;
    }
    gFile_Buffer[--loc] = '\n';
    file_write(file, gFile_Buffer);
  }

  if (inversion_) {
    sprintf(gFile_Buffer, "UnShuffled: ");
    loc = 12;
    for (i=0; i<gene_vec_.size(); i++) {
      sprintf(gFile_Buffer+loc, "%d", gene_vec_[Bit_Order[location_[i]]]);
      //    while (gFile_Buffer[++loc] != ' ')
      //      ;
      loc++;
    }
    gFile_Buffer[--loc] = '\n';
    file_write(file, gFile_Buffer);

  } else if (Bit_Order) {
    sprintf(gFile_Buffer, "UnShuffled: ");
    loc = 12;
    for (i=0; i<gene_vec_.size(); i++) {
      sprintf(gFile_Buffer+loc, "%d", gene_vec_[Bit_Order[i]]);
      //    while (gFile_Buffer[++loc] != ' ')
      //      ;
      loc++;
    }
    gFile_Buffer[--loc] = '\n';
    file_write(file, gFile_Buffer);
  }


  loc = 0;
  for (i=0; i<gene_vec_.size(); i++) {
    sprintf(gFile_Buffer+loc, "%d ", location_[i]);
    while (gFile_Buffer[++loc] != ' ')
      ;
    loc++;
  }
  gFile_Buffer[--loc] = '\n';
  file_write(file, gFile_Buffer);

  if (Bit_Order) {
    loc = 0;
    for (i=0; i<gene_vec_.size(); i++) {
      sprintf(gFile_Buffer+loc, "%d ", Bit_Order[i]);
      while (gFile_Buffer[++loc] != ' ')
	;
      loc++;
    }
    gFile_Buffer[--loc] = '\n';
    file_write(file, gFile_Buffer);
  }
#endif

  return ostr;
}


bool Individ_Bits :: write(const char *fname)
{
  ofstream file(fname);
  write(file);
  return true;
}



} // namespace alps

/********************************************************/
