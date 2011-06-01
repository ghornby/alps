/*******************************************************

  File: individual.cpp
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


#include "alps_utils.h"
using namespace alps_utils;

#include "alps_individ.h"


#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
using namespace std;


namespace alps {


const int INDIVIDUAL_VERSION = 1;

const int NUM_VARIATION_OPS = 10;


int gNum_Fitness = 5;
int gNum_Operators = NUM_VARIATION_OPS;
double gRec_delta = 1.0;
double gMutate_size = 1.0;


int gNum_Mutations = 0;
int gNum_Mutations_Equal = 0;
int gNum_Mutations_Better = 0;
int gNum_Recombinations = 0;
int gNum_Recombinations_Equal = 0;
int gNum_Recombinations_Better = 0;


int gOp_Used[NUM_HISTORY_OPS];
int gOp_Used_Hist[NUM_HISTORY_OPS];

const int OP_STACK_SIZE = 20;
int gOp_Stack[OP_STACK_SIZE];
int gOp_Stack_Ptr = 0;



/***********************************************************/

const int Age_Take_Older = 1;
const int Age_Take_Average = 2;
const int Age_Take_Younger = 3;


// Individual class functions.

Individual :: Individual()
  : assign_age_type_(Age_Take_Older)
{
}


Individual :: Individual(Individual *sample)
{
  assign_age_type_ = sample->assign_age_type_;
  printf("individual :: individual(individual* sample)\n");
}


Individual :: ~Individual()
{
}


Individual* Individual :: new_instance()
{
  return new Individual();
}


// Initialize variables.
void Individual :: initialize()
{
  history_self_.initialize();
  history_parent1_.initialize();
  history_parent2_.initialize();

  evaluations_ = 0;
  age_move_ = 0;
}


bool Individual :: valid()
{
  return true;
}


// fitness_hist_diff() - returns the different in fitness
// between this individual and its parent: (this fitness - parent fitness).
// parent is taken as the first parent.
double Individual :: fitness_hist_diff(int index)
{
  if ((index < 0) || (index >= gNum_Fitness)) {
    cout << "individual :: fitness_hist_diff() - error on fitness index: "
	 << index << "(0 to " << gNum_Fitness << ").\n";
    while (1) ;
  }

  return history_self_.get_fitness(index)
    - history_parent1_.get_fitness(index);
}

double Individual :: fitness_hist_diff()
{
  return (history_self_.get_fitness() - history_parent1_.get_fitness());
}



// Returns 1 if better fitness than parent(s)
// returns 0 if equal/inbetween parent(s)
// returns -1 if strictly less than parent(s).
int Individual :: better(bool maximize)
{
  if (maximize) {
    if ((history_self_.get_fitness() > history_parent1_.get_fitness()) &&
        (history_self_.get_fitness() > history_parent2_.get_fitness()))
      return 1;
    else if ((history_self_.get_fitness() == history_parent1_.get_fitness()) ||
             (history_self_.get_fitness() == history_parent2_.get_fitness()))
      return 0;
    else if ((history_self_.get_fitness() > history_parent1_.get_fitness()) ||
             (history_self_.get_fitness() > history_parent2_.get_fitness()))
      return 0;
    else
      return -1;
  }

    if ((history_self_.get_fitness() < history_parent1_.get_fitness()) &&
        (history_self_.get_fitness() < history_parent2_.get_fitness()))
      return 1;
    else if ((history_self_.get_fitness() == history_parent1_.get_fitness()) ||
             (history_self_.get_fitness() == history_parent2_.get_fitness()))
      return 0;
    else if ((history_self_.get_fitness() < history_parent1_.get_fitness()) ||
             (history_self_.get_fitness() < history_parent2_.get_fitness()))
      return 0;
    else
      return -1;
}


// Conditions for keeping this individual versus re-applying
// variation to the parent(s).
bool Individual :: keep(bool maximize)
{
  if (maximize) {
    if (history_self_.get_how_created() == CREATE_RANDOM) {
      if (history_self_.get_fitness() > 0.0) {
	/*
	  printf("individual :: keep() - init fitness %f > 0.0, keep.\n",
	  history_self_.get_fitness());
	*/
	return true;
      }
      /*
	printf("individual :: keep() - init fitness %f <= 0.0, don't keep.\n",
	history_self_.get_fitness());
      */
      return false;
    }
    
    if (history_self_.get_fitness() >= history_parent1_.get_fitness()*0.1) {
      /*
	printf("individual :: keep() - fitness %f > 0.1* %f keep.  (Created=%d).\n",
	history_self_.get_fitness(), history_parent1_.get_fitness(),
	history_self_.get_how_created());
      */
      return true;
    }

  } else {
    if (history_self_.get_how_created() == CREATE_RANDOM) {
      if (history_self_.get_fitness() < Worst_Fitness) {
	/*
	  printf("individual :: keep() - init fitness %f > 0.0, keep.\n",
	  history_self_.get_fitness());
	*/
	return true;
      }
      /*
	printf("individual :: keep() - init fitness %f <= 0.0, don't keep.\n",
	history_self_.get_fitness());
      */
      return false;
    }
    
    if (history_self_.get_fitness() <= history_parent1_.get_fitness()*100.0) {
      /*
	printf("individual :: keep() - fitness %f > 0.1* %f keep.  (Created=%d).\n",
	history_self_.get_fitness(), history_parent1_.get_fitness(),
	history_self_.get_how_created());
      */
      return true;
    }
  }


  /*
  printf("individual :: keep() - fitness %f < 0.1x %f, don't keep.\n",
	 history_self_.get_fitness(), history_parent1_.get_fitness());
  */

  return false;
}



// Returns 1 if this is better fitness than individ2
// returns 0 if this is equal to individ2
// returns -1 if this is strictly less than individ2.
int Individual :: compare_fitness(bool maximizing, Individual *individ2)
{
  //  double *fitness2 = individ2->history_self_.fitness;

  if (maximizing) {
    if (history_self_.get_fitness() > 
	individ2->history_self_.get_fitness()) {
      return 1; // This individual is better.
    }

  } else {
    if (history_self_.get_fitness() < 
	individ2->history_self_.get_fitness()) {
      return 1; // This individual is better.
    }
  }

  if (history_self_.get_fitness() == 
      individ2->history_self_.get_fitness()) {
    return 0; // Both are equal.
  }

  return -1; // The other individual is better.
}


//bool Individual :: same(Individual *individ2)
bool Individual :: same(Individual *)
{
  // To be implemented by subclass.
  return true;
}


bool Individual :: created_by_variation()
{
  if ((history_self_.get_how_created() != CREATE_MUTATE) &&
      (history_self_.get_how_created() != CREATE_RECOMBINE)) {
    return false;
  }

  return (evaluations_ == 1);
}


double Individual :: get_fitness()
{
  return history_self_.get_fitness();
}


void Individual :: get_fitness(vector<double>& fitness_vec)
{
  fitness_vec = history_self_.get_fitness_vec();
}

vector<double>& Individual :: get_fitness_vec()
{
  return history_self_.get_fitness_vec();
}



void Individual :: set_fitness(bool maximize, vector<double>& fitness_vec)
{
  int result;

  /*
  cout << "alps_individ() :: set_fit() " << fitness_vec[0]
       << " evals=" << evaluations_ << endl;
  */

  history_self_.backup_fitness();

  if (evaluations_ >= 1) {
    // Average fitness values.
    evaluations_++;
    history_self_.avg_fitness(fitness_vec);

  } else {
    evaluations_ = 1;
    history_self_.set_fitness(fitness_vec);
  }

  if (history_self_.get_how_created() == CREATE_MUTATE) {
    gNum_Mutations++;
    result = better(maximize);
    if (result == 1) {
      gNum_Mutations_Better++;
    } else if (result == 0) {
      gNum_Mutations_Equal++;
    }

  } else if (history_self_.get_how_created() == CREATE_RECOMBINE) {
    gNum_Recombinations++;
    result = better(maximize);
    if (result == 1) {
      gNum_Recombinations_Better++;
    } else if (result == 0) {
      gNum_Recombinations_Equal++;
    }
  }

  //  cout << "Op count size: " << history_self_.get_op_count()->size() << endl;
  Op_update(*history_self_.get_op_count());
}

void Individual :: set_fitness(vector<double>& fitness_vec)
{
  set_fitness(true, fitness_vec);
}


void Individual :: set_fitness(double fitness)
{
  vector<double> v(1, fitness);
  set_fitness(v);
}


int Individual :: get_num_evaluations()
{
  return evaluations_;
}

void Individual :: incr_num_evaluations()
{
  evaluations_++;
}


int Individual :: get_age() const
{
  return history_self_.get_age();
}

int Individual :: get_age(int evals, int num_individs) const
{
  return history_self_.get_age(evals, num_individs);
}

void Individual :: set_age(int age)
{
  history_self_.set_age(age);
}

void Individual :: increase_age()
{
  history_self_.incr_age();
}


int Individual :: get_age_move() const
{
  return age_move_;
}

void Individual :: set_age_move(int age)
{
  age_move_ = age;
}


void Individual :: set_take_age_older()
{
  assign_age_type_ = Age_Take_Older;
}
void Individual :: set_take_age_average()
{
  assign_age_type_ = Age_Take_Average;
}
void Individual :: set_take_age_younger()
{
  assign_age_type_ = Age_Take_Younger;
}




void Individual :: set_creation(int type)
{
  history_self_.set_how_created(type);
}

int Individual :: get_creation()
{
  return history_self_.get_how_created();
}


void Individual :: print_history()
{
  history_self_.print(cout);
}



void Individual :: make_random()
{
  evaluations_ = 0;
  history_self_.initialize();
  history_self_.set_how_created(CREATE_RANDOM);

  history_parent1_.initialize();
  history_parent2_.initialize();
}

void Individual :: duplicate_settings(Individual *individ2)
{
  assign_age_type_ = individ2->assign_age_type_;
}


void Individual :: duplicate(Individual *ind)
{
  evaluations_ = ind->evaluations_;
  history_self_ = ind->history_self_;
  assign_age_type_ = ind->assign_age_type_;

  history_parent1_ = ind->history_parent1_;
  history_parent2_ = ind->history_parent2_;
}



bool Individual :: mutate()
{
  evaluations_ = 0;

  history_parent1_ = history_self_;
  history_self_.incr_num_mutate();
  history_self_.set_how_created(CREATE_MUTATE);

  return true;
}

bool Individual :: mutate(vector<Individual*>&)
{
  return mutate();
}



bool Individual :: recombine(Individual *parent2)
{
  evaluations_ = 0;

  // History for this individual.
  history_parent1_ = history_self_;
  history_parent2_ = parent2->history_self_;
  history_self_.incr_num_recombine();
  history_self_.set_how_created(CREATE_RECOMBINE);  

  if (assign_age_type_ == Age_Take_Older) {
    history_self_.take_older(parent2->history_self_);
  } else if (assign_age_type_ == Age_Take_Younger) {
    history_self_.take_younger(parent2->history_self_);
  } else if (assign_age_type_ == Age_Take_Average) {
    history_self_.take_average(parent2->history_self_);
  } else {
    cout << "individual :: error, invalid take age type: "
	 << assign_age_type_ << "\n";
    while (1) ;
  }

  return true;
}

bool Individual :: recombine_rand2(Individual *parent2)
{
  return recombine(parent2);
}

bool Individual :: recombine(vector<Individual*>& ind2s)
{
  if (ind2s.size() == 0) {
    return false;
  }

  return recombine(ind2s[0]);
}



/*************************************************/

//int Individual :: phenotype_distance(Individual* ind2)
int Individual :: phenotype_distance(Individual*)
{
  // Should be implented in subclass (if appropriate).
  return 0;
}


//bool Individual :: make_phenotype(void* arg)
bool Individual :: make_phenotype(void*)
{
  // Should be implented in subclass (if appropriate).
  return true;
}

bool Individual :: make_phenotype()
{
  return true;
}


/*************************************************/




///////////////////////////////////////////////////////////


void Individual :: write_log(char *string)
{
  int loc = 0;

  for (int i=0; i<gNum_Fitness; i++) {
    sprintf(&string[loc], "%.2f %.2f ", history_self_.get_fitness(i),
	    fitness_hist_diff(i));
    while (string[loc] != 0) {
      ++loc;
    }
  }

  if (history_self_.get_how_created() == CREATE_INIT) {
    sprintf(&string[loc], "3\n");
  } else if (history_self_.get_how_created() == CREATE_MUTATE) {
    sprintf(&string[loc], "1\n");
  } else if (history_self_.get_how_created() == CREATE_RECOMBINE) {
    sprintf(&string[loc], "2\n");
  } else if (history_self_.get_how_created() == CREATE_RANDOM) {
    sprintf(&string[loc], "4\n");
  } else {
    sprintf(&string[loc], "-1\n");
  }
}


ostream& Individual :: write_log(ostream& log_stream)
{
  //  printf("write log\n");

  for (int i=0; i<gNum_Fitness; i++) {
    log_stream << history_self_.get_fitness(i) << " ";
  }

  log_stream << get_age() << " ";

  for (int i=0; i<HIST_NUM_COMPLEXITY; i++) {
    log_stream << history_self_.get_complexity(i) << " ";
  }

  if (history_self_.get_how_created() == CREATE_INIT) {
    log_stream << "3";
  } else if (history_self_.get_how_created() == CREATE_MUTATE) {
    log_stream << "1";
  } else if (history_self_.get_how_created() == CREATE_RECOMBINE) {
    log_stream << "2";
  } else if (history_self_.get_how_created() == CREATE_RANDOM) {
    log_stream << "4";
  } else {
    log_stream << "-1";
  }

  if ((history_self_.get_how_created() == CREATE_MUTATE) ||
      (history_self_.get_how_created() == CREATE_RECOMBINE)) {
    for (int i=0; i<gNum_Fitness; i++) {
      log_stream << " " << history_parent1_.get_fitness(i);
    }

    for (int i=0; i<HIST_NUM_COMPLEXITY; i++) {
      log_stream << " " << history_parent1_.get_complexity(i);
    }

  } else {
    for (int i=0; i<gNum_Fitness; i++) {
      log_stream << " 0";
    }

    for (int i=0; i<HIST_NUM_COMPLEXITY; i++) {
      log_stream << " 0";
    }
  }

  log_stream << "\n";

  return log_stream;
}



istream& Individual :: read(istream& istr)
{
  int version = 0;
  string s;
  istr >> s >> version;
  if ((s != "IndV") || (version != INDIVIDUAL_VERSION)) {
    throw AlpsError("individual :: read() - error, bad version.");
  }

  istr >> history_self_;
  istr >> history_parent1_;
  istr >> history_parent2_;
  
  return istr;
}


bool Individual :: read(const char *fname)
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


ostream& Individual :: write(ostream& ostr)
{
  ostr << "IndV " << INDIVIDUAL_VERSION << endl;
  history_self_.write(ostr);
  ostr << history_parent1_;
  ostr << history_parent2_;
  return ostr;
}


bool Individual :: write(const char *fname)
{
  ofstream file(fname);
  write(fname);
  file.close();
  return true;
}



/*************************************************/

void Individual :: print_history_ops()
{
  history_self_.print_ops(cout);
}


void Individual :: print_results()
{
  if (history_self_.get_how_created() == CREATE_MUTATE) {
    if (history_self_.get_fitness() > history_parent1_.get_fitness())
      printf("  Higher!");
    else if (history_self_.get_fitness() == history_parent1_.get_fitness())
      printf("  Same.");
    /*
    else
      printf("  lower. %f < %f", history_self_.get_fitness(),
	     history_parent1_.get_fitness());
    */
 } else if (history_self_.get_how_created() == CREATE_RECOMBINE) {
    if ((history_self_.get_fitness() > history_parent1_.get_fitness()) &&
	(history_self_.get_fitness() > history_parent2_.get_fitness()))
      printf("  Highest!");
    else if ((history_self_.get_fitness() == history_parent1_.get_fitness()) ||
	     (history_self_.get_fitness() == history_parent2_.get_fitness()))
      printf("  Same.");
    else if ((history_self_.get_fitness() > history_parent1_.get_fitness()) ||
	     (history_self_.get_fitness() > history_parent2_.get_fitness()))
      printf("  Middle.");
    /*
    else
      printf("  lower. %f < %f & %f", history_self_.get_fitness(),
	     history_parent1_.get_fitness(), history_parent2_.get_fitness());
    */
  }
}


void Individual :: print_history_full()
{
  streamsize prec = cout.precision();

  if (evaluations_ > 1) {
    double diff, val;

    cout << "Evl(" << evaluations_ << ")                     ";
    val = get_fitness();
    diff = val - history_self_.get_fitness_prev();
    history_parent1_.print_f1(cout);
    cout << " -> ";

    cout << setprecision(6);
    cout << val + diff;
    cout << " diff: " << 2.0*diff;
    //    printf("%6.0f", val + diff);
    //    printf(" diff: %6.0f ", 2.0f*diff);
    return;

  } else if (history_self_.get_how_created() == CREATE_INIT) {
    cout << "Init ";
    //    printf("                     ");
    history_self_.print(cout);
    return;

  } else if (history_self_.get_how_created() == CREATE_RANDOM) {
    cout << "Rand ";
    //    printf("                     ");
    history_self_.print(cout);
    return;

  } else if (history_self_.get_how_created() == CREATE_MUTATE) {
    cout << "Mut ";
    cout << "                         ";
    history_parent1_.print_f1(cout);

  } else if (history_self_.get_how_created() == CREATE_RECOMBINE) {
    cout << "Rec ";
    history_parent1_.print_f1(cout);
    cout << " + ";
    history_parent2_.print_f1(cout);

  }

  cout << " -> ";
  history_self_.print_f1(cout);

  print_results();

  cout.precision(prec);
}


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

// Returns:
//  true: if ind1 is not as good as ind2
//  false: otherwise
bool individ_compare_max(void const *elt1, void const *elt2)
{
  Individual* ind1 = (Individual*)elt1;
  Individual* ind2 = (Individual*)elt2;

  if (ind1->compare_fitness(true, ind2) != 1) {
    return true;
  }
  return false;
}

bool individ_compare_min(const void *elt1, const void *elt2)
{
  Individual* ind1 = (Individual*)elt1;
  Individual* ind2 = (Individual*)elt2;

  if (ind1->compare_fitness(false, ind2) != 1) {
    return true;
  }
  return false;
}


/***********************************************************/


void individual_print_op_history()
{
  if (gNum_Operators == 0)
    return;

  cout << "m:" << gNum_Mutations_Better << ":" << gNum_Mutations_Equal
       << ":" << gNum_Mutations
       << ", r:" << gNum_Recombinations_Better << ":"
       << gNum_Recombinations_Equal << ":" << gNum_Recombinations;
}


void individual_print_op_history(char *string)
{
  int i, loc;

  if (gNum_Operators == 0)
    return;

  sprintf(string, "[");
  loc = 1;
  for (i=0; i<gNum_Operators-1; i++) {
    sprintf(&string[loc], "%2d ]", gOp_Used_Hist[i]);
    while (string[loc] != ']')
      loc++;
  }
  sprintf(&string[loc], "%2d] ", gOp_Used_Hist[i]);

  while (string[loc] != ' ')
    loc++;

  sprintf(&string[loc], "  m:%d:%d:%d r:%d:%d:%d\n",
	  gNum_Mutations_Better, gNum_Mutations_Equal, gNum_Mutations,
	  gNum_Recombinations_Better, gNum_Recombinations_Equal, gNum_Recombinations);
}


void individual_write_op_history(char *string)
{
  int i, loc;

  sprintf(string, "%d\n", gNum_Operators);
  if (gNum_Operators == 0)
    return;

  loc = 1;
  for (i=0; i<gNum_Operators; i++) {
    sprintf(&string[loc], " %d\n", gOp_Used_Hist[i]);
    while (string[loc] != '\n')
      loc++;
  }

  sprintf(&string[loc], "  m:%d:%d:%d r:%d:%d:%d\n",
	  gNum_Mutations_Better, gNum_Mutations_Equal,
	  gNum_Mutations, gNum_Recombinations_Better,
	  gNum_Recombinations_Equal, gNum_Recombinations);
}


void individual_read_op_history(const char *string)
{
  int i, loc, num;


  num = read_integer(string);
  if (num != gNum_Operators) {
    printf("individual_read_op_history :: # operators does not match: %d != %d.\n",
	   num, gNum_Operators);
  }

  loc = 0;
  for (i=0; i<num; i++) {
    while (string[loc] != ' ') {
      loc++;
      if (loc > 1000) {
	printf("individual_read_op_history :: error reading op hist %d:%s\n",
	       loc, string);
	return;
      }
    }
    loc++;
    gOp_Used_Hist[i] = read_integer(&string[loc]);
  }

  while (string[loc] != 'm') {
    loc++;
    if (loc > 1000) {
      printf("individual_read_op_history :: error finding varation history %d:%s\n",
	     loc, string);
      return;
    }
  }

  num = sscanf(&string[loc], "m:%d:%d:%d r:%d:%d:%d\n",
	       &gNum_Mutations_Better, &gNum_Mutations_Equal,
	       &gNum_Mutations, &gNum_Recombinations_Better,
	       &gNum_Recombinations_Equal, &gNum_Recombinations);
  if (num != 6) {
    printf("individual_read_op_history :: error read %d:  m:%d:%d:%d r:%d:%d:%d\n",
	   num, gNum_Mutations_Better, gNum_Mutations_Equal,
	   gNum_Mutations, gNum_Recombinations_Better,
	   gNum_Recombinations_Equal, gNum_Recombinations);
  }
}


}
