/************************************************************

  File: alps_utils.cpp
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

************************************************************/

#include "alps_utils.h"
#include "alps_random_mt.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <sys/time.h>


using namespace alps_random;
using namespace std;





namespace alps_utils {

const int Max_Line_Length = 255; // max line allowed in configuration file.


inline double reduce(double val)
{
  return (floor((val) * 1.0E+8 + 0.5) / 1.0E+8);
}

inline double reduce2(double val)
{
  return (floor((val) * 1.0E+5 + 0.5) / 1.0E+5);
}


int Option[50];

int Max_operator;
int Last_op;
int Last_op2;
int Op_better[100];
int Operation_hist[100];
int Op_better_hist[100];

const int OP_STACK_SIZE = 20;
int Op_stack[OP_STACK_SIZE];
int Op_stack_ptr = 0;
int Num_fit_dims = 1;


double Rec_delta = 1.0;


template<typename T> inline T absv(const T a)
{
  return ((a >= 0.0) ? a : -a);
}



void print_vec2(char *string, double vec[2])
{
  printf("%s %f %f\n", string, vec[0], vec[1]);
}


double recombine_double(double val1, double val2)
{
  double diff;

  diff = Rec_delta * (val1 - val2);
  return reduce(val1 + diff * (random_double() * 2.0 - 1.0));
}


double recombine(double val1, double val2)
{
  if (val1 == 0.0) {
    if (random_int(2))
      return val1;
    else
      return val2;
  } else if (val2 == 0.0) {
    if (random_int(2))
      return val1;
    else
      return val2;
  }

  double diff;

  diff = Rec_delta * (val1 - val2);
  return reduce(val1 + diff * (random_double() * 2.0f - 1.0f));
}



void recombine_elliptic(vector<double>& parent1, vector<double>& parent2,
			vector<double>& child1)
{
  cerr << "alps_util :: recombine_elliptic() - not implemented.\n";
  while (1) ;
  
   /* child1 - the random vector in the circle */
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     child1[i] = random_norm();
   }


   /* normalize child1 to be a unit vector */
   double val = 0.0;
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     val += child1[i] * child1[i];
   }
   val = sqrt(val);

   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     child1[i] /= val;
   }

   /* give it a random length */
   val = pow(random_double(), 1.0/(double)child1.size());
   for (vector<double>::size_type i=0; i<child1.size(); i++)
     child1[i] *= val;

   cout << "C1: ";
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     cout << child1[i] << " ";
   }  
   cout << "\n";

   
   /* Set up I, a and b, and d */
   
   /* d - the vector from p1 to p2 */
   //   Individual d = pop_string_temp(string_pop(parent1));
   //   Individual uu = pop_string_temp2(string_pop(parent1));

   vector<double> d(parent2);
   for (vector<double>::size_type i=0; i<parent1.size(); i++) {
     d[i] -= parent1[i];
   }

   val = 0.0;
   for (vector<double>::size_type i=0; i<parent1.size(); i++) {
     val += d[i] * d[i];
   }
   double I = sqrt(val);

   //   double a = I * 1.0 / sqrt((double)child1.size());
   double a = I * 1.0;
   double b = I * 0.3 / sqrt((double)child1.size());

   /* uu is projection of child1/u onto d */
   double sum = 0.0;
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     sum += child1[i] * d[i];
   }
	 
   val = sum / val;
   vector<double> uu(d);
   for (vector<double>::size_type i=0; i<parent1.size(); i++) {
     uu[i] *= val;
   }


   /* Now stretch child1 to proper shape and translate it */
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     child1[i] = b * (child1[i] - uu[i])
       + a * (uu[i] + parent1[i]);
   }
   
   
   cout << "P1: ";
   for (vector<double>::size_type i=0; i<parent1.size(); i++) {
     cout << parent1[i] << " ";
   }  
   cout << "\n";

   cout << "P2: ";
   for (vector<double>::size_type i=0; i<parent2.size(); i++) {
     cout << parent2[i] << " ";
   }  
   cout << "\n";



   cout << "d:  ";
   for (vector<double>::size_type i=0; i<d.size(); i++) {
     cout << d[i] << " ";
   }  
   cout << "\n";

   cout << "UU: ";
   for (vector<double>::size_type i=0; i<uu.size(); i++) {
     cout << setw(7) << uu[i] << " ";
   }  
   cout << "\n";



   cout << "C1: ";
   for (vector<double>::size_type i=0; i<child1.size(); i++) {
     cout << child1[i] << " ";
   }  
   cout << "\n";
}




const double INTEGER_RANGE = 100.0;
const int INT_RANGE = 100;


int recombine_int(int val1, int val2)
{
  
  if (val1 == val2)
    return val1;

  if (val1 < val2) {
    int diff = val2 - val1;
    diff = random_int(2 * diff + 1);
//    printf("general :: recombine_int() - a. %d:%d -> diff=%d -> %d => ",
//    	   val1, val2, (val2-val1), diff);

    val2 -= diff;
    //    printf("%d.\n", val2);
    return val2;
  }

  int diff = val1 - val2;
  diff = random_int(2 * diff + 1);
  //  printf("general :: recombine_int() - b. %d:%d -> diff=%d -> %d => ",
  //	 val1, val2, (val1-val2), diff);
  
  val2 += diff;
  //  printf("%d,\n", val1);
  return val2;
}

double recombine_double_as_int(double val1, double val2)
{
  
  if (val1 == val2)
    return val1;

  int diff = int(val1 - val2);
  diff = random_int(2 * diff + 1) - diff;
  //  diff = rand()%(2 * diff + 1) - diff;

  val1 += double(diff);
  if (val1 > INTEGER_RANGE)
    val1 = INTEGER_RANGE;
  else if (val1 < 0.0)
    val1 = 0.0;
  return val1;
}

bool Op_used(int op_num)
{
  if (Op_stack_ptr >= OP_STACK_SIZE)
    return false;

  Op_stack[Op_stack_ptr++] = op_num;
  return true;
}

int Op_last(void)
{
  if (Op_stack_ptr > 0)
    return Op_stack[Op_stack_ptr-1];
  return -1;
}


void Op_reset(void)
{
  Op_stack_ptr = 0;
}

void Op_update(vector<int>& operation)
{
  int i;

  //  cout << "alps_util :: Op_update() - #" << Op_stack_ptr << endl;
  for (i=0; i<Op_stack_ptr; i++) {
    //    cout << "alps_util :: op_update(): " << Op_stack[i] << endl;
    if (Op_stack[i] >= (int)operation.size()) {
      operation.resize(Op_stack[i]+1);
    }
    operation[Op_stack[i]]++;
  }
  /*
  if (better) {
    for (i=0; i<Op_stack_ptr; i++)
      Op_better[Op_stack[i]]++;
  }
  */
  Op_stack_ptr = 0;
}


int next_field(const char *line)
{
  int i = 0;

  while ((line[i] != ' ') && (line[i] != ':') &&
	 (line[i] != ',') && (line[i] != 0)) {
    i++;
  }

  if (line[i] != 0)
    return i+1;
  return i;
}


int read_integer(const char *line)
{
  int i = 0;
  int val;

  char tmp[255];
  while ((line[i] != ' ') && (line[i] != ',')
	 && (line[i] != ')') && (line[i] != '>') &&
	 (line[i] != 0)) {
    tmp[i] = line[i];
    if (++i>254) {
      // Should throw an error or something.
      cout << "alps_util :: read_integer() error reading line: " << line << endl;
      while (1) ;
    }
  }

  tmp[i] = 0;
  val = atoi(tmp);

  return val;
}


double read_double(const char *line)
{
  double val;

  if ((line[0] != '-') &&
      ((line[0] < '0') || (line[0] > '9'))) {
    // not a number.
    printf("read_double() - error reading value at: %s.\n", line);
    return -9999.9;
  }

  int i = 0;
  char tmp[255];
  while ((line[i] != ' ') && (line[i] != ':') &&
	 (line[i] != 0)) {
    i++;
    tmp[i] = line[i];
    if (i>254) {
      // Should throw an error or something.
      cout << "alps_util :: read_double() error reading line: " << line << endl;
      while (1) ;
    }
  }

  tmp[i] = 0;
  val = atof(tmp);

  return val;
}


bool get_int(const char *string, int *val)
{
  int i;

  i = find_field(string);
  *val = atoi(string+i);
  return true;
}


bool get_double(const char *string, double *val)
{
  int i;

  i = find_field(string);
  *val = atof(string+i);
  return true;
}



// Returns location of input field, 0 if not found.
int find_field(const char *string)
{
  int i=0;
  // Find separator: ':'.
  while (string[i] != ':') {
    i++;
    if (i > Max_Line_Length) {
      printf("Pop :: parsing configuration file - error finding separator<%s>.\n", string);
      return 0;
    }
  }
  i++;

  // Move to end of whitespace.
  while (string[i] == ' ') {
    i++;
    if (i >= Max_Line_Length) {
      printf("Pop :: parsing configuration file - error finding field.\n");
      return 0;
    }
  }

  return i;
}



void print_time()
{
  int sec, min, hour;
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);

  tp.tv_sec -= 7*3600; // change to PST.
  hour = (tp.tv_sec/3600) % 24;
  min = (tp.tv_sec/60) % 60;
  sec = tp.tv_sec % 60;

  cout << "(" << hour << ":" << min << ":" << sec << ")";
}



AlpsError :: AlpsError(const char* msg)
  : data(msg)
{
}



}
