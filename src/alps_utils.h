/************************************************************
  Copyright Gregory S. Hornby.
  All Rights Reserved.
  File: alps_general.h
************************************************************/

#ifndef ALPS_UTIL_HEADER_FILE
#define ALPS_UTIL_HEADER_FILE

#include <cmath>
#include <vector>


namespace alps_utils {


extern const int Max_Line_Length;
extern int Max_operator;

class AlpsError {
  const char* const data;
public:
  AlpsError(const char* msg = 0);
};


bool Op_used(int op_num);
int Op_last(void);
void Op_reset(void);
void Op_update(std::vector<int>& operation);



int recombine_int(int val1, int val2);
double recombine_double_as_int(double val1, double val2);
double recombine_double(double val1, double val2);
double recombine(double val1, double val2);
void recombine_elliptic(std::vector<double>& p1, std::vector<double>& p2,
			std::vector<double>& c1);

int next_field(const char *line);
int read_integer(const char *line);
double read_double(const char *line);

bool get_int(const char *string, int *val);
bool get_double(const char *string, double *val);
int find_field(const char *string);

void print_time(void);

}

#endif
