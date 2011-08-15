/*******************************************************

  File: evoreal.cpp
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

#include "alps/alps_random_mt.h"
using namespace alps_random;


#include "alps/alps_utils.h"
using namespace alps_utils;


#include "alps/alps.h"
#include "alps/alps_gen.h"
#include "alps/alps_sstate.h"
#include "alps/alps_layer.h"
#include "alps/alps_individ_real.h"
using namespace alps;


#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#include <sys/time.h>


/***********************************/

int Fitness_Func = 0;
int Num_Evaluations = 0;

bool Save_Log = false;
ofstream Log_File;


/***********************************/


bool Finished = false;

char Seed_Individual[255];
char Population_Name[255];
bool gSet_Population_Name = false;

string Logfile_Name;
bool Set_Logfile_Name = false;

int Population_Size = -1;
int Number_Layers = -1;
int Number_Generations = -1;
int Num_Runs = -1;
int Population_Elitism = -1;
int Run_Number = -1;
int Number_Genes = -1;
double Population_Pressure = -1.0;
double Population_Prob_Recomb = -1.0;

Alps *Population = 0;

bool gReplay = false;
char *gFilename1 = 0;

int Print[50];

bool configure_pop(Alps* pop);



int get_time_sec(void)
{
  struct timeval tp;
  struct timezone tzp;

  gettimeofday(&tp, &tzp);
  return tp.tv_sec;
}


int str_eq(const char *string1, const char *string2)
{
  int i = 0;

  while (string1[i] != 0) {
    if (string1[i] != string2[i])
      return false;
    i++;
    if (i > 500) {
      cout << "Pop :: parsing configuration file - max line length exceeded.\n";
      return false;
    }
  }
  return true;
}



double f102(double x, double y)
{
  // F102 or Rana Function.
  // f102(x,y) = xsin(sqrt(|y+1-x|))cos(sqrt|x+y+1|)
  //             + (y+1)cos(sqrt|y+1-x|)sin(sqrt|x+y+1|)
  // Note this implementation is NOT shifted or rotated.
  // An improved test function would do both.

  double sum1 = y+1.0-x;
  sum1 = sqrt(abs(sum1));

  double sum2 = x+y+1.0;
  sum2 = sqrt(abs(sum2));

  double fitness = x*sin(sum1)*cos(sum2);
  fitness += (y+1.0)*cos(sum1)*sin(sum2);
  //  fitness += 511.707735; // To put the minimum at 0.
  return fitness;
}

double f102_func(const vector<double> params)
{
  double num = 1.0;
  double fitness = f102(params.back(), params[0]);

  for (vector<double>::size_type i=1; i<params.size(); i++) {
    num += 1.0;
    double val = f102(params[i-1], params[i]);
    fitness += val;
  }

  fitness /= num;

  return fitness;  
}



bool evaluate_indreal(vector<double>& fitness, const vector<double>& genes)
{
  fitness[0] = f102_func(genes);

  return true;
}

bool evaluate_individ(vector<double>& fitness, Individual* individ)
{
  if (individ == NULL) {
    cout << "evo_real :: evaluate_individ() - error, null individ sent.\n";
    while (1) ;
    return false;
  }

  bool res = evaluate_indreal(fitness, ((Individ_Real*)individ)->get_genes());

  return res;
}


int reevaluate_individ(const char *filename)
{
  Individ_Real *individ1 = new Individ_Real;
  if (!individ1->read(filename)) {
    printf("evoreal :: error reading %s\n", filename);
    return false;
  }

  if (!individ1->read(filename)) {
    printf("evo_real :: error reading %s\n", filename);
    return false;
  }


  vector<double> fitness;
  fitness.resize(1);
  int result = evaluate_individ(fitness, individ1);

  cout << "Re-eval fit: " << fitness[0] << " : ";
  individ1->print_history();
  individ1->print_history_ops();
  cout << endl;

  individ1->write("read.ind");

  return result;
}



bool configure_pop(Alps* pop)
{
  if (Num_Runs > 0) {
    pop->set_num_runs(Num_Runs);
  }

  if (Run_Number > -1) {
    pop->set_run_number(Run_Number);
  }

  if (gSet_Population_Name) {
    pop->set_name(Population_Name);
  }


  if (Set_Logfile_Name) {
    pop->set_logfile(Logfile_Name.c_str());
  }

  return true;
}



void setup_pop_gen(Individual* individ_config, AlpsGen* pop)
{
  configure_pop(pop);
  pop->set_save_best(true);


  AlpsLayer layer_def;
  int type = 2;
  int age_gap = 5;
  int age_scheme = ALPS_AGING_FIBONACCI1;
  if (type == 1) {
    // Configuration for a regular EA/GA:
    Number_Layers = 1;
    layer_def.set_select_type(ALPS_SELECT_TOURN);
    layer_def.set_size(400);
    layer_def.set_elitism(5);
    layer_def.set_tourn_size(5);
    pop->set_recomb_prob(0.5);
    pop->set_rec_rand2_prob(1.0); // 1.0
    pop->set_print_results_rate(400); // 400

  } else if (type == 2) {
    Number_Layers = 10;
    age_gap = 3; //4
    age_scheme = ALPS_AGING_EXP3;
    Number_Layers = 10;
    layer_def.set_select_type(ALPS_SELECT_TOURN);
    //    layer_def.set_select_type(ALPS_SELECT_DC);
    layer_def.set_size(40);
    layer_def.set_elitism(5);
    layer_def.set_tourn_size(5);
    layer_def.set_prob_select_prev(0.2);
    pop->set_recomb_prob(0.5);
    pop->set_rec_rand2_prob(1.0);
    pop->set_print_results_rate(1);

  } else if (type == 3) {
    Number_Layers = 10;
    age_gap = 3; //4
    age_scheme = ALPS_AGING_EXP3;
    Number_Layers = 10;
    //    layer_def.set_select_type(ALPS_SELECT_TOURN);
    layer_def.set_select_type(ALPS_SELECT_DC);
    layer_def.set_size(40);
    layer_def.set_elitism(0);
    layer_def.set_tourn_size(5);
    layer_def.set_prob_select_prev(0.25);
    pop->set_recomb_prob(0.5);
    pop->set_rec_rand2_prob(1.0);
    pop->set_print_results_rate(400);

  } else {
    cout << "evoreal :: setup_pop_gen() - error, invalid EA type: "
	 << type << "\n";
    return;
  }

  pop->config_layers_same(age_scheme, age_gap,
			  Number_Layers, layer_def);
  pop->print_layers();
  pop->set_num_runs(1);
  pop->set_max_gen(2525);
}



void setup_pop_sstate(Individual* individ_config, AlpsSState* pop)
{
  configure_pop(pop);
  pop->set_save_best(true);

  AlpsLayer layer_def;
  int type = 2;
  int age_gap = 5;
  int age_scheme = ALPS_AGING_FIBONACCI1;
  if (type == 1) {
    // Configuration for a regular EA/GA:
    Number_Layers = 1;
    layer_def.set_select_type(ALPS_SELECT_TOURN);
    layer_def.set_size(400);
    layer_def.set_elitism(5);
    layer_def.set_tourn_size(5);
    pop->set_recomb_prob(0.5);
    pop->set_rec_rand2_prob(1.0); // 1.0
    pop->set_print_results_rate(400); // 400

  } else if (type == 2) {
    // Configuration for a 10-layer ALPS run:
    //  Note: targeted for ~1million evaluations.
    age_gap = 3; //4
    age_scheme = ALPS_AGING_EXP3;
    Number_Layers = 10;
    layer_def.set_select_type(ALPS_SELECT_TOURN);
    layer_def.set_size(40);
    layer_def.set_elitism(5);
    layer_def.set_tourn_size(5);
    layer_def.set_prob_select_prev(0.25);
    pop->set_recomb_prob(0.5);
    pop->set_rec_rand2_prob(1.0);
    pop->set_print_results_rate(400);

  } else {
    cout << "evoreal :: setup_pop_sstate() - error, invalid EA type: "
	 << type << "\n";
    return;
  }

  pop->config_layers_same(age_scheme, age_gap,
			  Number_Layers, layer_def);
  pop->print_layers();
  pop->set_max_evals(1010000);
}




void *ea_engine(void *arg1)
{
  cout << "EA engine started.\n";

  if (Number_Genes <= 0) {
    Number_Genes = 20;
  }

  Individ_Real *individ_config = new Individ_Real(Number_Genes);
    //    individ_config->set_init_minmax(-2.048, 2.047);
    //    individ_config->set_minmax(-2.048, 2.047);
  individ_config->set_init_minmax(-512.0, 511.0);
  individ_config->set_minmax(-512.0, 511.0);
  //  individ_config->set_init_minmax(-100.0, 100.0);
  //  individ_config->set_minmax(-100.0, 100.0);


  vector<double> fitness;
  fitness.resize(1);


  // Configure a generational ALPS population:
  Population = new AlpsGen("real", individ_config);
  setup_pop_gen(individ_config, (AlpsGen*)Population);

  /*
  // Configure a steady-state ALPS population:
  Population = new AlpsSState("real", individ_config);
  setup_pop_sstate(individ_config, (AlpsSState*)Population);
  */

  //  Population->set_print_debug(true);
  Population->set_minimize();

  if (Save_Log) {
    Log_File << "# Func: " << Fitness_Func << "\n";
    Population->write_header(Log_File);
  }

  while (!Population->is_finished() && !Finished) {
    int index;
    Individual* individ;
    int res = Population->get_next_individ(index, individ);
    if (res == -1) {
      continue; // Get another index / Try again.

    } else if (res == -2) {
      // Evolution is over.
      break;
    }


    vector<double> fitness;
    fitness.resize(1);
    int result = evaluate_individ(fitness, individ);
    if (result == false) {
      // Error evaluating this individual.
      Population->evaluate_error(index, individ);
    } else {
      // Evaluated successfully.
      Population->insert_evaluated(fitness, index, individ, Print[1]);
    }
  }

  printf("EA engine ended.\n");

  if (Save_Log) {
    Log_File.close();
  }

  return 0;
}


int process_args(int argc, char **argv)
{
  int i, j, index;

  for (i=1; i<argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'f':
	if (++i >= argc) {
	  cout << "Missing argument for -f." << endl;
	  return false;
	}
	
	index = atoi(argv[i]);
	if ((index < 0) || (index > 5)) {
	  cout << "Invalid function # for -f: " << index << endl;
	  return false;
	}
	Fitness_Func = index;
	printf("Using fitness function #%d.\n", Fitness_Func);
	break;

      case 'l':
	if (++i >= argc) {
	  cout << "Missing argument for -l." << endl;
	  return false;
	}

	Set_Logfile_Name = true;
	Logfile_Name.assign(argv[i]);
	break;

      case 'L':
        // Save log to this file
        if (++i >= argc) {
          cout << "Missing argument for -l." << endl;
          return false;
        }
        
        {
          string fname(argv[i]);
          fname += ".log";
          Log_File.open(fname.c_str(), ios::app);
          if ( !Log_File.is_open() ) {
            // The file could not be opened
            cout << "Error :: Log file could not be opened." << endl;
            return false;
          }
          Log_File << "# ALPS\n";
          Save_Log = true;
          cout << "Saving log to: " << fname << "\n";
        }
        break;
   
      case 'm':
	if (++i < argc) {
	  gReplay = false;
	  gFilename1 = argv[i];
	} else {
	  cout << "Missing argument for -m." << endl;
	  exit(-1);
	}
	break;

      case 'p':
	if (++i >= argc) {
	  cout << "Missing argument for -p." << endl;
	  return false;
	}
	
	index = atoi(argv[i]);
	if ((index < 0) || (index > 49)) {
	  cout << "Invalid argument for -p: " << index << endl;
	  return false;
	}
	Print[index] = false;
	break;
	
      case 'P':
	if (++i >= argc) {
	  cout << "Missing argument for -P." << endl;
	  return false;
	}
	
	if (argv[i][0] == 'A') {
	  for (j=0; j<10; j++)
	    Print[j] = true;
	} else if (argv[i][0] == 'a') {
	  for (j=0; j<10; j++)
	    Print[j] = false;
	} else {
	  index = atoi(argv[i]);
	  if ((index < 0) || (index > 50)) {
	    cout << "Missing argument for -P." << endl;
	    return false;
	  }
	  Print[index] = true;
	}
	break;

      case 'r':
	if (++i < argc) {
	  gReplay = true;
	  gFilename1 = argv[i];
	} else {
	  cout << "Missing argument for -r." << endl;
	  exit(-1);
	}

      case 's':
	if (++i >= argc) {
	  cout << "Missing argument for -s." << endl;
	  return false;
	}
	
	int seed = atoi(argv[i]);
	seed_random(seed);
	break;
      }
	
    } else {
      printf("Error in command line >%c<.\n", argv[i][1]);
      return false;
    }
  }

  return true;
}



const int Max_Line_Len = 255;

// Parses a line from the configuration file,
// returns false if a problem, otherwise true.
bool parse_config_line(char *config_line, int verbose)
{
  int i, j, i_val, loc;
  double f_val;
  
  for (i=0; i<Max_Line_Len; i++) {
    if (config_line[i] == '\n')
      return true; // Empty line, no error in parsing.
    if (config_line[i] == '#')
      break; // Comment, skip rest of line.
    if (config_line[i] != ' ')
      break;
  }

  // Population size.
  if (str_eq("pop size", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 1) {
      printf("evo_real :: parse_config_file() - Error, number of individuals (%d) must be > 0.\n",
	     i_val);
      return false;
    }

    Population_Size = i_val;
    if (verbose) {
      printf("evo_real :: parse_config_file() - setting number of individuals to: %d.\n",
	     Population_Size);
    }
    return true;
  }

  // Number of layers:
  if (str_eq("layers", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 1) {
      printf("evo_real :: parse_config_file() - Error, number of layers (%d) must be > 0.\n",
	     i_val);
      return false;
    }

    Number_Layers = i_val;
    if (verbose) {
      printf("evo_real :: parse_config_file() - setting number of layers to: %d.\n",
	     Number_Layers);
    }
    return true;
  }

  // Number of generations
  if (str_eq("num generation", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 1) {
      printf("evo_real :: parse_config_file() - Error, number of generations must be 0 < %d.\n",
	     i_val);
      return false;
    }

    Number_Generations = i_val;
    if (verbose) {
      printf("evo_real :: parse_config_file() - setting number of generations to: %d.\n",
	     Number_Generations);
    }
    return true;
  }

  // Number of runs
  if (str_eq("num runs", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 1) {
      printf("evo_real :: parse_config_file() - Error, number of runs must be 0 < %d.\n",
	     i_val);
      return false;
    }

    Num_Runs = i_val;
    if (verbose) {
      printf("evo_real :: parse_config_file() - setting number of runs to: %d.\n", Num_Runs);
    }
    return true;
  }

  // Elitism
  if (str_eq("elitism", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 0) {
      printf("evo_real :: parse_config_file() - Error, elitism must be >= %d.\n",
	     i_val);
      return false;
    }

    Population_Elitism = i_val;
    if (verbose) {
      printf("evo_real :: parse_config_file() - setting elitism to: %d.\n", 
	     Population_Elitism);
    }
    return true;
  }

  // Starting run number
  if (str_eq("starting run", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 1) {
      printf("evo_real :: parse_config_file() - Error, starting number, %d,  must be >0.\n",
	     i_val);
      return false;
    }

    Run_Number = i_val;
    if (verbose) {
      printf("evo_real :: parse_config_file() - setting run number to: %d.\n", Run_Number);
    }
    return true;
  }

  // Selective Pressure.
  if (str_eq("pressure", &config_line[i])) {
    if (!get_double(&config_line[i], &f_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if ((f_val <= 0.0f) || (f_val >= 1.0f)) {
      printf("evo_real :: parse_config_file() - Error, pressure must be 0.0 < %.2f < 1.0.\n",
	     f_val);
      return false;
    }

    Population_Pressure = f_val;
    if (verbose) {
      printf("evo_real :: parse_config_file() - setting pressure to: %f.\n",
	     Population_Pressure);
    }
    return true;
  }


  // Recombine rate.
  if (str_eq("recomb rate", &config_line[i])) {
    if (!get_double(&config_line[i], &f_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if ((f_val < 0.0f) || (f_val > 1.0f)) {
      printf("evo_real :: parse_config_file() - Error, recombination rate must be 0.0 < %.2f < 1.0.\n",
	     f_val);
      return false;
    }

    Population_Prob_Recomb = f_val;
    if (verbose) {
      printf("evo_real :: parse_config_file() - setting recombination rate to: %f.\n",
	     Population_Prob_Recomb);
    }
    return true;
  }


  // Random number seed.
  if (str_eq("random seed", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 0) {
      printf("evo_real :: parse_config_file() - Error, random seed, %d < 0.\n",
	     i_val);
      return false;
    }

    if (verbose) {
      printf("evo_real :: parse_config_file() - setting random seed %d.\n", i_val);
    }
    seed_random(i_val);
    srand(i_val);
    while (i_val-- > 0) {
      f_val = random_double();
      i = random_int(100);
    }
    return true;
  }

  // Log file name.
  if (str_eq("log file", &config_line[i])) {
    loc = find_field(&config_line[i]);
    j = loc;
    while (config_line[j] != '\n') {
      j++;
      if (j >= Max_Line_Len) {
	config_line[j-1] = 0;
	printf("evo_real :: parse_config_file() - parsing - run name, max line length exceeded: %s.\n",
	       config_line);
	return false;
      }      
    }

    if (config_line[j] == '\n') {
      config_line[j] = 0;
    }

    loc += i;
    Logfile_Name.assign((char*)config_line[loc]);
    Set_Logfile_Name = true;
    cout << "Saving log file: " << Logfile_Name << "\n";

    return true;
  }


  ////////// Individual settings.

  // Number of genes.
  if (str_eq("num genes", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if (i_val < 2) {
      printf("evo_real :: parse_config_file() - Error, number of genes must be %d > 1.\n",
	     i_val);
      return false;
    }

    Number_Genes = i_val;
    if (verbose) {
      printf("evo_real :: parse_config_file() - setting number of genes to: %d.\n",
	     Number_Genes);
    }
    return true;
  }


  ////////////////////////////
  // Print commands.
  if (str_eq("Print", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if ((i_val < 0) || (i_val >= 50)) {
      printf("evo_real :: parse_config_file() - Error, 0<= Print[%d] < 50.\n", i_val);
      return false;
    }
    
    Print[i_val] = 1;
    return true;
  }
  if (str_eq("print", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("evo_real :: parse_config_file() - Error parsing: %s\n", config_line);
      return false;
    }

    if ((i_val < 0) || (i_val >= 50)) {
      printf("evo_real :: parse_config_file() - Error, 0<= Print[%d] < 50.\n", i_val);
      return false;
    }
    
    Print[i_val] = 0;
    return true;
  }

  printf("evo_real :: parse_config_file() - Error with: %s", config_line);
  return false; // Not a recognized line type.
}


bool parse_configure_file(const char *fname, int verbose)
{
  ifstream ifile(fname);

  /*
  while (get a line from file)
    if (line[0] == '#') {
      // Ignore comment lines.
      continue;
    }

    parse_config_line(line, verbose);
  }
  */
  
  ifile.close();

  return true;
}


int main(int argc, char **argv) {
  cout << "ALPS example file.  Copyright (C) [2009] [G. S. Hornby].\n";
  cout << "This program is free software: you can redistribute it and/or modify\n";
  cout << "it under the terms of the GNU General Public License as published by\n";
  cout << "the Free Software Foundation, either version 3 of the License, or\n";
  cout << "(at your option) any later version.\n";
  cout << "\n";
  cout << "This program is distributed in the hope that it will be useful,\n";
  cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
  cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
  cout << "GNU General Public License for more details.\n";
  cout << "\n";
  cout << "You should have received a copy of the GNU General Public License\n";
  cout << "along with this program.  If not, see <http://www.gnu.org/licenses/>.\n";
  cout << "\n";


  seed_random(0); // Input 0 to seed from clock, otherwise uses the input number.

  parse_configure_file(".evo_real", true);
  if (!process_args(argc, argv)) {
    return 0;
  }

  ea_engine(0);

  return 0;
}


