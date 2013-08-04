/**********************************************
  Copyright Gregory S. Hornby.
  All Rights Reserved.
  File: main.cpp
**********************************************/


#include <sys/types.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <ostream>
#include <iostream>
#include <string>
#include <malloc.h>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include <math.h>


#include <mcheck.h> /* Header file to include mtrace related functions */



#include "general.h"

#include "alps.h"
#include "individual.h"
#include "individ_bits.h"

using namespace std;



char Seed_Individual[255];
char Population_Name[255];
char Population_Restart[255];
char Logfile_Name[255];
int Set_Logfile_Name = FALSE;
int Max_Repeats = 4;
int Do_Population_Restart = FALSE;
int Population_Size = -1;
int Number_Layers = -1;
int Number_Generations = -1;
int Num_Runs = -1;
int Population_Elitism = -1;
int Run_Number = -1;
int Number_Genes = -1;
int Number_Productions = -1;
int Number_Bodies = -1;
int Keep_Stats = TRUE;
Real Population_Pressure = -1.0;
Real Population_Prob_Recomb = -1.0;

int Changed_Prob_HTS = FALSE;
int Changed_Prob_Inversion = FALSE;
int Changed_Prob_Invert_HTS = FALSE;
int* Bit_Order = 0;
int Shuffled_Hiff = FALSE;
int Randomized_Hiff = FALSE;
int* Randomized_Value = 0;


Individ_bits *Individ_Config = 0;
Individ_bits *individ1 = 0;
Alps *Population = 0;

int gSeed_Individual = FALSE;
int gReplay = FALSE;
int gMutate = FALSE;
int gSet_Population_Name = FALSE;
char *gFilename1 = 0;
char *gFilename2 = 0;

int (*Translator_func)(char*);



int configure_pop(Alps* pop);



int get_time_sec(void)
{
  struct timeval tp;
  struct timezone tzp;

  gettimeofday(&tp, &tzp);
  return tp.tv_sec;
}




int test_pop_alps(void)
{
  Individ_Config = new Individ_bits();
  Individ_Config->set_num_bits(16);
  Individ_Config->create();
  
  Population = new Alps(Individ_Config);

  if (!Population->created()) {
    Population->create();
  }

  Layer_Def layer_def;
  layer_def.select_type = POP_SELECT_TOURN;
  layer_def.size = 100;
  layer_def.prob_select_prev = 0.25;
  Population->config_layers_tree(POP_AGING_EXP, 10, 8, 2, 2, &layer_def);
  //  Population->print_layers();

  int gen = 80;
  int layer1 = Population->calc_max_layer(gen);
  Population->Generation = gen;
  int layer2 = Population->calc_max_layer();

  printf("Generation: %2d,  max layer: %2d  %2d.\n",
	 gen, layer1, layer2);

  
  

  return TRUE;
}




#define NUL -1
#define ZERO 0
#define ONE 1

#define hiff(a, b) (a == NUL) ? NUL : ((a == b) ? a : NUL)


void print_bit(char bit)
{
  if (bit == NUL)
    printf("-");
  else
    printf("%d", (int)bit);
}

void print_ind(int size, char* values)
{
  int i;
  for (i=0; i<size; i++) {
    print_bit(values[i]);
  }
}


void print_ind(Individ_bits* indbits)
{
  int len = indbits->get_num_bits();
  char values[len];
  indbits->copy_genes(values);
  print_ind(len, values);
}
  


int hiff_function(int len, char *values)
{
  int i, blocks, block_size;
  int lvl, ind[len];
  int fitness = 0, blk_val;

  block_size = 2;
  blocks = len / block_size;
  len = blocks * block_size;

  for (i=0; i<len; i++) {
    ind[i] = values[i];
  }

  lvl = 0;
  blk_val = 1;
  while (len > 0) {
    //    printf("%2d : ", len);
    for (i=0; i<len; i++) {
      if (ind[i] == ONE) {
	fitness += blk_val;
      } else if (ind[i] == ZERO) {
	fitness += blk_val;
      }
      //      print_bit(ind[i]);
    }
    //    printf(" : %2d\n", fitness);

    if (len > 1) {
      for (i=0; i<len; i+=2) {
	ind[i/2] = hiff(ind[i], ind[i+1]);
      }
    }

    lvl++;
    len = len/2;
    blk_val = blk_val * 2;
  }

  return fitness;
}



int evaluate_indbits(Real* fitness, Individ_bits* indbits)
{
  int len = indbits->get_num_bits();
  char values[len];
  char values2[len];
  indbits->copy_genes(values);

  if (Shuffled_Hiff) {
    for (int i=0; i<len; i++) {
      values2[i] = values[Bit_Order[i]];
    }

    if (Randomized_Hiff) {
      for (int i=0; i<len; i++) {
	values2[i] = values2[i] ^ Randomized_Value[i];
      }
    }

    fitness[0] = (Real)hiff_function(len, values2);
    return TRUE;
  }

  if (Randomized_Hiff) {
    for (int i=0; i<len; i++) {
      //      values2[i] = values[i];
      values[i] = values[i] ^ Randomized_Value[i];
      //      printf("%d(%d%d) ", values[i], values2[i], Randomized_Value[i]);
    }
  }
  //  printf("\n");

  fitness[0] = (Real)hiff_function(len, values);
  return TRUE;
}

int evaluate_index(Real* fitness, int index)
{
  Individ_bits *indbits;

  indbits = (Individ_bits*)Population->get_individual(index);
  if (indbits == NULL) {
    printf("evo_bits :: evaluate_index() - error getting individ: %d\n", index);
    return FALSE;
  }

  return evaluate_indbits(fitness, indbits);
}


int reevaluate_individ(char *filename)
{
  Real fitness[10];
  Individual *individ1;

  Translator_func = 0;

  individ1 = new Individual;
  if (!individ1->read(filename, Translator_func)) {
    printf("evo_bits :: error reading %s\n", filename);
    return FALSE;
  }

  Population = new Alps(individ1);
  Population->Run_Sample = TRUE;

   
  if (!Population->created()) {
    Population->create();
  }

  individ1 = Population->get_individual(0);
  if (individ1 == NULL) {
    printf("evo_bits :: ea_test_voxels2() - error getting individ: %d\n", 0);
    return FALSE;
  }

  if (!individ1->read(filename, Translator_func)) {
    printf("evo_bits :: error reading %s\n", filename);
    return FALSE;
  }

  int result = evaluate_index(fitness, 0);

  individ1->print_history();
  individ1->print_history_ops();

  individ1->write("read.ind");

  return result;
}



int mutate_individual(char *filename)
{
  int result;
  Real fitness[10];
  Individual *individ1;

  //  ea_setup_voxels();
  Translator_func = 0;

  individ1 = new Individual;
  if (!individ1->read(filename, Translator_func)) {
    printf("evo_bits :: error reading %s\n", filename);
    return FALSE;
  }

  Population = new Alps(individ1);
  Population->Run_Sample = TRUE;

  /*
  if (!Population->configure(".evolve", FALSE, Translator_func)) {
    printf("evolve :: ea_test_dist() - error configuring population.");
    return 0;
  }
  */

   
  if (!Population->created()) {
    Population->create();
  }

  individ1 = Population->get_individual(0);
  if (individ1 == NULL) {
    printf("evo_bits :: ea_test_voxels2() - error getting individ: %d\n", 0);
    return FALSE;
  }

  if (!individ1->read(filename, Translator_func)) {
    printf("evo_bits :: error reading %s\n", filename);
    return FALSE;
  }

  result = evaluate_index(fitness, 0);
  individ1->print_history();
  individ1->print_history_ops();

  result = individ1->mutate();
  if (!result) {
    printf("evo_bits :: mutate_individual() - error mutating.\n");

  } else {
    result = evaluate_index(fitness, 0);
    individ1->print_history();
    individ1->print_history_ops();
  }

  individ1->write("mutant.ind");

  return result;
}





void test_bits1(void)
{
  printf("evo_bits :: test_bits1()\n");

  Individ_Config = new Individ_bits();
  Individ_Config->set_num_bits(16);
  Individ_Config->create();
  Individ_Config->make_random();
  //  Individ_Config->write(stdout);

  Population = new Alps("bits", Individ_Config);
  configure_pop(Population);

  Layer_Def layer_def;
  //  layer_def.select_type = POP_SELECT_ROULETTE;
  layer_def.select_type = POP_SELECT_TOURN;
  layer_def.size = 10;
  layer_def.elitism = 2;
  layer_def.prob_select_prev = 0.0;
  //  Population->config_layer(&layer_def);
  //  Population->config_layer(20, POP_SELECT_TOURN);
  Population->config_layers_same(POP_AGING_EXP, 25, 3, &layer_def);


  Population->print_layers();
  if (!Population->created()) {
    Population->create();
  }


  int index = Population->get_next_index();
  printf("evo_bits :: test_bits1() - got index: %d.\n", index);

  Real fitness[10];
  int result = evaluate_index(fitness, index);
  if (result) {
    Population->set_fitness(index, fitness, Print[1], TRUE);
  } else {
    Population->evaluate_error(index);
  }


  Individ_bits* indbits = (Individ_bits*)Population->get_individual(index);
  print_ind(indbits);
  printf(" : %2d\n", (int)fitness[0]);
}

void shuffle_bit_order(void)
{
  for (int i=0; i<Number_Genes; i++) {
    Bit_Order[i] = i;
  }

  for (int i=0; i<Number_Genes; i++) {
    int loc = random_int(Number_Genes-1);
    if (loc == i) {
      continue;
    }

    int tmp = Bit_Order[i];
    Bit_Order[i] = Bit_Order[loc];
    Bit_Order[loc] = tmp;
  }

  /*
  for (int i=0; i<Number_Genes; i++) {
    printf("%2d ", Bit_Order[i]);
  }
  printf("\n");
  */
}


void randomize_bit_value(void)
{
  for (int i=0; i<Number_Genes; i++) {
    Randomized_Value[i] = random_integer()%2;
  }

  for (int i=0; i<Number_Genes; i++) {
    printf("%d", Randomized_Value[i]);
  }
  printf("\n");
}


void *ea_engine(void *arg1)
{
  int index;

  //  printf("EA engine started.\n");

  // 128   =>   1024
  // 256   =>   2304
  // 512   =>   5120
  // 1024  =>  11264
  // 2048  =>  24576
  // 4096  =>  53248
  // 8192  => 114688
  //16384  => 245760


  Individ_Config = new Individ_bits();
  //  Individ_Config->set_num_bits(16);
  //  Individ_Config->set_num_bits(512);
  //  Individ_Config->set_num_bits(1024);
  //  Individ_Config->set_num_bits(2048);
  if (Number_Genes > 0) {
    Individ_Config->set_num_bits(Number_Genes);
  } else {
    Individ_Config->set_num_bits(512);
    Number_Genes = 512;
  }


  if (Randomized_Hiff) {
    Randomized_Value = (int*)malloc(Number_Genes*sizeof(int));
    randomize_bit_value();
  }

  if (Shuffled_Hiff) {
    Individ_Config->Inversion = TRUE;
    Bit_Order = (int*)malloc(Number_Genes*sizeof(int));
    for (int i=0; i<Number_Genes; i++) {
      Bit_Order[i] = i;
    }
    if (Shuffled_Hiff) {
      shuffle_bit_order();
    }
  }


  //  Individ_Config->Inversion = TRUE;
  Individ_Config->create();
  Individ_Config->zero_bits();


  Real fitness[10];
  int tmp_rv = Randomized_Hiff;
  Randomized_Hiff = FALSE;
  evaluate_indbits(fitness, Individ_Config);
  Randomized_Hiff = tmp_rv;

  Real Max_Fitness = fitness[0];
  printf("Optimal fitness: %2d.  #bits: %d\n", (int)Max_Fitness,
	 Individ_Config->Num_Genes);


  Individ_Config->Bit_Order = Bit_Order;
  Individ_Config->Bit_Value = Randomized_Value;
  Population = new Alps("bits", Individ_Config);
  configure_pop(Population);
  Population->Save_Best = TRUE;

  if (!Changed_Prob_HTS) {
    Prob_HTS = 0.5;
  }

  if (!Changed_Prob_Inversion) {
    Prob_Inversion = 0.0;
  }

  if (!Changed_Prob_Invert_HTS) {
    Prob_Invert_HTS = 0.0;
  }


  Population->Tournament_Size = 5;

  Layer_Def layer_def;
  layer_def.max_age = -1;
  layer_def.start_index = 0;
  layer_def.parent_layer = -1;
  layer_def.num_prev_layers = 0;
  layer_def.prob_select_prev = 0.0;

  if (Population_Size == 1) {
    layer_def.select_type = POP_SELECT_DC;
    layer_def.size = 1;
    layer_def.elitism = 0;
    Population->Tournament_Size = 1;

  } else if (Population_Size > 0) {
    //    layer_def.select_type = POP_SELECT_TOURN;
    layer_def.select_type = POP_SELECT_DC;
    layer_def.size = Population_Size;
    layer_def.elitism = 2;
    Population->Tournament_Size = 5;

  } else {
    layer_def.select_type = POP_SELECT_TOURN;
    layer_def.size = 200;
    layer_def.elitism = 2;
    Population->Tournament_Size = 5;
  }

  /*
  layer_def.select_type = POP_SELECT_TOURN;
  layer_def.size = 100; // 1
  layer_def.elitism = 2; // 2
  layer_def.prob_select_prev = 0.35;
  */

  //  Population->config_layers_tree(POP_AGING_EXP, 10, 8, 2, 2, &layer_def);
  //  Population->config_layers_same(POP_AGING_EXP, 15, 20, &layer_def);
  //  Population->config_layers_same(POP_AGING_POLY1, 15, 5, &layer_def);
  //  Population->config_layers_same(POP_AGING_LINEAR1, 25, 10, &layer_def);

  // For single layer evolution:
  //  Population->config_layers_same(POP_AGING_FIBONACCI1, -1, 1, &layer_def);



  //  Population->config_layers_same(POP_AGING_EXP, 100, Number_Layers, &layer_def);
  //  Population->config_layers_same(POP_AGING_POLY1, 20, Number_Layers, &layer_def);
  Population->config_layers_same(POP_AGING_FIBONACCI1, 10, Number_Layers, &layer_def);
  //  Population->config_layers_same(POP_AGING_LINEAR1, 250, Number_Layers, &layer_def);
  //  Population->config_layer(&layer_def);
  Population->print_layers();
  Population->Print_Results_Rate = 10000; // 400

  /*
  if (!Population->configure(".evolve", TRUE, Translator_func)) {
    printf("evolve :: ea_engine() - error configuring population.");
    return 0;
  }
  */


  if (!Population->created()) {
    Population->create();
  }

  //  sem_wait(&gSem_Slave);

  //  printf("EA started\n");


  int max_evaluations = 10000000;
  int num_failures = 0;
  int num_successes = 0;
  Real total_evaluations = 0.0;

  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  int start_time_sec = tp.tv_sec;
  int start_time_usec = tp.tv_usec;
  int total_time_sec = 0;
  int total_time_usec = 0;


  while (!Population->finished()) {
    index = Population->get_next_index();
    if (index == -1) {
      continue; // Get another index / Try again.

    } else if (index == -2) {
      // Evolution is over.
      break;
    }

    Real fitness[10];
    int result = evaluate_index(fitness, index);
    if (result == FALSE) {
      // Error evaluating this individual.
      Population->evaluate_error(index);
    } else {
      Population->set_fitness(index, fitness, Print[1], TRUE);
      //    Population->set_fitness(index, &fitness, FALSE, TRUE);
    }

    if (fitness[0] == Max_Fitness) {
      //      printf("evo_bits() - success!  (%.0f)\n", fitness[0]);
      num_successes++;
      total_evaluations += Population->Num_Evaluations;
      Population->set_finished();
      printf("evo_bits() - success! (%.0f:%d)  Evals:%d:%.0f \n",
	     fitness[0], Population->get_individual(index)->get_age(),
	     Population->Num_Evaluations, total_evaluations);

      gettimeofday(&tp, &tzp);
      total_time_sec += tp.tv_sec - start_time_sec;
      total_time_usec += tp.tv_usec - start_time_usec;

      start_time_sec = tp.tv_sec;
      start_time_usec = tp.tv_usec;

      if (Shuffled_Hiff) {
	shuffle_bit_order();
      }
      if (Randomized_Hiff) {
	randomize_bit_value();
      }

    } else if (Population->Num_Evaluations >= max_evaluations) {
      //      printf("Failed after %d evaluations.\n", max_evaluations);
      num_failures++;
      Population->set_finished();

      gettimeofday(&tp, &tzp);
      total_time_sec += tp.tv_sec - start_time_sec;
      total_time_usec += tp.tv_usec - start_time_usec;

      start_time_sec = tp.tv_sec;
      start_time_usec = tp.tv_usec;
    }
  }

  printf("EA engine ended.\n");


  Real avg_time = (Real)total_time_sec + (Real)total_time_usec/1000000.0;
  avg_time /= (Real)num_successes;

  printf("Trials: %d,  success:%d,  fail:%d,  succ-rate: %.2f,  avg-evals: %.2f,  avg-time: %.2f.\n",
	 num_successes+num_failures, num_successes, num_failures,
	 (Real)num_successes/total_evaluations,
	 total_evaluations/(Real)num_successes, avg_time);

  return 0;
}



int test_individual1(void)
{
  int index;

  Individ_Config = new Individ_bits();
  Individ_Config->set_num_bits(16);
  Individ_Config->create();

  Population = new Alps(Individ_Config);

  /*
  if (!Population->configure(".evolve", TRUE, Translator_func)) {
    printf("evolve :: test_individual1() - error configuring population.");
    return 0;
  }
  */

   
  if (!Population->created()) {
    Population->create();
  }

  index = Population->get_next_index();
  if (index != 0) {
    printf("evo_bits :: test_individual1() - error, didn't get individual 0 -> %d.\n",
	   index);
    while (1) ;
  }

  individ1 = (Individ_bits*)Population->get_individual(0);
  if (individ1 == NULL) {
    printf("evo_bits :: test_individual1() - error getting individ: %d\n", 0);
    return FALSE;
  }
  

  individ1->write("individ1.ind");

  Individ_bits* individ2 = (Individ_bits*)Population->get_individual(1);
  individ2->read("individ1.ind", Translator_func);
  individ2->write("individ2.ind");

  return TRUE;
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
	  return FALSE;
	}
	
	index = atoi(argv[i]);
	if ((index < 1) || (index > 49)) {
	  cout << "Invalid function # for -f: " << index << endl;
	  return FALSE;
	}
	GParam.Fitness_func = index;
	printf("Using fitness function #%d.\n", GParam.Fitness_func);
	break;
	
      case 'm':
	if (++i < argc) {
	  gMutate = TRUE;
	  gReplay = FALSE;
	  gFilename1 = argv[i];
	} else {
	  cout << "Missing argument for -m." << endl;
	  exit(-1);
	}
	break;

      case 'p':
	if (++i >= argc) {
	  cout << "Missing argument for -p." << endl;
	  return FALSE;
	}
	
	index = atoi(argv[i]);
	if ((index < 0) || (index > 49)) {
	  cout << "Invalid argument for -p: " << index << endl;
	  return FALSE;
	}
	Print[index] = FALSE;
	break;
	
      case 'P':
	if (++i >= argc) {
	  cout << "Missing argument for -P." << endl;
	  return FALSE;
	}
	
	if (argv[i][0] == 'A') {
	  for (j=0; j<10; j++)
	    Print[j] = TRUE;
	} else if (argv[i][0] == 'a') {
	  for (j=0; j<10; j++)
	    Print[j] = FALSE;
	} else {
	  index = atoi(argv[i]);
	  if ((index < 0) || (index > 50)) {
	    cout << "Missing argument for -P." << endl;
	    return FALSE;
	  }
	  Print[index] = TRUE;
	}
	break;

      case 'o':
      case 'O':
	if (++i+1 >= argc) {
	  cout << "Missing arguments for -o." << endl;
	  return FALSE;
	}
	
	index = atoi(argv[i]);
	if ((index < 0) || (index > 49)) {
	  printf("Invalid index %d for -o.\n", index);
	  return FALSE;
	}
	Option[index] = atoi(argv[++i]);
	printf("Option[%d] = %d\n", index, Option[index]);
	break;

      case 'r':
	if (++i < argc) {
	  gReplay = TRUE;
	  gMutate = FALSE;
	  gFilename1 = argv[i];
	} else {
	  cout << "Missing argument for -r." << endl;
	  exit(-1);
	}

	if (++i < argc) {
	  if (argv[i][0] != '-') {
	    gFilename2 = argv[i];
	  } else
	    i--;
	}
	break;

      case 's':
	if (++i >= argc) {
	  cout << "Missing argument for -s." << endl;
	  return FALSE;
	}
	
	int seed = atoi(argv[i]);
	seed_random(seed);
	break;
      }
	
    } else {
      printf("Error in command line >%c<.\n", argv[i][1]);
      return FALSE;
    }
  }

  return TRUE;
}


int configure_pop(Alps* pop)
{
  if (Population_Size > 0) {
    pop->config_layer(Population_Size, POP_SELECT_TOURN, 2);
  }

  if (Number_Generations > 0) {
    pop->Max_Generations = Number_Generations;
  }

  if (Num_Runs > 0) {
    pop->Num_Runs = Num_Runs;
  }

  if (Run_Number > -1) {
    pop->Run_Number = Run_Number;
  }

  if (Population_Pressure > 0.0) {
    pop->Pressure = Population_Pressure;
  }

  if (Population_Prob_Recomb > -1.0) {
    pop->Prob_Recomb = Population_Prob_Recomb;
  }

  if (gSet_Population_Name) {
    pop->set_name(Population_Name);
  }

  if (Set_Logfile_Name) {
    pop->set_logfile(Logfile_Name);
  }

  return TRUE;
}


// Parses a line from the configuration file,
// returns FALSE if a problem, otherwise TRUE.
int parse_config_line(char *config_line, int verbose)
{
  int i, j, i_val, loc;
  Real f_val;
  
  for (i=0; i<MAX_LINE_LENGTH; i++) {
    if (config_line[i] == '\n')
      return TRUE; // Empty line, no error in parsing.
    if (config_line[i] == '#')
      break; // Comment, skip rest of line.
    if (config_line[i] != ' ')
      break;
  }

  // Population size.
  if (str_eq("pop size", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if (i_val < 1) {
      printf("evo_bits :: parse_config_file() - Error, number of individuals (%d) must be > 0.\n",
	     i_val);
      return FALSE;
    }

    Population_Size = i_val;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting number of individuals to: %d.\n",
	     Population_Size);
    }
    return TRUE;
  }

  // Number of layers:
  if (str_eq("layers", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if (i_val < 1) {
      printf("evo_bits :: parse_config_file() - Error, number of layers (%d) must be > 0.\n",
	     i_val);
      return FALSE;
    }

    Number_Layers = i_val;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting number of layers to: %d.\n",
	     Number_Layers);
    }
    return TRUE;
  }

  // Number of generations
  if (str_eq("num generation", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if (i_val < 1) {
      printf("evo_bits :: parse_config_file() - Error, number of generations must be 0 < %d.\n",
	     i_val);
      return FALSE;
    }

    Number_Generations = i_val;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting number of generations to: %d.\n",
	     Number_Generations);
    }
    return TRUE;
  }

  // Number of runs
  if (str_eq("num runs", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if (i_val < 1) {
      printf("evo_bits :: parse_config_file() - Error, number of runs must be 0 < %d.\n",
	     i_val);
      return FALSE;
    }

    Num_Runs = i_val;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting number of runs to: %d.\n", Num_Runs);
    }
    return TRUE;
  }

  // Elitism
  if (str_eq("elitism", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if (i_val < 0) {
      printf("evo_bits :: parse_config_file() - Error, elitism must be >= %d.\n",
	     i_val);
      return FALSE;
    }

    Population_Elitism = i_val;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting elitism to: %d.\n", 
	     Population_Elitism);
    }
    return TRUE;
  }

  // Starting run number
  if (str_eq("starting run", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if (i_val < 1) {
      printf("evo_bits :: parse_config_file() - Error, starting number, %d,  must be >0.\n",
	     i_val);
      return FALSE;
    }

    Run_Number = i_val;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting run number to: %d.\n", Run_Number);
    }
    return TRUE;
  }

  // Selective Pressure.
  if (str_eq("pressure", &config_line[i])) {
    if (!get_Real(&config_line[i], &f_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if ((f_val <= 0.0f) || (f_val >= 1.0f)) {
      printf("evo_bits :: parse_config_file() - Error, pressure must be 0.0 < %.2f < 1.0.\n",
	     f_val);
      return FALSE;
    }

    Population_Pressure = f_val;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting pressure to: %f.\n",
	     Population_Pressure);
    }
    return TRUE;
  }


  // Recombine rate.
  if (str_eq("recomb rate", &config_line[i])) {
    if (!get_Real(&config_line[i], &f_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if ((f_val < 0.0f) || (f_val > 1.0f)) {
      printf("evo_bits :: parse_config_file() - Error, recombination rate must be 0.0 < %.2f < 1.0.\n",
	     f_val);
      return FALSE;
    }

    Population_Prob_Recomb = f_val;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting recombination rate to: %f.\n",
	     Population_Prob_Recomb);
    }
    return TRUE;
  }


  // Random number seed.
  if (str_eq("random seed", &config_line[i])) {
    if (!get_int(&config_line[i], &i_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if (i_val < 0) {
      printf("evo_bits :: parse_config_file() - Error, random seed, %d < 0.\n",
	     i_val);
      return FALSE;
    }

    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting random seed %d.\n", i_val);
    }
    seed_random(i_val);
    srand(i_val);
    while (i_val-- > 0) {
      f_val = random_real();
      i = random_int(100);
    }
    return TRUE;
  }


  // Seed individual.
  if (str_eq("seed ind", &config_line[i])) {
    loc = find_field(&config_line[i]);
    j = i+loc;
    while (config_line[j] != '\n') {
      j++;
      if (j > MAX_LINE_LENGTH) {
	printf("evo_bits :: parse_config_file() - parsing - seed ind, max line length exceeded.\n");
	return FALSE;
      }
    }
    config_line[j] = 0;

    if (verbose) {
      printf("Seed individual:  %s\n", &(config_line[i+loc]));
    }

    int x = 0;
    do {
      Seed_Individual[x] = config_line[i+loc+x];
    } while (config_line[i+loc+x++] != 0);
    gSeed_Individual = TRUE;

    return TRUE;
  }


  // Run individual.
  if (str_eq("run", &config_line[i])) {
    loc = find_field(&config_line[i]);
    j = i+loc;
    while (config_line[j] != '\n') {
      j++;
      if (j > MAX_LINE_LENGTH) {
	printf("evo_bits :: parse_config_file() - parsing - run, max line length exceeded.\n");
	return FALSE;
      }
    }
    config_line[j] = 0;


    int x = 0;
    do {
      gFilename1[x] = config_line[i+loc+x];
    } while (config_line[i+loc+x++] != 0);
    gReplay = TRUE;

    return TRUE;
  }


  // Run name.
  if (str_eq("name", &config_line[i])) {
    loc = find_field(&config_line[i]);
    j = loc;
    while (config_line[j] != '\n') {
      j++;
      if (j >= MAX_LINE_LENGTH) {
	config_line[j-1] = 0;
	printf("evo_bits :: parse_config_file() - parsing - run name, max line length exceeded: %s.\n",
	       config_line);
	return FALSE;
      }      
    }

    if (config_line[j] == '\n') {
      config_line[j] = 0;
    }

    int x = 0;
    do {
      Population_Name[x] = config_line[i+loc+x];
    } while (config_line[i+loc+x++] != 0);
    gSet_Population_Name = TRUE;

    return TRUE;
  }


  // Log file name.
  if (str_eq("log file", &config_line[i])) {
    loc = find_field(&config_line[i]);
    j = loc;
    while (config_line[j] != '\n') {
      j++;
      if (j >= MAX_LINE_LENGTH) {
	config_line[j-1] = 0;
	printf("evo_bits :: parse_config_file() - parsing - run name, max line length exceeded: %s.\n",
	       config_line);
	return FALSE;
      }      
    }

    if (config_line[j] == '\n') {
      config_line[j] = 0;
    }

    int x = 0;
    do {
      Logfile_Name[x] = config_line[i+loc+x];
    } while (config_line[i+loc+x++] != 0);
    Set_Logfile_Name = TRUE;
    //    printf("log file: %s. (%d)\n", Logfile_Name, x);

    return TRUE;
  }


  // Restart population
  if (str_eq("restart pop", &config_line[i])) {
    loc = i+find_field(&config_line[i]);
    i = loc;
    while (config_line[i] != '\n') {
      i++;
      if (i > MAX_LINE_LENGTH) {
	printf("evo_bits :: parse_config_file() - parsing - restart pop, max line length exceeded.\n");
	return FALSE;
      }
    }
    config_line[i] = 0;

    int x = 0;
    do {
      Population_Restart[x] = config_line[loc+x];
    } while (config_line[loc+x++] != 0);
    Do_Population_Restart = TRUE;

    if (verbose) {
      printf("Restarting population: %d:%d  %s\n",
	     loc, i, &(config_line[loc]));
    }

    return TRUE;
  }

  ////////// Individual settings.

  // Number of genes.
  if (str_eq("num genes", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if (i_val < 2) {
      printf("evo_bits :: parse_config_file() - Error, number of genes must be %d > 1.\n",
	     i_val);
      return FALSE;
    }

    Number_Genes = i_val;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting number of genes to: %d.\n",
	     Number_Genes);
    }
    return TRUE;
  }


  // Probability of hierarchical tree shaped variation.
  if (str_eq("hts rate", &config_line[i])) {
    if (!get_Real(&config_line[i], &f_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if ((f_val < 0.0) || (f_val > 1.0)) {
      printf("evo_bits :: parse_config_file() - Error, HTS rate must be 0.0 < %.2f < 1.0.\n",
	     f_val);
      return FALSE;
    }

    Prob_HTS = f_val;
    Changed_Prob_HTS = TRUE;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting HTS rate to: %f.\n",
	     Prob_HTS);
    }
    return TRUE;
  }


  // Probability of individuals do inversion.
  if (str_eq("prob inversion", &config_line[i])) {
    if (!get_Real(&config_line[i], &f_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if ((f_val < 0.0) || (f_val > 1.0)) {
      printf("evo_bits :: parse_config_file() - Error, prob of doing inversion must be 0.0 < %.2f < 1.0.\n",
	     f_val);
      return FALSE;
    }

    Prob_Inversion = f_val;
    Changed_Prob_Inversion = TRUE;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting inversion probability to: %f.\n",
	     Prob_Inversion);
    }
    return TRUE;
  }


  // Probability of doing hierarchical tree shaped inversion.
  if (str_eq("prob hts inversion", &config_line[i])) {
    if (!get_Real(&config_line[i], &f_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if ((f_val < 0.0) || (f_val > 1.0)) {
      printf("evo_bits :: parse_config_file() - Error, prob of doing HTS inversion must be 0.0 < %.2f < 1.0.\n",
	     f_val);
      return FALSE;
    }

    Prob_Invert_HTS = f_val;
    Changed_Prob_Invert_HTS = TRUE;
    if (verbose) {
      printf("evo_bits :: parse_config_file() - setting HTS inversion probability to: %f.\n",
	     Prob_Invert_HTS);
    }
    return TRUE;
  }


  // Shuffled HIFF
  if (str_eq("shuffled-hiff", &config_line[i])) {
    loc = find_field(&config_line[i]);
    j = i+loc;

    if (str_eq("on", &config_line[j])) {
      Shuffled_Hiff = TRUE;
      printf("evo_bits :: parse_config_file() - shuffled hiff is on.\n");

    } else if (str_eq("off", &config_line[j])) {
      Shuffled_Hiff = FALSE;
      printf("evo_bits :: parse_config_file() - shuffled hiff is off.\n");
    } else {
      printf("evo_bits :: parse_config_file() - Error parsing: %s. (shuffled hiff)\n", config_line);
      return FALSE;
    }
    return TRUE;
  }


  // Randomized HIFF
  if (str_eq("randomized-hiff", &config_line[i])) {
    loc = find_field(&config_line[i]);
    j = i+loc;

    if (str_eq("on", &config_line[j])) {
      Randomized_Hiff = TRUE;
      printf("evo_bits :: parse_config_file() - randomized hiff is on.\n");

    } else if (str_eq("off", &config_line[j])) {
      Randomized_Hiff = FALSE;
      printf("evo_bits :: parse_config_file() - randomized hiff is off.\n");
    } else {
      printf("evo_bits :: parse_config_file() - Error parsing: %s. (randomized hiff)\n", config_line);
      return FALSE;
    }
    return TRUE;
  }
  


  ////////////////////////////
  // Stats
  if (str_eq("statistics", &config_line[i])) {
    loc = find_field(&config_line[i]);
    j = i+loc;

    if (str_eq("on", &config_line[j])) {
      Keep_Stats = TRUE;
      printf("evo_bits :: parse_config_file() - statistics kept.\n");

    } else if (str_eq("off", &config_line[j])) {
      Keep_Stats = FALSE;
      printf("evo_bits :: parse_config_file() - statistics not kept.\n");
    } else {
      printf("evo_bits :: parse_config_file() - Error parsing: %s. (statistics)\n", config_line);
      return FALSE;
    }
    return TRUE;
  }


  ////////////////////////////
  // Print commands.
  if (str_eq("Print", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if ((i_val < 0) || (i_val >= 50)) {
      printf("evo_bits :: parse_config_file() - Error, 0<= Print[%d] < 50.\n", i_val);
      return FALSE;
    }
    
    Print[i_val] = 1;
    return TRUE;
  }
  if (str_eq("print", config_line)) {
    if (!get_int(config_line, &i_val)) {
      printf("evo_bits :: parse_config_file() - Error parsing: %s\n", config_line);
      return FALSE;
    }

    if ((i_val < 0) || (i_val >= 50)) {
      printf("evo_bits :: parse_config_file() - Error, 0<= Print[%d] < 50.\n", i_val);
      return FALSE;
    }
    
    Print[i_val] = 0;
    return TRUE;
  }

  printf("evo_bits :: parse_config_file() - Error with: %s", config_line);
  return FALSE; // Not a recognized line type.
}


int parse_configure_file(char *fname, int verbose)
{
  int len;
  FILE_HANDLE file;

  file = file_open(fname, "r");
  if (file <= 0) {
    return FALSE;
  }
  
  while ((len = file_read(file, gFile_Buffer, gFILE_BUFFER_LEN))) {
    if (gFile_Buffer[0] == '#') {
      // Ignore comment lines.
      continue;
    }

    parse_config_line(gFile_Buffer, verbose);
  }
  
  file_close(file);

  return TRUE;
}



int main(int argc, char **argv) {
  //  mtrace(); /* This starts memory tracing. */


  Print[1] = TRUE;
  //  Print[1] = FALSE;


  //  seed_random(0);
  seed_random(1);

  gNum_Fitness = 1;

  parse_configure_file(".evo_bits", TRUE);
  if (!process_args(argc, argv)) {
    return 0;
  }


  //  test_bits1();
  //  exit(0);


  if (gReplay) {
    printf("genre :: main() - doing replay of : %s.\n", gFilename1);
    reevaluate_individ(gFilename1);
  } else if (gMutate) {
    mutate_individual(gFilename1);
  } else {
    ea_engine(0);
  }    

  return 0;
}


