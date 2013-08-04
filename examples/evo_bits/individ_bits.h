/************************************************************
  Copyright Gregory S. Hornby.
  All Rights Reserved.
************************************************************/
#ifndef INDIVID_BITS_HEADER_FILE
#define INDIVID_BITS_HEADER_FILE

#include "general.h"
#include "individual.h"


#define MAX_CONNECTIONS 7

class Vec;

extern Real Prob_Inversion;
extern Real Prob_HTS;
extern Real Prob_Invert_HTS;

typedef struct Gene_Link_Struct {
  int index;
  int strength;
};

class Individ_bits : public Individual
{
public:
  int Num_Genes;
  char *Gene_Vec;
  int Num_HTS_Points;

  int Inversion;
  int* Location;
  int* Bit_Order;
  int* Bit_Value;

  Gene_Link_Struct Graph_Root[MAX_CONNECTIONS];
  Gene_Link_Struct **Connections;

  /********** Methods **********/

  Individ_bits(void); // perhaps input data struct
  Individ_bits(Individual *sample);
  ~Individ_bits();

  /*** Interface Methods ***/

  Individ_bits* new_instance(void);
  void initialize(void);
  int create(void);
  int uncreate(void);
  int same(Individ_bits *individ2);

  void zero_bits(void);
  int get_num_bits(void);
  int set_num_bits(int num_bits);
  void copy_genes(char* values);

  int make_random(void);
  int duplicate_settings(Individ_bits *individ2);
  int duplicate_settings(Individual *individ2);
  int duplicate(Individ_bits *ind);
  int duplicate(Individual *ind);

  int mutate_virtual_tree(void);
  int mutate_inversion(void);
  int mutate(void);
  int recombine_graph(Individ_bits *parent2);
  int recombine_tree(Individ_bits *parent2);
  int recombine(Individ_bits *parent2);
  int recombine(Individual *parent2);

  int make_phenotype(void* arg);
  int make_phenotype(void);
  int phenotype_distance(Individ_bits* ind2);
  int phenotype_distance(Individual* ind2);

  int read(FILE_HANDLE file, int (*translator_func)(char*));
  int read(char *fname, int (*translator_func)(char*));
  void write(FILE_HANDLE file);
  int write(char *fname);
};

#endif
