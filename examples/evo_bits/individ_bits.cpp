/*******************************************************
  Copyright Gregory S. Hornby.
  All Rights Reserved.
********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <malloc.h>
#include "individ_bits.h"
#include "individual.h"

#define INDIVID_BITS_VERSION 1

Real Prob_Inversion = 0.0;
Real Prob_HTS = 0.0;
Real Prob_Invert_HTS = 0.0;

/***********************************************************/

// Individual class functions.

Individ_bits :: Individ_bits(void)
{
  Num_Genes = 0;
  Gene_Vec = 0;
  Connections = 0;
  Num_HTS_Points = 0;
  Inversion = FALSE;
  Location = NULL;
  Bit_Order = 0;
  Bit_Value = 0;
}

Individ_bits :: Individ_bits(Individual *sample)
{
  printf("individ_bits :: individ_bits(individual* sample)\n");
  Num_Genes = ((Individ_bits*)sample)->Num_Genes;
  Gene_Vec = NULL;
  Connections = 0;
  Bit_Order = ((Individ_bits*)sample)->Bit_Order;
  Bit_Value = ((Individ_bits*)sample)->Bit_Value;

  Inversion = FALSE;
  Location = NULL;

  //  initialize();

  if (sample != NULL) {
    Num_Genes = ((Individ_bits*)sample)->Num_Genes;
    Inversion = ((Individ_bits*)sample)->Inversion;
    printf("Inversion = %d (%d).\n", Inversion, 
	   ((Individ_bits*)sample)->Inversion);
  }
}


Individ_bits :: ~Individ_bits()
{
  if (Gene_Vec) {
    free(Gene_Vec);
    Gene_Vec = 0;
  }

  if (Connections) {
    for (int i=0; i<Num_Genes; i++) {
      free(Connections[i]);
    }
    free(Connections);
    Connections = 0;
  }

  if (Location) {
    free(Location);
  }
}


Individ_bits* Individ_bits :: new_instance(void)
{
  //  printf("individ_bits :: new_instance()\n");
  return new Individ_bits();
}


// Initialize variables.
void Individ_bits :: initialize(void)
{
  Individual::initialize();
}



// Malloc space.
int Individ_bits :: create(void)
{
  if (Created) {
    Created = FALSE;
    if (Gene_Vec) {
      free(Gene_Vec);
    }
  }

  if (Num_Genes > 0) {
    Gene_Vec = (char *)malloc(Num_Genes*sizeof(char));
    if (Gene_Vec == 0) {
      printf("individual :: create() - error malloc'ing space.\n");
      return FALSE;
    }

    Location = (int *)malloc(Num_Genes*sizeof(int));
    if (Location == 0) {
      printf("individual :: create() - error malloc'ing space.\n");
      return FALSE;
    }

    /*
    Connections = (Gene_Link_Struct**)malloc(Num_Genes*sizeof(int*));
    for (int i=0; i<Num_Genes; i++) {
      Connections[i] = (Gene_Link_Struct*)malloc(MAX_CONNECTIONS*sizeof(Gene_Link_Struct));
    }
    */

    int num_levels = 1;
    int level_pts = Num_Genes;
    Num_HTS_Points = level_pts;

    // Actually, if Num_Genes is a power of 2 this is 2*Num_Genes - 1.
    while (level_pts > 1) {
      level_pts = level_pts >> 1;
      Num_HTS_Points += level_pts;
      num_levels++;
    }
  }

  if (!Individual::create()) {
    return FALSE;
  }

  Created = TRUE;
  return TRUE;
}


int Individ_bits :: uncreate(void)
{
  if (!Created) {
    return TRUE;
  }

  Created = FALSE;
  free(Gene_Vec);

  Individual::uncreate();

  return TRUE;
}



int Individ_bits :: same(Individ_bits *individ2)
{
  int result = TRUE;

  if (!Created) {
    return FALSE;
  }

  if (!individ2->Created) {
    return FALSE;
  }

  if (Num_Genes != individ2->Num_Genes) {
    printf("individual :: same() - number of genes different %d != %d.\n",
	   Num_Genes, individ2->Num_Genes);
    result = FALSE;
  }

  if (!result) {
    while (1) ;
  }

  return result;
}


void Individ_bits :: zero_bits(void)
{
  if (!Created) {
    return;
  }

  for (int i=0; i<Num_Genes; i++) {
    Gene_Vec[i] = 0;
  }
}  



int Individ_bits :: get_num_bits(void)
{
  return Num_Genes;
}


int Individ_bits :: set_num_bits(int num_bits)
{
  if (num_bits <= 0) {
    return FALSE;
  }

  uncreate();
  Num_Genes = num_bits;
  return create();
}

void Individ_bits :: copy_genes(char* values)
{
  if (!Created) {
    return;
  }

  if (!Inversion) {
    for (int i=0; i<Num_Genes; i++) {
      values[i] = Gene_Vec[i];
    }

  } else {
    for (int i=0; i<Num_Genes; i++) {
      //      values[Location[i]] = Gene_Vec[i];
      values[i] = Gene_Vec[Location[i]];
    }
  }
}



int Individ_bits :: make_random(void)
{
  if (!Created) {
    printf("individual :: make_random() - Error, not created!\n");
    return FALSE;
  }

  Individual::make_random();
  for (int i=0; i<Num_Genes; i++) {
    Gene_Vec[i] = random_int(2);
  }

  /*
  for (int i=0; i<MAX_CONNECTIONS; i++) {
    Graph_Root[i].strength = random_int(3)+1;
    Graph_Root[i].index = random_int(Num_Genes);
  }

  for (int i=0; i<Num_Genes; i++) {
    for (int j=0; j<MAX_CONNECTIONS; j++) {
      Connections[i][j].strength = random_int(3)+1;
      Connections[i][j].index = random_int(Num_Genes-1);
      if (Connections[i][j].index == i) {
	Connections[i][j].index = Num_Genes-1;
      }
    }
  }
  */


  if (Inversion) {
    for (int i=0; i<Num_Genes; i++) {
      Location[i] = i;
    }

    for (int i=0; i<Num_Genes; i++) {
      int loc = random_int(Num_Genes);
      int tmp = Location[i];
      Location[i] = Location[loc];
      Location[loc] = tmp;
    }
  }

  /*
  for (int i=0; i<Num_Genes; i++) {
    printf("%3d ", Location[i]);
  }
  printf("\n");
  */

  return TRUE;
}

int Individ_bits :: duplicate_settings(Individ_bits *individ2)
{
  if (Created) {
    uncreate();
  }
  Num_Genes = individ2->Num_Genes;
  Inversion = individ2->Inversion;
  Bit_Order = individ2->Bit_Order;
  Bit_Value = individ2->Bit_Value;

  return TRUE;
}
int Individ_bits :: duplicate_settings(Individual *individ2)
{
  return duplicate_settings((Individ_bits*)individ2);
}



int Individ_bits :: duplicate(Individ_bits *ind)
{
  if (Num_Genes != ind->Num_Genes) {
    uncreate();
  }
  if (!Created) {
    Num_Genes = ind->Num_Genes;
    if (!create()) {
      printf("Individ_bits :: duplicate() - error creating.\n");
      return FALSE;
    }
  }

  Individual::duplicate(ind);

  Inversion = ind->Inversion;

  for (int i=0; i<Num_Genes; i++) {
    Gene_Vec[i] = ind->Gene_Vec[i];
    Location[i] = ind->Location[i];
  }

  /*
  for (int i=0; i<Num_Genes; i++) {
    for (int j=0; j<MAX_CONNECTIONS; j++) {
      Connections[i][j].index = ind->Connections[i][j].index;
      Connections[i][j].strength = ind->Connections[i][j].strength;
    }
  }
  */

  return TRUE;
}
int Individ_bits :: duplicate(Individual *ind)
{
  return duplicate((Individ_bits*)ind);
}


int Individ_bits :: mutate_virtual_tree(void)
{
  // Virtual Sub-Tree Mutation.
  // - works by picking a random subtree and flipping the bits.
  int tree_point = random_int(Num_HTS_Points);

  int level_pts = Num_Genes;
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
  if (stop_pt > Num_Genes) {
    stop_pt = Num_Genes;
  }

  if (!Inversion) {
    for (int i=start_pt; i<stop_pt; i++) {
      if (Gene_Vec[i] == 1) {
	Gene_Vec[i] = 0;
      } else {
	Gene_Vec[i] = 1;
      }
    }

  } else {
    for (int i=start_pt; i<stop_pt; i++) {
      if (Gene_Vec[Location[i]] == 1) {
	Gene_Vec[Location[i]] = 0;
      } else {
	Gene_Vec[Location[i]] = 1;
      }
    }
  }

  return TRUE;
}

int Individ_bits :: mutate_inversion(void)
{
  if ((random_real() >= Prob_Invert_HTS) || (Num_HTS_Points < 3)) {
    int num = random_int(4) + 1;
    for (int i=0; i<num; i++) {
      int index1 = random_int(Num_Genes);
      int index2 = random_int(Num_Genes-1);
      if (index1 == index2) {
	index2 = Num_Genes-1;
      }

      int tmp = Location[index1];
      Location[index1] = Location[index2];
      Location[index2] = tmp;
    }

  } else {
    int hts_point = random_int(Num_HTS_Points-2);

    int level_pts = Num_Genes;
    int level = 1;
    int level_size = 1;
    while (hts_point >= level_pts) {
      hts_point -= level_pts;
      level_pts = level_pts >> 1;
      level++;
      level_size = level_size << 1;
    }

    int hts_point2;
    if (level_pts == 2) {
      if (hts_point == 0) {
	hts_point2 = 1;
      } else {
	hts_point = 0;
      }

    } else if (level_pts > 2) {
      hts_point2 = random_int(level_pts-2);
      if (hts_point2 == hts_point) {
	hts_point2 = level_pts - 1;
      } else {
	if (hts_point % 2) {
	  if (hts_point2 == hts_point-1) {
	    hts_point2 = level_pts - 2;
	  } else if (hts_point2 == hts_point+1) {
	    hts_point2 = level_pts - 2;
	  }
	}
      }

    } else {
      printf("individ_bits :: mutate_invert() - error, invalid number of level points: %d.\n",
	     level_pts);
      while (1) ;
    }
  }

  return TRUE;
}


int Individ_bits :: mutate(void)
{
  int mutated = FALSE;

  //  write("mutating.ind");

  //  printf("ibits :: mutate()\n");

  if (Gene_Vec == NULL) {
    return FALSE;
  }

  if (Inversion) {
    if (random_real() < Prob_Inversion) {
      mutate_inversion();
    }
  }


  if (random_real() >= Prob_HTS) {
    int num = random_int(4) + 1;
    for (int i=0; i<num; i++) {
      int index = random_int(Num_Genes);
      if (Gene_Vec[index]) {
	Gene_Vec[index] = 0;
      } else {
	Gene_Vec[index] = 1;
      }
      mutated = TRUE;
    }

  } else {
    mutated = mutate_virtual_tree();
  }

  if (mutated) {
    Individual::mutate();
  }

  return mutated;
}





int Individ_bits :: recombine_graph(Individ_bits *parent2)
{
  // Part I: pick random point to recombine.
  int num_levels = 1;
  int level_pts = Num_Genes;
  int num_points = level_pts;

  // Actually, if Num_Genes is a power of 2 this is 2*Num_Genes - 1.
  while (level_pts > 1) {
    level_pts = level_pts >> 1;
    num_points += level_pts;
    num_levels++;
  }

  int recomb_point = random_int(num_points);


  // Part II: pick root of graph.
  // First module:
  int sum_total = 0;
  for (int i=0; i<MAX_CONNECTIONS; i++) {
    sum_total += Graph_Root[i].strength;
  }

  int sum = random_int(sum_total);
  int root1 = -1;
  for (int i=0; MAX_CONNECTIONS; i++) {
    if (sum < Graph_Root[i].strength) {
      root1 = i;
      break;
    }
    sum -= Graph_Root[i].strength;
  }

  // Move index to use to position 0.
  // Moving the index used to the front helps in recombining the weights
  // so we know which are worth keeping.
  root1 = Graph_Root[root1].index;
  if (root1 > 0) {
    int tmp_index = Graph_Root[0].index;
    int tmp_strength = Graph_Root[0].strength;
    for (int i=0; i<root1; i++) {
      Graph_Root[i].index = Graph_Root[i+1].index;
      Graph_Root[i].strength = Graph_Root[i+1].strength;
    }
    Graph_Root[root1].index = tmp_index;
    Graph_Root[root1].strength = tmp_strength;
  }


  // Second module:
  sum_total = 0;
  for (int i=1; i<MAX_CONNECTIONS; i++) {
    sum_total += Graph_Root[i].strength;
  }

  sum = random_int(sum_total);
  int root2 = -1;
  for (int i=1; MAX_CONNECTIONS; i++) {
    if (sum < Graph_Root[i].strength) {
      root2 = i;
      break;
    }
    sum -= Graph_Root[i].strength;
  }

  // Move index to use to position 1.
  root2 = Graph_Root[root2].index;
  if (root2 > 1) {
    int tmp_index = Graph_Root[1].index;
    int tmp_strength = Graph_Root[1].strength;
    for (int i=1; i<root2; i++) {
      Graph_Root[i].index = Graph_Root[i+1].index;
      Graph_Root[i].strength = Graph_Root[i+1].strength;
    }
    Graph_Root[root2].index = tmp_index;
    Graph_Root[root2].strength = tmp_strength;
  }

  


  // Part III: Traverse down the graph, randomly selecting the branches.
  int level = 1;
  int level_size = Num_Genes;
  level_pts = 1;
  while (recomb_point >= level_pts) {
    recomb_point -= level_pts;
    level_pts = level_pts << 1;
    level++;
    level_size = level_size >> 1;
  }

  int start_pt = level_size*recomb_point;
  int stop_pt = start_pt + level_size;
  if (stop_pt > Num_Genes) {
    stop_pt = Num_Genes;
  }


  
  



  // Part IV: Recombine the connection strengths.
  // First index cannot decrease.
  //  other indices, can they increase, or only decrease?
  //  perhaps fall in range of p1 & p2.


  return TRUE;
}



int Individ_bits :: recombine_tree(Individ_bits *parent2)
{
  int recomb_point = random_int(Num_HTS_Points);
  //  int recomb_pt2 = recomb_point;

  int level_pts = Num_Genes;
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
  if (stop_pt > Num_Genes) {
    stop_pt = Num_Genes;
  }

  /*
  printf("recombine_tree() - #Genes:%d, #levels:%d, #pts:%d;  recomb-pt:%d(%d), level:%d, level-size:%d, start:%d, stop:%d.\n",
	 Num_Genes, num_levels, num_points, recomb_pt2, recomb_point,
	 level, level_size, start_pt, stop_pt);
  */

  if (!Inversion) {
    for (int i=start_pt; i<stop_pt; i++) {
      Gene_Vec[i] = parent2->Gene_Vec[i];
    }

  } else {
    for (int i=start_pt; i<stop_pt; i++) {
      Gene_Vec[Location[i]] = parent2->Gene_Vec[Location[i]];
    }
  }

  /*
  for (int i=start_pt; i<stop_pt; i++) {
    printf("%d ", (int)Gene_Vec[i]);
  }
  printf("\n");

  for (int i=start_pt; i<stop_pt; i++) {
    printf("%d ", (int)parent2->Gene_Vec[i]);
  }
  printf("\n\n");
  */  

  return TRUE;
}



int Individ_bits :: recombine(Individ_bits *parent2)
{
  if (Num_Genes < 2) {
    return FALSE;
  }

  Individual::recombine(parent2);

  //  recombine_tree(parent2);

  if (Inversion) {
    if (random_real() < Prob_Inversion) {
      mutate_inversion();
    }
  }


  if (random_real() >= Prob_HTS) {
    if (Num_Genes == 2) {
      if (random_int(2)) {
	Gene_Vec[0] = parent2->Gene_Vec[0];
      }
      if (random_int(2)) {
	Gene_Vec[1] = parent2->Gene_Vec[1];
      }

    } else {
      int loc1 = random_int(Num_Genes);
      int loc2 = random_int(Num_Genes-1);

      if (loc2 == loc1) {
	loc2 = Num_Genes - 1;
      }

      if (!Inversion) {
	if (loc1 < loc2) {
	  for (int i=loc1; i<loc2; i++) {
	    Gene_Vec[i] = parent2->Gene_Vec[i];
	  }

	} else {
	  for (int i=loc2; i<loc1; i++) {
	    Gene_Vec[i] = parent2->Gene_Vec[i];
	  }
	}
	
      } else {
	if (loc1 < loc2) {
	  for (int i=loc1; i<loc2; i++) {
	    Gene_Vec[Location[i]] = parent2->Gene_Vec[Location[i]];
	  }

	} else {
	  for (int i=loc2; i<loc1; i++) {
	    Gene_Vec[Location[i]] = parent2->Gene_Vec[Location[i]];
	  }
	}
      }
    }

  } else {
    recombine_tree(parent2);
  }

  return TRUE;
}

int Individ_bits :: recombine(Individual *parent2)
{
  return recombine((Individ_bits*)parent2);
}



int Individ_bits :: phenotype_distance(Individ_bits* ind2)
{
  if ((Created && !ind2->Created) ||
      (!Created && ind2->Created)) {
    return Num_Genes;
  }

  int num_genes = minv(Num_Genes, ind2->Num_Genes);

  int dist = maxv(Num_Genes, ind2->Num_Genes) - num_genes;
  for (int i=0; i<num_genes; i++) {
    if (Gene_Vec[i] != ind2->Gene_Vec[i]) {
      num_genes++;
    }
  }

  return dist;
}
int Individ_bits :: phenotype_distance(Individual* ind2)
{
  return phenotype_distance((Individ_bits*)ind2);
}


int Individ_bits :: make_phenotype(void* arg)
{
  return TRUE;
}

int Individ_bits :: make_phenotype(void)
{
  return TRUE;
}


int Individ_bits :: read(FILE_HANDLE file, int (*translator_func)(char*))
{
  int i, j, num, version;

  // if (Created) {
  //    uncreate();
  // }

  Individual::read(file, translator_func);


  file_read(file, gFile_Buffer, gFILE_BUFFER_LEN);
  num = sscanf(gFile_Buffer, "Individ_bits Version: %d\n", &version);
  if (num != 1) {
    printf("Individ_bits: error in reading version\n");
    return FALSE;
  } else if (version != INDIVID_BITS_VERSION) {
    printf("Individ_bits: Error reading unsupported version #%d\n", version);
    return FALSE;
  }    

  Evaluated = 2;
  
  file_read(file, gFile_Buffer, gFILE_BUFFER_LEN);  
  num = sscanf(gFile_Buffer, "%d\n", &Num_Genes);
  if (num != 1) {
    printf("Individ_bits: error in reading Num_Genes\n");
    printf("<%s>\n", gFile_Buffer);
    return FALSE;
  }

  if (Num_Genes < 1) {
    printf("Individ_bits: error in reading Num_Genes: %d < 1.\n",
	   Num_Genes);
    printf("<%s>\n", gFile_Buffer);
    return FALSE;
  }    

  create();


  if (!file_read(file, gFile_Buffer, gFILE_BUFFER_LEN)) {
    printf("Individ_bits :: read() - error, empty line.");
    return FALSE;
  }

  j = 0;
  for (i=0; i<Num_Genes; i++) {
    Gene_Vec[i] = read_integer(&(gFile_Buffer[j]));
    while (gFile_Buffer[j] != ' ')
      j++;
    j++;
  }

  if (!file_read(file, gFile_Buffer, gFILE_BUFFER_LEN)) {
    // end of data file.
    return TRUE;
  }

  return TRUE;
}


int Individ_bits :: read(char *fname, int (*translator_func)(char*))
{
  return Individual::read(fname, translator_func);
}


void Individ_bits :: write(FILE_HANDLE file)
{
  int i, loc;

  Individual::write(file);

  sprintf(gFile_Buffer, "Individ_bits Version: %d\n", INDIVID_BITS_VERSION);
  file_write(file, gFile_Buffer);

  sprintf(gFile_Buffer, "%d\n", Num_Genes);
  file_write(file, gFile_Buffer);

  if (Num_Genes <= 0) {
    printf("Ind::write()  Error - size is %d.\n", Num_Genes);
    return;
  }

  loc = 0;
  for (i=0; i<Num_Genes; i++) {
    //    sprintf(gFile_Buffer+loc, "%d ", Gene_Vec[i]);
    //    while (gFile_Buffer[++loc] != ' ')
    //      ;
    sprintf(gFile_Buffer+loc, "%d", Gene_Vec[i]);
    loc++;
  }
  gFile_Buffer[loc-1] = '\n';
  gFile_Buffer[loc] = '\n';
  gFile_Buffer[loc+1] = 0;
  file_write(file, gFile_Buffer);

  loc = 0;
  for (i=0; i<Num_Genes; i++) {
    sprintf(gFile_Buffer+loc, "%d", Gene_Vec[Location[i]]);
    //    while (gFile_Buffer[++loc] != ' ')
    //      ;
    loc++;
  }
  gFile_Buffer[--loc] = '\n';
  file_write(file, gFile_Buffer);


  if (Bit_Value) {
    sprintf(gFile_Buffer, "DeRandomized: ");
    loc = 14;
    for (i=0; i<Num_Genes; i++) {
      sprintf(gFile_Buffer+loc, "%d", Gene_Vec[i] ^ Bit_Value[i]);
      //    while (gFile_Buffer[++loc] != ' ')
      //      ;
      loc++;
    }
    gFile_Buffer[--loc] = '\n';
    file_write(file, gFile_Buffer);
  }

  if (Inversion) {
    sprintf(gFile_Buffer, "UnShuffled: ");
    loc = 12;
    for (i=0; i<Num_Genes; i++) {
      sprintf(gFile_Buffer+loc, "%d", Gene_Vec[Bit_Order[Location[i]]]);
      //    while (gFile_Buffer[++loc] != ' ')
      //      ;
      loc++;
    }
    gFile_Buffer[--loc] = '\n';
    file_write(file, gFile_Buffer);

  } else if (Bit_Order) {
    sprintf(gFile_Buffer, "UnShuffled: ");
    loc = 12;
    for (i=0; i<Num_Genes; i++) {
      sprintf(gFile_Buffer+loc, "%d", Gene_Vec[Bit_Order[i]]);
      //    while (gFile_Buffer[++loc] != ' ')
      //      ;
      loc++;
    }
    gFile_Buffer[--loc] = '\n';
    file_write(file, gFile_Buffer);
  }


  loc = 0;
  for (i=0; i<Num_Genes; i++) {
    sprintf(gFile_Buffer+loc, "%d ", Location[i]);
    while (gFile_Buffer[++loc] != ' ')
      ;
    loc++;
  }
  gFile_Buffer[--loc] = '\n';
  file_write(file, gFile_Buffer);

  if (Bit_Order) {
    loc = 0;
    for (i=0; i<Num_Genes; i++) {
      sprintf(gFile_Buffer+loc, "%d ", Bit_Order[i]);
      while (gFile_Buffer[++loc] != ' ')
	;
      loc++;
    }
    gFile_Buffer[--loc] = '\n';
    file_write(file, gFile_Buffer);
  }
}


int Individ_bits :: write(char *fname)
{
  return Individual::write(fname);
}

