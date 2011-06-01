/************************************************************
  Copyright Gregory S. Hornby.

************************************************************/
#ifndef INDIVID_BITS_HEADER_FILE
#define INDIVID_BITS_HEADER_FILE

#include "alps_individ.h"


namespace alps {


extern double Prob_Inversion;
extern double Prob_HTS;
extern double Prob_Invert_HTS;


class Individ_Bits : public alps::Individual
{
 public:
  /********** Methods **********/

  Individ_Bits(int num_genes = 0); // perhaps input data struct
  Individ_Bits(Individual *sample);
  ~Individ_Bits();

  /*** Interface Methods ***/
  Individ_Bits* new_instance();
  void initialize();
  int same(Individ_Bits *individ2);

  void zero_genes(void);
  int get_num_genes(void);
  bool set_num_genes(int num_genes);
  const std::vector<char>& get_genes();

  void shuffle();
  void make_random();
  void duplicate_settings(Individ_Bits *individ2);
  void duplicate_settings(Individual *individ2);
  void duplicate(Individ_Bits *ind);
  void duplicate(Individual *ind);

  bool mutate_virtual_tree(void);
  bool mutate(void);

  bool recombine_tree(Individ_Bits *parent2);
  bool recombine(Individ_Bits *parent2);
  bool recombine(Individual *parent2);
  bool recombine(std::vector<Individual*>& p2s);

  bool make_phenotype(void* arg);
  bool make_phenotype(void);
  int phenotype_distance(Individ_Bits* ind2);
  int phenotype_distance(Individual* ind2);

  std::istream& read(std::istream& file);
  bool read(const char *fname);
  std::ostream& write(std::ostream& file);
  bool write(const char *fname);


 private:
  int num_genes_;
  std::vector<char> gene_vec_;

  int num_hts_points_;
  bool is_shuffled_;
  std::vector<int> gene_loc_; // Shuffling should probably be handled outside class.

  std::vector<int> location_;
  std::vector<int> bit_order_;
  std::vector<char> bit_value_;
 };

} // Namespace ALPS.


#endif
