/************************************************************
  Copyright Gregory S. Hornby.

************************************************************/
#ifndef INDIVID_REAL_HEADER_FILE
#define INDIVID_REAL_HEADER_FILE


#include <vector>

#include "alps_individ.h"


namespace alps {


class Individ_Real : public Individual
{
public:
  /********** Methods **********/

  Individ_Real(int num_genes = 0); // perhaps input data struct
  Individ_Real(Individual *sample);
  ~Individ_Real();

  /*** Interface Methods ***/
  Individ_Real* new_instance(void);
  void initialize(void);
  int same(Individ_Real *individ2);

  void set_minmax(std::vector<double>& min,
		  std::vector<double>& max);
  void set_init_minmax(std::vector<double>& min,
		       std::vector<double>& max);
  void set_minmax(double min, double max);
  void set_init_minmax(double min, double max);
  void zero_genes(void);
  int get_num_genes(void);
  bool set_num_genes(int num_genes);
  const std::vector<double>& get_genes();

  void make_random(void);
  void duplicate_settings(Individ_Real *individ2);
  void duplicate_settings(Individual *individ2);
  void duplicate(Individ_Real *ind);
  void duplicate(Individual *ind);

  bool mutate(double scale);
  bool mutate();
  void get_unique(std::vector<Individual*>& inds2,
		  std::vector<Individ_Real*>& indreals);
  bool mutate(std::vector<Individual*>& inds2);
  void recombine_linear_clip(double dist, std::vector<double>& parent1,
			     std::vector<double>& parent2);
  void recombine_linear(double dist, std::vector<double>& parent1,
			std::vector<double>& parent2);
  bool recombine_rand2(Individ_Real *parent2);
  bool recombine_rand2(Individual *parent2);
  bool recombine(Individ_Real *parent2);
  bool recombine(Individual *parent2);
  bool recombine(std::vector<Individual*>& inds2);

  bool make_phenotype(void* arg);
  bool make_phenotype(void);
  int phenotype_distance(Individ_Real* ind2);
  int phenotype_distance(Individual* ind2);

  std::istream& read(std::istream& file);
  bool read(const char *fname);
  std::ostream& write(std::ostream& file);
  bool write(const char *fname);


 private:
  int num_genes_;
  std::vector<double> gene_vec_;
  std::vector<double> max_val_;
  std::vector<double> min_val_;
  std::vector<double> init_max_;
  std::vector<double> init_min_;
};

}

#endif
