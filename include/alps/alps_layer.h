/*******************************************************
  Copyright Gregory S. Hornby.
  All Rights Reserved.
********************************************************/

#ifndef ALPS_LAYER_HEADER_FILE
#define ALPS_LAYER_HEADER_FILE

#include <iostream>
#include <vector>


namespace alps {

class AlpsLayer
{
  friend std::ostream& operator<<(std::ostream&, AlpsLayer&);
  friend std::istream& operator>>(std::istream&, AlpsLayer&);

 public:
  AlpsLayer();
  AlpsLayer(int s);
  ~AlpsLayer();

  void set_elitism(unsigned int v);
  void set_parent_layer(int p);
  void set_max_age(int a);
  void set_prev_layers(int n, std::vector<int>& l);
  void set_select_type(int t);
  void set_size(int s);
  void set_start_index(int i);
  void set_prob_select_prev(double p);
  void set_tourn_size(unsigned int t);
  unsigned int get_tourn_size() const;

  unsigned int get_elitism() const;
  int get_max_age() const;
  int get_num_prev_layers() const;
  int get_parent_layer() const;
  int get_prev_layer(unsigned int l) const;
  double get_prob_select_prev() const;
  unsigned int get_size() const;
  int get_select_type() const;
  unsigned int get_start_index() const;
  int get_last_index() const;
  unsigned int get_stop_index() const;

  void print() const;
  void print_info() const;

  std::istream& read(std::istream& file);
  std::ostream& write(std::ostream& file);
  std::ostream& write_header(std::ostream& file);

 private:
  unsigned int elitism_;
  int max_age_;
  int parent_layer_;
  int num_prev_layers_;
  int select_type_;
  unsigned int size_;
  unsigned int start_index_;
  int last_index_; // Last actual index.
  unsigned int stop_index_; // One past last actual index.
  unsigned int tourn_size_;
  double prob_select_prev_;
  std::vector<int> prev_layers_;
};

}

#endif
