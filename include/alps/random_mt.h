#ifndef RANDOM_MT_HEADER_FILE
#define RANDOM_MT_HEADER_FILE

namespace alps_random {

void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
double genrand_double1(void);
double genrand_double2(void);

bool is_random_seed_set();
unsigned long get_random_seed();
void seed_random(unsigned long seed);
unsigned int random_int(int max);
unsigned int random_integer();
double random_double();
double random_norm();



}

#endif
