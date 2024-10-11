#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED

/* output normal random number using Box-Muller method */
double box_muller(double mu, double sigma);

void *my_malloc(size_t size);

void my_free(void *ptr);

/* calculation avarage of x[num] */
double calc_avg(double x[], int num);

/* calculation variance of x[num] */
double calc_var(double x[], int num);

/* calculation 3rd moment (skewness) of x[num] */
double calc_3rd_moment(double x[], int num);

/* calculation 4th moment (kurtosis) of x[num] */
double calc_4th_moment(double x[], int num);

/* Mersenne twister */

/* Period parameters */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length);

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void);

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void);

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void);

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void);

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void);

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void);

#endif /* ZENBUNOSE11_H_INCLUDED */
