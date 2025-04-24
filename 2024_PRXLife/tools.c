#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tools.h"

/* return normal random number using Box-Muller method */
double box_muller(double mu, double sigma)
{
  static int flag = 0;
  static double save = 0.0;

  if (flag == 0)
  {
    double u_1 = 1.0 - genrand_real2();
    double u_2 = genrand_real2();
    double v_1 = mu + sigma * sqrt(-2.0 * log(u_1)) * sin(2.0 * M_PI * u_2);
    double v_2 = mu + sigma * sqrt(-2.0 * log(u_1)) * cos(2.0 * M_PI * u_2);
    save = v_2;
    flag = 1;
    return v_1;
  }
  else
  {
    flag = 0;
    return save;
  }
}

void *my_malloc(size_t size)
{
  void *p;

  p = malloc(size);
  if (p == NULL)
  {
    fprintf(stderr, "malloc fault\n");
    exit(1);
  }
  return p;
}

void my_free(void *ptr)
{
  free(ptr);
  ptr = NULL;
}

/* calculation avarage of x[num] */
double calc_avg(double x[], int num)
{
  int i;
  double avg = 0.0;
  for (i = 0; i < num; i++)
    avg += x[i];
  avg /= (double)num;
  return avg;
}

/* calculation variance of x[num] */
double calc_var(double x[], int num)
{
  int i;
  double avg = 0.0;
  double var = 0.0;
  avg = calc_avg(x, num);
  for (i = 0; i < num; i++)
    var += (x[i] - avg) * (x[i] - avg);
  var /= (double)num;
  return var;
}

/* calculation 3rd moment (skewness) of x[num] */
double calc_3rd_moment(double x[], int num)
{
  int i;
  double avg = 0.0;
  double skew = 0.0;
  avg = calc_avg(x, num);
  for (i = 0; i < num; i++)
    skew += (x[i] - avg) * (x[i] - avg) * (x[i] - avg);
  skew /= (double)num;
  return skew;
}

/* calculation 4th moment (kurtosis) of x[num] */
double calc_4th_moment(double x[], int num)
{
  int i;
  double avg = 0.0;
  double kurt = 0.0;
  avg = calc_avg(x, num);
  for (i = 0; i < num; i++)
    kurt += (x[i] - avg) * (x[i] - avg) * (x[i] - avg) * (x[i] - avg);
  kurt /= (double)num;
  return kurt;
}

/* Mersenne twister */

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti = N + 1;     /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
  mt[0] = s & 0xffffffffUL;
  for (mti = 1; mti < N; mti++)
  {
    mt[mti] =
        (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    mt[mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
  int i, j, k;
  init_genrand(19650218UL);
  i = 1;
  j = 0;
  k = (N > key_length ? N : key_length);
  for (; k; k--)
  {
    mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL)) + init_key[j] + j; /* non linear */
    mt[i] &= 0xffffffffUL;                                                             /* for WORDSIZE > 32 machines */
    i++;
    j++;
    if (i >= N)
    {
      mt[0] = mt[N - 1];
      i = 1;
    }
    if (j >= key_length)
      j = 0;
  }
  for (k = N - 1; k; k--)
  {
    mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) - i; /* non linear */
    mt[i] &= 0xffffffffUL;                                                  /* for WORDSIZE > 32 machines */
    i++;
    if (i >= N)
    {
      mt[0] = mt[N - 1];
      i = 1;
    }
  }

  mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
  unsigned long y;
  static unsigned long mag01[2] = {0x0UL, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */

  if (mti >= N)
  { /* generate N words at one time */
    int kk;

    if (mti == N + 1)       /* if init_genrand() has not been called, */
      init_genrand(5489UL); /* a default initial seed is used */

    for (kk = 0; kk < N - M; kk++)
    {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (; kk < N - 1; kk++)
    {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    mti = 0;
  }

  y = mt[mti++];

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
  return (long)(genrand_int32() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
  return genrand_int32() * (1.0 / 4294967295.0);
  /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
  return genrand_int32() * (1.0 / 4294967296.0);
  /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
  return (((double)genrand_int32()) + 0.5) * (1.0 / 4294967296.0);
  /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
  unsigned long a = genrand_int32() >> 5, b = genrand_int32() >> 6;
  return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */
