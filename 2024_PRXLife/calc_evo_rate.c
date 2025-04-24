#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tools.h"

#define N 1000         // total number of cells
#define GENOME_NUM 100 // maximam limit of the number of genome
#define GENERATION 100000

typedef struct cell
{
    double fitness;
    double pheno;
    int genome_num;
    int next_genome_num;
    double genome[GENOME_NUM];
    double var_genome;
    double skew_genome;
    double kurt_genome;
} cell_t;

double calc_pheno(const cell_t *indiv)
{
    int i;
    double pheno;
    pheno = 0.0;
    for (i = 0; i < indiv->genome_num; i++)
        pheno += indiv->genome[i];
    pheno /= (double)indiv->genome_num;
    return pheno;
}

double calc_fitness(const cell_t *indiv)
{
    // potential = -sin(x / 5) + 2 sin(2 x / 5)
    // initial (local minimum) = -3.441212
    // global minimum = 11.374583
    // local maximum = 3.441212
    // global maximum = -11.374583 (20.041344)
    return -sin(indiv->pheno / 5.0) + 2.0 * sin(2.0 * indiv->pheno / 5.0);
}

double calc_genome_num(const cell_t *indiv, int init_genome_num)
{
    return init_genome_num;
}

int compare_fitness(const void *a, const void *b)
{
    if (((cell_t *)a)->fitness > ((cell_t *)b)->fitness)
        return -1;
    else
        return 1;
}

int compare_phenotype(const void *a, const void *b)
{
    if (((cell_t *)a)->pheno > ((cell_t *)b)->pheno)
        return -1;
    else
        return 1;
}

void circle_avg(double *avg, double *var, const double data[], int num)
{
    int i;
    double x, y;
    x = 0.0;
    y = 0.0;
    for (i = 0; i < num; i++)
    {
        x += cos(data[i]);
        y += sin(data[i]);
    }
    x /= (double)num;
    y /= (double)num;
    *avg = atan2(y, x);
    *var = 1.0 - sqrt(x * x + y * y);
    return;
}

int main(int argc, char *argv[])
{
    int i, j, k;
    int generation;
    int init_genome_num = 1;
    cell_t *indiv;
    cell_t *next_indiv;
    cell_t next_indiv_candidate;
    double fitness_array[N], growth_array[N], pheno_array[N], var_genome_array[N], skew_genome_array[N], kurt_genome_array[N];
    double avg_fitness;
    double avg_pheno;
    double var_pheno;
    double avg_var_genome;
    double var_var_genome;
    double avg_skew_genome;
    double var_skew_genome;
    double avg_kurt_genome;
    double var_kurt_genome;
    double select_beta = 0.25;
    double growth_sum;
    double mutate_variance = 1.0;
    double plateau_width = 2.0;
    double hole_width = 1.0;
    int seed;
    double rand_num;
    char filepath[1024], filepath_2[1024];
    int init_flag;
    int final_generation = GENERATION;

    if (argc > 1)
    {
        init_genome_num = atoi(argv[1]);
    }
    else
    {
        fprintf(stderr, "Num of argment error\n");
        exit(1);
    }

    FILE *fp, *fp_2;

    indiv = (cell_t *)my_malloc(N * sizeof(cell_t));
    next_indiv = (cell_t *)my_malloc(N * sizeof(cell_t));

    sprintf(filepath_2, "first_evoltion_time_%d.data", init_genome_num);
    fp_2 = fopen(filepath_2, "w");

    for (seed = 0; seed < 10000; seed++)
    {

        init_genrand(seed);

        generation = 0;
        init_flag = 1;

        for (i = 0; i < N; i++)
        {
            indiv[i].genome_num = init_genome_num;
            for (j = 0; j < indiv[i].genome_num; j++)
                indiv[i].genome[j] = -3.441212;
            // indiv[i].pheno = calc_pheno(indiv + i);
            indiv[i].pheno = calc_avg(indiv[i].genome, indiv[i].genome_num);
            indiv[i].var_genome = calc_var(indiv[i].genome, indiv[i].genome_num);
            indiv[i].skew_genome = calc_3rd_moment(indiv[i].genome, indiv[i].genome_num);
            indiv[i].skew_genome /= pow(indiv[i].var_genome, 1.5);
            indiv[i].kurt_genome = calc_4th_moment(indiv[i].genome, indiv[i].genome_num);
            indiv[i].kurt_genome /= indiv[i].var_genome * indiv[i].var_genome;
            indiv[i].kurt_genome -= 3.0;
            indiv[i].fitness = -calc_fitness(indiv + i);
            indiv[i].next_genome_num = calc_genome_num(indiv + i, init_genome_num);
            pheno_array[i] = indiv[i].pheno;
            fitness_array[i] = indiv[i].fitness;
            growth_array[i] = exp(select_beta * indiv[i].fitness);
            var_genome_array[i] = indiv[i].var_genome;
            skew_genome_array[i] = indiv[i].skew_genome;
            kurt_genome_array[i] = indiv[i].kurt_genome;
        }
        growth_sum = 0.0;
        for (i = 0; i < N; i++)
            growth_sum += growth_array[i];
        for (i = 0; i < N; i++)
            growth_array[i] /= growth_sum;
        avg_fitness = calc_avg(fitness_array, N);
        avg_pheno = calc_avg(pheno_array, N);
        var_pheno = calc_var(pheno_array, N);
        avg_var_genome = calc_avg(var_genome_array, N);

        for (generation = 1; generation <= GENERATION; generation++)
        {
            k = 0;
            while (1)
            {
                rand_num = genrand_real2();
                i = 0;
                while (1)
                {
                    if (rand_num < growth_array[i])
                        break;
                    rand_num -= growth_array[i];
                    i++;
                }
                next_indiv_candidate.genome_num = indiv[i].next_genome_num;
                for (j = 0; j < next_indiv_candidate.genome_num; j++)
                {
                    rand_num = indiv[i].genome_num * genrand_real2();
                    next_indiv_candidate.genome[j] = indiv[i].genome[(int)rand_num] + box_muller(0.0, mutate_variance);
                }
                next_indiv_candidate.pheno = calc_avg(next_indiv_candidate.genome, next_indiv_candidate.genome_num);
                next_indiv[k] = next_indiv_candidate;
                next_indiv[k].var_genome = calc_var(next_indiv[k].genome, next_indiv[k].genome_num);
                next_indiv[k].skew_genome = calc_3rd_moment(next_indiv[k].genome, next_indiv[k].genome_num);
                next_indiv[k].skew_genome /= pow(next_indiv[k].var_genome, 1.5);
                next_indiv[k].kurt_genome = calc_4th_moment(next_indiv[k].genome, next_indiv[k].genome_num);
                next_indiv[k].kurt_genome /= next_indiv[k].var_genome * next_indiv[k].var_genome;
                next_indiv[k].kurt_genome -= 3.0;
                next_indiv[k].fitness = -calc_fitness(next_indiv + k);
                next_indiv[k].next_genome_num = calc_genome_num(next_indiv + k, init_genome_num);
                k++;
                if (k >= N)
                    break;
            }

            qsort(next_indiv, N, sizeof(cell_t), compare_fitness);

            for (i = 0; i < N; i++)
            {
                indiv[i] = next_indiv[i];
                pheno_array[i] = indiv[i].pheno;
                fitness_array[i] = indiv[i].fitness;
                growth_array[i] = exp(select_beta * indiv[i].fitness);
                var_genome_array[i] = indiv[i].var_genome;
                skew_genome_array[i] = indiv[i].skew_genome;
                kurt_genome_array[i] = indiv[i].kurt_genome;
            }
            growth_sum = 0.0;
            for (i = 0; i < N; i++)
                growth_sum += growth_array[i];
            for (i = 0; i < N; i++)
                growth_array[i] /= growth_sum;
            avg_fitness = calc_avg(fitness_array, N);
            avg_pheno = calc_avg(pheno_array, N);
            var_pheno = calc_var(pheno_array, N);
            avg_var_genome = calc_avg(var_genome_array, N);
            var_var_genome = calc_var(var_genome_array, N);
            avg_skew_genome = calc_avg(skew_genome_array, N);
            var_skew_genome = calc_var(skew_genome_array, N);
            avg_kurt_genome = calc_avg(kurt_genome_array, N);
            var_kurt_genome = calc_var(kurt_genome_array, N);

            if (init_flag == 1)
            {
                if (avg_fitness > 1.5)
                {
                    fprintf(fp_2, "%d %d %f\n", seed, generation, 1.0 / (double)(generation));
                    fflush(fp_2);
                    init_flag = 0;
                    break;
                }
            }

            if (generation >= final_generation)
            {
                if (init_flag == 1)
                {
                    fprintf(fp_2, "%d -1 -1\n", seed);
                    fflush(fp_2);
                }
                break;
            }
        }
    }

    fclose(fp_2);

    return 0;
}
