 #include "parameter.h"

/**************************************************************** 
 * struct Parameter()
 * {
 *     //Params for mean-shift clustering
 *     long dim, log_dim, dim_up;
 *     long num_of_dusts, log_num_of_dusts, num_of_dusts_up;
 *     long num_of_points, log_num_of_points, num_of_points_up;
 *     long num_of_slots, log_num_of_slots;
 *     
 *     long zeta1, zeta3;
 *     long gamma1, gamma2, gamma3;
 *     long num_of_shiftings;
 * 
 *     //Params for multi-threads
 *     long num_of_pvecs, num_of_dvecs;
 *     long dusts_per_cipher, points_per_cipher;
 * 
 *     //Params for HEAAN
 *     long logp, logc, logq;
 *     long num_of_threads;
 * 
 *     //Params for Bootstrapping
 *     long radix;
 * }
 ****************************************************************/

void print_parameter(Parameter &params)
{
    printf("Params for mean shift clustering:\n");
    printf("* Dimension: %ld (2**%ld = %ld)\n", params.dim, params.log_dim, params.dim_up);
    printf("* Num of dusts: %ld (2**%ld = %ld)\n", params.num_of_dusts, params.log_num_of_dusts, params.num_of_dusts_up);
    printf("* Num of points: %ld (2**%ld = %ld)\n", params.num_of_points, params.log_num_of_points, params.num_of_points_up);
    printf("* Num of slots: %ld (2**%ld)\n", params.num_of_slots, params.log_num_of_slots);
    printf("* ZETA: %ld %ld\n", params.zeta1, params.zeta3);
    printf("* gamma: %ld %ld %ld\n", params.gamma1, params.gamma2, params.gamma3);
    
    printf("Params for packing\n");
    printf("* %ld pvecs, %ld(2**%ld) points in each\n", params.num_of_pvecs, params.points_per_cipher, params.log_points_per_cipher);
    printf("* %ld dvecs, %ld dusts in each\n", params.num_of_dvecs, params.dusts_per_cipher);

    printf("Params for HEAAN\n");
    printf("* logp = %ld, logc = %ld, logq = %ld\n", params.logp, params.logc, params.logq);
    printf("* logN = %ld, logQ = %ld\n", logN, logQ);
    printf("* Num of threads: %ld\n", params.num_of_threads);

    printf("Params for FFT-Bootstrapping\n");
    printf("* Radix: %ld\n", params.radix);
}

void set_parameter(Parameter &params, long dim, long num_of_dusts, long num_of_points, long num_of_slots, long zeta1, long zeta3, long gamma1, long gamma2, long gamma3, long num_of_shiftings, long logp, long logc, long logq, long num_of_threads, long radix)
{
    params.dim = dim;
    params.log_dim = log2up(dim);
    params.dim_up = 1 << params.log_dim;

    params.num_of_dusts = num_of_dusts;
    params.log_num_of_dusts = log2up(num_of_dusts);
    params.num_of_dusts_up = 1 << params.log_num_of_dusts;

    params.num_of_points = num_of_points;
    params.log_num_of_points = log2up(num_of_points);
    params.num_of_points_up = 1 <<params.log_num_of_points;

    params.num_of_slots = num_of_slots;
    params.log_num_of_slots = log2up(num_of_slots);

    params.zeta1 = zeta1;
    params.zeta3 = zeta3;
    params.gamma1 = gamma1;
    params.gamma2 = gamma2;
    params.gamma3 = gamma3;
    params.num_of_shiftings = num_of_shiftings;

    params.points_per_cipher = num_of_slots / params.dim_up;
    if(params.points_per_cipher > params.num_of_points_up)
        params.points_per_cipher = params.num_of_points_up;
    params.log_points_per_cipher = log2up(params.points_per_cipher);
    params.num_of_pvecs = (num_of_points - 1) / params.points_per_cipher + 1;

    params.dusts_per_cipher = num_of_slots / params.dim_up / params.points_per_cipher;
    if(params.dusts_per_cipher > params.num_of_dusts_up)
        params.dusts_per_cipher = params.num_of_dusts_up;
    params.log_dusts_per_cipher = log2up(params.dusts_per_cipher);
    params.num_of_dvecs = (num_of_dusts - 1) / params.dusts_per_cipher + 1;

    params.logp = logp;
    params.logc = logc;
    params.logq = logq;
    params.num_of_threads = num_of_threads;
    
    params.radix = radix;

    params.cost_meanshift = 550;
}


