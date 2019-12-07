#pragma once
#include <stdio.h>
#include "common.h"
#include "heaansrc.h"

#define DPC     params.dusts_per_cipher
#define LDPC    params.log_dusts_per_cipher
#define PPC     params.points_per_cipher
#define LPPC    params.log_points_per_cipher

#define ZETA1   params.zeta1
#define ZETA3   params.zeta3
#define GAMMA1  params.gamma1
#define GAMMA2  params.gamma2
#define GAMMA3  params.gamma3

#define DIM     params.dim
#define LDIM    params.log_dim
#define UDIM    params.dim_up

#define NOD     params.num_of_dusts
#define LNOD    params.log_num_of_dusts
#define UNOD    params.num_of_dusts_up

#define NOP     params.num_of_points
#define LNOP    params.log_num_of_points
#define UNOP    params.num_of_points_up

#define NOS     params.num_of_slots
#define LNOS    params.log_num_of_slots

#define NOSH    params.num_of_shiftings

#define PVEC    params.num_of_pvecs
#define DVEC    params.num_of_dvecs

#define NOTH    params.num_of_threads

#define _logp   params.logp
#define _logc   params.logc
#define _logq   params.logq

#define RADIX   params.radix
  
#define COST_MS params.cost_meanshift

struct Parameter
{
    //Params for mean-shift clustering
    long dim, log_dim, dim_up;
    long num_of_dusts, log_num_of_dusts, num_of_dusts_up;
    long num_of_points, log_num_of_points, num_of_points_up;
    long num_of_slots, log_num_of_slots;
    
    long zeta1, zeta3;
    long gamma1, gamma2, gamma3;
    long num_of_shiftings;

    //Params for multi-threads
    long num_of_pvecs, num_of_dvecs;
    long dusts_per_cipher, points_per_cipher;
    long log_dusts_per_cipher, log_points_per_cipher;

    //Params for HEAAN
    long logp, logc, logq;
    long num_of_threads;

    //Params for Bootstrapping
    long radix;

    long cost_meanshift;
};


void set_parameter(Parameter &params, long dim, long num_of_dusts, long num_of_points, long num_of_slots, long zeta1, long zeta3, long gamma1, long gamma2, long gamma3, long num_of_shiftings, long logp, long logc, long logq, long num_of_threads, long radix);
void print_parameter(Parameter &params);
