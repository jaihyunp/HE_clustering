#pragma once
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "typedefs.h"
#include "heaansrc.h"
#include "parameter.h"

#define P_DUMMY ((PVEC * PPC) - NOP)
#define D_DUMMY ((DVEC * DPC) - NOD)
#define NEED_NOP_ADJUSTMENT (P_DUMMY != 0)
#define NEED_NOD_ADJUSTMENT (D_DUMMY != 0)
#define logm _logc


//void left_rotate(Ciphertext &res, Ciphertext &val, long idx, Scheme &scheme);
//void right_rotate(Ciphertext &res, Ciphertext &val, long idx, Scheme &scheme);


/************************************************************************
 * Elementary functions of HEAAN
 ************************************************************************/

void left_rotate_sum(Ciphertext & res, Ciphertext &val, long len_of_seg, long log_num_of_iter, Scheme &scheme);

void right_rotate_sum(Ciphertext &res, Ciphertext &val, long len_of_seg, long log_num_of_iter, Scheme &scheme);

void set_mask(Ciphertext &mask, long len_of_cycle, long fill_starts, long fill_ends, long n, long logp, long logq, Scheme &scheme);

void mult_by_mask(Ciphertext &out, Ciphertext &in, double *vals, long vals_size, long logp, Scheme &scheme); //TODO mask

void goldschmidt(Ciphertext &res, Ciphertext &val, double given_max, long zeta, double logp, double logc, Scheme &scheme);



/************************************************************************
 * Complex functions for Mean shift Clustering in HEAAN
 ************************************************************************/
//Kernel
void kernel(Ciphertext &res, Ciphertext &c1, Ciphertext &c2, double max, double GAMMA, long log_dim, long logp, long logc, Scheme &scheme);

void sample_dusts(CVec &D, CVec &P, Parameter &params, Scheme &scheme);
void meanshift(CVec &D, CVec &P, Parameter &params, Scheme &scheme);
void rearrange_dusts(Ciphertext &RD, CVec &D, CVec &TMP, Parameter &params, Scheme &scheme);
void rerearrange_dusts(CVec &D, Ciphertext &RD, Parameter &params, Scheme &scheme);
void num_of_nbhds(CVec &NBHD, CVec &D, Ciphertext &RD, Parameter &params, Scheme &scheme);
void labeling_points(std::vector<CVec> &C, CVec &P, CVec &NBHD, CVec &D, Parameter &params, Scheme &scheme);

void cluster(std::vector<CVec> &C, CVec &P, Parameter &params, Scheme &scheme);
void cluster_verbose(std::vector<CVec> &C, CVec &P, Parameter &params, Scheme &scheme, SecretKey &sk);
