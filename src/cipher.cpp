#include "cipher.h"

void left_rotate_sum(Ciphertext & res, Ciphertext &val, long len_of_seg, long log_num_of_iter, Scheme &scheme)
{
    Ciphertext tmp, _res;

    _res.copy(val);

    for(long i = 0; i < log_num_of_iter; i ++){
        scheme.leftRotateFast(tmp, _res, (1 << i) * len_of_seg);
        scheme.addAndEqual(_res, tmp);
    }

    res.copy(_res);
}

void right_rotate_sum(Ciphertext &res, Ciphertext &val, long len_of_seg, long log_num_of_iter, Scheme &scheme)
{
    Ciphertext tmp, _res;

    _res.copy(val);

    for(long i = 0; i < log_num_of_iter; i ++){
    
        scheme.rightRotateFast(tmp, _res, (1 << i) * len_of_seg);
        scheme.addAndEqual(_res, tmp);
        
    }
    
    res.copy(_res);
}

void set_mask(Ciphertext &mask, long len_of_cycle, long fill_starts, long fill_ends, long n, long logp, long logq, Scheme &scheme)
{
    double *vals = (double *) calloc(len_of_cycle, sizeof(double));
	Plaintext plain;

    for(long j = fill_starts; j < fill_ends; j ++)
        vals[j] = 1;

    scheme.encode(plain, vals, len_of_cycle, logp, logq);
    scheme.encryptMsg(mask, plain);
    mask.n = n;
	
	free(vals);
}

void mult_by_mask(Ciphertext &out, Ciphertext &in, double *vals, long vals_size, long logc, Scheme &scheme)
{
    ZZ *mask = new ZZ[N];
    scheme.ring.encode(mask, vals, vals_size, logc);
    scheme.multByPoly(out, in, mask, logc);
}

void kernel(Ciphertext &res, Ciphertext &c1, Ciphertext &c2, double max, double GAMMA, long log_dim, long logp, long logc, Scheme &scheme)
{
    Ciphertext _c1, _c2, tmp;


    if(c1.logq < c2.logq) {
        _c1.copy(c1);
        scheme.modDownTo(_c2, c2, c1.logq);
    } else {
        scheme.modDownTo(_c1, c1, c2.logq);
        _c2.copy(c2);
    }
    
    if(_c1.n > _c2.n)
        scheme.sub(tmp, _c1, _c2);
    else
        scheme.sub(tmp, _c2, _c1);    //tmp <- +-(c1 - c2)

    scheme.multAndEqual(tmp, tmp);
    scheme.reScaleByAndEqual(tmp, logp);   //tmp <- (c1 - c2) ^ 2
    
    left_rotate_sum(tmp, tmp, 1, log_dim, scheme);    //tmp <- distance ^2 [1 / dim]
    
    
    scheme.multByConst(tmp, tmp, - 1.0 / max, logc);
    scheme.reScaleByAndEqual(tmp, logc);
    scheme.addConstAndEqual(tmp, 1, -1);  //tmp <- 1 - distance ^2 / max [1 / dim]

    for(long i = 0; i < GAMMA; i ++){
        scheme.multAndEqual(tmp, tmp);
        scheme.reScaleByAndEqual(tmp, logp);   //tmp <- (1 - distance ^ 2 / max) ^ (2 ^ (i + 1)) [1 / dim]
    }
    
    
    res.copy(tmp); //res <- (1 - distance ^ 2 / max) ^ (2 ^ GAMMA) [1 / dim]
}


void goldschmidt(Ciphertext &res, Ciphertext &val, double given_max, long zeta, double logp, double logc, Scheme &scheme)
{
    Ciphertext base, prod, tmp;

    scheme.multByConst(base, val, -2.0 / given_max, logc);
    scheme.reScaleByAndEqual(base, logc);
    scheme.addConstAndEqual(base, 1, -1);   //base <- 1 - 2 / max * val

    scheme.addConst(prod, base, 1, -1);
    scheme.modDownByAndEqual(prod, logp);   //prod <- 1 + base

    
    for(long i = 0; i < (zeta - 1); i ++) {
        scheme.multAndEqual(base, base);
        scheme.reScaleByAndEqual(base, logp);   //base <- (1 - 2 / max * val) ^ (2 ^ (i + 1))
        
        scheme.addConst(tmp, base, 1, -1);
        scheme.multAndEqual(prod, tmp);
        scheme.reScaleByAndEqual(prod, logp);  //prod <- PROD from j = 0 to (i + 1) {(1 - 2 / max * val) ^ (2 ^ j)}
    }
    
    scheme.multByConstAndEqual(prod, 2.0 / given_max, logc);
    scheme.reScaleBy(res, prod, logc);  //res <- 1 / max * PROD from i = 0 to (zeta - 1) {(1 - 2 / max * val) ^ (2 ^ i)}
}


void sample_dusts_single(CVec &D, CVec &P, std::vector<long> idcs, Parameter &params, Scheme &scheme, long idx)
{
    long idxx, idxy;
    Ciphertext tmp;
    double *vals = (double *) calloc(UDIM * PPC, sizeof(double));
    double *vals2 = (double *) calloc(UDIM * PPC * DPC, sizeof(double));

    for(long i = 0; i < DPC; i ++) {
        if((DPC * idx + i) < NOD) {
            //Extract idxx-th element of P[idxy]
            idxx = idcs[idx * DPC + i] % PPC;
            idxy = idcs[idx * DPC + i] / PPC;

            for(long j = 0; j < DIM; j ++)
                vals[idxx * UDIM + j] = 1;
            mult_by_mask(tmp, P[idxy], vals, UDIM * PPC, logm, scheme);
            scheme.reScaleByAndEqual(tmp, logm);
            for(long j = 0; j < DIM; j ++)
                vals[idxx * UDIM + j] = 0;

            left_rotate_sum(tmp, tmp, UDIM, LPPC, scheme);

            //Sum up
            for(long j = 0; j < UDIM * PPC; j ++){
                vals2[i * UDIM * PPC + j] = 1;
            }
            
            mult_by_mask(tmp, tmp, vals2, UDIM * PPC * DPC, logm, scheme);
            
            for(long j = 0; j < UDIM * PPC; j ++)
                vals2[i * UDIM * PPC + j] = 0;

            if(i)
                scheme.addAndEqual(D[idx], tmp);
            else
                D[idx].copy(tmp);
        }
    }
    scheme.reScaleByAndEqual(D[idx], logm);

    free(vals);
    free(vals2);
}

void sample_dusts(CVec &D, CVec &P, Parameter &params, Scheme &scheme)
{
    long num_threads_here = (NOTH > DVEC) ? DVEC : NOTH;
    long unit = (DVEC - 1) / NOTH + 1;

    std::vector<long> idcs = sample_indices(NOD, NOP);
    D.resize(DVEC);

    NTL_EXEC_RANGE(num_threads_here, first, last);
    for(long i = first; i < last; i ++){
        long idx;
        for(long j = 0; j < unit; j ++){

            idx = i * unit + j;
            if(idx >= DVEC)
                break;
            
            sample_dusts_single(D, P, idcs, params, scheme, idx);

        }
    }
    NTL_EXEC_RANGE_END;
}


void meanshift_single(CVec &D, CVec &P, Parameter &params, Scheme &scheme, long idx)
{
    double *mask = (double *) calloc(UDIM * PPC, sizeof(double));
    Ciphertext ker, ker_all, normalizer, tmp;
    
    for(long i = 0; i < PPC; i ++)
        mask[i * UDIM] = 1;

    for(long i = 0; i < PVEC; i ++){
        kernel(ker, D[idx], P[i], DIM, GAMMA1, LDIM, _logp, _logc, scheme);

        if(NEED_NOP_ADJUSTMENT && (i == (PVEC - 1))) {
            for(long j = PPC - P_DUMMY; j < PPC; j ++)
                mask[j * UDIM] = 0;
        }
        mult_by_mask(ker, ker, mask, UDIM * PPC, logm, scheme);
        scheme.reScaleByAndEqual(ker, logm);
        if(NEED_NOP_ADJUSTMENT && (i == (PVEC - 1))) {
            for(long j = PPC - P_DUMMY; j < PPC; j ++)
                mask[j * UDIM] = 1;
        }

        right_rotate_sum(ker, ker, 1, LDIM, scheme);

        if(i)
            scheme.addAndEqual(normalizer, ker);
        else
            normalizer.copy(ker);

        scheme.modDownTo(tmp, P[i], ker.logq);
        scheme.multAndEqual(ker, tmp);
        scheme.reScaleByAndEqual(ker, _logp);

        if(i)
            scheme.addAndEqual(ker_all, ker);
        else
            ker_all.copy(ker);            
    }

    left_rotate_sum(normalizer, normalizer, UDIM, LPPC, scheme);    
    goldschmidt(normalizer, normalizer, NOP, ZETA1, _logp, _logc, scheme);

    left_rotate_sum(ker_all, ker_all, UDIM, LPPC, scheme);

    scheme.modDownToAndEqual(ker_all, normalizer.logq);
    scheme.multAndEqual(ker_all, normalizer);
    scheme.reScaleByAndEqual(ker_all, _logp);

    for(long i = 0; i < PPC; i ++)
        mask[i * UDIM] = 0;
    for(long i = 0; i < DIM; i ++)
        mask[i] = 1;
    mult_by_mask(ker_all, ker_all, mask, UDIM * PPC, logm, scheme);
    scheme.reScaleByAndEqual(ker_all, logm);
    free(mask);

    right_rotate_sum(D[idx], ker_all, UDIM, LPPC, scheme);
}

void meanshift(CVec &D, CVec &P, Parameter &params, Scheme &scheme)
{
    long num_threads_here = (NOTH > DVEC) ? DVEC : NOTH;
    long unit = (DVEC - 1) / NOTH + 1;

    NTL_EXEC_RANGE(num_threads_here, first, last);
    for(long i = first; i < last; i ++){
        long idx;
        for(long j = 0; j < unit; j ++){
            
            idx = i * unit + j;
            if(idx >= DVEC)
                break;

            meanshift_single(D, P, params, scheme, idx);

        }
    }
    NTL_EXEC_RANGE_END;    
}

void rearrange_dusts_single(CVec &TMP, CVec &D, Parameter &params, Scheme &scheme, long idx)
{
    double *mask = (double *) calloc(UDIM * PPC * DPC, sizeof(double));

    if(NEED_NOD_ADJUSTMENT && (idx == (DVEC - 1))){
        
        for(long i = 0; i < DPC - D_DUMMY; i ++)
            for(long j = 0; j < DIM; j ++)
                mask[idx * DPC * UDIM + (i * PPC  + i) * UDIM + j] = 1;
     
    } else {
    
        for(long i = 0; i < DPC; i ++)
            for(long j = 0; j < DIM; j ++)
                mask[idx * DPC * UDIM + (i * PPC + i) * UDIM + j] = 1;
    
    }

    mult_by_mask(TMP[idx], D[idx], mask, UDIM * PPC * DPC, logm, scheme);
    free(mask);
    left_rotate_sum(TMP[idx], TMP[idx], PPC * UDIM, LDPC, scheme);
}

void rearrange_dusts(Ciphertext &RD, CVec &D, CVec &TMP, Parameter &params, Scheme &scheme)
{
    long num_threads_here = (NOTH > DVEC) ? DVEC : NOTH;
    long unit = (DVEC - 1) / NOTH + 1;

    TMP.resize(DVEC);
    
    NTL_EXEC_RANGE(num_threads_here, first, last);
    for(long i = first; i < last; i ++){
        long idx;
        for(long j = 0; j < unit; j ++){

            idx = i * unit + j;
            if(idx >= DVEC)
                break;
            
            rearrange_dusts_single(TMP, D, params, scheme, idx);

        }
    }
    NTL_EXEC_RANGE_END;

    RD.copy(TMP[0]);
    for(long i = 1; i < DVEC; i ++)
        scheme.addAndEqual(RD, TMP[i]);

    scheme.reScaleByAndEqual(RD, logm);
}

void rerearrange_dusts_single(CVec &D, Ciphertext &RD, Parameter &params, Scheme &scheme, long idx)
{
    long idxx;
    double *vals1 = (double *) calloc(UDIM * PPC * DPC, sizeof(double));
    double *vals2 = (double *) calloc(UDIM * PPC, sizeof(double));

    for(long i = 0; i < DPC; i ++) {
        idxx = idx * DPC + i;
        if(idxx >= NOD)
            break;

        for(long j = 0; j < DIM; j ++)
            vals1[UDIM * PPC * i + idxx * UDIM + j] = 1;
    }

    mult_by_mask(D[idx], RD, vals1, UDIM * PPC * DPC, logm, scheme);
    free(vals1);
    scheme.reScaleByAndEqual(D[idx], logm);

    left_rotate_sum(D[idx], D[idx], UDIM, LNOD, scheme);

    for(long i = 0; i < DIM; i ++)
        vals2[i] = 1;

    mult_by_mask(D[idx], D[idx], vals2, UDIM * PPC, logm, scheme);
    free(vals2);
    scheme.reScaleByAndEqual(D[idx], logm);

    right_rotate_sum(D[idx], D[idx], UDIM, LPPC, scheme);
}

void rerearrange_dusts(CVec &D, Ciphertext &RD, Parameter &params, Scheme &scheme)
{
    long num_threads_here = (NOTH > DVEC) ? DVEC : NOTH;
    long unit = (DVEC - 1) / NOTH + 1;

    D.resize(DVEC);

    NTL_EXEC_RANGE(num_threads_here, first, last);
    for(long i = first; i < last; i ++){
        long idx;
        for(long j = 0; j < unit; j ++){

            idx = i * unit + j;
            if(idx >= DVEC)
                break;
            
            rerearrange_dusts_single(D, RD, params, scheme, idx);

        }
    }
    NTL_EXEC_RANGE_END;
}

void num_of_nbhds_single(CVec &NBHD, CVec &D, Ciphertext &RD, double *mask1, double *mask2, Parameter &params, Scheme &scheme, long idx)
{
    kernel(NBHD[idx], D[idx], RD, DIM, GAMMA2, LDIM, _logp, _logc, scheme);

    mult_by_mask(NBHD[idx], NBHD[idx], mask1, UDIM * PPC, logm, scheme);
    scheme.reScaleByAndEqual(NBHD[idx], logm);

    left_rotate_sum(NBHD[idx], NBHD[idx], UDIM, LNOD, scheme);

    mult_by_mask(NBHD[idx], NBHD[idx], mask2, UDIM * PPC, logm, scheme);
    scheme.reScaleByAndEqual(NBHD[idx], logm);

    right_rotate_sum(NBHD[idx], NBHD[idx], UDIM, LPPC, scheme);
}


void num_of_nbhds(CVec &NBHD, CVec &D, Ciphertext &RD, Parameter &params, Scheme &scheme)
{
    double *mask1 = (double *) calloc(UDIM * PPC, sizeof(double));
    double *mask2 = (double *) calloc(UDIM * PPC, sizeof(double));
    long num_threads_here = (NOTH > DVEC) ? DVEC : NOTH;
    long unit = (DVEC - 1) / NOTH + 1;

    for(long i = 0; i < NOD; i ++)
        for(long j = 0; j < DIM; j ++)
            mask1[i * UDIM + j] = 1;
    mask2[0] = 1;
     
    NTL_EXEC_RANGE(num_threads_here, first, last);
    for(long i = first; i < last; i ++){
        long idx;
        for(long j = 0; j < unit; j ++){

            idx = i * unit + j;
            if(idx >= DVEC)
                break;
            
            num_of_nbhds_single(NBHD, D, RD, mask1, mask2, params, scheme, idx);

        }
    }
    NTL_EXEC_RANGE_END;

    free(mask1);
    free(mask2);
}

void labeling_points_kernel(std::vector<CVec> &C, CVec &P, CVec &D, double *mask1, Parameter &params, Scheme &scheme, long idx)
{
    for(long i = 0; i < PVEC; i ++) {

        kernel(C[i][idx], D[idx], P[i], DIM, GAMMA3, LDIM, _logp, _logc, scheme);

        if(NEED_NOD_ADJUSTMENT) {
            if(idx == (DVEC - 1)) {
                mult_by_mask(C[i][idx], C[i][idx], mask1, UDIM * PPC * DPC, logm, scheme);
                scheme.reScaleByAndEqual(C[i][idx], logm);
            } else {
                scheme.modDownByAndEqual(C[i][idx], logm);
            }
        }
    }
}

void labeling_points_weight(std::vector<CVec> &C, CVec &NBHD, CVec &normalizer, Parameter &params, Scheme &scheme, long idx)
{
    for(long i = 0; i < PVEC; i ++){
        if(NBHD[idx].logq > C[i][idx].logq)
            scheme.modDownToAndEqual(NBHD[idx], C[i][idx].logq);
        else
            scheme.modDownToAndEqual(C[i][idx], NBHD[idx].logq);
        
        scheme.multAndEqual(C[i][idx], NBHD[idx]);
        scheme.reScaleByAndEqual(C[i][idx], _logp);

        scheme.modDownToAndEqual(C[i][idx], normalizer[i].logq);
        scheme.multAndEqual(C[i][idx], normalizer[i]);
        scheme.reScaleByAndEqual(C[i][idx], _logp);
    }
}


void labeling_points(std::vector<CVec> &C, CVec &P, CVec &NBHD, CVec &D, Parameter &params, Scheme &scheme)
{
    long num_threads_here = (NOTH > DVEC) ? DVEC : NOTH;
    long unit = (DVEC - 1) / NOTH + 1;
    double *mask1 = NULL, *mask2 = NULL;
    CVec normalizer;
    
    normalizer.resize(PVEC);
    C.resize(PVEC);
    for(long i = 0; i < PVEC; i ++)
        C[i].resize(DVEC);

    if(NEED_NOD_ADJUSTMENT) {
        mask1 = (double *) calloc(UDIM * PPC * DPC, sizeof(double));
        for(long i = 0; i < (DPC - D_DUMMY); i ++)
            for(long j = 0; j < PPC; j ++)
                mask1[i * UDIM * PPC + j * UDIM] = 1;
    }

    if(DPC > 1) {
        mask2 = (double *) calloc(UDIM * PPC * DPC, sizeof(double));
        for(long i = 0; i < (UDIM * PPC); i ++)
            mask2[i] = 1;
    }
 

    NTL_EXEC_RANGE(num_threads_here, first, last);
    for(long i = first; i < last; i ++) {
        long idx;
        for(long j = 0; j < unit; j ++) {
            idx = i * unit + j;
            if(idx >= DVEC)
                break;
            labeling_points_kernel(C, P, D, mask1, params, scheme, idx);
        }
    }
    NTL_EXEC_RANGE_END;    

    long num_threads_here2 = (NOTH > PVEC) ? PVEC : NOTH;
    long unit2 = (PVEC - 1) / NOTH + 1;    
    NTL_EXEC_RANGE(num_threads_here2, first, last);
    for(long i = first; i < last; i ++){
        long idx;
        for(long j = 0; j < unit2; j ++){
            idx = i * unit2 + j;
            if(idx >= PVEC)
                break;
            normalizer[idx].copy(C[idx][0]);
            for(long k = 1; k < DVEC; k ++)
                scheme.addAndEqual(normalizer[idx], C[idx][k]);
        }
    }
    NTL_EXEC_RANGE_END;

    if(DPC > 1) {//it means PVEC == 1, so does not need multi-threading
        for(long i = 0; i < PVEC; i ++){
            left_rotate_sum(normalizer[i], normalizer[i], UDIM * PPC, LDPC, scheme);
            mult_by_mask(normalizer[i], normalizer[i], mask2, UDIM * PPC * DPC, logm, scheme);
            scheme.reScaleByAndEqual(normalizer[i], logm);
            right_rotate_sum(normalizer[i], normalizer[i], UDIM * PPC, LDPC, scheme);
        }
    }
 
    NTL_EXEC_RANGE(num_threads_here2, first, last);
    for(long i = first; i < last; i ++){
        long idx;
        for(long j = 0; j < unit2; j ++){
            idx = i * unit2 + j;
            if(idx >= PVEC)
                break;
            goldschmidt(normalizer[idx], normalizer[idx], NOD, ZETA3, _logp, _logc, scheme);
        }
    }
    NTL_EXEC_RANGE_END;
    

    NTL_EXEC_RANGE(num_threads_here, first, last);
    for(long i = first; i < last; i ++) {
        long idx;
        for(long j = 0; j < unit; j ++) {
            idx = i * unit + j;
            if(idx >= DVEC)
                break;
            labeling_points_weight(C, NBHD, normalizer, params, scheme, idx);
        }
    }
    NTL_EXEC_RANGE_END;    


    if(NEED_NOD_ADJUSTMENT)
        free(mask1);
    if(DPC > 1)
        free(mask2);
    for(long i = 0; i < PVEC; i ++)
        normalizer[i].free();
}

void cluster(std::vector<CVec> &C, CVec &P, Parameter &params, Scheme &scheme)
{
    Ciphertext RD;
    CVec NBHD, D;
    BootHelper boothelper(LDIM + LNOD, RADIX, _logc + 10, scheme, scheme.ring);
    TimeUtils time;
    double time_all = 0;

    //Dusts sampling
    time.start("Sampling dusts\n");
    sample_dusts(D, P, params, scheme);
    time.stop("Sampling dusts\n");
    time_all += time.timeElapsed;

    //Mode seeking
    for(long shift = 0; shift < NOSH; shift ++) {
        time.start("Mean shift\n");
        meanshift(D, P, params, scheme);
        time.stop("Mean shift\n");
        time_all += time.timeElapsed;
        printf("meanshift %ld done %ld\n", shift, D[0].logq);

        time.start("Preprocess for Bootstrapping");
        rearrange_dusts(RD, D, NBHD, params, scheme);
        right_rotate_sum(RD, RD, UDIM * UNOD, LPPC - LNOD, scheme);
        RD.n = RADIX;
        time.stop("Preprocess for Bootstrapping");
        time_all += time.timeElapsed;

        time.start("Bootstrapping");
        boothelper.bootstrapping(RD, 45, logQ, 4);
        time.stop("Bootstrapping");
        time_all += time.timeElapsed;
        
        time.start("Postprocess for Bootstrapping");
        RD.n = NOS;
        rerearrange_dusts(D, RD, params, scheme);
        time.stop("Postprocess for Bootstrapping");
        time_all += time.timeElapsed;
    }

    printf("MEANSHIFT done %ld\n", D[0].logq);

    //#NBHD
    time.start("Counting neighberhoods");
    num_of_nbhds(NBHD, D, RD, params, scheme);
    time.stop("Counting neighberhoods");
    time_all += time.timeElapsed;
    printf("nbhd calculated %ld\n", NBHD[0].logq);

    //Labeling points
    time.start("Labeling points");
    labeling_points(C, P, NBHD, D, params, scheme);
    time.stop("Labeling points");
    time_all += time.timeElapsed;
    printf("labeling done %ld\n", C[0][0].logq);

    //Free memories
    for(long i = 0; i < DVEC; i ++) {
        NBHD[i].free();
        D[i].free();
    }
    printf("Total time cost: %f ms\n", time_all);
}


void cluster_verbose(std::vector<CVec> &C, CVec &P, Parameter &params, Scheme &scheme, SecretKey &sk)
{
    //Initialize instances
    Ciphertext RD;
    CVec NBHD, D;
    BootHelper boothelper(LDIM + LNOD, RADIX, _logc + 10, scheme, scheme.ring);
    print_parameter(params);
    complex<double> *msg;
    TimeUtils time;
    double time_all = 0;

    for(long i = 0; i < PVEC; i ++){
        msg = scheme.decrypt(sk, P[i]);
        for(long j = 0; j < PPC; j ++){
            long idx = i * PPC + j;
            if(idx >= NOP)
                break;
            printf("%5ld: ", idx);
            for(long k = 0; k < DIM; k ++){
                printf("%f", msg[j * UDIM + k].real());
            }
            printf("\n");
        }
        free(msg);
    }


    //Dusts sampling
    time.start("Sampling dusts\n");
    sample_dusts(D, P, params, scheme);
    time.stop("Sampling dusts\n");
    time_all += time.timeElapsed;
    
    printf("dust set %ld\n", D[0].logq);
    for(long i = 0; i < DVEC; i ++){
        msg = scheme.decrypt(sk, D[i]);
        for(long j = 0; j < DPC; j ++){
            long idx = i * DPC + j;
            printf("%5ld: ", idx);
            for(long k = 0; k < DIM; k ++){
                printf("%f ", msg[j * UDIM * PPC + k].real());
            }
            printf("\n");
        }
        free(msg);
    }


    //Mode seeking
    for(long shift = 0; shift < NOSH; shift ++) {

        time.start("Mean shift\n");
        meanshift(D, P, params, scheme);
        time.stop("Mean shift\n");
        time_all += time.timeElapsed;

        printf("meanshift %ld done %ld\n", shift, D[0].logq);
        for(long i = 0; i < DVEC; i ++){
            msg = scheme.decrypt(sk, D[i]);
            for(long j = 0; j < DPC; j ++){
                long idx = i * DPC + j;
                printf("%5ld: ", idx);
                for(long k = 0; k < DIM; k ++){
                    printf("%f ", msg[j * UDIM * PPC + k].real());
                }
                printf("\n");
            }
            free(msg);
        }


        time.start("Preprocess for Bootstrapping");
        rearrange_dusts(RD, D, NBHD, params, scheme);
        right_rotate_sum(RD, RD, UDIM * UNOD, LPPC - LNOD, scheme);
        RD.n = RADIX;
        time.stop("Preprocess for Bootstrapping");
        time_all += time.timeElapsed;

        printf("before bootstrapping %ld\n", RD.logq);
        msg = scheme.decrypt(sk, RD);
        for(long i = 0; i < NOD; i ++){
            printf("%5ld: ", i);
            for(long j = 0; j < DIM; j ++)
                 printf("%f ", msg[i * UDIM + j].real());
             printf("\n");
        }
        free(msg);


        time.start("Bootstrapping");
        boothelper.bootstrapping(RD, 45, logQ, 4);
        time.stop("Bootstrapping");
        time_all += time.timeElapsed;

        printf("after bootstrapping %ld\n", RD.logq);
        msg = scheme.decrypt(sk, RD);
        for(long i = 0; i < NOD; i ++){
            printf("%5ld: ", i);
            for(long j = 0; j < DIM; j ++)
                 printf("%f ", msg[i * UDIM + j].real());
             printf("\n");
        }
        free(msg);
        

        time.start("Postprocess for Bootstrapping");
        RD.n = NOS;
        rerearrange_dusts(D, RD, params, scheme);
        time.stop("Postprocess for Bootstrapping");
        time_all += time.timeElapsed;

        printf("MEANSHIFT done %ld\n", D[0].logq);
        for(long i = 0; i < DVEC; i ++){
            msg = scheme.decrypt(sk, D[i]);
            for(long j = 0; j < DPC; j ++){
                long idx = i * DPC + j;
                printf("%5ld: ", idx);
                for(long k = 0; k < DIM; k ++){
                    printf("%f ", msg[j * UDIM * PPC + k].real());
                }
                printf("\n");
            }
            free(msg);
        }
    }


    //#NBHD
    time.start("Counting neighberhoods");
    num_of_nbhds(NBHD, D, RD, params, scheme);
    time.stop("Counting neighberhoods");
    time_all += time.timeElapsed;
    
    printf("nbhd calculated %ld\n", NBHD[0].logq);
    for(long i = 0; i < DVEC; i ++){
        msg = scheme.decrypt(sk, NBHD[i]);
        for(long j = 0; j < DPC; j ++){
            long idx = i * DPC + j;
            if(idx >= NOD)
                break;
            printf("%5ld: ", idx);
            for(long k = 0; k < DIM; k ++){
                printf("%f ", msg[j * UDIM * PPC + k].real());
            }
            printf("\n");
        }
        free(msg);
    }


    //Labeling points
    time.start("Labeling points");
    labeling_points(C, P, NBHD, D, params, scheme);
    time.stop("Labeling points");
    time_all += time.timeElapsed;
    printf("labelnig done %ld\n", C[0][0].logq);

    //Free memories
    for(long i = 0; i < DVEC; i ++) {
        NBHD[i].free();
        D[i].free();
    }
    printf("Total time cost: %f ms\n", time_all);
}
