#include "ui.h"

void ui_gen_points(std::string path, long seed)
{
    long num, num_of_clusters = 0, dim;
    double tmp;
    char yn;
    std::vector<long> num_of_elements;
    Points means, vars;

    cout << "Start of gen_points\n";
    cout << "Enter dimension of points:";
    cin >> dim;
    
    while(1){
        cout << "More cluster?(y/n):";
        cin >> yn;

        if(yn == 'n') {
            break;
        } else if (yn != 'y') {
            continue;
        } else {
            cout << "Enter parameters for " << (++ num_of_clusters) << "th cluster.\n";
            
            cout << "  Number of elements:";
            cin >> num;
            
            Point mean, var;
            cout << "  Mean value of cluster\n";
            for(long i = 0; i < dim; i ++){
                cout << "    " << (i + 1) << "th coordinate of mean:";
                cin >> tmp;
                mean.push_back(tmp);
            }

            cout << "  Variance of cluster\n";
            for(long i = 0; i < dim; i ++){
                cout << "    " << (i + 1) << "th coordinate of variance:";
                cin >> tmp;
                var.push_back(tmp);
            }

            num_of_elements.push_back(num);
            means.push_back(mean);
            vars.push_back(var);
        }   
    }

    srand(seed);
    generate_points(path, num_of_elements, means, vars, num_of_clusters, dim);
    printf("End of gen_plain\n");
}


void ui_gen_scheme(std::string path_sk, std::string path_pk, long seed, Ring &ring)
{
    printf("Start of gen_scheme\n");

    srand(seed);
    SecretKey sk(ring);

    generate_secretKey(path_sk, seed, ring);

    read_secretKey(sk, path_sk + ".key");
    generate_scheme(path_pk, sk, ring);

    printf("End of gen_scheme\n");
}


void ui_gen_cvecs(std::string path_enc, std::string path_dat, std::string path_pk, long logp, long logq, long seed, Ring &ring)
{
    std::vector<long> rots, params;
    Points points;
    Scheme scheme(ring, false);
    long dim, num_of_points, num_of_slots = Nh;

    printf("Start of gen_cvecs\n");

    read_params(params, 2, path_dat + ".params");
    dim = params[0];
    num_of_points = params[1];

    read_params(params, 2, path_pk + ".params");
    if((params[0] != logN) || (params[1] != logQ)) {
        printf("ui_gen error: HEAAN params not matches\n");
        return;
    }

    points = read_points(path_dat + ".dat");

    read_scheme(scheme, rots, path_pk);

    srand(seed);
    generate_cvec(path_enc, dim, num_of_points, num_of_slots, logp, logq, points, scheme);

    printf("End of gen_cvecs\n");
}

void ui_cluster_plain(std::string path_dat, long num_of_dusts, long num_of_shiftings, long zeta_ms, long zeta_mi, long gamma_ms, long gamma_nn, long gamma_mi, long seed)
{
    long dim;
    std::vector<long> params;
    Points points, c;

    printf("Start of plain-cluster\n");

    read_params(params, 2, path_dat + ".params");
    dim = params[0];
    points = read_points(path_dat + ".dat");

    srand(seed);
    c = cluster(points, dim, num_of_dusts, num_of_shiftings, zeta_ms, zeta_mi, gamma_ms, gamma_nn, gamma_mi);
    print_plain_result(c, points);

    printf("End of plain-cluster\n");
}

void ui_gaussian_cluster_plain(std::string path_dat, long num_of_dusts, long num_of_shiftings, long zeta_ms, long zeta_mi, double sigma_ms, long gamma_nn, long gamma_mi, long seed)
{
    long dim;
    std::vector<long> params;
    Points points, c;

    printf("Start of plain-gaussian-cluster\n");

    read_params(params, 2, path_dat + ".params");
    dim = params[0];
    points = read_points(path_dat + ".dat");

    srand(seed);
    c = gaussian_cluster(points, dim, num_of_dusts, num_of_shiftings, zeta_ms, zeta_mi, sigma_ms, gamma_nn, gamma_mi);
    print_plain_result(c, points);

    printf("End of plain-gaussian-cluster\n");
}

void ui_fk_cluster_plain(std::string path_dat, long num_of_shiftings, long zeta_ms, long gamma_ms, long seed)
{
    long dim;
    std::vector<long> params;
    Points points;

    printf("Start of plain-fk-cluster\n");

    read_params(params, 2, path_dat + ".params");
    dim = params[0];
    points = read_points(path_dat + ".dat");

    for(long num_of_dusts = 8; num_of_dusts <= 256; num_of_dusts = num_of_dusts << 1){
        printf("Freedman with %ld dusts\n", num_of_dusts);
        srand(seed);
        freedman_cluster(points, dim, num_of_dusts, num_of_shiftings, zeta_ms, gamma_ms);
    }

    printf("End of plain-fk-cluster\n");
}

void ui_cluster_cipher(std::string path_enc, std::string path_pk, std::string path_sk, long num_of_dusts, long num_of_shiftings, long zeta1, long zeta3, long gamma1, long gamma2, long gamma3, long logc, long num_of_threads, long seed, Ring &ring)
{
    Parameter params;
    Scheme scheme(ring, false);
    SecretKey sk(ring);
    CVec P;
    std::vector<CVec> C;
    std::vector<long> pars, rots;
    long dim, num_of_points, logp, logq, logradix;
    long ppc, pvec, unop, udim;

    printf("Start of cipher-cluster\n");

    //Set Parameters
    read_params(pars, 7, path_enc + ".params");
    dim = pars[0];
    num_of_points = pars[1];
    unop = 1 << log2up(num_of_points);
    udim = 1 << log2up(dim);
    logp = pars[3];
    logq = pars[4];
    if((pars[2] != Nh) || (pars[5] != logN) || (pars[6] != logQ)){
        printf("error: ui_cluster path_enc: HEAAN\n");
        return;
    }

    read_params(pars, 2, path_pk + ".params");
    if((pars[0] != logN) || (pars[1] != logQ)){
        printf("error: ui_cluster path_pk\n");
        return;
    }

    logradix = log2up(dim) + log2up(num_of_dusts);
    rotkeys(rots, logradix);
    read_scheme(scheme, rots, path_pk);

    ppc = Nh / udim;
    if(ppc > unop)
        ppc = unop;
    pvec = (num_of_points - 1) / ppc + 1;
    read_cvecs(P, pvec, path_enc);

    set_parameter(params, dim, num_of_dusts, num_of_points, Nh, zeta1, zeta3, gamma1, gamma2, gamma3, num_of_shiftings, logp, logc, logq, num_of_threads, 1 << logradix);
    printf("parameter set\n");


    //Clustering
    srand(seed);
    cluster(C, P, params, scheme); //Note that sk is currently empty
    printf("clustering done\n");


    //Load SecretKey
    read_params(pars, 2, path_sk + ".params");
    if((pars[0] != logN) || (pars[1] != logQ)){
        printf("error: ui_cluster path_sk\n");
        return;
    }
    read_secretKey(sk, path_sk + ".key");
    

    //Decrypt the result
    Point clus;
    for(long i = 0; i < PVEC; i ++){
        for(long j = 0; j < DVEC; j ++){
            Plaintext plain;
            complex<double> *vals;
            scheme.decryptMsg(plain, sk, C[i][j]);
            vals = scheme.decode(plain);
            for(long k = 0; k < UDIM * DPC * PPC; k ++){
                clus.push_back(vals[k].real());
            }
        }
    }

    //Print the result
    double val;
    long idx, idx2;
    for(long i = 0; i < PVEC; i ++){
        for(long j = 0; j < PPC; j ++){
            idx = i * PPC + j;
            if(idx >= NOP)
                break;
            printf("%5ld: ", idx);
            for(long k = 0; k < DVEC; k ++){
                for(long l = 0; l < DPC; l ++){
                    idx2 = k * DPC + l;
                    if(idx2 >= NOD)
                        break;
                    val = clus[i * NOS * DVEC + k * NOS + l * PPC * UDIM + j * UDIM];
                    printf("%f ", (val > 0)? val: 0);
                } 
            }
            printf("\n");
        }
    }

    printf("End of cipher-cluster\n");
}

void ui_cluster_cipher_verbose(std::string path_enc, std::string path_pk, std::string path_sk, long num_of_dusts, long num_of_shiftings, long zeta1, long zeta3, long gamma1, long gamma2, long gamma3, long logc, long num_of_threads, long seed, Ring &ring)
{
    Parameter params;
    Scheme scheme(ring, false);
    SecretKey sk(ring);
    CVec P;
    std::vector<CVec> C;
    std::vector<long> pars, rots;
    long dim, num_of_points, logp, logq, logradix;
    long ppc, pvec, unop, udim;

    printf("Start of cipher-cluster\n");

    //Set Parameters
    read_params(pars, 7, path_enc + ".params");
    dim = pars[0];
    num_of_points = pars[1];
    unop = 1 << log2up(num_of_points);
    udim = 1 << log2up(dim);
    logp = pars[3];
    logq = pars[4];
    if((pars[2] != Nh) || (pars[5] != logN) || (pars[6] != logQ)){
        printf("error: ui_cluster path_enc: HEAAN\n");
        return;
    }

    read_params(pars, 2, path_pk + ".params");
    if((pars[0] != logN) || (pars[1] != logQ)){
        printf("error: ui_cluster path_pk\n");
        return;
    }

    read_params(pars, 2, path_sk + ".params");
    if((pars[0] != logN) || (pars[1] != logQ)){
        printf("error: ui_cluster path_sk\n");
        return;
    }
    read_secretKey(sk, path_sk + ".key");

    logradix = log2up(dim) + log2up(num_of_dusts);
    rotkeys(rots, logradix);
    read_scheme(scheme, rots, path_pk);

    ppc = Nh / udim;
    if(ppc > unop)
        ppc = unop;
    pvec = (num_of_points - 1) / ppc + 1;
    read_cvecs(P, pvec, path_enc);

    set_parameter(params, dim, num_of_dusts, num_of_points, Nh, zeta1, zeta3, gamma1, gamma2, gamma3, num_of_shiftings, logp, logc, logq, num_of_threads, 1 << logradix);
    printf("parameter set\n");


    //Clustering
    srand(seed);
    cluster_verbose(C, P, params, scheme, sk);
    printf("clustering done\n");


    //Print the results
    Point clus;
    for(long i = 0; i < PVEC; i ++){
        for(long j = 0; j < DVEC; j ++){
            Plaintext plain;
            complex<double> *vals;
            scheme.decryptMsg(plain, sk, C[i][j]);
            vals = scheme.decode(plain);
            for(long k = 0; k < UDIM * DPC * PPC; k ++){
                clus.push_back(vals[k].real());
            }
        }
    }

    double val;
    long idx, idx2;
    for(long i = 0; i < PVEC; i ++){
        for(long j = 0; j < PPC; j ++){
            idx = i * PPC + j;
            if(idx >= NOP)
                break;
            printf("%5ld: ", idx);
            for(long k = 0; k < DVEC; k ++){
                for(long l = 0; l < DPC; l ++){
                    idx2 = k * DPC + l;
                    if(idx2 >= NOD)
                        break;
                    val = clus[i * NOS * DVEC + k * NOS + l * PPC * UDIM + j * UDIM];
                    printf("%f ", (val > 0) ? val: 0);
                }                
            }
            printf("\n");
        }
    }

    printf("End of cipher-cluster\n");
}

