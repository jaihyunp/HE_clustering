#include "ui.h"
#include <MemoryUsage.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    std::string path_enc = "", path_sk = "", path_pk = "", path_dat = "";
    std::string root_data = "data/", dataset = "", root_scheme = "data/SCHEME/";
    long num_of_shiftings = 5, num_of_dusts = 8, logp = 38, logc = 25, num_of_threads = 8;
    long zeta_ms = 5, zeta_mi = 8;
    long gamma_ms = 6, gamma_nn = 7, gamma_mi = 7;
    long seed = 7;
    long mode = 0;
    long verbose = 0;
    double sigma_ms = 0.08;
    
    for(int i = 0; i < argc; i ++){
        if((!strcmp(argv[i], "--path-plain")) && (i < (argc - 1))) {

            path_dat = argv[++ i];

        } else if((!strcmp(argv[i], "--path-cipher")) && (i < (argc - 1))) {

            path_enc = argv[++ i];

        } else if((!strcmp(argv[i], "--path-sk")) && (i < (argc - 1))) {

            path_sk = argv[++ i];

        } else if((!strcmp(argv[i], "--path-pk")) && (i < (argc - 1))) {

            path_pk = argv[++ i];

        } else if((!strcmp(argv[i], "--dataset")) && (i < (argc - 1))) {

            dataset = argv[++ i];
            path_dat = root_data + "plain/" + dataset + "/";
            path_enc = root_data + "cipher/" + dataset + "/cipher";
            path_sk = root_scheme + "sk";
            path_pk = root_scheme + "pk/";

        } else if((!strcmp(argv[i], "--params")) && (i < (argc - 1))) {

            num_of_shiftings = atol(argv[++ i]);
            gamma_ms = atol(argv[++ i]);
            gamma_nn = atol(argv[++ i]);
            zeta_ms = atol(argv[++ i]);
            zeta_mi = atol(argv[++ i]);
            gamma_mi = atol(argv[++ i]);

        } else if((!strcmp(argv[i], "--num-of-shiftings")) && (i < (argc - 1))) {

            num_of_shiftings = atol(argv[++ i]);

        } else if((!strcmp(argv[i], "--num-of-dusts")) && (i < (argc - 1))) {

            num_of_dusts = atol(argv[++ i]);

        } else if((!strcmp(argv[i], "--num-of-threads")) && (i < (argc - 1))) {

            num_of_threads = atol(argv[++ i]);

        } else if((!strcmp(argv[i], "--logp")) && (i < (argc - 1))) {

            logp = atol(argv[++ i]);

        } else if((!strcmp(argv[i], "--logc")) && (i < (argc - 1))) {

            logc = atol(argv[++ i]);

        } else if((!strcmp(argv[i], "--zeta")) && (i < (argc - 2))) {

            zeta_ms = atol(argv[++ i]);
            zeta_mi = atol(argv[++ i]);

        } else if((!strcmp(argv[i], "--gamma")) && (i < (argc - 3))) {

            gamma_ms = atol(argv[++ i]);
            gamma_nn = atol(argv[++ i]);
            gamma_mi = atol(argv[++ i]);

        } else if((!strcmp(argv[i], "--seed")) && (i < (argc - 1))) {

            seed = atol(argv[++ i]);

        } else if(!strcmp(argv[i], "--verbose")) {

            verbose = 1;

        } else if((!strcmp(argv[i], "--mode")) && (i < (argc - 1))) {

            mode = atol(argv[++ i]);

        } else if((!strcmp(argv[i], "--sigma")) && (i < (argc - 1))) {

            sigma_ms = atof(argv[++ i]);

        } else {}
    }


    if(mode & 0x1)
        ui_gen_points(path_dat, seed & 0x10048);
    if(mode & 0x2)
        ui_cluster_plain(path_dat, num_of_dusts, num_of_shiftings, zeta_ms, zeta_mi, gamma_ms, gamma_nn, gamma_mi, seed);
    if(mode & 0x20)
        ui_gaussian_cluster_plain(path_dat, num_of_dusts, num_of_shiftings, zeta_ms, zeta_mi, sigma_ms, gamma_nn, gamma_mi, seed);
    if(mode & 0x40)
        ui_fk_cluster_plain(path_dat, num_of_shiftings, zeta_ms, gamma_ms, seed);

    if(mode & 0x1C) {
        
        Ring ring = Ring();
        SetNumThreads(num_of_threads);

        if(mode & 0x4)
            ui_gen_scheme(path_sk, path_pk, seed & 1310754, ring);

        if(mode & 0x8)
            ui_gen_cvecs(path_enc, path_dat, path_pk, logp, logQ, seed, ring);

        if(mode & 0x10){
            if(verbose)
                ui_cluster_cipher_verbose(path_enc, path_pk, path_sk, num_of_dusts, num_of_shiftings, zeta_ms, zeta_mi, gamma_ms, gamma_nn, gamma_mi, logc, num_of_threads, seed, ring);
            else
                ui_cluster_cipher(path_enc, path_pk, path_sk, num_of_dusts, num_of_shiftings, zeta_ms, zeta_mi, gamma_ms, gamma_nn, gamma_mi, logc, num_of_threads, seed, ring);
        }

    }

    size_t peakAfterSchemeSize = getPeakRSS() >> 20;
    cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

    printf("end of script\n");
    return 1;
}
