#include "common.h"
std::vector<long> sample_indices(int num_of_dusts, int num_of_points)
{
    unsigned char *checks = (unsigned char *) calloc(num_of_points, sizeof(unsigned char)), check;
    long idx;
    std::vector<long> idcs;
    
    idcs.resize(num_of_dusts);
    
    for(long i = 0; i < num_of_dusts; i ++){
        while(1){
            idx = (long) (num_of_points * ((double) random() / RAND_MAX));
             if(!checks[idx]){
                checks[idx] = 1;
                idcs[i] = idx;
                break;
            }

            check = 0;
            for(long l = 0; l < num_of_points; l ++){
                if(!checks[l])
                    check = 1;
            }

            if(!check){
                for(long l = 0; l < num_of_points; l ++)
                    checks[l] = 0;
            }
        }
    }

    free(checks);
    return idcs;
}


