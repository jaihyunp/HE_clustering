#include "plain.h"

std::vector<long> freedman_cluster(std::vector<Point> &points, long dim, long num_of_dusts, long num_of_shiftings, long zeta_ms, long gamma_ms)
{
    long num_of_points = points.size();
    double max = (double) dim, tmp, dist;
    std::vector<Point> dusts, dusts2;
    std::vector<long> closest, modes, clus;

    //Sample dusts
    dusts = sample_points(num_of_dusts, points, num_of_points);
    for(long i = 0; i < num_of_dusts; i ++){
        Point dust;
        dust.resize(dim);
        for(long j = 0; j < dim; j ++)
            dust[j] = dusts[i][j];
        dusts2.push_back(dust);
    }
    
    //Map backwards
    closest.resize(num_of_points);
    for(long i = 0; i < num_of_points; i ++){
        closest[i] = 0;
        tmp = 9999;
        for(long j = 1; j < num_of_dusts; j ++){
            dist = euclidean_distance_sqr(points[i], dusts[j]);
            if(tmp > dist){
                tmp = dist;
                closest[i] = j;
            }
        }
    }

    //Mean-shift
    for(long i = 0; i < num_of_shiftings; i ++)
        dusts = meanshift(dusts, dusts2, max, zeta_ms, gamma_ms);

    //Monitor the dusts
    printf("Detected modes:\n");
    for(long i = 0; i < num_of_dusts; i ++){
        printf("%3ld: ( ", i);
        for(long j = 0; j < dim; j ++)
            printf("%f ", dusts2[i][j]);
        printf(")->(");
        for(long j = 0; j < dim; j ++)
            printf("%f ", dusts[i][j]);
        printf(")\n");
    }

    //Set eps roughly
    double eps = 0;
    for(long i = 0; i < num_of_dusts; i ++){
        for(long j = 0; j < num_of_dusts; j ++){
            eps += euclidean_distance_sqr(dusts[i], dusts[j]);
        }
    }
    eps /= num_of_dusts * num_of_dusts * 8;
    printf("Epsilon: %f\n", eps);

    //Generate labels of dusts using eps nbhd
    modes.resize(num_of_dusts);
    for(long i = 0; i < num_of_dusts; i ++){
        modes[i] = i;
        for(long j = 0; j < i; j ++){
            if(euclidean_distance_sqr(dusts[i], dusts[j]) < eps){
                modes[i] = modes[j];
                break;
            }
        }
        printf("%ld: %ld\n", i, modes[i]);
    }


    //Print the labels of the data
    clus.resize(num_of_points);
    for(long i = 0; i < num_of_points; i ++){
        printf("%ld: %ld\n", i, clus[i] = modes[closest[i]]);
    }


    //Calculate accuracy for the Tetra dataset
    int candidates[400], cnt = 0, cand;
    for(long i = 0; i < 4; i ++){
        cand = 0;
        for(long j = 0; j < 400; j ++)
            candidates[j] = 0;

        for(long j = 0; j < 100; j ++){
            candidates[modes[closest[100 * i + j]]] ++;
        }

        for(long j = 1; j < 400; j ++){
            if(candidates[cand] < candidates[j])
                cand = j;
        }
        
        for(long j = 0; j < 100; j ++){
            if(modes[closest[100 * i + j]] != modes[cand]){
                cnt ++;
            }
        }
        printf("after %ld cluster (cand: %d): %d\n", i, cand, cnt);
    }
    printf("%d / 400 errors\n", cnt);

    return clus;
}

