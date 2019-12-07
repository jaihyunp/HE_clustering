#pragma once
#include <vector>
#include <stdlib.h>

std::vector<long> sample_indices(int num_of_dusts, int num_of_points);
inline long log2up(long val)
{
    long res = 0, prod = 1;
    while(prod < val){
        res ++;
        prod = prod << 1;
    }
    return res;
};
