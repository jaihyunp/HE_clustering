#include "plain.h"

double gaussian_kernel(double val, double max, double sigma)
{
    return exp(- val / max / 2.0 / sigma / sigma);
}


double gaussian_kernel(Point &point_a, Point &point_b, double max, double sigma)
{
    return gaussian_kernel(euclidean_distance_sqr(point_a, point_b), max, sigma);
}


std::vector<Point> gaussian_meanshift(std::vector<Point> &dusts, std::vector<Point> &points, double max, long zeta, double sigma)
{
    long dim = (int) dusts[0].size(), num_of_points = points.size(), num_of_dusts = dusts.size();
    double sum, ker;
    std::vector<Point> res;
    
    res.resize(num_of_dusts);
        
    for(long i = 0; i < num_of_dusts; i ++){
    
        sum = 0;
        for(long k = 0; k < dim; k ++)
            res[i].push_back(0);

        for(long j = 0; j < num_of_points; j ++){

            ker = gaussian_kernel(dusts[i], points[j], max, sigma);

            for(long k  = 0; k < dim; k ++)
                res[i][k] += ker * points[j][k];

            sum += ker;
        }
        
        sum = goldschmidt(sum, num_of_points, zeta);

        for(long k = 0; k < dim; k ++)
            res[i][k] *= sum;

    }
    
    return res;
}


std::vector<Point> gaussian_cluster(std::vector<Point> &points, long dim, long num_of_dusts, long num_of_shiftings, long zeta_ms, long zeta_sm, double sigma_ms, long gamma_sm, long gamma_mm)
{
    long num_of_points = (long) points.size();
    double max = (double) dim;
    std::vector<Point> dusts, clus;

    dusts = sample_points(num_of_dusts, points, num_of_points);

    for(long i = 0; i < num_of_shiftings; i ++)
        dusts = gaussian_meanshift(dusts, points, max, zeta_ms, sigma_ms);

    clus = index_numbering(dusts, points, max, zeta_sm, gamma_sm, gamma_mm);

    return clus;
}
