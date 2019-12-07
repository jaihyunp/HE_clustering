#include "plain.h"


std::vector<Point> sample_points(int num_of_dusts, std::vector<Point> &points, int num_of_points)
{
    std::vector<long> idcs = sample_indices(num_of_dusts, num_of_points);
    std::vector<Point> dusts;

    dusts.resize(num_of_dusts);

    for(long i = 0; i < num_of_dusts; i ++){
        Point dust;
        for(long j = 0; j < (long) points[idcs[i]].size(); j ++)
            dust.push_back(points[idcs[i]][j]);

        dusts[i] = dust;

    }

    return dusts;
}


double euclidean_distance_sqr(Point &point_a, Point &point_b)
{
    double res = 0, temp;

    for(int i = 0; i < (int) point_a.size(); i ++){
    
        temp = (point_a[i] - point_b[i]);
        res += temp * temp;
        
    }
    
    return res;
}


double goldschmidt(double val, double given_max, long zeta)
{
    double base = 1.0 - 2.0 * val / given_max;
    double res = 2.0 / given_max;
    
    for(int i = 0; i < zeta; i ++) {
    
        if(i)
            base = base * base;
            
        res *= (1 + base);
        
    }
    
    return res;
}


double kernel(double val, double max, long gamma)
{
    double base = 1 - val / max;

    for(long i = 0; i < gamma; i ++)
        base *= base;
        
    return base;
}

double kernel(Point &point_a, Point &point_b, double max, long gamma)
{
    return kernel(euclidean_distance_sqr(point_a, point_b), max, gamma);
}


std::vector<Point> meanshift(std::vector<Point> &dusts, std::vector<Point> &points, double max, long zeta, long gamma)
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

            ker = kernel(dusts[i], points[j], max, gamma);

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


std::vector<Point> minidx(std::vector<Point> &dusts, std::vector<Point> &points, double max, long zeta, long gamma)
{
    double sum;
    long num_of_points = points.size(), num_of_dusts = dusts.size();
    std::vector<Point> res;

    res.resize(num_of_points);

    for(long i = 0; i < num_of_points; i ++){

        res[i].resize(num_of_dusts);
        sum = 0;

        for(long j = 0; j < num_of_dusts; j ++){
            res[i][j] = kernel(dusts[j], points[i], max, gamma);
            sum += res[i][j];
        }

        sum = goldschmidt(sum, num_of_dusts, zeta);

        for(long j = 0; j < num_of_dusts; j ++) 
            res[i][j] *= sum;        
    }
    
    return res;
}


std::vector<double> nbhd(std::vector<Point> &dusts, double max, long gamma)
{
    long num_of_dusts = dusts.size();
    std::vector<double> res;

    res.resize(num_of_dusts);

    for(long i = 0; i < num_of_dusts; i ++){
        res[i] = 0;
        for(long j = 0; j < num_of_dusts; j ++){
            res[i] += kernel(dusts[i], dusts[j], max, gamma);
        }
    }

    return res;
}


std::vector<Point> index_numbering(std::vector<Point> &dusts, std::vector<Point> &points, double max, long zeta_mi, long gamma_nn, long gamma_mi)
{
    std::vector<Point> mi, res;
    std::vector<double> nn;
    long num_of_points = points.size(), num_of_dusts = dusts.size();

    mi = minidx(dusts, points, max, zeta_mi, gamma_mi);
    nn = nbhd(dusts, max, gamma_nn);

    for(long i = 0; i < num_of_points; i ++) {

        for(long j = 0; j < num_of_dusts; j ++) {

            mi[i][j] *= nn[j];

        }
    }
    
    return mi;
}


std::vector<Point> cluster(std::vector<Point> &points, long dim, long num_of_dusts, long num_of_shiftings, long zeta_ms, long zeta_sm, long gamma_ms, long gamma_sm, long gamma_mm)
{
    long num_of_points = (long) points.size();
    double max = (double) dim;
    std::vector<Point> dusts, clus;

    dusts = sample_points(num_of_dusts, points, num_of_points);

    for(long i = 0; i < num_of_shiftings; i ++)
        dusts = meanshift(dusts, points, max, zeta_ms, gamma_ms);

    clus = index_numbering(dusts, points, max, zeta_sm, gamma_sm, gamma_mm);

    return clus;
}


void print_plain_result(std::vector<Point> &clus, std::vector<Point> &points)
{
    long num_of_points = clus.size(), num_of_dusts = clus[0].size(), dim = points[0].size();
    
    printf("Clustering results:\n");
    for(long i = 0; i < num_of_points; i ++){

        printf("%5ld(", i);
        for(long j = 0; j < dim; j ++){
            if(j)
                printf(" %f", points[i][j]);
            else
                printf("%f", points[i][j]);
        }
        printf(") ");

        for(long j = 0; j < num_of_dusts; j ++){

            if(clus[i][j] > 1.1){
                printf("XXXXXXXX ");
            } else if(clus[i][j] > 0.6) {
                printf("++++++++ ");
            } else if (clus[i][j] < 0.4) {
                printf("-------- ");
            } else {
                printf("%f ", clus[i][j]);
            }

        }
        printf("\n");
    }
}

