#pragma once
#include "common.h"
#include <vector>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "typedefs.h"

std::vector<Point> sample_points(int num_of_dusts, std::vector<Point> &points, int num_of_points);
double euclidean_distance_sqr(Point &point_a, Point &point_b);
double goldschmidt(double val, double given_max, long zeta);


std::vector<Point> index_numbering(std::vector<Point> &dusts, std::vector<Point> &points, double max, long zeta_mi, long gamma_nn, long gamma_mi);


double kernel(Point &point_a, Point &point_b, double max, long gamma);
std::vector<Point> meanshift(std::vector<Point> &dusts, std::vector<Point> &points, double max, long zeta, long gamma);
std::vector<Point> cluster(std::vector<Point> &points, long dim, long num_of_dusts, long num_of_shiftings, long zeta_ms, long zeta_sm, long gamma_ms, long gamma_sm, long gamma_mm);


double gaussian_kernel(Point &point_a, Point &point_b, double max, double sigma);
std::vector<Point> gaussian_meanshift(std::vector<Point> &dusts, std::vector<Point> &points, double max, long zeta, double sigma);
std::vector<Point> gaussian_cluster(std::vector<Point> &points, long dim, long num_of_dusts, long num_of_shiftings, long zeta_ms, long zeta_sm, double sigma_ms, long gamma_sm, long gamma_mm);


std::vector<long> freedman_cluster(std::vector<Point> &points, long dim, long num_of_dusts, long num_of_shiftings, long zeta_ms, long gamma_ms);


void print_plain_result(std::vector<Point> &clus, std::vector<Point> &points);
