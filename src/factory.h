#pragma once
#include "common.h"
#include <vector>
#include "io.h"
#include "typedefs.h"
#include "heaansrc.h"
#include "parameter.h"

void generate_points(std::string path, std::vector<long> num_of_elements, Points means, Points vars, long num_of_clusters, long dim);

void generate_secretKey(std::string path, long seed, Ring &ring);
void generate_scheme(std::string path, SecretKey &sk, Ring &ring);

void generate_cvec(std::string path, long dim, long num_of_points, long num_of_slots, long logp, long logq, Points &points, Scheme &scheme);

void rotkeys(std::vector<long> &rots, long bootstrapping_slots);
