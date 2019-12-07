#pragma once
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include "typedefs.h"
#include "heaansrc.h"
#include <complex>
#include <dirent.h>

Points read_points(std::string path);
void write_points(Points points, std::string path);
void print_points(Points points);

void write_secretKey(SecretKey &sk, std::string path);
int read_secretKey(SecretKey &sk, std::string path);

void write_key(Key *key, std::string path);
int read_key(Key *key, std::string path);

void write_ciphertext(Ciphertext& cipher, std::string path);
void read_ciphertext(Ciphertext &cipher, std::string path);

void write_cvecs(CVec &vecs, long num_of_vecs, std::string path);
void read_cvecs(CVec &vecs, long num_of_vecs, std::string path);

void read_scheme(Scheme &scheme, std::vector<long> rots, std::string path);
void write_scheme(Scheme &scheme, std::string path);

void read_params(std::vector<long> &params, long num_of_params, std::string path);
void write_params(std::vector<long> &params, long num_of_params, std::string path);
