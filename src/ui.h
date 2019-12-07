#pragma once
#include "heaansrc.h"
#include <stdio.h>
#include "io.h"
#include "factory.h"
#include "plain.h"
#include "typedefs.h"
#include "cipher.h"
#include "parameter.h"
#include <string>
#include "common.h"

void ui_gen_points(std::string path, long seed);


void ui_gen_scheme(std::string path_sk, std::string path_pk, long seed, Ring &ring);


void ui_gen_cvecs(std::string path_enc, std::string path_dat, std::string path_pk, long logp, long logq, long seed, Ring &ring);


void ui_cluster_plain(std::string path_dat, long num_of_dusts, long num_of_shiftings, long zeta_ms, long zeta_mi, long gamma_ms, long gamma_nn, long gamma_mi, long seed);

void ui_gaussian_cluster_plain(std::string path_dat, long num_of_dusts, long num_of_shiftings, long zeta_ms, long zeta_mi, double sigma_ms, long gamma_nn, long gamma_mi, long seed);
void ui_fk_cluster_plain(std::string path_dat, long num_of_shiftings, long zeta_ms, long gamma_ms, long seed);


void ui_cluster_cipher(std::string path_enc, std::string path_pk, std::string path_sk, long num_of_dusts, long num_of_shiftings, long zeta1, long zeta3, long gamma1, long gamma2, long gamma3, long logc, long num_of_threads, long seed, Ring &ring);
void ui_cluster_cipher_verbose(std::string path_enc, std::string path_pk, std::string path_sk, long num_of_dusts, long num_of_shiftings, long zeta1, long zeta3, long gamma1, long gamma2, long gamma3, long logc, long num_of_threads, long seed, Ring &ring);
