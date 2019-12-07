#include "factory.h"


double standard_normal_sampling()
{
    double res = 0;
    long sqrt_num = 8;

    for(long i = 0; i < sqrt_num * sqrt_num; i ++)
        res += ((double) random() / RAND_MAX);
    res = (res / sqrt_num - 0.5 * sqrt_num) * 12;

    return res;
}

void gaussian_sampling(Point &point, Point mean, Point var, long dim)
{
    point.resize(dim);
    for(long i = 0; i < dim; i ++){
        while(1){
            point[i] = mean[i] + var[i] * standard_normal_sampling();
            if((point[i] < 1) && (point[i] > 0))
                break;
        }
    }
}

void shuffle_index(std::vector<long> &deck, long num_of_cards)
{
    long idx, j;
    unsigned char *check = (unsigned char *) calloc(num_of_cards, sizeof(unsigned char));

    deck.resize(num_of_cards);

    for(long i = 0; i < num_of_cards; i ++){

        idx = (num_of_cards - i) * (((double) random() / RAND_MAX)) + 1;

        for(j = 0; j < num_of_cards; j ++){
            if(!check[j])
                idx --;
            if(!idx) {
                check[j] = 1;
                deck[i] = j;
                break;
            }
        }
        if(j == num_of_cards){
            printf("Error on shuffle_index\n");
        }
    }

    free(check);
}

void generate_points(std::string path, std::vector<long> num_of_elements, Points means, Points vars, long num_of_clusters, long dim) //dim, nop
{
    Points cluster, points, answers;
    std::vector<long> deck;
    std::vector<long> i_answers;
    long num_of_points = 0, ctr = 0;

    for(long i = 0; i < num_of_clusters; i ++)
        num_of_points += num_of_elements[i];

    //shuffle_index(deck, num_of_points);

    deck.resize(num_of_points);
    for(long i = 0; i < num_of_points; i ++)
        deck[i] = i;

    points.resize(num_of_points);
    answers.resize(num_of_points);
            
    for(long i = 0; i < num_of_clusters; i ++) {
        for(long j = 0; j < num_of_elements[i]; j ++){ 

            gaussian_sampling(points[deck[ctr]], means[i], vars[i], dim);
            answers[deck[ctr]].push_back(i);
            ctr ++;

        }
    }
    
    write_points(points, path + ".dat");
    write_points(answers, path + ".ans");

    std::vector<long> params;
    params.push_back(dim);
    params.push_back(num_of_points);
    write_params(params, 2, path + ".params");
}


void generate_secretKey(std::string path, long seed, Ring &ring)//logN, logQ
{
    srand(seed);
    SecretKey sk(ring);
    write_secretKey(sk, path + ".key");

    std::vector<long> params;
    params.push_back(logN);
    params.push_back(logQ);
    write_params(params, 2, path + ".params");
}


void add_if_not(std::vector<long> &rots, long idx)
{
    unsigned char check = 0;
    for(long i = 0; i < (long) rots.size(); i ++){
        if(rots[i] == idx){
            check = 1;
            break;
        }
    }
    if(!check)
        rots.push_back(idx);
}

void add_boots_rots(std::vector<long> &rots, long logSlots, long radix)
{
	// parameters		
	long slots = 1 << logSlots;
	long log2r = log2(radix);
	long logrSlots = logSlots / log2r;
	
	// make public keys for left rotations
	for(long i = 0; i < logrSlots; i++) {
		long gap = slots / (radix * pow(radix, i));
		if(i == 0) {
			if(radix > 4) {
				long bs = 1 << (long)floor((log2(radix)) / 2.);
				long gs = (radix) / bs;
				for(long j = 1; j < bs; j++) {
					long idx = j * gap;
					idx %= slots;
					if(idx != 0)
                        add_if_not(rots, idx);
				}
				for(long j = 1; j < gs; j++) {
					long idx = j * gap * bs;
					idx %= slots;
                    if(idx != 0)
                        add_if_not(rots, idx);
				}
			} else {
				for(long j = 1; j < radix; j++) {
					long idx = j * gap;
					idx %= slots;
					if(idx != 0)
                        add_if_not(rots, idx);
				}
			}
		} else {
            add_if_not(rots, Nh - (radix - 1) * gap);
			if(radix > 2) {
				long bs = 1 << (long)floor((log2(radix) + 1) / 2.);
				long gs = (2 * radix) / bs;
				for(long j = 1; j < bs; j++) {
					long idx = j * gap;
					idx %= slots;
					if(idx != 0)
                        add_if_not(rots, idx);
				}
				for(long j = 1; j < gs; j++) {
					long idx = j * gap * bs;
					idx %= slots;
					if(idx != 0)
                        add_if_not(rots, idx);
				}
			} else {
				for(long j = 1; j < 2 * radix - 1; j++) {
					long idx = j * gap;
                    idx %= slots;
                    if(idx != 0)
                        add_if_not(rots, idx);
                }
            }
        }
    }
}


void generate_scheme(std::string path, SecretKey &sk, Ring &ring)
{
    Key *key = new Key();
    ZZ *ax = new ZZ[N];
    ZZ *bx = new ZZ[N];
    ZZ *spow = new ZZ[N];
    std::vector<long> rots;
    long np = ceil((1 + logQQ + logN + 2) / 59.0);

    Scheme scheme(sk, ring, false);
    scheme.addConjKey(sk);

    for(long i = 0; i < logN - 1; i ++){
        add_if_not(rots, 1 << i);
        add_if_not(rots, Nh - (1 << i));
        printf("adding rotkeys..(%ld)\n", rots.size());
    }
    
    for(long slots = 1; slots < 8; slots ++){
        add_boots_rots(rots, slots, 1 << slots);
        printf("adding rotkeys..(%ld)\n", rots.size());
    }
   
    printf("Need %ld rotkeys:\n", rots.size());

    for(long i = 0; i < (long) rots.size(); i ++){
        ring.sampleUniform2(ax, logQQ);
        ring.mult(bx, sk.sx, ax, np, QQ);
        ring.subFromGaussAndEqual(bx, QQ);

        ring.leftRotate(spow, sk.sx, rots[i]);
        ring.leftShiftAndEqual(spow, logQ, QQ);
        ring.addAndEqual(bx, spow, QQ);

        ring.CRT(key->rax, ax, nprimes);
        ring.CRT(key->rbx, bx, nprimes);

        write_key(key, path + "ROTATION" + to_string(rots[i]) + ".key");
    }

    
    write_key(scheme.keyMap.at(0), path + "ENCRYPTION.key");
    write_key(scheme.keyMap.at(1), path + "MULTIPLICATION.key");
    write_key(scheme.keyMap.at(2), path + "CONJUGATION.key");

    delete[] ax;
    delete[] bx;
    delete[] spow;

    std::vector<long> params;
    params.push_back(logN);
    params.push_back(logQ);
    write_params(params, 2, path + ".params");
}


void generate_cvec(std::string path, long dim, long num_of_points, long num_of_slots, long logp, long logq, Points &points, Scheme &scheme)
{
    Ciphertext cipher;
    long udim = 1 << log2up(dim);
    long unop = 1 << log2up(num_of_points);
    long ppc = num_of_slots / udim;
    if(ppc > unop)
        ppc = unop;
    long pvec = (num_of_points - 1) / ppc + 1;

    double *vals = (double *) calloc(udim * ppc, sizeof(double));
    Plaintext plain;

    for(long i = 0; i < pvec; i ++){

        for(long j = 0; j < ppc; j ++){
            if((i * ppc + j) >= num_of_points){
                for(long k = 0; k < dim; k ++)
                    vals[j * udim + k] = 0;
            } else {
                for(long k = 0; k < dim; k ++)
                    vals[j * udim + k] = points[i * ppc + j][k];
            }
        }
        scheme.encode(plain, vals, udim * ppc, logp, logq);
        scheme.encryptMsg(cipher, plain);
        cipher.n = num_of_slots;
        write_ciphertext(cipher, path + to_string(i) + ".enc");
    }

    free(vals);

    std::vector<long> params;
    params.push_back(dim);
    params.push_back(num_of_points);
    params.push_back(num_of_slots);
    params.push_back(logp);
    params.push_back(logq);
    params.push_back(logN);
    params.push_back(logQ);
    write_params(params, 7, path + ".params");

}

void rotkeys(std::vector<long> &rots, long bootstrapping_slots)
{
    rots.resize(0);
    for(long i = 0; i < logN - 1; i ++){
        add_if_not(rots, 1 << i);
        add_if_not(rots, Nh - (1 << i));
    }
    add_boots_rots(rots, bootstrapping_slots, 1 << bootstrapping_slots);
}


