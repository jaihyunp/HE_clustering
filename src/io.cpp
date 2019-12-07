#include "io.h"


std::vector<std::vector<double> > read_points(string path)
{
    const char *filename = path.c_str();
    std::vector<std::vector<double> > points;
    FILE *fp = fopen(filename, "r");
    char line[1000], temp[1000];
    while (fgets(line, sizeof(line), fp) != NULL) {
        double val;
        std::vector<double> point;
        int flag = 0, i, j;
        while (flag < (int) strlen(line)) {
            i = 0;
            while((line[flag + i] != ',') && (line[flag + i] != '\0'))
                i ++;
            for(j = 0; j < i; j ++)
                temp[j] = line[flag + j];
            temp[j] = 0;
            val = atof(temp);
            point.push_back(val);
            flag += i + 1;
        }
        points.push_back(point);
    }
    fclose(fp);
 
    return points;
}


void print_points(std::vector<std::vector<double> > points)
{
    for(int i = 0; i < (int) points.size(); i ++){
        printf("%d: ", i);
        for(int dim = 0; dim < (int) points[i].size(); dim ++) {
            printf("%f ", points[i][dim]);
        }
        printf("\n");
    }
}


void write_points(std::vector<std::vector<double> > points, std::string path)
{
    const char *filename = path.c_str();
    FILE *fp = fopen(filename, "w");
    for(int i = 0; i < (int) points.size(); i ++) {    
        for(int j = 0; j < (int) points[i].size(); j ++){
            if(j == ((int) points[i].size() - 1))
                fprintf(fp, "%f\n", points[i][j]);
            else   
                fprintf(fp, "%f,", points[i][j]);
        }
    }
    fclose(fp);
}


void write_secretKey(SecretKey &sk, std::string path)
{
    long sig;
	unsigned char* bytes = new unsigned char[1];
	fstream fout;
    long n = N;

	fout.open(path, ios::binary|ios::out);

	fout.write(reinterpret_cast<char*>(&n), sizeof(long));

    for(int i = 0; i < N; i ++) {
        sig = sign(sk.sx[i]);
        
        if(sig == 1) {
            bytes[0] = 0x0F;
        } else if (sig == -1) {
            bytes[0] = 0xFF;
        } else {
            bytes[0] = 0x00;
        }

        fout.write(reinterpret_cast<char*>(bytes), 1);
    }
    
	fout.close();
}


int read_secretKey(SecretKey &sk, std::string path)
{
    int stat = 1;
    long n;
	fstream fin;
	unsigned char *bytes = new unsigned char[1];

	fin.open(path, ios::binary|ios::in);
	fin.read(reinterpret_cast<char*>(&n), sizeof(long));
	
	if(n != N){
	
	    stat = 0;
	    
	} else {
	    
	    for(int i = 0; i < N; i ++){
    		fin.read(reinterpret_cast<char*>(bytes), 1);
    		if(bytes[0] == 0x0F) {
                sk.sx[i] = ZZ(1);
    		} else if (bytes[0] == 0xFF) {
                sk.sx[i] = ZZ(-1);
    		} else{
                sk.sx[i] = ZZ(0);
    		}
	    }
	    
	}
	fin.close();
    return stat;
}


void write_ciphertext(Ciphertext& cipher, std::string path)
{
	fstream fout;
	fout.open(path, ios::binary|ios::out);
	long n = cipher.n;
	long logp = cipher.logp;
	long logq = cipher.logq;
	fout.write(reinterpret_cast<char*>(&n), sizeof(long));
	fout.write(reinterpret_cast<char*>(&logp), sizeof(long));
	fout.write(reinterpret_cast<char*>(&logq), sizeof(long));

	long np = ceil(((double)logq + 1)/8);
	ZZ q = conv<ZZ>(1) << logq;
	unsigned char* bytes = new unsigned char[np];
	for (long i = 0; i < N; ++i) {
		cipher.ax[i] %= q;
		BytesFromZZ(bytes, cipher.ax[i], np);
		fout.write(reinterpret_cast<char*>(bytes), np);
	}
	for (long i = 0; i < N; ++i) {
		cipher.bx[i] %= q;
		BytesFromZZ(bytes, cipher.bx[i], np);
		fout.write(reinterpret_cast<char*>(bytes), np);
	}
	fout.close();
}


void read_ciphertext(Ciphertext &cipher, std::string path)
{
	long n, logp, logq;
	fstream fin;

	fin.open(path, ios::binary|ios::in);

	fin.read(reinterpret_cast<char*>(&n), sizeof(long));
	fin.read(reinterpret_cast<char*>(&logp), sizeof(long));
	fin.read(reinterpret_cast<char*>(&logq), sizeof(long));
    cipher.logp = logp;
    cipher.logq = logq;
    cipher.n = n;

	long np = ceil(((double)logq + 1)/8);
	unsigned char* bytes = new unsigned char[np];
    
	for (long i = 0; i < N; ++i) {
		fin.read(reinterpret_cast<char*>(bytes), np);
		ZZFromBytes(cipher.ax[i], bytes, np);
	}

	for (long i = 0; i < N; ++i) {
		fin.read(reinterpret_cast<char*>(bytes), np);
		ZZFromBytes(cipher.bx[i], bytes, np);
	}

	fin.close();
}


void write_cvecs(CVec &vecs, long num_of_vecs, std::string path)
{
    for(long i = 0; i < num_of_vecs; i ++){
        write_ciphertext(vecs[i], path + to_string(i) + ".enc");
    }
}

void read_cvecs(CVec &vecs, long num_of_vecs, std::string path)
{
    vecs.resize(num_of_vecs);
    for(long i = 0; i < num_of_vecs; i ++){
        read_ciphertext(vecs[i], path + to_string(i) + ".enc");
    }
}


void write_key(Key* key, string path)
{
    FILE *fout = fopen(path.c_str(), "w");
    if(fout == NULL){
        printf("Invalid path: %s\n", path.c_str());
        return;
    }
    fwrite(key->rax, sizeof(uint64_t), Nnprimes, fout);
    fwrite(key->rbx, sizeof(uint64_t), Nnprimes, fout);
    fclose(fout);
}


int read_key(Key *key, string path)
{
    int stat = 1;
    FILE *fin = fopen(path.c_str(), "r");

    if(fin == NULL){
        printf("Invalid path: %s\n", path.c_str());
        return 0;
    }

    if(fread(key->rax, sizeof(uint64_t), Nnprimes, fin) < Nnprimes){
        stat = 0;
    }
    if(fread(key->rbx, sizeof(uint64_t), Nnprimes, fin) < Nnprimes){
        stat = 0;
    }
    fclose(fin);
    
    return stat;
}


void read_params(std::vector<long> &params, long num_of_params, std::string path)
{
    FILE *fp = fopen(path.c_str(), "r");
    char line[500];
  
    params.resize(num_of_params); 
    for(long i = 0; i < num_of_params; i ++){
        if(fgets(line, sizeof(line), fp) == NULL){
            printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            break;
        }
        params[i] = atol(line);
    }
    fclose(fp);
}

void write_params(std::vector<long> &params, long num_of_params, std::string path)
{
    FILE *fp = fopen(path.c_str(), "w");
    for(long i = 0; i < num_of_params; i ++)
        fprintf(fp, "%ld\n", params[i]);
    fclose(fp);
}


void read_scheme(Scheme &scheme, std::vector<long> rots, std::string path)
{
    scheme.isSerialized = false;
    
    Key *multkey = new Key();
    read_key(multkey, path + "MULTIPLICATION.key");
    scheme.keyMap.insert(pair<long, Key*>(1, multkey));

    Key *enckey = new Key();
    read_key(enckey, path + "ENCRYPTION.key");
    scheme.keyMap.insert(pair<long, Key*>(0, enckey));

    Key *conKey = new Key();
    read_key(conKey, path + "CONJUGATION.key");
    scheme.keyMap.insert(pair<long, Key*>(2, conKey));

    for(long i = 0; i < (long) rots.size(); i ++) {
        Key *rotkey = new Key();
        read_key(rotkey, path + "ROTATION" + to_string(rots[i]) + ".key");
        scheme.leftRotKeyMap.insert(pair<long, Key*>(rots[i], rotkey));
    }
}

void write_scheme(Scheme &scheme, std::string ipath)
{
    if(scheme.isSerialized) {
        printf("Not implemented!\n");
        return;
    }

    std::string path(ipath);
    if(ipath[ipath.length() - 1] != '/')
        path += "/";
        
    std::string rpath = path + "ENCRYPTION.key";
    write_key(scheme.keyMap.at(0), rpath);
    
    rpath = path + "MULTIPLICATION.key";
    write_key(scheme.keyMap.at(1), rpath);
    
    
    for(long idx = 1; idx < N; idx ++) {
		if(scheme.leftRotKeyMap.find(idx) != scheme.leftRotKeyMap.end()) {
            rpath = path + "ROTATION" + to_string(idx) + ".key";
            write_key(scheme.leftRotKeyMap.at(idx), rpath);
		}
	}
}
