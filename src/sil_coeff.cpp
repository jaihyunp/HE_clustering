#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>

using namespace std;

int read_data(vector<vector<double> > &points, vector<int> &idcs, const char *path)
{
    FILE *fp = fopen(path, "r");
    char line[500], temp[500];
    int num_of_clus = 0;
    while(fgets(line, sizeof(line), fp) != NULL) {
        vector<double> point;
        int flag = 0, i, j;
        while(flag < (int) strlen(line)) {
            i = 0;
            if(flag)
                point.push_back(atof(temp));
            while((line[flag + i] != ',') && (line[flag + i] != '\0'))
                i ++;
            for(j = 0; j < i; j ++)
                temp[j] = line[flag + j];
            temp[j] = 0;
            flag += i + 1;
        }
        points.push_back(point);
        if(atoi(temp) > num_of_clus)
            num_of_clus = atoi(temp);
        idcs.push_back(atoi(temp));
    }
    fclose(fp);
    return num_of_clus;
}

void classify(vector<vector<vector<double> > > &clus, vector<vector<double> > &points, vector<int> &idcs, int num_of_clus)
{
    int num_of_points = points.size(), dim = points[0].size();
    clus.resize(num_of_clus);
    for(int i = 0; i < num_of_clus; i ++){
        clus[i].resize(0);
    }
    
    for(int i = 0; i < num_of_points; i ++){
        vector<double> point;
        point.resize(dim);
        for(int j = 0; j < dim; j ++)
            point[j] = points[i][j];
        clus[idcs[i] - 1].push_back(point);
    }
}

void print_data(vector<vector<vector<double> > > &clus)
{
    int num_of_clus = clus.size(), dim = clus[0][0].size();

    for(int i = 0; i < num_of_clus; i ++){
        printf("[%d]\n", i);
        for(int j = 0; j < (int)clus[i].size(); j ++){
            printf("%d: ", j);
            for(int k = 0; k < dim; k ++){
                printf("%f ", clus[i][j][k]);
            }
            printf("\n");
        }
        printf("=====%d\n", (int) clus[i].size());
    }
}

double dist(vector<double> &p1, vector<double> &p2)
{
    double d = 0;
    int dim = p1.size();

    for(int i = 0; i < dim; i ++)
        d += (p1[i] - p2[i]) * (p1[i] - p2[i]);

    return sqrt(d);
}

double calc_a(vector<vector<double> > &clu, vector<double> point)
{
    double sum = 0;

    for(int i = 0; i < (int)clu.size(); i ++)
        sum += dist(clu[i], point);

    if(clu.size() > 1)
        sum /= (clu.size() - 1);

    return sum;
}

double calc_b(vector<vector<double> > &clu, vector<double> point)
{
    double sum = 0;

    for(int i = 0; i < (int)clu.size(); i ++)
        sum += dist(clu[i], point);

    return sum / clu.size();
}

double silhouette_coefficient_single(vector<vector<vector<double> > > &clus, vector<double> point, int idx)
{
    double a = 0, b = 999, sil, tmp;

    for(int i = 0; i < (int)clus.size(); i ++){
        
        if(i == idx) {
            a = calc_a(clus[i], point);
        } else {
            tmp = calc_b(clus[i], point);
            if(tmp < b)
                b = tmp;
        }

    }

    if(clus[idx].size() == 1) {
        sil = 0;
    } else {
        if(a < b)
            sil = 1 - a / b;
        else if (a > b)
            sil = b / a - 1;
        else //if (a == b)
            sil = 0;
    }
    return sil;
}


double silhouette_coefficient(vector<vector<vector<double> > > &clus)
{
    int num_of_clus = clus.size(), num_of_points = 0;
    double sum = 0, val;

    for(int i = 0; i < num_of_clus; i ++){
        printf("%d th cluster:\n", i);
        for(int j = 0; j < (int)clus[i].size(); j ++){
            val = silhouette_coefficient_single(clus, clus[i][j], i);
            sum += val;
            num_of_points ++;
            printf("- %d th point: %f\n", j, val);
        }
    }
    return sum / num_of_points;
}

int main(int argc, char **argv)
{
    vector<vector<double> > points;
    vector<int> idcs;
    vector<vector<vector<double> > > clus;

    if(argc < 2)
        return 0;
    int num_of_clus = read_data(points, idcs, argv[1]);
    printf("%d clusters\n", num_of_clus);
    classify(clus, points, idcs, num_of_clus);
    print_data(clus);
    double sil = silhouette_coefficient(clus);
    printf("Total sil val = %f\n", sil);
    return 0;
}
