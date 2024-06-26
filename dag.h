#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define INIT_VAL 1e9
#define INIT_NODE (nodeArrLen + 100)
#define M_PI 3.14159265358979323846

struct Node {
    struct Node **deps;
    int depsLen;
    int index;
    struct Node **children;
    int childrenLen;
};

struct DAGConfig {
    int nodeCost;
    int commCost;
    int numProc;
    int populationSize;
    int iterarions;
};

struct Whale {
    double *pos;
    int lenPos;

    int *seq;
    int obj;
};

struct Ocean {
    struct Whale *whales;
    int lenWhales;

    int *globBestSeq;
    int lenGlobBestSeq;
    int globBestObj;
};

struct Node *nodeArr = NULL;
int nodeArrLen = 0;

struct Ocean ocean;

const struct DAGConfig GJE_CONFIG = {40, 100, 4, 72, 50};

void buildNodes(int len);
void buildGJE(int topLevel);
void buildLU(int topLevel);
void initOcean();
void topoSortSeq(int ind);

int *procState, *procAlloc, *nodeEndTime;
void calcMakespan(int ind);

void runOptimization();

void printDAG();
void printWhales();

double *tempSeqGenerator;
int compareDoubles(const void *a, const void *b);