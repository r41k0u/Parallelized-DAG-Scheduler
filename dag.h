#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

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
};

struct Node *nodeArr = NULL;
int nodeArrLen = 0;

const struct DAGConfig GJE_CONFIG = {40, 100, 4, 30};

void buildNodes(int len);
void buildGJE(int topLevel);
void buildLU(int topLevel);
void printDAG();