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

struct Node *nodeArr = NULL;
int nodeArrLen = 0;

void buildNodes(int len);
void buildGJE(int topLevel);
void buildLU(int topLevel);
void printDAG();