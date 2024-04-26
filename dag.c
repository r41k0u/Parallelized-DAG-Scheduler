#include "dag.h"

// Initialize All nodes
void buildNodes(int len) {
    nodeArr = (struct Node *)malloc(len * sizeof(struct Node));
    if (nodeArr == NULL) {
        fprintf(stderr, "Node array malloc failed, %s\n", strerror(errno));
        return;
    }

    for (int i = 0; i < len; i++) {
        nodeArr[i].index = i;
        nodeArr[i].childrenLen = nodeArr[i].depsLen = 0;
        nodeArr[i].deps = nodeArr[i].children = NULL;
    }

    nodeArrLen = len;
}

// Build a GJE adjacency list by populating the dependenices and children
void buildGJE(int topLevel) {
    buildNodes((topLevel * (topLevel + 1)) / 2);
    int prefix;

    for (int i = 0; i < topLevel; i++) {
        // Exit node
        if (i == (topLevel - 1)) break;
        prefix = ((topLevel * i) - ((i * (i - 1)) / 2)); 
        for (int j = 0; j < (topLevel - i); j++) {
            // First node in a layer
            if (j == 0) {
                nodeArr[prefix + j].children = (struct Node **)malloc((topLevel - i - 1) * sizeof(struct Node *));
                if (nodeArr[prefix + j].children == NULL) {
                    fprintf(stderr, "Node array child malloc failed, %s\n", strerror(errno));
                    return;
                }
                nodeArr[prefix + j].childrenLen = (topLevel - i - 1);

                for (int k = 0; k < (topLevel - i - 1); k++) {
                    nodeArr[prefix + (topLevel - i) + k].deps = (struct Node **)malloc(2 * sizeof(struct Node *));
                    if (nodeArr[prefix + (topLevel - i) + k].deps == NULL) {
                        fprintf(stderr, "Node array deps malloc failed, %s\n", strerror(errno));
                        return;
                    }
                    nodeArr[prefix + (topLevel - i) + k].depsLen = 2;

                    nodeArr[prefix + (topLevel - i) + k].deps[0] = nodeArr + prefix + j;
                    nodeArr[prefix + j].children[k] = nodeArr + prefix + (topLevel - i) + k; 
                }
                continue;
            }

            nodeArr[prefix + j].children = (struct Node **)malloc(1 * sizeof(struct Node *));
            if (nodeArr[prefix + j].children == NULL) {
                fprintf(stderr, "Node array child malloc failed, %s\n", strerror(errno));
                return;
            }
            nodeArr[prefix + j].childrenLen = 1;
            nodeArr[prefix + j].children[0] = nodeArr + prefix + (topLevel - i - 1) + j;
            nodeArr[prefix + (topLevel - i - 1) + j].deps[1] = nodeArr + prefix + j;
        }
    }
}

void buildLU(int topLevel) {
    buildNodes((topLevel * (topLevel + 3)) / 2);

    // screw it I will write this later
}

void initOcean() {
    tempSeqGenerator = NULL;
    ocean.whales = (struct Whale *)malloc(GJE_CONFIG.populationSize * sizeof(struct Whale));
    if (ocean.whales == NULL) {
        fprintf(stderr, "Ocean malloc failed, %s\n", strerror(errno));
        return;
    }

    ocean.lenWhales = GJE_CONFIG.populationSize;
    ocean.globBestSeq = NULL;
    ocean.lenGlobBestSeq = 0;
    ocean.globBestObj = INIT_VAL;

    procState = (int *)malloc(GJE_CONFIG.populationSize * sizeof(int));
    if (procState == NULL) {
        fprintf(stderr, "procState malloc failed, %s\n", strerror(errno));
        return;
    }
    procAlloc = (int *)malloc(nodeArrLen * sizeof(int));
    if (procAlloc == NULL) {
        fprintf(stderr, "procAlloc malloc failed, %s\n", strerror(errno));
        return;
    }
    nodeEndTime = (int *)malloc(nodeArrLen * sizeof(int));
    if (nodeEndTime == NULL) {
        fprintf(stderr, "nodeEndTime malloc failed, %s\n", strerror(errno));
        return;
    }
    
    // Seed the random number generator
    srand(time(NULL));

    for (int i = 0; i < ocean.lenWhales; i++) {
        ocean.whales[i].pos = (double *)malloc(nodeArrLen * sizeof(double));
        ocean.whales[i].seq = (int *)malloc(nodeArrLen * sizeof(int));
        if (ocean.whales[i].pos == NULL || ocean.whales[i].seq == NULL) {
            fprintf(stderr, "Whale malloc failed, %s\n", strerror(errno));
            return;   
        }
        ocean.whales[i].lenPos = nodeArrLen;

        // Random position of whale
        for (int j = 0; j < nodeArrLen; j++) {
            ocean.whales[i].pos[j] = ((double)rand() / RAND_MAX) * 20.0 - 10.0;
            ocean.whales[i].seq[j] = j;
        }

        // Sort the sequence array wrt the values of pos array
        tempSeqGenerator = ocean.whales[i].pos;
        qsort(ocean.whales[i].seq, nodeArrLen, sizeof(int), compareDoubles);
        
        // Change the pos and seq array to comply with the toposort
        // toposort for GJE and LU is simply the integral order
        // Just check if all the dependencies of a node appear before it or not sequentially
        // If not, swap with the last one
        topoSortSeq(i);
        calcMakespan(i);
    }

    tempSeqGenerator = NULL;
}

void topoSortSeq(int ind) {
    int *visitedArr = (int *)malloc(nodeArrLen * sizeof(int));
    if (visitedArr == NULL) {
        fprintf(stderr, "Whale vistedArr malloc failed, %s\n", strerror(errno));
        return;  
    }
    int tempSeq;
    double tempPos;
    for (int i = 0; i < nodeArrLen; i++) visitedArr[i] = INIT_NODE;
    int iter = 0, childrenResolved, minChain;
    while (iter < nodeArrLen) {
        if (nodeArr[ocean.whales[ind].seq[iter]].childrenLen == 0) {
            visitedArr[ocean.whales[ind].seq[iter]] = iter;
            iter++;
            continue;
        }
        
        childrenResolved = 0;
        while (!childrenResolved) {
            // This just finds the first occurence of a child
            minChain = INIT_NODE;
            for (int i = 0; i < nodeArr[ocean.whales[ind].seq[iter]].childrenLen; i++)
                minChain = (visitedArr[nodeArr[ocean.whales[ind].seq[iter]].children[i]->index] < minChain) ? visitedArr[nodeArr[ocean.whales[ind].seq[iter]].children[i]->index] : minChain;

            if (minChain < iter)
                childrenResolved = 0;
            else
                childrenResolved = 1;
            if (!childrenResolved) {
                tempSeq = ocean.whales[ind].seq[iter];
                ocean.whales[ind].seq[iter] = ocean.whales[ind].seq[minChain];
                ocean.whales[ind].seq[minChain] = tempSeq;
                tempPos = ocean.whales[ind].pos[ocean.whales[ind].seq[iter]];
                ocean.whales[ind].pos[ocean.whales[ind].seq[iter]] = ocean.whales[ind].pos[ocean.whales[ind].seq[minChain]];
                ocean.whales[ind].pos[ocean.whales[ind].seq[minChain]] = tempPos;
                visitedArr[ocean.whales[ind].seq[minChain]] = minChain;
                visitedArr[ocean.whales[ind].seq[iter]] = INIT_NODE;
            }
        }
        visitedArr[ocean.whales[ind].seq[iter]] = iter;
        iter++;
    }
}

void calcMakespan(int ind) {
    for (int i = 0; i < GJE_CONFIG.numProc; i++)
        procState[i] = 0;
    for (int i = 0; i < nodeArrLen; i++) {
        procAlloc[ocean.whales[ind].seq[i]] = (i % GJE_CONFIG.numProc);
        nodeEndTime[i] = 0;
    }

    int maxEndTime;
    for (int i = 0; i < nodeArrLen; i++) {
        if (nodeArr[ocean.whales[ind].seq[i]].depsLen == 0) {
            procState[procAlloc[ocean.whales[ind].seq[i]]] += GJE_CONFIG.nodeCost;
            nodeEndTime[ocean.whales[ind].seq[i]] = procState[procAlloc[ocean.whales[ind].seq[i]]];
            continue;
        }

        maxEndTime = procState[procAlloc[ocean.whales[ind].seq[i]]];
        for (int j = 0; j < nodeArr[ocean.whales[ind].seq[i]].depsLen; j++) {
            if (procAlloc[ocean.whales[ind].seq[i]] != procAlloc[nodeArr[ocean.whales[ind].seq[i]].deps[j]->index])
                maxEndTime = ((nodeEndTime[nodeArr[ocean.whales[ind].seq[i]].deps[j]->index] + GJE_CONFIG.commCost) > maxEndTime) ? (nodeEndTime[nodeArr[ocean.whales[ind].seq[i]].deps[j]->index] + GJE_CONFIG.commCost) : maxEndTime;
        }

        procState[procAlloc[ocean.whales[ind].seq[i]]] = (maxEndTime + GJE_CONFIG.nodeCost);
        nodeEndTime[ocean.whales[ind].seq[i]] = procState[procAlloc[ocean.whales[ind].seq[i]]];
    }

    ocean.whales[ind].obj = procState[0];
    for (int i = 0; i < GJE_CONFIG.numProc; i++)
        ocean.whales[ind].obj = (procState[i] > ocean.whales[ind].obj) ? procState[i] : ocean.whales[ind].obj;
}

void printDAG() {
    for (int i = 0; i < nodeArrLen; i++) {
        printf("Node %d:\nDeps:\n", nodeArr[i].index);
        for (int j = 0; j < nodeArr[i].depsLen; j++)
            printf("%d ", nodeArr[i].deps[j]->index);
        printf("\nChildren:\n");
        for (int j = 0; j < nodeArr[i].childrenLen; j++)
            printf("%d ", nodeArr[i].children[j]->index);
        printf("\n\n\n");
    }
}

void printWhales() {
    printf("Whales:\n\n");
    for (int i = 0; i < ocean.lenWhales; i++) {
        printf("Whale %d:\nPos Array:\n", i);
        for (int j = 0; j < ocean.whales[i].lenPos; j++) {
            printf("%f ", ocean.whales[i].pos[j]);
        }
        printf("\nSeq Array:\n");
        for (int j = 0; j < ocean.whales[i].lenPos; j++) {
            printf("%d ", ocean.whales[i].seq[j]);
        }
        printf("\nObj: %d\n\n\n", ocean.whales[i].obj);
    }
}

// Comparison function for qsort
int compareDoubles(const void *a, const void *b) {
    double diff = tempSeqGenerator[*(int*)a] - tempSeqGenerator[*(int*)b];
    return (diff > 0) - (diff < 0);
}

int main() {
    buildGJE(5);
    initOcean();
    printWhales();
    return 0;
}