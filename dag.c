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

int main() {
    return 0;
}