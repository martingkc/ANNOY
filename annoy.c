#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


#define K 16
#define LEN_UUID 37


typedef struct ScorePair {
    float key;
    float *value;
    char uuid[LEN_UUID];
} ScorePair;

typedef struct VectorPair {
    char uuid[LEN_UUID];
    float *value;
} VectorPair;

typedef struct Node {
    VectorPair **data; // list of VectorPair
    int size;
    struct Node *left;
    struct Node *right;
    struct Node *parent;
    float *normal; // normal vector for the sector
    float indexedMedian; // decision val
} Node;

typedef struct List {
    VectorPair **data;
    int size;
    struct List *next;
} List;

void generateUUID(char uuid[LEN_UUID]) {
    unsigned char b[16];

    for (int i = 0; i < 16; i++) {
        b[i] = (unsigned char) (rand() & 0xFF);
    }

    b[6] = (b[6] & 0x0F) | 0x40;

    b[8] = (b[8] & 0x3F) | 0x80;

    sprintf(uuid,
            "%02x%02x%02x%02x-"
            "%02x%02x-"
            "%02x%02x-"
            "%02x%02x-"
            "%02x%02x%02x%02x%02x%02x",
            b[0], b[1], b[2], b[3],
            b[4], b[5],
            b[6], b[7],
            b[8], b[9],
            b[10], b[11], b[12], b[13], b[14], b[15]);
}

int compPairs(const void *p1, const void *p2) {
    const ScorePair *x = p1;
    const ScorePair *y = p2;
    if (x->key < y->key)
        return 1;
    if (x->key > y->key)
        return -1;
    return 0;
}

int comp(const void *elem1, const void *elem2) {
    float f = *((float *) elem1);
    float s = *((float *) elem2);
    if (f > s)
        return 1;
    if (f < s)
        return -1;
    return 0;
}

void deleteList(List *head) {
    List *tmp;
    while (head) {
        tmp = head->next;
        free(head);
        head = tmp;
    }
}

Node *createNode(VectorPair **data, int size) {
    Node *newNode = (Node *) malloc(sizeof(Node));
    if (newNode == NULL)
        return NULL;
    newNode->data = data;
    newNode->size = size;
    newNode->left = NULL;
    newNode->right = NULL;
    newNode->parent = NULL;
    return newNode;
}

VectorPair *createVectorPair(float *data) {
    VectorPair *newNode = (VectorPair *) malloc(sizeof(VectorPair));
    if (newNode == NULL)
        return NULL;
    generateUUID(newNode->uuid);
    newNode->value = data;
    return newNode;
}

List *createListNode(VectorPair **data, int size) {
    List *newList = (List *) malloc(sizeof(List));
    if (newList == NULL)
        return NULL;
    newList->data = data;
    newList->size = size;
    newList->next = NULL;
    return newList;
}

void pushToList(List *head, VectorPair **data, int size) {
    List *tmp;
    List *newItem;
    if (head->data == NULL || head->size == -1) {
        head->data = data;
        head->size = size;
    } else {
        tmp = head;
        newItem = createListNode(data, size);
        while (tmp->next != NULL) {
            tmp = tmp->next;
        }

        tmp->next = newItem;
        newItem->next = NULL;
    }

    return;
}

float dotProd(float *a, float *b, int size) {
    float sum;
    sum = 0.0f;

    for (int i = 0; i < size; i++) {
        sum += a[i] * b[i];
    }

    return sum;
}

float computeProjection(float *normal, float *vector, int sizeVectors) {
    float nDotv;
    float nDotn;
    float factor;
    float tempRes;
    float result;

    result = 0;

    nDotv = dotProd(normal, vector, sizeVectors);
    nDotn = dotProd(normal, normal, sizeVectors);

    factor = nDotv / nDotn;

    for (int i = 0; i < sizeVectors; i++) {
        tempRes = normal[i] * factor;
        result += tempRes * tempRes;
    }

    return sqrt(result);
}

float cosineSimilarity(float *a, float *b, int dataSize) {
    float aDotb;
    float aNormal;
    float bNormal;

    aNormal = 0;
    bNormal = 0;

    for (int i = 0; i < dataSize; i++) {
        aNormal += a[i] * a[i];
        bNormal += b[i] * b[i];
    }

    aNormal = (float) sqrt(aNormal);
    bNormal = (float) sqrt(bNormal);

    aDotb = dotProd(a, b, dataSize);

    return aDotb / (aNormal + bNormal);
}

float *determineNormal(float *a, float *b, int size) {
    float *normal;

    normal = malloc(size * sizeof(float));

    if (normal == NULL) {
        return NULL;
    }

    for (int i = 0; i < size; i++) {
        normal[i] = a[i] - b[i];
    }

    return normal;
}

float calculateMedian(float *projections, int size) {
    float median;
    median = 0.0f;
    qsort(projections, size, sizeof(float), comp);

    if (size % 2 == 0) {
        median += projections[size / 2];
        median += projections[size / 2 - 1];
        median /= 2;
    } else {
        median += projections[(size + 1) / 2 - 1];
    }

    return median;
}

void addSplit(Node *self, int dataSize) {
    float *a;
    float *b;
    float *normal;
    float *projections;
    VectorPair **dataL;
    VectorPair **dataR;
    VectorPair **data;
    float median;
    int idxA;
    int idxB;
    int size;
    int sizeL;
    int sizeR;
    int cnt;
    int idxL, idxR;
    idxL = 0;
    idxR = 0;
    size = self->size;
    data = self->data;
    cnt = 0;


    projections = malloc(sizeof(float) * size);

    if (projections == NULL) {
        return;
    }

    idxA = rand() % size;
    idxB = rand() % size;

    if (idxA == idxB) {
        idxB = rand() % size; // if we get 2 times the same values o zmn sictik
    }

    a = data[idxA]->value;
    b = data[idxB]->value;

    normal = determineNormal(a, b, dataSize);
    if (normal == NULL) return;

    for (int i = 0; i < size; i++) {
        projections[i] = computeProjection(normal, data[i]->value, dataSize);
    }

    median = calculateMedian(projections, size);


    for (int i = 0; i < size; i++) {
        if (projections[i] <= median) {
            cnt += 1;
        }
    }

    sizeL = cnt;
    sizeR = size - cnt;

    if (sizeL == 0 || sizeR == 0) {
        free(dataL);
        free(dataR);
        free(projections);
        return;
    }


    dataL = malloc(sizeL * sizeof(*dataL));
    dataR = malloc(sizeR * sizeof(*dataR));

    if (dataL == NULL || dataR == NULL) {
        free(dataL);
        free(dataR);
        free(projections);
        return;
    }

    for (int i = 0; i < size; i++) {
        if (projections[i] <= median) {
            dataL[idxL] = data[i];
            idxL += 1;
        } else {
            dataR[idxR] = data[i];
            idxR += 1;
        }
    }

    self->normal = normal;
    self->indexedMedian = median;

    self->left = createNode(dataL, sizeL);
    self->right = createNode(dataR, sizeR);

    if (self->left == NULL || self->right == NULL) {
        return;
    }

    self->right->parent = self;
    self->left->parent = self->right->parent;

    // If the size of a node is greater than K it gets split so it means that they arent leafs
    if (sizeL > K + 1) {
        addSplit(self->left, dataSize);
    }
    if (sizeR > K + 1) {
        addSplit(self->right, dataSize);
    }

    // Nodes that aren't leafs don't get the size or data fields filled.

    free(self->data);
    free(projections);

    self->data = NULL;
    self->size = -1;
}

void addVector(Node *self, float *vector, int dataSize) {
    VectorPair **tmpData;
    float projection;

    // this would be better with a check on the leafs
    if (self->data == NULL) {
        projection = computeProjection(self->normal, vector, dataSize);
        if (projection > self->indexedMedian) {
            addVector(self->right, vector, dataSize);
        } else {
            addVector(self->left, vector, dataSize);
        }
    } else {
        tmpData = malloc((self->size + 1) * sizeof(VectorPair *));
        if (tmpData == NULL)
            return;
        for (int i = 0; i < (self->size + 1); i++) {
            if (i != (self->size)) {
                tmpData[i] = self->data[i];
            } else {
                tmpData[i] = createVectorPair(vector);
            }
        }

        free(self->data);
        self->data = tmpData;

        self->size = self->size + 1;
        if (self->size > K) {
            addSplit(self, dataSize);
        }
    }
}

/*
Add serialization and deserialization for the tree.
Add search
*/

/*
go to the leaf that you belong then traverse the adjacent leafs until you reach topK dataPoints

possible solutions:
	- go to the adjacent (topK/K) * 2 nodes then gather them to calculate cosine or eucledian similarity.
	- for now i've included parent node pointer in the node def. I should try also adding a trasverse stack to see if it performs better
*/

Node *recursiveNodeSearch(Node *self, float *vector, int dataSize) {
    if (!self) return NULL;

    /* internal node → decide which branch to follow */
    if (self->data == NULL) {
        float proj = computeProjection(self->normal, vector, dataSize);
        Node *next = (proj > self->indexedMedian) ? self->right : self->left;

        /* If that branch is missing, treat *this* internal node as the leaf. */
        return next ? recursiveNodeSearch(next, vector, dataSize) : self;
    }

    /* leaf node */
    return self;
}

void collectLevel(Node *self, List *head) {
    if (self == NULL) /* <- NEW: safe-guard against NULL child   */
        return;

    if (self->data == NULL) {
        /* internal node – recurse on children     */
        collectLevel(self->right, head);
        collectLevel(self->left, head);
    } else {
        /* leaf – add its vectors to the list      */
        pushToList(head, self->data, self->size);
    }
}

VectorPair **getVectorsList(List *head, int dataSize, int *retSize) {
    List *tmp;
    VectorPair **vectors;
    int totSize;
    int pos;
    int i;

    totSize = 0;
    tmp = head;
    while (tmp) {
        totSize += tmp->size;
        tmp = tmp->next;
    }

    vectors = malloc(totSize * sizeof(*vectors));
    if (vectors == NULL) {
        return NULL;
    }

    tmp = head;
    pos = 0;
    while (tmp) {
        for (i = 0; i < tmp->size; i++) {
            vectors[pos++] = tmp->data[i];
        }
        tmp = tmp->next;
    }

    *retSize = totSize;

    return vectors;
}

ScorePair *searchTopK(Node *self, float *vector, int topK, int dataSize, int *size) {
    /**
     * Node self -> head of the ANN tree
     * float* vector -> vector to do similarity search
     * int topK -> num of results to fetch
     * int dataSize -> size of vectors
     * int *size -> variable where the size of the return vector gets written tbh its unnecessary since topK is set but still havent added a delimiter and its not a given that there will be topK results maybe we can add dummy data or NULL vals to the list to complete it.
     *
     *
     * TODOS:
     * 	1 - Check the mem allocation I think there might be a mem leak somewhere
     *  2 - Add an enum as an input so that you can choose the similarity function to use.
     */

    List *head; // list of vectors holding the results
    Node *startNode; // starting node from which to start trasversing
    VectorPair **vectors;
    float *similarity;
    ScorePair *arr; // a map style return type key = simScore // value = vector returns ordered
    ScorePair *topKResults;
    int nodesToFind;
    head = createListNode((VectorPair **) NULL, -1);

    startNode = recursiveNodeSearch(self, vector, dataSize);
    if (startNode == NULL) {
        *size = 0;
        return NULL;
    }
    nodesToFind = (int) sqrt((topK * 2 / (float) K) * 1.25f) + 2;

    if (startNode->parent == NULL)
        return NULL;

    do {
        startNode = startNode->parent;
        nodesToFind -= 1;
    } while (nodesToFind >= 0 && startNode->parent != NULL);

    collectLevel(startNode, head);


    vectors = getVectorsList(head, dataSize, size);
    deleteList(head);
    similarity = malloc(sizeof(float) * (*size));
    arr = malloc(sizeof(ScorePair) * (*size));

    if (similarity == NULL || arr == NULL)
        return NULL;

    for (int i = 0; i < (*size); i++) {
        similarity[i] = cosineSimilarity(vectors[i]->value, vector, dataSize);
    }

    for (int i = 0; i < (*size); i++) {
        arr[i].key = similarity[i];
        arr[i].value = vectors[i]->value;
        strcpy(arr[i].uuid, vectors[i]->uuid);
    }
    free(vectors);
    free(similarity);

    qsort(arr, *(size), sizeof *arr, compPairs);
    if (*size > topK) {
        topKResults = malloc(sizeof(*topKResults) * topK);
        if (topKResults == NULL) { return arr; }
        *size = topK;
        for (int i = 0; i < topK; i++) {
            topKResults[i] = arr[i];
        }
        free(arr);
        arr = topKResults;
    }
    return arr;
}

void freeTree(Node *n, int dataSize) {
    if (!n)
        return;
    if (n->data) {
        for (int i = 0; i < n->size; i++) {
            free(n->data[i]); // free each VectorPair
        }
        free(n->data);
    } else {
        free(n->normal);
    }
    freeTree(n->left, dataSize);
    freeTree(n->right, dataSize);
    // only for internal nodes
    free(n);
}

void serializeTree(Node *root, FILE *f, int dataSize) {
    if (!root) return;

    if (root->data) {
        // Leaf node
        int marker = 0;
        fwrite(&marker, sizeof(int), 1, f);
        fwrite(&root->size, sizeof(int), 1, f);
        for (int i = 0; i < root->size; i++) {
            fwrite(root->data[i]->uuid, sizeof(char), LEN_UUID, f);
            fwrite(root->data[i]->value, sizeof(float), dataSize, f);
        }
    } else {
        // Internal node
        int marker = 1;
        fwrite(&marker, sizeof(int), 1, f);
        fwrite(&root->indexedMedian, sizeof(float), 1, f);
        fwrite(root->normal, sizeof(float), dataSize, f);
        serializeTree(root->left, f, dataSize);
        serializeTree(root->right, f, dataSize);
    }
}

Node *deserializeTree(FILE *f, int dataSize) {
    int marker;
    if (fread(&marker, sizeof(int), 1, f) != 1) return NULL;

    if (marker == 0) {
        // Leaf node
        int size;
        if (fread(&size, sizeof(int), 1, f) != 1) return NULL;

        VectorPair **data = malloc(size * sizeof(VectorPair *));
        if (!data) return NULL;

        for (int i = 0; i < size; i++) {
            VectorPair *vp = malloc(sizeof(VectorPair));
            if (!vp) {
                while (i-- > 0) free(data[i]);
                free(data);
                return NULL;
            }

            if (fread(vp->uuid, sizeof(char), LEN_UUID, f) != LEN_UUID ||
                !(vp->value = malloc(dataSize * sizeof(float)))) {
                free(vp);
                while (i-- > 0) free(data[i]);
                free(data);
                return NULL;
            }

            if (fread(vp->value, sizeof(float), dataSize, f) != dataSize) {
                free(vp->value);
                free(vp);
                while (i-- > 0) free(data[i]);
                free(data);
                return NULL;
            }
            data[i] = vp;
        }

        Node *node = createNode(data, size);
        if (!node) {
            for (int i = 0; i < size; i++) free(data[i]);
            free(data);
        }
        return node;
    } else if (marker == 1) {
        // Internal node
        Node *node = malloc(sizeof(Node));
        if (!node) return NULL;

        if (fread(&node->indexedMedian, sizeof(float), 1, f) != 1 ||
            !(node->normal = malloc(dataSize * sizeof(float)))) {
            free(node);
            return NULL;
        }

        if (fread(node->normal, sizeof(float), dataSize, f) != dataSize) {
            free(node->normal);
            free(node);
            return NULL;
        }

        node->left = deserializeTree(f, dataSize);
        node->right = deserializeTree(f, dataSize);
        node->parent = NULL;
        node->data = NULL;
        node->size = -1;

        if (node->left) node->left->parent = node;
        if (node->right) node->right->parent = node;
        return node;
    }
    return NULL; // Invalid marker
}

void saveTree(Node *root, const char *path, int dataSize) {
    FILE *f = fopen(path, "wb");
    if (!f) {
        perror("fopen");
        return;
    }
    serializeTree(root, f, dataSize);
    fclose(f);
}

Node *loadTree(const char *path, int dataSize) {
    FILE *f = fopen(path, "rb");
    if (!f) {
        perror("fopen");
        return NULL;
    }
    Node *root = deserializeTree(f, dataSize);
    fclose(f);
    return root;
}

int main(int argc, char **argv) {
    /**
     * `--create [db name] -f [file path] -s [size of vectors] -m [metadata in json]` creates vectordb named db name and uses the file containing the vectors on file path
     * `--load [db name]` loads db with db name
     * `--save [db name]` saves current db
     * `--search [db name] -v [vector path or vector] -k [num of neighbors] [opt] -similarity "cosine"||"eucledian"` sear
     * `--add [db name] -v [vector path or vector] -m [metadata in json]
     *
     * NOTES:
     * 1 - Serialization and Deserialization dont work
     *
     *
     * TODO:
     * 1. add collection metadata serialization ie metadata fields, datasize, collection name, paths, user defined constraints.
     * 2. add sqlite to store vector uuid, vector and metadata.
     * 3. add CRUD methods for metadata and vectors.
     * 4. add other similarity functions.
     * 5. define python and dart wrappers.
     */

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <datafile>\n", argv[0]);
        return -1;
    }
    const char *path = argv[1];
    FILE *fp = fopen(path, "r");
    if (fp == NULL) {
        perror("fopen");
        return -1;
    }

    long count = 0;
    int ch;
    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\n') {
            count++;
        }
    }
    if (ferror(fp)) {
        perror("reading for line count");
        fclose(fp);
        return -1;
    }

    const int DATASIZE = 512;
    float **vectors = malloc(sizeof *vectors * count);
    if (vectors == NULL) {
        perror("malloc vectors array");
        fclose(fp);
        return -1;
    }

    for (long j = 0; j < count; j++) {
        vectors[j] = malloc(sizeof *vectors[j] * DATASIZE);
        if (vectors[j] == NULL) {
            perror("malloc vectors[j]");
            while (j-- > 0) free(vectors[j]);
            free(vectors);
            fclose(fp);
            return -1;
        }
    }

    rewind(fp);
    char token[DATASIZE];
    int pos = 0;
    long i = 0, j = 0;
    while ((ch = fgetc(fp)) != EOF && j < count) {
        if (isspace(ch)) {
            if (pos > 0) {
                token[pos] = '\0';
                char *endptr;
                float val = strtof(token, &endptr);
                if (*endptr == '\0') {
                    if (i >= DATASIZE) {
                        fprintf(stderr, "ERROR: vector %ld too many elements\n", j);
                        for (long idx = 0; idx < count; idx++) {
                            free(vectors[idx]);
                        }
                        free(vectors);
                    }
                    vectors[j][i++] = val;
                }
                pos = 0;
            }
            if (ch == '\n') {
                j++;
                i = 0;
            }
        } else if (pos < DATASIZE - 1) {
            token[pos++] = ch;
        }
    }

    if (ferror(fp)) {
        perror("parsing floats");
        for (long idx = 0; idx < count; idx++) {
            free(vectors[idx]);
        }
        free(vectors);
    }
    fclose(fp);

    VectorPair **values = malloc(sizeof *values * count);
    if (values == NULL) {
        perror("malloc values");
        for (long idx = 0; idx < count; idx++) {
            free(vectors[idx]);
        }
        free(vectors);
    }
    for (long idx = 0; idx < count; idx++) {
        values[idx] = createVectorPair(vectors[idx]);
    }


    Node *head = createNode(values, (int) count);
    if (head == NULL) {
        fprintf(stderr, "Failed to create root node\n");
        free(values);
        for (long idx = 0; idx < count; idx++) {
            free(vectors[idx]);
        }
        free(vectors);
    }

    addSplit(head, DATASIZE);
    saveTree(head, "tree.bin", DATASIZE);
    ScorePair *results, *resultsLoad;
    int size;
    results = searchTopK(head, vectors[0], 10, DATASIZE, &size);
    for (int i = 0; i < size; i++) {
        printf("%s\n", results[i].uuid);
    }
    freeTree(head, DATASIZE);

    head = loadTree("tree.bin", DATASIZE);
    resultsLoad = searchTopK(head, vectors[0], 10, DATASIZE, &size);
    for (int i = 0; i < size; i++) {
        int same = (strcmp(resultsLoad[i].uuid,
                           results[i].uuid) == 0);
        printf("loaded[%d]=%s  original[%d]=%s  same=%d\n",
               i, resultsLoad[i].uuid,
               i, results[i].uuid,
               same);
    }

    free(resultsLoad);
    freeTree(head, DATASIZE);

    free(results);



    for (long idx = 0; idx < count; idx++) {
        free(vectors[idx]);
    }
    free(vectors);

    printf("Processed %ld vectors (lines).\n", count);
    return 0;
}
