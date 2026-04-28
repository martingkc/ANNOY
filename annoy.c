#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <sqlite3.h>



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

// Like addVector but takes a pre-created VectorPair so the caller knows the UUID upfront.
void addVectorPair(Node *self, VectorPair *vp, int dataSize) {
    VectorPair **tmpData;
    float projection;

    if (self->data == NULL) {
        projection = computeProjection(self->normal, vp->value, dataSize);
        if (projection > self->indexedMedian) {
            addVectorPair(self->right, vp, dataSize);
        } else {
            addVectorPair(self->left, vp, dataSize);
        }
    } else {
        tmpData = malloc((self->size + 1) * sizeof(VectorPair *));
        if (tmpData == NULL) return;
        for (int i = 0; i < self->size; i++) {
            tmpData[i] = self->data[i];
        }
        tmpData[self->size] = vp;
        free(self->data);
        self->data = tmpData;
        self->size++;
        if (self->size > K) {
            addSplit(self, dataSize);
        }
    }
}

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
    if (self == NULL)
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

// ============================================================
// SQLite metadata layer
// ============================================================

static char *buildPath(const char *name, const char *ext) {
    size_t len = strlen(name) + strlen(ext) + 1;
    char *path = malloc(len);
    if (path) snprintf(path, len, "%s%s", name, ext);
    return path;
}

sqlite3 *openMetadataDB(const char *collection_name) {
    char *path = buildPath(collection_name, ".db");
    if (!path) return NULL;

    sqlite3 *db;
    if (sqlite3_open(path, &db) != SQLITE_OK) {
        fprintf(stderr, "Cannot open DB %s: %s\n", path, sqlite3_errmsg(db));
        sqlite3_close(db);
        free(path);
        return NULL;
    }
    free(path);
    return db;
}

void closeMetadataDB(sqlite3 *db) {
    if (db) sqlite3_close(db);
}

int initCollectionDB(sqlite3 *db, const char *name, int dataSize) {
    const char *sql =
        "CREATE TABLE IF NOT EXISTS collection_info ("
        "  key   TEXT PRIMARY KEY,"
        "  value TEXT NOT NULL"
        ");"
        "CREATE TABLE IF NOT EXISTS vectors ("
        "  uuid     TEXT PRIMARY KEY,"
        "  metadata TEXT"
        ");";

    char *err = NULL;
    if (sqlite3_exec(db, sql, NULL, NULL, &err) != SQLITE_OK) {
        fprintf(stderr, "SQL error (init): %s\n", err);
        sqlite3_free(err);
        return -1;
    }

    const char *infoSql =
        "INSERT OR REPLACE INTO collection_info (key, value) VALUES (?, ?);";
    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db, infoSql, -1, &stmt, NULL) != SQLITE_OK) return -1;

    sqlite3_bind_text(stmt, 1, "name", -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 2, name, -1, SQLITE_STATIC);
    sqlite3_step(stmt);
    sqlite3_reset(stmt);

    char sizeStr[32];
    snprintf(sizeStr, sizeof(sizeStr), "%d", dataSize);
    sqlite3_bind_text(stmt, 1, "data_size", -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 2, sizeStr, -1, SQLITE_STATIC);
    sqlite3_step(stmt);

    sqlite3_finalize(stmt);
    return 0;
}

int getCollectionDataSize(sqlite3 *db, int *dataSize) {
    const char *sql = "SELECT value FROM collection_info WHERE key = 'data_size';";
    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL) != SQLITE_OK) return -1;

    int rc = -1;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        *dataSize = atoi((const char *) sqlite3_column_text(stmt, 0));
        rc = 0;
    }
    sqlite3_finalize(stmt);
    return rc;
}

int insertMetadata(sqlite3 *db, const char *uuid, const char *metadata_json) {
    const char *sql =
        "INSERT OR REPLACE INTO vectors (uuid, metadata) VALUES (?, ?);";
    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL) != SQLITE_OK) return -1;

    sqlite3_bind_text(stmt, 1, uuid, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 2, metadata_json ? metadata_json : "null",
                     -1, SQLITE_STATIC);
    int rc = sqlite3_step(stmt) == SQLITE_DONE ? 0 : -1;
    sqlite3_finalize(stmt);
    return rc;
}

// Returns a malloc'd string; caller must free. Returns NULL if not found.
char *getMetadata(sqlite3 *db, const char *uuid) {
    const char *sql = "SELECT metadata FROM vectors WHERE uuid = ?;";
    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL) != SQLITE_OK) return NULL;

    sqlite3_bind_text(stmt, 1, uuid, -1, SQLITE_STATIC);
    char *result = NULL;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        const char *text = (const char *) sqlite3_column_text(stmt, 0);
        if (text) result = strdup(text);
    }
    sqlite3_finalize(stmt);
    return result;
}

int deleteMetadata(sqlite3 *db, const char *uuid) {
    const char *sql = "DELETE FROM vectors WHERE uuid = ?;";
    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL) != SQLITE_OK) return -1;

    sqlite3_bind_text(stmt, 1, uuid, -1, SQLITE_STATIC);
    int rc = sqlite3_step(stmt) == SQLITE_DONE ? 0 : -1;
    sqlite3_finalize(stmt);
    return rc;
}

int updateMetadata(sqlite3 *db, const char *uuid, const char *metadata_json) {
    const char *sql = "UPDATE vectors SET metadata = ? WHERE uuid = ?;";
    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL) != SQLITE_OK) return -1;

    sqlite3_bind_text(stmt, 1, metadata_json, -1, SQLITE_STATIC);
    sqlite3_bind_text(stmt, 2, uuid, -1, SQLITE_STATIC);
    int rc = sqlite3_step(stmt) == SQLITE_DONE ? 0 : -1;
    sqlite3_finalize(stmt);
    return rc;
}

// ============================================================
// File helpers
// ============================================================

// Reads space/newline separated floats; one vector per line.
static float **readVectorsFromFile(const char *path, int dataSize, long *outCount) {
    FILE *fp = fopen(path, "r");
    if (!fp) { perror("fopen"); return NULL; }

    long count = 0;
    int ch;
    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\n') count++;
    }
    if (count == 0) { fclose(fp); *outCount = 0; return NULL; }
    rewind(fp);

    float **vectors = malloc(sizeof *vectors * count);
    if (!vectors) { fclose(fp); return NULL; }

    for (long j = 0; j < count; j++) {
        vectors[j] = malloc(sizeof(float) * dataSize);
        if (!vectors[j]) {
            while (j-- > 0) free(vectors[j]);
            free(vectors);
            fclose(fp);
            return NULL;
        }
    }

    char token[64];
    int pos = 0;
    long i = 0, j = 0;
    while ((ch = fgetc(fp)) != EOF && j < count) {
        if (isspace(ch)) {
            if (pos > 0) {
                token[pos] = '\0';
                char *endptr;
                float val = strtof(token, &endptr);
                if (*endptr == '\0' && i < dataSize) {
                    vectors[j][i++] = val;
                }
                pos = 0;
            }
            if (ch == '\n') { j++; i = 0; }
        } else if (pos < (int) sizeof(token) - 1) {
            token[pos++] = ch;
        }
    }
    fclose(fp);

    *outCount = count;
    return vectors;
}

// Reads one JSON metadata string per line (matches vector file line-for-line).
static char **readMetadataFromFile(const char *path, long count) {
    FILE *fp = fopen(path, "r");
    if (!fp) { perror("fopen"); return NULL; }

    char **metadata = calloc(count, sizeof(char *));
    if (!metadata) { fclose(fp); return NULL; }

    char line[4096];
    long i = 0;
    while (i < count && fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        if (len > 0 && line[len - 1] == '\n') line[len - 1] = '\0';
        metadata[i++] = strdup(line);
    }
    fclose(fp);
    return metadata;
}

static void freeVectors(float **vectors, long count) {
    for (long i = 0; i < count; i++) free(vectors[i]);
    free(vectors);
}

static void freeMetadata(char **metadata, long count) {
    if (!metadata) return;
    for (long i = 0; i < count; i++) free(metadata[i]);
    free(metadata);
}

// ============================================================
// CLI entry point
// ============================================================

int main(int argc, char **argv) {
    /*
     * --create <collection> <vector_file> <metadata_file> -d <dim>
     *     Build a new ANN index from vector_file.  One JSON string per line
     *     in metadata_file is stored in SQLite keyed by vector UUID.
     *
     * --search <collection> <vector_file> -k <topK>
     *     Load index, run ANN search with the first vector in vector_file,
     *     and print results with their metadata.
     *
     * --add <collection> <vector_file> -m <metadata_file>
     *     Insert new vectors (with metadata) into an existing collection.
     *
     * --remove <collection> <uuid>
     *     Delete a vector's metadata entry from the DB by UUID.
     *
     * --update <collection> <uuid> <metadata_json>
     *     Replace the metadata for an existing UUID.
     */

    if (argc < 2) {
        fprintf(stderr,
                "Usage:\n"
                "  --create <collection> <vector_file> <metadata_file> -d <dim>\n"
                "  --search <collection> <vector_file> -k <topK>\n"
                "  --add    <collection> <vector_file> -m <metadata_file>\n"
                "  --remove <collection> <uuid>\n"
                "  --update <collection> <uuid> <metadata_json>\n");
        return 1;
    }

    // --create <collection> <vector_file> <metadata_file> -d <dim>
    if (strcmp(argv[1], "--create") == 0) {
        if (argc != 7 || strcmp(argv[5], "-d") != 0) {
            fprintf(stderr,
                    "Usage: --create <collection> <vector_file> <metadata_file> -d <dim>\n");
            return 1;
        }
        const char *collection = argv[2];
        const char *vecFile    = argv[3];
        const char *metaFile   = argv[4];
        int dataSize = atoi(argv[6]);
        if (dataSize <= 0) { fprintf(stderr, "Invalid dim\n"); return 1; }

        long count = 0;
        float **vectors = readVectorsFromFile(vecFile, dataSize, &count);
        if (!vectors) return 1;

        char **metadata = readMetadataFromFile(metaFile, count);

        VectorPair **values = malloc(sizeof *values * count);
        if (!values) {
            freeMetadata(metadata, count);
            freeVectors(vectors, count);
            return 1;
        }
        for (long idx = 0; idx < count; idx++) {
            values[idx] = createVectorPair(vectors[idx]);
        }

        sqlite3 *db = openMetadataDB(collection);
        if (!db || initCollectionDB(db, collection, dataSize) != 0) {
            free(values);
            freeMetadata(metadata, count);
            freeVectors(vectors, count);
            return 1;
        }

        // Insert metadata before addSplit (addSplit frees the values[] array on root).
        sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, NULL);
        for (long idx = 0; idx < count; idx++) {
            const char *meta = (metadata && metadata[idx]) ? metadata[idx] : "null";
            insertMetadata(db, values[idx]->uuid, meta);
        }
        sqlite3_exec(db, "COMMIT;", NULL, NULL, NULL);

        Node *root = createNode(values, (int) count);
        if (!root) {
            closeMetadataDB(db);
            freeMetadata(metadata, count);
            freeVectors(vectors, count);
            return 1;
        }
        addSplit(root, dataSize);

        char *tp = buildPath(collection, ".bin");
        saveTree(root, tp, dataSize);
        free(tp);

        freeTree(root, dataSize);
        closeMetadataDB(db);
        freeMetadata(metadata, count);
        freeVectors(vectors, count);

        printf("Created collection '%s' with %ld vectors (dim=%d).\n",
               collection, count, dataSize);
        return 0;
    }

    // --search <collection> <vector_file> -k <topK>
    if (strcmp(argv[1], "--search") == 0) {
        if (argc != 6 || strcmp(argv[4], "-k") != 0) {
            fprintf(stderr, "Usage: --search <collection> <vector_file> -k <topK>\n");
            return 1;
        }
        const char *collection = argv[2];
        const char *vecFile    = argv[3];
        int topK = atoi(argv[5]);
        if (topK <= 0) { fprintf(stderr, "Invalid topK\n"); return 1; }

        sqlite3 *db = openMetadataDB(collection);
        if (!db) return 1;

        int dataSize = 0;
        if (getCollectionDataSize(db, &dataSize) != 0) {
            fprintf(stderr, "Cannot read collection info\n");
            closeMetadataDB(db);
            return 1;
        }

        long count = 0;
        float **vectors = readVectorsFromFile(vecFile, dataSize, &count);
        if (!vectors || count == 0) { closeMetadataDB(db); return 1; }

        char *tp = buildPath(collection, ".bin");
        Node *root = loadTree(tp, dataSize);
        free(tp);
        if (!root) {
            closeMetadataDB(db);
            freeVectors(vectors, count);
            return 1;
        }

        int size = 0;
        ScorePair *results = searchTopK(root, vectors[0], topK, dataSize, &size);

        for (int i = 0; i < size; i++) {
            char *meta = getMetadata(db, results[i].uuid);
            printf("[%d] uuid=%s  score=%.4f  metadata=%s\n",
                   i + 1, results[i].uuid, results[i].key,
                   meta ? meta : "null");
            free(meta);
        }

        free(results);
        freeTree(root, dataSize);
        closeMetadataDB(db);
        freeVectors(vectors, count);
        return 0;
    }

    // --add <collection> <vector_file> -m <metadata_file>
    if (strcmp(argv[1], "--add") == 0) {
        if (argc != 6 || strcmp(argv[4], "-m") != 0) {
            fprintf(stderr, "Usage: --add <collection> <vector_file> -m <metadata_file>\n");
            return 1;
        }
        const char *collection = argv[2];
        const char *vecFile    = argv[3];
        const char *metaFile   = argv[5];

        sqlite3 *db = openMetadataDB(collection);
        if (!db) return 1;

        int dataSize = 0;
        if (getCollectionDataSize(db, &dataSize) != 0) {
            fprintf(stderr, "Cannot read collection info\n");
            closeMetadataDB(db);
            return 1;
        }

        long count = 0;
        float **vectors = readVectorsFromFile(vecFile, dataSize, &count);
        if (!vectors) { closeMetadataDB(db); return 1; }

        char **metadata = readMetadataFromFile(metaFile, count);

        char *tp = buildPath(collection, ".bin");
        Node *root = loadTree(tp, dataSize);
        if (!root) {
            free(tp);
            closeMetadataDB(db);
            freeMetadata(metadata, count);
            freeVectors(vectors, count);
            return 1;
        }

        sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, NULL);
        for (long idx = 0; idx < count; idx++) {
            VectorPair *vp = createVectorPair(vectors[idx]);
            if (!vp) continue;
            const char *meta = (metadata && metadata[idx]) ? metadata[idx] : "null";
            insertMetadata(db, vp->uuid, meta);
            addVectorPair(root, vp, dataSize);
        }
        sqlite3_exec(db, "COMMIT;", NULL, NULL, NULL);

        saveTree(root, tp, dataSize);
        free(tp);

        freeTree(root, dataSize);
        closeMetadataDB(db);
        freeMetadata(metadata, count);
        freeVectors(vectors, count);

        printf("Added %ld vector(s) to '%s'.\n", count, collection);
        return 0;
    }

    // --remove <collection> <uuid>
    if (strcmp(argv[1], "--remove") == 0) {
        if (argc != 4) {
            fprintf(stderr, "Usage: --remove <collection> <uuid>\n");
            return 1;
        }
        const char *collection = argv[2];
        const char *uuid       = argv[3];

        sqlite3 *db = openMetadataDB(collection);
        if (!db) return 1;

        if (deleteMetadata(db, uuid) == 0) {
            printf("Removed UUID %s from '%s'.\n", uuid, collection);
        } else {
            fprintf(stderr, "Failed to remove UUID %s\n", uuid);
            closeMetadataDB(db);
            return 1;
        }
        closeMetadataDB(db);
        return 0;
    }

    // --update <collection> <uuid> <metadata_json>
    if (strcmp(argv[1], "--update") == 0) {
        if (argc != 5) {
            fprintf(stderr, "Usage: --update <collection> <uuid> <metadata_json>\n");
            return 1;
        }
        const char *collection    = argv[2];
        const char *uuid          = argv[3];
        const char *metadata_json = argv[4];

        sqlite3 *db = openMetadataDB(collection);
        if (!db) return 1;

        if (updateMetadata(db, uuid, metadata_json) == 0) {
            printf("Updated metadata for UUID %s in '%s'.\n", uuid, collection);
        } else {
            fprintf(stderr, "Failed to update UUID %s\n", uuid);
            closeMetadataDB(db);
            return 1;
        }
        closeMetadataDB(db);
        return 0;
    }

    fprintf(stderr, "Unknown command: %s\n", argv[1]);
    return 1;
}
