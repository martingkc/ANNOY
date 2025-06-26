#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define K 10

typedef struct  Node
{ 
float** data; // list of vectors 
int size;
struct Node *left; 
struct Node *right;
struct Node *parent;
float* normal; // normal vector for the sector 
float indexedMedian; // decision val
}Node;

typedef struct List{
	float** data; 
	int size;
	struct List *next; 
} List;

int comp (const void * elem1, const void * elem2) 
{
    float f = *((float*)elem1);
    float s = *((float*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
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

Node* createNode(float** data, int size)
{
    Node* newNode = (Node*)malloc(sizeof(Node));
    newNode->data = data;
    newNode->size = size;
    newNode->left = NULL;
    newNode->right = NULL;
    newNode->parent = NULL; 
    return newNode;
}

List* createListNode(float** data, int size)
{
    List* newList = (List*)malloc(sizeof(List));
    newList->data = data;
    newList->size = size;
    newList->next = NULL;
    return newList;
}

void pushToList(List* head,float** data, int size){
	List *tmp;
	List *newItem;
	if(head->data == NULL || head->size == -1){
		head->data = data;
		head->size = size;
	}else{
		tmp = head;
		newItem = createListNode(data, size);
		while(tmp->next != NULL){
			tmp = tmp->next;
		}

		tmp->next = newItem;
		newItem->next = NULL;
	}

	return;
}
 
float dotProd(float* a, float* b, int size){

	float sum; 
	sum = 0.0f;

	for(int i = 0; i<size; i++){
		sum += a[i]*b[i];
	}

	return sum; 
}



float computeProjection(float* normal, float* vector, int sizeVectors){
	
	float nDotv; 
	float nDotn; 
	float factor; 
	float tempRes; 
	float result; 

	result = 0;
	
	nDotv = dotProd(normal, vector, sizeVectors);
	nDotn = dotProd(normal, normal, sizeVectors); 

	factor = nDotv / nDotn; 

	for(int i = 0; i<sizeVectors; i++){
		tempRes = normal[i]*factor;
		result += tempRes* tempRes;
	}

	return sqrt(result);

}


float* determineNormal(float *a, float *b, int size){
	float* normal;

	normal = malloc(size*sizeof(float)); 

	if(normal == NULL){
		return NULL; 
	}

	for(int i=0; i<size; i++){
		normal[i] = a[i]-b[i]; 
	}

	return normal;
}


float calculateMedian(float *projections, int size){
	float median; 
	median = 0.0f;
	qsort(projections, size, sizeof(float), comp);

	if(size%2 == 0){
		median += projections[size/2]; 
		median += projections[size/2 - 1];
		median /= 2;
	}else{
		median += projections[(size + 1)/2 - 1]; 
	}

	return median; 
}




void addSplit(Node* self, int dataSize){
	float* a; 
	float* b;
	float* normal; 
	float* projections; 
	float** dataL; 
	float** dataR;
	float** data;
	float median;
	int idxA; 
	int idxB; 
	int size; 
	int sizeL;
	int sizeR;
	int cnt;
	int idxL, idxR;
	idxL=0;
	idxR=0;
	size = self->size; 
	data = self->data;
	cnt = 0;

	projections = malloc(sizeof(float)*size);

	if(projections == NULL){
		return;
	}

	idxA = rand()%size;
	idxB = rand()%size; 

	if(idxA == idxB){
		idxB = rand()%size; // if we get 2 times the same values o zmn sictik 
	}

	a = data[idxA]; 
	b = data[idxB]; 

	normal = determineNormal(a, b, dataSize); 

	for(int i=0; i<size; i++){
		projections[i] = computeProjection(normal, data[i], dataSize);
	}

	median = calculateMedian(projections, size); 

	if(size%2 == 0){
		sizeL = size/2; // Left block 
		sizeR = size/2; // Right block
	}else{
		for(int i = 0; i<size; i++){
			if(projections[i]<= median){
				cnt++;
			}
		}

		sizeL = cnt;
		sizeR = size - cnt;
	}

	dataL = malloc(sizeL * sizeof(float*));
	dataR = malloc(sizeR * sizeof(float*));

	if(dataL == NULL || dataR == NULL){
		free(dataL);
		free(dataR);
		return;
	}



	for(int i = 0; i<size; i++){
		if(projections[i]<= median){
			dataL[idxL] = data[i];
			idxL += 1;
		}else{
			dataR[idxR] = data[i];
			idxR += 1; 
		}
	}

	self->normal = normal;
	self->indexedMedian = median; 

	self->left = createNode(dataL, sizeL); 
	self->right = createNode(dataR, sizeR);

	if(self->left == NULL || self->right == NULL){
		return; 
	}

	self->right->parent = self;
	self->left->parent  = self->right->parent;

	// If the size of a node is greater than K it gets split so it means that they arent leafs
	if(sizeL > K ){
		addSplit(self->left, dataSize);
	}
	if(sizeR > K ){
		addSplit(self->right, dataSize);
	}



	// Nodes that aren't leafs don't get the size or data fields filled. 
	
	free(self->data);
	free(projections);

	self->data = NULL;
	self->size = -1;

}

void addVector(Node* self, float* vector, int dataSize){
	float** tmpData; 
	float projection;

	// this would be better with a check on the leafs 
	if(self->data == NULL){
		projection = computeProjection(self->normal, vector, dataSize);
		if (projection > self->indexedMedian){
			addVector(self->right, vector, dataSize);
		}else{
			addVector(self->left, vector, dataSize);			
		}
	}else{
		tmpData = malloc((self->size+1)*sizeof(float*)); 
		if(tmpData == NULL) return;
		for(int i = 0; i<(self->size + 1); i++){
			if(i !=(self->size)){
				tmpData[i] = self->data[i];
			}else{
				tmpData[i] = vector;
			}
		}

		free(self->data); 
		self->data = tmpData;

		self->size = self->size + 1; 
		if(self->size > K){
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

Node* recursiveNodeSearch(Node* self, float* vector, int dataSize){
	float** tmpData; 
	float projection;

	// this would be better with a check on the leafs 
	if(self->data == NULL){
		projection = computeProjection(self->normal, vector, dataSize);
		if (projection > self->indexedMedian){
			return recursiveNodeSearch(self->right, vector, dataSize);
		}else{
			return recursiveNodeSearch(self->left, vector, dataSize);		
		}
	}else{
		return self; 
	}
}


void collectLevel(Node* self, List* head){
	if(self->data == NULL){
		 collectLevel(self->right, head);
		 collectLevel(self->left, head);		
	}else{
		pushToList(head, self->data,self->size);
	}
}

float** getVectorsList(List* head, int dataSize, int *retSize) {
    List* tmp;
    float **vectors;
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

float** searchTopK(Node* self, float* vector, int topK, int dataSize, int *size){
	List* head; // list of vectors holding the results 

	Node* startNode; // starting node from which to start trasversing
	float** vectors;
	int nodesToFind; 
	head = createListNode(NULL, -1); 

	startNode = recursiveNodeSearch(self, vector, dataSize); 
	nodesToFind = (int)sqrt((topK/(float)K)*1.25f) + 1;

	if (startNode->parent == NULL) return NULL;

	do{
		startNode = startNode->parent; 
		nodesToFind -= 1;
	}while(nodesToFind>=0 && startNode->parent != NULL);

	collectLevel(startNode, head);

	deleteList(head); 

	return getVectorsList(head, dataSize, size);

}









int main() {
	printf("hello world \n");
	return 0;
}