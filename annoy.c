#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define K 10

struct Node
{ 
float** data; // list of vectors 
int** size;
struct node *left; 
struct node *right;
struct node *parent;
float* normal; // normal vector for the sector 
float indexedMedian; // decision val
};

int comp (const void * elem1, const void * elem2) 
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
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
	
	nDotv = dotProd(normal, vector);
	nDotn = dotProd(normal, normal); 

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
		return NULL;
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
		return NULL;
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
		return NULL; 
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
		if(tmpData == NULL) return NULL;
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
float** search(Node* self, float* vector, int topK, int dataSize){

	if(self->data == NULL){
		projection = computeProjection(self->normal, vector, dataSize);
		if (projection > self->indexedMedian){
			search(self->right, vector,topK, dataSize);
		}else{
			search(self->left, vector, topK, dataSize);			
		}
	}else{

	}
}









int main() {
	printf("hello world \n");
	return 0;
}