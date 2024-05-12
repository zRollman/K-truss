//
// Created by Administrator on 2023/4/17 0017.
//

#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#define MAX_ADJ 10000
#define MAX_EDGE 300000
#define MAX_MX 10


static int kpgn[MAX_EDGE][MAX_MX];
static int npgn[MAX_EDGE][MAX_MX];
static int ncount;
static int adjUsed[MAX_ADJ];
static int kusdpgn[MAX_ADJ][MAX_ADJ];



static int **adj_list;

static int isSyc = 0;

static unsigned kcount = 0;
static int trow ;
static int omp_n;
static int tk;
static int kused = 0;

void QuickSort(int* a, int left, int right);
int findVal(int* arr, int val, int m);
int matchArrs(int* a, int* b, int m);

void findkPolygonCount();
void printkPolygon();

void run(int triangle_k,int row, int sy,int nz,int **triples,int len,int ompn);