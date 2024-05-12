//
// Created by Administrator on 2023/4/17 0017.
//

#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#define MAX_ADJ 10000 ///vertex count
#define MAX_EDGE 300000 ///edge count
#define MAX_MX 10 ///k-truss count

typedef struct {
    int index;
    int selected_so_far[MAX_MX];
    int selected_count;
} State;

//int col = 56;

static int kpgn[MAX_EDGE][MAX_MX];
static int npgn[MAX_EDGE][MAX_MX];
static int PGN[2*MAX_EDGE][MAX_MX];

static unsigned ncount = 0;
static unsigned kcount = 0;
static unsigned pcount = 0;

static int **adj_list;

static int isSyc = 0;


static int trow ;
static int omp_n;
static int tk;

void findkPolygonCount();

void run(int triangle_k,int row, int sy,int nz,int **triples,int len,int ompn);
void BFS(int arr[],int n, int k);
int matchAdjs(int temp[],int k,int begin);
int matchArrs(int arr[]);
int IsRepeatKtruss(int *temp,int count, int k);
int IsRepeatNtruss(int *temp,int count, int k);
void DelIncludePgn();
