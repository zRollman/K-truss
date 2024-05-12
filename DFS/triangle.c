
#include "triangle.h"
#include "mmio.h"
#include "omp.h"
typedef struct bool bool;

omp_lock_t lock;

// 交换元素
void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

int maxValue(int a,int b) {
    return a>b?a:b;
}
void swapArr(int *a,int *b,int len){
//    printf("len %d:\n",len);
    for (int k = 0; k <= len; k++) {
//        printf("a %d,b %d\n",a[k],b[k]);
        swap(&a[k], &b[k]);
    }
}

void quickSort(int arr[][MAX_MX], int left, int right) {
    if (left >= right) {
        return;
    }

    int pivot = arr[(left + right) / 2][0];
//    printf("pivot %d\n",pivot);
    int i = left;
    int j = right;

    while (i <= j) {
        while (arr[i][0] < pivot) {
            i++;
        }
        while (arr[j][0] > pivot) {
            j--;
        }
        if (i <= j) {

            swapArr(arr[i],arr[j], maxValue(arr[i][0],arr[j][0]));
            i++;
            j--;
        }
    }

    quickSort(arr, left, j);
    quickSort(arr, i, right);
}

int singleCmp(const void* a, const void* b)
{
    return  *(int*)a - *(int*)b ;
}

int doubleCmp(const void* a, const void* b){
    return((int*)a)[0]-((int*)b)[0];
}

int findVal(int* arr, int val, int m) {////Find the same value
    for (int i = 0; i < m; ++i) {
        if (arr[i] == val) return 1;
    }
    return 0;
}

int matchArrs(int* a, int* b, int m) { /// the array appear
    for (int i = 0; i < m; ++i) {
        if (a[i] != b[i]) return 0;
    }
    return 1;
}
int nmatchArrs(int *a,int *b,int m){/// the array appear
    for(int i=1;i<=m;++i){
        if(a[i]!=b[i-1]) return 0;
    }
    return 1;
}

int matchIncludeArrs(int *a,int *b,int m,int i) {
    //a 为 pgn  , b为 tmp
    int j = 0;
    for (; i <= tk && j <= m;) {
        if (a[i] == b[j]) {
            ++i;
            ++j;
        } else i++;
    }
    if (j > m) return 1;
    return 0;
}



void initAdjList(int **triples,int len) {
    printf("%d,%d\n",trow,len);

    for(int i=0;i<trow+1;++i){
        if(isSyc) {
            adj_list[i] = (int*) malloc((5+2*trow)*sizeof (int));
            memset(adj_list[i],0,(5+2*trow)*sizeof (int));
        }
        else {
            adj_list[i] = (int*) malloc((5+trow)*sizeof (int));
            memset(adj_list[i],0,(5+trow)*sizeof (int));
        }
        if(adj_list == NULL) {
            printf("malloc error!!\n");
            exit(1);
        }
    }

    int tmpLen = len;
    for(int i=0;i<len;++i){
        int adjx = triples[i][0];
        int adjy = triples[i][1];
        int loc = adj_list[adjx][0]+1;
        ///Judging duplicate vertices
        int repeatAdj = 0;
        for(int j=1;j<loc;++j){
            if(adjy == adj_list[adjx][j]) {
                repeatAdj = 1;
                tmpLen--;
                break;
            }
        }
        if(!repeatAdj){
            adj_list[adjx][loc] = adjy;
            adj_list[adjx][0]= loc;
        }

    }
    len = tmpLen;

    if(isSyc){
        ///Turn a directed graph into an undirected graph.
        ///traverse the matrix
        for (int i = 0; i < trow+1; ++i) {
            ///Iterate over all vertices connected to vertex i and join vertex i at the connected vertices
            for(int j=1;j<=adj_list[i][0];++j){
                if(i!=adj_list[i][j]){
                    int v = adj_list[i][j];
                    ///Determine if a v-vertex also contains an i-vertex
                    int cnt = 0;
                    for(int k = 1;k<=adj_list[v][0];++k){
                        if(adj_list[v][k]==i){
                            cnt = 1;
                            break;
                        }
                    }
                    if(!cnt){
                        adj_list[v][++adj_list[v][0]] = i;
                    }
                }
            }
        }
    }


    ///output
    for (int i = 1; i < trow; ++i) {
        printf("%d,%d : ", i,adj_list[i][0]);

        for (int j = 1; j < adj_list[i][0] + 1; ++j) {
            printf("%d-", adj_list[i][j]);
        }
        printf("\n");
    }

}


///
int dfs(int j, int adj[], int m,int visited[],int second) {

    if (m == tk) {
        qsort(adj, m, sizeof(int), singleCmp);
        int cnt = 1;
//        omp_set_lock(&lock);
        ///Excluding the same
        if(!second){
            for (int i = 0; i < kcount; ++i) {
                if (matchArrs(kpgn[i], adj, m)) cnt = 0;
                if (!cnt) break;
            }
            if (cnt) {
                for (int i = 0; i < m; ++i) {
                    kpgn[kcount][i] = adj[i];
                printf("%d,",adj[i]);
                }
                kcount++;
            printf("k-count : %d\n", kcount);
            }
        }else{
            for(int i=0;i<kused;++i){
                if(matchIncludeArrs(kusdpgn[i],adj,m,0)) cnt =0;
                if(!cnt) break;
            }

            if(cnt){
                for(int i=0;i<ncount;++i){
                    if(npgn[i][0] == m){
                        if(nmatchArrs(npgn[i],adj,m)) cnt = 0;
                        if(!cnt) break;
                    }
                }
                if(cnt){
                    //去重 -- kused-truss
                    for(int k = 0;k<kused;++k){
                        if(matchIncludeArrs(kusdpgn[k],adj,m,0)){
                            cnt =0;
                            break;
                        }
                    }
                    if(cnt){
                        npgn[ncount][0] = m;
//                        printf("npgn :\n");
                        for(int i=0;i<m;++i){
                            npgn[ncount][i+1] = adj[i];
//                            printf("%d,",adj[i]);
                        }
//                        printf("\n");
                        ncount++;
                    }
                }

            }
        }

//        omp_unset_lock(&lock);
        return 1;
    }
    int temp[MAX_MX];
    int w = 0;///Determine if dfs is ready to go deeper
    for (int i = 0; i < m; ++i) {
        temp[i] = adj[i];
    }
    for (int i = 1; i < adj_list[j][0] + 1; ++i) {
        int v = adj_list[j][i];
//        printf("j:%d v:%d\n",j,v);

//        if (adjUsed[v]) continue;

        if (visited[v] || findVal(adj, v, m)) continue;
        int acount = 0;///Calculate the number of vertices in v that are connected to the vertices in adj
        for (int ii = 1; ii < adj_list[v][0] + 1; ++ii) {
            ///Whether v is connected to a vertex in adj
            if (findVal(adj, adj_list[v][ii], m)) acount++;
        }
        if (acount == m) {
            adj[m] = v;
        } else continue;

        int cnt = 1;
//        omp_set_lock(&lock);
        ///Excluding the same

        w = 1; // One more layer than the last
        visited[v] = 1;
        int d = dfs(v, adj, m + 1, visited, second);

        int tmp[MAX_MX];
        for (int ii = 0; ii < m; ++ii) {
            adj[ii] = temp[ii];
            tmp[ii] = temp[ii];
        }
        tmp[m] = v;
        visited[v] = 0;
        if (d) continue;

        if (m >= 1) {///When the absence of the next edge deposits,at least two vertices

            qsort(tmp, m + 1, sizeof(int), singleCmp);
            int ncnt = 1, kcnt = 1;
            for (int n = 0; n < maxValue(kcount, ncount); ++n) {
                if (0 == second)///if it is first dfs
                    if (n < kcount && matchIncludeArrs(kpgn[n], tmp, m, 0)) kcnt = 0;
                if (n < ncount && matchIncludeArrs(npgn[n], tmp, m, 1)) ncnt = 0;
                if (!ncnt || !kcnt) break;
            }

            if (ncnt && kcnt) {
                int k = 0;
                for (int n = 1; n <= m + 1; ++n) {
                    npgn[ncount][n] = tmp[n - 1];
//                    printf("tmp %d,",tmp[n-1]);
                    if (tmp[n - 1] != 0) k++;
                }
//                printf("k %d\n",k);
                npgn[ncount][0] = k;
//                printf("ncount %d\n",ncount);
                ncount++;
            }
        }

    }
    return w;
//    free(temp);
}



int edglink(int *adj,int v,int m){
    int cnt = 0;
    for(int i=1;i<adj_list[v][0]+1;++i){
        if(findVal(adj,adj_list[v][i],m)) cnt++;
    }
    if(cnt == m) return 1;
    return 0;
}

////Iterative recursion
void findkPolygonCount() {
    printf("begin\n");
//#pragma omp parallel for
    for (int i = 1; i < trow; ++i) {

        int adj[MAX_MX],visited[MAX_ADJ];
        adj[0] = i;
        visited[i] = 1;
        dfs(i, adj, 1,visited,0);
    }
}

void printkPolygon() {
    printf("k-truss count %d :\n", kcount);

    printf("%d-truss subgraph :\n", tk);

    ///(1) Output mutually exclusive k-truss
    int i = 1;
    for (int j = 0; j < kcount; ++j) {
        int w = 0;
        for (int jj = 0; jj < tk; ++jj) {
            if(adjUsed[kpgn[j][jj]]){
                w = 1;
                break;
            }
        }
        if(!w){
            printf("\t%d : ", i++);
            for (int jj = 0; jj < tk; ++jj) {
                adjUsed[kpgn[j][jj]] = 1;
                kusdpgn[kused][jj] = kpgn[j][jj];
                printf("%d,", kpgn[j][jj]);
            }
            kused++;
            printf("\n");
        }
    }
    tk--;
    ///(2) Search for relationships between the remaining vertices
    int unadjUsed[MAX_ADJ] ={0};
    int ud =0;
    ///Store all unused vertices
    for (int ii = 1; ii < trow; ++ii) {
        if(!adjUsed[ii]) unadjUsed[ud++] = ii;
    }
    for(int ii=0;ii<ud;++ii){
        int j = unadjUsed[ii];
        int adj[MAX_MX],visited[MAX_ADJ];
        adj[0] = j;
        visited[j] = 1;
        dfs(j, adj, 1,visited,1);
    }

    printf("n-truss count : %d \n",ncount);

    int nct = 0;
    i=-1;
    int twoT = -1;

    if(ncount == 0){
        i = 1;
    }
    twoT = ncount;
    for(int ii=0;ii<ud;++ii){
        int j = unadjUsed[ii];
        for(int jj=1;jj<adj_list[j][0];++jj){
            int v = adj_list[j][jj];
            if(!adjUsed[v]) continue;

            ///de-repeat -- n-truss
            int cnt =0;
            int tmp[2];
            if(j>v){
                tmp[0]=v;
                tmp[1]=j;
            }else{
                tmp[0]=j;
                tmp[1]=v;
            }
            for(int k = 0;k<ncount;++k){
                if(matchIncludeArrs(npgn[k],tmp,1,1)) {
                    cnt = 1;
                    break;
                }
            }
            if(cnt) continue;
            ///de-repeat- kused-truss
            for(int k = 0;k<kused;++k){
                if(matchIncludeArrs(kusdpgn[k],tmp,1,0)){
                    cnt =1;
                    break;
                }
            }
            if(cnt) continue;
            npgn[ncount][0]=2;
            if(j>v) {
                npgn[ncount][1] = v;
                npgn[ncount][2] = j;
            }else{
                npgn[ncount][1] = j;
                npgn[ncount][2] = v;
            }
            ncount++;

        }
    }

    ///Use unoutput edges between points
    for(int ii=0;ii<kcount;++ii){
        for(int j=0;j<tk-1;++j){
            int v1 = kpgn[ii][j];
            for(int jj=j+1;jj<tk;++jj){
                int v2 = kpgn[ii][jj];
                int isused = 0;
                ///Match (whether this edge is used)
                for(int k=0;k<kused;++k) {
                    int cnt = 0;
                    for(int kk=0;kk<tk;++kk){
                        if(kusdpgn[k][kk] == v1) cnt++;
                        if(kusdpgn[k][kk] == v2) cnt++;
                        if(cnt == 2) break;
                    }
                    ///the dege is used
                    if(cnt == 2) {
                        isused = 1;
                        break;
                    };
                }
                if(!isused){
                    //de-repeat
                    int cnt =0;
                    int tmp[2]={0};
                    if(v1>v2){
                        tmp[0]=v2;
                        tmp[1]=v1;
                    }else{
                        tmp[0]=v1;
                        tmp[1]=v2;
                    }
                    for(int k = 0;k<ncount;++k){
                        if(matchIncludeArrs(npgn[k],tmp,1,1)) {
                            cnt = 1;
                            break;
                        }
                    }
                    if(cnt) continue;
                    //de-repeat -- kused-truss
                    for(int k = 0;k<kused;++k){
                        if(matchIncludeArrs(kusdpgn[k],tmp,1,0)){
                            cnt =1;
                            break;
                        }
                    }
                    if(cnt) continue;
                    npgn[ncount][0]=2;
                    if(v1>v2) {
                        npgn[ncount][1] = v2;
                        npgn[ncount][2] = v1;
                    }else{
                        npgn[ncount][1] = v1;
                        npgn[ncount][2] = v2;
                    }
                    ncount++;
                }
            }
        }
    }
    ///output 2-truss
    quickSort(npgn,0,ncount-1);
    for (int j = 0; j < ncount; ++j) {
        if(nct != npgn[j][0]) {
            printf("%d-truss subgraph :\n",npgn[j][0] );
            nct = npgn[j][0];
            if(twoT==-1 && nct == 2)  twoT = j;
            i=1;
        }
        printf("\t%d : ",i++);
        for (int jj = 1; jj <= npgn[j][0]; ++jj) {

            printf("%d,", npgn[j][jj]);
        }
        printf("\n");
    }

}
void run(int triangle_k,int r,int sy,int nz,int **triples,int len,int ompn) {

    printf("run !!!\n");
    double starttime = omp_get_wtime();
    isSyc = sy;
    trow = r;
    omp_n = ompn;///thread count
    tk = triangle_k;
    kcount = 0,ncount=0;
    adj_list = (int **) malloc((10+r)*sizeof (int*));

    initAdjList(triples,len);

    findkPolygonCount();

    printkPolygon();
    double endtime = omp_get_wtime();
    printf("\ntime :%lfs\n",endtime-starttime);
    omp_destroy_lock(&lock);
    free(adj_list);
}