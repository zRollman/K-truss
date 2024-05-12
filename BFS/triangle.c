/**
 * Requirements .
 * 1. Find all k-edge types (mutually exclusive).
 * 2. Find all m(2<=m<k) edge types (de-duplication).
 * 3. Do not contain duplicates (e.g. 4-truss has 1,2,3,4; 3-truss can not contain the related ones inside, 1,2,3; 1,2,4, etc., and so on)
 * **/

/**
 * Methods.
 * 1. All possible permutations (2~k), choose from m (up to k vertices).
 * 2. Determine if all permutations can form a polygon (all points must be connected).
 * 3. Sort and de-weight the m-edge types.
 * 4. De-emphasise k-edge exclusivity
 */

#include "triangle.h"
#include "mmio.h"
#include "omp.h"
typedef struct bool bool;

omp_lock_t lock;

void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

int maxValue(int a,int b) {
    return a>b?a:b;
}
void swapArr(int *a,int *b,int len){

    for (int k = 0; k <= len; k++) {

        swap(&a[k], &b[k]);
    }
}

int singleCmp(const void* a, const void* b)
{
    return  *(int*)a - *(int*)b ;
}

int doubleCmp(const void* a, const void* b){
    return((int*)a)[0]-((int*)b)[0];
}

int charCmp(const void* p1, const void* p2)
{
    //先用（char**）强转，再用*解引用一层
    return strcmp(*(char**)p1,*(char**)p2);
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
//            printf("i %d,j %d\n",i,j);
/*            for (int k = 0; k < tk; k++) {
                swap(&arr[i][k], &arr[j][k]);
            }*/
            swapArr(arr[i],arr[j], maxValue(arr[i][0],arr[j][0]));
            i++;
            j--;
        }
    }

    quickSort(arr, left, j);
    quickSort(arr, i, right);
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

        //repeat point
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

            for(int j=1;j<=adj_list[i][0];++j){
                if(i!=adj_list[i][j]){
                    int v = adj_list[i][j];

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

    for (int i = 1; i < trow; ++i) {
        printf("%d,%d : ", i,adj_list[i][0]);


        for (int j = 1; j < adj_list[i][0] + 1; ++j) {
            printf("%d-", adj_list[i][j]);
        }
        printf("\n");
    }

}
///Determine whether the vertices are connected
int matchAdjs(int arr[],int k,int begin){
    ///arr[i], arr[i+1~k) deterministic connectivity
    for(int i=begin;i<k;++i){
        int n = adj_list[arr[i]][0];
        int *temp = (int *) malloc((n+2)*sizeof (int));
        temp[0] = arr[i];
        for(int j = 1;j<n+1;++j){
            temp[j] = adj_list[arr[i]][j];
        }
        qsort(temp, n+1, sizeof(int), singleCmp);
        int j=1,jj = i+1;
        while(j<n+1){
            ///Both are sequential
            if(arr[jj] ==temp[j]){
                j++;
                jj++;
            }else j++;
            if(jj>=k) break;
        }
        free(temp);
        if(jj < k){
            ///arr[i] 与 arr[i+1 ~ k ) not connectivity
            return 0;
        }

    }
    return 1;
}
int matchArrs(int arr[]){
    for(int i=0;i<pcount;++i){
        if(arr[0]<=PGN[i][0]){
            ///big include small
            int j = 1;
            int jj = 1;
            while(j<=PGN[i][0]){
                if(arr[jj] == PGN[i][j]){
                    j++;
                    jj++;
                }else j++;
            }
            if(jj > arr[0]){
                ///included
                return 1;
            }

        }
    }
    return 0;
}
int IsRepeatKtruss(int *temp,int count, int k){

    for(int i=0;i<count;++i){
        int flag = 0;
        for(int j=0;j<k;++j){
            if(kpgn[i][j]!=temp[j]) {
                flag = 1;
                break;
            }
        }
        if(!flag){
            ///Find the same
            return 1;
        }
    }
    return 0;
}
int IsRepeatNtruss(int *temp,int count, int m){

    for(int i=0;i<count;++i){

        if(npgn[i][0] == m){
            int flag = 0;
            for(int j=1;j<=m;++j){
                if(npgn[i][j]!=temp[j-1]){
                    flag = 1;
                    break;
                }
            }
            if(!flag){
                ///找到相同的
                return 1;
            }
        }

    }
    return 0;
}
//// Select k numbers from n numbers
//// Tree - Branching
void BFS(int arr[], int n, int k) {
    /// Initialising the queue
    State *queue = (State *)malloc(MAX_ADJ * MAX_ADJ * sizeof(State));
    if (queue == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    int front = 0, rear = 0;

    /// initial state
    State initial_state = {0, {0}, 0};
    queue[rear++] = initial_state;

    while (front < rear) {
        /// out queue
        State current_state = queue[front++];

        /// If k numbers have been selected, output the result
        if (current_state.selected_count == k ) {
            /// k-truss;
            int temp[MAX_MX] = {0};
            for (int i = 0; i < k; ++i) {
                temp[i] = arr[current_state.selected_so_far[i]];
            }
            ///Judging k be connected
            if(matchAdjs(temp,k,0)){
                ///linked
                ///de-repeat

                if(!IsRepeatKtruss(temp,kcount,k)){
                /// new
                    for (int i = 0; i < k; ++i) {
                        kpgn[kcount][i] = temp[i];
                    }
                    kcount++;
                }
            }

            continue;
        }
        else if(current_state.selected_count < k && current_state.selected_count >=2)
        {
            /// m-truss
            int temp[MAX_MX] = {0};
            int m = current_state.selected_count;
            for (int i = 0; i < m; ++i) {
                temp[i] = arr[current_state.selected_so_far[i]];
            }


            if(matchAdjs(temp,m,0)){

                if(!IsRepeatNtruss(temp,ncount,m)){

                    npgn[ncount][0] = m;
                    for (int i = 1; i <= m; ++i) {

                        npgn[ncount][i] = temp[i-1];
                    }
                    ncount++;

                }

            }
        }

        /// Add the next element to the queue
        for (int next_index = current_state.index; next_index < n; ++next_index) {
            State next_state = current_state;
            next_state.selected_so_far[next_state.selected_count++] = next_index;
            next_state.index = next_index + 1;

            /// join queue
            queue[rear++] = next_state;
        }
    }
    free(queue);
}

void findkPolygonCount() {
    printf("begin\n");
    #pragma omp parallel for
    for (int i = 1; i < trow; ++i) {

        int arr[MAX_ADJ+5] = {0};
        arr[0] = i;

        for(int j=1;j<adj_list[i][0]+1;++j){
            arr[j] = adj_list[i][j];

        }
        qsort(arr, adj_list[i][0]+1, sizeof(int), singleCmp);
        BFS(arr,adj_list[i][0]+1,tk);
    }
    quickSort(npgn,0,ncount-1);

}

void DelIncludePgn(){
    ///PGN = part(kpgn,npgn)
    int vis[MAX_ADJ] = {0};
    int trussNumber[MAX_MX] = {0};
    for(int i=0;i<kcount;++i){
        ///互斥
        if(pcount == 0){
            for(int j=0;j<tk;++j){
                vis[kpgn[i][j]] = 1;
                PGN[pcount][j+1] = kpgn[i][j];
            }
            trussNumber[tk]++;
            PGN[pcount++][0] = tk;
        }else{
            int flag =0 ;
            for(int j=0;j<tk;++j){
                if(vis[kpgn[i][j]]) {
                    flag = 1;
                    break;
                }
            }
            if(!flag){
                for(int j=0;j<tk;++j){
                    vis[kpgn[i][j]] = 1;
                    PGN[pcount][j+1] = kpgn[i][j];
                }
                trussNumber[tk]++;
                PGN[pcount++][0] = tk;
            }
        }
    }

    ///De-include repeat
    for(int i=ncount-1;i>=0;--i){
        ///back-to-front search --- PGN,npgn
        if(!matchArrs(npgn[i])){
            for(int j =0 ;j<=npgn[i][0];++j){
                PGN[pcount][j] = npgn[i][j];
            }
            trussNumber[npgn[i][0]]++;
            pcount++;
        }
    }
    printf("PGN pcount == %d\n",pcount);
    printf("%d-truss is %d :\n",tk,trussNumber[tk]);
    int truss = tk;
    for(int i=0;i<pcount;++i){
        int m = PGN[i][0];
        if(m != truss) {
            truss = m;
            printf("%d-truss is %d :\n",truss,trussNumber[truss]);
        }
        for(int j=1;j<=m;++j){
            printf("%d,",PGN[i][j]);
        }
        printf("\n");
    }
}
void run(int triangle_k,int r,int sy,int nz,int **triples,int len,int ompn) {
    printf("run !!!\n");
    double starttime = omp_get_wtime();
    isSyc = sy;
    trow = r;
    omp_n = ompn;
    tk = triangle_k;
    kcount = 0,ncount=0;
    adj_list = (int **) malloc((10+r)*sizeof (int*));
    initAdjList(triples,len);

    findkPolygonCount();

    DelIncludePgn();

    double endtime = omp_get_wtime();
    printf("\ntime :%lfs\n",endtime-starttime);
    omp_destroy_lock(&lock);
    free(adj_list);
}
