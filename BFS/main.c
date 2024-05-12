#include <omp.h>
#include "triangle.h"
#include "mmio.h"
#include "stdlib.h"

//input
int triangle_k = 4; /// k-truss
/*int rowPtr[] ={0,8,15,24,33,41,50,59,68,77,84};
int colInd[] ={2,3,4,5,6,7,8,9,
                2,3,4,5,6,7,8,
                0,1,3,4,5,6,7,8,9,
                0,1,2,4,5,6,7,8,9,
                0,1,2,3,5,6,7,8,
                0,1,2,3,4,6,7,8,9,
                0,1,2,3,4,5,7,8,9,
                0,1,2,3,4,5,6,8,9,
                0,1,2,3,4,5,6,7,9,
                0,2,3,5,6,7,8};*/
int *rowPtr ;// = { 0,5,12,19,26,33,38,42,47,52,56 };
int *colInd;//  { 1,2,3,4,5,  0,2,3,4,5,6,7,
//                          0,1,3,4,6,7,8, 0,1,2,4,
//                          5,8,9,  0,1,2,3,7,8,9,
//                          0,1,3,6,9, 1,2,5,7,
//                          1,2,4,6,8, 2,3,4,7,9,
//                          3,4,5,8 };


int row = 0;



int main(int argc, char** argv) {
    //char* filename = "matrix.mtx";
    char* filename = "matrix.mtx";
    if (argc > 1) {
        filename = argv[1];
        if(argc >2 ){
            triangle_k = atoi(argv[2]);
        }
    }
    else {
//        printf("file not found!\n");
//        exit(0);
        filename = "test.mtx";
//        filename = "testDemo2.mtx";
//        filename = "testDemo3.mtx";
    }
    printf("filename: %s\n", filename);
    // read matrix from mtx file
    //printf("read matrix from mtx file.\n");
    int ret_code;
    int m, n,nnzA_mtx_report;//rows,cols,not Zero number
    MM_typecode matcode;
    FILE* f;
    int *row_ind;
    int **triples;
    int **edgeVisited;////Used of edges
    int isInteger = 0, isReal = 0, isPattern = 0, isSymmetric = 0;/


    int maxThreads = omp_get_max_threads();

    int omp_n = maxThreads>>1;
    printf("max threads : %d \n",maxThreads);

    omp_set_num_threads(omp_n);


//    int omp_size = ;


    // load matrix
    if ((f = fopen(filename, "r")) == NULL)
        return -1;

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n" );
        return -2;
    }

    if (mm_is_complex(matcode))
    {
        printf("Sorry, data type 'COMPLEX' is not supported.\n");
        return -3;
    }

    if (mm_is_pattern(matcode)) { isPattern = 1; /*cout << "type = Pattern" << endl;*/ }
    if (mm_is_real(matcode)) { isReal = 1; /*cout << "type = real" << endl;*/ }
    if (mm_is_integer(matcode)) { isInteger = 1; /*cout << "type = integer" << endl;*/ }

    /* find out size of sparse matrix .... */
    ret_code = mm_read_mtx_crd_size(f, &m, &n, &nnzA_mtx_report);
    if (ret_code != 0)
        return -4;
    printf("%d,%d,%d\n",m,n,nnzA_mtx_report);


    if (mm_is_symmetric(matcode) || mm_is_hermitian(matcode))
    {
        isSymmetric = 1;
        triples = (int**) malloc(2*(10+nnzA_mtx_report)*sizeof (int*));
        //cout << "symmetric = true" << endl;
    }
    else
    {
        //cout << "symmetric = false" << endl;
    }

    if(n>row) row = n;
    if(m>row) row = m;
    row++;
    row_ind = (int*)malloc((5+row)*sizeof(int));
    rowPtr = (int*) malloc((5+row)*sizeof (int));
    if(isSymmetric)
        colInd = (int*) malloc(2*(5+nnzA_mtx_report)*sizeof (int));
    else
        colInd = (int*) malloc((5+nnzA_mtx_report)*sizeof (int));

    if(row_ind == NULL || rowPtr == NULL || colInd == NULL) {
        printf("error malloc!!!");
        exit(0);
    }
    memset(row_ind,0,(5+row)*sizeof (int));
    memset(rowPtr,0,(5+row)*sizeof (int));
    memset(colInd,0,(5+nnzA_mtx_report)*sizeof (int));



    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    // 0,5,12,19,26,33,38,42,47,52,56,
    // 1,2,3,4,5,0,2,3,4,5,6,7,0,1,3,4,6,7,8,0,1,2,4,5,8,9,0,1,2,3,7,8,9,0,1,3,6,9,1,2,5,7,1,2,4,6,8,2,3,4,7,9,3,4,5,8,

    int maxadj = 1;
    int j = 0;
    for (int i = 0; i < nnzA_mtx_report; i++)
    {
        int idxi, idxj;//row col
        double fval;
        int ival;

        if (isReal)
            fscanf(f, "%d %d %lg\n", &idxi, &idxj, &fval);
        else if (isInteger)//~U~U
        {
            fscanf(f, "%d %d %d\n", &idxi, &idxj, &ival);
            fval = ival;
        }
        else if (isPattern)
        {
            fscanf(f, "%d %d\n", &idxi, &idxj);
            fval = 1.0;
        }
        if(idxi == idxj) continue;

        triples[j] = (int *) malloc(3*sizeof (int));
        triples[j][0] = idxi;
        triples[j][1] = idxj;
        if(maxadj < idxi) maxadj = idxi;
        if(maxadj < idxj) maxadj = idxj;
        j++;
    }
    printf("maxadj == %d\n",maxadj);
    if (f != stdin)
        fclose(f);

    run(triangle_k,row,isSymmetric,nnzA_mtx_report,triples,j,omp_n);

    free(row_ind);
    free(rowPtr);
    free(colInd);

    return 0;
}