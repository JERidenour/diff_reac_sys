#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void printarr(int **data, int n, char *str);
int **allocarray(int n);

int main(int argc, char **argv) {

    // arrays
    int **bigarray = NULL;

    /* array sizes */
    const int bigsize =10;
    const int subsize =8;

    /* communications parameters */
    const int sender  =0;
    const int receiver=1;
    const int ourtag  =2;

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < receiver+1) {
        if (rank == 0)
            fprintf(stderr,"%s: Needs at least %d  processors.\n", argv[0], receiver+1);
        MPI_Finalize();
        return 1;
    }

    if (rank == sender) {
        bigarray = allocarray(bigsize);
        for (int i=0; i<bigsize; i++)
            for (int j=0; j<bigsize; j++)
                bigarray[i][j] = i*bigsize+j;


        printarr(bigarray, bigsize, " Sender: Big array ");

        MPI_Datatype sendSubarray;

        int starts[2] = {1,bigsize-1};          // E
        int subsizes[2]  = {subsize,1};         // vert stripe
        int bigsizes[2]  = {bigsize, bigsize};
        MPI_Type_create_subarray(2, bigsizes, subsizes, starts,
                                 MPI_ORDER_C, MPI_INT, &sendSubarray);
        MPI_Type_commit(&sendSubarray);

        MPI_Send(&(bigarray[0][0]), 1, sendSubarray, receiver, ourtag, MPI_COMM_WORLD);
        MPI_Type_free(&sendSubarray);

        free(bigarray[0]);
        free(bigarray);

    } else if (rank == receiver) {

        bigarray = allocarray(bigsize);
        for (int i=0; i<bigsize; i++)
            for (int j=0; j<bigsize; j++)
                bigarray[i][j] = i*bigsize+j;


        //printarr(bigarray, bigsize, " Receiver: Big array ");

        MPI_Datatype recSubarray;

        int starts[2] = {1,0};          // W
        int subsizes[2]  = {subsize,1};         // vert stripe
        int bigsizes[2]  = {bigsize, bigsize};
        MPI_Type_create_subarray(2, bigsizes, subsizes, starts,
                                 MPI_ORDER_C, MPI_INT, &recSubarray);
        MPI_Type_commit(&recSubarray);

        // int **subarray = allocarray(subsize);

        // for (int i=0; i<subsize; i++)
        //     for (int j=0; j<subsize; j++)
        //         subarray[i][j] = 0;

        MPI_Recv(&(bigarray[0][0]), subsize*subsize, recSubarray, sender, ourtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // printarr(subarray, subsize, " Receiver: Subarray -- after receive");
        printarr(bigarray, bigsize, " Receiver: Bigarray -- after receive");

        MPI_Type_free(&recSubarray);

        // free(subarray[0]);
        // free(subarray);

        free(bigarray[0]);
        free(bigarray);
    }

    MPI_Finalize();
    return 0;
}

void printarr(int **data, int n, char *str) {    
    printf("-- %s --\n", str);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            printf("%3d ", data[i][j]);
        }
        printf("\n");
    }
}

int **allocarray(int n) {
    int *data = malloc(n*n*sizeof(int));
    int **arr = malloc(n*sizeof(int *));
    for (int i=0; i<n; i++)
        arr[i] = &(data[i*n]);

    return arr;
}