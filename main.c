#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<errno.h>
#include<mpi.h>

int monte_carlo_sampling(double* monte_carlo_result, int NumOfEpoch, int rank);
int main(int argc, char* argv[])
{
	int rank, nproc, root_node;
	root_node = 0;// The root node
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);// Number of process
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);// The index of current process
	int NumOfEpoch = 1e2;
	int BufferSize;
	BufferSize = NumOfEpoch;
	double* send_buff = (double*)malloc(BufferSize * sizeof(double));
	double* recv_buff = (double*)malloc(BufferSize * (double)nproc * sizeof(double));
	monte_carlo_sampling(send_buff, NumOfEpoch, rank);
	MPI_Gather(send_buff, BufferSize, MPI_DOUBLE, recv_buff, BufferSize, MPI_DOUBLE, root_node, MPI_COMM_WORLD);
	int index_epoch, index_proc;
	MPI_Barrier(MPI_COMM_WORLD);
	FILE* fp;
	errno_t err;
	if (rank == 0)
	{   
		err = fopen_s(&fp,"MCdata.dat", "w");
		for (index_epoch = 0; index_epoch < BufferSize; index_epoch++)
		{
			for (index_proc = 0; index_proc < nproc; index_proc++)
			{
				printf("%4.3lf\t", recv_buff[index_proc * BufferSize + index_epoch]);
				fprintf(fp,"%4.3lf\t", recv_buff[index_proc * BufferSize + index_epoch]);
			};
			printf("\n");
			fprintf(fp,"\n");
		};
		fclose(fp);
	};
	free(send_buff);
	free(recv_buff);
	MPI_Finalize();
	return 0;
}

int monte_carlo_sampling(double* monte_carlo_result, int NumOfEpoch, int rank)
{

	//Setting Up the MonteCarlo parameters
	int NumInEdge = 10;
	int NumOfWarm = 5000;
	double beta = 0.2;
	srand(rank + (int)time(0));
	int indexer;//iterator(i)
	int indexer_x, indexer_y;
	int x_index, y_index;
	int val_left, val_right, val_up, val_down;
	double Ene_old, Ene_new;
	double sum = 0.0;
	int** ising_zone;
	ising_zone = (int**)malloc(sizeof(int*) * NumInEdge);
	for (indexer = 0; indexer < NumInEdge; indexer++)//Creates an 2d array
	{
		ising_zone[indexer] = (int*)malloc(sizeof(int) * NumInEdge);
	};

	for (indexer_x = 0; indexer_x < NumInEdge; indexer_x++)
	{
		for (indexer_y = 0; indexer_y < NumInEdge; indexer_y++)
		{
			if (((double)rand() )/ RAND_MAX < 0.5)
			{
				ising_zone[indexer_x][indexer_y] = 1;
			}
			else
			{
				ising_zone[indexer_x][indexer_y] = -1;
			};
		};
	};
	for (indexer = 0; indexer < NumOfWarm; indexer++)
	{
		
		x_index = rand() % NumInEdge;
		y_index = rand() % NumInEdge;

		val_left = ising_zone[(x_index - 1 + NumInEdge) % NumInEdge][y_index];
		val_right = ising_zone[(x_index + 1 + NumInEdge) % NumInEdge][y_index];
		val_up = ising_zone[x_index][(y_index - 1 + NumInEdge) % NumInEdge];
		val_down = ising_zone[x_index][(y_index + 1 + NumInEdge) % NumInEdge];
		Ene_old = -ising_zone[x_index][y_index] * (val_left + val_right + val_up + val_down);
		Ene_new = -Ene_old;
		if (((double)rand() )/ RAND_MAX < exp(-beta * (Ene_new - Ene_old)))
		{
			ising_zone[x_index][y_index] = -ising_zone[x_index][y_index];
		};
	}

	for (indexer = 0; indexer < NumOfEpoch; indexer++)
	{
		x_index = rand() % NumInEdge;
		y_index = rand() % NumInEdge;

		val_left = ising_zone[(x_index - 1 + NumInEdge) % NumInEdge][y_index];
		val_right = ising_zone[(x_index + 1 + NumInEdge) % NumInEdge][y_index];
		val_up = ising_zone[x_index][(y_index - 1 + NumInEdge) % NumInEdge];
		val_down = ising_zone[x_index][(y_index + 1 + NumInEdge) % NumInEdge];
		Ene_old = -ising_zone[x_index][y_index] * (val_left + val_right + val_up + val_down);
		Ene_new = -Ene_old;
		if (((double)rand()) / RAND_MAX < exp(-beta * (Ene_new - Ene_old)))
		{
			ising_zone[x_index][y_index] = -ising_zone[x_index][y_index];
		};
		sum = 0.0;
		for (indexer_x = 0; indexer_x < NumInEdge; indexer_x++)
		{
			for (indexer_y = 0; indexer_y < NumInEdge; indexer_y++)
			{
				sum = sum + ising_zone[indexer_x][indexer_y];
			};
		};
		monte_carlo_result[indexer] = sum / (NumInEdge * NumInEdge);
	}

	return 0;

}