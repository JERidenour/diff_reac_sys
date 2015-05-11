if (!(meshRank%2))
	{
		err = MPI_Send(&(localData[0][0]), n, borders[NORTH], RANK_NORTH, 666, MPI_COMM_WORLD);
		// err = MPI_Send(&(localData[0][0]), n, borders[EAST], RANK_EAST, 666, MPI_COMM_WORLD);
		// err = MPI_Send(&(localData[0][0]), n, borders[SOUTH], RANK_SOUTH, 666, MPI_COMM_WORLD);
		// err = MPI_Send(&(localData[0][0]), n, borders[WEST], RANK_WEST, 666, MPI_COMM_WORLD);

		err = MPI_Recv(&(localData[0][0]), n, ghosts[NORTH], RANK_NORTH, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// err = MPI_Recv(&(localData[0][0]), n, ghosts[EAST], RANK_EAST, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// err = MPI_Recv(&(localData[0][0]), n, ghosts[SOUTH], RANK_SOUTH, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// err = MPI_Recv(&(localData[0][0]), n, ghosts[WEST], RANK_WEST, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else
	{
		err = MPI_Recv(&(localData[0][0]), n, ghosts[SOUTH], RANK_SOUTH, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// err = MPI_Recv(&(localData[0][0]), n, ghosts[WEST], RANK_WEST, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// err = MPI_Recv(&(localData[0][0]), n, ghosts[NORTH], RANK_NORTH, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// err = MPI_Recv(&(localData[0][0]), n, ghosts[EAST], RANK_EAST, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		err = MPI_Send(&(localData[0][0]), n, borders[SOUTH], RANK_SOUTH, 666, MPI_COMM_WORLD);
		// err = MPI_Send(&(localData[0][0]), n, borders[WEST], RANK_WEST, 666, MPI_COMM_WORLD);
		// err = MPI_Send(&(localData[0][0]), n, borders[NORTH], RANK_NORTH, 666, MPI_COMM_WORLD);
		// err = MPI_Send(&(localData[0][0]), n, borders[EAST], RANK_EAST, 666, MPI_COMM_WORLD);
	}

	if(meshRank == RANK_CENTER)
	{
		MPI_Sendrecv(	&(localData[0][0]),
					n,
					borders[NORTH],
					RANK_NORTH,
					666,
					&(localData[0][0]),
					n,
					ghosts[NORTH],
					RANK_NORTH,
					777,
					processMesh,
					MPI_STATUS_IGNORE
					);
	}
	else if(meshRank == RANK_NORTH)
	{
		MPI_Sendrecv(	&(localData[0][0]),
					n,
					borders[SOUTH],
					RANK_SOUTH,
					666,
					&(localData[0][0]),
					n,
					ghosts[SOUTH],
					RANK_SOUTH,
					777,
					processMesh,
					MPI_STATUS_IGNORE
					);
	}

	//printf("MeshRank %d %d communicating with %d.\n", meshRank, RANK_CENTER, RANK_NORTH );
		// MPI_Sendrecv(	&(localData[0][0]),
		//                 n,
		//                 borders[NORTH],
		//                 RANK_NORTH,
		//                 666,
		//                 &(localData[0][0]),
		//                 n,
		//                 ghosts[NORTH],
		//                 RANK_NORTH,
		//                 666,
		//                 processMesh,
		//                 MPI_STATUS_IGNORE
		//             );


	//printf("MeshRank %d %d communicating with %d.\n", meshRank, RANK_CENTER, RANK_SOUTH );
		// MPI_Sendrecv(	&(localData[0][0]),
		//                 n,
		//                 borders[SOUTH],
		//                 RANK_SOUTH,
		//                 666,
		//                 &(localData[0][0]),
		//                 n,
		//                 ghosts[SOUTH],
		//                 RANK_SOUTH,
		//                 666,
		//                 processMesh,
		//                 MPI_STATUS_IGNORE
		//             );
		//printf("Communication step done.\n");


		// if (meshRank == 0)
	// {
	// 	printf("Rank: %d, meshRank: %d, i, j: %d,%d, hood: %d %d %d %d %d %d %d %d %d\n",
	// 	       rank, meshRank, meshCoords[0], meshCoords[1],
	// 	       hood[0][0], hood[0][1], hood[0][2],
	// 	       hood[1][0], hood[1][1], hood[1][2],
	// 	       hood[2][0], hood[2][1], hood[2][2]
	// 	      );
	// }


		/*int ghostStarts[] = {	0, 1,
	                        1, n + 1,		//	east starts
	                        n + 1, 1,		// 	south starts
	                        1, 0,			//	west starts
	                        0, n + 1,		// 	NE starts
	                        n + 1, n + 1, 	// 	SE starts
	                        n + 1, 0,		// 	SW starts
	                        0, 0			// 	NW starts
	                    };*/


		/*int borderStarts[] = {	1, 1,			// 	north starts
	                        1, n,			//	east starts
	                        n, 1,			// 	south starts
	                        1, 1,			//	west starts
	                        1, n,			// 	NE starts
	                        n, n,		 	// 	SE starts
	                        n, 1,			// 	SW starts
	                        1, 1,			// 	NW starts
	                     };*/                        


// err = MPI_Send(&(localData[0][1]), 1, border[NORTH], 2, 666, MPI_COMM_WORLD);
		/*err = MPI_Send(&(DATA_BORDER_NORTH), 1, border[NORTH], neighbour[NORTH], 666, MPI_COMM_WORLD);
		err = MPI_Send(&(DATA_BORDER_EAST), n, border[EAST], neighbour[NORTH], 666, MPI_COMM_WORLD);
		err = MPI_Send(&(DATA_BORDER_SOUTH), 1, border[SOUTH], neighbour[NORTH], 666, MPI_COMM_WORLD);
		err = MPI_Send(&(DATA_BORDER_WEST), n, border[WEST], neighbour[NORTH], 666, MPI_COMM_WORLD);

		err = MPI_Send(&(DATA_BORDER_NORTHEAST), 1, MPI_DOUBLE, neighbour[NORTH], 666, MPI_COMM_WORLD);
		err = MPI_Send(&(DATA_BORDER_SOUTHEAST), 1, MPI_DOUBLE, neighbour[NORTH], 666, MPI_COMM_WORLD);
		err = MPI_Send(&(DATA_BORDER_SOUTHWEST), 1, MPI_DOUBLE, neighbour[NORTH], 666, MPI_COMM_WORLD);
		err = MPI_Send(&(DATA_BORDER_NORTHWEST), 1, MPI_DOUBLE, neighbour[NORTH], 666, MPI_COMM_WORLD);
		*/

/*err = MPI_Recv(&(DATA_GHOST_SOUTH), 1, ghost[SOUTH], neighbour[SOUTH], 666, MPI_COMM_WORLD, &status);
		err = MPI_Recv(&(DATA_GHOST_WEST), n, ghost[WEST], neighbour[SOUTH], 666, MPI_COMM_WORLD, &status);
		err = MPI_Recv(&(DATA_GHOST_NORTH), 1, ghost[NORTH], neighbour[SOUTH], 666, MPI_COMM_WORLD, &status);
		err = MPI_Recv(&(DATA_GHOST_EAST), n, ghost[WEST], neighbour[SOUTH], 666, MPI_COMM_WORLD, &status);
		
		err = MPI_Recv(&(DATA_GHOST_SOUTHWEST), n, MPI_DOUBLE, neighbour[SOUTH], 666, MPI_COMM_WORLD, &status);
		err = MPI_Recv(&(DATA_GHOST_NORTHWEST), n, MPI_DOUBLE, neighbour[SOUTH], 666, MPI_COMM_WORLD, &status);
		err = MPI_Recv(&(DATA_GHOST_NORTHEAST), n, MPI_DOUBLE, neighbour[SOUTH], 666, MPI_COMM_WORLD, &status);
		err = MPI_Recv(&(DATA_GHOST_SOUTHEAST), n, MPI_DOUBLE, neighbour[SOUTH], 666, MPI_COMM_WORLD, &status);
		*/



		if (meshRank == 0)
	{
		//printf("rank %d sending to %d.\n", rank, 2 );
		
		err = MPI_Sendrecv(	&(DATA_BORDER_NORTH),
					1,
					border[NORTH],
					neighbour[NORTH],
					666,
					&(DATA_GHOST_NORTH),
					1,
					ghost[NORTH],
					neighbour[NORTH],
					666,
					processMesh,
					&status
					);

		MPI_Get_count(&status, ghost[NORTH], &rcount);
		printf("0 received: %d from NORTH.\n", rcount);

		err = MPI_Sendrecv(	&(DATA_BORDER_EAST),
					n,
					border[EAST],
					neighbour[EAST],
					777,
					&(DATA_GHOST_EAST),
					n,
					ghost[EAST],
					neighbour[EAST],
					777,
					processMesh,
					&status
					);

		MPI_Get_count(&status, ghost[EAST], &rcount);
		printf("0 received: %d from EAST.\n", rcount);

		// printf("Send step done, err: %d.\n", err);

		// printf("cerr %d terr %d\n", cerr, terr);
	}
	if (meshRank == 2)
	{

		//printf("MeshRank %d receiving from %d.\n", rank, 0 );

		err = MPI_Sendrecv(	&(DATA_BORDER_SOUTH),
					1,
					border[SOUTH],
					neighbour[SOUTH],
					666,
					&(DATA_GHOST_SOUTH),
					1,
					ghost[SOUTH],
					neighbour[SOUTH],
					666,
					processMesh,
					&status
					);

		MPI_Get_count(&status, ghost[SOUTH], &rcount);
		printf("2 received: %d from SOUTH.\n", rcount);

	}
	if(meshRank == 1)
	{
		err = MPI_Sendrecv(	&(DATA_BORDER_WEST),
					n,
					border[WEST],
					0,
					//neighbour[WEST],
					777,
					&(DATA_GHOST_WEST),
					n,
					ghost[WEST],
					0,
					//neighbour[WEST],
					777,
					processMesh,
					&status
					);

		MPI_Get_count(&status, ghost[WEST], &rcount);
		printf("2 received: %d from WEST.\n", rcount);
	}