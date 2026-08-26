#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "time.h"
// must for c language
// program to calculate the bond percolation threshold for a system of two square lattices (with bond probability p)
// To find pc, vary the probabilities until the upper cumulative of the bins, multiplied by s^(tau-2), is a constant for large s,
// where s = 2^n.

void    randinit(int seed);            // about random number
#define SEED 123
#define M 16383                        // used for RZ random number generator
#define NewRandomInteger (++nd,ra[nd&M] = ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])
int 	ra[M+1], nd;

#define HEIGHT 16384    // must be a power of 2 for periodic b.c., system size: HEIGHT*WIDTH*DEPTH
#define WIDTH HEIGHT
#define H (HEIGHT-1)   // these allow for periodic b.c. by using "&"
#define W (WIDTH-1)    // not a array, is a value 1023, 1111111111
// 1(512), 1(256), 1(128), 1(64), 1(32), 1(16), 1(8), 1(4), 1(2), 1(1)
int     lat[HEIGHT][WIDTH];      // the lattice

#define S 65535   // used for list length power of 2, minus 1, 1111111111111111
// 1(32768), 1(16384), 1(8192), 1(4096), 1(2048), 1(1024), 1(512), 1(256), 1(128), 1(64), 1(32), 1(16), 1(8), 1(4), 1(2), 1(1)
#define GetFromStack(X,Y) {X = xlist[gptr & S]; Y = ylist[gptr & S]; ++gptr;}  // get lsit from stack
#define PutOnStack(X,Y)	{xlist[pptr & S] = X; ylist[pptr & S] = Y; ++pptr;}  // put the list into stack
int 	xlist[S+1], ylist[S+1];       // the stack is used to temporary storage of data: xlist[S+1], ylist[S+1], zlist[S+1]

#define RUNSMAX  2147483647     // number of runs, -2147483648(-2^31),...,-1,0,1,2,...,2147483647(2^31-1), related to the bits of computer 
#define PROB  0.258        // bond probability in the plane
#define MAX 65536//1048576/2             // size cutoff -- make a power of 2; maximun cutoff size of the cluster, 2^20, related to bin[]     
#define PRINTFREQ 65535          // power of 2, minus 1, 111111111

int main(void)         // main function
{
	int dx[14]={-2,-2,-2,  -1,-1,-1,   0,0,   1,1,1,   2,2,2};
	int dy[14]={-1, 0, 1,  -1, 0, 1,  -1,1,  -1,0,1,  -1,0,1};       // direction vectors, x, y and z
	int x, y, xo, yo, xp, yp, dir, prob, gptr, pptr, runs, i, bin[32];
    long  nocctot;

    prob=(int) (2147483648.0 * PROB); // convert probabilities to integers, random number range [1, 2147483648]
	randinit(SEED);  // call function
	srand((unsigned int)time(NULL));

	for (x = 0; x < WIDTH; ++x)
        for (y = 0; y < HEIGHT; ++y)
            lat[x][y] = 0;         // clear the lattice

	for (i = 0; i < 32; ++i) 
		bin[i] = 0;                       // clear the bin[], bin[] is used to record the number of bins

    nocctot = 0;   // clear nocctot, to keep track of the total number of occupied sites generated
    
	for (runs = 1; runs < RUNSMAX; ++runs)
	{
        xo = yo = 0;   // start at origin point
		// don't need to clear the lattice each time using "runs" this way
        lat[xo][yo] = runs;  // make first site occupied

        gptr = pptr = 0;      // stack pointers 
        PutOnStack(xo,yo)  // put the information of the origin point to stack
		// #define PutOnStack(X,Y,Z)	{xlist[pptr & S] = X; ylist[pptr & S] = Y; zlist[pptr & S] = Z; ++pptr;}
		
		

        do
        {
            GetFromStack(x,y)  // stack is really a "queue", get the information of list from the stack 
			// #define GetFromStack(X,Y,Z) {X = xlist[gptr & S]; Y = ylist[gptr & S]; Z = zlist[gptr & S]; ++gptr;}
            
			int ret = rand()%10 +1;
		//	printf("%d\n", ret);
						

			if (ret < 6)      // horizontal
			{
				for (dir = 0; dir < 14; ++dir)
				{
					xp = x + dx[dir];
					yp = y + dy[dir];
					if (lat[xp & W][yp & H] < runs)
					if (NewRandomInteger < prob)  // then make occupied, two conditions: runs and prbobility p
					{
						lat[xp & W][yp & H] = runs;
						PutOnStack(xp,yp)
						// #define PutOnStack(X,Y,Z)	{xlist[pptr & S] = X; ylist[pptr & S] = Y; zlist[pptr & S] = Z; ++pptr;}
					}
				}
			}
			else    // vertical
			{
				for (dir = 0; dir < 14; ++dir)
				{
					xp = x + dy[dir];
					yp = y + dx[dir];
					if (lat[xp & W][yp & H] < runs)
					if (NewRandomInteger < prob)  // then make occupied, two conditions: runs and prbobility p
					{
						lat[xp & W][yp & H] = runs;
						PutOnStack(xp,yp)
						// #define PutOnStack(X,Y,Z)	{xlist[pptr & S] = X; ylist[pptr & S] = Y; zlist[pptr & S] = Z; ++pptr;}
					}
				}
			}



        }
        while ((gptr != pptr) && (pptr < MAX));    // when gptr=pptr, the cluster stop growing
        
        ++bin[(int)(log(pptr)/log(1.9999999))];    // pptr is the cluster size, put in log2 bins
        nocctot += pptr;                           // to keep track of the total number of occupied sites generated.
		if ((runs & PRINTFREQ) == 0)               // output the results every 1024 runs
		{
            FILE *fpWrite=fopen("data.txt","a");   // open a file,	
			//FILE *fpWrite=fopen("data.txt","a"); // if you don't want to overwrite the original data
			printf("%12.8f%13d%13ld\n", PROB, runs, nocctot);
			fprintf(fpWrite, "%12.8f%13d%13ld\n", PROB, runs, nocctot);   // write the data to a file
            for (i = 0; i < 32; ++i)
				if (bin[i])
				{
					printf("%10d%18.10e\n", i, bin[i]*1.0/runs);
					fprintf(fpWrite,"%10d%18.10e\n",i, bin[i]*1.0/runs);  // write the data to a file
				}
			fclose (fpWrite);
		} // if runs	
    } // for runs
	
}         // end of the main function

void randinit(int seed)  //seed the XOR random number generator with a congruential generator.
{	
	double a, ee = -1 + 1/2147483648.0;
	int i;
	extern int nd, ra[M+1];
	
	a = seed/2147483648.0;
	for (nd = 0; nd <= M; nd++)
	{
		a *= 16807;
		a += ee * (int)(a);
		if (a >= 1) a += ee;
		ra[nd] = (int) (2147483648.0 * a);
	}
	nd = M;
	for(i = 0; i<100001; i++) NewRandomInteger;  // warmup generator
    for (i = 0; i < 10; ++i) printf("rng %20d\n",ra[i]); // print out to show it's working (not all zeros!!!)
}