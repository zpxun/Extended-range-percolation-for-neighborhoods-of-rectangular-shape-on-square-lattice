#include "stdlib.h"
#include "stdio.h"
#include "math.h"

// Model II of Xun, Hao, Ziff 2026 with connections in both directions

void randinit(long seed);

#define HEIGHT 8192 //power of two
#define WIDTH HEIGHT
#define L 1 //neighborhood is (2L+1)(4L+1)
#define RUNSMAX  2147483647
#define PRINTFREQ  1023
#define PROB 0.22203 //0.21
#define H (HEIGHT-1)
#define W (WIDTH-1)
#define OUTFILE "output"
#define S 65535
#define SEED 123
#define DIRMAX ((2*L+1)*(4*L+1)-1)
#define M  16383
#define GetFromStack(X,Y,Z) {X = xlist[gptr & S]; Y = ylist[gptr & S], Z = zlist[gptr & S]; ++gptr;}
#define PutOnStack(X,Y,Z)    {xlist[pptr & S]=X; ylist[pptr & S]=Y; zlist[pptr & S]=lat[X&W][Y&W]; ++pptr;}
#define NewRandomInteger (++nd,ra[nd&M] = ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])

long     ra[M+1], nd;
long    lat[WIDTH][HEIGHT];
long     xlist[S+1], ylist[S+1], zlist[S+1];
FILE    *fp1;

int main(void)
{
    long    dx[4] = {1, 0, -1, 0}, dy[4] = {0, 1, 0, -1};
    long     x, y, z, xo, yo, xp, yp, dir,  prob, gptr, pptr, nocc, runs, runs2, i, j, chemdist, hit, nocctot, nc, nctot, big;
    long  bin[32], sbin[32];
    long  delx[DIRMAX],dely[DIRMAX];
    int  index;

    FILE   *fp1;

    prob = (long) (2147483648.0 * PROB);
    randinit(SEED);
    
    index = 0;
   for (x = -2*L; x <= 2*L; ++x)
    for (y = -L; y <= L; ++y)
    if ((x != 0) || (y != 0))
         
    {
        delx[index] = x;
        dely[index] = y;
        printf(" %d %ld %ld\n", index, delx[index], dely[index]);
        ++index;
    }
    printf("index %d  %d \n", index, DIRMAX);
    
  
    //clear the lattice
    
  
    for (i = 0; i < 32; ++i)
        bin[i] = sbin[i] = 0;
    
    nctot = 0;
        
    for (runs = 1; runs <= RUNSMAX; runs+=1)
    {
        if ((runs & 65535) == 0) printf("%10ld\n", runs);
        for (x = 0; x < WIDTH; ++x)
        for (y = 0; y < HEIGHT; ++y)
        if (NewRandomInteger < prob)
           // lat[x][y] = 1;  for model I
          lat[x][y] = 1 + (NewRandomInteger & 1);   //note, (NewRandomInteger & 1) gives 1 or 0 with equal probability
        else
            lat[x][y] = 0;
     
        nocctot = 0;

        for (xo = 0; xo < WIDTH; ++xo)
        for (yo = 0; yo < HEIGHT; ++yo)
        if (lat[xo][yo] > 0)
        {
            //printf("%10ld%10ld%10ld\n", xo, yo, lat[xo][yo]);
            gptr = pptr = 0;
            ++nctot;
            PutOnStack(xo,yo,lat[xo][yo])
            lat[xo][yo] = 0;
            
            do
            {
                GetFromStack(x,y,z)
                //z = orientation
                
                //printf("a %10ld%10ld%10ld%10ld\n", x, y, z, lat[x&W][y&H]);
                
                if(z == 1)
                {
                    for (dir = 0; dir < DIRMAX; ++dir)
                    {
                        xp = x + delx[dir]; //not switched
                        yp = y + dely[dir];
                        
                        if (lat[xp & W][yp & H])
                        {
                            PutOnStack(xp,yp,lat[xp & W][yp & H])
                            lat[xp & W][yp & H] = 0;
                        }
                    }
                    
                    for (dir = 0; dir < DIRMAX; ++dir)
                    {
                        xp = x + dely[dir];  ///switched
                        yp = y + delx[dir];
                        
                        if (lat[xp & W][yp & H] == 2)
                        {
                            PutOnStack(xp,yp,lat[xp & W][yp & H])
                            lat[xp & W][yp & H] = 0;
                        }
                    }
                }
                else //z=2
                {
                    for (dir = 0; dir < DIRMAX; ++dir)
                    {
                        xp = x + dely[dir]; //switched
                        yp = y + delx[dir];
                        
                        if (lat[xp & W][yp & H])
                        {
                            PutOnStack(xp,yp,lat[xp & W][yp & H])
                            lat[xp & W][yp & H] = 0;
                        }
                    }
                    for (dir = 0; dir < DIRMAX; ++dir)
                    {
                        xp = x + delx[dir];  ///not switched
                        yp = y + dely[dir];
                        
                        if (lat[xp & W][yp & H] == 1)
                        {
                            PutOnStack(xp,yp,lat[xp & W][yp & H])
                            lat[xp & W][yp & H] = 0;
                        }
                    }
                    //printf("here %10ld%10ld%10ld\n", xo, yo, lat[xo][yo]);
                
                }
            }
            while (gptr != pptr);
            
//printf("%10ld%10ld%10ld%10ld\n", nocc, xo, yo, zo);

            ++bin[(int)(log(pptr)/log(1.9999999))];
            sbin[(int)(log(pptr)/log(1.9999999))]+=pptr;
            nocctot += pptr;
        }
       // printf("%10ld\n", nocctot);
        
        if ((runs & PRINTFREQ) == 0)
        {
            fp1 = fopen(OUTFILE, "w");
            printf("%10d%14.6f%12ld%10d\n", HEIGHT, PROB, runs, SEED);
            fprintf(fp1,"%14.6f%12ld%10d\n", PROB, runs, SEED);
            printf("%10ld%15.8f%10ld\n", nctot, nctot*1.0/(runs*HEIGHT*WIDTH), nocctot);
            
            for (i = 1; i < 32; ++i)
            {
                if(bin[i]) printf("%10ld%15ld%18.8f%18.8f\n", i, bin[i],bin[i]*1.0/runs,sbin[i]*1.0/runs);
                fprintf(fp1,"%10ld%18.10e\n", i, bin[i]*1.0/runs);
            }
            fclose(fp1);
        } // if runs
   } // for runs
}

void randinit(long seed)
{
    double a, ee = -1 + 1/2147483648.0;
    long i;
    extern long nd, ra[M+1];
    
    a = seed/2147483648.0;
    for (nd = 0; nd <= M; nd++)
    {
        a *= 16807;
        a += ee * (long)(a);
        if (a >= 1) a += ee;
        ra[nd] = (long) (2147483648.0 * a);
    }
    nd = M;
    for(i = 0; i<100001; i++)
        NewRandomInteger;
}
