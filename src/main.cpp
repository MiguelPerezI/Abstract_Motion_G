#include <iostream>
#include <cstdlib>
#include "SimplexAbstract.hpp"
#include <math.h>
#include <list>
#include <time.h>
#include "SimplexAlpha.hpp"
#include "Covering_CORE.hpp"
#include "OptimizedCovering.hpp"
//#include "Lex.hpp"
#include <mpi.h>
#include <stdio.h>
#include <cstdio>
#include "LocalSearch.hpp"
#include <functional>

using namespace std;

LocalSearch local;
Complex komplex;
SimplicialMap map0, map1;
ComplexProduct KxK;
SubComplexJ L;

Covering covering;
OptimizedCovering O;

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main(int argc, char **argv) {

	struct timespec start, end;
    	clock_gettime(CLOCK_MONOTONIC, &start);
    	ios_base::sync_with_stdio(false);

        srand(time(NULL));
	
   	char *ptr;
   	double ret;

   	ret = strtod(argv[3], &ptr);


        int M = atoi(argv[2]);
        double r = ret;
        int Ns = atoi(argv[4]);

	int num = atoi(argv[1]);	

        int maxInt = 3;
        int numOfMaxSimplex = 5;
        komplex.initComplex(numOfMaxSimplex, maxInt + 1);

        komplex.K.A[0][0].initSimplex(2, 0, 1);
        komplex.K.A[0][1].initSimplex(2, 1, 2);
        komplex.K.A[0][2].initSimplex(2, 0, 2);
	komplex.K.A[0][3].initSimplex(2, 1, 3);
	komplex.K.A[0][4].initSimplex(2, 0, 3);

        KxK.initComplexProduct(komplex);
        KxK.escMaximalSimplices();
        map1.projection1(KxK);
        map0.projection2(KxK);

        L.initSubComplexJ(50);
        int counting = 0;

        for (int i = 0; i < KxK.listOfFacets.m; i++) {
                for (int j = 0; j < KxK.listOfFacets.rowLength.getA(0, i); j++) {
                        L.initA(counting, KxK.listOfFacets.A[i][j]);
                        counting += 1;
                }
        }
	
	
	 //L.initA(0, KxK.listOfFacets.A[6][1]);
	 //L.initA(1, KxK.listOfFacets.A[4][1]);

	 //L.initA(2, KxK.listOfFacets.A[5][1]);
	 //L.initA(3, KxK.listOfFacets.A[6][0]);

	 //L.initA(4, KxK.listOfFacets.A[0][1]);
	 //L.initA(5, KxK.listOfFacets.A[1][0]);
	 //
	 //L.initA(6, KxK.listOfFacets.A[1][1]);
	 //L.initA(7, KxK.listOfFacets.A[8][0]);


	 L.initZero_Skeleton();

	 L.escSubComplexJ();
	
         komplex.initAdjMat();
         komplex.graph.addWeight(0, 1, 1);
         komplex.graph.addWeight(0, 2, 1);
         komplex.graph.addWeight(1, 2, 1);
	 komplex.graph.addWeight(0, 3, 1);
	 komplex.graph.addWeight(1, 3, 1);
	
	//covering.initCovering(L, komplex, map1, map0, M, r, argc, argv);	
	////covering.runCovering(L, komplex, map1, map0, M, r, argc, argv);
	////covering.runRCC(L, komplex, map0, map1, M, r, argc, argv);
	//covering.plusFacet(L, KxK.listOfFacets.A[5][0], komplex, map1, map0, M, r, argc, argv);
	////covering.J.escSubComplexJ();
	//	
	////covering.plusSimplex(0, KxK.listOfFacets.A[5][1], L, komplex, map0, map1, M, r, argc, argv);

	////L.escSubComplexJ();
	//local.initLocalSearch(covering.J, komplex, map1, map0, M, r);

	//for (int i = 0; i < 10; i++)	
	//	if (local.ghost != 420)
	//		local.updateLocalSearch(covering.J, komplex, map1, map0, M, r);	
	//
	//local.escReducedChainResults(covering.J, map1);

	//covering.endCovering();
	//	

	O.initOptimizedCovering_CORE(L, komplex, map0, map1, M, r,  argc, argv);
	O.runOptimizedCovering_CORE(num, L, komplex, map0, map1, M, Ns,r, argc, argv);
	
	 //cout << "\n\n-->Argumento " << num << "\n";
	O.endOptimizedCovering();
        cout << "\n\n\n";


 	clock_gettime(CLOCK_MONOTONIC, &end);
    	double time_taken;
    	time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    	time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;

    	cout << "El tiempo que el programa tomó es: " << fixed
         << time_taken << setprecision(9);
    	cout << " sec" << endl;

    return 0;
}
