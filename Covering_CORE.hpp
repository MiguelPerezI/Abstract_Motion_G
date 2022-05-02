#include <stdarg.h>

#ifndef COVERINGCORE
#define COVERINGCORE

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <list>
#include <iterator>
#include "SimplexAbstract.hpp"
#include "SimplexAlpha.hpp"
#include "LocalSearch.hpp"
#include "Lex.hpp"
#include <stdio.h>
#include <mpi.h>


using namespace std;

int getAcomp(list<int> _list, int _i){
    list<int>::iterator it = _list.begin();
    for(int i=0; i<_i; i++){
        ++it;
    }
    return *it;
}

int get(list<int> _list, int _i){
    list<int>::iterator it = _list.begin();
    for(int i=0; i<_i; i++){
        ++it;
    }
    return *it;
}

class Covering {

	private:

	LocalSearch suchen, suchenPrime;
	SimplexAlpha sigma, sigmaPrime;
     

	int message[1000], message2[1000];
	int bigOsize, bigOsize2; 
	list <int> O, Oprime, cpu;
	int success;

	int thread, core, coreS;
	MPI_Status status;

	int sendList[10000];
	int sendRandom[10000];
	int sendState;
	int receiveState;
	int receiveList[10000];
	int sendSimplex[10000];
	int receiveSimplex[10000];

	
    	public:

        //RCC rcc;
	SubComplexJ J, Jprime;
        MatrixSimplexAlpha AA;
        MatrixSimplexAlpha P;
        SubComplexJ I;
        list <int> Acomp;

        TensorSimplexAlpha p, pOrder;
        TensorSimplicialMap pMap, pOrderMap;


        void initAcomp() {
            for (int i = 0; i < AA.getN(); i++) 
                Acomp.push_back(i);
        }

        TensorSimplicialMap getPOrderMap() {return pOrderMap;}
        TensorSimplicialMap getPMap() {return pMap;}



        void initCovering(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {
            cout << "\n\n\n<<<<<<<< Begining Covering >>>>>>>>>> " << endl;
            
            //L.escSubComplexJ();
            P.initMatrixSimplexAlpha(1, 1);
            P.initAByCopyOf(0, 0, L.listOfFacets.A[0][0]);
            P.n = 0;

            AA.initMatrixSimplexAlpha(1, 1);
            AA.initAByCopyOf(0, 0, L.listOfFacets.A[0][0]);

            //std :: cout << "-->We begin by copying the facets of L into A:";
            for(int i = 1; i < L.numSimplex; i++){
                AA.push(L.listOfFacets.A[0][i]);
            }

            
            //for (int i = 0; i < AA.getN(); i++) 
            //    Acomp.push_back(i);

            //cout << "\n-->Initializing auxilary list";
            //a.escListInt(Acomp);
            //cout << endl;


            cout << "\n--> L List Given By :\n";
            AA.escMatrixSimplexAlpha();

            //std :: cout << "\n-->Now we define an empty P = {}";
            p.initTensorSimplexAlpha(1, 1);
            pOrder.initTensorSimplexAlpha(1, 1);
            pMap.initTensorSimplicialMap();
            pOrderMap.initTensorSimplicialMap();
            //p.push(L.listOfFacets);
            //p.push(L.listOfFacets);
            //p.A[0][1].A[0][1].updateSimplexAlpha(L.listOfFacets.A[0][0]);
            //p.escTensorSimplexAlpha();

            //std :: cout << "\n-->Initializing a subComplex I that inherits all of A's facets";
            I.initSubComplexJ(1);
            I.listOfFacets.initAByCopyOf(0, 0, AA.A[0][0]);

            for(int i = 1; i < AA.getN(); i++)
                I.listOfFacets.push(AA.A[0][i]);

            I.initZero_Skeleton();
            
            //I.escSubComplexJ();
	    //
	    
		thread = MPI_Init(&argc, &argv);
                thread = MPI_Comm_rank(MPI_COMM_WORLD, &core);
                thread = MPI_Comm_size(MPI_COMM_WORLD, &coreS);
                srand(time(NULL));
	
        	string file = to_string(core);
        	string file_name = "core" + file + ".txt";


                //printf("\n\n<<<<< RCC >>>>>>\n\n");
                J.initSubComplexJ(1);
                Jprime.initSubComplexJ(1);
                //printf("-->Fetching random SimplexAlpha from L and saving to J.\n");
                J.initA(0, L.getRandomJ());
                J.initZero_Skeleton();
                
		Jprime.initA(0, L.getRandomJ());
                Jprime.initZero_Skeleton();
		sigma.makeCopySimplexAlpha(L.listOfFacets.A[0][0]);
		sigmaPrime.makeCopySimplexAlpha(L.listOfFacets.A[0][0]);

		srand(time(NULL) + core);
		for (int i = 0; i < coreS; i++) {
			sendList[i] = 0;
			receiveList[i] = 0;
		}
                suchen.initLocalSearch(L, K, a, b, M, r);    
	
        }


        void setPAuxilary(int ccc) {
            
            //cout << "\n-->Check for facets in J that intersect with facets in L";
            //compareSimplexAlpha(SimplexAlpha)
            if (ccc == 1) {
                //p.push(rcc.getJ_Facets());
                for (int iterJ = 0; iterJ < J.listOfFacets.getN(); iterJ++)
                    P.pushZero(J.listOfFacets.A[0][iterJ]);
                p.initA(0, 0, P);
                pOrder.initA(0, 0, P);
                //p.escTensorSimplexAlpha();
                //cout << "\n-->Tensor of size 1, pushing J directly";
            }
            else {
                //cout << "\nccc = " << ccc << endl;
                //cout << "\n-->Tensor of size "<< p.n <<"       --> Building P";
                //cout << "\n-->Number of Partitions " << p.n << endl;
                //p.escTensorSimplexAlpha();
                //cout << "\nnewJFacets := \n";
                //rcc.getJ_Facets().escMatrixSimplexAlpha();



                for (int iterJ = 0; iterJ < J.numSimplex; iterJ += 1) {
                    int bit = 0;
                    for (int i = 0; i < p.getN(); i++) {
                        //printf("\n-->Partition := %d of size %d", i, p.sizeOfPartition(i)); //could be zero
                        

                        int bb = 0;
                        while (bb < p.sizeOfPartition(i)) {

                                    if (J.listOfFacets.A[0][iterJ].compareSimplexAlpha(p.getPartitionSimplex(i, bb)) == 1) {
                                        bit = 1;
                                        i = p.getN() + 1;
                                        break;
                                    }

                            bb += 1;
                        }
                    }

                    if (bit == 0) {
                        //cout << "\n-->No Match, Adding Facet";
                        P.pushZero(J.listOfFacets.A[0][iterJ]);
                    }
                }


                //cout << "\n-->newP := ";
                //P.escMatrixSimplexAlpha();
    
                //cout << "\n-->Pushing P into Tensor";
                p.push(P);
                p.n += 1;

                pOrder.push(P);
                pOrder.n += 1;

                //cout << "size >>>>>>>>>:= " << p.n;
                //cout << "\n-->Tensor newSize " << p.getN();
            }


            
            //cout << "\n\n\n";
        }

        void resetSubComplexI_FromA() {

            I.resetSubComplexJ(AA.A[0][getAcomp(Acomp, 0)]);
            if (AA.getN() > 0)
                for (int i = 1; i < Acomp.size(); i++)
                    I.pushSimplexAlpha(AA.A[0][getAcomp(Acomp, i)]);
        }

        void eliminateP_From_A() {

            for (int i = 0; i < P.n; i++) {
                //cout << "-----------------------------------> " << i << endl;
                int count = 0;
                while (count < Acomp.size()) {

                    int RR = getAcomp(Acomp, count);
                    //cout << "\n-->Comparing \n";
                    //printf("-->A(%d) := ", RR); AA.getA(0, RR).escSimplexAlpha(); cout << "\n";
                    //printf("-->P(%d) := ", i);P.getA(0, i).escSimplexAlpha(); cout << "\n";
                    if (AA.getA(0, RR).compareSimplexAlpha(P.getA(0, i)) == 1) {
                        //cout << ">>>>>Difference FOUND<<<<< Removing From A:\n";
                        //printf("-->A(%d) := ", RR); AA.getA(0, RR).escSimplexAlpha(); cout << "\n";
                        Acomp.remove(RR);
                        break;
                    }

                    count += 1;
                }
            }

            //cout << "\n>>>>>>>>Elimination Finish<<<<<<<<\n";
        }

	void communication(int send) {
                        
			sendState = send;
                        for (int an_id = 0; an_id < coreS; an_id++)///////////TELL THE REST
                        	if (an_id != core)
                                	thread = MPI_Send(&sendState, 1, MPI_INT,
                                                        an_id, 2001, MPI_COMM_WORLD     );

                        for (int iter = 0; iter < coreS; iter++)/////////////////////RECEIVING NEWS FROM THE REST AND UPDATING
                                if (core != iter) {
                                        thread = MPI_Recv(&receiveState, 1, MPI_INT, iter, 2001, MPI_COMM_WORLD, &status);
                                        sendList[iter] = receiveState;
                        }

	}

	int chooseIndex() {
		
		int i = 0;
		int ret = -11;
		while (i < coreS) {
			if (sendList[i] > - 1) {
				ret = sendList[i];
				break;
			}
			i++;
		}

		return ret;
	}

	int resetIndex() {
		
		for (int i = 0; i < coreS; i++)
			sendList[i] = -11;
	
		return 0;
	}	

	

	void runRCC(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {
		
		int count  = 0;
		while (count < L.numSimplex) {
			sendList[core] = simplexScrutiny(count, L, K, a, b, M, r, argc, argv);
			communication(sendList[core]);
			int index = chooseIndex();
			J.popSimplexAlpha();
			if (index > -1) {
				sigma.updateSimplexAlpha(L.listOfFacets.A[0][index]);
                        	J.pushSimplexAlpha(sigma);
			}

			resetIndex();

			if (J.numSimplex == L.numSimplex)
				break;

			count += 1;
		}

		string cc = to_string(core);
		string ss = "core" + cc + ".txt";
		SimplexAlphaWriter sim = SimplexAlphaWriter(ss);
		sim.write(J.stringSubComplexJ());
		sim.closeSimplexAlphaWriter();

	}

	int plusFacet(SubComplexJ L, SimplexAlpha alpha, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {

		int ret = 0;
		J.resetSubComplexJ(L.listOfFacets.A[0][0]);
                for (int i = 1; i < L.numSimplex; i++)
                        J.pushSimplexAlpha(L.listOfFacets.A[0][i]);
		
		Jprime.resetSubComplexJ(alpha);

                int count  = 0;
                while (count < L.numSimplex) {
                        sendList[core] = simplexScrutiny(count, Jprime, K, a, b, M, r, argc, argv);
                        communication(sendList[core]);
                        int index = chooseIndex();
                        J.popSimplexAlpha();
                        if (index > -1) {
                                sigma.updateSimplexAlpha(L.listOfFacets.A[0][index]);
                                J.pushSimplexAlpha(sigma);
                        }

                        resetIndex();

                        if (J.numSimplex == L.numSimplex + 1) {
                                ret = 1;
				break;
			}

                        count += 1;
                }

                string cc = to_string(core);
                string ss = "plusSimplex" + cc + ".txt";
                SimplexAlphaWriter sim = SimplexAlphaWriter(ss);
                sim.write(alpha.stringSimplexAlpha());
		sim.write("\n\n");
		sim.write(J.stringSubComplexJ());
		if (ret == 1) sim.write("\n\nADDEDSIMPLEX");
                sim.closeSimplexAlphaWriter();
		return ret;
        }

	void runAddFacets(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {

                int count  = 0;
                while (count < L.numSimplex) {
                        sendList[core] = simplexScrutinyGeneral(L, K, a, b, M, r, argc, argv);
                        communication(sendList[core]);
                        int index = chooseIndex();
                        Jprime.popSimplexAlpha();
                        if (index > -1) {
                                sigma.updateSimplexAlpha(L.listOfFacets.A[0][index]);
                                Jprime.pushSimplexAlpha(sigma);
                        }

                        resetIndex();

                        if (Jprime.numSimplex == L.numSimplex)
                                break;

                        count += 1;
                }

                string cc = to_string(core);
                string ss = "core" + cc + ".txt";
                SimplexAlphaWriter sim = SimplexAlphaWriter(ss);
                sim.write(Jprime.stringSubComplexJ());
                sim.closeSimplexAlphaWriter();

        }

	//void runAddFacets(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {

        //        int count  = 0;
        //        while (count < L.numSimplex) {

        //                simplexScrutinyGeneral(L, K, a, b, M, r, argc, argv);
        //                if (Jprime.numSimplex == L.numSimplex)
        //                        break;

	//		count += 1;
        //        }
        //}
	
	void plusSimplex(int stb, const SimplexAlpha& s, SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {

                int cycleRCC = 0;
                while (  cycleRCC < 4) {

			 int RR0;
                        if (O.size() < coreS)
                                RR0 = get(O, core);
                        else {
                                RR0 = get(O, core%O.size());
                        }

                        //////////////////////////////////////////////////////////////Unión of J and {random(O)}
                        sigma.updateSimplexAlpha(s);
                        L.pushSimplexAlpha(sigma);
			L.escSubComplexJ();
                        //////////////////////////////////////////////////////////////FINDING A MAP OVER JU{random(O)}
                        suchen.ghostLocalSearch(L, K, a, b, M, r);
                        int callack = 0;
                        if (suchen.ghost == 420) {////////////////////////////////////MAP FOUND!! DO:
                                cout << "\nSIMPLEX FOUND BY CPU("<< core << ")\n";
                                //cout << "\n\n";
                                sendState = RR0;

                                for (int an_id = 0; an_id < coreS; an_id++)///////////TELL THE REST
                                        if (an_id != core)
                                        thread = MPI_Send(&sendState, 1, MPI_INT,
                                                        an_id, 2001, MPI_COMM_WORLD     );
                        } else {//////////////////////////////////////////////////////NO MAP FOUND
                                //cout << "\nSSSSSHHHHHHHHHHHHHH\n\n";

                                L.popSimplexAlpha();

                                sendState = -RR0 - 1;
                                for (int an_id = 0; an_id < coreS; an_id++)//////////TELL THE OVERS
                                        if (an_id != core)
                                        thread = MPI_Send(&sendState, 1, MPI_INT,
                                                        an_id, 2001, MPI_COMM_WORLD     );

                        }

                        sendList[core] = sendState;

                        for (int iter = 0; iter < coreS; iter++)/////////////////////RECEIVING NEWS FROM THE REST AND UPDATING
                                if (core != iter) {
                                        thread = MPI_Recv(&receiveState, 1, MPI_INT, iter, 2001, MPI_COMM_WORLD, &status);
                                        sendList[iter] = receiveState;
                        }

                        for (int iter = 0; iter < coreS; iter++){

                                if (sendList[iter] > -1)//////////////WE MUST CHOOSE ONLY ONE SIMPLEX
                                {
                                //      cout << "\n==============================================>> SIGMA := " << sendList[iter] << "\n\n";
                                        if (sendList[core] > -1) J.popSimplexAlpha();//SO WE REMOVE WHAT WE UNITED

                                        sigma.updateSimplexAlpha(s);//AND WE KEEP ONLY ONE
                                        L.pushSimplexAlpha(sigma);
                                        cycleRCC += 20;
                                        iter += 10000000000;

                                }
                        }

			cycleRCC++;
                }

        }
	
			
	int simplexScrutiny(int stb, SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {
							//////RUN RCC/////
		int b1 = 0;
                int nn1 = 0;
			for (int i = 0; i < L.numSimplex; i++) {
                	                for (int j = 0; j < J.numSimplex; j++)
                	                        b1 += L.listOfFacets.A[0][i].compareSimplexAlpha(J.listOfFacets.A[0][j]);
                	                if (b1 < 1) {
                	                        if (nn1 == 0)
                	                                O.push_back(i);
                	                        else
                	                                O.push_back(i);
                	                        nn1 += 1;
                	                        }
                	                b1 = 0;
                	}

			cout << "\n";
			a.escListInt(O);
			cout << "\n";
		int cycleRCC = 0;
		int count0 = -11;
		//while (  cycleRCC < 1) {
			//////////////////////////////////////////////////////////////Choosing random O element
                        
			int i0 = core%(O.size());
			int RR0;
		            RR0 = get(O, i0);

				//////////////////////////////////////////////////////////////Unión of J and {random(O)}
			sigma.updateSimplexAlpha(L.listOfFacets.A[0][RR0]);
			J.pushSimplexAlpha(sigma);
				//////////////////////////////////////////////////////////////FINDING A MAP OVER JU{random(O)}
			suchen.ghostLocalSearch(J, K, a, b, M, r);
				if (suchen.ghost == 420) {////////////////////////////////////MAP FOUND!! DO:
					//cout << "core(" << core << ") found simplex at index := " << RR0 << endl;
                			O.clear();
					return RR0;
				} else {//////////////////////////////////////////////////////NO MAP FOUND
					//J.popSimplexAlpha();

					//cout << "core(" << core << ") found simplex at index := " << -11 << endl;
					O.clear();
					return -11;
				}
	}




	int simplexScrutinyGeneral(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {
                                                        //////RUN RCC/////
                int b1 = 0;
                int nn1 = 0;
                        for (int i = 0; i < L.numSimplex; i++) {
                                        for (int j = 0; j < Jprime.numSimplex; j++)
                                                b1 += L.listOfFacets.A[0][i].compareSimplexAlpha(Jprime.listOfFacets.A[0][j]);
                                        if (b1 < 1) {
                                                if (nn1 == 0)
                                                        Oprime.push_back(i);
                                                else
                                                        Oprime.push_back(i);
                                                nn1 += 1;
                                                }
                                        b1 = 0;
                        }

                        cout << "\n";
                        a.escListInt(Oprime);
                        cout << "\n";
                int cycleRCC = 0;
                int count0 = -11;
                //while (  cycleRCC < 1) {
                        //////////////////////////////////////////////////////////////Choosing random O element

                        int i0 = core%(Oprime.size());
                        int RR0;
                            RR0 = get(Oprime, i0);

                                //////////////////////////////////////////////////////////////Unión of J and {random(O)}
                        sigma.updateSimplexAlpha(L.listOfFacets.A[0][RR0]);
                        Jprime.pushSimplexAlpha(sigma);
                                //////////////////////////////////////////////////////////////FINDING A MAP OVER JU{random(O)}
                        suchen.ghostLocalSearch(Jprime, K, a, b, M, r);
                                if (suchen.ghost == 420) {////////////////////////////////////MAP FOUND!! DO:
                                        //cout << "core(" << core << ") found simplex at index := " << RR0 << endl;
                                        Oprime.clear();
                                        return RR0;
                                } else {//////////////////////////////////////////////////////NO MAP FOUND
                                        //J.popSimplexAlpha();

                                        //cout << "core(" << core << ") found simplex at index := " << -11 << endl;
                                        Oprime.clear();
                                        return -11;
                                }
        }








	//void simplexScrutinyGeneral(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {
	//	///////////////////////////////////////////////////////////////
	//			///////////////////////////////////////////////////////////////
	//						//////////////////
	//						//////////////////
	//						//////RUN RCC/////
	//						//////////////////
	//						//////////////////
	//			


	//	int b1 = 0;
        //        int nn1 = 0;
        //        //printf("\n<<<<<<<<ADDFACET>>>>>>>>>"); First we Serach for all the simplicies not in J
        //        for (int i = 0; i < L.numSimplex; i++) {
        //                        for (int j = 0; j < Jprime.numSimplex; j++)
        //                                b1 += L.listOfFacets.A[0][i].compareSimplexAlpha(Jprime.listOfFacets.A[0][j]);
        //                        if (b1 == 0) {
        //                                if (nn1 == 0)
        //                                        Oprime.push_back(i);
        //                                else
        //                                        Oprime.push_back(i);
        //                                nn1 += 1;
        //                                }
        //                        b1 = 0;
        //        }
	//	int cycleRCC = 0;

	//	while (  Oprime.size() > 0) {


	//		//cout << "\nSearching for Simplicies in CORE 0\n\n";

	//		//////////////////////////////////////////////////////////////Choosing random O element
	//		int i0 = rand()%(Oprime.size());
	//		int RR0 = get(Oprime, core%Oprime.size());

	//		//////////////////////////////////////////////////////////////Unión of J and {random(O)}
	//		sigmaPrime.updateSimplexAlpha(L.listOfFacets.A[0][RR0]);
	//		Jprime.pushSimplexAlpha(sigmaPrime);


	//		//////////////////////////////////////////////////////////////FINDING A MAP OVER JU{random(O)}
	//		suchen.ghostLocalSearch(Jprime, K, a, b, M, r);

	//		int callack = 0;

	//		if (suchen.ghost == 420) {////////////////////////////////////MAP FOUND!! DO:
	//			//cout << "\n\n";
	//			sendState = RR0;

	//			for (int an_id = 0; an_id < coreS; an_id++)///////////TELL THE REST
	//				if (an_id != core)
	//				thread = MPI_Send(&sendState, 1, MPI_INT,
	//					        an_id, 2001, MPI_COMM_WORLD	);
	//		} else {//////////////////////////////////////////////////////NO MAP FOUND
	//			//cout << "\nNO SIMPLEX FOUND IN CORE 0, POPING SIMPLEX FROM J\n\n";
	//			Jprime.popSimplexAlpha();

	//			sendState = -RR0 - 1;
	//			for (int an_id = 0; an_id < coreS; an_id++)//////////TELL THE OVERS
	//				if (an_id != core)
        //                                thread = MPI_Send(&sendState, 1, MPI_INT,
        //                                                an_id, 2001, MPI_COMM_WORLD     );
	//		}

	//		sendList[core] = sendState;

	//		for (int iter = 0; iter < coreS; iter++)/////////////////////RECEIVING NEWS FROM THE REST AND UPDATING
	//			if (core != iter) {
	//				thread = MPI_Recv(&receiveState, 1, MPI_INT, iter, 2001, MPI_COMM_WORLD, &status);
	//				sendList[iter] = receiveState;
	//		}


	//		int check = 0;
	//		for (int iter = 0; iter < coreS; iter++){

	//			if (sendList[iter] <= -1) {check += 1;}//////////////WE MUST CHOOSE ONLY ONE SIMPLEX
	//			else {
	//				//cout << "\n==============================================>> SIGMA := " << sendList[iter] << "\n\n";
        //                                if (sendList[core] > -1) Jprime.popSimplexAlpha();//SO WE REMOVE WHAT WE UNITED

        //                                sigmaPrime.updateSimplexAlpha(L.listOfFacets.A[0][sendList[iter]]);//AND WE KEEP ONLY ONE
        //                                Jprime.pushSimplexAlpha(sigmaPrime);
        //                                Oprime.remove(sendList[iter]);//////////////////////MAKING SURE THAT THE O-LIST's REMAIN PARALLEL
        //                                iter += 10000000000;
	//			}


	//		}

	//			for (int iter = 0; iter < coreS; iter++) {////////////REMOVING all FAILED ATEMPS
	//				if (sendList[iter] <= -1)
	//					Oprime.remove(-1 * (sendList[iter] + 1));
	//			}


	//		//cout << "\n\n<<<<<<<<<<<<<<<<<<\n";
        //        	//a.escListInt(O);/////////////////////////////////////////////WE PRINT OUR O LIST
        //        	//cout << "\n<<<<<<<<<<<<<<<<<<\n\n";


	//		cycleRCC++;
	//	}



	//	Oprime.clear();
	//	//Jprime.escSubComplexJ();
	//			////////////////////////////////////////////////////////////////
	//		 	////////////END OF RCC//////////////////////////////////////////
	//		    	////////////////////////////////////////////////////////////////
	//		    	////////////////////////////////////////////////////////////////
	//		    	////////////////////////////////////////////////////////////////
	//
	//}







	void resetRCCSubComplexJ(SubComplexJ L) {

                int RR = rand()%(L.listOfFacets.n);
                sendState = RR;
                sendRandom[core] = RR;
                for (int an_id = 0; an_id < coreS; an_id++)//////////TELL THE OVERS
                        if (an_id != core)
                                thread = MPI_Send(&sendState, 1, MPI_INT, an_id, 2001, MPI_COMM_WORLD);

                for (int iter = 0; iter < coreS; iter++)/////////////////////RECEIVING NEWS FROM THE REST AND UPDATING
                                if (core != iter) {
                                        thread = MPI_Recv(&receiveState, 1, MPI_INT, iter, 2001, MPI_COMM_WORLD, &status);
                                        sendRandom[iter] = receiveState;
                                }
                J.resetSubComplexJ(L.listOfFacets.A[0][sendRandom[0]]);
        }


        void resetAddFacetsSubComplexJ(SubComplexJ JsubL) {

                Jprime.resetSubComplexJ(JsubL.listOfFacets.A[0][0]);
                for (int i = 1; i < JsubL.numSimplex; i++)
                        Jprime.pushSimplexAlpha(JsubL.listOfFacets.A[0][i]);
        }

	void resetI(MatrixSimplexAlpha M) {

                I.resetSubComplexJ(M.A[0][0]);
                for (int i = 1; i < M.getM(); i++)
                        I.pushSimplexAlpha(M.A[0][i]);
        }

	void AUX_CORE(SubComplexJ L, SubComplexJ * JD) {

		int RR = rand()%(L.listOfFacets.n);
		sendState = RR;
                sendRandom[core] = RR;
                for (int an_id = 0; an_id < coreS; an_id++)//////////TELL THE OVERS
                        if (an_id != core)
                                thread = MPI_Send(&sendState, 1, MPI_INT, an_id, 2001, MPI_COMM_WORLD);

                for (int iter = 0; iter < coreS; iter++)/////////////////////RECEIVING NEWS FROM THE REST AND UPDATING
                                if (core != iter) {
                                        thread = MPI_Recv(&receiveState, 1, MPI_INT, iter, 2001, MPI_COMM_WORLD, &status);
                                        sendRandom[iter] = receiveState;
                                }


                JD->resetSubComplexJ(L.listOfFacets.A[0][sendRandom[0]]);
                JD->listOfFacets.n = 0;
                JD->escSubComplexJ();
	}




        void runCovering(SubComplexJ L, Complex K, SimplicialMap a, SimplicialMap b, int M, double r, int argc, char **argv) {

            cout << "\n\n--> Covering Running";
            cout << "\n-->Running Main Loop";
            int cycle = 1;

            initAcomp();

            while (Acomp.size() != 0) {

                //cout << "\n\n\n -->Cycle " << cycle;

                //cout << "\n debug (0)---------------------------------------\n";
                resetSubComplexI_FromA();

                //cout << "\n\n--> reset of Subcomplex I := \n";
                //I.escSubComplexJ(); cout << "\n\n";
                //cout << "\n debug (1)---------------------------------------\n";
                if (Acomp.size() > 1) { 

                 //           cout << "\n debug (2)---------------------------------------\n";
                            
			    this->runRCC( I, K, a, b, M, r, argc, argv);
                            setPAuxilary(cycle);
                            eliminateP_From_A();
                            P.resetPCovering();
                            if (Acomp.size() > 1)
                            	resetSubcomplexJ2(Acomp, AA, cycle);
                            
                } else {

			if (core == 0)cout << "\n----------------------------------------------------------------------->ONE LEFT\n\n";
			int bb = 0;	
			for (int i = 0; i < pOrder.n; i++) {
				
				resetI(pOrder[0][i]);				
				bb = plusFacet(I, AA.getA(0, getAcomp(Acomp, 0)), K, a, b, M, r, argc, argv);
				if (bb == 1) {
					p.pushSimplexAlpha(0, i, AA.getA(0, getAcomp(Acomp, 0)));
					pOrder.pushSimplexAlpha(0, i, AA.getA(0, getAcomp(Acomp, 0)));
					 P.pushZero(AA.getA(0, getAcomp(Acomp, 0)));
					 eliminateP_From_A();
					 P.resetPCovering();
					i += (2*pOrder.n);
				}

			}

			if (bb == 0) {

				//for (int iterJ = 0; iterJ < rcc.getJSize(); iterJ++)
                        //    cout << "\n debug (13)---------------------------------------\n";
                            P.pushZero(AA.getA(0, getAcomp(Acomp, 0)));
                        //    cout << "\n debug (14)---------------------------------------\n";
                            p.push(P);
                        //    cout << "\n debug (15)---------------------------------------\n";
                            p.n += 1;

                        //    cout << "\n debug (16)---------------------------------------\n";
                            pOrder.push(P);
                        //    cout << "\n debug (17)---------------------------------------\n";
                            pOrder.n += 1;

                            //break;
                        //    cout << "\n debug (18)---------------------------------------\n";
                            eliminateP_From_A();
                        //    cout << "\n debug (19)---------------------------------------\n";
                            //cout << "\n-->new A list := ";
                            //a.escListInt(Acomp);
                            P.resetPCovering();
                        //    cout << "\n debug (20)---------------------------------------\n";
			}
                }//

                //Acomp.pop_back();
                //cout << "\n\n";
                //I.escSubComplexJ();
                //I.escZero_Skeleton();

                cycle += 1;
                //cout << "\n---------------------------------------------------------------------------------------------------------->\n";

            }


            //I.escSubComplexJ();

            //cout << "\n\n<<<<<<<<<<<End of Covering>>>>>>>>>>>>>>>\n\n";
            
            //p.escTensorSimplexAlpha();
	    string cc = to_string(core);
            string ss = "core" + cc + ".txt";
            SimplexAlphaWriter sim = SimplexAlphaWriter(ss);
            sim.write(p.stringTensorSimplexAlpha());
            sim.closeSimplexAlphaWriter();

            int si = 0;
            for (int i = 0; i < p.n; i++) {
                si += p.sizeOfPartition(i);
            }

            if (si != L.numSimplex) {
		    thread = MPI_Finalize();
                throw runtime_error("WROOOOOOOOOOOOONG");
		cout << "\n==================>> " << si << "\n\n";
	    }
		
   	    cout << "\n\n<<<<<<<<<<<End of Covering>>>>>>>>>>>>>>>\n\n";	    
	    //thread = MPI_Finalize();
	    return;
        }

	void endCovering() {
		//cout << "\n\n===================>>Covering Finalized\n\n";
		//rcc.endRCC_CORE();
		thread = MPI_Finalize();
	}


    void orderPartition(SimplicialMap a) {

        int count = 0;
        for (int i = 0; i < p.n; i++) {
                if (count == 0)
                    pOrder.initA(0, 0, p.A[0][i]);
                else {
                    pOrder.push(p.A[0][i]);
                    pOrder.n += 1;
                }
                count += 1;
            }

        for (int i = 0; i < p.n; i++) {
            pOrder.replaceAt(i, p.A[0][i]);
        }

        //cout << "\n\n-->ordered tensor should be:\n";
        pOrder.escTensorSimplexAlpha();
    }

    int getPM() {
        return pOrder.getM();
    }

    int getPN() {
        return pOrder.getN();
    }

	void resetSubcomplexJ2(list <int> aaa0, MatrixSimplexAlpha AA0, int ccc) {
		
		int ran = rand()%(aaa0.size());
		int RR = get(aaa0, ran);
		sendState = RR;
		sendRandom[core] = RR;
		for (int an_id = 0; an_id < coreS; an_id++)//////////TELL THE OVERS
                	if (an_id != core)
                        	thread = MPI_Send(&sendState, 1, MPI_INT, an_id, 2001, MPI_COMM_WORLD);
		
		for (int iter = 0; iter < coreS; iter++)/////////////////////RECEIVING NEWS FROM THE REST AND UPDATING
                                if (core != iter) {
                                        thread = MPI_Recv(&receiveState, 1, MPI_INT, iter, 2001, MPI_COMM_WORLD, &status);
					sendRandom[iter] = receiveState;
				}
		
		J.resetSubComplexJ(AA0.A[0][sendRandom[0]]);
//		cout << "\n\n\n				XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX CORE  := " << core  << "    cycle :=  " << ccc << "      " << sendRandom[0] << "\n\n";

	}





    int getPartitionSize(int i) {
        return  pOrder.sizeOfPartition(i);
    }



    // { [m0, n0], [m1, n1],  [m2, n2],  [m3, n3]}
    //[0, ni] = sigma(0, ni)

    SimplexAlpha getPartitionSimplex(int i, int j) {
        return pOrder.getPartitionSimplex(i, j);
    }

    void replaceAt_i(int i, MatrixSimplexAlpha M) {
        pOrder.replaceAt(i, M);
    }
	
    void escPOrder() {
        pOrder.escTensorSimplexAlpha();
    }







//	void escRCCResults() {
//
//		cout << "\n\n-->RCC found \n";
//		//escMapeos(J);
//		cout << "\n Defined over the subcomplex J :\n";
//		J.escSubComplexJ();
//	}
//
//	SubComplexJ getJ() {
//		return J;
//	}
//
//	MatrixSimplexAlpha getJListOfFacets() {
//		return J.listOfFacets;
//	}
//
//	MatrixSimplexAlpha getJ_Facets() {
//		return J.listOfFacets;
//	}
//
//	//for listOfFacets oneDimensional
//	int getJSize() {
//		return J.listOfFacets.getN();
//	}
//
//	//for listOfFacets oneDimensional
//	SimplexAlpha getJ_Facet(int k) {
//		return J.listOfFacets.A[0][k];
//	}
};

#endif
