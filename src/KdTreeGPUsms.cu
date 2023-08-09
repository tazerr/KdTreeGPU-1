//
//  KdTreeGPUsms.cu
//
//  Created by John Robinson on 7/15/15.
//  Copyright (c) 2015 John Robinson. All rights reserved.
/*
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSEARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
//

/*
 * Copyright (c) 2015, Russell A. Brown
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSEARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* @(#)kdTreeSingleThread.cc	1.61 04/13/15 */

/*
 * The k-d tree was described by Jon Bentley in "Multidimensional Binary Search Trees
 * Used for Associative Searching", CACM 18(9): 509-517, 1975.  For k dimensions and
 * n elements of data, a balanced k-d tree is built in O(kn log n) + O((k+1)n log n)
 * time by first sorting the data in each of k dimensions, then building the k-d tree
 * in a manner that preserves the order of the k sorts while recursively partitioning
 * the data at each level of the k-d tree. No further sorting is necessary.  Moreover,
 * it is possible to replace the O((k+1)n log n) term with a O((k-1)n log n) term but
 * this approach sacrifices the generality of building the k-d tree for points of any
 * number of dimensions.
 */

#include <stdbool.h>
#include <stdlib.h>
#include <vector>
#include <list>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>

using std::setprecision;

using namespace std;

#include "Gpu.h"
#include "KdNode.h"


//#if __cplusplus != 201103L
#if 0

#include <chrono>
#define TIMER_DECLARATION()						\
		auto startTime = std::chrono::high_resolution_clock::now();		\
		auto endTime = <std::chrono::high_resolution_clock::now();
#define TIMER_START()							\
		startTime = std::chrono::high_resolution_clock::now(); // high_resolution_clock::is_steady
#define TIMER_STOP(__TIMED)						\
		endTime = std::chrono::high_resolution_clock::now();			\
		__TIMED = (std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - startTime).count())/1000.0

#elif defined(MACH)

#define TIMER_DECLARATION()				\
		struct timespec startTime, endTime;
#define TIMER_START()						\
		mach_gettime(CLOCK_REALTIME, &startTime);
#define TIMER_STOP(__TIMED)					\
		clock_gettime(CLOCK_REALTIME, &endTime);			\
		__TIMED = (endTime.tv_sec - startTime.tv_sec) +			\
		1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec))

#else

#define TIMER_DECLARATION()				\
		struct timespec startTime, endTime;
#define TIMER_START()						\
		clock_gettime(CLOCK_REALTIME, &startTime);
#define TIMER_STOP(__TIMED)					\
		clock_gettime(CLOCK_REALTIME, &endTime);			\
		__TIMED = (endTime.tv_sec - startTime.tv_sec) +			\
		1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec))

#endif

Gpu *gpu;

/*
 * The superKeyCompare method compares two sint arrays in all k dimensions,
 * and uses the sorting or partition coordinate as the most significant dimension.
 *
 * calling parameters:
 *
 * a - a int*
 * b - a int*
 * p - the most significant dimension
 * dim - the number of dimensions
 *
 * returns: +1, 0 or -1 as the result of comparing two sint arrays
 */
KdCoord KdNode::superKeyCompare(const KdCoord *a, const KdCoord *b, const sint p, const sint dim)
{
	KdCoord diff = 0;
	for (sint i = 0; i < dim; i++) {
		sint r = i + p;
		r = (r < dim) ? r : r - dim;
		diff = a[r] - b[r];
		if (diff != 0) {
			break;
		}
	}
	return diff;
}

/*
 * Walk the k-d tree and check that the children of a node are in the correct branch of that node.
 *
 * calling parameters:
 *
 * dim - the number of dimensions
 * depth - the depth in the k-d tree
 *
 * returns: a count of the number of kdNodes in the k-d tree
 */
sint KdNode::verifyKdTree( const KdNode kdNodes[], const KdCoord coords[], const sint dim, const sint depth) const
{
	sint count = 1 ;

	// The partition cycles as x, y, z, w...
	sint axis = depth % dim;

	if (ltChild != -1) {
		if (superKeyCompare(coords+kdNodes[ltChild].tuple*dim, coords+tuple*dim, axis, dim) >= 0) {
			cout << "At Depth " << depth << " LT child is > node on axis " << axis << "!" << endl;
			printTuple(coords+tuple*dim, dim);
			cout << " < [" << ltChild << "]";
			printTuple(coords+kdNodes[ltChild].tuple*dim, dim);
			cout << endl;
			exit(1);
		}
		count += kdNodes[ltChild].verifyKdTree(kdNodes, coords, dim, depth + 1);
	}
	if (gtChild != -1) {
		if (superKeyCompare(coords+kdNodes[gtChild].tuple*dim, coords+tuple*dim, axis, dim) <= 0) {
			cout << "At Depth " << depth << " GT child is < node on axis " << axis << "!" << endl;
			printTuple(coords+tuple*dim, dim);
			cout << " > [" << gtChild << "]";
			printTuple(coords+kdNodes[gtChild].tuple*dim, dim);
			cout << endl;
			exit(1);
		}
		count += kdNodes[gtChild].verifyKdTree(kdNodes, coords, dim, depth + 1);
	}
	return count;
}

/*
 * The createKdTree function performs the necessary initialization then calls the buildKdTree function.
 *
 * calling parameters:
 *
 * coordinates - a vector<int*> of references to each of the (x, y, z, w...) tuples
 * numDimensions - the number of dimensions
 *
 * returns: a KdNode pointer to the root of the k-d tree
 */
KdNode* KdNode::createKdTree(KdNode kdNodes[], KdCoord coordinates[],  const sint numDimensions, const sint numTuples)
{

	TIMER_DECLARATION();

	TIMER_START();
	Gpu::initializeKdNodesArray(coordinates, numTuples, numDimensions);
	cudaDeviceSynchronize();
	TIMER_STOP (double initTime);

	// Sort the reference array using multiple threads if possible.

	TIMER_START();
	sint end[numDimensions]; // Array used to collect results of the remove duplicates function
	Gpu::mergeSort(end, numTuples, numDimensions);
	TIMER_STOP (double sortTime);

	// Check that the same number of references was removed from each reference array.
	for (sint i = 0; i < numDimensions-1; i++) {
		if (end[i] < 0) {
			cout << "removeDuplicates failed on dimension " << i << endl;
			cout << end[0];
			for (sint k = 1;  k<numDimensions; k++) cout << ", " << end[k] ;
			cout << endl;
			exit(1);
		}
		for (sint j = i + 1; j < numDimensions; j++) {
			if ( end[i] != end[j] ) {
				cout << "Duplicate removal error" << endl;
				cout << end[0];
				for (sint k = 1;  k<numDimensions; k++) cout << ", " << end[k] ;
				cout << endl;
				exit(1);
			}
		}
	}
	cout << numTuples-end[0] << " equal nodes removed. "<< endl;

	// Build the k-d tree.
	TIMER_START();
	//  refIdx_t root = gpu->startBuildKdTree(kdNodes, end[0], numDimensions);
	refIdx_t root = Gpu::buildKdTree(kdNodes, end[0], numDimensions);
	TIMER_STOP (double kdTime);

	// Verify the k-d tree and report the number of KdNodes.
	TIMER_START();
	sint numberOfNodes = Gpu::verifyKdTree(kdNodes, root, numDimensions, numTuples);
	// sint numberOfNodes = kdNodes[root].verifyKdTree( kdNodes, coordinates, numDimensions, 0);
	cout <<  "Number of nodes = " << numberOfNodes << endl;
	TIMER_STOP (double verifyTime);

	cout << "totalTime = " << fixed << setprecision(4) << initTime + sortTime + kdTime + verifyTime
			<< "  initTime = " << initTime << "  sortTime + removeDuplicatesTime = " << sortTime
			<< "  kdTime = " << kdTime << "  verifyTime = " << verifyTime << endl << endl;

	// Return the pointer to the root of the k-d tree.
	return &kdNodes[root];
}

// Define a comparator for the priority queue
struct QueueComparator {
    bool operator() (const pair<double, const KdNode*>& lhs, const pair<double, const KdNode*>& rhs) const {
        return lhs.first < rhs.first;
    }
};

// Define a function to calculate the distance between two points in 3D space
double euclideanDistance(KdNode* p1, const KdCoord* query, const KdCoord coords[], const sint dim) {
    double dx = coords[p1->tuple*dim + 0] - query[0];
    double dy = coords[p1->tuple*dim + 1] - query[1];
    double dz = coords[p1->tuple*dim + 2] - query[2];
    return (dx*dx + dy*dy + dz*dz);
}

/*
 * The csearchKdTree function is a resursive CPU based function which finds the nearest neighbours 
 * of a query points and stores the results in a priority queu. It requires a loop later
 * , after its call has been finished to write the obtained points to the results array that will be 
 * sent to fortran
 *
 * calling parameters:
 *
 * node - a reference of the node of the kdtree which is stored in the kdNodes array
 * kdNodes - pointer to the start of the array holding the KD Nodes. This array is copied from the GPU to CPU
 * coords - pointer to the start of the array holding the coordinates. This array is on the CPU
 * query - pointer to the array holding the query point
 * numResults - number of nearest neighbour required
 * dim - the number of dimensions
 * i - variable to keep track of the depth/axis of comparison
 * pq - pritority queue
 *
 * returns: a KdNode pointer to the root of the k-d tree
 */
void csearchKdTree(KdNode *node, KdNode kdNodes[], const KdCoord coords[], const KdCoord* query, const sint numResults,
		const sint dim, sint i, priority_queue<pair<double, KdNode*>, vector<pair<double, KdNode*>>, QueueComparator> &pq) {

	double curr_dist = euclideanDistance(node, query, coords, dim);

	if (pq.size() < numResults) {
		pq.push({curr_dist, node});
	}
	else if (curr_dist < pq.top().first) {
		pq.pop();
        pq.push({curr_dist, node});
	}

	double perp_;

	perp_ = coords[node->tuple*dim + i] - query[i];

	i = (i+1) % 3;

	if (pow(perp_,2)<pq.top().first) {

		if (node->ltChild!=-1)	{
			csearchKdTree(&kdNodes[node->ltChild], kdNodes, coords, query, numResults, dim, i, pq);
		}
		
		if (node->gtChild!=-1)	{
			csearchKdTree(&kdNodes[node->gtChild], kdNodes, coords, query, numResults, dim, i, pq);
		}
    }

	else {
        if (perp_<0) {
			if (node->gtChild!=-1)	{
            	csearchKdTree(&kdNodes[node->gtChild], kdNodes, coords, query, numResults, dim, i, pq);
			}	
        }
        else {
			if (node->ltChild!=-1)	{
				csearchKdTree(&kdNodes[node->ltChild], kdNodes, coords, query, numResults, dim, i, pq);
			}
		}
    }

    return;

}

/*
 * call the search the k-d tree  function and find the numResults number KdNodes 
 * that lie closest to the query point. It calls the GPU functions to search
 * closest points on gpu itseld
 * calling parameters:
 *
 * root - index of root node of the kdtree in the Kdnode array
 * query - pointer to the query point tuple
 * numResult - number of k nn rewuired
 * dim - the number of dimensions
 * pq - array that will hold the results
 */

/*
 * Print one tuple.
 *
 * calling parameters:
 *
 * tuple - the tuple to print
 * dim - the number of dimensions
 */
void KdNode::printTuple(const KdCoord* tuple, const sint dim)
{
	cout << "(" << tuple[0] << ","; //changed tuple[dim] from citation doce to tuple[0]
	for (sint i=1; i<dim-1; i++) cout << tuple[i] << ",";
	cout << tuple[dim-1] << ")";
}

/*
 * Print the k-d tree "sideways" with the root at the ltChild.
 *
 * calling parameters:
 *
 * dim - the number of dimensions
 * depth - the depth in the k-d tree
 */
void KdNode::printKdTree(KdNode kdNodes[], const KdCoord coords[], const sint dim, const sint depth) const
{
	if (gtChild != -1) {
		kdNodes[gtChild].printKdTree(kdNodes, coords, dim, depth+1);
	}
	for (sint i=0; i<depth; i++) cout << "       ";
	printTuple(coords+tuple*dim, dim);
	cout << endl;
	if (ltChild != -1) {
		kdNodes[ltChild].printKdTree(kdNodes, coords, dim, depth+1);
	}
}









void cprintResultwhole(pair_coord_dist* gpq, sint* numQuerys, const sint* numResults, const KdCoord* coords, const sint* dim, float* mltip) {

	for(int q=0;q<*numQuerys; q++){
		printf("\n");
		pair_coord_dist* pq = &gpq[q*(*numResults)];
		for(int i=0; i<*numResults;i++){
			for(int j=0;j<*dim;j++){
				printf(" %f ", coords[pq[i].tpl*(*dim)+j]/(*mltip));
			}
			printf("\n");
		}
	}

}

void cdistCalc(pair_coord_dist* nd_pair, const KdCoord* query, const KdCoord* coords, const sint* dim, float* mltip) {

	double dx = coords[(nd_pair->tpl)*(*dim) + 0]/(*mltip) - query[0];
    double dy = coords[(nd_pair->tpl)*(*dim) + 1]/(*mltip) - query[1];
    double dz = coords[(nd_pair->tpl)*(*dim) + 2]/(*mltip) - query[2];

	nd_pair->dist = (dx*dx + dy*dy + dz*dz);

}

void cinsertSort(pair_coord_dist* pq, sint* counter, const sint* numResults, pair_coord_dist* nd_pair) {

	if (*counter<(*numResults)) {
		pq[*counter]=*nd_pair;
		*counter = (*counter) +1;
		for (int i=1; i<(*counter); i++) {
			pair_coord_dist key = pq[i];
			int j=i-1;

			while (j >= 0 && pq[j].dist > key.dist) {
				pq[j+1] = pq[j];
				j=j-1;
			}

			pq[j+1] = key;
		}
	}
	else if (nd_pair->dist < pq[(*numResults)-1].dist) {
		pq[(*numResults)-1] = *nd_pair;
		for (int i=1; i<(*numResults); i++) {
			pair_coord_dist key = pq[i];
			int j=i-1;

			while (j >= 0 && pq[j].dist > key.dist) {
				pq[j+1] = pq[j];
				j=j-1;
			}

			pq[j+1] = key;
		}
	}

}

void cinsertInList(KdNode* node, litem** ch, const sint* dim) {
	//litem *li = (litem*) malloc(sizeof(litem));
	litem* li;
	li = (litem*)malloc(sizeof(litem));
	if(li==nullptr) printf("ERRORR!!!!!\n\n");
    
    li->data = node;
    li->axis = (ch[0]->axis+1)%(*dim);
    li->next = nullptr;
 
    ch[1]->next = li;
    ch[1] = li;
}

void cIterSearchKdTree(KdNode *root, KdNode kdNodes[], const KdCoord coords[], const KdCoord* query, const sint* numResults,
		const sint* dim, sint* counter, float* mltip, sint* numQuerys, pair_coord_dist* gpq) {

	int tid;
	for(tid=0; tid<*numQuerys; tid++) {
		if(tid<*numQuerys) {
			
			const KdCoord* qu = &query[tid*(*dim)];
			pair_coord_dist* pq = &gpq[tid*(*numResults)];
			sint dcounter=0;
			
			litem *ls;
			ls = (litem*)malloc(sizeof(litem));
			ls->data = root;
			ls->axis=0;
			ls->next=nullptr;
			litem* ch[2];
			ch[0] = ls;
			ch[1] = ls;

			//printf("%d\n", tid);

			printf("Query Point %d :: %f %f %f\n", tid, qu[0], qu[1], qu[2]);

			int numIter=0;
			
			//printf("%d\n", tid);
			
			while(ch[0]!=ch[1] || ch[0]->data==root) {
				
				//if(ch[0]->data==nullptr) {
					//ch[0] = ch[0]->next;
					//continue;
				//}
				numIter++;
				//printf("tid: %d\n",tid);
				//printf("%f %f %f \n", qu[0], qu[1], qu[2]);
				pair_coord_dist nd_pair;

				//printf("Assigned ndcoord\n");
				nd_pair.tpl = ch[0]->data->tuple;

				//printf("Assigned ndcoord tpl\n");
				cdistCalc(&nd_pair, qu, coords, dim, mltip);
				//printf("Calculated Distance\n");

				cinsertSort(pq, &dcounter, numResults, &nd_pair);
				//printf("Inserted the pair in pq\n");

				double perp_ = coords[ch[0]->data->tuple*(*dim) + ch[0]->axis]/(*mltip) - qu[ch[0]->axis];

				if (pow(perp_,2)<pq[dcounter-1].dist){

					if (ch[0]->data->ltChild!=-1)
					{
						
						//printf("Inside 1st if %d %d\n", ch[0]->data->ltChild, sizeof(kdNodes));
						cinsertInList(&kdNodes[ch[0]->data->ltChild], ch, dim);
						//printf("Exit 1st if\n");
					}
					
					if (ch[0]->data->gtChild!=-1)
					{
						//printf("Inside 2ns if\n");
						cinsertInList(&kdNodes[ch[0]->data->gtChild], ch, dim);
						//printf("Exit 2nd if\n");
					}
				}

				else 
				{
					if (perp_<0) 
					{
						if (ch[0]->data->gtChild!=-1)
						{
							//printf("Inside 3rd if\n");
							cinsertInList(&kdNodes[ch[0]->data->gtChild], ch, dim);
							//printf("Exit 3rd if\n");
						}	
					}
					else 
					{
						if (ch[0]->data->ltChild!=-1)
						{
							//printf("Inside 4th if\n");
							cinsertInList(&kdNodes[ch[0]->data->ltChild], ch, dim);
							//printf("Exit 4th if\n");
						}
					}
				}
				litem* temp = ch[0];
				ch[0] = ch[0]->next;
				free(temp);
				//printf("changed the current pointer\n");

			}
			

			//printResult<<<1,1>>>(pq, numResults, coords, dim, mltip);
			printf("tid: %d :::: Iter:%d\n\n",tid, numIter);

		}
	
	}
	

}

void cSearchKdtreeCPU(KdCoord* coordinates, KdNode kdNodes[], refIdx_t root, const KdCoord* query, const sint numResults, 
	const sint dim, KdCoord* results, sint numQuerys, pair_coord_dist** pqRefs, float mltip) {

		//Declare variables pointers for device
	float* d_mltip = &mltip;
	const sint* d_numResults= &numResults;
	const sint* d_dim = &dim;
	
	pair_coord_dist* d_pq;
	d_pq = (pair_coord_dist*) malloc(sizeof(pair_coord_dist)*numQuerys*numResults);

	sint counter = 0;
	sint*d_counter = &counter;
	
	sint* d_numQuerys = &numQuerys;
	pair_coord_dist** d_pqRefs;
	litem* temp;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start);
	

	cIterSearchKdTree(kdNodes+root, kdNodes, coordinates, query, d_numResults, d_dim, d_counter, d_mltip, d_numQuerys, d_pq);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float searchSeconds = 0;
	cudaEventElapsedTime(&searchSeconds, start, stop);
	searchSeconds = searchSeconds/1000;
	
	cudaEventRecord(start);
	cprintResultwhole(d_pq, d_numQuerys, d_numResults, coordinates, d_dim, d_mltip);
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float printSeconds;
	cudaEventElapsedTime(&printSeconds, start, stop);
	printSeconds = printSeconds/1000;
	printf("cuIterSearchKdTree Kernel execution time: %f s\n", searchSeconds);
	printf("printResultwhole Kernel execution time: %f s\n", printSeconds);


	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	//checkCudaErrors(cudaDeviceSynchronize());
	pqRefs[0] = d_pq;
	}

//Declaration for the custom function so that it can be called by fortran
extern "C" {
    void custom_funct(KdCoord* coordinates, sint numPoints, sint numQuerys, KdCoord* query, KdCoord* results, sint numResults, sint* rootIdx);
}

/*
 * A custom function which acts as the main function 
 * This is the function that is called by fortran
 * calling parameters:
 *
 * 
 * 
*/
void custom_funct( KdCoord* coordinates, sint numPoints, sint numQuerys, KdCoord* query, KdCoord* results, sint numResults, sint* rootIdx)
{
	sint extraPoints = 100;
	sint numDimensions = 3;
	sint numThreads = 512;
	sint numBlocks = 16;
	sint maximumNumberOfNodesToPrint = 5;
	
	//sint  i = maximumNumberOfNodesToPrint + numDimensions + extraPoints;
	
	
	//GPU
	Gpu::gpuSetup(2, numThreads,numBlocks,numDimensions);
	if (Gpu::getNumThreads() == 0 || Gpu::getNumBlocks() == 0) {
		cout << "KdNode Tree cannot be built with " << numThreads << " threads or " << numBlocks << " blocks." << endl;
		exit(1);
	}
	cout << "Points = " << numPoints << " dimensions = " << numDimensions << ", threads = " << numThreads << ", blocks = " << numBlocks << endl;

	// Create the k-d tree.  First copy the data to a tuple in its kdNode.
	// also null out the gt and lt references
	// create and initialize the kdNodes array
	KdNode *kdNodes = new KdNode[numPoints];
	if (kdNodes == NULL) {
		printf("Can't allocate %d kdNodes\n", numPoints);
		exit (1);
	}
	
	//MAKE KDTREE
	KdNode *root = KdNode::createKdTree(kdNodes, coordinates, numDimensions, numPoints);
	
	cout << endl;

	cout<<"KD TREE MADE and stored on GPU"<<endl;


	// rootIdx stores the index of the root node on the kdnodes array
	// this value is later passed to the searchKdTree function  
	*rootIdx = root-kdNodes;

	// TIMER_DECLARATION();

	// TIMER_START();
	// Gpu::getKdTreeResults( kdNodes,  coordinates, numPoints, numDimensions);
	// TIMER_STOP(double copyTreeTime);

	// TIMER_START();
	// // SEARCH FOR Nearest Neighbour
	// pair_coord_dist* pqRefs[numQuerys];
	// cSearchKdtreeCPU(coordinates, kdNodes, *rootIdx, query, numResults, numDimensions, results, numQuerys, pqRefs, 100000);
	// TIMER_STOP(double searchTime);

	// printf("time to copy tree %f\n", copyTreeTime);
	// cout << "Total Search Time = " << searchTime << endl;
	
}


//Declaration for the custom function so that it can be called by fortran
extern "C" {
    void search_funct(KdCoord* coordinates, sint numPoints, sint numDimensions, sint numQuerys, KdCoord* query, KdCoord* results, sint numResults, sint rootIdx, float mltip, sint* gindices, double* dists);
}

/*
 * A custom function which acts as the main function 
 * This is the function that is called by fortran
 * calling parameters:
 *
 * 
 * 
*/
void search_funct( KdCoord* coordinates, sint numPoints, sint numDimensions, sint numQuerys, KdCoord* query, KdCoord* results, sint numResults, sint rootIdx, float mltip, sint* gindices, double* dists)
{
	// rootIdx stores the index of the root node on the kdnodes array
	// this value is later passed to the searchKdTree function  
	KdNode *kdNodes = new KdNode[numPoints];
	if (kdNodes == NULL) {
		printf("Can't allocate %d kdNodes\n", numPoints);
		exit (1);
	}

	TIMER_DECLARATION();

	TIMER_START();
	Gpu::getKdTreeResults( kdNodes,  coordinates, numPoints, numDimensions);
	TIMER_STOP(double copyTreeTime);

	printf("time to copy tree %f\n", copyTreeTime);
	printf("HERE %d\n", kdNodes[rootIdx].tuple);

	TIMER_START();
	// SEARCH FOR Nearest Neighbour
	pair_coord_dist* pqRefs[numQuerys];
	cSearchKdtreeCPU(coordinates, kdNodes, rootIdx, query, numResults, numDimensions, results, numQuerys, pqRefs, mltip);
	TIMER_STOP(double searchTime);
	//cout << "Total Search Time = " << searchTime << endl;


	// TIMER_DECLARATION();

	// TIMER_START();
	// // SEARCH FOR Nearest Neighbour
	// pair_coord_dist* pqRefs[numQuerys];
	// Gpu::searchKdTree(coordinates, rootIdx, query, numResults, numDimensions, results, numQuerys, pqRefs, mltip);
	// TIMER_STOP(double searchTime);
	// //cout << "Total Search Time = " << searchTime << endl;

	// TIMER_START();
	// Gpu::getSearchResults(pqRefs, coordinates, numResults, numDimensions, results, numQuerys, mltip, gindices, dists);
	// TIMER_STOP(double getSearchResultsTime);
	// printf("Total time taken to get the results back cudaMemcpy (including cudaMemcpy (Gpu::getSearchResultsGPU) time): %f s\n", getSearchResultsTime);
	

	
	
}
