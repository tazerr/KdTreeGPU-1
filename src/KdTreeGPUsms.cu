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


//Declaration for the custom function so that it can be called by fortran
extern "C" {
    void custom_funct(KdCoord* coordinates, sint numPoints, sint numQuerys, KdCoord* query, KdCoord* results, sint numResults, sint* rootIdx);
}

/*
 * A custom function which acts as the main function 
 * This is the function that is called by fortran
 * calling parameters:
 * coordinates: Array containing the grig points
 * numPoints: Total number pf points in the grid
 * rootIdx: Variable that will hold the index of the root of the kdTree
*/
void custom_funct( KdCoord* coordinates, sint numPoints, sint numQuerys, KdCoord* query, KdCoord* results, sint numResults, sint* rootIdx)
{
	sint numDimensions = 3;
	sint numThreads = 512;
	sint numBlocks = 16;
	
	
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
	
}


//Declaration for the search_funct function so that it can be called by fortran
extern "C" {
    void search_funct(KdCoord* coordinates, sint numPoints, sint numDimensions, sint numQuerys, KdCoord* query, KdCoord* results, sint numResults, sint rootIdx,
		 float mltip, sint* gindices, double* dists, double* u, double* v, double* uinterp, double* vinterp);
}

/*
 * A custom function which acts as the main function 
 * This is the function that is called by fortran
 * calling parameters:
 *
 * 
 * 
*/
void search_funct( KdCoord* coordinates, sint numPoints, sint numDimensions, sint numQuerys, KdCoord* query, KdCoord* results, sint numResults, sint rootIdx, 
	float mltip, sint* gindices, double* dists, double* u, double* v, double* uinter, double* vinter)
{
	
	TIMER_DECLARATION();

	TIMER_START();
	// SEARCH FOR Nearest Neighbour
	double* interpRefs[4];
	Gpu::searchKdTree(coordinates, numPoints, rootIdx, query, numResults, numDimensions, results, numQuerys, interpRefs, mltip);

	//Interpolation
	Gpu::getInterpResults(interpRefs, numPoints, numQuerys, u, v, uinter, vinter);
	Gpu::interpolation(coordinates, numPoints, numDimensions, numQuerys, query, numResults, gindices, dists, mltip, u, v, uinter, vinter);
	
	TIMER_STOP(double searchTime);
	cout << "Total Search+Interpolation Time = " << searchTime << endl;
	
}

