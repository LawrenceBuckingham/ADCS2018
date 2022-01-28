// ------------------------------------------------------------------
// Copyright (C) 2018 Lawrence Buckingham.
//
// This file is part of the supplementary material which accompanies 
// the paper:
//	Lawrence Buckingham, Shlomo Geva, and James M. Hogan. 2018.
//	Protein database search using compressed k-mer vocabularies.
//	In 23rd Australasian Document Computing Symposium (ADCS '18), 
//	December 11--12, 2018, Dunedin, New Zealand
//	https://doi.org/10.1145/3291992.3291997
//
// This file is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3, or (at your 
// option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License 
// along with this program; see the file "license.txt".  If not, see
// <http://www.gnu.org/licenses/>.
// ------------------------------------------------------------------

#include "AlphabetHelper.hpp"
#include "Args.hpp"
#include "Assert.hpp"
#include "Delegates.hpp"
#include "KmerDistanceCache.hpp"
#include "FastaSequence.hpp"
#include "HBRandom.hpp"
#include "Kmer.hpp"
#include "KmerCluster.hpp"
#include "TestFramework.h"
#include "OmpTimer.h"
#include "KmerClusterPrototype.hpp"
#include "FileUtil.hpp"

#include <bitset>
#include <cstdio>
#include <omp.h>
#include <set>

using namespace QutBio;
using namespace std;

// Address of singleton argument table.
Args *arguments;

struct AAClust {
public:
	static int Run() {
		bool ok = true;

		OMP_TIMER_DECLARE(loadTime);
		OMP_TIMER_DECLARE(clusterTime);

		string protoIn;
		string protoOut;
		string fastaFile;
		int numThreads = 7;
		int wordLength = 32;
		int threshold;
		int  seed;
		int idIndex;
		string clusterOut;
		int increment;
		int clusterMode = 1;

		if (arguments->IsDefined("help")) {
			vector<string> text{
				"AAClust: Greedy clustering of Amino Acid kmers by substitution matrix.",
				"--help	Gets this text.",
				"--fastaFile	Required. A list of one or more file paths. Each file will be parsed as a FASTA file which contains DNA sequences to be clustered.",
				"--idIndex	Required. The 0-origin position of the sequence ID field in the pipe-separated definition line.",
				"--generateEdges	Optional boolean, default value false. If true, edges for a multiple alignment will be generated.",
				"--clusterOut	Required. The name the output file produced by the program.",
				"--increment	Required. The number of new clusters to add on each pass. Make this smaller to minimise the chance of a prototype belonging to a cluster whose centroid is outside its basin of attraction.",
				"--threshold	Required. Threshold for assignment of points to clusters. Distance less than or equal to the threshold corresponds to cluster membership.",
				"--numThreads	Optional; default value = 7. The number of OpenMP threads to use in parallel regions.",
				"--wordLength	Optional; default value = 32. The word length used for kmer tiling.",
				"--seed		Required. The random number seed.",
				"--merge	Required. The merge mode used to combine overlapping HSKP (Highly Significant Kmer Pairs). Valid values are:",
				"		none: Do not merge;",
				"		consecutive: Merge new edge onto existing edge when both endpoints of the new edge are exact continuations of the previous edge.",
				"		overlapping: Merge new edge onto existing edge when both endpoints of the new edge are within the current extent of the corresponding intervals represented by the existing edge. Some kmers covered by an edge built this way may not be HSKPs.",
				"--matrixId	Optional, default = 62, but you need either --matrixId or --matrixFile. BLOSUM Matrix ID, one of { 35, 40, 45, 50, 62, 80, 100 }. This is ignored if a custom similarity matrix file is specified.",
				"--matrixFile	Optional, but you need either --matrixId or --matrixFile. File name for custom similarity matrix. Use this to specify some matrix other than BLOSUM, or if a custom alphabet is in use.",
				"--isCaseSensitive	Optional, default = true. Should symbols be treated as case-sensitive.",
				"--clusterMode	Optional [1, 2], default = 1. Experimental clustering mode.",
				"		1: Use Insertion-sort inspired modification to reduce worst case complexity by average factor of at least 2.",
				"		2: Use banded version of 1 to partition work to threads ahead of time (which in the end slows things down).",
			};

			for (auto s : text) {
				cerr << s << "\n";
			}

			return 0;
		}

		if (!arguments->Get("protoOut", protoOut)) {
			cerr << arguments->ProgName() << ": error - required argument '--protoOut' not supplied.\n";
			ok = false;
		}

		if (!arguments->Get("fastaFile", fastaFile)) {
			cerr << arguments->ProgName() << ": error - required argument '--fastaFile' not supplied.\n";
			ok = false;
		}

		if (!arguments->Get("idIndex", idIndex)) {
			cerr << arguments->ProgName() << ": error - required argument '--idIndex' not supplied.\n";
			ok = false;
		}

		if (!arguments->Get("seed", seed)) {
			cerr << arguments->ProgName() << ": error - required argument '--seed' not supplied.\n";
			ok = false;
		}

		if (!arguments->Get("threshold", threshold)) {
			cerr << arguments->ProgName() << ": error - required argument '--threshold' not provided.\n";
			ok = false;
		}

		if (!arguments->Get("increment", increment)) {
			cerr << arguments->ProgName() << ": error - required argument '--increment' not provided.\n";
			ok = false;
		}

		if (!arguments->Get("clusterOut", clusterOut)) {
			cerr << arguments->ProgName() << ": error - required argument '--clusterOut' not provided.\n";
			ok = false;
		}

		int matrixId;

		if (arguments->IsDefined("matrixId")) {
			if (!(arguments->Get("matrixId", matrixId))) {
				cerr << arguments->ProgName() << ": error - argument 'matrixId' not valid." << endl;
				ok = false;
			}

			vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

			bool found = false;

			for (auto x : matrices) {
				if (x == matrixId) { found = true; }
			}

			if (!found) {
				cerr << arguments->ProgName() << ": error - matrix id not recognised." << endl;
				ok = false;
			}
		}

		if (!arguments->Get("numThreads", numThreads)) {
			cerr << arguments->ProgName() << ": note  - optional argument '--numThreads' not set; running with default value " << numThreads << ".\n";
		}

		if (!arguments->Get("wordLength", wordLength)) {
			cerr << arguments->ProgName() << ": note  - optional argument '--wordLength' not set; running with default value " << wordLength << ".\n";
		}

		if (!arguments->Get("protoIn", protoIn)) {
			cerr << arguments->ProgName() << ": note  - optional argument '--protoIn' not supplied.\n";
			cerr << arguments->ProgName() << ": note  - (new prototypes will be generated from kmers in the current input dataset)\n";
		}

		if (arguments->Get("clusterMode", clusterMode)) {
			if (clusterMode < 1 || clusterMode > 2) {
				cerr << arguments->ProgName() << ": error - invalid value for '--clusterMode'.\n";
				ok = false;
			}
		}

		if (!ok) {
			cerr << "Invalid command line arguments supplied. For help, run: AAClust --help\n";
			return 1;
		}

		string matrixFile;
		DistanceType * distanceType = DistanceType::BlosumDistance();

		if (arguments->IsDefined("matrixFile")) {
			arguments->Get("matrixFile", matrixFile);
			distanceType = DistanceType::Custom();
			matrixId = -1;
		}

		bool isCaseSensitive = true;

		if (arguments->IsDefined("isCaseSensitive")) {
			if (!arguments->Get("isCaseSensitive", isCaseSensitive)) {
				cerr << arguments->ProgName() << ": error - Invalid data for argument 'isCaseSensitive'." << "\n";
				ok = false;
			}
		}

		SimilarityMatrix * matrix = SimilarityMatrix::GetMatrix(distanceType, matrixId, matrixFile, isCaseSensitive);

		if (!matrix) {
			cerr << arguments->ProgName() << ": error - unable to construct similarity matrix.\n";
			cerr << arguments->ProgName() << ": error - you need to supply either matrixId or matrixFile arguments.\n";
			ok = false;
		}

		if (ok == false || !matrix) {
			cerr << arguments->ProgName() << ": error - unable to construct similarity matrix. For help, run: AAClust --help\n";
			return 1;
		}

		Alphabet * alphabet = new Alphabet(matrix);
		BlosumDifferenceFunction rawDistanceFunction(matrix);
		KmerDistanceCache2 distanceFunction(alphabet, &rawDistanceFunction);
		size_t charsPerWord = distanceFunction.CharsPerWord();
		UniformRealRandom rand(seed);

		using DistanceFunction = KmerDistanceCache2;
		using Cluster = KmerCluster<DistanceFunction, Kmer>;

		vector<Cluster *> clusters;

		omp_set_num_threads(numThreads);

		OMP_TIMER_START(loadTime);
		PointerList<EncodedFastaSequence> protos;

		if (protoIn.length() > 0) {
			fstream protoStream(protoIn);
			EncodedFastaSequence::ReadSequences(protos, protoStream, 0, -1, alphabet, wordLength, charsPerWord, 'x', KmerClusterPrototype::DefaultFactory);

			if ( ! isCaseSensitive ) {
				for ( auto p: protos ) {
					String::ToLowerInPlace( p->Sequence() );
				}
			}
			
			Cluster::InitialiseClusters(protos, wordLength, distanceFunction, clusters);
			(cerr << "AAClust: " << protos.Length() << " prototypes loaded.\n").flush();
		}

		PointerList<EncodedFastaSequence> db;

		{
			fstream fasta(fastaFile);
			EncodedFastaSequence::ReadSequences(db, fasta, idIndex, -1, alphabet, wordLength, distanceFunction.CharsPerWord(), 'x', EncodedFastaSequence::DefaultFactory);

			if ( ! isCaseSensitive ) {
				for ( auto p: db ) {
					String::ToLowerInPlace( p->Sequence() );
				}
			}
		}

		cerr << "AAClust: " << db.Length() << " sequences loaded.\n";

		KmerIndex kmerIndex(db.Items(), wordLength);

		OMP_TIMER_END(loadTime);
		OMP_TIMER_START(clusterTime);

		auto createPrototype = [=, &protos](Kmer *kmer) {
			KmerClusterPrototype * protoSeq = new KmerClusterPrototype(kmer->Word(), alphabet, wordLength, charsPerWord);
			protos.Add([protoSeq]() { return protoSeq; });
			return protoSeq;
		};

		if (clusterMode == 2) {
			// Banded version: doesn't work as fast as default.
			Cluster::DoExhaustiveIncrementalClustering(
				kmerIndex,
				wordLength,
				threshold,
				alphabet->Size(),
				distanceFunction,
				rand,
				increment,
				createPrototype,
				clusters,
				numThreads
			);
		}
		else {
			// Default: uses the "first-fit" criterion to assign k-mers to 
			// cluster, and remove from consideration.
			Cluster::DoExhaustiveIncrementalClustering(
				kmerIndex,
				wordLength,
				threshold,
				alphabet->Size(),
				distanceFunction,
				rand,
				increment,
				createPrototype,
				clusters
			);
		}
#if REMOVE_TINY_CLUSTERS
		int i;

		for (i = clusters.size() - 1; i >= 0 && (clusters[i]->kmers.size() < 2; i--) {
			delete clusters[i];
				clusters[i] = 0;
		}

		clusters.resize(i + 1);
#endif
			OMP_TIMER_END(clusterTime);

			// Update prototype sizes.
		{
			for (auto cluster : clusters) {
				auto proto = (pKmerClusterPrototype)cluster->prototype.Sequence();
				proto->Size(proto->Size() + cluster->InstanceCount());
			}
		}

		// Save the prototypes.
		{
			ofstream protoFile(protoOut);

			for (auto proto_ : protos) {
				auto proto = (pKmerClusterPrototype)proto_;

				if (proto->Size() > 0) {
					protoFile << *proto;
				}
			}
		}

		ofstream cOut( clusterOut );
		for ( auto c: clusters) cOut << (*c);

		cerr << "Elapsed time loading: " << OMP_TIMER(loadTime) << "\n";
		cerr << "Elapsed time clustering: " << OMP_TIMER(clusterTime) << "\n";

		return 0;
	}
};

mutex QutBio::DistanceType::m;

int main(int argc, char *argv[]) {
	try {
		Args args(argc, argv);

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = AAClust::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << (end_time - start_time) << "s" << endl;

		return retCode;
	}
	catch (Exception ex) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}


