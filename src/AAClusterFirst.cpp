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
#include "KmerCodebook.hpp"
#include "TestFramework.h"
#include "OmpTimer.h"
#include "FileUtil.hpp"

#include <bitset>
#include <cstdio>
#include <omp.h>
#include <set>

using namespace QutBio;
using namespace std;

// Address of singleton argument table.
Args *arguments;

struct AAClusterFirst {

	using DistanceFunction = KmerDistanceCache2;
	using Cluster = KmerCluster<DistanceFunction, Kmer>;
	using Codebook = KmerCodebook<DistanceFunction, Kmer>;

public:
	static int Run() {
		bool ok = true;

		string fastaFile;
		string clusterIn;
		string protoIn;
		string clusterOut;
		string protoOut;
		size_t numThreads = 7;
		size_t wordLength = 32;
		int idIndex = 0;
		size_t numClusters = 0;

		if ( arguments->IsDefined( "help" ) ) {
			vector<string> text{
				"AAClusterFirst: Gets the ${numClusters} largest clusters from a kmer codebook.",
				"--help	Gets this text.",
				"--fastaFile    Required. A list of one or more file paths. Each file will be parsed as a FASTA file which contains DNA sequences that have been clustered.",
				"--clusterIn    Required. The name of a file that contains a list of k-mer clusters.",
				"--protoIn      Required. The name of a file containing the prototypes."
				"--numClusters  Required. The number of clusters to select fromt he codebook.",
				"--clusterOut   Required. The name of the cluster output file.",
				"--protoOut     Required. The name of the prototype output file.",
				"--idIndex      Required. The 0-origin position of the sequence ID field in the pipe-separated definition line.",
				"--numThreads   Optional; default value = 7. The number of OpenMP threads to use in parallel regions.",
				"--wordLength   Optional; default value = 32. The word length used for kmer tiling.",
				"--matrixId     Optional, default = 62. BLOSUM Matrix ID, one of { 35, 40, 45, 50, 62, 80, 100 }. This is ignored if a custom similarity matrix file is specified.",
				"--matrixFile   Optional. File name for custom similarity matrix. Use this to specify some matrix other than BLOSUM, or if a custom alphabet is in use.",
			};

			for ( auto s : text ) {
				cerr << s << "\n";
			}

			return 0;
		}

		if ( !arguments->Get( "fastaFile", fastaFile ) ) {
			cerr << arguments->ProgName() << ": error - required argument '--fastaFile' not supplied.\n";
			ok = false;
		}

		if ( !arguments->Get( "clusterIn", clusterIn ) ) {
			cerr << arguments->ProgName() << ": error - required argument '--clusterIn' not supplied.\n";
			ok = false;
		}

		if ( !arguments->Get( "protoIn", protoIn ) ) {
			cerr << arguments->ProgName() << ": error - required argument '--protoIn' not supplied.\n";
			ok = false;
		}

		if ( !arguments->Get( "idIndex", idIndex ) ) {
			cerr << arguments->ProgName() << ": error - required argument '--idIndex' not supplied.\n";
			ok = false;
		}

		if ( !arguments->Get( "numClusters", numClusters ) ) {
			cerr << arguments->ProgName() << ": error - required argument '--numClusters' not supplied.\n";
			ok = false;
		}

		if ( !arguments->Get( "numThreads", numThreads ) ) {
			cerr << arguments->ProgName() << ": note - optional argument '--numThreads' not set; running with default value " << numThreads << ".\n";
		}

		if ( !arguments->Get( "wordLength", wordLength ) ) {
			cerr << arguments->ProgName() << ": note - optional argument '--wordLength' not set; running with default value " << wordLength << ".\n";
		}

		if ( !arguments->Get( "clusterOut", clusterOut ) ) {
			cerr << arguments->ProgName() << ": note - required argument '--clusterOut' not provided.\n";
			ok = false;
		}

		if ( !arguments->Get( "protoOut", protoOut ) ) {
			cerr << arguments->ProgName() << ": note - required argument '--protoOut' not provided.\n";
			ok = false;
		}

		int matrixId;

		if ( arguments->IsDefined( "matrixId" ) ) {
			if ( !( arguments->Get( "matrixId", matrixId ) ) ) {
				cerr << arguments->ProgName() << ": error - argument 'matrixId' not valid." << endl;
				ok = false;
			}

			vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

			bool found = false;

			for ( auto x : matrices ) {
				if ( x == matrixId ) { found = true; }
			}

			if ( !found ) {
				cerr << arguments->ProgName() << ": error - matrix id not recognised." << endl;
				ok = false;
			}
		}

		if ( !ok ) {
			cerr << "Invalid command line arguments supplied. For help, run: " << arguments->ProgName() << " --help\n";
			return 1;
		}

		string matrixFile;
		DistanceType * distanceType = DistanceType::BlosumDistance();

		if ( arguments->IsDefined( "matrixFile" ) ) {
			arguments->Get( "matrixFile", matrixFile );
			distanceType = DistanceType::Custom();
			matrixId = -1;
		}

		bool isCaseSensitive = false;

		if ( arguments->IsDefined( "isCaseSensitive" ) ) {
			if ( !arguments->Get( "isCaseSensitive", isCaseSensitive ) ) {
				cerr << "Invalid data for argument 'isCaseSensitive'." << "\n";
				ok = false;
			}
		}

		SimilarityMatrix * matrix = SimilarityMatrix::GetMatrix( distanceType, matrixId, matrixFile, isCaseSensitive );

		if ( ok == false || !matrix ) {
			cerr << "Unable to construct similarity matrix. For help, run: AAClusterFirst --help\n";
			return 1;
		}

		if ( clusterOut == clusterIn || clusterOut == fastaFile || clusterOut == protoIn ) {
			cerr << "AASample: Output file " << clusterOut << " will overwrite one of your input files.\n";
			return 1;
		}

		if ( protoOut == clusterIn || protoOut == fastaFile || protoOut == protoIn ) {
			cerr << "AASample: Output file " << protoOut << " will overwrite one of your input files.\n";
			return 1;
		}

		Alphabet * alphabet = new Alphabet( matrix );
		BlosumDifferenceFunction rawDistanceFunction( matrix );
		using DistanceFunction = KmerDistanceCache2;
		DistanceFunction distanceFunction( alphabet, &rawDistanceFunction );

		omp_set_num_threads( numThreads );

		PointerList<EncodedFastaSequence> db;

		{
			fstream fasta( fastaFile );
			EncodedFastaSequence::ReadSequences(
				db,
				fasta,
				idIndex,
				-1,
				alphabet,
				wordLength,
				distanceFunction.CharsPerWord(),
				'x',
				EncodedFastaSequence::DefaultFactory
			);
		}

		cerr << "AAClusterFirst: " << db.Length() << " sequences loaded.\n";

		EncodedFastaSequence::Index seqIndex( db.Items() );
		KmerIndex kmerIndex( db.Items(), wordLength );

		PointerList<EncodedFastaSequence> protos;
		ifstream protoStream( protoIn );
		EncodedFastaSequence::ReadSequences(
			protos,
			protoStream,
			0,
			-1,
			alphabet,
			wordLength,
			distanceFunction.CharsPerWord(),
			'x',
			KmerClusterPrototype::DefaultFactory
		);
		EncodedFastaSequence::Index protoIndex( protos.Items() );

		//for (auto & i : protoIndex) {
		//	cerr << i.first << "\n";
		//	cerr << (*i.second) << "\n";
		//}

		Codebook * codebook = 0;

		FILE * f = fopen( clusterIn.c_str(), "r" );
		if ( f ) {
			codebook = new Codebook( alphabet, distanceFunction, distanceFunction.CharsPerWord(), wordLength, seqIndex, protoIndex, kmerIndex, f );
			fclose( f );

			if ( codebook->Size() == 0 ) {
				cerr << "Cluster dataset contains no entries; run terminated.\n";
				exit( 1 );
			}
		}
		else {
			( cerr << "Cluster dataset " << clusterIn << " cannot be opened for reading.\n" ).flush();
			exit( 1 );
		}

		using pCluster = Cluster * ;
		vector<pCluster> &clusters{ codebook->Codebook() };

		auto descendingClusterSize = []( const pCluster & lhs, const pCluster & rhs ) {
			return lhs->InstanceCount() > rhs->InstanceCount();
		};

		std::sort( clusters.begin(), clusters.end(), descendingClusterSize );

		( cerr << "selecting largest " << numClusters << " clusters from " << clusters.size() << "\n" ).flush();

		vector<pCluster> clusterSubset;
		vector<KmerClusterPrototype *> protoSubset;

		for ( auto cluster : clusters ) {
			if ( clusterSubset.size() >= numClusters ) break;

			clusterSubset.push_back( cluster );
			protoSubset.push_back( (KmerClusterPrototype *) cluster->prototype.Sequence() );
		}

		ofstream cOut( clusterOut );
		for ( auto c: clusterSubset ) cOut << (*c);

		ofstream pOut( protoOut );
		for ( auto p: protoSubset ) pOut << (*p);

		return 0;
	}
};

mutex QutBio::DistanceType::m;

int main( int argc, char *argv[] ) {
	try {
		Args args( argc, argv );

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = AAClusterFirst::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << ( end_time - start_time ) << "s" << endl;

		return retCode;
	}
	catch ( Exception ex ) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}


