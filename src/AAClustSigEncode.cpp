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
#include "KmerCodebook.hpp"
#include "KmerDistanceCache.hpp"
#include "FastaSequence.hpp"
#include "HBRandom.hpp"
#include "Kmer.hpp"
#include "KmerCluster.hpp"
#include "TestFramework.h"
#include "OmpTimer.h"
#include "BitSet.hpp"
#include "kNearestNeighbours.hpp"
#include "Ranking.hpp"

#include <cstdio>
#include <stdlib.h>
#include <omp.h>
#include <set>
#include <vector>

using namespace QutBio;
using namespace std;

// Singletons.
Args *arguments;
mutex QutBio::DistanceType::m;

struct AAClustSig {
public:
	using DistanceFunction = KmerDistanceCache2;
	using Cluster = KmerCluster<DistanceFunction, Kmer>;
	using pCluster = Cluster * ;

	struct Params {
	public:
		string seqFile;
		string protoFile;
		string outFile;
		size_t numThreads = 7;
		size_t wordLength = 0;
		int idIndex = 0;
		int classIndex = 0;
		bool ok = true;
		DistanceType *distanceType = DistanceType::BlosumDistance();
		SimilarityMatrix *matrix;
		Distance threshold;
		bool assignNearest = false;

		Params() {

			if ( arguments->IsDefined( "help" ) ) {
				vector<string> text{
					"AAClustSigEncode: Generates a sparse binary signature for each input sequence based on",
					"                  kmer proximity to one of a list of prototypes. Bit i is set in the",
					"                  signature iff there exists kmer k in sequence for which dist(k,p_i) <= T",
					"                  (p_i is prototype kmer; T is a real-valued threshold > 0).",
					"",
					"--help         Gets this text.",
					"--seqFile      Required. A file path. The file will be parsed as a FASTA file which contains ",
					"                         amino acid sequences that have been clustered.",
					"--protoFile    Required. The name of a file containing the prototypes."
					"--outFile      Required. The name of a file which will be overwritten with signatures.",
					"--idIndex      Required. The 0-origin position of the sequence ID field in the pipe-separated",
					"                         definition line.",
					"--classIndex   Required. The 0-origin position of the sequence class label field in the pipe-",
					"                         separated definition line.",
					"                         Class labels are a semicolon-separated list of arbitrary strings (no",
					"                         embedded semicolons!)",
					"--wordLength   Required; The word length used for kmer tiling.",
					"--threshold    Required. Positive integer specifying the distance cutoff for assignment of ",
					"                         kmers to clusters. A kmer is considered to be a member of the cluster ",
					"                         if the distance from kmer to cluster centroid is equal to or less than ",
					"                         the threshold distance. The threshold should match that used when the ",
					"                         codebook was constructed.",
					"--numThreads   Optional; default value = 7. The number of OpenMP threads to use in parallel regions.",
					"--matrixId     Optional, default = 62. BLOSUM Matrix ID, one of { 35, 40, 45, 50, 62, 80, 100 }.",
					"                         This is ignored if a custom similarity matrix file is specified.",
					"                         {Why do we need this? The clusters depend on a distanceFunction function for",
					"                         membership, and I currently need to know the distanceFunction function to form",
					"                         a cluster. This may be cleaned up at some point in the future.}",
					"--matrixFile   Optional. File name for custom similarity matrix. Use this to specify some matrix ",
					"                         other than BLOSUM, or if a custom alphabet is in use.",
					"--assignNearest Opt.     Boolean, default = false. Assign k-mers to only one cluster instead of all that fall ",
					"                         within threshold."
				};

				for ( auto s : text ) {
					cerr << s << "\n";
				}
			}

			if ( !arguments->Get( "seqFile", seqFile ) ) {
				cerr << arguments->ProgName() << ": error - required argument '--seqFile' not supplied.\n";
				ok = false;
			}

			if ( !arguments->Get( "protoFile", protoFile ) ) {
				cerr << arguments->ProgName() << ": error - required argument '--protoFile' not "
					"supplied.\n";
				ok = false;
			}

			if ( !arguments->Get( "idIndex", idIndex ) ) {
				cerr << arguments->ProgName() << ": error - required argument '--idIndex' not "
					"supplied.\n";
				ok = false;
			}

			if ( !arguments->Get( "classIndex", classIndex ) ) {
				cerr << arguments->ProgName() << ": error - required argument '--classIndex' not "
					"supplied.\n";
				ok = false;
			}

			if ( !arguments->Get( "numThreads", numThreads ) ) {
				cerr << arguments->ProgName() << ": note - optional argument '--numThreads' not set"
					"; running with default value "
					<< numThreads << ".\n";
			}

			if ( !arguments->Get( "wordLength", wordLength ) ) {
				cerr << arguments->ProgName() << ": error: required argument '--wordLength' not supplied.\n";
				ok = false;
			}

			if ( !arguments->Get( "outFile", outFile ) ) {
				cerr << arguments->ProgName() << ": note - required argument '--outFile' not provided.\n";
				ok = false;
			}

			if ( !arguments->Get( "threshold", threshold ) ) {
				cerr << arguments->ProgName() << ": Error - required argument '--threshold' not provided.\n";
				ok = false;
			}

			if ( arguments->IsDefined( "assignNearest" ) && !arguments->Get( "assignNearest", assignNearest ) ) {
				cerr << arguments->ProgName() << ": Error - invalid boolean data for argument '--assignNearest'.\n";
				ok = false;
			}

			string error;
			if ( !arguments->Get( matrix, error ) ) {
				cerr << error << '\n';
				ok = false;
			}

			if ( outFile == seqFile || outFile == protoFile ) {
				cerr << arguments->ProgName() << ": Output file " << outFile << " will overwrite one of your input files.\n";
				ok = false;
			}
		}
	};

	using DF = KmerDistanceCache2;
	using Codebook = KmerCodebook<DF, Kmer>;

	static int Run() {
		Params parms;

		if ( !parms.ok ) {
			return 1;
		}

		Alphabet *alphabet = new Alphabet( parms.matrix );
		BlosumDifferenceFunction rawDistanceFunction( parms.matrix );
		DistanceFunction distanceFunction( alphabet, &rawDistanceFunction );

		omp_set_num_threads( parms.numThreads );

		PointerList<EncodedFastaSequence> db;
		EncodedFastaSequence::ReadSequences( db, parms.seqFile, parms.idIndex, parms.classIndex, alphabet, parms.wordLength, distanceFunction.CharsPerWord(), alphabet->DefaultSymbol(), EncodedFastaSequence::DefaultFactory );
		cerr << arguments->ProgName() << ": " << db.Length() << " reference sequences loaded from " << parms.seqFile << ".\n";

		PointerList<KmerClusterPrototype> protos;
		EncodedFastaSequence::ReadSequences<KmerClusterPrototype>( protos, parms.protoFile, 0, -1, alphabet, parms.wordLength, distanceFunction.CharsPerWord() );
		cerr << arguments->ProgName() << ": " << protos.Length() << " prototypes loaded from " << parms.protoFile << ".\n";

		OMP_TIMER_DECLARE( encodeDb );
		OMP_TIMER_START( encodeDb );
		Encode( db, protos, distanceFunction, parms.wordLength, parms.threshold, parms.assignNearest, parms.outFile );
		OMP_TIMER_END( encodeDb );

		cerr << "Database encoded in " << OMP_TIMER( encodeDb ) << "s.\n";

		// SaveSignatures(db, parms.outFile);
		return 0;
	}

	static void Encode(
		PointerList<EncodedFastaSequence> &sequences,
		PointerList<KmerClusterPrototype> &protos,
		DistanceFunction &distanceFunction,
		uint K,
		Distance threshold,
		bool assignNearest,
		string &outFile //
	) {
		if ( assignNearest ) {
			EncodeNearest( sequences, protos, distanceFunction, K, threshold, outFile );
		}
		else {
			EncodeAny( sequences, protos, distanceFunction, K, threshold, outFile );
		}
	}

	static void EncodeNearest(
		PointerList<EncodedFastaSequence> &sequences,
		PointerList<KmerClusterPrototype> &protos,
		DistanceFunction &distanceFunction,
		uint K,
		Distance threshold,
		string &outFile //
	) {
		const uint Q = sequences.Length();
		const uint C = protos.Length();

#define INTERLEAVE 1
#if INTERLEAVE
		ofstream str( outFile );
#else
		vector<BitSet> signatures;

		for ( uint i = 0; i < Q; i++ ) {
			signatures.emplace_back( C );
		}
#endif

		// With no schedule, encodes 500000 sequences against 10000 prototypes in 1751s
		// With schedule(guided), takes 1509s
		// With the output interleaved with calculation, takes:
#pragma omp parallel
		{
#if INTERLEAVE
			BitSet signature( C );
#endif
#pragma omp for schedule(guided)
			for ( uint q = 0; q < Q; q++ ) {
				auto seq = sequences[q];
				uint M = seq->KmerCount( K );
#if INTERLEAVE
				signature.Clear();
#else
				BitSet & signature = signatures[q];
#endif

				for ( uint m = 0; m < M; m++ ) {
					EncodedKmer kmerCode = seq->GetEncodedKmer( m );
					Distance nearestDistance = numeric_limits<Distance>::max();
					uint nearestIndex = 0;

					for ( uint c = 0; c < C; c++ ) {
						EncodedKmer centroidCode = protos[c]->SingletonKmer().PackedEncoding();
						auto dist = distanceFunction( centroidCode, kmerCode, K );

						if ( dist <= threshold && dist < nearestDistance ) {
							nearestIndex = c;
							nearestDistance = dist;
						}
					}

					if ( nearestDistance < numeric_limits<Distance>::max() ) {
						signature.Insert( nearestIndex );
					}
				}

#if INTERLEAVE
#pragma omp critical
				{
					str << sequences[q]->Id() << " " << signature << "\n";
				}
#endif
			}
		}

#if !INTERLEAVE
		ofstream str( outFile );

		for ( uint q = 0; q < Q; q++ ) {
			str << sequences[q]->Id() << " " << signatures[q] << "\n";
		}
#endif
	}


	static void EncodeAny(
		PointerList<EncodedFastaSequence> &sequences,
		PointerList<KmerClusterPrototype> &protos,
		DistanceFunction &distanceFunction,
		uint K,
		Distance threshold,
		string &outFile //
	) {
		const uint Q = sequences.Length();
		const uint C = protos.Length();

#define INTERLEAVE 1
#if INTERLEAVE
		ofstream str( outFile );
#else
		vector<BitSet> signatures;

		for ( uint i = 0; i < Q; i++ ) {
			signatures.emplace_back( C );
		}
#endif

#pragma omp parallel
		{
#if INTERLEAVE
			BitSet signature( C );
#endif
#pragma omp for schedule(guided)
			for ( uint q = 0; q < Q; q++ ) {
				auto seq = sequences[q];
				uint M = seq->KmerCount( K );
#if INTERLEAVE
				signature.Clear();
#else
				BitSet & signature = signatures[q];
#endif

				for ( uint c = 0; c < C; c++ ) {
					EncodedKmer centroidCode = protos[c]->PackedEncoding();

					for ( uint m = 0; m < M; m++ ) {
						EncodedKmer kmerCode = seq->GetEncodedKmer( m );
						auto dist = distanceFunction( centroidCode, kmerCode, K );

						if ( dist <= threshold ) {
							signature.Insert( c );
							break;
						}
					}
				}

#if INTERLEAVE
#pragma omp critical
				{
					str << sequences[q]->Id() << " " << signature << "\n";
				}
#endif
			}
		}

#if !INTERLEAVE
		ofstream str( outFile );

		for ( uint q = 0; q < Q; q++ ) {
			str << sequences[q]->Id() << " " << signatures[q] << "\n";
		}
#endif
	}
};

int main( int argc, char *argv[] ) {
	try {
		Args args( argc, argv );

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = AAClustSig::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << ( end_time - start_time ) << "s" << endl;

		return retCode;
	}
	catch ( Exception ex ) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}
