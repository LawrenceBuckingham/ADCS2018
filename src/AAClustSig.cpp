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
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <utmpx.h>

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

	struct Signature {
		string id;
		BitSet signature;
		vector<uint> indices;

		Signature( const string &id, uint sigLength ) : id( id ), signature( sigLength ) {}

		~Signature() {}
	};

	static int Run() {
		Params parms;

		if ( !parms.ok ) {
			return 1;
		}

		if ( parms.numThreads > 0 ) {
			omp_set_num_threads( parms.numThreads );
		}

		vector<Signature *> querySigs, dbSigs;

		ReadSignatures( parms.querySigs, querySigs, parms.sigLength );
		ReadSignatures( parms.dbSigs, dbSigs, parms.sigLength );

		vector<vector<uint>> dbIndex( parms.sigLength );
		CreateIndex( dbSigs, dbIndex );

		OMP_TIMER_DECLARE( rank );
		OMP_TIMER_START( rank );
		if ( parms.mode == "merge" ) {
			RankMerge( querySigs, dbSigs, dbIndex, parms.maxResults, parms.outFile );
		}
		else {
			RankBits( querySigs, dbSigs, dbIndex, parms.maxResults, parms.outFile );
		}
		OMP_TIMER_END( rank );
		return 0;
	}

	static void CreateIndex(
		const vector<Signature *> &dbSigs,
		vector<vector<uint>> &index
		//
	) {
		const uint D = dbSigs.size();

		for ( uint d = 0; d < D; d++ ) {
			for ( auto i : dbSigs[d]->indices ) {
				index[i].push_back( d );
			}
		}
	}

	static void RankBits(
		const vector<Signature *> &queries,
		const vector<Signature *> &database,
		const vector<vector<uint>> &dbIndex,
		uint maxResults,
		string &outFile //
	) {
		uint Q = queries.size();

#define INTERLEAVE 1

#if INTERLEAVE
		ofstream out( outFile );
#else
		KnnVector<size_t, double> exemplar( maxResults );
		vector<KnnVector<size_t, double>> allRankings( Q, exemplar );
#endif

#pragma omp parallel
		{

#if INTERLEAVE
			KnnVector<size_t, double> rankings( maxResults );
#endif
			BitSet processed( database.size() );

#pragma omp for
			for ( uint q = 0; q < Q; q++ ) {
#if ! INTERLEAVE
				auto & rankings = allRankings[q];
#endif
				const BitSet &querySignature = queries[q]->signature;

				rankings.clear();
				processed.Clear();

				for ( uint c : queries[q]->indices ) {
					for ( uint d : dbIndex[c] ) {
						if ( !processed.Contains( d ) ) {
							processed.Insert( d );
							const BitSet &dbSignature = database[d]->signature;
							double distance = 1.0 - querySignature.Similarity( dbSignature );

							if ( rankings.canPush( distance ) ) {
								rankings.push( d, distance );
							}
						}
					}
				}

				rankings.sort();

#if INTERLEAVE
#pragma omp critical
				{
					out << queries[q]->id;

					for ( auto & ranking : rankings ) {
						out << " " << database[ranking.second]->id << " " << ( -ranking.first );
					}

					out << " ___eol___ -100000\n";
				}
#endif
			}

		}
#if ! INTERLEAVE
		{
			ofstream out( outFile );

			for ( uint q = 0; q < Q; q++ ) {
				out << queries[q]->id;

				for ( auto & ranking : allRankings[q] ) {
					out << " " << database[ranking.second]->id << " " << ( -ranking.first );
				}

				out << " ___eol___ -100000\n";
			}
		}
#endif
	}

	static void RankMerge(
		const vector<Signature *> &queries,
		const vector<Signature *> &database,
		const vector<vector<uint>> &dbIndex,
		uint maxResults,
		string &outFile //
	) {
		cerr << "RankMerge\n";

		uint Q = queries.size();

#define INTERLEAVE 1

#if INTERLEAVE
		ofstream out( outFile );
#else
		KnnVector<size_t, double> exemplar( maxResults );
		vector<KnnVector<size_t, double>> allRankings( Q, exemplar );
#endif

#pragma omp parallel
		{

#if INTERLEAVE
			KnnVector<size_t, double> rankings( maxResults );
#endif
			BitSet processed( database.size() );

#pragma omp for
			for ( uint q = 0; q < Q; q++ ) {
#if ! INTERLEAVE
				auto & rankings = allRankings[q];
#endif
				auto &queryIndices = queries[q]->indices;

				rankings.clear();
				processed.Clear();

				for ( uint c : queryIndices ) {
					for ( uint d : dbIndex[c] ) {
						if ( !processed.Contains( d ) ) {
							processed.Insert( d );
							auto &dbIndices = database[d]->indices;
							double distance = 1.0 - Jaccard( queryIndices, dbIndices );

							if ( rankings.canPush( distance ) ) {
								rankings.push( d, distance );
							}
						}
					}
				}

				rankings.sort();

#if INTERLEAVE
#pragma omp critical
				{
					out << queries[q]->id;

					for ( auto & ranking : rankings ) {
						out << " " << database[ranking.second]->id << " " << ( -ranking.first );
					}

					out << " ___eol___ -100000\n";
				}
#endif
			}

		}
#if ! INTERLEAVE
		{
			ofstream out( outFile );

			for ( uint q = 0; q < Q; q++ ) {
				out << queries[q]->id;

				for ( auto & ranking : allRankings[q] ) {
					out << " " << database[ranking.second]->id << " " << ( -ranking.first );
				}

				out << " ___eol___ -100000\n";
			}
		}
#endif
	}

	static double Jaccard( const vector<uint> & a, const vector<uint> & b ) {
		uint m = a.size();
		uint n = b.size();
		uint i = 0, j = 0, intersect = 0, union_ = 0;

		while ( i < m && j < n ) {
			uint x = a[i], y = b[j];
			
			union_++;
			
			if ( x < y ) {
				i++;
			}
			else if ( y < x ) {
				j++;
			}
			else {
				intersect ++;
				i++;
				j++;
			}
		}

		union_ += m + n - i - j;

		return (double) intersect / union_;
	}

	static void ReadSignatures(
		string &sigFile,
		vector<Signature *> &signatures,
		uint sigLength //
	) {
		ifstream sigStream( sigFile );

		if ( sigStream.fail() ) {
			cerr << "File " << sigFile << " did not open properly\n";
			throw Exception( "Error reading file " + sigFile, FileAndLine );
		}

		while ( !sigStream.eof() ) {
			string seqId;
			sigStream >> seqId;

			if ( seqId.length() == 0 ) {
				break;
			}

			Signature * sig = new Signature( seqId, sigLength );
			signatures.push_back( sig );
			sigStream >> sig->signature;

			sig->signature.Foreach( [&]( size_t index ) {
				sig->indices.push_back( index );
			} );
		}
	}

	struct Params {
	public:
		string dbSigs;
		string querySigs;
		string outFile;
		size_t numThreads = 8;
		bool ok = true;
		uint maxResults = 1000;
		uint sigLength = 0;
		string mode = "merge";

		Params() {

			if ( arguments->IsDefined( "help" ) ) {
				vector<string> text{
"AAClustSig:  Ranks the top K reference signatures for each sequence in a ",
"             query dataset. This program uses the Jaccard Similarity between ",
"             signatures of reference and query as a proxy for biological ",
"             sequence similarity. Requires precomputed signatures (use ",
"             AAClustSigEncode) but no fasta files, prototypes, or  clusters.",
"",
"--help       Gets this text.",
"",
"--sigLength  Required. The number of bits per signature.",
"",
"--dbSigs     Required. The name of the file which contains signatures for the ",
"             reference sequences.",
"",
"--querySigs  Required. The name of the file which contains signatures for the ",
"             query sequences.",
"",
"--outFile    Required. The name of the output file. This will be a CSV ",
"             document with records containing two fields: the prototype ",
"             sequence ID and information gain.",
"",
"--numThreads Optional; default value = '# cores'. The number of OpenMP ",
"             threads to use in parallel regions.",
"",
"--mode       Optional; default value = 'merge'. The low-level bit counting ",
"             method used to compute cardinality of intersection and union of ",
"             the binary signatures of query and reference sequences. Valid ",
"             values are 'merge', and 'bits'. Merge uses an ordered merge of ",
"             the bit indices (suitable for sparse signatures), while bits ",
"             uses a packed array of boolean together with bitwise operators ",
"             (suitable for dense signatures).",
"",
				};

				for ( auto s : text ) {
					cerr << s << "\n";
				}
			}

			if ( !arguments->Get( "numThreads", numThreads ) ) {
				cerr << arguments->ProgName() << ": note - optional argument '--numThreads' not set"
					"; running with default value "
					<< numThreads << ".\n";
			}

			if ( !arguments->Get( "sigLength", sigLength ) ) {
				cerr << arguments->ProgName() << ": error - required argument '--sigLength' not set.\n";
				ok = false;
			}

			if ( !arguments->Get( "dbSigs", dbSigs ) ) {
				cerr << arguments->ProgName() << ": error - required argument '--dbSigs' not set.\n";
				ok = false;
			}

			if ( !arguments->Get( "querySigs", querySigs ) ) {
				cerr << arguments->ProgName() << ": error - required argument '--querySigs' not set.\n";
				ok = false;
			}

			if ( !arguments->Get( "outFile", outFile ) ) {
				cerr << arguments->ProgName() << ": note - required argument '--outFile' not provided.\n";
				ok = false;
			}

			if ( !arguments->Get( "maxResults", maxResults ) ) {
				cerr << arguments->ProgName() << ": note - optional argument '--maxResults' not set"
					"; running with default value "
					<< maxResults << ".\n";
			}

			if ( !arguments->Get( "mode", mode ) ) {
				cerr << arguments->ProgName() << ": note - optional argument '--mode' not set"
					"; running with default value "
					<< mode << ".\n";
			}

			if ( mode != "bits" && mode != "merge" ) {
				cerr << arguments->ProgName() << ": Mode " << outFile << " is not valid. Use 'merge' or 'bits'.\n";
				ok = false;
			}

			if ( outFile == dbSigs || outFile == querySigs ) {
				cerr << arguments->ProgName() << ": Output file " << outFile << " will overwrite one of your input files.\n";
				ok = false;
			}
		}
	};
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
