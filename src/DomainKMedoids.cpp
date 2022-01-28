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

#include <map>
#include <string>
#include <vector>
#include <omp.h>
#include <unordered_set>

#include <Args.hpp>
#include <Domain.hpp>
#include <FastaSequence.hpp>
#include <KmerCodebook.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <FileUtil.hpp>
#include <KMedoids.hpp>

using namespace QutBio;
using namespace std;

namespace AdHoc {
	struct DomainKMedoids {
		using DistanceFunction = KmerDistanceCache2;
		using KM = KMedoids<DistanceFunction, Kmer>;
		using Cluster = KM::Cluster;
		using SortMode = KM::SortMode;
		using SelectMode = KM::SelectMode;
		using MedoidMode = KM::MedoidMode;

		static void Run( int argc, char** argv ) {
			Args args( argc, argv );
			Params parms( args );
			Alphabet alphabet( parms.matrix );
			BlosumDifferenceFunction rawDist( parms.matrix );
			DistanceFunction distance( &alphabet, &rawDist );

			omp_set_num_threads( parms.numThreads );

			map<string, Domain> domains;
			LoadDomains( parms.domains, domains );

			PointerList<EncodedFastaSequence> db;
			LoadSequences( parms.db, parms.idIndex, parms.classIndex, db, parms.isCaseSensitive );
			EncodedFastaSequence::Index dbIdx( db.Items() );

			vector<const Domain*> domainList;
			for ( auto & p : domains ) {
				if ( parms.wantedDomains.size() == 0 || parms.wantedDomains.find( p.second.pfamId ) != parms.wantedDomains.end() ) {
					domainList.push_back( &( p.second ) );
				}
			}

//			cerr << "domainList.size() = " << domainList.size() << "\n";

			ofstream protoFile( parms.protos );
			ofstream clusterFile( parms.clusters );
			uint clusterCount = 0;

#pragma omp parallel for schedule(dynamic)
			for ( uint i = 0; i < domainList.size(); i++ ) {
				auto dom = domainList[i];
				vector<Subsequence> domainInstances;
				dom->GetInstances( domainInstances, dbIdx );

//#pragma omp critical
//				{
//					cerr << dom->pfamId << ": domainInstances.size() = " << domainInstances.size() << "\n";
//				}

				if ( domainInstances.size() == 0 ) continue;

				vector<Kmer *> clusterProtos;
				vector<Cluster *> clusters;

				KM::Partition( domainInstances, clusterProtos, clusters, parms.kmerLength, parms.threshold, parms.seed, alphabet, distance );

#pragma omp critical
				{
					// cerr << dom->pfamId << ": clusters.size() = " << clusters.size() << "\n";
					for ( uint i = 0; i < clusters.size(); i++ ) {
						Cluster *c = clusters[i];
						Kmer * p = clusterProtos[i];

						// Yuck. KmerClusterPrototype is not thread-safe!!!
						string id = "proto_" + std::to_string( clusterCount++ );
						string defline = id + "|" + dom->pfamId + "|size=" + std::to_string( c->Size() );
						EncodedFastaSequence seq( id, "", defline, p->Word(), &alphabet, parms.kmerLength, distance.CharsPerWord(), alphabet.DefaultSymbol() );
						Kmer newProto( &seq, 0, parms.kmerLength );
						c->prototype = newProto;

						( protoFile << seq ).flush();
						( clusterFile << ( *c ) ).flush();

						delete c;
						delete p;
					}
				}
			}
		}

		static void LoadDomains( const string & domFileName, map<string, Domain> & domains ) {
			OMP_TIMER_DECLARE( domLoad );
			OMP_TIMER_START( domLoad );

			ifstream domFile( domFileName );
			Domain::Load( domFile, domains );
			domFile.close();

			OMP_TIMER_END( domLoad );
			cerr << domains.size() << " domains loaded from " << domFileName << " in " << OMP_TIMER( domLoad ) << "s\n";
		}

		static void LoadSequences(
			const string & fileName,
			int idIndex,
			int classIndex,
			PointerList<EncodedFastaSequence> & seqs,
			bool isCaseSensitive
			//
		) {
			OMP_TIMER_DECLARE( load );
			OMP_TIMER_START( load );
			EncodedFastaSequence::ReadSequences( seqs, fileName, idIndex, classIndex, Alphabet::AA(), 30 );

			if ( !isCaseSensitive ) {
				for ( auto seq : seqs ) {
					String::ToLowerInPlace( seq->Sequence() );
				}
			}
			OMP_TIMER_END( load );
			cerr << seqs.Length() << " sequences loaded from " << fileName << " in " << OMP_TIMER( load ) << "s\n";
		}

		struct Params {
			bool ok;
			string domains, db, protos, clusters;
			int idIndex, classIndex, kmerLength, seed;
			SimilarityMatrix *matrix;
			bool isCaseSensitive;
			Distance threshold;
			int numThreads;
			unordered_set<string> wantedDomains;

			Params( Args & args ) {
				ok = true;

				if ( !args.Get( "domains", domains ) ) {
					cerr << "Argument 'domains' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "db", db ) ) {
					cerr << "Argument 'db' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "protos", protos ) ) {
					cerr << "Argument 'protos' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "clusters", clusters ) ) {
					cerr << "Argument 'clusters' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "kmerLength", kmerLength ) ) {
					cerr << "Argument 'kmerLength' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "idIndex", idIndex ) ) {
					cerr << "Argument 'idIndex' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "classIndex", classIndex ) ) {
					cerr << "Argument 'classIndex' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "isCaseSensitive", isCaseSensitive ) ) {
					cerr << "Argument 'isCaseSensitive' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "threshold", threshold ) ) {
					cerr << "Argument 'threshold' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "seed", seed ) ) {
					cerr << "Argument 'seed' not supplied.\n";
					ok = false;
				}
				if ( !args.Get( "numThreads", numThreads ) ) {
					cerr << "Argument 'numThreads' not supplied.\n";
					ok = false;
				}

				args.Get( "wantedDomains", wantedDomains );

				string matrixError;

				if ( !args.Get( matrix, matrixError ) ) {
					ok = false;
					cerr << matrixError;
				}

				if ( !ok ) {
					cerr << "Example:\nDomainKMedoids.exe --wantedDomains PF00001 PF00002 --domains swissprot.domains --db sp500000.faa --protos sp500000.domain.protos --clusters sp500000.domain.clusters --idIndex 2 --classIndex 3 --matrixId 62 --isCaseSensitive false --kmerLength 30 --threshold 305" << '\n';
					throw Exception( "Invalid arguments.", FileAndLine );
				}
			}

		};

	};
}

using namespace AdHoc;

int main( int argc, char** argv ) {
	try {
		DomainKMedoids::Run( argc, argv );
	}
	catch ( Exception ex ) {
		cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")" << endl;
		return 1;
	}
	catch ( runtime_error & err ) {
		cerr << "Unhandled exception:" << endl << err.what() << endl;
		return 1;
	}
	return 0;
}

mutex QutBio::DistanceType::m;
