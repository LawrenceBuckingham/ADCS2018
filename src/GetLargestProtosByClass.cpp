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

/*
 *	Extracts a subset of a set of clusters based on class labels.
 *	The program reads prototypes and clusters and selects the K
 *	largest clusters (based on prototype size property) from each
 *	class and saves them as a new file.
 *	Prototypes are expected to have id in field 0 and class label
 *	in field 1 of the FASTA defline.
 */

#include <map>
#include <string>
#include <vector>
#include <cstdio>

#include <Args.hpp>
#include <Alphabet.hpp>
#include <Domain.hpp>
#include <FastaSequence.hpp>
#include <KmerCodebook.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <SimilarityMatrix.hpp>
#include <FileUtil.hpp>

using namespace std;
using namespace QutBio;

struct GetLargestProtosByClass {
	using DistanceFunction = KmerDistanceCache2;
	using Cluster = KmerCluster<DistanceFunction, Kmer>;
	using pCluster = Cluster * ;
	using Codebook = KmerCodebook<DistanceFunction, Kmer>;

	static void Run( int argc, char** argv ) {
		Args args( argc, argv );
		Params parms( args );

		//	The alphabet and distance function are dummies required to satisfy
		//	the rather top-heavy API for reading a codebook. That will be fixed 
		//	one day.
		Alphabet *alphabet = Alphabet::AA();
		BlosumDifferenceFunction dist( SimilarityMatrix::Blosum62() );
		DistanceFunction distanceFunction( alphabet, &dist );

		PointerList<EncodedFastaSequence> db;

		try {
			EncodedFastaSequence::ReadSequences( db, parms.db, parms.idIndex, parms.classIndex, alphabet, parms.kmerLength );
		}
		catch ( Exception ex ) {
			cerr << "Error reading from '" << parms.db << "'\n";
			throw;
		}

		EncodedFastaSequence::Index dbIndex( db.Items() );
		KmerIndex kmerIndex( db.AsPointersTo<EncodedFastaSequence>(), db.Length(), parms.kmerLength );

		PointerList<KmerClusterPrototype> protos;

		try {
			EncodedFastaSequence::ReadSequences( protos, parms.protosIn, 0, 1, alphabet, parms.kmerLength );
		}
		catch ( Exception ex ) {
			cerr << "Error reading from '" << parms.protosIn << "'\n";
			throw;
		}

		EncodedFastaSequence::Index protoIndex( protos.AsPointersTo<EncodedFastaSequence>(), protos.Length() );

		Codebook *codebook;
		FILE * f = fopen( parms.clustersIn.c_str(), "r" );

		if ( f ) {
			codebook = new Codebook(
				alphabet, distanceFunction, distanceFunction.CharsPerWord(), parms.kmerLength, dbIndex,
				protoIndex, kmerIndex, f
			);
			fclose( f );
		}
		else {
			cerr << "Error reading codebook from '" << parms.clustersIn << "'\n";
			throw Exception( "Unable to read codebook.", FileAndLine );
		}

		map<string, vector<pKmerClusterPrototype>> seqFamilies;

		for ( auto p : protos ) {
			seqFamilies[p->ClassLabel()].push_back( p );
		}

		map<string, pCluster> clusters;

		for ( auto c : codebook->codebook ) {
			clusters[c->prototype.Sequence()->Id()] = c;
		}

		ofstream pOut( parms.protosOut );
		ofstream cOut( parms.clustersOut );

		for ( auto & pair : seqFamilies ) {
			auto & protosPerFamily = pair.second;

			std::sort( protosPerFamily.begin(), protosPerFamily.end(),
				[]( const pKmerClusterPrototype & lhs, const pKmerClusterPrototype & rhs ) {
				return lhs->Size() > rhs->Size();
			}
			);

			for ( uint i = 0; i < parms.protosPerClass && i < protosPerFamily.size(); i++ ) {
				auto p = protosPerFamily[i];
				pOut << *p;
				cOut << *clusters[p->Id()];
			}
		}
	}

	struct Params {
		bool ok;
		string protosIn, clustersIn, protosOut, clustersOut, db;
		int idIndex, classIndex, protosPerClass, kmerLength;

		Params( Args & args ) {
			ok = true;

			if ( !args.Get( "db", db ) ) {
				cerr << "Argument 'db' not supplied.\n";
				ok = false;
			}
			if ( !args.Get( "protosIn", protosIn ) ) {
				cerr << "Argument 'protosIn' not supplied.\n";
				ok = false;
			}
			if ( !args.Get( "protosOut", protosOut ) ) {
				cerr << "Argument 'protosOut' not supplied.\n";
				ok = false;
			}
			if ( !args.Get( "clustersOut", clustersOut ) ) {
				cerr << "Argument 'clustersOut' not supplied.\n";
				ok = false;
			}
			if ( !args.Get( "clustersIn", clustersIn ) ) {
				cerr << "Argument 'clustersIn' not supplied.\n";
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
			if ( !args.Get( "protosPerClass", protosPerClass ) ) {
				cerr << "Argument 'protosPerClass' not supplied.\n";
				ok = false;
			}
			if ( !args.Get( "kmerLength", kmerLength ) ) {
				cerr << "Argument 'kmerLength' not supplied.\n";
				ok = false;
			}

			if ( !ok ) {
				throw Exception( "Invalid arguments.", FileAndLine );
			}
		}

	};

};

int main( int argc, char** argv ) {
	try {
		GetLargestProtosByClass::Run( argc, argv );
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
