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


#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include <set>

#include <Args.hpp>
#include <Domain.hpp>
#include <FastaSequence.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <HBRandom.hpp>

using namespace std;
using namespace QutBio;
// a:\Phd\Notes\SplitFastaHomologs.html

struct SplitFastaHomologs {
	static void Run( int argc, char** argv ) {
		Args args( argc, argv );
		Params parms( args );
		UniformIntRandom<size_t> rand( parms.seed, 1, parms.parts );

		vector<FastaSequence> seqs;
		FastaSequence::ReadSequences( parms.fasta, parms.idIndex, seqs );

		unordered_map<string, size_t> seqIds;

		for ( size_t i = 0; i < seqs.size(); i++ ) {
			seqIds[seqs[i].Id()] = i;
		}

		// parallel array of homolog ids. 
		vector<vector<size_t>> homologs( seqIds.size() );
		ReadHomologs( parms.homologs, seqIds, homologs );

		for ( auto & vec : homologs ) {
			sort( vec.begin(), vec.end() );
			removeDuplicates( vec );
		}

		vector<size_t> partNumbers( seqIds.size() );

		for ( size_t seqId = 0; seqId < seqIds.size(); seqId++ ) {
			partNumbers[seqId] = rand();
		}

		//	We now have allocated each sequence to a random partition.
		//	For each part, must write the sequences to corresponding 
		//	FASTA file, and in the corresponding homologs file, write
		//	a row for each sequence in the part, but which only includes
		//	homologs ids which belong to OTHER parts.
		//	And the OTHER sequences, will be written to a per-partition 
		//	database file.

		for ( size_t part = 1; part <= parms.parts; part++ ) {
			char partBuff[128];
			sprintf( partBuff, "%02zu", part );
			ofstream testFile( parms.outStub + "." + partBuff + ".test.faa" );
			ofstream trainFile( parms.outStub + "." + partBuff + ".train.faa" );
			ofstream homologFile( parms.outStub + "." + partBuff + ".homologs" );

			for ( size_t seqId = 0; seqId < seqs.size(); seqId++ ) {
				if ( partNumbers[seqId] == part ) {
					testFile << seqs[seqId];
					homologFile << seqs[seqId].Id();

					for ( size_t homologId : homologs[seqId] ) {
						if ( partNumbers[homologId] != part ) {
							homologFile << " " << seqs[homologId].Id();
						}
					}

					homologFile << "\n";
				}
				else {
					trainFile << seqs[seqId];
				}
			}
		}
	}

	// https://www.geeksforgeeks.org/remove-duplicates-sorted-array/
	static void removeDuplicates( vector<size_t> & arr ) {
		const size_t n = arr.size();

		if ( n == 0 || n == 1 )
			return;

		int j = 0;

		for ( int i = 0; i < n - 1; i++ ) {
			if ( arr[i] != arr[i + 1] ) {
				arr[j++] = arr[i];
			}
		}

		arr[j++] = arr[n - 1];

		arr.resize(j);
	}


	static void ReadHomologs(
		const string & fileName,
		const unordered_map<string, size_t> & seqIds,
		vector<vector<size_t>> & homologs
		//
	) {
		FILE *fp;
		fp = fopen( fileName.c_str(), "r" );

		if ( !fp ) {
			fprintf( stderr, "homologs file %s does not exist.\n", fileName.c_str() );
			exit( 1 );
		}

		const size_t MISSING = string::npos;

		//size_t n = 0;

		for ( ;; ) {
			char topic_[256];
			char delimiter;

			int items = fscanf( fp, "%256s%c", topic_, &delimiter );

			if ( items < 2 ) {
				break;
			}

			string topic( topic_ );
			size_t topicId = GetTopicId( topic, seqIds );

			//if ( n++ ) {
			//	( cerr << "\r" << n << ": Current topic: " << topic << "                    " ).flush();
			//}

			while ( delimiter == ' ' ) {
				char doc_[256];

				items = fscanf( fp, "%256s%c", doc_, &delimiter );

				if ( items < 2 ) break;

				string doc( doc_ );
				size_t docId = GetTopicId( doc, seqIds );

				if ( topicId != MISSING && docId != MISSING ) {
					homologs[topicId].push_back( docId );
				}
			}
		}

		fclose( fp );
	}

	static size_t GetTopicId(
		const string & topic,
		const unordered_map<string, size_t> & topicIds
	) {
		auto topicPos = topicIds.find( topic );
		size_t topicId;

		if ( topicPos == topicIds.end() ) {
			topicId = string::npos;
		}
		else {
			topicId = topicPos->second;
		}

		return topicId;
	}

	struct Params {
		bool ok;
		string fasta, homologs, outStub;
		int idIndex, seed, parts;

		Params( Args & args ) {
			ok = true;

			if ( !args.Get( "fasta", fasta ) ) {
				cerr << "Argument 'fasta' not supplied.\n";
				ok = false;
			}

			if ( !args.Get( "homologs", homologs ) ) {
				cerr << "Argument 'homologs' not supplied.\n";
				ok = false;
			}

			if ( !args.Get( "outStub", outStub ) ) {
				cerr << "Argument 'outStub' not supplied.\n";
				ok = false;
			}

			if ( !args.Get( "idIndex", idIndex ) ) {
				cerr << "Argument 'idIndex' not supplied.\n";
				ok = false;
			}

			if ( !args.Get( "seed", seed ) ) {
				cerr << "Argument 'seed' not supplied.\n";
				ok = false;
			}

			if ( !args.Get( "parts", parts ) ) {
				cerr << "Argument 'parts' not supplied.\n";
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
		SplitFastaHomologs::Run( argc, argv );
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

