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

#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201103L
#endif

#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <cstdio>

#include "CharMap.hpp"
#include "PrecisionRecallRecord.hpp"
#include "PointerList.hpp"
#include "Selector.hpp"
#include "Util.hpp"
#include "TrecEvalRecord.hpp"
#include "EncodedKmer.hpp"
#include "Alphabet.hpp"

namespace QutBio {
	using EncodingMatrix = vector<vector<KmerWord>>;

	typedef class EncodedFastaSequence *pEncodedFastaSequence;

	class FastaSequence {
	protected:
		string sequence;
		vector<string> metadata;
		size_t idIndex;

	public:
		const string & Id() const {
			if ( idIndex >= metadata.size() ) {
				throw Exception( "Index Out Of Bounds: idIndex", FileAndLine );
			}

			return metadata[idIndex];
		}

		const string & Sequence() const {
			return sequence;
		}

		void SetSequence( const string &value ) {
			sequence.clear();

			for ( auto ch : value ) {
				if ( ch == '-' || isspace( ch ) )
					continue;

				sequence.push_back( ch );
			}
		}

		/**
		**	<summary>
		**	Get the definition line.
		**	</summary>
		*/
		string DefLine() const {
			return String::Join( metadata, "|" );
		}

		/**
		**	<summary>
		**	Set the definition line.
		**	<para>
		**	This does not re-parse the definition line, but it probably should.
		**	</para>
		**	</summary>
		*/

		void SetDefLine( const string &defLine ) {
			metadata = String::Split( defLine, '|' );

			if ( metadata.size() > 0 && metadata[0][0] == '>' ) {
				metadata[0].erase( 0, 1 );
			}
		}



		FastaSequence(
			const string &defLine,
			const string &sequence,
			size_t idIndex
		) : idIndex( idIndex ) {
			SetSequence( sequence );
			SetDefLine( defLine );
		}

		virtual ~FastaSequence() {}

		/// <summary>
		/// Parses a FASTA file and extracts a set of sequences. The fields of the
		/// definition line are assumed to be separated by '|' symbols, with an identifier
		/// in one field and the class label in another.
		/// <para>
		///		The id number may be a GI number, or something else such as a binding site ID or a
		///		application specific id number. This may be any string except one containing a pipe '|'
		///		symbol.
		/// </para>
		/// <para>
		///		The class label is an arbitrary string, which once again, cannot contain embedded pipe
		///		'|' symbols.
		/// </para>
		/// </summary>
		/// <param name="sequenceFile">The name of the file containing the sequences to be loaded.</param>
		/// <param name="idIndex">The index of the pipe-separated field (counting from 0) that contains the id number.</param>
		/// <param name="classIndex">The index of the pipe-separated field (counting from 0) that contains the class label.</param>
		/// <param name="sequences">A list onto which the sequences from the file will be appended.</returns>
		static void ReadSequences(
			istream &reader,
			int idIndex,
			vector<FastaSequence> & db
			//
		) {
			string currentLine;
			string currentDefLine;
			string currentSequence;

			auto update = [&]() {
				if ( currentSequence.size() > 0 ) {
					db.emplace_back( currentDefLine, currentSequence, idIndex );
				}
			};

			auto reset = [&]() {
				currentSequence.clear();
				currentDefLine = currentLine.substr( 1 );
			};

			while ( !( reader.fail() || reader.bad() || reader.eof() ) ) {
				getline( reader, currentLine );
				String::TrimInPlace( currentLine );

				if ( currentLine[0] == '>' ) {
					update();
					reset();
				}
				else {
					currentSequence += String::Trim( currentLine );
				}
			}

			if ( !reader.eof() ) {
				throw Exception( "Error reading from stream.", __FILE__, __LINE__ );
			}

			update();

			// Console.Error( String.Format( "{0} definition lines read.", sequences.Count ) );
		}



		/// <summary>
		/// Parses a FASTA file and extracts a set of sequences. The fields of the
		/// definition line are assumed to be separated by '|' symbols, with an identifier
		/// in one field and the class label in another.
		/// <para>
		///		The id number may be a GI number, or something else such as a binding site ID or a
		///		application specific id number. This may be any string except one containing a pipe '|'
		///		symbol.
		/// </para>
		/// <para>
		///		The class label is an arbitrary string, which once again, cannot contain embedded pipe
		///		'|' symbols.
		/// </para>
		/// </summary>
		/// <param name="sequenceFile">The name of the file containing the sequences to be loaded.</param>
		/// <param name="idIndex">The index of the pipe-separated field (counting from 0) that contains the id number.</param>
		/// <param name="classIndex">The index of the pipe-separated field (counting from 0) that contains the class label.</param>
		/// <returns>A dictionary containing the sequences from the file, indexed by their id numbers.</returns>

		static void ReadSequences(
			const string &fileName,
			int idIndex,
			vector<FastaSequence> & db
		) {
			ifstream reader( fileName );

			if ( reader.fail() ) {
				cerr << "Unable to read from '" << fileName << "'.\n";
			}
			else {
				ReadSequences( reader, idIndex, db );
				reader.close();
			}
		}

		/// <summary> Serialises this FastaSequence to text.
		/// </summary>
		/// <returns></returns>

		string ToString() {
			ostringstream out;
			out << '>' << DefLine() << endl
				<< sequence << endl;
			return out.str();
		}

		/**
		*	<summary>
		*		Gets a normalised histogram with the counts of each symbol.
		*	</summary>
		*	<param name="db">
		*		Reference to a pointer-list of FastaSequence objects in which symbol occurrences are to be calculated.
		*	</param>
		*	<param name="histogram">
		*		Reference to a histogram that will be overwritten with the symbol probabilities observed in the database.
		*	</param>
		*/

		template<typename Collection>
		static Histogram<char> GetSymbolHistogram( const Collection &db ) {
			Histogram<char> histogram;

			for ( auto *seq : db ) {
				histogram.AddRange( seq->Sequence() );
			}

			histogram.Normalise();

			return histogram;
		}

		/**
		*	<summary>
		*		Gets a normalised histogram with the counts of each symbol.
		*	</summary>
		*	<param name="db">
		*		Reference to a pointer-list of FastaSequence objects in which symbol occurrences are to be calculated.
		*	</param>
		*	<param name="histogram">
		*		Reference to a histogram that will be overwritten with the symbol probabilities observed in the database.
		*	</param>
		*/

		/***
		*	<summary>
		*		Uses the supplied selector to choose zero or more kmers from this sequence and pass them back for processing.
		*	</summary>
		*	<param name="kmerLength">The size of the desired kmers.</param>
		*	<param name="selector">A reference to a Selector that determines which, if any, kmers are processed.</param>
		*	<param name="process">A function that is invoked once per selected kmer to process the selection.</param>
		***/

		void SelectKmers( size_t kmerLength, Selector &selector, function<void( FastaSequence &seq, size_t pos, size_t length )> process ) {
			size_t length = sequence.size();

			if ( length < kmerLength )
				return;

			size_t n = length - kmerLength + 1;

			for ( size_t i = 0; i < n; i++ ) {
				if ( selector.SelectThis() ) {
					process( *this, i, kmerLength );
				}
			}
		}

		/***
		*	<summary>
		*		Iterates over the kmers in this sequence and passes them back for processing.
		*	</summary>
		*	<param name="kmerLength">The size of the desired kmers.</param>
		*	<param name="process">A function that is invoked once per selected kmer to process the selection.</param>
		***/

		void SelectKmers( size_t kmerLength, function<void( FastaSequence & seq, size_t pos, size_t length )> process ) {
			size_t n = KmerCount( kmerLength );

			for ( size_t i = 0; i < n; i++ ) {
				process( *this, i, kmerLength );
			}
		}

		class Index : public unordered_map<string, const FastaSequence *> {
		public:
			/**
			*	<summary>
			*		Gets an index which can be used to look up sequences by Id.
			*	</summary>
			*	<param name="db">
			*		Reference to a pointer-list of FastaSequence objects in which symbol occurrences are to be calculated.
			*	</param>
			*	<param name="index">
			*		Reference to a hash table that will be populated with the lookup table. Any previous contents in this
			*		data structure will be erased.
			*	</param>
			*/
			Index( const vector<FastaSequence> & dataset ) {
				for ( auto & seq : dataset ) {
					( *this )[seq.Id()] = &seq;
					// cerr << "Added (" << seq.Id() << "," << (&seq) << ")\n";
				}
			}
		};

		/**
		*	<summary>
		*		Returns the total number of kmers in a dataset.
		*	</summary>
		*	<param name="db">A list of sequences.</param>
		*	<param name="kmerLength">The kmer length.</param>
		*/
		static size_t GetTotalKmerCount( const vector<FastaSequence> & db, size_t kmerLength ) {
			size_t totalKmerCount = 0;

			for ( auto & seq : db ) {
				totalKmerCount += seq.KmerCount( kmerLength );
			}

			return totalKmerCount;
		}

		/**
		*	<summary>
		*	Returns the number of complete kmers of length K in the sequence.
		*	</summary>
		*/
		size_t KmerCount( size_t K ) const {
			auto length = sequence.size();
			return length >= K ? length + 1 - K : 0;
		}

		size_t Length() const {
			return sequence.size();
		}

		friend ostream &operator<<( ostream &out, const FastaSequence &seq ) {
			out << ">" << seq.DefLine() << "\n"
				<< seq.sequence << "\n";
			return out;
		}

		void fprint( FILE * out ) {
			fprintf( out, ">%s\n%s\n", DefLine().c_str(), sequence.c_str() );
		}
	};

	class EncodedFastaSequence {
	protected:
		string id;
		string classLabel;
		string defLine;
		string sequence;

		vector<uint64_t> embedding;
		size_t length;
		pAlphabet alphabet;
		size_t charsPerWord = 0;
		size_t kmerLength;
		char defaultSymbol;

	public:
		// I know, this is messy as anything, but pending a total rewrite, this is what we have for now.
		// Position, rowMinima, colMinima, fragDistAuxData and hit are used by database sequences to store
		// information used in the collection-wide kmer sweep.
		size_t position;
		vector<Distance> rowMinima;
		vector<Distance> colMinima;
		void *fragDistAuxData = 0;
		function<void( void *p )> deleteFragDistAuxData;
		vector<pEncodedFastaSequence> homologs;
		vector<int> classNumbers;

		// When I'm using the cached 2-mer and/or 3-mer score tables, this is the packed image of the kmers in the collection.
		// Update (hacky!) When switching to the centroid-based codebook, I need to access one-char and two-char encodings.
		EncodingMatrix encoding1, encoding2;

		const string &Id() const { return id; }
		const string &ClassLabel() const { return classLabel; }
		const vector<uint64_t> Embedding() { return embedding; }

		string &Sequence() { return sequence; }

		void Sequence( const string &value ) {
			sequence.clear();

			for ( auto ch : value ) {
				if ( ch == '-' || isspace( ch ) )
					continue;

				sequence.push_back( ch );
			}

			length = sequence.size();
		}

		/**
		**	<summary>
		**	Get the definition line.
		**	</summary>
		*/
		const string &DefLine() const { return defLine; }

		/**
		**	<summary>
		**	Set the definition line.
		**	<para>
		**	This does not re-parse the definition line, but it probably should.
		**	</para>
		**	</summary>
		*/

		void SetDefLine( const string &defLine ) {
			// TODO: parse definition line to make sure all properties are consistent.
			this->defLine = defLine;
		}

		static unordered_map<string, int> &ClassNumberRegistry() {
			static unordered_map<string, int> classNumberRegistry;
			return classNumberRegistry;
		}

		static vector<string> &ClassNameRegistry() {
			static vector<string> classNameRegistry;
			return classNameRegistry;
		}

		static int GetClassId( const string &classLabel ) {
			auto &&classNumberRegistry = ClassNumberRegistry();
			int classNumber = -1;

			auto existingRecord = classNumberRegistry.find( classLabel );

			if ( existingRecord == classNumberRegistry.end() ) {
				classNumber = (int) classNumberRegistry.size();
				classNumberRegistry.emplace( classLabel, classNumber );
				auto &&classNameRegistry = ClassNameRegistry();
				classNameRegistry.push_back( classLabel );
			}
			else {
				classNumber = existingRecord->second;
			}

			return classNumber;
		}

		using Factory = function<EncodedFastaSequence *(
			const string &id,
			const string &classLabel,
			const string &defLine,
			const string &sequence,
			pAlphabet alphabet,
			size_t kmerLength,
			size_t charsPerWord,
			char defaultSymbol )>;

		static pEncodedFastaSequence DefaultFactory(
			const string &id,
			const string &classLabel,
			const string &defLine,
			const string &sequence,
			pAlphabet alphabet,
			size_t kmerLength,
			size_t charsPerWord,
			char defaultSymbol ) {
			return new EncodedFastaSequence( id, classLabel, defLine, sequence, alphabet, kmerLength, charsPerWord, defaultSymbol );
		}

		EncodedFastaSequence(
			const string &id,
			const string &classLabel,
			const string &defLine,
			const string &sequence,
			pAlphabet alphabet,
			size_t kmerLength,
			size_t charsPerWord,
			char defaultSymbol ) : id( String::Trim( id ) ),
			classLabel( String::Trim( classLabel ) ),
			defLine( String::Trim( defLine ) )
			//
		{
			Sequence( sequence );

			if ( classLabel.size() > 0 ) {
				auto splitClassLabels = String::Split( classLabel, ';' );

				for ( string &classLabel : splitClassLabels ) {
					int classNumber = GetClassId( classLabel );
					classNumbers.push_back( classNumber );
				}
			}

			Encode( alphabet, kmerLength, charsPerWord, defaultSymbol );
		}

		virtual ~EncodedFastaSequence() {
			if ( fragDistAuxData )
				deleteFragDistAuxData( fragDistAuxData );
		}

		/// <summary>
		/// Parses a FASTA file and extracts a set of sequences. The fields of the
		/// definition line are assumed to be separated by '|' symbols, with an identifier
		/// in one field and the class label in another.
		/// <para>
		///		The id number may be a GI number, or something else such as a binding site ID or a
		///		application specific id number. This may be any string except one containing a pipe '|'
		///		symbol.
		/// </para>
		/// <para>
		///		The class label is an arbitrary string, which once again, cannot contain embedded pipe
		///		'|' symbols.
		/// </para>
		/// </summary>
		/// <param name="sequenceFile">The name of the file containing the sequences to be loaded.</param>
		/// <param name="idIndex">The index of the pipe-separated field (counting from 0) that contains the id number.</param>
		/// <param name="classIndex">The index of the pipe-separated field (counting from 0) that contains the class label.</param>
		/// <param name="sequences">A list onto which the sequences from the file will be appended.</returns>
		static void ReadSequences(
			istream &reader,
			int idIndex,
			int classIndex,
			pAlphabet alphabet,
			size_t kmerLength,
			size_t charsPerWord,
			char defaultSymbol,
			Factory sequenceFactory //
		) {
			string currentLine;
			string currentDefLine;
			string currentSequence;

			auto update = [&]() {
				if ( currentSequence.size() > 0 ) {
					vector<string> parts = String::Split( currentDefLine, "|" );

					if ( idIndex >= (int) parts.size() ) {
						throw Exception( "Index Out Of Bounds: idIndex", FileAndLine );
					}

					string &gi = parts[idIndex];

					if ( classIndex >= 0 && classIndex >= (int) parts.size() ) {
						throw Exception( "Index Out Of Bounds: classIndex", FileAndLine );
					}

					static string empty = "";
					string &classLabel = classIndex >= 0 ? parts[classIndex] : empty;
					sequenceFactory(
						gi,
						classLabel,
						currentDefLine,
						currentSequence,
						alphabet,
						kmerLength,
						charsPerWord,
						defaultSymbol );
				}
			};

			auto reset = [&]() {
				currentSequence = "";
				currentDefLine = currentLine.substr( 1 );
			};

			while ( !( reader.fail() || reader.bad() || reader.eof() ) ) {
				getline( reader, currentLine );
				String::TrimInPlace( currentLine );

				if ( currentLine[0] == '>' ) {
					update();
					reset();
				}
				else {
					currentSequence += String::Trim( currentLine );
				}
			}

			if ( !reader.eof() ) {
				throw Exception( "Error reading from stream.", __FILE__, __LINE__ );
			}

			update();

			// Console.Error( String.Format( "{0} definition lines read.", sequences.Count ) );
		}

		static void ReadSequences(
			PointerList<EncodedFastaSequence> &sequences,
			istream &reader,
			int idIndex,
			int classIndex,
			pAlphabet alphabet,
			size_t kmerLength,
			size_t charsPerWord,
			char defaultSymbol,
			Factory factory ) {
			Factory nested = [&](
				const string &id,
				const string &classLabel,
				const string &defLine,
				const string &sequence,
				pAlphabet alphabet,
				size_t kmerLength,
				size_t charsPerWord,
				char defaultSymbol ) {
				pEncodedFastaSequence seq = factory( id, classLabel, defLine, sequence, alphabet, kmerLength, charsPerWord, defaultSymbol );
				sequences.Add( [seq]() { return seq; } );
				return seq;
			};

			ReadSequences(
				reader,
				idIndex,
				classIndex,
				alphabet,
				kmerLength,
				charsPerWord,
				defaultSymbol,
				nested );
		}

		/// <summary>
		/// Parses a FASTA file and extracts a set of sequences. The fields of the
		/// definition line are assumed to be separated by '|' symbols, with an identifier
		/// in one field and the class label in another.
		/// <para>
		///		The id number may be a GI number, or something else such as a binding site ID or a
		///		application specific id number. This may be any string except one containing a pipe '|'
		///		symbol.
		/// </para>
		/// <para>
		///		The class label is an arbitrary string, which once again, cannot contain embedded pipe
		///		'|' symbols.
		/// </para>
		/// </summary>
		/// <param name="sequenceFile">The name of the file containing the sequences to be loaded.</param>
		/// <param name="idIndex">The index of the pipe-separated field (counting from 0) that contains the id number.</param>
		/// <param name="classIndex">The index of the pipe-separated field (counting from 0) that contains the class label.</param>
		/// <returns>A dictionary containing the sequences from the file, indexed by their id numbers.</returns>

		static void ReadSequences(
			PointerList<EncodedFastaSequence> &sequences,
			const string &fileName,
			int idIndex,
			int classIndex,
			pAlphabet alphabet,
			size_t kmerLength,
			size_t charsPerWord,
			char defaultSymbol,
			Factory factory ) {
			ifstream reader( fileName );

			if ( reader.fail() ) {
				cerr << "Unable to read from '" << fileName << "'.\n";
			}
			else {
				ReadSequences( sequences, reader, idIndex, classIndex, alphabet, kmerLength, charsPerWord, defaultSymbol, factory );
				reader.close();
			}
		}

		/**
		*	<summary>
		*		Pads the sequence out to the designated minimum length
		*	</summary>
		*	<param name="minLength">The required minimum length.</param>
		*	<param name="padding">The symbol used to pad the sequence.</param>
		*/

		void Pad( size_t minLength, char padding ) {
			while ( sequence.size() < minLength ) {
				sequence.push_back( padding );
			}

			length = sequence.size();
		}

		/// <summary> Serialises this FastaSequence to text.
		/// </summary>
		/// <returns></returns>

		string ToString() {
			ostringstream out;
			out << defLine << endl
				<< sequence << endl;
			return out.str();
		}

		/**
		*	<summary>
		*		Returns true iff a designated sequence appears in the
		*		list of homologs of this sequence.
		*		(Needless to say, this is not to be used in search/classification.)
		*	</summary>
		*	<param name="other">A candidate homolog.</param>
		*/
		bool IsHomolog( const EncodedFastaSequence *other ) {
			if ( homologs.size() > 0 ) {
				return find( homologs.begin(), homologs.end(), other ) != homologs.end();
			}
			else {
				for ( auto i : classNumbers ) {
					for ( auto j : other->classNumbers ) {
						if ( i == j )
							return true;
					}
				}
			}

			return false;
		}

		bool IsHomologAny( const vector<const EncodedFastaSequence *> &others ) {
			for ( auto other : others ) {
				if ( IsHomolog( other ) ) {
					return true;
				}
			}

			return false;
		}

		void SetEmbedding( const CharMap &charMap ) {
			embedding.resize( sequence.size() );

			for ( size_t i = 0; i < sequence.size(); i++ ) {
				embedding[i] = charMap.bits[(uint8_t) sequence[i]].lo;
			}
		}

		/***
		*	<summary>
		*		Uses the supplied selector to choose zero or more kmers from this sequence and pass them back for processing.
		*	</summary>
		*	<param name="kmerLength">The size of the desired kmers.</param>
		*	<param name="selector">A reference to a Selector that determines which, if any, kmers are processed.</param>
		*	<param name="process">A function that is invoked once per selected kmer to process the selection.</param>
		***/

		void SelectKmers( size_t kmerLength, Selector &selector, function<void( EncodedFastaSequence *seq, size_t pos, size_t length )> process ) {
			if ( length < kmerLength )
				return;

			size_t n = length - kmerLength + 1;

			for ( size_t i = 0; i < n; i++ ) {
				if ( selector.SelectThis() ) {
					process( this, i, kmerLength );
				}
			}
		}

		/***
		*	<summary>
		*		Iterates over the kmers in this sequence and passes them back for processing.
		*	</summary>
		*	<param name="kmerLength">The size of the desired kmers.</param>
		*	<param name="process">A function that is invoked once per selected kmer to process the selection.</param>
		***/

		void SelectKmers( size_t kmerLength, function<void( EncodedFastaSequence *seq, size_t pos, size_t length )> process ) {
			size_t n = KmerCount( kmerLength );

			for ( size_t i = 0; i < n; i++ ) {
				process( this, i, kmerLength );
			}
		}

		class Index : public unordered_map<string, EncodedFastaSequence *> {
		public:
			/**
			*	<summary>
			*		Gets an index which can be used to look up sequences by Id.
			*	</summary>
			*	<param name="db">
			*		Reference to a pointer-list of FastaSequence objects in which symbol occurrences are to be calculated.
			*	</param>
			*	<param name="index">
			*		Reference to a hash table that will be populated with the lookup table. Any previous contents in this
			*		data structure will be erased.
			*	</param>
			*/
			Index(
				EncodedFastaSequence **dataset,
				size_t length //
			) {
				for ( size_t i = 0; i < length; i++ ) {
					auto *seq = dataset[i];
					( *this )[seq->Id()] = seq;
				}

				if ( length != this->size() ) {
					cerr << "dataset length does not match expected value (this->size)\n"
						<< "length = " << length << "\n"
						<< "this->size() = " << this->size() << "\n";
				}

				assert_equal( length, this->size() );
			}

			Index(
				const vector<EncodedFastaSequence *> dataset					 //
			) : Index( (EncodedFastaSequence **) ( dataset.data() ), dataset.size() ) //
			{}

		}; // namespace QutBio

		   /**
		   *	<summary>
		   *		Returns the total number of kmers in a dataset.
		   *	</summary>
		   *	<param name="db">A list of sequences.</param>
		   *	<param name="kmerLength">The kmer length.</param>
		   */
		static size_t GetTotalKmerCount( PointerList<EncodedFastaSequence> &db, size_t kmerLength ) {
			size_t totalKmerCount = 0;

			for ( auto seq : db ) {
				totalKmerCount += seq->KmerCount( kmerLength );
			}

			return totalKmerCount;
		}

		/**
		*	<summary>
		*		Packs each kmer in the present sequence into a list of KmerWord arrays with the
		*		designated length and density. Each kmer is then represented as one row in the
		*		encoding member of this sequence.
		*
		*		As of January 2017, I need to have simultaneous access to both 1- and 2-letter
		*		numeric codes (1 for the KmerCLusterAL centroid distance calc, and 2 for the
		*		Kmer-vs-Kmer distance calc. So I create two encodings.
		*
		*	</summary>
		*/

		void Encode( Alphabet *alphabet, size_t kmerLength, size_t charsPerWord, char defaultSymbol = 'x' ) {
			this->charsPerWord = charsPerWord;
			this->kmerLength = kmerLength;
			Pad( kmerLength, defaultSymbol );
			alphabet->Encode( sequence.c_str(), length, kmerLength, 1, encoding1 );

			if ( charsPerWord > 1 ) {
				alphabet->Encode( sequence.c_str(), length, kmerLength, charsPerWord, encoding2 );
			}
		}

		/**
		*	<summary>
		*	Returns the number of complete kmers of length K in the sequence.
		*	</summary>
		*/
		size_t KmerCount( size_t K ) {
			return length >= K ? length + 1 - K : 0;
		}

		EncodedKmer GetEncodedKmer( size_t pos ) {
			return charsPerWord == 0 ? GetEncodedKmerError( pos ) : charsPerWord == 1 ? GetEncodedKmer1( pos ) : charsPerWord == 2 ? GetEncodedKmer2( pos ) : charsPerWord == 3 ? GetEncodedKmer3( pos ) : GetEncodedKmerGeneral( pos );
		}

		EncodedKmer GetEncodedKmerGeneral( size_t pos ) {
			return kmerLength <= charsPerWord
				? &encoding2[0][pos]
				: &encoding2[pos % charsPerWord][pos / charsPerWord];
		}

		EncodedKmer GetEncodedKmer1( size_t pos ) {
			return &encoding1[0][pos];
		}

		EncodedKmer GetEncodedKmer2( size_t pos ) {
			return &encoding2[pos % 2][pos / 2];
		}

		EncodedKmer GetEncodedKmer3( size_t pos ) {
			return &encoding2[pos % 3][pos / 3];
		}

		EncodedKmer GetEncodedKmerError( size_t pos ) {
			( cerr << "GetEncodedKmerError -- this function should never be called.\n" ).flush();
			return 0;
		}

		size_t Length() const {
			return length;
		}

		friend ostream &operator<<( ostream &out, const EncodedFastaSequence &seq ) {
			out << ">" << seq.defLine << "\n"
				<< seq.sequence << "\n";
			return out;
		}

		void fprint( FILE * out ) {
			fprintf( out, ">%s\n%s\n", defLine.c_str(), sequence.c_str() );
		}

		template <typename S>
		static void ReadSequences(
			PointerList<S> &sequences,
			const string &dbFileName,
			int idIndex,
			int classIndex,
			Alphabet *alphabet,
			uint wordLength,
			uint encodedCharsPerWord = 2,
			char paddingChar = 'x'
			//
		) {
			EncodedFastaSequence::Factory fact = [&sequences](
				const string &id,
				const string &classLabel,
				const string &defLine,
				const string &sequence,
				pAlphabet alphabet,
				size_t kmerLength,
				size_t charsPerWord,
				char defaultSymbol //
				) {
				// (cerr << "Creating sequence " << id << '\n').flush();
				S *seq = new S( id, classLabel, defLine, sequence, alphabet, kmerLength, charsPerWord, defaultSymbol );
				sequences.Add( [seq]() { return seq; } );
				return seq;
			};

			fstream fastaStream( dbFileName );
			ReadSequences(
				fastaStream,
				idIndex,
				classIndex,
				alphabet,
				wordLength,
				encodedCharsPerWord,
				paddingChar,
				fact );
		}

	};

	struct pFastaHash {
		size_t operator()( const pEncodedFastaSequence &__val ) const noexcept {
			auto &s = __val->Id();
			return _Hash_impl::hash( (void *) s.c_str(), s.length() );
		}
	};

	struct Subsequence {
		EncodedFastaSequence * source;
		size_t start;
		size_t length;
	};

} // namespace QutBio
