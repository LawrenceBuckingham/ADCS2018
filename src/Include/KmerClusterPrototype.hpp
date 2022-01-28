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

#include "FastaSequence.hpp"

namespace QutBio {
	/**
	** <summary>
	**	An abuse of the FASTA format which allows a sequence
	**	which consists of a single kmer to be managed, along
	**	with additional properties, notably the size of a cluster.
	**	I might add other things later as needed.
	**
	**	<para>
	**	A prototype represents the centre of a kmer cluster that contains
	**	elements accumulated across sequences in multiple files. An
	**	individual cluster file will index kmers from multiple sequences
	**	in a single file. Prototypes unite clusters across files.
	**	</para>
	**
	** </summary>
	*/
	class KmerClusterPrototype : public EncodedFastaSequence {
	protected:
		size_t size = 0;
		Kmer thisKmer;
		size_t serialNumber = 0;

		const string & idPrefix() {
			static string s = "proto_";
			return s;
		}

		static size_t largestSerialNumber(size_t latest = 0) {
			static size_t value = 0;

#pragma omp critical
			{
				if (latest > value) {
					value = latest;
				}
			}

			return value;
		}

	public:
		typedef KmerClusterPrototype Type;
		typedef Type *Pointer;

		static KmerClusterPrototype * DefaultFactory(
			const string & id,
			const string & classLabel,
			const string & defLine,
			const string & sequence,
			pAlphabet alphabet,
			size_t kmerLength,
			size_t charsPerWord,
			char defaultSymbol
		) {
			return new KmerClusterPrototype(id, classLabel, defLine, sequence, alphabet, kmerLength, charsPerWord, defaultSymbol);
		}

		KmerClusterPrototype(
			const string & id,
			const string & classLabel,
			const string & defLine,
			const string & sequence,
			Alphabet * alphabet,
			size_t wordLength,
			size_t charsPerWord,
			char defaultSymbol = 'x'
		) :
			EncodedFastaSequence(id, classLabel, defLine, sequence, alphabet, wordLength, charsPerWord, defaultSymbol),
			size(0),
			thisKmer( Substring( Sequence().c_str(), 0, Sequence().length() ))
			//
		{
			istringstream idStream(id.substr(idPrefix().length()));
			idStream >> serialNumber;

			largestSerialNumber(serialNumber);

			auto metadata = String::Split(defLine, "|;");
			bool done = false;

			for (auto & meta : metadata) {
				if (done) break;

				if (meta.find("size=") != string::npos) {
					auto properties = String::Split(meta, "=");
					bool sizeFound = false;

					for (auto & p : properties) {
						if (p == "size") {
							sizeFound = true;
						}
						else if (sizeFound) {
							istringstream(p) >> size;
							done = true;
							break;
						}
					}
				}
			}

			thisKmer.Add( this, 0 );
		}

		KmerClusterPrototype(
			size_t serialNumber,
			const string & kmerWord,
			Alphabet * alphabet,
			size_t wordLength,
			size_t charsPerWord,
			char defaultSymbol = 'x'
		) :
			EncodedFastaSequence(GetId(serialNumber), "", "", kmerWord, alphabet, wordLength, charsPerWord, defaultSymbol),
			size(0),
			thisKmer( Substring( Sequence().c_str(), 0, Sequence().length() ) ),
			serialNumber(serialNumber)
			//
		{
			Encode(alphabet, wordLength, charsPerWord, defaultSymbol);
			UpdateDefLine();
			largestSerialNumber(serialNumber);
		}

		KmerClusterPrototype(
			const string & kmerWord,
			Alphabet * alphabet,
			size_t wordLength,
			size_t charsPerWord,
			char defaultSymbol = 'x'
		) :
			EncodedFastaSequence(GetId(largestSerialNumber(largestSerialNumber() + 1)), "", "", kmerWord, alphabet, wordLength, charsPerWord, defaultSymbol),
			size(0),
			thisKmer( Substring( Sequence().c_str(), 0, Sequence().length() ) ),
			serialNumber(largestSerialNumber())
			//
		{
			Encode(alphabet, wordLength, charsPerWord, defaultSymbol);
			UpdateDefLine();
		}

		/// <summary>Get the (total) size of the cluster (s) represented by this prototype.</summary>
		size_t Size() const {
			return size;
		}

		/// <summary>Set the size of the combined collection of indexed kmers.</summary>
		KmerClusterPrototype & Size(size_t size_) {
			this->size = size_;
			UpdateDefLine();
			return *this;
		}

		/// <summary>Make the definition line consistent with the id and size.</summary>
		void UpdateDefLine() {
			ostringstream s;
			s << idPrefix() << serialNumber << "|size=" << size;
			SetDefLine(s.str());
		}

		/**
		**	<summary>
		**	Get the singleton kmer represented by this prototype.
		**	</summary>
		*/

		Kmer & SingletonKmer() {
			return thisKmer;
		}

		EncodedKmer PackedEncoding() {
			return GetEncodedKmer(0);
		}

		/**
		**	<summary>
		**	Get an ID string for use when persisting to a FASTA file.
		**	</summary>
		*/

		string GetId(size_t serialNumber) {
			return idPrefix() + std::to_string(serialNumber);
		}
	};

	typedef KmerClusterPrototype * pKmerClusterPrototype;

}
