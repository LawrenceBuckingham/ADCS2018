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

#include <ios>
#include <unordered_map>
#include <unordered_set>

#include "FastaSequence.hpp"
#include "EncodedKmer.hpp"
#include "Substring.hpp"
#include "DistanceType.hpp"

namespace QutBio {

	using pKmer = class Kmer *;

	class Kmer {
	public:
		/*
		** <summary>The location of an instance of the Kmer.</summary>
		*/
		typedef struct Instance {
			// The Id of the containing sequence.
			pEncodedFastaSequence sequence;

			// The offset of the first character of the kmer from the start of the sequence.
			size_t kmerPosition;

			Instance(pEncodedFastaSequence sequence, size_t kmerPosition) : sequence(sequence), kmerPosition(kmerPosition) {}

			friend ostream & operator<<(ostream & str, const Instance & instance) {
				str << instance.sequence->Id() << ":" << instance.kmerPosition;
				return str;
			}

			EncodedKmer PackedEncoding() {
				return sequence->GetEncodedKmer(kmerPosition);
			}

			EncodedKmer UnpackedEncoding() const {
				return sequence->GetEncodedKmer1(kmerPosition);
			}

			static Instance & Zero() {
				static Instance zero(0, 0);
				return zero;
			}

			friend bool operator==(const Instance & lhs, const Instance & rhs) {
				return lhs.sequence == rhs.sequence && lhs.kmerPosition == rhs.kmerPosition;
			}

			friend bool operator!=(const Instance & lhs, const Instance & rhs) {
				return lhs.sequence != rhs.sequence || lhs.kmerPosition != rhs.kmerPosition;
			}
		} Instance;

	private:
		// The string representation of the kmer.
		Substring substring;

		// The locations of all instances of the current kmer.
		vector<Instance> instances;

		// the address of an encoded form of the kmer.
		EncodedKmer encoding = 0;

		// Certain applications like to know how far a kmer is from a distinguished
		// point such as the centroid of a cluster.
		Distance distanceFromPrototype = numeric_limits<Distance>::max();

		// And certain applications like to know the probability of occurrence

	public:
		/**
		**	Summary:
		**		Construct a kmer belonging to a sequence.
		**
		**	Parameters:
		**		charData: A sub-string containing the character pattern represented by the Kmer;
		*/
		Kmer(const Substring & charData) : substring(charData) {}

		/**
		**	Summary:
		**		Construct a kmer belonging to a sequence.
		**
		**	Parameters:
		**		seq:          a sequence containing a prototypical instance of the kmer;
		**		kmerPosition: the offset of the kmer from the start of the sequence;
		**		kmerLength:   the length of the kmer;
		**		dist:         a distance value associated with the kmer, used when clustering
		**			          Arguably this should not be present in the base class.
		*/
		Kmer(
			pEncodedFastaSequence seq,
			size_t kmerPosition,
			size_t kmerLength,
			Distance dist = numeric_limits<Distance>::max()
		) :
			substring(seq->Sequence().c_str(), kmerPosition, kmerLength),
			distanceFromPrototype(dist)
			//
		{
			Add(seq, kmerPosition);
		}

		virtual ~Kmer() {

		}

		/**
		**	Summary:
		**		Constructs an instance of the base Kmer class.
		**
		**	Parameters:
		**		charData: the defining substring of the kmer.
		*/
		static pKmer DefaultFactory(const Substring & charData) {
			return new Kmer(charData);
		}

		/**
		**	Summary:
		**		Adds a new instance to the current kmer.
		**
		**	Parameters:
		**		seq:          a sequence containing a prototypical instance of the kmer;
		**		kmerPosition: the offset of the kmer from the start of the sequence;
		**		kmerLength:   the length of the kmer;
		**		dist:         a distance value associated with the kmer, used when clustering
		**			          Arguably this should not be present in the base class.
		*/
		void Add(pEncodedFastaSequence seq, size_t kmerPosition, Distance dist = numeric_limits<Distance>::max()) {
			distanceFromPrototype = dist;
			instances.emplace_back(seq, kmerPosition);

			if (instances.size() == 1) {
				encoding = seq->GetEncodedKmer(kmerPosition);
			}
		}

		void Add(const vector<Instance> & other) {
			for (auto & i : other) {
				Add(i.sequence, i.kmerPosition);
			}
		}

		const Substring & Substr() const {
			return substring;
		}

		/// <summary>
		/// Returns a string containing a copy of the kmer.
		/// </summary>

		const string Word() const {
			return string(substring.Chars(), substring.Length());
		}

		const vector<Instance> & Instances() const {
			return instances;
		}

		// Gets the address of the first word in the packed numerically encoded kmer
		// array.
		//	*	This is currently a unit16_t array containing 1, 2, or 3 symbols
		//		packed into a single number.
		//	*	If I switch to a 64-bit dual embedding this will become a uint64_t
		//		array.
		//	*	I may have to keep both formats alive side-by-side for a while.

		EncodedKmer PackedEncoding() const {
			return encoding;
		}

		// Gets the address of the first word in the unpacked numerically encoded kmer
		// array.
		//	*	This is a unit16_t array containing 1 symbol in each element of the
		//		EncodedKmer.

		EncodedKmer UnpackedEncoding() {
			return instances.size() == 0 ? 0 : instances[0].UnpackedEncoding();
		}

		friend ostream & operator << (ostream & str, const Kmer & kmer) {
			for (auto & instance : kmer.instances) {
				str << instance << ";";
			}
			return str;
		}

		Kmer & operator=(const Kmer & other) {
			this->substring = other.substring;
			this->instances = other.instances;
			this->encoding = other.encoding;
			return *this;
		}

		friend bool operator==(const Kmer & lhs, const Kmer & rhs) {
			return lhs.substring == rhs.substring;
		}

		friend bool operator!=(const Kmer & lhs, const Kmer & rhs) {
			return lhs.substring != rhs.substring;
		}

		friend bool operator<(const Kmer & lhs, const Kmer & rhs) {
			return lhs.substring < rhs.substring;
		}

		ostream & write(ostream & stream) const {
			return stream << substring;
		}

		Distance DistanceFromPrototype() const { return distanceFromPrototype; }

		Kmer & DistanceFromPrototype(Distance d) {
			distanceFromPrototype = d;
			return *this;
		}

		/**
		*	<summary>
		*		Returns the number of kmers required to tile the longest sequence in a dataset.
		*	</summary>
		*	<param name="db">A list of sequences.</param>
		*	<param name="kmerLength">The kmer length.</param>
		*/
		static size_t GetMaxKmerCount(PointerList<EncodedFastaSequence> & db, size_t kmerLength) {
			size_t maxKmerCount = 0;

			db.ForEach([&maxKmerCount, kmerLength]( EncodedFastaSequence * seq) {
				size_t K = seq->KmerCount(kmerLength);

				if (K > maxKmerCount) {
					maxKmerCount = K;
				}
			});

			return maxKmerCount;
		}

		pEncodedFastaSequence Sequence() {
			return instances.size() > 0 ? instances[0].sequence : 0;
		}

		size_t KmerPosition() {
			return instances.size() > 0 ? instances[0].kmerPosition : 0;
		}

		Instance & FirstInstance() {
			return instances.size() > 0 ? instances[0] : Instance::Zero();
		}

		size_t Length() const {
			return substring.Length();
		}
	};
}

