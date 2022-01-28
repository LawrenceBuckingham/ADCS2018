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
#include <ctime>
#include <cfloat>
#include <functional>

#include "Alphabet.hpp"
#include "Console.hpp"
#include "Delegates.hpp"
#include "EncodedKmer.hpp"
#include "Util.hpp"
#include "Array.hpp"


namespace QutBio {

	/**
	 *	<summary>
	 *	Base class for delegate type which knows how to calculate a raw kmer
	 *	distance.
	 *	</summary>
	 */
	class RawKmerDistanceFunction {
	public:
		/**
		 *	<summary>
		 *	Override this function to implement the low-level kmer distance calc
		 *	between two strings.
		 *	</summary>
		 *	<param name="x">One kmer.</param>
		 *	<param name="y">Another kmer.</param>
		 *	<param name="k">The k in kmer -- length of the sequences to be compared.</param>
		 */
		virtual Distance operator()(const char * x, const char * y, uint k) = 0;
	};

	/**
	 *	<summary>
	 *	Implement the BLOSUM distance used by Halperin et al.:
	 *		d(x,y) = b(x,x) + b(y,y) - 2 b(x,y).
	 *	</summary>
	 */
	class HalperinBlosumDistanceFunction : public RawKmerDistanceFunction {
	private:
		SimilarityMatrix * matrix;
	public:
		/**
		 *	<summary>
		 *	initialise the calculator, saving the similarity matrix for later use.
		 *	</summary>
		 *	<param name="matrix">
		 *		A non-null pointer to a similarity matrix that will provide the values.
		 *	</param>
		 */
		HalperinBlosumDistanceFunction(SimilarityMatrix * matrix) : matrix(matrix) {}

		/**
		 *	<summary>
		 *	Override which implements the BLOSUM distance: d(x,y) = max { b(a,a) | a \in alphabet } - b(x,y).
		 *	</summary>
		 *	<param name="x">One kmer.</param>
		 *	<param name="y">Another kmer.</param>
		 *	<param name="k">The k in kmer -- length of the sequences to be compared.</param>
		 */
		Distance operator () (const char * x, const char * y, uint length) {
			Distance d = matrix->HalperinDistance(x, y, length);
			return d;
		}
	};

	/**
	 *	<summary>
	 *	Implement the BLOSUM distance: d(x,y) = max { b(a,a) | a \in alphabet } - b(x,y).
	 *	</summary>
	 */
	class BlosumDifferenceFunction : public RawKmerDistanceFunction {
	private:
		SimilarityMatrix * matrix;
	public:
		/**
		 *	<summary>
		 *	Override this function to implement the low-level kmer distance calc
		 *	between two strings.
		 *	</summary>
		 *	<param name="matrix">A non-null pointer to a similarity matrix that will provide the values.</param>
		 */
		BlosumDifferenceFunction(SimilarityMatrix * matrix) : matrix(matrix) {}

		/**
		 *	<summary>
		 *	Override which implements the BLOSUM distance: d(x,y) = max { b(a,a) | a \in alphabet } - b(x,y).
		 *	</summary>
		 *	<param name="x">One kmer.</param>
		 *	<param name="y">Another kmer.</param>
		 *	<param name="k">The k in kmer -- length of the sequences to be compared.</param>
		 */
		Distance operator () (const char * x, const char * y, uint length) {
			Distance d = length * matrix->maxValue - matrix->Similarity(x, y, length);
			return d;
		}
	};

	/**
	 *	<summary>
	 *	UngappedEditDistanceCalculator.
	 *	</summary>
	 */
	class UngappedEditDistanceFunction : public RawKmerDistanceFunction {
	public:

		/**
		 *	<summary>
		 *	Override which implements the BLOSUM distance: d(x,y) = max { b(a,a) | a \in alphabet } - b(x,y).
		 *	</summary>
		 *	<param name="x">One kmer.</param>
		 *	<param name="y">Another kmer.</param>
		 *	<param name="k">The k in kmer -- length of the sequences to be compared.</param>
		 */
		Distance operator () (const char * x, const char * y, uint length) {
			int diffs = 0;

			for ( uint i = 0; i < length; i++ ) {
				if ( x[i] != y[i] ) diffs++;
			}

			return diffs;
		}
	};

	class RawKmerDistanceFunctionFactory {
	public:
		static RawKmerDistanceFunction * Factory(DistanceType * dist, int matrixId) {
			if ( dist == DistanceType::UngappedEdit() ) {
				return new UngappedEditDistanceFunction();
			}

			SimilarityMatrix * matrix = SimilarityMatrix::GetBlosum(matrixId);

			if ( dist == DistanceType::HalperinEtAl() )
				return new HalperinBlosumDistanceFunction(matrix);
			else
				return new BlosumDifferenceFunction(matrix);
		}
	};

	class KmerDistanceCache {
	protected:
		Alphabet * alphabet;
		RawKmerDistanceFunction & dist;

		typedef int8_t CacheType;
		typedef CacheType * pCacheType;
		typedef pCacheType * ppCacheType;

	public:
		KmerDistanceCache(Alphabet * alphabet, RawKmerDistanceFunction * dist) 
			: alphabet(alphabet), dist(*dist)
		{}

		virtual ~KmerDistanceCache() {}

	protected:
		void PrecomputeDistances(uint charsPerWord, pCacheType &kmerDistanceTable, uint & vocabSize_) {
			int len = alphabet->Size();

			uint vocabSize = len;

			for ( uint i = 1; i < charsPerWord; i++ ) { vocabSize *= len; }

			vocabSize_ = vocabSize;

			assert_true(vocabSize <= numeric_limits<KmerWord>::max());

			vector<char *> vocab(vocabSize);

			KmerWord code;

			for ( uint i = 0; i < vocabSize; i++ ) {
				code = i;
				char * x = new char[charsPerWord + 1];
				x[charsPerWord] = 0;
				alphabet->Decode(&code, charsPerWord, charsPerWord, x);
				vocab[i] = x;
			}

			kmerDistanceTable = new CacheType[vocabSize*vocabSize];

			for ( uint i = 0; i < vocabSize; i++ ) {
				char * x = vocab[i];

				// Console.Log( String.Format( "Adding entries for {0} - {1} out of {2}", x, i * vocabSize, vocabSize * vocabSize ) );

				for ( uint j = 0; j <= i; j++ ) {
					char * y = vocab[j];
					Distance t = dist(x, y, charsPerWord);
					kmerDistanceTable[i * vocabSize + j] = (CacheType) t;
					kmerDistanceTable[j * vocabSize + i] = (CacheType) t;
				}
			}

			for ( uint i = 0; i < vocabSize; i++ ) {
				delete[] vocab[i];
			}
		}

	};

	/**
	*	<summary>
	*		Pre-computed kmer distance tables for k == 1, 2, 3, and
	*		kmer distance function which uses them instead of looking
	*		up the matrix to save a lookups and/or loops.
	*	</summary>
	*/
	class KmerDistanceCache3 : public KmerDistanceCache {
	public:
		typedef Distance(*KmerDistanceFunction)(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength);

	protected:
		pCacheType kmerDistances1;
		uint vocabSize1;

		pCacheType kmerDistances2;
		uint vocabSize2;

		pCacheType kmerDistances3;
		uint vocabSize3;

		KmerDistanceFunction kmerDistanceFunctions[18];

	public:

		KmerDistanceCache3(Alphabet * alphabet, RawKmerDistanceFunction * dist) : KmerDistanceCache(alphabet, dist) {
			PrecomputeDistances();
			int i = 0;
			kmerDistanceFunctions[i++] = GetKmerDistances_1;
			kmerDistanceFunctions[i++] = GetKmerDistances_2;
			kmerDistanceFunctions[i++] = GetKmerDistances_3;
			kmerDistanceFunctions[i++] = GetKmerDistances_4;
			kmerDistanceFunctions[i++] = GetKmerDistances_5;
			kmerDistanceFunctions[i++] = GetKmerDistances_6;
			kmerDistanceFunctions[i++] = GetKmerDistances_7;
			kmerDistanceFunctions[i++] = GetKmerDistances_8;
			kmerDistanceFunctions[i++] = GetKmerDistances_9;
			kmerDistanceFunctions[i++] = GetKmerDistances_10;
			kmerDistanceFunctions[i++] = GetKmerDistances_11;
			kmerDistanceFunctions[i++] = GetKmerDistances_12;
			kmerDistanceFunctions[i++] = GetKmerDistances_13;
			kmerDistanceFunctions[i++] = GetKmerDistances_14;
			kmerDistanceFunctions[i++] = GetKmerDistances_15;
			kmerDistanceFunctions[i++] = GetKmerDistances_16;
			kmerDistanceFunctions[i++] = GetKmerDistances_17;
			kmerDistanceFunctions[i++] = GetKmerDistances_18;
		}

		KmerDistanceCache3( const KmerDistanceCache3 & other ) = delete;

		KmerDistanceCache3 & operator=( const KmerDistanceCache3 & other ) = delete;

		virtual ~KmerDistanceCache3() {
			delete[] kmerDistances1;
			delete[] kmerDistances2;
			delete[] kmerDistances3;
		}

		const KmerDistanceFunction * KmerDistanceFunctions(void) {
			return kmerDistanceFunctions;
		}

		size_t CharsPerWord() {
			return 3;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 3 or less.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_Any(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			Distance dist = 0;
			KmerWord * s = sKmerCode;
			KmerWord * t = tKmerCode;
			uint remaining = kmerLength;

			while ( remaining > 18 ) {
				dist += GetKmerDistances_18(instance, s, t, 18);
				s += 6;
				t += 6;
				remaining -= 18;
			}

			auto f = instance.kmerDistanceFunctions[remaining - 1];
			dist += f(instance, s, t, remaining);

			return dist;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 3 or less.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_3(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			size_t s1 = sKmerCode[0];
			size_t t1 = tKmerCode[0];

			Distance result = instance.kmerDistances3[s1*instance.vocabSize3 + t1];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 3 or less.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_1(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord t1 = tKmerCode[0];

			Distance result = instance.kmerDistances1[s1 * instance.vocabSize1 + t1];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 3 or less.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_2(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord t1 = tKmerCode[0];

			Distance result = instance.kmerDistances2[s1 * instance.vocabSize2 + t1];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 4 to 6 inclusive.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_4(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord s2 = sKmerCode[1];

			KmerWord t1 = tKmerCode[0];
			KmerWord t2 = tKmerCode[1];

			Distance result = instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances1[s2 * instance.vocabSize1 + t2];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 4 to 6 inclusive.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_5(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord s2 = sKmerCode[1];

			KmerWord t1 = tKmerCode[0];
			KmerWord t2 = tKmerCode[1];

			Distance result = instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances2[s2 * instance.vocabSize2 + t2];
			return result;
		}


		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 4 to 6 inclusive.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_6(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord s2 = sKmerCode[1];

			KmerWord t1 = tKmerCode[0];
			KmerWord t2 = tKmerCode[1];

			Distance result = instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 7 to 9 inclusive.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_7(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord s2 = sKmerCode[1];
			KmerWord s3 = sKmerCode[2];

			KmerWord t1 = tKmerCode[0];
			KmerWord t2 = tKmerCode[1];
			KmerWord t3 = tKmerCode[2];

			Distance result = instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances1[s3 * instance.vocabSize1 + t3];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 7 to 9 inclusive.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_8(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord s2 = sKmerCode[1];
			KmerWord s3 = sKmerCode[2];

			KmerWord t1 = tKmerCode[0];
			KmerWord t2 = tKmerCode[1];
			KmerWord t3 = tKmerCode[2];

			Distance result = instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances2[s3 * instance.vocabSize2 + t3];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 7 to 9 inclusive.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_9(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord s2 = sKmerCode[1];
			KmerWord s3 = sKmerCode[2];

			KmerWord t1 = tKmerCode[0];
			KmerWord t2 = tKmerCode[1];
			KmerWord t3 = tKmerCode[2];

			Distance result = instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances3[s3 * instance.vocabSize3 + t3];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 10.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_10(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord s2 = sKmerCode[1];
			KmerWord s3 = sKmerCode[2];
			KmerWord s4 = sKmerCode[3];

			KmerWord t1 = tKmerCode[0];
			KmerWord t2 = tKmerCode[1];
			KmerWord t3 = tKmerCode[2];
			KmerWord t4 = tKmerCode[3];

			Distance result = instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances3[s3 * instance.vocabSize3 + t3];
			/*   */ result += instance.kmerDistances1[s4 * instance.vocabSize1 + t4];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 10.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_11(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord s2 = sKmerCode[1];
			KmerWord s3 = sKmerCode[2];
			KmerWord s4 = sKmerCode[3];

			KmerWord t1 = tKmerCode[0];
			KmerWord t2 = tKmerCode[1];
			KmerWord t3 = tKmerCode[2];
			KmerWord t4 = tKmerCode[3];

			Distance result = instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances3[s3 * instance.vocabSize3 + t3];
			/*   */ result += instance.kmerDistances2[s4 * instance.vocabSize2 + t4];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 12.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_12(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord s2 = sKmerCode[1];
			KmerWord s3 = sKmerCode[2];
			KmerWord s4 = sKmerCode[3];

			KmerWord t1 = tKmerCode[0];
			KmerWord t2 = tKmerCode[1];
			KmerWord t3 = tKmerCode[2];
			KmerWord t4 = tKmerCode[3];

			Distance result = instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances3[s3 * instance.vocabSize3 + t3];
			/*   */ result += instance.kmerDistances3[s4 * instance.vocabSize3 + t4];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 13.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_13(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s1 = sKmerCode[0];
			KmerWord s2 = sKmerCode[1];
			KmerWord s3 = sKmerCode[2];
			KmerWord s4 = sKmerCode[3];
			KmerWord s5 = sKmerCode[4];

			KmerWord t1 = tKmerCode[0];
			KmerWord t2 = tKmerCode[1];
			KmerWord t3 = tKmerCode[2];
			KmerWord t4 = tKmerCode[3];
			KmerWord t5 = tKmerCode[4];

			Distance result = instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances3[s3 * instance.vocabSize3 + t3];
			/*   */ result += instance.kmerDistances3[s4 * instance.vocabSize3 + t4];
			/*   */ result += instance.kmerDistances1[s5 * instance.vocabSize1 + t5];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 14.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_14(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s0 = sKmerCode[0];
			KmerWord s1 = sKmerCode[1];
			KmerWord s2 = sKmerCode[2];
			KmerWord s3 = sKmerCode[3];
			KmerWord s4 = sKmerCode[4];

			KmerWord t0 = tKmerCode[0];
			KmerWord t1 = tKmerCode[1];
			KmerWord t2 = tKmerCode[2];
			KmerWord t3 = tKmerCode[3];
			KmerWord t4 = tKmerCode[4];

			Distance result = instance.kmerDistances3[s0 * instance.vocabSize3 + t0];
			/*   */ result += instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances3[s3 * instance.vocabSize3 + t3];
			/*   */ result += instance.kmerDistances2[s4 * instance.vocabSize2 + t4];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 15.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_15(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s0 = sKmerCode[0];
			KmerWord s1 = sKmerCode[1];
			KmerWord s2 = sKmerCode[2];
			KmerWord s3 = sKmerCode[3];
			KmerWord s4 = sKmerCode[4];

			KmerWord t0 = tKmerCode[0];
			KmerWord t1 = tKmerCode[1];
			KmerWord t2 = tKmerCode[2];
			KmerWord t3 = tKmerCode[3];
			KmerWord t4 = tKmerCode[4];

			Distance result = instance.kmerDistances3[s0 * instance.vocabSize3 + t0];
			/*   */ result += instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances3[s3 * instance.vocabSize3 + t3];
			/*   */ result += instance.kmerDistances3[s4 * instance.vocabSize3 + t4];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 16.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_16(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s0 = sKmerCode[0];
			KmerWord s1 = sKmerCode[1];
			KmerWord s2 = sKmerCode[2];
			KmerWord s3 = sKmerCode[3];
			KmerWord s4 = sKmerCode[4];
			KmerWord s5 = sKmerCode[5];

			KmerWord t0 = tKmerCode[0];
			KmerWord t1 = tKmerCode[1];
			KmerWord t2 = tKmerCode[2];
			KmerWord t3 = tKmerCode[3];
			KmerWord t4 = tKmerCode[4];
			KmerWord t5 = tKmerCode[5];

			Distance result = instance.kmerDistances3[s0 * instance.vocabSize3 + t0];
			/*   */ result += instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances3[s3 * instance.vocabSize3 + t3];
			/*   */ result += instance.kmerDistances3[s4 * instance.vocabSize3 + t4];
			/*   */ result += instance.kmerDistances1[s5 * instance.vocabSize1 + t5];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 16.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_17(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s0 = sKmerCode[0];
			KmerWord s1 = sKmerCode[1];
			KmerWord s2 = sKmerCode[2];
			KmerWord s3 = sKmerCode[3];
			KmerWord s4 = sKmerCode[4];
			KmerWord s5 = sKmerCode[5];

			KmerWord t0 = tKmerCode[0];
			KmerWord t1 = tKmerCode[1];
			KmerWord t2 = tKmerCode[2];
			KmerWord t3 = tKmerCode[3];
			KmerWord t4 = tKmerCode[4];
			KmerWord t5 = tKmerCode[5];

			Distance result = instance.kmerDistances3[s0 * instance.vocabSize3 + t0];
			/*   */ result += instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances3[s3 * instance.vocabSize3 + t3];
			/*   */ result += instance.kmerDistances3[s4 * instance.vocabSize3 + t4];
			/*   */ result += instance.kmerDistances2[s5 * instance.vocabSize2 + t5];
			return result;
		}

		/// <summary> Computes distance between two kmers, s and t, 
		///		assuming that they are length 16.
		/// </summary>
		/// <param name="sCode">The numeric code of kmer s.</param>
		/// <param name="tCode">The numeric code of kmer t.</param>
		/// <returns></returns>

		static Distance GetKmerDistances_18(KmerDistanceCache3 & instance, KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			KmerWord s0 = sKmerCode[0];
			KmerWord s1 = sKmerCode[1];
			KmerWord s2 = sKmerCode[2];
			KmerWord s3 = sKmerCode[3];
			KmerWord s4 = sKmerCode[4];
			KmerWord s5 = sKmerCode[5];

			KmerWord t0 = tKmerCode[0];
			KmerWord t1 = tKmerCode[1];
			KmerWord t2 = tKmerCode[2];
			KmerWord t3 = tKmerCode[3];
			KmerWord t4 = tKmerCode[4];
			KmerWord t5 = tKmerCode[5];

			Distance result = instance.kmerDistances3[s0 * instance.vocabSize3 + t0];
			/*   */ result += instance.kmerDistances3[s1 * instance.vocabSize3 + t1];
			/*   */ result += instance.kmerDistances3[s2 * instance.vocabSize3 + t2];
			/*   */ result += instance.kmerDistances3[s3 * instance.vocabSize3 + t3];
			/*   */ result += instance.kmerDistances3[s4 * instance.vocabSize3 + t4];
			/*   */ result += instance.kmerDistances3[s5 * instance.vocabSize3 + t5];
			return result;
		}

		Distance GetDistance(KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			uint numThrees = kmerLength / 3;
			uint rem = kmerLength % 3;
			Distance dist = 0;

			uint i = 0;

			for ( ; i < numThrees; i++ ) {
				dist += kmerDistances3[sKmerCode[i] * vocabSize3 + tKmerCode[i]];
			}

			if ( rem > 0 ) {
				dist += rem == 1
					? kmerDistances1[sKmerCode[i] * vocabSize1 + tKmerCode[i]]
					: kmerDistances2[sKmerCode[i] * vocabSize2 + tKmerCode[i]];
			}

			return dist;
		}

		/**
		*	<summary>
		*		Computes the distance between two kmers, on the assumption that all terms of the
		*		sum will be strictly non-negative. If the distance is greater than the designated
		*		threshold, function returns false. Otherwise, if the distance is equal to or less
		*		than the threshold, returns true and as a side effect, stores the distance in the
		*		variable referenced by result.
		*	</summary>
		*/

		bool IsWithin(KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength, Distance threshold, Distance & result) {
			uint numThrees = kmerLength / 3;
			uint rem = kmerLength % 3;
			Distance dist = 0;

			uint i = 0;

			for ( ; i < numThrees; i++ ) {
				dist += kmerDistances3[sKmerCode[i] * vocabSize3 + tKmerCode[i]];

				if ( dist > threshold ) {
					return false;
				}
			}

			if ( rem > 0 ) {
				dist += rem == 1
					? kmerDistances1[sKmerCode[i] * vocabSize1 + tKmerCode[i]]
					: kmerDistances2[sKmerCode[i] * vocabSize2 + tKmerCode[i]];

				if ( dist > threshold ) {
					return false;
				}
			}

			result = dist;
			return true;
		}

	private:
		void PrecomputeDistances() {
			KmerDistanceCache::PrecomputeDistances(1, kmerDistances1, vocabSize1);
			KmerDistanceCache::PrecomputeDistances(2, kmerDistances2, vocabSize2);
			KmerDistanceCache::PrecomputeDistances(3, kmerDistances3, vocabSize3);
		}
	};

	/**
	*	<summary>
	*		Pre-computed kmer distance tables for k == 1, 2, and
	*		kmer distance function which uses them instead of looking
	*		up the matrix to save a lookups and/or loops.
	*	</summary>
	*/
	class KmerDistanceCache2: KmerDistanceCache {
		// TODO: Define and implement IKmerWordDistance, offering the GetDistance function.
	protected:
		CacheType * kmerDistances1;
		uint vocabSize1;

		CacheType * kmerDistances2;
		uint vocabSize2;

	public:
		KmerDistanceCache2(Alphabet * alphabet, RawKmerDistanceFunction * dist) : KmerDistanceCache(alphabet, dist) {
			PrecomputeDistances();
		}

		KmerDistanceCache2( const KmerDistanceCache2 & other ) = delete;

		KmerDistanceCache2 & operator=( const KmerDistanceCache2 & other ) = delete;


		virtual ~KmerDistanceCache2() {
			delete[] kmerDistances1;
			delete[] kmerDistances2;
		}

		size_t CharsPerWord() const {
			return 2;
		}

		//Distance GetDistance( KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) const {
		//	uint numTwos = kmerLength >> 1;
		//	uint rem = kmerLength & 1;
		//	Distance dist = 0;

		//	uint i;

		//	for ( i = 0; i < numTwos; i++ ) {
		//		dist += kmerDistances2[(*sKmerCode++) * vocabSize2 + (*tKmerCode++)];
		//	}

		//	if ( rem > 0 ) {
		//		dist += kmerDistances1[(*sKmerCode++) * vocabSize1 + (*tKmerCode++)];
		//	}

		//	return dist;
		//}

		Distance operator()( KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) const {
			uint numTwos = kmerLength >> 1;
			uint rem = kmerLength & 1;
			Distance dist = 0;

			uint i;

			for ( i = 0; i < numTwos; i++ ) {
				dist += kmerDistances2[(*sKmerCode++) * vocabSize2 + (*tKmerCode++)];
			}

			if ( rem > 0 ) {
				dist += kmerDistances1[(*sKmerCode++) * vocabSize1 + (*tKmerCode++)];
			}

			return dist;
		}

		Distance GetDistance1(KmerWord x, KmerWord y) const {
			return kmerDistances1[x * vocabSize1 + y];
		}

		Distance GetDistance2(KmerWord x, KmerWord y) const {
			return kmerDistances2[x * vocabSize2 + y];
		}

		/**
		*	<summary>
		*		Computes the distance between two kmers, on the assumption that all terms of the
		*		sum will be strictly non-negative. If the distance is greater than the designated
		*		threshold, function returns false. Otherwise, if the distance is equal to or less
		*		than the threshold, returns true and as a side effect, stores the distance in the
		*		variable referenced by result.
		*	</summary>
		*/

		bool IsWithin(const KmerWord * sKmerCode, const KmerWord * tKmerCode, uint kmerLength, Distance threshold, Distance & result) {
			uint numTwos = kmerLength / 2;
			uint rem = kmerLength % 2;
			Distance dist = 0;

			uint i = 0;

			for ( ; i < numTwos; i++ ) {
				dist += kmerDistances2[(*sKmerCode++)*vocabSize2 + (*tKmerCode++)];

				if ( dist > threshold ) {
					return false;
				}
			}

			if ( rem > 0 ) {
				dist += kmerDistances1[(*sKmerCode++)*vocabSize1 + (*tKmerCode++)];

				if ( dist > threshold ) {
					return false;
				}
			}

			result = dist;
			return true;
		}

	private:
		void PrecomputeDistances() {
			KmerDistanceCache::PrecomputeDistances(1, kmerDistances1, vocabSize1);
			KmerDistanceCache::PrecomputeDistances(2, kmerDistances2, vocabSize2);
		}
	};

	/**
	*	<summary>
	*		Pre-computed kmer distance tables for k == 1, and
	*		kmer distance function which uses them instead of looking
	*		up the matrix to save a lookups and/or loops.
	*	</summary>
	*/
	class KmerDistanceCache1: public KmerDistanceCache {
		// TODO: Define and implement IKmerWordDistance, offering the GetDistance function.
	protected:
		CacheType * kmerDistances1;
		uint vocabSize1;

	public:
		KmerDistanceCache1(Alphabet * alphabet, RawKmerDistanceFunction * dist) : KmerDistanceCache(alphabet, dist) {
			PrecomputeDistances();
		}

		KmerDistanceCache1( const KmerDistanceCache1 & other ) = delete;

		KmerDistanceCache1 & operator=( const KmerDistanceCache1 & other ) = delete;

		virtual ~KmerDistanceCache1() {
			delete[] kmerDistances1;
		}

		size_t CharsPerWord() const {
			return 1;
		}

		Distance GetDistance(KmerWord * sKmerCode, KmerWord * tKmerCode, uint kmerLength) {
			Distance dist = 0;

			uint i;

			for ( i = 0; i < kmerLength; i++ ) {
				dist += kmerDistances1[(*sKmerCode++) * vocabSize1 + (*tKmerCode++)];
			}

			return dist;
		}

		Distance GetDistance1(KmerWord x, KmerWord y) {
			return kmerDistances1[x * vocabSize1 + y];
		}

	private:
		void PrecomputeDistances() {
			KmerDistanceCache::PrecomputeDistances(1, kmerDistances1, vocabSize1);
		}
	};
}
