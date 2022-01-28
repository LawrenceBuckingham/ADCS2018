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

#include <string>
#include <cstdint>
#include <cstddef>

#include "EncodedKmer.hpp"
#include "EnumBase.hpp"
#include "SimilarityMatrix.hpp"
#include "Array.hpp"
#include "BitSet.hpp"

using namespace QutBio;

namespace QutBio {
	typedef class Alphabet *pAlphabet;

	class Alphabet : public EnumBase {
	private:
		string symbols;
		unsigned char inverse[128];

		Alphabet(string literal, int value)
			: EnumBase(literal, value) {
			symbols = value == 0 ? SimilarityMatrix::Blosum62()->Symbols() : string("acgt");

			std::fill(inverse, inverse + sizeof(inverse), 0);

			for (size_t i = 0; i < symbols.size(); i++) {
				inverse[tolower(symbols[i])] = (unsigned char)i;
				inverse[toupper(symbols[i])] = (unsigned char)i;
			}
		}

	public:
		Alphabet(SimilarityMatrix * matrix) : EnumBase("Custom", 2) {
			symbols = matrix->Symbols();

			std::fill(inverse, inverse + sizeof(inverse), 0);

			for (size_t i = 0; i < symbols.size(); i++) {
				inverse[tolower(symbols[i])] = (unsigned char)i;
				inverse[toupper(symbols[i])] = (unsigned char)i;
			}
		}

		static Alphabet * AA() {
			static Alphabet value("AA", 0);
			return &value;
		}

		static Alphabet * DNA() {
			static Alphabet value("DNA", 1);
			return &value;
		}

		static vector<EnumBase *> & Values() {
			static vector<EnumBase *> result(2);
			result[0] = AA();
			result[1] = DNA();
			return result;
		}

		static size_t WordsPerKmer(size_t kmerLength, size_t charsPerWord) {
			return (kmerLength + charsPerWord - 1) / charsPerWord;
		}

		size_t BitsPerSymbol() {
			size_t bitsPerSymbol = 1;

			while (((size_t)1 << bitsPerSymbol) < symbols.size()) {
				bitsPerSymbol++;
			}

			return bitsPerSymbol;
		}

		/// <summary> Encodes a single kmer as a sequence of zero-origin numeric values.
		/// </summary>
		/// <param name="s"></param>
		/// <param name="charsPerWord">The number of encoded characters to pack into each word.
		///		This is independent of any particular kmer tiling scheme. It has to do with the
		///		maximum size of a pre-computed kmer similarity lookup table that can be created.
		///		<para>
		///			For my current experiments, I use charPerWord == 2 or sometimes 3 for amino acid, with the
		///			BLOSUM alphabet, this yields tables of approx 1/4GB.
		///		</para>
		/// </param>
		/// <returns></returns>

		void Encode(const char * s, size_t kmerLength, size_t charsPerWord, EncodedKmer code) {
			size_t words = WordsPerKmer(kmerLength, charsPerWord);
			size_t size = symbols.size();
			size_t wordIndex = 0;
			code[wordIndex] = 0;

			for (size_t i = 0; i < kmerLength; i++) {
				int j = inverse[(uint8_t)s[i]];
				code[wordIndex] = code[wordIndex] * (int)size + j;

				if (i > 0 && i % charsPerWord == (charsPerWord - 1)) {
					wordIndex++;

					if (wordIndex < words) {
						code[wordIndex] = 0;
					}
				}
			}
		}

		/// <summary> Encodes the designated string by computing a codeword for each tuple of 
		///		charsPerWord contiguous symbols. Multiple lists of codewords are calculated,
		///		one for each offset in 0, ..., (charsPerWord-1), and each of these lists becomes 
		///		a row of the resulting output dataset, code. 
		///		
		///		After processing:
		///		*	code[offset] contains the codewords of tuples starting at positions
		///			i*charsPerWord + offset for i = 1..(len / charsPerWord).
		/// </summary>
		/// <param name="s">The character sequence to encode.</param>
		/// <param name="charsPerWord">The number of encoded characters to pack into each word.
		///		This is independent of any particular kmer tiling scheme. It has to do with the
		///		maximum size of a pre-computed kmer similarity lookup table that can be created.
		///		<para>
		///			For my current experiments, I use charsPerWord \in {2,3} for amino acid, with the
		///			BLOSUM alphabet, this yields in-memory distance tables of approx 1/4GB.
		///		</para>
		///		<para>
		///			To encode DNA, I place one kmer per codeword, which introduces a constraint 
		///			on the maximum length, which becomes CHAR_BIT * sizeof(KmerWord) / 2.
		///		</para>
		/// </param>
		/// <param name="len">The number of characters to encode.</param>
		/// <param name="kmerLength">The number of characters in each kmer.</param>
		///	<param name="code">
		///		A reference to a list of arrays which will hold the codes 
		///		starting at each possible offset: 0 .. (charsPerWord - 1).
		///		<para>
		///			The staggered arrangement allows a sequence of successive kmers to be read from
		///			contiguous memory locations, so that code[0] is a list of partial kmers which start 
		///			at {0, charsPerWord, 2*charsPerWord, ...}; code[1] is the list of partial k-mers that
		///			start at {1, 1+charsPerWord, 1+2*charsPerWord, ...}.
		///		</para>
		///		<para>
		///			For DNA, where a KmerWord contains a whole k-mer, code will contain just a single
		///			row. That is, code[0] will contain the list of packed kmers.
		///		</para>
		///	</param>
		/// <returns></returns>

		void Encode(
			const char * s,
			size_t len,
			size_t kmerLength,
			size_t charsPerWord,
			vector<vector<KmerWord>> & code
		) {
			if (len + 1 < kmerLength ) {
				stringstream str;
				str << "Alphabet::Encode: string to encode must contain at least one k-mer: kmerLength = "
					<< kmerLength
					<< ", len = "
					<< len
					<< "\n";
				throw Exception(str.str(), FileAndLine);
			}

			if (kmerLength > charsPerWord) {
				if (kmerLength % charsPerWord != 0) {
					stringstream str;
					str << "Alphabet::Encode: kmerLength must be divisible by charsPerWord: kmerLength = "
						<< kmerLength
						<< ", charsPerWord = "
						<< charsPerWord
						<< "\n";
					throw Exception(str.str(), FileAndLine);
				}

				code.resize(charsPerWord);

				for (size_t i = 0; i < charsPerWord; i++) {
					code[i].clear();
					code[i].reserve(len / charsPerWord);
				}

				for (size_t i = 0; i < len - charsPerWord + 1; i++) {
					KmerWord codeWord;
					Encode(s + i, charsPerWord, charsPerWord, &codeWord);
					code[i % charsPerWord].push_back(codeWord);
				}
			}
			else {
				// This is primarily for DNA;
				code.resize(1);
				code[0].clear();

				for (size_t i = 0; i < len - kmerLength + 1; i++) {
					KmerWord codeWord;
					Encode(s + i, kmerLength, charsPerWord, &codeWord);
					code[0].push_back(codeWord);
				}
			}
		}

		/// <summary> Decodes a sequence of zero-origin numeric values into a string.
		/// </summary>
		/// <param name="code">The sequence of numeric code values.</param>
		/// <param name="k">The number of symbols required in the decoded string.</param>
		/// <param name="charsPerWord">The number of characters packed into each 
		///		word. See Encode for details.
		///	</param>
		/// <returns></returns>

		void Decode(const KmerWord * code, size_t k, size_t charsPerWord, char * charBuffer) const {
			size_t words = (k + charsPerWord - 1) / charsPerWord;
			size_t size = symbols.size();
			string * s = new string();
			vector<string *> t = {};
			t.push_back(s);
			size_t wordIndex = 0;
			KmerWord n = code[wordIndex];

			for (size_t i = 0; i < k; i++) {
				size_t j = n % size;
				n = (n - j) / size;
				s->insert(0, 1, symbols[j]);

				if (i > 0 && i % charsPerWord == (charsPerWord - 1)) {
					wordIndex++;

					if (wordIndex < words) {
						s = new string();
						t.push_back(s);
						n = code[wordIndex];
					}
				}
			}

			int buffPos = 0;

			for (size_t i = 0, iMax = t.size(); i < iMax; i++) {
				for (size_t j = 0, jMax = t[i]->size(); j < jMax; j++) {
					charBuffer[buffPos++] = (*t[i])[j];
				}
			}
		}

		/// <summary> Decodes a sequence of zero-origin numeric values into a string.
		/// </summary>
		/// <param name="code">The sequence of numeric code values.</param>
		/// <param name="len">The number of symbols required in the decoded string.</param>
		/// <param name="charsPerWord">The number of characters packed into each 
		///		word. See Pack for details.
		///	</param>
		/// <returns></returns>

		void Decode(const vector<vector<KmerWord>> & code, size_t len, size_t kmerLength, size_t charsPerWord, char * charBuffer) const {
			if (kmerLength > charsPerWord) {
				for (size_t i = 0; i < len - charsPerWord + 1; i++) {
					Decode(&code[i%charsPerWord][i / charsPerWord], charsPerWord, charsPerWord, charBuffer + i);
				}
			}
			else {
				for (size_t i = 0; i < len; i++) {
					Decode(&code[0][i], kmerLength, charsPerWord, charBuffer + i);
				}

			}
			charBuffer[len] = 0;
		}

		/// <summary> Gets the number of symbols in this alphabet.
		/// </summary>

		int Size() const {
			return (int)symbols.size();
		}

		/// <summary> Gets the symbols defined in this alphabet.
		/// </summary>

		string Symbols() const {
			return symbols;
		}

		/// <summary>Gets a symbol which can be used to pad strings.</summary>

		char DefaultSymbol() {
			auto idx = symbols.find_first_of("xX");

			return idx = string::npos ? symbols.front() : symbols[idx];
		}

		const uint8_t * Inverse() const {
			return inverse;
		}
	};

	typedef Alphabet * pAlphabet;
}
