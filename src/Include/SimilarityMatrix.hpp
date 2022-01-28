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

#include <cmath>
#include <climits>
#include <fstream>
#include <cstdint>

#include "Assert.hpp"
#include "Delegates.hpp"
#include "DistanceType.hpp"
#include "Histogram.hpp"
#include "IntegerDistribution.hpp"
#include "Types.hpp"

namespace QutBio {

#if USE_DOUBLE_DIST
	typedef double Distance;
#define BAD_DIST NAN
#define IS_BAD_DIST(x) std::isnan(x) 
#else
	typedef uint16_t Distance;
#define BAD_DIST numeric_limits<uint16_t>::min()
#define MAX_DIST numeric_limits<uint16_t>::max()
#define IS_BAD_DIST(x) (x == BAD_DIST) 
#endif

	typedef Distance * pDistance;
	typedef pDistance * ppDistance;

	typedef class SimilarityMatrix *pSimilarityMatrix;

	struct SimilarityMatrix {
		int8_t dict[128][128];
		bool isDefined[128];
		string symbols;
		int8_t maxValue = numeric_limits<int8_t>::min();
		int8_t minValue = numeric_limits<int8_t>::max();
		bool isCaseSensitive = false;
		bool isCustom = false;

		void SetSimilarity(unsigned char s, unsigned char t, int8_t value) {
			assert_true(s < 128 && t < 128);

			if (isCaseSensitive) {
				dict[s][t] = value;
				isDefined[s] = true;
			}
			else {
				dict[tolower(s)][tolower(t)] = value;
				dict[toupper(s)][tolower(t)] = value;
				dict[tolower(s)][toupper(t)] = value;
				dict[toupper(s)][toupper(t)] = value;
				isDefined[tolower(s)] = true;
				isDefined[toupper(s)] = true;
			}

			if (maxValue < value) {
				maxValue = value;
			}

			if (minValue > value) {
				minValue = value;
			}

		}

		void Parse(istream & reader) {
			symbols.clear();

			string currentLine;
			int rowCounter = 0;

			int8_t worst = 127;

			while (!(reader.bad() || reader.eof())) {
				getline(reader, currentLine);
				currentLine = String::Trim(currentLine);

				if (currentLine.size() == 0) break;

				if (currentLine[0] == '#') continue;

				auto parts = String::Split(currentLine, ' ');
				size_t len = parts.size();

				if (isalpha(parts[0][0])) {
					// process the heading row

					for (uint i = 0; i < len; i++) {
						char ch = isCaseSensitive ? parts[i][0] : (char)tolower(parts[i][0]);
						symbols += ch;
						isDefined[(uint8_t) ch] = true;
						isDefined[(uint8_t) toupper(ch)] = true;
					}
				}
				else {
					// process a data row.
					for (uint colCounter = 0; colCounter < len; colCounter++) {
						int8_t d = (int8_t)Double::Parse(parts[colCounter]);
						SetSimilarity(symbols[rowCounter], symbols[colCounter], d);

						if (d < worst) {
							worst = d;
						}
					}

					rowCounter++;
				}
			}

			for (int i = 0; i < 128; i++) {
				for (int j = 0; j < 128; j++) {
					if (dict[i][j] == BAD_DIST) {
						dict[i][j] = worst;
					}
				}
			}
		}

		SimilarityMatrix() {
			for (int i = 0; i < 128; i++) {
				isDefined[i] = false;

				for (int j = 0; j < 128; j++) {
					dict[i][j] = BAD_DIST;
				}
			}
		}

		bool IsDefined(unsigned char ch) const {
			return isDefined[ch];
		}

		int8_t Similarity(unsigned char s, unsigned char t) const {
			return dict[s][t];
		}

		int8_t MaxValue(void) const {
			return maxValue;
		}

		int8_t MinValue(void) const {
			return minValue;
		}

		void Foreach(Action3<int, int, Distance> action) const {
			for (int i = 0; i < 128; i++) {
				for (int j = 0; j < 128; j++) {
					if (!IS_BAD_DIST(dict[i][j])) {
						action(i, j, dict[i][j]);
					}
				}
			}
		}

		string Symbols() const {
			return symbols;
		}

		Distance Similarity(const char * x, const char * y, uint length) const {
			Distance score = 0;

			for (uint i = 0; i < length; i++) {
				score += dict[(unsigned char)x[i]][(unsigned char)y[i]];
			}

			return score;
		}

		Distance Similarity(const char * x, uint length) const {
			Distance score = 0;

			for (uint i = 0; i < length; i++) {
				uint t = x[i];
				score += dict[(unsigned char)t][(unsigned char)t];
			}

			return score;
		}

		Distance HalperinDistance(const char * x, const char * y, int length) const {
			Distance d = Similarity(x, length);
			d += Similarity(y, length);
			d -= 2 * Similarity(x, y, length);
			return d;
		}

		Distance Difference(const char * x, const char * y, int length) const {
			Distance d = length * maxValue - Similarity(x, y, length);
			return d;
		}

		Distance Difference(const char x, const char y) const {
			Distance d = maxValue - Similarity(x, y);
			return d;
		}

		bool IsWithin(const char * x, const char * y, int length, Distance threshold, Distance & result) {
			Distance dist = 0;

			for (int i = 0; i < length; i++) {
				dist += maxValue - dict[(uint8_t) x[i]][(uint8_t) y[i]];

				if (dist > threshold) {
					return false;
				}
			}

			result = dist;
			return true;
		}

		static SimilarityMatrix * Blosum100() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum100.iij							\n"
				"#  * column uses minimum score 										\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 100											\n"
				"#  Entropy =   1.4516, Expected =  -1.0948								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -2 -2 -3 -1 -1 -2 -1 -3 -3 -3 -2 -2 -4 -1  1 -1 -4 -4 -1 -3 -2 -1 -7\n"
				"-2  7 -1 -3 -5  0 -2 -4 -1 -4 -4  2 -2 -4 -3 -2 -2 -4 -3 -4 -2 -1 -2 -7\n"
				"-2 -1  7  1 -4 -1 -1 -2  0 -5 -5 -1 -4 -5 -4  0 -1 -6 -3 -4  4 -1 -2 -7\n"
				"-3 -3  1  7 -5 -2  1 -3 -2 -6 -6 -2 -5 -5 -3 -1 -2 -7 -5 -5  4  0 -3 -7\n"
				"-1 -5 -4 -5  9 -5 -6 -5 -5 -2 -3 -5 -3 -3 -5 -2 -2 -5 -4 -2 -5 -6 -3 -7\n"
				"-1  0 -1 -2 -5  7  1 -3  0 -4 -3  1 -1 -4 -2 -1 -2 -3 -3 -3 -1  3 -2 -7\n"
				"-2 -2 -1  1 -6  1  6 -4 -1 -5 -5  0 -4 -5 -3 -1 -2 -5 -4 -3  0  5 -2 -7\n"
				"-1 -4 -2 -3 -5 -3 -4  6 -4 -6 -5 -3 -5 -5 -4 -1 -3 -5 -6 -5 -2 -4 -3 -7\n"
				"-3 -1  0 -2 -5  0 -1 -4  9 -5 -4 -2 -3 -2 -3 -2 -3 -3  1 -5 -1 -1 -2 -7\n"
				"-3 -4 -5 -6 -2 -4 -5 -6 -5  5  1 -4  1 -1 -4 -4 -2 -4 -3  2 -5 -4 -2 -7\n"
				"-3 -4 -5 -6 -3 -3 -5 -5 -4  1  5 -4  2  0 -4 -4 -3 -4 -3  0 -5 -4 -2 -7\n"
				"-2  2 -1 -2 -5  1  0 -3 -2 -4 -4  6 -2 -4 -2 -1 -2 -5 -4 -4 -1  0 -2 -7\n"
				"-2 -2 -4 -5 -3 -1 -4 -5 -3  1  2 -2  8 -1 -4 -3 -2 -3 -3  0 -4 -3 -2 -7\n"
				"-4 -4 -5 -5 -3 -4 -5 -5 -2 -1  0 -4 -1  7 -5 -3 -3  0  3 -2 -5 -5 -3 -7\n"
				"-1 -3 -4 -3 -5 -2 -3 -4 -3 -4 -4 -2 -4 -5  8 -2 -3 -6 -5 -4 -3 -3 -3 -7\n"
				" 1 -2  0 -1 -2 -1 -1 -1 -2 -4 -4 -1 -3 -3 -2  6  1 -4 -3 -3 -1 -1 -1 -7\n"
				"-1 -2 -1 -2 -2 -2 -2 -3 -3 -2 -3 -2 -2 -3 -3  1  6 -5 -3 -1 -2 -2 -1 -7\n"
				"-4 -4 -6 -7 -5 -3 -5 -5 -3 -4 -4 -5 -3  0 -6 -4 -5 11  1 -4 -6 -4 -4 -7\n"
				"-4 -3 -3 -5 -4 -3 -4 -6  1 -3 -3 -4 -3  3 -5 -3 -3  1  8 -3 -4 -4 -3 -7\n"
				"-1 -4 -4 -5 -2 -3 -3 -5 -5  2  0 -4  0 -2 -4 -3 -1 -4 -3  5 -5 -3 -2 -7\n"
				"-3 -2  4  4 -5 -1  0 -2 -1 -5 -5 -1 -4 -5 -3 -1 -2 -6 -4 -5  4  0 -2 -7\n"
				"-2 -1 -1  0 -6  3  5 -4 -1 -4 -4  0 -3 -5 -3 -1 -2 -4 -4 -3  0  4 -2 -7\n"
				"-1 -2 -2 -3 -3 -2 -2 -3 -2 -2 -2 -2 -2 -3 -3 -1 -1 -4 -3 -2 -2 -2 -2 -7\n"
				"-7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7  1";
#pragma endregion
			static SimilarityMatrix matrix;

			if (IS_BAD_DIST(matrix.dict['a']['a'])) {
				istringstream reader(data);
				matrix.Parse(reader);
			}

			return &matrix;
		}

		static SimilarityMatrix * Blosum80() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum80.iij							\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 80											\n"
				"#  Entropy =   0.9868, Expected =  -0.7442								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -2 -2 -2 -1 -1 -1  0 -2 -2 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -1 -1 -6\n"
				"-2  6 -1 -2 -4  1 -1 -3  0 -3 -3  2 -2 -4 -2 -1 -1 -4 -3 -3 -2  0 -1 -6\n"
				"-2 -1  6  1 -3  0 -1 -1  0 -4 -4  0 -3 -4 -3  0  0 -4 -3 -4  4  0 -1 -6\n"
				"-2 -2  1  6 -4 -1  1 -2 -2 -4 -5 -1 -4 -4 -2 -1 -1 -6 -4 -4  4  1 -2 -6\n"
				"-1 -4 -3 -4  9 -4 -5 -4 -4 -2 -2 -4 -2 -3 -4 -2 -1 -3 -3 -1 -4 -4 -3 -6\n"
				"-1  1  0 -1 -4  6  2 -2  1 -3 -3  1  0 -4 -2  0 -1 -3 -2 -3  0  3 -1 -6\n"
				"-1 -1 -1  1 -5  2  6 -3  0 -4 -4  1 -2 -4 -2  0 -1 -4 -3 -3  1  4 -1 -6\n"
				" 0 -3 -1 -2 -4 -2 -3  6 -3 -5 -4 -2 -4 -4 -3 -1 -2 -4 -4 -4 -1 -3 -2 -6\n"
				"-2  0  0 -2 -4  1  0 -3  8 -4 -3 -1 -2 -2 -3 -1 -2 -3  2 -4 -1  0 -2 -6\n"
				"-2 -3 -4 -4 -2 -3 -4 -5 -4  5  1 -3  1 -1 -4 -3 -1 -3 -2  3 -4 -4 -2 -6\n"
				"-2 -3 -4 -5 -2 -3 -4 -4 -3  1  4 -3  2  0 -3 -3 -2 -2 -2  1 -4 -3 -2 -6\n"
				"-1  2  0 -1 -4  1  1 -2 -1 -3 -3  5 -2 -4 -1 -1 -1 -4 -3 -3 -1  1 -1 -6\n"
				"-1 -2 -3 -4 -2  0 -2 -4 -2  1  2 -2  6  0 -3 -2 -1 -2 -2  1 -3 -2 -1 -6\n"
				"-3 -4 -4 -4 -3 -4 -4 -4 -2 -1  0 -4  0  6 -4 -3 -2  0  3 -1 -4 -4 -2 -6\n"
				"-1 -2 -3 -2 -4 -2 -2 -3 -3 -4 -3 -1 -3 -4  8 -1 -2 -5 -4 -3 -2 -2 -2 -6\n"
				" 1 -1  0 -1 -2  0  0 -1 -1 -3 -3 -1 -2 -3 -1  5  1 -4 -2 -2  0  0 -1 -6\n"
				" 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -2 -1 -1 -2 -2  1  5 -4 -2  0 -1 -1 -1 -6\n"
				"-3 -4 -4 -6 -3 -3 -4 -4 -3 -3 -2 -4 -2  0 -5 -4 -4 11  2 -3 -5 -4 -3 -6\n"
				"-2 -3 -3 -4 -3 -2 -3 -4  2 -2 -2 -3 -2  3 -4 -2 -2  2  7 -2 -3 -3 -2 -6\n"
				" 0 -3 -4 -4 -1 -3 -3 -4 -4  3  1 -3  1 -1 -3 -2  0 -3 -2  4 -4 -3 -1 -6\n"
				"-2 -2  4  4 -4  0  1 -1 -1 -4 -4 -1 -3 -4 -2  0 -1 -5 -3 -4  4  0 -2 -6\n"
				"-1  0  0  1 -4  3  4 -3  0 -4 -3  1 -2 -4 -2  0 -1 -4 -3 -3  0  4 -1 -6\n"
				"-1 -1 -1 -2 -3 -1 -1 -2 -2 -2 -2 -1 -1 -2 -2 -1 -1 -3 -2 -1 -2 -1 -1 -6\n"
				"-6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6  1";
#pragma endregion
			static SimilarityMatrix matrix;

			if (IS_BAD_DIST(matrix.dict['a']['a'])) {
				istringstream reader(data);
				matrix.Parse(reader);
			}

			return &matrix;
		}

		static SimilarityMatrix * Blosum62() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum62.iij							\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 62											\n"
				"#  Entropy =   0.6979, Expected =  -0.5209								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n"
				"-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n"
				"-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n"
				"-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n"
				" 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n"
				"-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n"
				"-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n"
				" 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n"
				"-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n"
				"-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n"
				"-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n"
				"-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n"
				"-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n"
				"-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n"
				"-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n"
				" 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n"
				" 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n"
				"-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n"
				"-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n"
				" 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n"
				"-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  0 -1 -4\n"
				"-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  0  4 -1 -4\n"
				" 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n"
				"-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1";
#pragma endregion
			static SimilarityMatrix matrix;

			if (IS_BAD_DIST(matrix.dict['a']['a'])) {
				istringstream reader(data);
				matrix.Parse(reader);
			}

			return &matrix;
		}

		static SimilarityMatrix * Blosum50() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum50.iij							\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/3 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 50											\n"
				"#  Entropy =   0.4808, Expected =  -0.3573								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -2 -1 -2 -1 -1 -1  0 -2 -1 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -1 -1 -5\n"
				"-2  7 -1 -2 -4  1  0 -3  0 -4 -3  3 -2 -3 -3 -1 -1 -3 -1 -3 -1  0 -1 -5\n"
				"-1 -1  7  2 -2  0  0  0  1 -3 -4  0 -2 -4 -2  1  0 -4 -2 -3  4  0 -1 -5\n"
				"-2 -2  2  8 -4  0  2 -1 -1 -4 -4 -1 -4 -5 -1  0 -1 -5 -3 -4  5  1 -1 -5\n"
				"-1 -4 -2 -4 13 -3 -3 -3 -3 -2 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -3 -3 -2 -5\n"
				"-1  1  0  0 -3  7  2 -2  1 -3 -2  2  0 -4 -1  0 -1 -1 -1 -3  0  4 -1 -5\n"
				"-1  0  0  2 -3  2  6 -3  0 -4 -3  1 -2 -3 -1 -1 -1 -3 -2 -3  1  5 -1 -5\n"
				" 0 -3  0 -1 -3 -2 -3  8 -2 -4 -4 -2 -3 -4 -2  0 -2 -3 -3 -4 -1 -2 -2 -5\n"
				"-2  0  1 -1 -3  1  0 -2 10 -4 -3  0 -1 -1 -2 -1 -2 -3  2 -4  0  0 -1 -5\n"
				"-1 -4 -3 -4 -2 -3 -4 -4 -4  5  2 -3  2  0 -3 -3 -1 -3 -1  4 -4 -3 -1 -5\n"
				"-2 -3 -4 -4 -2 -2 -3 -4 -3  2  5 -3  3  1 -4 -3 -1 -2 -1  1 -4 -3 -1 -5\n"
				"-1  3  0 -1 -3  2  1 -2  0 -3 -3  6 -2 -4 -1  0 -1 -3 -2 -3  0  1 -1 -5\n"
				"-1 -2 -2 -4 -2  0 -2 -3 -1  2  3 -2  7  0 -3 -2 -1 -1  0  1 -3 -1 -1 -5\n"
				"-3 -3 -4 -5 -2 -4 -3 -4 -1  0  1 -4  0  8 -4 -3 -2  1  4 -1 -4 -4 -2 -5\n"
				"-1 -3 -2 -1 -4 -1 -1 -2 -2 -3 -4 -1 -3 -4 10 -1 -1 -4 -3 -3 -2 -1 -2 -5\n"
				" 1 -1  1  0 -1  0 -1  0 -1 -3 -3  0 -2 -3 -1  5  2 -4 -2 -2  0  0 -1 -5\n"
				" 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  2  5 -3 -2  0  0 -1  0 -5\n"
				"-3 -3 -4 -5 -5 -1 -3 -3 -3 -3 -2 -3 -1  1 -4 -4 -3 15  2 -3 -5 -2 -3 -5\n"
				"-2 -1 -2 -3 -3 -1 -2 -3  2 -1 -1 -2  0  4 -3 -2 -2  2  8 -1 -3 -2 -1 -5\n"
				" 0 -3 -3 -4 -1 -3 -3 -4 -4  4  1 -3  1 -1 -3 -2  0 -3 -1  5 -4 -3 -1 -5\n"
				"-2 -1  4  5 -3  0  1 -1  0 -4 -4  0 -3 -4 -2  0  0 -5 -3 -4  5  0 -1 -5\n"
				"-1  0  0  1 -3  4  5 -2  0 -3 -3  1 -1 -4 -1  0 -1 -2 -2 -3  0  5 -1 -5\n"
				"-1 -1 -1 -1 -2 -1 -1 -2 -1 -1 -1 -1 -1 -2 -2 -1  0 -3 -1 -1 -1 -1 -1 -5\n"
				"-5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1";
#pragma endregion
			static SimilarityMatrix matrix;

			if (IS_BAD_DIST(matrix.dict['a']['a'])) {
				istringstream reader(data);
				matrix.Parse(reader);
			}

			return &matrix;
		}

		static SimilarityMatrix * Blosum45() {
#pragma region Data
			static const string data = "#  Matrix made by matblas from blosum45.iij\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/3 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 45											\n"
				"#  Entropy =   0.3795, Expected =  -0.2789								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -2 -1 -2 -1 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -2 -2  0 -1 -1  0 -5\n"
				"-2  7  0 -1 -3  1  0 -2  0 -3 -2  3 -1 -2 -2 -1 -1 -2 -1 -2 -1  0 -1 -5\n"
				"-1  0  6  2 -2  0  0  0  1 -2 -3  0 -2 -2 -2  1  0 -4 -2 -3  4  0 -1 -5\n"
				"-2 -1  2  7 -3  0  2 -1  0 -4 -3  0 -3 -4 -1  0 -1 -4 -2 -3  5  1 -1 -5\n"
				"-1 -3 -2 -3 12 -3 -3 -3 -3 -3 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -2 -3 -2 -5\n"
				"-1  1  0  0 -3  6  2 -2  1 -2 -2  1  0 -4 -1  0 -1 -2 -1 -3  0  4 -1 -5\n"
				"-1  0  0  2 -3  2  6 -2  0 -3 -2  1 -2 -3  0  0 -1 -3 -2 -3  1  4 -1 -5\n"
				" 0 -2  0 -1 -3 -2 -2  7 -2 -4 -3 -2 -2 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -5\n"
				"-2  0  1  0 -3  1  0 -2 10 -3 -2 -1  0 -2 -2 -1 -2 -3  2 -3  0  0 -1 -5\n"
				"-1 -3 -2 -4 -3 -2 -3 -4 -3  5  2 -3  2  0 -2 -2 -1 -2  0  3 -3 -3 -1 -5\n"
				"-1 -2 -3 -3 -2 -2 -2 -3 -2  2  5 -3  2  1 -3 -3 -1 -2  0  1 -3 -2 -1 -5\n"
				"-1  3  0  0 -3  1  1 -2 -1 -3 -3  5 -1 -3 -1 -1 -1 -2 -1 -2  0  1 -1 -5\n"
				"-1 -1 -2 -3 -2  0 -2 -2  0  2  2 -1  6  0 -2 -2 -1 -2  0  1 -2 -1 -1 -5\n"
				"-2 -2 -2 -4 -2 -4 -3 -3 -2  0  1 -3  0  8 -3 -2 -1  1  3  0 -3 -3 -1 -5\n"
				"-1 -2 -2 -1 -4 -1  0 -2 -2 -2 -3 -1 -2 -3  9 -1 -1 -3 -3 -3 -2 -1 -1 -5\n"
				" 1 -1  1  0 -1  0  0  0 -1 -2 -3 -1 -2 -2 -1  4  2 -4 -2 -1  0  0  0 -5\n"
				" 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -1 -1  2  5 -3 -1  0  0 -1  0 -5\n"
				"-2 -2 -4 -4 -5 -2 -3 -2 -3 -2 -2 -2 -2  1 -3 -4 -3 15  3 -3 -4 -2 -2 -5\n"
				"-2 -1 -2 -2 -3 -1 -2 -3  2  0  0 -1  0  3 -3 -2 -1  3  8 -1 -2 -2 -1 -5\n"
				" 0 -2 -3 -3 -1 -3 -3 -3 -3  3  1 -2  1  0 -3 -1  0 -3 -1  5 -3 -3 -1 -5\n"
				"-1 -1  4  5 -2  0  1 -1  0 -3 -3  0 -2 -3 -2  0  0 -4 -2 -3  4  0 -1 -5\n"
				"-1  0  0  1 -3  4  4 -2  0 -3 -2  1 -1 -3 -1  0 -1 -2 -2 -3  0  4 -1 -5\n"
				" 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -2 -1 -1 -1 -1 -1 -5\n"
				"-5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1";
#pragma endregion
			static SimilarityMatrix matrix;

			if (IS_BAD_DIST(matrix.dict['a']['a'])) {
				istringstream reader(data);
				matrix.Parse(reader);
			}

			return &matrix;
		}

		static SimilarityMatrix * Blosum40() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum40.iij							\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/4 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 40											\n"
				"#  Entropy =   0.2851, Expected =  -0.2090								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -2 -1 -1 -2  0 -1  1 -2 -1 -2 -1 -1 -3 -2  1  0 -3 -2  0 -1 -1  0 -6\n"
				"-2  9  0 -1 -3  2 -1 -3  0 -3 -2  3 -1 -2 -3 -1 -2 -2 -1 -2 -1  0 -1 -6\n"
				"-1  0  8  2 -2  1 -1  0  1 -2 -3  0 -2 -3 -2  1  0 -4 -2 -3  4  0 -1 -6\n"
				"-1 -1  2  9 -2 -1  2 -2  0 -4 -3  0 -3 -4 -2  0 -1 -5 -3 -3  6  1 -1 -6\n"
				"-2 -3 -2 -2 16 -4 -2 -3 -4 -4 -2 -3 -3 -2 -5 -1 -1 -6 -4 -2 -2 -3 -2 -6\n"
				" 0  2  1 -1 -4  8  2 -2  0 -3 -2  1 -1 -4 -2  1 -1 -1 -1 -3  0  4 -1 -6\n"
				"-1 -1 -1  2 -2  2  7 -3  0 -4 -2  1 -2 -3  0  0 -1 -2 -2 -3  1  5 -1 -6\n"
				" 1 -3  0 -2 -3 -2 -3  8 -2 -4 -4 -2 -2 -3 -1  0 -2 -2 -3 -4 -1 -2 -1 -6\n"
				"-2  0  1  0 -4  0  0 -2 13 -3 -2 -1  1 -2 -2 -1 -2 -5  2 -4  0  0 -1 -6\n"
				"-1 -3 -2 -4 -4 -3 -4 -4 -3  6  2 -3  1  1 -2 -2 -1 -3  0  4 -3 -4 -1 -6\n"
				"-2 -2 -3 -3 -2 -2 -2 -4 -2  2  6 -2  3  2 -4 -3 -1 -1  0  2 -3 -2 -1 -6\n"
				"-1  3  0  0 -3  1  1 -2 -1 -3 -2  6 -1 -3 -1  0  0 -2 -1 -2  0  1 -1 -6\n"
				"-1 -1 -2 -3 -3 -1 -2 -2  1  1  3 -1  7  0 -2 -2 -1 -2  1  1 -3 -2  0 -6\n"
				"-3 -2 -3 -4 -2 -4 -3 -3 -2  1  2 -3  0  9 -4 -2 -1  1  4  0 -3 -4 -1 -6\n"
				"-2 -3 -2 -2 -5 -2  0 -1 -2 -2 -4 -1 -2 -4 11 -1  0 -4 -3 -3 -2 -1 -2 -6\n"
				" 1 -1  1  0 -1  1  0  0 -1 -2 -3  0 -2 -2 -1  5  2 -5 -2 -1  0  0  0 -6\n"
				" 0 -2  0 -1 -1 -1 -1 -2 -2 -1 -1  0 -1 -1  0  2  6 -4 -1  1  0 -1  0 -6\n"
				"-3 -2 -4 -5 -6 -1 -2 -2 -5 -3 -1 -2 -2  1 -4 -5 -4 19  3 -3 -4 -2 -2 -6\n"
				"-2 -1 -2 -3 -4 -1 -2 -3  2  0  0 -1  1  4 -3 -2 -1  3  9 -1 -3 -2 -1 -6\n"
				" 0 -2 -3 -3 -2 -3 -3 -4 -4  4  2 -2  1  0 -3 -1  1 -3 -1  5 -3 -3 -1 -6\n"
				"-1 -1  4  6 -2  0  1 -1  0 -3 -3  0 -3 -3 -2  0  0 -4 -3 -3  5  0 -1 -6\n"
				"-1  0  0  1 -3  4  5 -2  0 -4 -2  1 -2 -4 -1  0 -1 -2 -2 -3  0  5 -1 -6\n"
				" 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1  0 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -6\n"
				"-6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6  1";
#pragma endregion
			static SimilarityMatrix matrix;

			if (IS_BAD_DIST(matrix.dict['a']['a'])) {
				istringstream reader(data);
				matrix.Parse(reader);
			}

			return &matrix;
		}

		static SimilarityMatrix * Blosum35() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum35.iij							\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/4 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 35											\n"
				"#  Entropy =   0.2111, Expected =  -0.1550								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -1 -1 -1 -2  0 -1  0 -2 -1 -2  0  0 -2 -2  1  0 -2 -1  0 -1 -1  0 -5\n"
				"-1  8 -1 -1 -3  2 -1 -2 -1 -3 -2  2  0 -1 -2 -1 -2  0  0 -1 -1  0 -1 -5\n"
				"-1 -1  7  1 -1  1 -1  1  1 -1 -2  0 -1 -1 -2  0  0 -2 -2 -2  4  0  0 -5\n"
				"-1 -1  1  8 -3 -1  2 -2  0 -3 -2 -1 -3 -3 -1 -1 -1 -3 -2 -2  5  1 -1 -5\n"
				"-2 -3 -1 -3 15 -3 -1 -3 -4 -4 -2 -2 -4 -4 -4 -3 -1 -5 -5 -2 -2 -2 -2 -5\n"
				" 0  2  1 -1 -3  7  2 -2 -1 -2 -2  0 -1 -4  0  0  0 -1  0 -3  0  4 -1 -5\n"
				"-1 -1 -1  2 -1  2  6 -2 -1 -3 -1  1 -2 -3  0  0 -1 -1 -1 -2  0  5 -1 -5\n"
				" 0 -2  1 -2 -3 -2 -2  7 -2 -3 -3 -1 -1 -3 -2  1 -2 -1 -2 -3  0 -2 -1 -5\n"
				"-2 -1  1  0 -4 -1 -1 -2 12 -3 -2 -2  1 -3 -1 -1 -2 -4  0 -4  0 -1 -1 -5\n"
				"-1 -3 -1 -3 -4 -2 -3 -3 -3  5  2 -2  1  1 -1 -2 -1 -1  0  4 -2 -3  0 -5\n"
				"-2 -2 -2 -2 -2 -2 -1 -3 -2  2  5 -2  3  2 -3 -2  0  0  0  2 -2 -2  0 -5\n"
				" 0  2  0 -1 -2  0  1 -1 -2 -2 -2  5  0 -1  0  0  0  0 -1 -2  0  1  0 -5\n"
				" 0  0 -1 -3 -4 -1 -2 -1  1  1  3  0  6  0 -3 -1  0  1  0  1 -2 -2  0 -5\n"
				"-2 -1 -1 -3 -4 -4 -3 -3 -3  1  2 -1  0  8 -4 -1 -1  1  3  1 -2 -3 -1 -5\n"
				"-2 -2 -2 -1 -4  0  0 -2 -1 -1 -3  0 -3 -4 10 -2  0 -4 -3 -3 -1  0 -1 -5\n"
				" 1 -1  0 -1 -3  0  0  1 -1 -2 -2  0 -1 -1 -2  4  2 -2 -1 -1  0  0  0 -5\n"
				" 0 -2  0 -1 -1  0 -1 -2 -2 -1  0  0  0 -1  0  2  5 -2 -2  1 -1 -1  0 -5\n"
				"-2  0 -2 -3 -5 -1 -1 -1 -4 -1  0  0  1  1 -4 -2 -2 16  3 -2 -3 -1 -1 -5\n"
				"-1  0 -2 -2 -5  0 -1 -2  0  0  0 -1  0  3 -3 -1 -2  3  8  0 -2 -1 -1 -5\n"
				" 0 -1 -2 -2 -2 -3 -2 -3 -4  4  2 -2  1  1 -3 -1  1 -2  0  5 -2 -2  0 -5\n"
				"-1 -1  4  5 -2  0  0  0  0 -2 -2  0 -2 -2 -1  0 -1 -3 -2 -2  5  0 -1 -5\n"
				"-1  0  0  1 -2  4  5 -2 -1 -3 -2  1 -2 -3  0  0 -1 -1 -1 -2  0  4  0 -5\n"
				" 0 -1  0 -1 -2 -1 -1 -1 -1  0  0  0  0 -1 -1  0  0 -1 -1  0 -1  0 -1 -5\n"
				"-5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1";
#pragma endregion
			static SimilarityMatrix matrix;

			if (IS_BAD_DIST(matrix.dict['a']['a'])) {
				istringstream reader(data);
				matrix.Parse(reader);
			}

			return &matrix;
		}

		static SimilarityMatrix * GetBlosum(int matrixId) {
			switch (matrixId) {
			case 100: return Blosum100();
			case 80: return Blosum80();
			case 62: return Blosum62();
			case 50: return Blosum50();
			case 45: return Blosum45();
			case 40: return Blosum40();
			case 35: return Blosum35();
			default: return 0;
			}
		}

		/**
		 *	<summary>
		 *		Computes the value of the Karlin-Altshul parameter lambda for stipulated
		 *		background distributions.
		 *	</summary>
		 *	<param name=p1>
		 *		A normalised histogram  which contains the distribution of query residues.
		 *	</summary>
		 *	<param name=p2>
		 *		A normalised histogram  which contains the distribution of subject residues.
		 *	</summary>
		 */
		double ComputeExtremeDistLambda(
			const Histogram<char> & p1,
			const Histogram<char> & p2
		) {
			auto f = [=](double lambda) {
				double sum = 0;

				for (auto x : p1.data) {
					for (auto y : p2.data) {
						double score = dict[(uint8_t) x.first][(uint8_t) y.first];
						sum += x.second * y.second * exp(lambda * score);
					}
				}

				return sum;
			};

			auto fPrime = [=](double lambda) {
				double sum = 0;

				for (auto x : p1.data) {
					for (auto y : p2.data) {
						double score = dict[(uint8_t) x.first][(uint8_t)y.first];
						sum += x.second * y.second * score * exp(lambda * score);
					}
				}

				return sum;
			};

			double lambda = 0;

			while (true) {
				double lambda2 = lambda - f(lambda) / fPrime(lambda);
				const double eps = 1e-10;

				if (fabs(lambda2 - lambda) < eps) break;

				lambda = lambda2;
			}

			return lambda;
		}

		/**
		 *	<summary>
		 *		Gets either an in-built or custom similarity matrix defined by the supplied arguments.
		 *
		 *		Beware: if Custom is used, a new instance is created via the new operator. You will have to free this yourself when finished with it.
		 *	</summary>
		 */
		static SimilarityMatrix * GetMatrix(
			DistanceType * dist,
			int id,
			const string & customFileName,
			bool isCaseSensitive
		) {
			if (dist == DistanceType::HalperinEtAl() || dist == DistanceType::BlosumDistance()) {
				return GetBlosum(id);
			}
			else if (dist == DistanceType::Custom()) {
				SimilarityMatrix * matrix = new SimilarityMatrix();
				matrix->isCaseSensitive = isCaseSensitive;
				matrix->isCustom = true;

				if (matrix == 0) {
					throw Exception("Out of memory!", FileAndLine);
				}

				ifstream f(customFileName);
				matrix->Parse(f);
				return matrix;
			}
			else {
				throw Exception("Unable to create instance of matrix type " + dist->Name(), FileAndLine);
			}
		}

		void GetDifferenceDistributions(
			Histogram<char> & symbolDistribution,
			int maxK,
			vector<pair<int, IntegerDistribution>> &distributions
		) {
			Histogram<int> h;
			h.GetOneMerHistogram<char, int>(
				symbolDistribution,
				[=](char x, char y) {return Difference(x, y); }
			);

			assert_intsEqual(0, h.data.begin()->first);
			assert_intsEqual(maxValue - minValue, h.data.rbegin()->first);

			int min = 0;
			int max = maxValue - minValue;

			vector<double> p(1 + max - min);

			for (auto kvp : h.data) {
				p[(int)kvp.first - min] = kvp.second;
			}

			IntegerDistribution oneMerDistances(min, max, p);
			IntegerDistribution *latestDistribution = &oneMerDistances;
			distributions.push_back(pair<int, IntegerDistribution>(1, oneMerDistances));

			for (int k = 2; k <= maxK; k++) {
				IntegerDistribution nextDistribution = oneMerDistances.Add(*latestDistribution);
				distributions.push_back(pair<int, IntegerDistribution>(k, nextDistribution));
				latestDistribution = &oneMerDistances;
			}
		}

		bool IsCustom() {
			return this->isCustom;
		}

		void PopulateDistanceTable( Distance lookup[128][128] ) {
			for ( int i = 0; i < 128; i++ ) {
				for ( int j = 0; j < 128; j++ ) {
					lookup[i][j] = maxValue - dict[i][j];
				}
			}
		}
	};

	typedef SimilarityMatrix * pSimilarityMatrix;
}
