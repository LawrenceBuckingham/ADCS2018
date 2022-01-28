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

#include <map>
#include <vector>
#include <iostream>

#include "CsvIO.hpp"

using std::map;

namespace QutBio {

	template < typename T >
	class Histogram {
	public:
		/// <summary>The actual histogram.</summary>
		/// TODO: protect this.
		map<T, double> data;

		template <class Q>
		void Initialise(Q & values) {
			data.clear();
			AddRange(values);
		}

		void Add(T x) {
			auto value = data.find(x);
			if ( value == data.end() ) {
				data[x] = 1;
			}
			else {
				value->second++;
			}
		}

		void Add(T x, double y) {
			auto value = data.find(x);
			if ( value == data.end() ) {
				data[x] = y;
			}
			else {
				value->second += y;
			}
		}

		template <class Q>
		void AddRange(Q & values) {
			for ( auto x : values ) {
				auto value = data.find(x);
				if ( value == data.end() ) {
					data[x] = 1;
				}
				else {
					value->second++;
				}
			}
		}

		void DoConvolution(
			const Histogram<T> & singleHistogram,
			Histogram<T> & newHistogram
			) {
			for ( auto i : data ) {
				double currentKey = i.first;
				double currentValue = i.second;

				for ( auto j : singleHistogram.data ) {
					double singleKey = j.first;
					double singleValue = j.second;
					double newKey = currentKey + singleKey;
					newHistogram.Add(newKey, currentValue * singleValue);
				}
			}
		}

		/**
		*	<summary>
		*		Populates this histogram with a normalised histogram (i.e. a probability mass function)
		*		of the values obtained by applying a pairwise operator (for example, similarity or
		*		difference) to the Cartesian product of a given symbol set with itself. That is, it is 
		*		the distribution of 1-mer pairwise values under the assumption that the symbols are
		*		drawn i.i.d. from a uniform distribution.
		*	<para>
		*		On completion, data[x] = count { a in alphabet, b in alphabet, f(a,b) = x } /
		*			count{ a in alphabet, b in alphabet}.
		*	</para>
		*	</summary>
		*	<param name="alphabet>The alphabet.</param>
		*	<param name="f">A function that maps a pair of symbols from the alphabet to a</param>
		*/
		template<typename C, typename U>
		void GetOneMerHistogram(
			const C & alphabet,
			Func2<U, U, T> f
			) {
			data.clear();

			for ( auto x : alphabet ) {
				for ( auto y : alphabet ) {
					double d = f(x, y);
					auto p = data.find(d);

					if ( p == data.end() ) {
						data[d] = 1;
					}
					else {
						p->second++;
					}
				}
			}

			Normalise();
		}

		/**
		 *	<summary>
		 *		Populates this histogram with a normalised histogram (i.e. a probability mass function)
		 *		of the values obtained by applying a pairwise operator (for example, similarity or
		 *		difference) to the weighted Cartesian product of a given symbol distribution with itself.
		 *	<para>
		 *		On completion, data[x] = sum { p(a)p(b) | a in alphabet, b in alphabet, f(a,b) = x } /
		 *			sum{ p(a)p(b) | a in alphabet, b in alphabet}.
		 *	</para>
		 *	</summary>
		 *	<param name="alphabetDistribution>The symbol distribution (a probability mass function).</param>
		 *	<param name="f">A function that maps a pair of symbols from the alphabet to a</param>
		 */
		template<typename C, typename U>
		void GetOneMerHistogram(
			const Histogram<C> & alphabetDistribution,
			Func2<U, U, T> f
			) {
			data.clear();

			for ( auto x : alphabetDistribution.data ) {
				for ( auto y : alphabetDistribution.data ) {
					T d = f(x.first, y.first);
					auto p = data.find(d);

					double weight = x.second * y.second;

					if ( p == data.end() ) {
						data[d] = weight;
					}
					else {
						p->second += weight;
					}
				}
			}

			Normalise();
		}

		double operator[](T t) const {
			auto loc = data.find( t );

			return loc == data.end() ? 0 : loc->second;
		}

		/// <summary>
		/// Emits the distribution to standard output in a tabular form.
		/// </summary>
		ostream & Print(ostream & out) {
			out << "x" << "\t" << "f" << endl;

			for ( auto k : data ) {
				T key = k.first;
				double f = k.second;
				out << key << "\t" << f << endl;
			}

			return out;
		}

		/// <summary>
		/// Emits the distribution to standard output in a tabular form.
		/// </summary>
		/// <param name="keyFormat">A function that converts the key to string format.</param>
		/// <param name="keyFormat">A function that converts the key to string format.</param>
		/// <returns>A function that converts the key to string format.</returns>

		ostream & Print(
			ostream & out,
			Func1<T, string> keyFormat,
			Func1<double, string> valFormat
			) {
			out << "x" << "," << "f" << endl;

			for ( auto k : data ) {
				T key = k.first;
				double f = k.second;
				out << keyFormat(key) << "," << valFormat(f) << endl;
			}

			return out;
		}

		/**
		 *	Parses a histogram from a stream in which the x-y pairs are presented as a
		 *	Nx2 matrix, with x in column 1 and y in column 2.
		 */

		void ParseCols( istream & in, char delimiter, Func1<string &, T> parser ) {
			data.clear();
			CsvReader csv( in, delimiter );
			vector<vector<string>> records;

			csv.Read( records );

			for ( auto record : records ) {
				if ( record.size() != 2 ) continue;

				if ( !isdigit( record[1][0] ) ) continue;

				T key = parser( record[0] );
				double value = atof( record[1].c_str() );
				data[key] = value;
			}
		}

		/**
		 *	Parses a histogram from a stream in which the x-y pairs are presented as a
		 *	2xN matrix, with x values in row 1 and y values in row 2.
		 */

		void ParseRows( istream & in, char delimiter, Func1<string &, T> parser ) {
			data.clear();
			CsvReader csv( in, delimiter );
			vector<vector<string>> records;

			csv.Read( records );

			vector<string> & xvals = records[0];
			vector<string> & yvals = records[1];

			for ( int i = 1; i < xvals.size() && i < yvals.size(); i++ ) {
				T key = parser( xvals[i] );
				double val = Double::Parse( yvals[i] );
				data[key] = val;
			}
		}

		/// <summary>
		///	Scales the total mass of the histogram to 1.
		/// </summary>
		void Normalise() {
			double total_weight = 0;

			for ( auto i : data ) {
				total_weight += i.second;
			}

			for ( auto i : data ) {
				data[i.first] /= total_weight;
			}
		}

		bool Equals(const Histogram<T> & other, double tolerance = 1e-10) {
			if ( data.size() != other.data.size() ) {
				return false;
			}

			//Histogram<char> hhh;
			//auto z = hhh.data.begin();

			auto this_iter = data.begin();
			auto other_iter = other.data.begin();

			while ( this_iter != data.end() && other_iter != other.data.end() ) {
				if ( this_iter->first != other_iter->first ) {
					return false;
				}
				if ( fabs(this_iter->second - other_iter->second) > tolerance ) {
					return false;
				}
				this_iter++;
				other_iter++;
			}

			return true;
		}


		/// <summary>
		///     Eliminates unnecessary entries for which p == 0.
		/// </summary>

		void Cleanup(Func2<const T &, double, bool> predicate) {
			vector<double> toRemove;

			for ( auto kvp : data ) {
				if ( predicate(kvp.first, kvp.second) ) toRemove.push_back(kvp.first);
			}

			for ( auto key : toRemove ) {
				data.erase(key);
			}
		}

		/**
		*	<summary>
		*		Repopulates a vector with an ordered list of keys belonging to this histogram.
		*	</summary>
		*/

		void GetKeys(vector<T> & keys) {
			keys.clear();

			for ( auto p : data ) {
				keys.push_back(p.first);
			}
		}

		/**
		*	<summary>
		*		Repopulates a vector with an ordered list of keys belonging to this histogram.
		*	</summary>
		*/

		void GetValues(vector<double> & values) {
			values.clear();

			for ( auto p : data ) {
				values.push_back(p.second);
			}
		}
	};
}

