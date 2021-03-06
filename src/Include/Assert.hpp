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
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <vector>
#include <iomanip>

#include "Exception.hpp"

using std::function;
using std::cerr;

namespace QutBio {

	class Assert {
	public:
		static void IsTrue(bool cond, const char * file, int line) {
			if (!cond) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				throw Exception("Condition is not true as expected.", file, line);
			}
		}

		static void IsTrue(bool cond, const char * file, int line, function<void(void)> showDetails) {
			if (!cond) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				showDetails();
				throw Exception("Condition is not true as expected.", file, line);
			}
		}

		static void IsFalse(bool cond, string file, int line) {
			if (cond) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				throw Exception("Condition is not false as expected.", file, line);
			}
		}

		static void IsFalse(bool cond, string file, int line, function<void(void)> showDetails) {
			if (cond) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				showDetails();
				throw Exception("Condition is not false as expected.", file, line);
			}
		}

		static void StringsEqual(const string & expected, const string & actual, const char * file, int line) {
			if (expected != actual) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				ostringstream out;
				out << "Expected value <" << expected << "> does not match actual value <" << actual << ">.";
				throw Exception(out.str(), file, line);
			}
		}

		//static void StringsEqual( string expected, string actual, const char * file, int line ) {
		//	if ( expected != actual ) {
		//		ostringstream out;
		//		out << "Expected value <" << expected << "> does not match actual value <" << actual << ">.";
		//		throw Exception( out.str(), file, line );
		//	}
		//}

		static void StringsEqual(const char * expected, const char * actual, const char * file, int line) {
			if (strcmp(expected, actual) != 0) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				ostringstream out;
				out << "Expected value <" << expected << "> does not match actual value <" << actual << ">.";
				throw Exception(out.str(), file, line);
			}
		}

		static void StringsEqual(const char * expected, const char * actual, const char * file, int line, function<void(void)> showDetails) {
			if (strcmp(expected, actual) != 0) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				ostringstream out;
				out << "Expected value <" << expected << "> does not match actual value <" << actual << ">.";
				showDetails();
				throw Exception(out.str(), file, line);
			}
		}

		template<typename T, typename U>
		static void IntsEqual(T expected, U actual, const char * file, int line) {
			if (expected != actual) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				ostringstream out;
				out << "Expected value <" << expected << "> does not match actual value <" << actual << ">.";
				throw Exception(out.str(), file, line);
			}
		}

		template<typename T, typename U>
		static void IntsEqual(T expected, U actual, const char * file, int line, function<void(void)> showDetails) {
			if (expected != actual) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				ostringstream out;
				out << "Expected value <" << expected << "> does not match actual value <" << actual << ">.";
				showDetails();
				throw Exception(out.str(), file, line);
			}
		}

		static void DoublesEqual(double expected, double actual, double epsilon, const char * file, int line) {
			if (fabs(expected - actual) >= epsilon) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				ostringstream out;
				out << "Expected value <" << std::setprecision(15) << std::setw(15) << expected << "> does not match actual value <" << actual << ">.";
				throw Exception(out.str(), file, line);
			}
		}

		static void DoublesEqual(double expected, double actual, double epsilon, const char * file, int line, function<void(void)> showDetails) {
			if (fabs(expected - actual) >= epsilon) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				ostringstream out;
				out << "Expected value <" << std::setprecision(15) << std::setw(15) << expected << "> does not match actual value <" << actual << ">.";
				showDetails();
				throw Exception(out.str(), file, line);
			}
		}

		template<typename T>
		static void VectorsEqual(const vector<T> & expected, const vector<T> & actual, const char * file, int line) {
			if (expected.size() != actual.size()) {
				cerr << file << "(" << line << ") - Assertion failed.\n";
				ostringstream buffer;
				buffer << "Vector lengths do not match: expected = " << expected.size() << ", actual = " << actual.size();
				throw Exception(buffer.str(), file, line);
			}

			for (size_t i = 0, len = expected.size(); i < len; i++) {
				if (expected[i] != actual[i]) {
					ostringstream buffer;
					cerr << file << "(" << line << ") - Assertion failed.\n";
					buffer << "Elements in list at position " << i << " do not match: expected = " << expected.size() << ", actual = " << actual.size();
					throw Exception(buffer.str(), file, line);
				}
			}
		}

		template<typename T>
		static void Equal(const T & expected, const T & actual, const char * file, int line) {
			if (expected != actual) {
					cerr << file << "(" << line << ") - Assertion failed.\n";
				ostringstream buffer;
				buffer << "Expected and actual values do not match: expected = " << expected << ", actual = " << actual;
				throw Exception(buffer.str(), file, line);
			}

		}
	};


#define assert_true( cond ) Assert::IsTrue( cond, __FILE__, __LINE__ )
#define assert_false( cond ) Assert::IsFalse( cond, __FILE__, __LINE__ )
#define assert_stringsEqual( s1, s2 ) Assert::StringsEqual(s1, s2, __FILE__, __LINE__ )
#define assert_intsEqual( i1, i2 ) Assert::IntsEqual(i1, i2, __FILE__, __LINE__ )
#define assert_doublesEqual( i1, i2, eps ) Assert::DoublesEqual(i1, i2, eps, __FILE__, __LINE__ )
#define assert_vectorsEqual( i1, i2 ) Assert::VectorsEqual(i1, i2, __FILE__, __LINE__ )
#define assert_equal( i1, i2 ) Assert::Equal(i1, i2, __FILE__, __LINE__ )

}
