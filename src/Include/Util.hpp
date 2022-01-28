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

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <cstdint>
#include <functional>
#include <cmath>
#include <vector>

#include "Types.hpp"
#include "Exception.hpp"

namespace QutBio {

	using string = std::string;
	using istringstream = std::istringstream;
	using ostringstream = std::ostringstream;
	using istream = std::istream;
	using ifstream = std::ifstream;
	using ostream = std::ostream;

	template <typename T>
	T min( T x, T y ) { return x < y ? x : y; }

	template <typename T>
	T max( T x, T y ) { return x > y ? x : y; }

	template<typename T>
	static T sum( typename vector<T>::iterator begin, typename vector<T>::iterator end, T initial ) {
		T total = initial;

		for ( auto iter = begin; iter != end; iter++ ) {
			total += *iter;
		}

		return total;
	}

	class Util {
	public:

		template<typename X, typename Y>
		static void LinFit( const X & x, const Y & y, size_t n, double & a, double &b ) {
			double sumX = 0, sumY = 0, sumXY = 0, sumXSquared = 0;

			for ( size_t i = 0; i < n; i++ ) {
				if ( !std::isfinite( y[i] ) ) continue;
				sumX += x[i];
				sumY += y[i];
				sumXY += x[i] * y[i];
				sumXSquared += x[i] * x[i];
			}

			double A = sumXSquared, B = sumX, C = sumX, D = n;
			double det = A * D - B * C;
			double A1 = D / det, B1 = -B / det, C1 = -C / det, D1 = A / det;
			a = A1 * sumXY + B1 * sumY;
			b = C1 * sumXY + D1 * sumY;
		}

		static double LogOnePlusX( double x ) {
			if ( fabs( x ) >= 1e-10 ) {
				return log( 1 + x );
			}

			double sum = 0;
			double currentTerm = 1;
			int sign = -1;

			for ( int i = 1;; i++ ) {
				sign = -sign;
				currentTerm *= x;

				double prev = sum;

				sum += sign * currentTerm / i;

				if ( sum == prev ) {
					break;
				}

				//std::cerr << "i = " << i << ", currentTerm = " << currentTerm << ", sum = " << sum << "\n";
			}

			//std::cerr << "\n";

			return sum;
		}

		static double OneMinusExpX( double x ) {
			if ( fabs( x ) >= 1e-10 ) {
				return 1 - exp( x );
			}

			double sum = 0;
			double currentTerm = 1;

			for ( int i = 1; ; i++ ) {
				currentTerm *= x / i;

				double prev = sum;

				sum += currentTerm;

				if ( sum == prev ) break;

				//std::cerr << "i = " << i << ", currentTerm = " << currentTerm << ", sum = " << sum << "\n";
			}

			//std::cerr << "\n";

			return -sum;
		}
	};

	class Int {
	public:
		static int Parse( const string & s ) {
			try {
				size_t read = 0;
				return std::stoi( s, &read );
			}
			catch ( std::invalid_argument ) {
				throw Exception( "Invalid integer data in string '" + s + "'", FileAndLine );
			}
		}

		static string ToString( int value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<int> & x, string & delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Uint {
	public:
		static uint Parse( const string & s ) {
			try {
				size_t read = 0;
				return (uint) std::stoul( s, &read );
			}
			catch ( std::invalid_argument ) {
				throw Exception( "Invalid unsigned integer data in string '" + s + "'", FileAndLine );
			}

		}

		static string ToString( uint value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<uint> & x, string & delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Uint64 {
	public:
		static uint64_t Parse( const string & s ) {
			try {
				size_t read = 0;
				return std::stoull( s, &read );
			}
			catch ( std::invalid_argument ) {
				throw Exception( "Invalid floating point data in string '" + s + "'", FileAndLine );
			}
		}

		static string ToString( uint64_t value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<uint64_t> & x, string & delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	template <typename T>
	class Convert {
	public:
		static T Parse( const string & s ) {
			istringstream str( s );
			T result;
			str >> result;
			return result;
		}

		static string ToString( T value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<T> & x, string & delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Ulong {
	public:
		static ulong Parse( const string & s ) {
			try {
				size_t read = 0;
				return std::stoul( s, &read );
			}
			catch ( std::invalid_argument ) {
				throw Exception( "Invalid unsigned long data in string '" + s + "'", FileAndLine );
			}
		}

		static string ToString( ulong value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<ulong> & x, string & delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Double {
		/**
		*	<summary>
		*	Scans a floating point value from a supplied string.
		*	</summary>
		*	<param name="s">The string whichis isupposed to contain text
		*	representation of a floating point value.</param>
		*	<returns>
		*	</returns>
		*/
	public:
		static double Parse( const string & s ) {
			try {
				size_t read = 0;
				return std::stod( s, &read );
			}
			catch ( std::invalid_argument ) {
				throw Exception( "Invalid floating point data in string '" + s + "'", FileAndLine );
			}
		}

		static string ToString( double value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<double> & x, string & delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Bool {
	public:
		static bool Parse( const string & s ) {
			const char
				* t = "true",
				*cs = s.c_str();

			for ( ; *cs && *t; cs++, t++ ) {
				if ( tolower( *cs ) != *t ) return false;
			}

			return *cs == 0 && *t == 0;
		}

		static string ToString( double value ) {
			ostringstream str;
			str << ( value ? "true" : "false" );
			return str.str();
		}

		static string Join( const vector<bool> & x, string & delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << ( *current ? 1 : 0 );
				}
			}

			return s.str();
		}
	};

	class File {
	public:
		static bool Exists( const string & fileName ) {
			FILE * f = fopen( fileName.c_str(), "r" );
			if ( f ) {
				fclose( f );
				return true;
			}
			return false;
		}

		static void ReadStrings( std::istream & str, std::function<void( const string & s )> action ) {
			string s;
			while ( std::getline( str, s ) ) {
				action( s );
			}
		}

		static void ReadStrings( const string & fileName, std::function<void( const string & s )> action ) {
			ifstream str( fileName );
			ReadStrings( str, action );
		}
	};
}
