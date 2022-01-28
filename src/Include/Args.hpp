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
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "Exception.hpp"
#include "EnumBase.hpp"
#include "SimilarityMatrix.hpp"
#include "String.hpp"
#include "Util.hpp"

using namespace std;

namespace QutBio {

	/// <summary> Main class contains the program entry point: Run().
	/// </summary>
	class Args {
		friend std::ostream & operator<< ( std::ostream & str, Args args );

	private:
		unordered_map<string, vector<string>> arguments;

	public:

		Args( int argc, char ** argv ) {
			ParseArgs( argc, argv, arguments );
		}

		Args( size_t argc, const char ** argv ) {
			ParseArgs( (int) argc, (char**) argv, arguments );
		}

		bool Contains( string & key ) {
			string key_ = String::ToLowerCase( key );
			return arguments.find( key_ ) != arguments.end();
		}

		bool IsDefined( const char * key ) {
			string key_ = String::ToLowerCase( key );
			return arguments.find( key_ ) != arguments.end();
		}

		/// <summary> Safely gets the value of a single-valued argument. If the 
		///		named argument is not present, returns undefined. Otherwise, returns 
		///		the first element of the named list of arguments.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, string & result ) {
			string key_( String::ToLowerCase( key ) );

			if ( arguments.find( key_ ) != arguments.end() ) {
				vector<string> & values( arguments[key_] );

				result = values.size() > 0 ? values[0] : "";
				return true;
			}

			return false;
		}

		/// <summary> Safely gets the value of a single-valued argument. If the 
		///		named argument is not present, returns false. Otherwise, returns 
		///		true and copies back by reference the first element of the 
		///		named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <param name="result"></param>
		/// <returns></returns>

		bool Get( const char * key, string & result ) {
			return Get( string( key ), result );
		}

		/// <summary> Replaces the value of an argument.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		void  Set( const string & key, const string & result ) {
			string key_( String::ToLowerCase(key) );
			auto & results = arguments[key_];
			results.push_back(result);
		}

		/// <summary> Safely gets the value of a vector-valued argument. If the 
		///		named argument is not present, returns false. Otherwise, returns 
		///		true and passes back (by reference) the named list of values.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( string & key, vector<string> & result ) {
			string key_( String::ToLowerCase( key ) );

			if ( arguments.find( key_ ) == arguments.end() ) return false;

			result.clear();

			for ( auto s : arguments[key_] ) {
				result.push_back( s );
			}

			return true;
		}

		/// <summary> Safely gets the value of a single-valued argument. If the 
		///		named argument is not present, returns undefined. Otherwise, returns 
		///		the first element of the named list of arguments.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const char * key, vector<string> & result ) {
			string k( key );
			return Get( k, result );
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a double. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, double & result ) {
			string key_( String::ToLowerCase( key ) );

			if ( arguments.find( key_ ) != arguments.end() ) {
				vector<string> & values( arguments[key_] );

				if ( values.size() > 0 ) {
					result = Double::Parse( values[0] );
					return true;
				}
			}

			return false;
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a double. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const char * key, double & result ) {
			return Get( string( key ), result );
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as an int. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, int & result ) {
			string key_( String::ToLowerCase( key ) );

			if ( arguments.find( key_ ) != arguments.end() ) {
				vector<string> & values( arguments[key_] );

				if ( values.size() > 0 ) {
					result = Int::Parse( values[0] );
					return true;
				}
			}

			return false;
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as an int. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, uint & result ) {
			string key_( String::ToLowerCase( key ) );

			if ( arguments.find( key_ ) != arguments.end() ) {
				vector<string> & values( arguments[key_] );

				if ( values.size() > 0 ) {
					result = Uint::Parse( values[0] );
					return true;
				}
			}

			return false;
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as an int. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, Distance & result ) {
			string key_( String::ToLowerCase( key ) );

			if ( arguments.find( key_ ) != arguments.end() ) {
				vector<string> & values( arguments[key_] );

				if ( values.size() > 0 ) {
					result = Int::Parse( values[0] );
					return true;
				}
			}

			return false;
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a size_t. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, size_t & result ) {
			string key_( String::ToLowerCase( key ) );

			if ( arguments.find( key_ ) != arguments.end() ) {
				vector<string> & values( arguments[key_] );

				if ( values.size() > 0 ) {
					result = Uint::Parse( values[0] );
					return true;
				}
			}

			return false;
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a double. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const char * key, int & result ) {
			return Get( string( key ), result );
		}


		/// <summary> Safely gets the value of a single-valued argument and parses it as a boolean
		///		value. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, bool & result ) {
			string key_( String::ToLowerCase( key ) );

			if ( arguments.find( key_ ) != arguments.end() ) {
				vector<string> & values( arguments[key_] );

				if ( values.size() > 0 ) {
					result = Bool::Parse( values[0] );
					return true;
				}
			}

			return false;
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a boolean. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const char * key, bool & result ) {
			return Get( string( key ), result );
		}


		/// <summary> Safely gets the value of a single-valued argument and parses it as a double. 
		///		If the named argument is not present, returns false. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key ) {
			string keyLC( String::ToLowerCase( key ) );

			if ( arguments.find( keyLC ) != arguments.end() ) {
				auto values = arguments[keyLC];
				return values.size() > 0 ? String::ToLowerCase( values[0] ) == "true" : true;
			}
			else {
				return false;
			}
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a double. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const char * key ) {
			return Get( string( key ) );
		}

		void Show() {
			cout << ( *this );
		}

		template<typename T>
		bool Get( const string & key, vector<EnumBase *> values, T *& result ) {
			string keyLC( String::ToLowerCase( key ) );

			if ( arguments.find( keyLC ) != arguments.end() ) {
				auto vals = arguments[keyLC];

				if ( vals.size() > 0 ) {
					result = EnumBase::Parse<T>( vals[0], values );
					return true;
				}
			}

			return false;
		}

		template<typename T>
		bool Get( const char *key, unordered_set<T> & values ) {
			vector<string> vals;

			if ( Get( key, vals ) ) {
				for ( auto s : vals ) {
					istringstream str( s );
					T value;
					str >> value;
					values.insert( value );
				}

				return true;
			}
			else {
				return false;
			}
		}

		template <typename T>
		bool GetOptionalArgument( const char * name, T & arg ) {
			return IsDefined( name ) ? Get( name, arg ) : true;
		}

		template <typename T>
		bool GetOptionalArgument( const char * name, T & arg, Action showHelp ) {
			if ( IsDefined( name ) ) {
				if ( !Get( name, arg ) ) {
					cerr << "Unable to parse " << name << "." << endl;
					showHelp();
					return false;
				}
			}

			return true;
		}

		template <typename T>
		bool Get( const char * name, vector<T> & values ) {
			vector<string> vals;

			if ( Get( name, vals ) ) {
				for ( auto s : vals ) {
					istringstream str( s );
					T value;
					str >> value;
					values.push_back( value );
				}

				return true;
			}
			else {
				return false;
			}
		}


		/**
		*	Gets a Similarity matrix based on an argument list, which will need to have
		*/

		bool Get( SimilarityMatrix *&matrix, string & error ) {
			int matrixId = 0;

			if ( IsDefined( "matrixId" ) ) {
				if ( !Get( "matrixId", matrixId ) ) {
					error = ProgName() + ": error - argument 'matrixId' not valid.";
					return false;
				}

				vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

				bool found = false;

				for ( auto x : matrices ) {
					if ( x == matrixId ) {
						found = true;
					}
				}

				if ( !found ) {
					error = ProgName() + ": error - matrix id not recognised.";
					return false;
				}
			}

			string matrixFile;
			DistanceType * distanceType = DistanceType::BlosumDistance();

			if ( IsDefined( "matrixFile" ) ) {
				Get( "matrixFile", matrixFile );
				distanceType = DistanceType::Custom();
				matrixId = -1;
			}

			bool isCaseSensitive = false;

			if ( IsDefined( "isCaseSensitive" ) ) {
				if ( !Get( "isCaseSensitive", isCaseSensitive ) ) {
					error = ProgName() + ": Invalid data for argument 'isCaseSensitive'.";
					return false;
				}
			}

			matrix = SimilarityMatrix::GetMatrix( distanceType, matrixId, matrixFile, isCaseSensitive );

			if ( !matrix ) {
				error = ProgName() + ": Unable to construct similarity matrix.\n"
					"For help, run: " + ProgName() + " --help\n";
				return false;
			}

			return true;
		}


		string ProgName() {
			string progName;
			Get( "", progName );
			return progName;
		}

		/// <summary> Parses the supplied command line arguments and returns
		///		a unordered_map containing the results. For each named argument 
		///		(specified by a leading '-') a list of values is returned.
		///	<para>
		///		Arguments which appear before any named arguments are returned 
		///		in the entry with key = string.Empty.
		/// </para>
		/// </summary>
		/// <param name="args">
		///		A list of string containing the arguments.
		/// </param>
		/// <returns>
		///		A unordered_map with the arguments
		/// </returns>

	private:

		static void ParseArgs( int argc, char ** argv, unordered_map<string, vector<string>> & arguments ) {
			string currentKey;
			vector<string> * currentValues = &arguments[currentKey];

			for ( int i = 0; i < argc; i++ ) {
				string arg( argv[i] );

				if ( arg.size() >= 2 && arg[0] == '-' && arg[1] == '-' ) {
					currentKey = arg;
					currentKey.erase( currentKey.begin() );
					currentKey.erase( currentKey.begin() );
					currentValues = &arguments[String::ToLowerCase( currentKey )];
				}
				else {
					currentValues->push_back( arg );
				}
			}
		}
	};

	/// <summary> Echoes the command line arguments, parsed into lists.
	/// </summary>
	/// <param name="arguments"></param>

	ostream & operator<< ( ostream & str, Args args ) {
		bool deja( false );

		for ( auto arg : args.arguments ) {
			if ( deja ) str << " \\" << endl;

			deja = true;

			str << "--" << arg.first;

			for ( auto value : arg.second ) {
				str << " " << value;
			}
		}

		str << endl;

		return str;
	}

}
