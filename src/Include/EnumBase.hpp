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
#include <sstream>

#include "String.hpp"
#include "Array.hpp"
#include "Util.hpp"

using std::string;
using std::cerr;

using namespace QutBio;

// TODO: This clumsy attempt at a portable enum type can be used without too much change in JavaScript ports of this code. See if it can be made better.

namespace QutBio {

	/// <summary> Enumeration of the available 
	/// </summary>

	class EnumBase {
		friend ostream & operator<<( ostream &, const EnumBase & );

	private:
		std::string name;
		int value;

	protected:
		EnumBase( std::string name, int value ) : name( String::ToLowerCase( name ) ), value( value ) {}

	public:
		const string & ToString() {
			return name;
		}

		template<typename T>
		static T * Parse( const std::string & s, std::vector<EnumBase *> & values ) {
			string t = String::ToLowerCase( s );

			for ( size_t i = 0; i < values.size(); i++ ) {
				// cerr << "Checking " << values[i]->name << "\n";

				if ( t == values[i]->name ) {
					return (T*) values[i];
				}
			}

			ostringstream str;
			str << "Format Exception. Enumerated value '" << s << "' not recognised.";
			( cerr << str.str() << "\n" ).flush();
			throw new Exception( str.str(), __FILE__, __LINE__ );
		}

		bool operator == ( const EnumBase & other ) const {
			return this == &other ? true : ( this->value == other.value && this->name == other.name );
		}

		bool operator != ( const EnumBase & other ) const {
			return !( this == &other ? true : ( this->value == other.value && this->name == other.name ) );
		}

		int Value() { return value; }

		const string & Name() { return name; }
	};

	ostream & operator << ( ostream & out, const EnumBase & item ) {
		out << item.name;
		return out;
	}
}
