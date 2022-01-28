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
#include <stdexcept>

using string = std::string;
using runtime_error = std::runtime_error;

namespace QutBio {

	class Exception : public std::runtime_error {
	private:
		string file;
		int line;

	public:
		Exception( string message_, string file, int line )
			: runtime_error( message_ ), file( file ), line( line ) {}

		const string & File() { return file; };

		int Line() { return line; }
	};

	class KeyNotFoundException : public Exception {
	private:
		string key;

	public:
		KeyNotFoundException( string message_, string file, int line )
			: Exception( message_, file, line ), key( "Key not found" ) {}

		KeyNotFoundException( string message_, string key_, string file, int line )
			: Exception( message_, file, line ), key( key_ ) {}

		string & Key() { return key; }
	};

	class NotImplementedException : public Exception {
	public:
		NotImplementedException( string file, int line )
			: Exception( "Not implemented.", file, line ) {}
	};
}

#define FileAndLine __FILE__, __LINE__

