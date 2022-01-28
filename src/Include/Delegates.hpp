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

#include <functional>

namespace QutBio {
	using Action = std::function<void( void )>;

	template <typename T>
	using Action1 = std::function<void( T t )>;

	template <typename T1, typename T2>
	using Action2 = std::function<void( T1  t1, T2  t2 )>;

	template <typename T1, typename T2, typename T3>
	using Action3 = std::function<void( T1 t1, T2 t2, T3 t3 )>;

	template <typename T>
	using Func = std::function<T( void )>;

	template <typename T, typename Out>
	using Func1 = std::function<Out( T t )>;

	template <typename T1, typename T2, typename Out>
	using Func2 = std::function<Out( T1  t1, T2  t2 )>;

	template <typename T1, typename T2, typename T3, typename Out>
	using Func3 = std::function<Out( T1 t1, T2 t2, T3 t3 )>;
}
