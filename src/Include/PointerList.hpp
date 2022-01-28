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

#include <vector>
#include <iterator>

#include "Assert.hpp"
#include "Delegates.hpp"

using namespace QutBio;

namespace QutBio
{
/**
	<summary>
		A lookup table mapping some key of type K to _pointer_to_ values of type V, which owns the stored pointers.
		The values are dynamically allocated when a value is added by a factory function which is responsible for
		correct initialisation.
		When the collection is deleted the pointers are destroyed.
	<summary>
		*/

template <typename V>
class PointerList
{
  private:
	vector<V *> vec;

  public:
	typedef V *Pointer;

	PointerList() {}

	~PointerList()
	{
		// This is just too darn slow!
		//	for ( int i = 0; i < vec.size(); i++ ) {
		//		delete vec[i];
		//	}
	}

	PointerList(const PointerList &other) = delete;

	PointerList &operator=(const PointerList &other) = delete;

	const vector<V *> &Items(void) { return vec; }

	template<typename T>
	T **AsPointersTo() { return (T**)(vec.data()); }

	/**
			Uses the supplied factory function to create a pointer to new dynamically
			allocated object. The pointer is appended to the current list. Later, the 
			pointer will be disposed via the delete operator, so
			the dynamic allocation needs to be essentially a single call to new,
			followed by initialisation.

			DO NOT USE ADD to append shallow copies of objects, or you will get a
			"dual call to free" event!.
		 */
	void Add(Func<V *> factory)
	{
		V *value = factory();
		vec.push_back(value);
	}

	size_t Length() const
	{
		return vec.size();
	}

	V *&operator[](size_t index) 
	{
		return vec[index];
	}

	void ForEach(Action1<V *> action) const
	{
		for (auto v : vec)
		{
			action(v);
		}
	}

	typename vector<Pointer>::iterator begin()
	{
		return vec.begin();
	}

	typename vector<Pointer>::iterator end()
	{
		return vec.end();
	}
};
} // namespace QutBio
