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


#include <algorithm>
#include <vector>
#include <iostream>
#include <string.h>
#include <cstddef>

#include "db.hpp"
#include "Delegates.hpp"
#include "Exception.hpp"

using std::vector;
using std::size_t;
using std::ostream;
using std::endl;
using std::ostringstream;
using std::output_iterator_tag;

namespace QutBio {
	template <typename T>
	struct RawMatrix {
		size_t rows, cols;
		vector<vector<T>> data;

		RawMatrix( size_t rows, size_t cols ) : rows( rows ), cols( cols ), data( rows ) {
			for ( size_t r = 0; r < rows; r++ ) {
				data[r].resize( cols );
			}
		}

		virtual ~RawMatrix() {}
	};

	template <typename T>
	class FlatMatrix {
		std::size_t rows_, cols_;
		std::vector<T> data;

		friend class KmerDistanceCache;
	public:
		FlatMatrix( size_t rows = 0, size_t cols = 0 ) : rows_( rows ), cols_( cols ), data( rows * cols ) {}

		FlatMatrix( size_t rows, size_t cols, T value ) : rows_( rows ), cols_( cols ), data( rows * cols ) {
			fill( value );
		}

		FlatMatrix( FlatMatrix & other ) :
			rows_( other.rows_ ),
			cols_( other.cols_ ),
			data( other.data ) {}

		FlatMatrix & operator=( FlatMatrix & other ) {
			rows_ = other.rows_;
			cols_ = other.cols_;
			data = other.data;
			return *this;
		}

		virtual ~FlatMatrix() {}

		void resize( size_t newRows, size_t newCols ) {
			if ( rows_ * cols_ != newRows * newCols ) {
				data.resize( newRows * newCols );
			}
			rows_ = newRows;
			cols_ = newCols;
		}

		void fill( const T & value ) {
			std::fill( data.begin(), data.end(), value );
		}

		/*
		**	Gets a reference to the matrix element element at
		**	position (r,c), where 0 <= r < rows and 0 <= c < cols.
		*/

		T & operator()( size_t r, size_t c ) {
			//if ( r >= rows_ || c >= cols_ ) {
			//	throw Exception( "Index out of bounds.", FileAndLine );
			//}
			return data[r * cols_ + c];
		}

		/*
		** Serialises the object to a stream.
		*/

		friend ostream & operator << ( ostream & out, const FlatMatrix<T> & matrix ) {
			T *ptr = matrix.data;
			for ( int i = 0; i < matrix.rows_; i++ ) {
				for ( int j = 0; j < matrix.cols_; j++ ) {
					if ( j > 0 ) out << ',';
					out << *ptr++;
				}
				out << endl;
			}

			return out;
		}

		size_t rows() const { return rows_; }

		size_t cols() const { return cols_; }

		T * buffer() { return data.data(); }

		T * row( size_t r ) {
			return data.data() + r * cols_;
		}

		/*
		** Adds another matrix to this one, element by element.
		*/
		void operator += ( const FlatMatrix & other ) {
			for ( size_t i = 0; i < rows_ * cols_; i++ ) {
				data[i] += other.data[i];
			}
		}

		/*
		** Subtracts another matrix from this one, element by element.
		*/
		void operator -= ( const FlatMatrix & other ) {
			for ( int i = 0; i < rows_ * cols_; i++ ) {
				data[i] -= other.data[i];
			}
		}

		/*
		** Multiplies the matrix by another, element by element.
		*/
		void operator *= ( const FlatMatrix & other ) {
			for ( int i = 0; i < rows_ * cols_; i++ ) {
				data[i] *= other.data[i];
			}
		}

		/*
		** Divides the matrix by another, element by element.
		*/
		void operator /= ( const FlatMatrix & other ) {
			for ( int i = 0; i < rows_ * cols_; i++ ) {
				data[i] /= other.data[i];
			}
		}

		/*
		** Multiplies the matrix by another, element by element.
		*/
		void operator *= ( T scalar ) {
			for ( int i = 0; i < rows_ * cols_; i++ ) {
				data[i] *= scalar;
			}
		}

		/*
		** Divides the matrix by another, element by element.
		*/
		void operator /= ( T scalar ) {
			for ( int i = 0; i < rows_ * cols_; i++ ) {
				data[i] /= scalar;
			}
		}

		/**
		**  This is an equivalence relation. Complexity is linear in the size of the
		**  matrices. Matrices are considered equivalent if their dimensions are equal,
		**  and if corresponding elements compare equal.
		*/
		bool operator==( const FlatMatrix<T> & y ) const {
			return ( rows_ == y.rows_
				&& cols_ == y.cols_
				&& memcmp( data, y.data, rows_ * cols_ * sizeof( T ) ) == 0
				);
		}

		/**
		**  This is an equivalence relation. Complexity is linear in the size of the
		**  matrices. Matrices are considered equivalent if their dimensions are equal,
		**  and if corresponding elements compare equal.
		*/
		bool operator!=( const FlatMatrix<T> & y ) const {
			return ( rows_ != y.rows_
				|| cols_ != y.cols_
				|| memcmp( data.data(), y.data.data(), rows_ * cols_ * sizeof( T ) ) != 0
				);
		}

		// TODO: other linear algebra, as required and makes sense.
	};

	/**
	 *	MatrixView lets us use a pre-allocated array as a matrix.
	 *	Items are stored in row-major order, with item addressed computed as:
	 *		&m(r,c) = m.data + r * m.cols_ + c;
	 *
	 */
	template <typename T>
	class MatrixView {
		T * data;
		size_t rows, cols;

	public:
		MatrixView( T* data, size_t rows = 0, size_t cols = 0 ) : data( data ), rows( rows ), cols( cols ) {}

		MatrixView( T* data, size_t rows, size_t cols, T initValue ) : data( data ), rows( rows ), cols( cols ) {
			Fill( initValue );
		}

		virtual ~MatrixView() {}

		// Use this with extreme caution: 
		void Reinterpret( size_t newRows, size_t newCols ) {
			rows = newRows;
			cols = newCols;
		}

		void Fill( const T & value ) {
			std::fill( data, data + rows, value );
		}

		/*
		**	Gets a reference to the matrix element element at
		**	position (r,c), where 0 <= r < rows and 0 <= c < cols.
		*/

		T & operator()( size_t r, size_t c ) {
			//if ( r >= rows || c >= cols ) {
			//	throw Exception( "Index out of bounds.", FileAndLine );
			//}
			return data[r * cols + c];
		}

		/*
		** Serialises the object to a stream.
		*/

		friend ostream & operator << ( ostream & out, const MatrixView<T> & matrix ) {
			T *ptr = matrix.data;

			for ( int i = 0; i < matrix.rows; i++ ) {
				for ( int j = 0; j < matrix.cols; j++ ) {
					if ( j > 0 ) out << ',';
					out << *ptr++;
				}
				out << endl;
			}

			return out;
		}

		size_t Rows() const { return rows; }

		size_t Cols() const { return cols; }

		T * Buffer() { return data; }

		T * Row( size_t r ) {
			return data + r * cols;
		}

		/*
		** Adds another matrix to this one, element by element.
		*/
		void operator += ( const MatrixView& other ) {
			for ( size_t i = 0; i < rows * cols; i++ ) {
				data[i] += other.data[i];
			}
		}

		/*
		** Subtracts another matrix from this one, element by element.
		*/
		void operator -= ( const MatrixView & other ) {
			for ( int i = 0; i < rows * cols; i++ ) {
				data[i] -= other.data[i];
			}
		}

		/*
		** Multiplies the matrix by another, element by element.
		*/
		void operator *= ( const MatrixView & other ) {
			for ( int i = 0; i < rows * cols; i++ ) {
				data[i] *= other.data[i];
			}
		}

		/*
		** Divides the matrix by another, element by element.
		*/
		void operator /= ( const MatrixView & other ) {
			for ( int i = 0; i < rows * cols; i++ ) {
				data[i] /= other.data[i];
			}
		}

		/*
		** Multiplies the matrix by another, element by element.
		*/
		void operator *= ( T scalar ) {
			for ( int i = 0; i < rows * cols; i++ ) {
				data[i] *= scalar;
			}
		}

		/*
		** Divides the matrix by another, element by element.
		*/
		void operator /= ( T scalar ) {
			for ( int i = 0; i < rows * cols; i++ ) {
				data[i] /= scalar;
			}
		}

		/**
		**  This is an equivalence relation. Complexity is linear in the size of the
		**  matrices. Matrices are considered equivalent if their dimensions are equal,
		**  and if corresponding elements compare equal.
		*/
		bool operator==( const FlatMatrix<T> & y ) const {
			return ( rows == y.rows
				&& cols == y.cols
				&& memcmp( data, y.data, rows * cols * sizeof( T ) ) == 0
				);
		}

		/**
		**  This is an equivalence relation. Complexity is linear in the size of the
		**  matrices. Matrices are considered equivalent if their dimensions are equal,
		**  and if corresponding elements compare equal.
		*/
		bool operator!=( const FlatMatrix<T> & y ) const {
			return ( rows != y.rows
				|| cols != y.cols
				|| memcmp( data(), y.data(), rows * cols * sizeof( T ) ) != 0
				);
		}

		// TODO: other linear algebra, as required and makes sense.
	};

	template<typename T>
	class SubVector {
	public:
		vector<T> * base;
		size_t offset;
		size_t length;

		SubVector() : base( 0 ), offset( 0 ), length( 0 ) {}

		SubVector( vector<T> * base, size_t offset, size_t length )
			: base( base ), offset( offset ), length( length ) {
			if ( offset + length > base->size() ) throw Exception( "Bad offset and length in SubVector", FileAndLine );
		}

		SubVector( const SubVector * other ) : base( other.base ), offset( other.offset ), length( other.length ) {}

		SubVector & operator=( const SubVector * other ) { base = other.base; offset = other.offset; length = other.length; }

		T & operator[]( size_t i ) {
#if defined(SUBVECTOR_BOUNDS)
			if ( i >= length ) throw Exception( "bad index in SubVector", FileAndLine );
#endif
			return ( *base )[offset + i];
		}

		T * data() {
			return base->data() + offset;
		}

		size_t size() { return length; }
	};
}
