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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   db.hpp
 * Author: Lawrence
 *
 * Created on 22 September 2017, 12:40 AM
 */

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201103L
#endif

#include <cstdio>

#ifdef TRON
#undef TRACE
#undef PR
#undef DB
#define TRACE {fprintf(stderr,"%s:%d\n",__FILE__, __LINE__); fflush(stderr);}
#define PR(symbol,format) {fprintf(stderr,"%s:" #format "\n", #symbol, symbol); fflush(stderr);}
#define DB(stuff) stuff
#else
#undef TRACE
#undef PR
#undef DB
#define TRACE 
#define PR(symbol,format) 
#define DB(stuff) 
#endif


