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

#include "Alphabet.hpp"
#include "Args.hpp"
#include "DistanceType.hpp"
#include "SimilarityMatrix.hpp"

namespace QutBio {
	class AlphabetHelper {
	public:

		/** Deciphers the arguments used to define an alphabet and distance measure.
		 * Possibilities are:
		 *
		 *	alphabet == "DNA":
		 *		the DNA alphabet will be returned.
		 *		matrix returned will be NULL, regardless of other values.
		 *		distance returned will be UngappedEdit regardless of other values.
		 *			If you want to do something more elaborate, set up a Similarity Maxtrix in a
		 *			file and use dist = Custom with no alphabet instead.
		 *
		 *	matrixId:
		 *		if this is present, we know we're using one of the BLOSUM matrices, so
		 *			the alphabet returned will be AA,
		 *			the appropriate BLOSUM matrix will be returned.
		 *
		 *			dist == "HalperinEtAl" or dist == "BlosumDistance" or dist == "UngappedEdit":
		 *				the value is recorded without any special further action.
		 *			dist not present: 
		 *				BlosumDistance will be returned.
		 *			
		 *		Otherwise:
		 *			the argument matrixFile must be present.
		 *			The custom similarity matrix is loaded from the file and will be returned.
		 *			The value of distance returned will be BlosumDistance. If you have actual
		 *				distances, negate them to produce similarities, so that BlosumDistance will
		 *				reverse the negation and poduce positive distances. It is what it is.
		 */

		static bool GetAlphabetAndMatrix(
			Args & arguments, 
			pAlphabet & alphabet, 
			pSimilarityMatrix & matrix,
			pDistanceType &distance,
			ostream & err
		) {
			bool OK = true;

			if (arguments.IsDefined("alphabet")) {
				if (!arguments.Get("alphabet", Alphabet::Values(), alphabet)) {
					err << "Unable to parse argument 'alphabet'." << endl;
					OK = false;
				}

				if (alphabet == Alphabet::DNA()) {
					matrix = 0;
					distance = DistanceType::UngappedEdit();
					return OK;
				}
			}
			
			if (arguments.IsDefined("matrixId")) {
				int matrixId = 0;

				if (!(arguments.Get("matrixId", matrixId))) {
					err << "Argument 'matrixId' not valid." << endl;
					OK = false;
				}

				vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

				bool found = false;

				for (auto x : matrices) {
					if (x == matrixId) { found = true; }
				}

				if (!found) {
					matrix = 0;
					err << "Matrix id not recognised." << endl;
					OK = false;
				}
				else {
					matrix = SimilarityMatrix::GetBlosum( matrixId );

					if ( arguments.IsDefined( "dist" ) ) {
						if (!arguments.Get("dist", DistanceType::Values(), distance)) {
							err << "Argument 'dist' is not valid." << endl;
							OK = false;
						}
					}
				}
			}

			else if (arguments.IsDefined("matrixFile")) {
				string matrixFile;
				bool isCaseSensitive = true;
				arguments.Get("matrixFile", matrixFile);
				
				if (!arguments.GetOptionalArgument("isCaseSensitive", isCaseSensitive)) {
					err << "Argument 'isCaseSensitive', if supplied, must be true or false." << endl;
					OK = false;
				}
				
				matrix = SimilarityMatrix::GetMatrix(DistanceType::Custom(), -1, matrixFile, isCaseSensitive);
				distance = DistanceType::BlosumDistance();
				alphabet = new Alphabet(matrix);
			}

			else {
				err << "Must have either 'alphabet', 'matrixId', or 'matrixFile' defined in arguments." << endl;
				OK = false;
			}

			return OK;
		}
	};
}
