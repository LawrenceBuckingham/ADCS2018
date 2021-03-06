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

/**
**	<summary>
**		Parses a tabulated discrete PDF expressed as a series <(x_i,p_i)|i=0..N>
**		from a file such as that produced by GetDnaHammingDistribution or
**		GetKmerTheoreticalDistanceDistributions.
**		Computes the inverse CDF of the model at a range of designated input
**		probabilities. Since this is a discrete model, for input probability P,
**		the value returned will be x_i* such that (sum_(i=0)^(i*) x_i) <= P
**		and (sum_(i=0)^(i*+1) x_i) > P.
**		Special cases:
**			P <= x_0 --> x_0 - 1;
**			P > x_N --> x_N.
**	</summary>
*/

#include <cstdlib>
#include <mutex>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Args.hpp"
#include "Delegates.hpp"
#include "Histogram.hpp"
#include "DiscreteDistribution.hpp"

using namespace QutBio;
using namespace std;

Args * args;

void Run(void) {
	string inFile;
	vector<double> pValues;
	bool ok = true;

	if (args->IsDefined("help")) {
		vector<string> text{
			"GetCdfInverse: Reports (to standard output) a list of inverse CDF values from a histogram file.",
			"Arguments:"
			"--help      : Gets this text.",
			"--inFile    : Required. The path to a file which contains a histogram such as that produced by\n"
			"              GetKmerTheoreticalDistanceDistributions.",
			"--pValues   : Required. A list of (floating point) probability thresholds for which the inverse\n"
			"              CDF is wanted.",
			"--numThreads: Optional; default value = 7. The number of OpenMP threads to use in parallel regions.",
		};

		for (auto s : text) {
			cerr << s << "\n\n";
		}
	}

	if (!args->Get("inFile", inFile)) {
		cerr << "Command line argument '--inFile' is required." << endl;
		ok = false;
	}

	if (!args->Get("pValues", pValues)) {
		cerr << "Command line argument '--pValues' is required." << endl;
		ok = false;
	}

	if ( ! ok ) {
		cerr << "Command line arguments not valid.\nFor help: GetCdfInverse --help\n\n";
		return;
	}

	ifstream inStream( inFile );

	if (inStream.fail()) {
		cerr << "Unable to read from '" << inFile << "'.\n";
		return;
	}

	Histogram<double> hist;

	hist.ParseRows( inStream, '\t', [](string & s){ 
		return atof(s.c_str());
	});
	
	DiscreteDistribution dist;
	dist.SetPmf(hist);

	cout << "p\tx\n";

	for ( double p: pValues ) {
		cout << p << "\t" << dist.InverseCdf(p) << "\n";
	}
}

int main(int argc, char** argv) {
	try {
		args = new Args(argc, argv);
		Run();
	}
	catch (Exception ex) {
		cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")\n";
	}
	catch (runtime_error & err) {
		cerr << "Unhandled exception:\n" << err.what() << endl;
	}
	return 0;
}

