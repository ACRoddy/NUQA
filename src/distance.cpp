/*  This file is part of NUQA.

    NUQA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    NUQA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with NUQA.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <fcntl.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


static const char *optString = "n:l:d:";
static const int BUFFER_SIZE = 64;// Size of the line buffer.


struct globalArgs_t {
    const char *normFileName;   /* -n option */
    FILE *normFile;
    const char *labFileName;    /* -l option */
    istream *labFile;
    const char *distance;       /* -d option */
    char **inputFiles;          /* input files */
    int numInputFiles;          /* # of input files */
} globalArgs;

//Find Jensen-Shannon Divergence values for each row vector passed
ArrayXXd JSDeigen(ArrayXd v, int N){
	static ArrayXXd P(N,N);
	static ArrayXXd M(N,N);
		
	for (size_t i = 0, size = P.size(); i < size; i++)
      	{
		int j = i%N;
		M(i) = (v[j] + v[i/N])/2.0;
		P(i) = v[j]*log2(v[j]/M(i));
	}

	P = P.isNaN().select(0,P);
	return (P+P.transpose())/2.0;

}


//Find Hellinger Distance values for each row vector passed
ArrayXXd HellingerEigen(ArrayXd v, int N){

	//make matrix out of vector*len, transpose, sqrt, subtract
	static ArrayXXd P(N,N);
	for (int i=0; i<N; i++) {
		P.row(i) = v;
	}
	
	return (P.sqrt() - (P.transpose()).sqrt()).pow(2);
}

//Print Eigen Array in a format that can be read in by Phylip
void printDM(ArrayXXd R, int N, vector<string> L){
	cout << '\t' << N;

	for (int i=0; i<N*N; i++){
		if(i%N==0){
			
			cout << '\n' << L[i/N];
		}else{ cout << '\t';}
		
		cout << R(i);
	}
}

int main(int argc, char *argv[]) {


	/*arguments
	-n --norm       input file containing total counts across each ffp
	-o --output     name of output file [use > instead?]
	-d --distance   metric to be used
	*/


	//set up
	std::setlocale(LC_ALL, "C");
	std::ios_base::sync_with_stdio(false);

	//init globalArgs
	int opt=0;
	globalArgs.normFileName = NULL;
	globalArgs.normFile = NULL;
	globalArgs.labFileName = NULL;
	globalArgs.labFile = NULL;
	globalArgs.distance = NULL;
	globalArgs.inputFiles = NULL;
	globalArgs.numInputFiles = 0;

	//Read parameters
	opt = getopt(argc, argv, optString);
	while (opt != -1){
		switch (opt){
			case 'n':
				globalArgs.normFileName = optarg;
				break;

			case 'l':
				globalArgs.labFileName = optarg;
				break;

			case 'd':
				globalArgs.distance = optarg;
				break;

			default:
				break;
		}
		opt = getopt(argc, argv, optString);
	}
    	globalArgs.inputFiles = argv + optind;
    	globalArgs.numInputFiles = argc - optind;

	//open and load total count values to normalise counts
	vector<double> norm;
    	globalArgs.normFile = fopen(globalArgs.normFileName, "r");
    	char buffer[BUFFER_SIZE];

    	while (fgets(buffer, BUFFER_SIZE, globalArgs.normFile)!=NULL){
    	norm.push_back(atof(buffer));
    	}
	fclose(globalArgs.normFile);
	int Nsize = norm.size();
	Map<ArrayXd> normVals(norm.data(),Nsize);	

	/********************************USING Eigen ***************************************/

	//open and initialise counts file
	ArrayXXd temp = ArrayXXd::Zero(Nsize, Nsize); //matrix to hold current sum of values
	ArrayXXd DistM = ArrayXXd::Zero(Nsize, Nsize); //matrix to hold final distances

	FILE *input = fopen(*globalArgs.inputFiles, "r");
	posix_fadvise(input->_fileno, 0,0, POSIX_FADV_SEQUENTIAL); //Read in counts file

	char buff[BUFFER_SIZE];
	ArrayXd kmerCount(Nsize);
	

	ArrayXXd (*DM)(ArrayXd, int);//pointer to function depending on parameter -d
	switch (*globalArgs.distance){
		case 'j':
			DM = &JSDeigen;
			break;

		case 'h':
			DM = &HellingerEigen;
			break;

		default:
			DM = &JSDeigen;
			break;
	}

	//read file row by row and find DM value
	while(fgets(buff, BUFFER_SIZE, input)!=NULL) {
		char *counts = strtok(buff, " ");
		int i=0;

		while (counts != NULL){
			kmerCount(i) = atof(counts);
			counts = strtok(NULL, " \n");
			i++;
		}

		
		temp += DM(kmerCount/normVals,Nsize);
		
	}
	//cerr << temp;
	//cerr << Nsize;

	//increment distance metric after checking cache (if?)

	//output distance metric

	switch (*globalArgs.distance){
		case 'j':
			DistM = temp;
			break;

		case 'h':
			DistM  = temp.sqrt()/sqrt(2);
			break;

		default:
			DistM = temp;
			break;
	}
	
	//Get labels based on file names
	vector<string> Labs;
	if(globalArgs.labFileName != NULL){
		//globalArgs.labFile = fopen(globalArgs.labFileName, "r");
		globalArgs.labFile = new ifstream(globalArgs.labFileName);
	    	char buffer2[BUFFER_SIZE];
		string line;
		string line1;
		while (getline(*globalArgs.labFile, line)){
			
			line1 = line.substr(0,8);
			line1.insert(line1.end(), 10-line1.size(), ' ');	
			Labs.push_back(line1);
		}
	}
	
	//Print final distance matrix to cout
	printDM(DistM,Nsize,Labs);

	return 0;
	

	/***********************************************************************************/


}
