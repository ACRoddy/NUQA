
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
#include <vector>
#include <algorithm>
#include <chrono>
#include <set>
#include <iterator>
#include <cstring>
#include <fcntl.h>

using namespace std;

static const int BUFFER_SIZE = 64; // Size of the line buffer.

static char placeholder = 'Z'; // An artificially high value to initially compare the next sequence to.

struct kmerCountFile
{
    FILE *file;
    int id;

    char *currentSequence;
    int count;

    // Valid is false if the file
    // is uninitialised or we have
    // reached then end.
    bool valid;
    bool target;

    char buffer[BUFFER_SIZE];

    kmerCountFile(){
        id=0;
        valid=false;
    }

    kmerCountFile(int i, std::string filename){
        id=i;
        file = fopen(filename.c_str(), "r");
        posix_fadvise(file->_fileno, 0,0, POSIX_FADV_SEQUENTIAL);

        //std::cerr << "File " << id << " opened." << std::endl;
        valid = false;
    }

    void getNextSequenceCount(){
        // Read the next line from the file into the buffer.
        if(fgets(this->buffer, BUFFER_SIZE, this->file)!=NULL) {
            // Null terminate the kmer sequence string.
            char *pos = strchr(buffer, ' ');
            (*pos) = '\0';
            this->currentSequence = buffer;

            // Parse the remainder of the string as the kmer count.
            this->count = atoi(pos+1);

            this->valid = true;
        }else{
            this->valid = false;
            //currentSequence = NULL; //this-> not used because it's a pointer/Null isn't stored in memory?
            this->count = 0;
        }
    }

};

bool any_valid(kmerCountFile **v, const int c){
    // Scan over the vector of kmer count files
    // to see if any still contain input.
    // Returns false if we are at the end of all input files.
    for(auto x=v;x<v+c;x++){
        if((*x)->valid){
            return true;
        }
    }
    return false;
}

// Determine the next sequence string to use across all the opened files.
const char * getNextMer(kmerCountFile **v, const int c){
    char *retval = &placeholder;

    for(auto x=v;x<v+c;x++){
        if(((*x)->valid) && (strcmp((*x)->currentSequence, retval) < 0)){
    		retval = (*x)->currentSequence;
        }
    }

    //cerr << "target" << retval << "\n";
    return retval;
}


int main(int argc, char *argv[]) {
    std::setlocale(LC_ALL, "C");
    std::ios_base::sync_with_stdio(false);

    int numberOfInputFiles = argc-1;
    kmerCountFile *inputFiles[numberOfInputFiles];

    int j=0;

    // Initialise the input file data structures.
    for(int x=1;x<argc;x++){
        inputFiles[x-1] = new kmerCountFile(j++, std::string(argv[x]));
        inputFiles[x-1]->getNextSequenceCount();
    }

    const char* targetMer; // The current target sequence being counted.

    // Iterate over all input files until they are all exhausted.
    while(any_valid(inputFiles, numberOfInputFiles)) {

        targetMer = getNextMer(inputFiles, numberOfInputFiles);
        //printf("%s ", targetMer); // Comment out to omit targetMer from output, if not required. Does have a speed advantage.

        for(int i=0;i<numberOfInputFiles;i++){
        	//cerr << "mer" << inputFiles[i]->currentSequence << targetMer << "\n";
            if(!(inputFiles[i]->valid) || strcoll(inputFiles[i]->currentSequence, targetMer)!=0) {
           		// targetMer does not appear in this file.
                printf("0 ");
                inputFiles[i]->target = false;
                continue;
            }else{
                printf("%d ", inputFiles[i]->count);
                inputFiles[i]->target = true;
            }
        }

        for(int i=0;i<numberOfInputFiles;i++){
			if((inputFiles[i]->valid) && inputFiles[i]->target==true ) {
				inputFiles[i]->getNextSequenceCount();
			}
		}


        printf("\n");


    }


    // Cleanup
    for(int i=0;i<numberOfInputFiles;i++){
        fclose(inputFiles[i]->file);
        //std::cerr << "File " << inputFiles[i]->id << " closed." << std::endl;
        delete inputFiles[i];
    }
    std::cerr << "Merging of counts vectors complete." << std::endl;
    return 0;
}





