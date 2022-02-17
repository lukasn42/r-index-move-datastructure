#include <iostream>

#include "internal/r_index_mds.hpp"
#include "internal/utils.hpp"

using namespace ri_mds;
using namespace std;

bool hyb=false;

void help(){
	cout << "ri-mds-space: breakdown of index space usage" << endl;
	cout << "Usage:       ri-mds-space <index>" << endl;
	//cout << "   -h        use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
	//cout << "             if the index was built with -h options enabled."<<endl;
	cout << "   <index>   index file (with extension .ri-mds)" << endl;
	exit(0);
}


void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	/*if(s.compare("-h")==0){

		hyb=true;

	}else*/{

		cout << "Error: unknown option " << s << endl;
		help();

	}

}

int main(int argc, char** argv){

	int ptr = 1;

	if(argc<2) help();

	while(ptr<argc-1)
		parse_args(argv, argc, ptr);

	if(hyb){

		r_index_mds<sparse_hyb_vector,rle_string_hyb> idx;
		idx.load_from_file(argv[ptr]);
		auto space = idx.print_space();
		cout << "\nTOT space: " << space << " Bytes" <<endl;

	}else{

		r_index_mds<> idx;
		idx.load_from_file(argv[ptr]);
		auto space = idx.print_space();
		cout << "\nTOT space: " << space << " Bytes" <<endl;

	}




}
