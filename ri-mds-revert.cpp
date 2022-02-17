#include <iostream>

#include "internal/r_index_mds.hpp"

#include "internal/utils.hpp"

using namespace ri_mds;
using namespace std;

string check = string();//check occurrences on this text
bool hyb=false;
string ofile;

void help(){
	cout << "ri-mds-revert: reconstruct the original file." << endl << endl;

	cout << "Usage: ri-mds-locate [options] <index> <patterns>" << endl;
	//cout << "   -h           use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
	//cout << "                if the index was built with -h options enabled."<<endl;
	cout << "   -o <ofile>   output file" << endl;
	cout << "   <index>      index file (with extension .ri-mds)" << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-o")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		ofile = string(argv[ptr]);
		ptr++;

	}/*else if(s.compare("-h")==0){

		hyb=true;

	}*/else{

		cout << "Error: unknown option " << s << endl;
		help();

	}

}


template<class idx_t>
void revert(std::ifstream& in){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;


    auto t1 = high_resolution_clock::now();

	ofstream out;

    if(ofile.compare(string())!=0){

    	out = ofstream(ofile);

    }

    idx_t idx;
	idx.load(in);
	string text_orig;
	text_orig.resize(idx.text_length()-1);

	auto t2 = high_resolution_clock::now();
	uint64_t load = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	cout << "Load time : " << load << " milliseconds" << endl;

	idx.revert(text_orig);

	auto t3 = high_resolution_clock::now();
	uint64_t revert = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
	cout << "Revert time : " << revert << " milliseconds" << endl;

	out << text_orig;
}

int main(int argc, char** argv){

	if(argc < 3)
		help();

	int ptr = 1;

	while(ptr<argc-2)
		parse_args(argv, argc, ptr);


	string idx_file(argv[ptr]);

	std::ifstream in(idx_file);

	bool fast;

	cout << "Loading r-index" << endl;

	if(hyb){

		//locate<r_index_mds<sparse_hyb_vector,rle_string_hyb> >(in, patt_file);

	}else{

		revert<r_index_mds<> >(in);

	}

	in.close();

}
