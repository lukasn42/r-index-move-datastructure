#include <iostream>

#include <omp.h>

#include <filesystem>

#include "utils.hpp"
#include "r_index_mds.hpp"

using namespace ri_mds;
using namespace std;

string check = string();//check occurrences on this text
bool hyb=false;
string ofile;
std::ofstream measurement_file;

void help(){
	cout << "ri-mds-revert: reconstruct the original file." << endl << endl;

	cout << "Usage: ri-mds-locate [options] <index> <patterns>" << endl;
	//cout << "   -h           use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
	//cout << "                if the index was built with -h options enabled."<<endl;
	cout << "   -l           file to write runtime data to"<<endl;
	cout << "   -o <ofile>   output file" << endl;
	cout << "   <index>      index file (with extension .ri-mds)" << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-l")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		measurement_file.open(argv[ptr],std::filesystem::exists(argv[ptr]) ? std::ios::app : std::ios::out);

		if (!measurement_file.good()) {
			cout << "Error: cannot open measurement file" << endl;
			help();
		}

		ptr++;

	} else if(s.compare("-o")==0){

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

	if (measurement_file.is_open()) {
		measurement_file << " a=" << idx.ret_a();
	}

	string text_orig;
	text_orig.resize(idx.text_length()-1);

	auto t2 = high_resolution_clock::now();
	uint64_t load = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	cout << "Load time : " << load << " milliseconds" << endl;

	idx.revert(text_orig);

	auto t3 = high_resolution_clock::now();
	uint64_t revert = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
	cout << "Revert time : " << revert << " milliseconds" << endl;

	if (measurement_file.is_open()) {
		measurement_file << " time=" << revert << endl;
	}

	out << text_orig;
}

int main(int argc, char** argv){

	if(argc < 3)
		help();

	int ptr = 1;

	while(ptr<argc-2)
		parse_args(argv, argc, ptr);


	string idx_file(argv[ptr]);

	string text_file_name = idx_file.substr(idx_file.find_last_of("/\\") + 1);
	if (measurement_file.is_open()) {
		measurement_file << "RESULT type=ri_mds name=" << text_file_name;
	}

	std::ifstream in(idx_file);

	bool is_64_bit;
	in >> is_64_bit;

	if (is_64_bit) {
		revert<r_index_mds<uint64_t>>(in);
	} else {
		revert<r_index_mds<uint32_t>>(in);
	}

	in.close();

}
