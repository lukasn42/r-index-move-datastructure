#include <iostream>
#include <filesystem>

#include <omp.h>

#include "utils.hpp"
#include "r_index_mds.hpp"

using namespace ri_mds;
using namespace std;

string check = string();//check occurrences on this text

bool hyb=false;
std::ofstream measurement_file;

void help(){
	cout << "ri-mds-count: number of occurrences of the input patterns." << endl << endl;

	cout << "Usage: ri-mds-count <index> <patterns>" << endl;
	//cout << "   -h           use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
	//cout << "                if the index was built with -h options enabled."<<endl;
	cout << "   -l           file to write runtime data to"<<endl;
	cout << "   <index>      index file (with extension .ri-mds)" << endl;
	cout << "   <patterns>   file in pizza&chili format containing the patterns." << endl;
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

	}/*if(s.compare("-h")==0){

		hyb=true;

	}else*/else {

		cout << "Error: unknown option " << s << endl;
		help();

	}

}


template<class idx_t>
void count(std::ifstream& in, string patterns){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    string text;
    bool c = false;

    if(check.compare(string()) != 0){

    	c = true;

		ifstream ifs1(check);
		stringstream ss;
		ss << ifs1.rdbuf();//read the file
		text = ss.str();

    }

    auto t1 = high_resolution_clock::now();

    idx_t idx;

	idx.load(in);

	if (measurement_file.is_open()) {
		measurement_file << " a=" << idx.ret_a();
	}

	auto t2 = high_resolution_clock::now();

	cout << "searching patterns ... " << endl;
	ifstream ifs(patterns);

	//read header of the pizza&chilli input file
	//header example:
	//# number=7 length=10 file=genome.fasta forbidden=\n\t
	string header;
	std::getline(ifs, header);

	ulint n = get_number_of_patterns(header);
	ulint m = get_patterns_length(header);

	uint last_perc = 0;

	ulint occ_tot=0;

	int64_t execution_time=0;

	//extract patterns from file and search them in the index
	for(ulint i=0;i<n;++i){

		uint perc = (100*i)/n;
		if(perc>last_perc){
			cout << perc << "% done ..." << endl;
			last_perc=perc;
		}

		string p = string();

		if (idx.ret_char_shift() != 0) {
			unsigned char char_shift = idx.ret_char_shift();
			for(ulint j=0;j<m;++j){
				char c;
				ifs.get(c);
				p+=c+char_shift;
			}
		} else {
			for(ulint j=0;j<m;++j){
				char c;
				ifs.get(c);
				p+=c;
			}
		}

		auto t_before = high_resolution_clock::now();

		occ_tot += idx.occ(p);

		auto t_after = high_resolution_clock::now();

		execution_time += std::chrono::duration_cast<std::chrono::microseconds>(t_after - t_before).count();
	}

	double occ_avg = (double)occ_tot / n;

	cout << endl << occ_avg << " average occurrences per pattern" << endl;

	ifs.close();

	auto t3 = high_resolution_clock::now();

	//printRSSstat();

	uint64_t load = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	cout << "Load time : " << load << " milliseconds" << endl;

	uint64_t search = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
	cout << "number of patterns n = " << n << endl;
	cout << "pattern length m = " << m << endl;
	cout << "total number of occurrences  occ_t = " << occ_tot << endl;

	cout << "Total time : " << search << " milliseconds" << endl;
	cout << "Search time : " << (double)search/n << " milliseconds/pattern (total: " << n << " patterns)" << endl;
	cout << "Search time : " << (double)search/occ_tot << " milliseconds/occurrence (total: " << occ_tot << " occurrences)" << endl;

	if (measurement_file.is_open()) {
		measurement_file << " time=" << execution_time/1000 << endl;
	}
}

int main(int argc, char** argv){

	if(argc < 3)
		help();

	int ptr = 1;

	while(ptr<argc-2)
		parse_args(argv, argc, ptr);

	string idx_file(argv[ptr]);
	string patt_file(argv[ptr+1]);

	std::ifstream in(idx_file);

	string text_file_name = idx_file.substr(idx_file.find_last_of("/\\") + 1);
	if (measurement_file.is_open()) {
		measurement_file << "RESULT type=ri_mds name=" << text_file_name;
	}

	bool fast;

	cout << "Loading r-index" << endl;

	bool is_64_bit;
	in >> is_64_bit;

	if (is_64_bit) {
		count<r_index_mds<int64_t>>(in, patt_file);
	} else {
		count<r_index_mds<int32_t>>(in, patt_file);
	}

	in.close();

}
