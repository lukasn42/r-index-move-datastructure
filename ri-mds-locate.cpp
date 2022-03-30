#include <iostream>
#include <filesystem>

#include <omp.h>

#include "utils.hpp"
#include "r_index_mds.hpp"

using namespace ri_mds;
using namespace std;

string check = string();//check occurrences on this text
bool hyb=false;
string ofile;
std::ofstream measurement_file;

void help(){
	cout << "ri-mds-locate: locate all occurrences of the input patterns." << endl << endl;

	cout << "Usage: ri-mds-locate [options] <index> <patterns>" << endl;
	cout << "   -c <text>    check correctness of each pattern occurrence on this text file (must be the same indexed)" << endl;
	//cout << "   -h           use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. -h is required "<<endl;
	//cout << "                if the index was built with -h options enabled."<<endl;
	cout << "   -l           file to write runtime data to"<<endl;
	cout << "   -o <ofile>   write pattern occurrences to this file (ASCII)" << endl;
	cout << "   <index>      index file (with extension .ri-mds)" << endl;
	cout << "   <patterns>   file in pizza&chili format containing the patterns." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-c")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -c option." << endl;
			help();
		}

		check = string(argv[ptr]);
		ptr++;

	} else if(s.compare("-l")==0){

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

	}else if(s.compare("-o")==0){

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
void locate(std::ifstream& in, string patterns){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    string text;
    bool c = false;

    ofstream out;

    if(ofile.compare(string())!=0){

    	out = ofstream(ofile);

    }

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

	int64_t execution_time=0;

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

	//extract patterns from file and search them in the index
	for(ulint i=0;i<n;++i){

		uint perc = (100*i)/n;
		if(perc>last_perc){
			cout << perc << "% done ..." << endl;
			last_perc=perc;
		}

		string p = string();

		if (idx.ret_chars_mapped()) {
			for(ulint j=0;j<m;++j){
				char c;
				ifs.get(c);
				p+=idx.ret_map_char(c);
			}
		} else {
			for(ulint j=0;j<m;++j){
				char c;
				ifs.get(c);
				p+=c;
			}
		}

		//cout << "locating " << idx.occ(p) << " occurrences of "<< p << " ... " << flush;

		auto t_before = high_resolution_clock::now();

		auto OCC = idx.locate_all(p);	//occurrences

		auto t_after = high_resolution_clock::now();

		execution_time += std::chrono::duration_cast<std::chrono::microseconds>(t_after - t_before).count();

		if(ofile.compare(string())!=0){

			sort(OCC.begin(),OCC.end());

			for(auto x:OCC) out << (int)x << endl;

		}

		occ_tot += OCC.size();

		if(c){//check occurrences

			//remove duplicates, if any (there shouldn't be!)
			sort(OCC.begin(),OCC.end());
			auto it = unique(OCC.begin(),OCC.end());
			OCC.resize(std::distance(OCC.begin(),it));

			if(OCC.size() != idx.occ(p)){

				cout << "Error: wrong number of located occurrences: " << OCC.size() << "/" << idx.occ(p) << endl;
				exit(0);

			}

			if (idx.ret_chars_mapped()) {
				for (ulint i=0; i<p.size(); i++) {
					p[i] = idx.ret_map_char_rev(p[i]);
				}
			}

			for(auto o:OCC){

				//for(auto c : p) cout << int(c) << " ";

				if(text.substr(o,p.size()).compare(p) != 0){

					cout << "Error: wrong occurrence: " << o << " ("  << occ_tot << " occurrences"  << ") "<< endl;
					for(auto c : text.substr(o,p.size())) cout << c << " ";
					cout << "  /  ";
					for(auto c : p) cout << c << " ";
					cout << endl;

					break;

					//exit(0);

				}

			}

		}

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
		locate<r_index_mds<int64_t>>(in, patt_file);
	} else {
		locate<r_index_mds<int32_t>>(in, patt_file);
	}

	in.close();

}
