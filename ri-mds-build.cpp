#include <iostream>
#include <filesystem>

#include <omp.h>

extern "C" {
	#include <malloc_count.h>
}

int64_t time_diff_ms(std::chrono::_V2::system_clock::time_point t1, std::chrono::_V2::system_clock::time_point t2) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
}

#include "utils.hpp"
#include "r_index_mds.hpp"

using namespace ri_mds;
using namespace std;

string out_basename=string();
string input_file=string();
int sa_rate = 512;
int p = omp_get_max_threads();
int a = 8;
int v = 3;
std::ofstream measurement_file;
ulint T = 0;//Build fast index with SA rate = T
bool fast = false;//build fast index
bool hyb = false; //use hybrid bitvectors instead of sd_vectors?

void help(){
	cout << "ri-mds-build: builds the r-index. Extension .ri-mds is automatically added to output index file" << endl << endl;
	cout << "Usage: ri-mds-build [options] <input_file_name>" << endl;
	cout << "   -o <basename>        use 'basename' as prefix for all index files. Default: basename is the specified input_file_name"<<endl;
	//cout << "   -h                   use hybrid bitvectors instead of elias-fano in both RLBWT and predecessor structures. Important: "<<endl;
	//cout << "                        if the index is built with -h, need to specify -h also when locating and counting (ri-mds-count/ri-mds-locte). "<<endl;
	//cout << "   -fast                build fast index (O(occ)-time locate, O(r log(n/r)) words of space). By default, "<<endl;
	//cout << "                        small index is built (O(occ*log(n/r))-time locate, O(r) words of space)"<<endl;
	//cout << "   -sa_rate <T>         T>0. if used, build the fast index (see option -fast) storing T SA samples before and after each"<<endl;
	//cout << "                        BWT equal-letter run. O(r*T) words of space, O(occ(log(n/r)/T) + log(n/r))-time locate. "<<endl;
	cout << "   -p                   number of threads to use (1 for v=1/2, 1<=p for v=3, 2<=p for v=4) (default: all threads)"<<endl;
	cout << "   -a                   balancing parameter, restricts size to O(r*(1+1/(a-1))^2) (default: 8)"<<endl;
	cout << "   -v                   balancing algorithm, (1/2/3/4) (default: 3)"<<endl;
	cout << "   -l                   file to write runtime and memory usage data to"<<endl;
	cout << "   <input_file_name>    input text file." << endl;
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

		out_basename = string(argv[ptr]);
		ptr++;

	}else if(s.compare("-p")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		p = atoi(argv[ptr]);
		ptr++;

	}else if(s.compare("-a")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		a = atoi(argv[ptr]);
		ptr++;

	}else if(s.compare("-v")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		v = atoi(argv[ptr]);
		ptr++;

	}else if(s.compare("-l")==0){

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

	}/*else if(s.compare("-h")==0){

		hyb=true;

	}else if(s.compare("-fast")==0){

		fast=true;

	}else if(s.compare("-T")==0){

		T = atoi(argv[ptr]);

		if(T<=0){
			cout << "Error: parameter T must be T>0" << endl;
			help();
		}

		ptr++;
		fast=true;

	}*/else{
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

	if (!(
		1 <= p && p <= omp_get_max_threads() &&
		((v == 1 && p == 1) || (v == 2 && p == 1) || v == 3 || (v == 4 && 2 <= p))
	)) {
		cout << "Error: incompatible combination of threads and balancing algorithm" << endl;
		help();
	}
}

int main(int argc, char** argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

	//parse options

    out_basename=string();
    input_file=string();
	int ptr = 1;

	if(argc<2) help();

	while(ptr<argc-1)
		parse_args(argv, argc, ptr);

	input_file = string(argv[ptr]);

	if(out_basename.compare("")==0)
		out_basename = string(input_file);

	string idx_file = out_basename;
	idx_file.append(".ri-mds");


	cout << "Building r-index of input file " << input_file << endl;
	cout << "Index will be saved to " << idx_file << endl;

	string T;
	int64_t n;

	malloc_count_reset_peak();

	{
		std::ifstream file(input_file);
		if (!file.good()) {
			std::cout << "invalid input: could not read textfile" << std::endl;
			return -1;
		}
		file.seekg(0,std::ios::end);
		n = file.tellg()+(std::streamsize)+1;
		file.seekg(0,std::ios::beg);
		if (n <= INT_MAX) {
			T.resize(n);
			file.read((char*)&T[0],n-1);
		} else {
			std::stringstream buffer;
			buffer << file.rdbuf();
			T.reserve(n);
			T = buffer.str();
		}
		file.close();
		T[n-1] = 1;
	}

	string path = string(out_basename).append(".ri-mds");
	std::ofstream out(path);
	string text_file_name = out_basename.substr(out_basename.find_last_of("/\\") + 1);

	if (measurement_file.is_open()) {
		measurement_file << "RESULT type=ri_mds name=" << text_file_name << " p=" << p << " v=" << v << " a=" << a;
	}
	
    auto t1 = high_resolution_clock::now();

	if(hyb){

		//auto idx = r_index_mds<sparse_hyb_vector,rle_string_hyb>(input,sais);
		//idx.serialize(out);

	}else{
		if (n <= INT_MAX) {
			out << false;
			auto idx = r_index_mds<int32_t>(T,(int32_t)n,p,v,(int32_t)a,true,measurement_file.is_open() ? &measurement_file : NULL);
			idx.serialize(out);
		} else {
			out << true;
			auto idx = r_index_mds<int64_t>(T,(int64_t)n,p,v,(int64_t)a,true,measurement_file.is_open() ? &measurement_file : NULL);
			idx.serialize(out);
		}

	}

	auto t2 = high_resolution_clock::now();
	ulint total = duration_cast<duration<double, std::ratio<1>>>(t2 - t1).count();
	cout << "Build time : " << get_time(total) << endl;

	if (measurement_file.is_open()) {
		measurement_file << " time_tot=" << time_diff_ms(t1,t2) << " file_size=" << std::filesystem::file_size(std::filesystem::path(path))/1000 << endl;
	}

	out.close();

}
