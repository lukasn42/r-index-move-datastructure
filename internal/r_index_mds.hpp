/*
 * r_index_mds.hpp
 *
 *  Created on: Apr 13, 2017
 *      Author: nico
 *
 * Small v of the r-index: O(r) words of space, O(log(n/r)) locate time per occurrence
 *
 */

#ifndef R_INDEX_MDS_S_H_
#define R_INDEX_MDS_S_H_

#include "definitions.hpp"
#include "rle_string.hpp"
#include "sparse_sd_vector.hpp"
#include "sparse_hyb_vector.hpp"
#include "utils.hpp"

#include <omp.h>

extern "C" {
	#include <malloc_count.h>
}

#include <mds.hpp>
#include <mds.cpp>

extern "C" {
	#include <libsais.h>
	#include <libsais64.h>
}

using namespace sdsl;

std::chrono::steady_clock::time_point now() {
    return std::chrono::steady_clock::now();
}

int64_t time_diff_ms(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
}

namespace ri_mds{

template<typename INT_T>
class r_index_mds{

public:

	using triple = std::tuple<range_t, ulint, ulint>;

	r_index_mds(){}

	/*
	 * Build index
	 */
	r_index_mds(string T, INT_T n, int p, int v = 3, INT_T a = 16, bool log = true, std::ofstream *measurement_file = NULL){
		this->n = n;
		this->a = a;
		omp_set_num_threads(p);

		chars_mapped = false;
		map_char = std::vector<unsigned char>(256,0);
		#pragma omp parallel for num_threads(p)
		for (INT_T i=0; i<n-1; i++) {
			map_char[(unsigned char) T[i]] = 1;
		}
		if (map_char[0] || map_char[1]) {
			uint8_t num_dist_chars = 0;
			for (uint16_t i=0; i<256; i++) {
				num_dist_chars += map_char[i];
			}
			if (num_dist_chars > 253) {
				cout << "Error: the input contains more than 253 distinct characters" << endl;
				return;
			}
			chars_mapped = true;
			map_char_rev = std::vector<unsigned char>(256,0);
			uint8_t j = 2;
			for (uint16_t i=0; i<256; i++) {
				if (map_char[i] != 0) {
					map_char[i] = j;
					j++;
				}
			}
			for (uint16_t i=0; i<256; i++) {
				if (map_char[i] != 0) {
					map_char_rev[map_char[i]] = i;
				}
			}
			#pragma omp parallel for num_threads(p)
			for (INT_T i=0; i<n-1; i++) {
				T[i] = map_char[(unsigned char) T[i]];
			}
		} else {
			map_char.clear();
		}

		auto time = now();

		if (log) cout << "building SA" << endl;
		std::vector<INT_T> SA(n);
		if (p > 1) {
			if (std::is_same<INT_T,int32_t>::value) {
				libsais_omp((uint8_t*)&T[0],(int32_t*)&SA[0],(int32_t)n,0,NULL,p);
			} else {
				libsais64_omp((uint8_t*)&T[0],(int64_t*)&SA[0],(int64_t)n,0,NULL,p);
			}
		} else {
			if (std::is_same<INT_T,int32_t>::value) {
				libsais((uint8_t*)&T[0],(int32_t*)&SA[0],(int32_t)n,0,NULL);
			} else {
				libsais64((uint8_t*)&T[0],(int64_t*)&SA[0],(int64_t)n,0,NULL);
			}
		}
		if (measurement_file != NULL) {
			*measurement_file << " phase_1=" << time_diff_ms(time,now());
			time = now();
		}

		if (log) cout << "building bwt" << endl;
		string bwt;
		bwt.resize(n);
		#pragma omp parallel for num_threads(p)
		for(INT_T i=0; i<n; i++) {
			bwt[i] = SA[i] == 0 ? T[n-1] : T[SA[i]-1];
		}
		T.clear();

		if (log) cout << "building C" << endl;
		std::vector<INT_T> C(256,0);
		for (uchar c : bwt) {
			C[c]++;
		}
		for (int i=255; i>0; i--) {
			C[i] = C[i-1];
		}
		C[0] = 0;
		for (int i=1; i<256; i++) {
			C[i] += C[i-1];
		}

		if (log) cout << "building LF-Array" << endl;
		std::vector<std::pair<INT_T,INT_T>> *I_LF = new std::vector<std::pair<INT_T,INT_T>>();
		I_LF->reserve(n/16);
		std::vector<INT_T> C_(256,0);
		INT_T l = 0;
		uchar c = bwt[0];
		I_LF->emplace_back(std::make_pair(0,C[c]));
		for (INT_T i=1; i<n; i++) {
			if (bwt[i] != c) {
				C_[c] += i-l;
				l = i;
				c = bwt[i];
				I_LF->emplace_back(std::make_pair(i,C[c]+C_[c]));
			}
		}
		r = I_LF->size();
		I_LF->shrink_to_fit();
		C.clear();
		C_.clear();
		if (measurement_file != NULL) {
			*measurement_file << " phase_2=" << time_diff_ms(time,now());
			time = now();
		}

		if (log) cout << "building phi-Array" << endl;
		std::vector<std::pair<INT_T,INT_T>> *I_phi = new std::vector<std::pair<INT_T,INT_T>>(r);
		I_phi->at(0) = std::make_pair(SA[0],SA[n-1]);
		for (INT_T i=1; i<r; i++) {
			I_phi->at(i) = std::make_pair(SA[I_LF->at(i).first],SA[I_LF->at(i).first-1]);
		}

		if (log) cout << "sorting phi-Array" << endl;
		auto comp = [](auto p1, auto p2){return p1.first < p2.first;};
		if (p > 1) {
			ips4o::parallel::sort(I_phi->begin(),I_phi->end(),comp);
		} else {
			ips4o::sort(I_phi->begin(),I_phi->end(),comp);
		}

		if (log) cout << "building move datastructure for LF" << endl;
		M_LF = mds<INT_T>(I_LF,n,a,p,v,measurement_file == NULL ? log : false);
		r_ = M_LF.intervals();

		if (log) cout << "SA-Samples and bwt runheads" << endl;
		SA_s = std::vector<INT_T>(r_);
		S_bwtr.resize(r_);
		SA_s[r_-1] = SA[n-1];
		S_bwtr[0] = bwt[0];
		#pragma omp parallel for num_threads(p)
		for (INT_T i=1; i<r_; i++) {
			SA_s[i-1] = SA[M_LF.pair(i).first-1];
			S_bwtr[i] = bwt[M_LF.pair(i).first];
		}
		bwt.clear();
		SA.clear();

		if (log) cout << "builing move datastructure for phi" << endl;
		M_phi = mds<INT_T>(I_phi,n,a,p,v,measurement_file == NULL ? log : false);
		if (measurement_file != NULL) {
			*measurement_file << " phase_3=" << time_diff_ms(time,now());
			time = now();
		}
		
		if (log) cout << "building SA-Sample indices" << endl;
		SA_x = std::vector<INT_T>(r_);
		#pragma omp parallel for num_threads(p)
		for (INT_T i=0; i<r_; i++) {
			INT_T b = 0;
			INT_T e = M_phi.intervals()-1;
			INT_T m;
			while (b != e) {
				m = (b+e)/2+1;
				if (M_phi.pair(m).first > SA_s[i]) {
					e = m-1;
				} else {
					b = m;
				}
			}
			SA_x[i] = b;
		}
		if (measurement_file != NULL) {
			*measurement_file << " phase_4=" << time_diff_ms(time,now());
			time = now();
		}

		if (log) cout << "buildung wavelet tree of bwt runheads" << endl;
		W_bwtr = huff_string(S_bwtr);
		if (measurement_file != NULL) {
			*measurement_file << " phase_5=" << time_diff_ms(time,now());
			time = now();
			*measurement_file << " memory_usage=" << malloc_count_peak()/1000000;
		}
	}

	void revert(string &text) {
		std::pair<INT_T,INT_T> mp = std::make_pair(0,0);
		text[n-2] = S_bwtr[mp.second];
		for (INT_T i=n-3; i>=0; i--) {
			M_LF.move(mp);
			text[i] = S_bwtr[mp.second];
		}
	}

	/*
	 * Return BWT range of pattern P
	 */
	range_t count(string &P){
		INT_T m = P.size();

		std::pair<INT_T,INT_T> mp_l(0,0);
		std::pair<INT_T,INT_T> mp_r(n-1,r_-1);

		for (INT_T i=m-1; i>=0; i--) {
			if (P[i] != S_bwtr[mp_l.second]) {
				mp_l.second = W_bwtr.select(W_bwtr.rank(mp_l.second,P[i])+1,P[i]);
				mp_l.first = M_LF.pair(mp_l.second).first;
			}
			if (P[i] != S_bwtr[mp_r.second]) {
				mp_r.second = W_bwtr.select(W_bwtr.rank(mp_r.second,P[i]),P[i]);
				mp_r.first = M_LF.pair(mp_r.second+1).first-1;
			}
			if (mp_l.first <= mp_r.first) {
				M_LF.move(mp_l);
				M_LF.move(mp_r);
			} else {
				return {1,0};
			}
		}

		return {mp_l.first,mp_r.first};
	}

	/*
	 * Return number of occurrences of P in the text
	 */
	ulint occ(string &P){
		auto rn = count(P);

		if (rn.second>=rn.first) {
			return (rn.second-rn.first)+1;
		} else {
			return 0;
		}
	}

	/*
	 * locate all occurrences of P and return them in an array
	 * (space consuming if result is big).
	 */
	vector<ulint> locate_all(string& P){
		INT_T m = P.size();

		std::pair<INT_T,INT_T> mp_l(0,0);
		std::pair<INT_T,INT_T> mp_r(n-1,r_-1);
		std::pair<INT_T,INT_T> mp_sa_r(0,r_-1);

		for (INT_T i=m-1; i>=0; i--) {
			if (P[i] != S_bwtr[mp_l.second]) {
				mp_l.second = W_bwtr.select(W_bwtr.rank(mp_l.second,P[i])+1,P[i]);
				mp_l.first = M_LF.pair(mp_l.second).first;
			}
			if (P[i] != S_bwtr[mp_r.second]) {
				mp_r.second = W_bwtr.select(W_bwtr.rank(mp_r.second,P[i]),P[i]);
				mp_r.first = M_LF.pair(mp_r.second+1).first-1;
				mp_sa_r = std::make_pair(0,mp_r.second);
			}
			if (mp_l.first <= mp_r.first) {
				M_LF.move(mp_l);
				M_LF.move(mp_r);
				mp_sa_r.first--;
			} else {
				return vector<ulint>(0);
			}
		}

		mp_sa_r.first += SA_s[mp_sa_r.second];
		mp_sa_r.second = SA_x[mp_sa_r.second];

		while (mp_sa_r.first < M_phi.pair(mp_sa_r.second).first) {
			mp_sa_r.second--;
		}

		vector<ulint> Occ(mp_r.first-mp_l.first+1);

		if (!Occ.empty()) {
			Occ[0] = mp_sa_r.first;
			for (INT_T i=1; i<Occ.size(); i++) {
				M_phi.move(mp_sa_r);
				Occ[i] = mp_sa_r.first;
			}
		}

		return Occ;
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){
		ulint w_bytes = 0;

		out.write((char*)&chars_mapped,sizeof(bool));
		w_bytes += sizeof(bool);

		if (chars_mapped) {
			out.write((char*)&map_char[0],256*sizeof(unsigned char));
			w_bytes += 256*sizeof(unsigned char);

			out.write((char*)&map_char_rev[0],256*sizeof(unsigned char));
			w_bytes += 256*sizeof(unsigned char);
		}

		out.write((char*)&n,sizeof(INT_T));
		w_bytes += sizeof(INT_T);

		out.write((char*)&r,sizeof(INT_T));
		w_bytes += sizeof(INT_T);

		out.write((char*)&r_,sizeof(INT_T));
		w_bytes += sizeof(INT_T);

		out.write((char*)&a,sizeof(INT_T));
		w_bytes += sizeof(INT_T);

		w_bytes += W_bwtr.serialize(out);

		out.write((char*)&S_bwtr[0],r_);
		w_bytes += r_;

		out.write((char*)&SA_s[0],r_*sizeof(INT_T));
		w_bytes += r_*sizeof(INT_T);

		out.write((char*)&SA_x[0],r_*sizeof(INT_T));
		w_bytes += r_*sizeof(INT_T);

		w_bytes += M_LF.serialize(out);
		
		w_bytes += M_phi.serialize(out);

		return w_bytes;
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {
		in.read((char*)&chars_mapped,sizeof(bool));

		if (chars_mapped) {
			map_char = std::vector<unsigned char>(256);
			in.read((char*)&map_char[0],256*sizeof(unsigned char));

			map_char_rev = std::vector<unsigned char>(256);
			in.read((char*)&map_char_rev[0],256*sizeof(unsigned char));
		}

		in.read((char*)&n,sizeof(INT_T));

		in.read((char*)&r,sizeof(INT_T));

		in.read((char*)&r_,sizeof(INT_T));

		in.read((char*)&a,sizeof(INT_T));

		W_bwtr.load(in);

		S_bwtr.resize(r_);
		in.read((char*)&S_bwtr[0],r_);

		SA_s = std::vector<INT_T>(r_);
		in.read((char*)&SA_s[0],r_*sizeof(INT_T));

		SA_x = std::vector<INT_T>(r_);
		in.read((char*)&SA_x[0],r_*sizeof(INT_T));

		M_LF = mds<INT_T>(in);

		M_phi = mds<INT_T>(in);
	}

	/*
	 * save the structure to the path specified.
	 * \param path_prefix prefix of the index files. suffix ".ri-mds" will be automatically added
	 */
	void save_to_file(string path_prefix){

		string path = string(path_prefix).append(".ri-mds");

		std::ofstream out(path);
		serialize(out);
		out.close();

	}

	/*
	 * load the structure from the path specified.
	 * \param path: full file name
	 */
	void load_from_file(string path){

		std::ifstream in(path);
		load(in);
		in.close();

	}

	INT_T text_length() {
		return n;
	}

	INT_T ret_a() {
		return this->a;
	}

	bool ret_chars_mapped() {
		return chars_mapped;
	}

	unsigned char ret_map_char(unsigned char c) {
		return map_char[c];
	}

	unsigned char ret_map_char_rev(unsigned char c) {
		return map_char_rev[c];
	}

private:

	huff_string W_bwtr;
	string S_bwtr;
	INT_T n;
	INT_T r;
	INT_T r_;
	INT_T a;
	std::vector<INT_T> SA_s;
	std::vector<INT_T> SA_x;
	mds<INT_T> M_LF;
	mds<INT_T> M_phi;
	bool chars_mapped;
	std::vector<unsigned char> map_char,map_char_rev;
};

}

#endif /* R_INDEX_MDS_S_H_ */
