/*
 * r_index_mds.hpp
 *
 *  Created on: Apr 13, 2017
 *      Author: nico
 *
 * Small version of the r-index: O(r) words of space, O(log(n/r)) locate time per occurrence
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

#include "../external/move-datastructure/include/mds.hpp"
#include "../external/move-datastructure/src/mds.cpp"

#include "../external/move-datastructure/extern/libsais/src/libsais.h"
#include "../external/move-datastructure/extern/libsais/src/libsais.c"

using namespace sdsl;

namespace ri_mds{

template	<	class sparse_bv_type = sparse_sd_vector,
				class rle_string_t = rle_string_sd
			>
class r_index_mds{

public:

	using triple = std::tuple<range_t, ulint, ulint>;

	r_index_mds(){}

	/*
	 * Build index
	 */
	r_index_mds(string T, int32_t length, int p, int version = 3, bool log = true){
		n = length;
		omp_set_num_threads(p);

		{
			if (log) cout << "calculating SA" << endl;
			std::vector<int32_t> SA(n);
			if (p > 1) {
				libsais_omp((uint8_t*)&T[0],(int32_t*)&SA[0],n,0,NULL,p);
			} else {
				libsais((uint8_t*)&T[0],(int32_t*)&SA[0],n,0,NULL);
			}

			if (log) cout << "calculating C" << endl;
			std::vector<uint32_t> C(256,0);
			#pragma omp parallel for num_threads(p)
			for(int i=0; i<n; i++) {
				#pragma omp atomic
				C[SA[i] == 0 ? T[n-1] : T[SA[i]-1]]++;
			}
			for(uint32_t i=255;i>0;--i) {
				C[i] = C[i-1];
			}
			C[0] = 0;
			for(uint32_t i=1;i<256;++i) {
				C[i] += C[i-1];
			}

			{
				if (log) cout << "calculating LF-Array" << endl;
				std::vector<std::pair<uint32_t,uint32_t>> *I_LF = new std::vector<std::pair<uint32_t,uint32_t>>();
				{
					std::vector<uint32_t> C_(256,0);
					uint32_t l = 0;
					uint8_t c = SA[0] == 0 ? T[n-1] : T[SA[0]-1];
					I_LF->emplace_back(std::make_pair(0,C[c]));
					for (uint32_t i=1; i<n; i++) {
						if (SA[i] == 0 ? T[n-1] : T[SA[i]-1] != c) {
							C_[c] += i-l;
							l = i;
							c = SA[i] == 0 ? T[n-1] : T[SA[i]-1];
							I_LF->emplace_back(std::make_pair(i,C[c]+C_[c]));
						}
					}
				}
				C.resize(0);

				if (log) cout << "building move datastructure for LF" << endl;
				mds_LF = mds<uint32_t>(I_LF,n,2,p,version,log);
				r = mds_LF.intervals();
			}

			{
				if (log) cout << "building phi-Array" << endl;
				std::vector<std::pair<uint32_t,uint32_t>> *I_phi = new std::vector<std::pair<uint32_t,uint32_t>>(r);
				I_phi->at(0) = std::make_pair(SA[0],SA[n-1]);
				#pragma omp parallel for num_threads(p)
				for (uint32_t i=1; i<r; i++) {
					I_phi->at(i) = std::make_pair(SA[mds_LF.pair(i).first],SA[mds_LF.pair(i).first-1]);
				}
				if (log) cout << "sorting phi-Array" << endl;
				auto comp = [](auto p1, auto p2){return p1.first < p2.first;};
				if (p > 1) {
					ips4o::parallel::sort(I_phi->begin(),I_phi->end(),comp);
				} else {
					ips4o::sort(I_phi->begin(),I_phi->end(),comp);
				}

				if (log) cout << "builing move datastructure for phi" << endl;
				mds_phi = mds<uint32_t>(I_phi,n,2,p,version,log);
			}

			if (log) cout << "building SA-Samples" << endl;
			SA_sampl = std::vector<uint32_t>(r);
			SA_sampl[r-1] = SA[n-1];
			#pragma omp parallel for num_threads(p)
			for (uint32_t i=1; i<r; i++) {
				SA_sampl[i-1] = SA[mds_LF.pair(i).first-1];
			}

			if (log) cout << "calculating SA-Sample indices" << endl;
			SA_sampl_idx = std::vector<uint32_t>(r);
			#pragma omp parallel for num_threads(p)
			for (int i=0; i<r; i++) {
				uint32_t b = 0;
				uint32_t e = mds_phi.intervals()-1;
				uint32_t m;
				while (b != e) {
					m = (b+e)/2+1;
					if (mds_phi.pair(m).first > SA_sampl[i]) {
						e = m-1;
					} else {
						b = m;
					}
				}
				SA_sampl_idx[i] = b;
			}

			if (log) cout << "run-length encoding bwt" << endl;
			bwt_rh_s.resize(r);
			#pragma omp parallel for num_threads(p)
			for (int i=0; i<r; i++) {
				bwt_rh_s[i] = SA[mds_LF.pair(i).first] == 0 ? T[n-1] : T[SA[mds_LF.pair(i).first]-1];
			}
		}
		bwt_rh = huff_string(bwt_rh_s);
	}

	void revert(string &text) {
		std::pair<uint32_t,uint32_t> mp = std::make_pair(0,0);
		for (int i=n-2; i>=0; i--) {
			text[i] = bwt_rh_s[mp.second];
			mds_LF.move(mp);
		}
	}

	/*
	 * Return BWT range of pattern P
	 */
	range_t count(string &P){
		uint32_t m = P.size();

		std::pair<uint32_t,uint32_t> mp_l(0,0);
		std::pair<uint32_t,uint32_t> mp_r(n-1,r-1);

		for (int i=m-1; i>=0; i--) {
			if (P[i] != bwt_rh_s[mp_l.second]) {
				mp_l.second = bwt_rh.select(bwt_rh.rank(mp_l.second,P[i])+1,P[i]);
				mp_l.first = mds_LF.pair(mp_l.second).first;
			}
			if (P[i] != bwt_rh_s[mp_r.second]) {
				mp_r.second = bwt_rh.select(bwt_rh.rank(mp_r.second,P[i]),P[i]);
				mp_r.first = mds_LF.pair(mp_r.second+1).first-1;
			}
			if (mp_l.first <= mp_r.first) {
				mds_LF.move(mp_l);
				mds_LF.move(mp_r);
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
		uint32_t m = P.size();

		std::pair<uint32_t,uint32_t> mp_l(0,0);
		std::pair<uint32_t,uint32_t> mp_r(n-1,r-1);
		std::pair<uint32_t,uint32_t> mp_sa_r(SA_sampl[r-1],SA_sampl_idx[r-1]);

		for (int i=m-1; i>=0; i--) {
			if (P[i] != bwt_rh_s[mp_l.second]) {
				mp_l.second = bwt_rh.select(bwt_rh.rank(mp_l.second,P[i])+1,P[i]);
				mp_l.first = mds_LF.pair(mp_l.second).first;
			}
			if (P[i] != bwt_rh_s[mp_r.second]) {
				mp_r.second = bwt_rh.select(bwt_rh.rank(mp_r.second,P[i]),P[i]);
				mp_r.first = mds_LF.pair(mp_r.second+1).first-1;
				mp_sa_r.first = SA_sampl[mp_r.second];
				mp_sa_r.second = SA_sampl_idx[mp_r.second];
			}
			if (mp_l.first <= mp_r.first) {
				mds_LF.move(mp_l);
				mds_LF.move(mp_r);
				mp_sa_r.first--;
				if (mp_sa_r.first < mds_phi.pair(mp_sa_r.second).first) {
					mp_sa_r.second--;
				}
			} else {
				return vector<ulint>(0);
			}
		}

		vector<ulint> occ(mp_r.first-mp_l.first+1);

		if (!occ.empty()) {
			occ[0] = mp_sa_r.first;
			for (int i=1; i<occ.size(); i++) {
				mds_phi.move(mp_sa_r);
				occ[i] = mp_sa_r.first;
			}
		}

		return occ;
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){
		ulint w_bytes = 0;

		out.write((char*)&n,sizeof(int32_t));
		w_bytes += sizeof(int32_t);

		out.write((char*)&r,sizeof(uint32_t));
		w_bytes += sizeof(uint32_t);

		w_bytes += bwt_rh.serialize(out);

		out.write((char*)&bwt_rh_s[0],r);
		w_bytes += r;

		out.write((char*)&SA_sampl[0],r*sizeof(uint32_t));
		w_bytes += r*sizeof(uint32_t);

		out.write((char*)&SA_sampl_idx[0],r*sizeof(uint32_t));
		w_bytes += r*sizeof(uint32_t);

		w_bytes += mds_LF.serialize(out);
		
		w_bytes += mds_phi.serialize(out);

		return w_bytes;
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {
		in.read((char*)&n,sizeof(int32_t));

		in.read((char*)&r,sizeof(uint32_t));

		bwt_rh.load(in);

		bwt_rh_s.resize(r);
		in.read((char*)&bwt_rh_s[0],r);

		SA_sampl = std::vector<uint32_t>(r);
		in.read((char*)&SA_sampl[0],r*sizeof(uint32_t));

		SA_sampl_idx = std::vector<uint32_t>(r);
		in.read((char*)&SA_sampl_idx[0],r*sizeof(uint32_t));

		mds_LF = mds<uint32_t>(in);

		mds_phi = mds<uint32_t>(in);
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

	int32_t text_length() {
		return n;
	}

private:

	static const uchar TERMINATOR = 1;

	huff_string bwt_rh;
	string bwt_rh_s;
	int32_t n;
	uint32_t r;
	std::vector<uint32_t> SA_sampl;
	std::vector<uint32_t> SA_sampl_idx;
	mds<uint32_t> mds_LF;
	mds<uint32_t> mds_phi;
};

}

#endif /* R_INDEX_MDS_S_H_ */
