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
		if (log) cout << endl << "n = " << n << endl;

		if (log) cout << "calculating SA" << endl;
		std::vector<int32_t> SA(n);
		if (p > 1) {
			libsais_omp((uint8_t*)&T[0],(int32_t*)&SA[0],n,0,NULL,p);
		} else {
			libsais((uint8_t*)&T[0],(int32_t*)&SA[0],n,0,NULL);
		}

		if (log) cout << "calculating bwt" << endl;
		string bwt_s;
		bwt_s.resize(n);
		#pragma omp parallel for num_threads(p)
		for (int i=0; i<n; i++) {
			bwt_s[i] = SA[i] == 0 ? T[n-1] : T[SA[i]-1];
		}

		if (log) cout << "calculating F" << endl;
		F = vector<ulint>(256,0);
		#pragma omp parallel for num_threads(p)
		for(uchar c : bwt_s) {
			#pragma omp atomic
			F[c]++;
		}
		for(ulint i=255;i>0;--i) {
			F[i] = F[i-1];
		}
		F[0] = 0;
		for(ulint i=1;i<256;++i) {
			F[i] += F[i-1];
		}
		#pragma omp parallel for num_threads(p)
		for(ulint i=0;i<bwt_s.size();++i) {
			if(bwt_s[i]==TERMINATOR) {
				terminator_position = i;
			}
		}

		std::vector<uint32_t> run_start_positions;
		std::vector<std::pair<uint32_t,uint32_t>> *I_phi;
		{
			if (log) cout << "calculating run start positions" << endl;
			run_start_positions.reserve(n);
			run_start_positions.push_back(0);
			for (int i=1; i<n; i++) {
				if (bwt_s[i] != bwt_s[i-1]) {
					run_start_positions.push_back(i);
				}
			}
			r = run_start_positions.size();
			if (log) cout << "r = " << r << endl;

			if (log) cout << "calculating LF-Array" << endl;
			std::vector<std::pair<uint32_t,uint32_t>> *I_LF = new std::vector<std::pair<uint32_t,uint32_t>>(r);
			{
				int_vector<64> F_(256,0);
				uchar c;
				for (uint32_t i=0; i<r; i++) {
					c = bwt_s[run_start_positions[i]];
					I_LF->at(i) = std::make_pair(run_start_positions[i],F[c]+F_[c]);
					F_[c] += run_start_positions[i+1] - run_start_positions[i];
				}
			}

			if (log) cout << "building move datastructure for LF" << endl;
			mds_LF = mds<uint32_t>(I_LF,n,4,2,p,version,log);
			I_LF = NULL;
			r = mds_LF.intervals();
			if (log) cout << "r = " << r << endl;

			if (log) cout << "adjusting run start positions" << endl;
			run_start_positions.resize(r);
			run_start_positions.shrink_to_fit();
			#pragma omp parallel for num_threads(p)
			for (uint32_t i=0; i<r; i++) {
				run_start_positions[i] = mds_LF.pair(i).first;
			}

			if (log) cout << "building phi-Array and SA-Samples" << endl;
			SA_sampl = int_vector<32>(r);
			I_phi = new std::vector<std::pair<uint32_t,uint32_t>>(r);
			SA_sampl[r-1] = SA[n-1];
			I_phi->at(0) = std::make_pair(SA[0],SA[n-1]);
			#pragma omp parallel for num_threads(p)
			for (uint32_t i=1; i<r; i++) {
				I_phi->at(i) = std::make_pair(SA[run_start_positions[i]],SA[run_start_positions[i]-1]);
				SA_sampl[i-1] = SA[run_start_positions[i]-1];
			}
		}

		if (log) cout << "sorting phi-Array" << endl;
		auto comp = [](auto p1, auto p2){return p1.first < p2.first;};
		if (p > 1) {
            ips4o::parallel::sort(I_phi->begin(),I_phi->end(),comp);
        } else {
            ips4o::sort(I_phi->begin(),I_phi->end(),comp);
        }

		if (log) cout << "builing move datastructure for phi" << endl;
		mds_phi = mds<uint32_t>(I_phi,n,4,2,p,version,log);
		I_phi = NULL;

		if (log) cout << "calculating SA-Sample indices" << endl;
		SA_sampl_idx = int_vector<32>(r);
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
			bwt_rh_s[i] = bwt_s[run_start_positions[i]];
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

		out.write((char*)&r,sizeof(ulint));
		w_bytes += sizeof(ulint);

		out.write((char*)&terminator_position,sizeof(ulint));
		w_bytes += sizeof(ulint);

		w_bytes += bwt_rh.serialize(out);

		out.write((char*)&bwt_rh_s[0],r);
		w_bytes += r;

		out.write((char*)F.data(),256*sizeof(ulint));
		w_bytes += 256*sizeof(ulint);

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

		in.read((char*)&r,sizeof(ulint));

		in.read((char*)&terminator_position,sizeof(ulint));

		bwt_rh.load(in);

		bwt_rh_s.resize(r);
		in.read((char*)&bwt_rh_s[0],r);

		F = vector<ulint>(256);
		in.read((char*)&F[0],256*sizeof(ulint));

		SA_sampl = int_vector<32>(r);
		in.read((char*)&SA_sampl[0],r*sizeof(uint32_t));

		SA_sampl_idx = int_vector<32>(r);
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
	ulint r;
	ulint terminator_position;
	vector<ulint> F;
	int_vector<32> SA_sampl;
	int_vector<32> SA_sampl_idx;
	mds<uint32_t> mds_LF;
	mds<uint32_t> mds_phi;
};

}

#endif /* R_INDEX_MDS_S_H_ */
