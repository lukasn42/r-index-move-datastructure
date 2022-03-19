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

#include <mds.hpp>
#include <mds.cpp>

extern "C" {
	#include <libsais.h>
	#include <libsais64.h>
}

using namespace sdsl;

namespace ri_mds{

template<typename INT_T>
class r_index_mds{

public:

	using triple = std::tuple<range_t, ulint, ulint>;

	r_index_mds(){}

	/*
	 * Build index
	 */
	r_index_mds(string T, INT_T n, int p, int version = 3, int a = 8, bool log = true){
		this->n = n;
		omp_set_num_threads(p);

		{
			if (log) cout << "building SA" << endl;
			std::vector<INT_T> SA(n);
			if (p > 1) {
				if (std::is_same<INT_T,uint32_t>::value) {
					libsais_omp((uint8_t*)&T[0],(int32_t*)&SA[0],(int32_t)n,0,NULL,p);
				} else {
					libsais64_omp((uint8_t*)&T[0],(int64_t*)&SA[0],(int64_t)n,0,NULL,p);
				}
			} else {
				if (std::is_same<INT_T,uint32_t>::value) {
					libsais((uint8_t*)&T[0],(int32_t*)&SA[0],(int32_t)n,0,NULL);
				} else {
					libsais64((uint8_t*)&T[0],(int64_t*)&SA[0],(int64_t)n,0,NULL);
				}
			}
			if (log) cout << "building bwt" << endl;
			string bwt;
			bwt.resize(n);
			#pragma omp parallel for num_threads(p)
			for(int64_t i=0; i<n; i++) {
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

			{
				if (log) cout << "building LF-Array" << endl;
				std::vector<std::pair<INT_T,INT_T>> *I_LF = new std::vector<std::pair<INT_T,INT_T>>();
				I_LF->reserve(n/16);
				{
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
				}
				I_LF->shrink_to_fit();
				C.clear();

				if (log) cout << "building move datastructure for LF" << endl;
				mds_LF = mds<INT_T>(I_LF,n,a,p,version,log);
				r = mds_LF.intervals();
			}

			{
				if (log) cout << "building phi-Array" << endl;
				std::vector<std::pair<INT_T,INT_T>> *I_phi = new std::vector<std::pair<INT_T,INT_T>>(r);
				I_phi->at(0) = std::make_pair(SA[0],SA[n-1]);
				#pragma omp parallel for num_threads(p)
				for (INT_T i=1; i<r; i++) {
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
				mds_phi = mds<INT_T>(I_phi,n,a,p,version,log);
			}

			if (log) cout << "building SA-Samples and run-length-bwt" << endl;
			SA_sampl = std::vector<INT_T>(r);
			bwt_rh_s.resize(r);
			bwt_rh_s[0] = bwt[0];
			SA_sampl[r-1] = SA[n-1];
			#pragma omp parallel for num_threads(p)
			for (INT_T i=1; i<r; i++) {
				SA_sampl[i-1] = SA[mds_LF.pair(i).first-1];
				bwt_rh_s[i] = bwt[mds_LF.pair(i).first];
			}
		}
		
		if (log) cout << "building SA-Sample indices" << endl;
		SA_sampl_idx = std::vector<INT_T>(r);
		#pragma omp parallel for num_threads(p)
		for (INT_T i=0; i<r; i++) {
			INT_T b = 0;
			INT_T e = mds_phi.intervals()-1;
			INT_T m;
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
		if (log) cout << "buildung wavelet tree of run-length-bwt" << endl;
		bwt_rh = huff_string(bwt_rh_s);
	}

	void revert(string &text) {
		std::pair<INT_T,INT_T> mp = std::make_pair(0,0);
		text[n-2] = bwt_rh_s[mp.second];
		for (int64_t i=n-3; i>=0; i--) {
			mds_LF.move(mp);
			text[i] = bwt_rh_s[mp.second];
		}
	}

	/*
	 * Return BWT range of pattern P
	 */
	range_t count(string &P){
		INT_T m = P.size();

		std::pair<INT_T,INT_T> mp_l(0,0);
		std::pair<INT_T,INT_T> mp_r(n-1,r-1);

		for (int64_t i=m-1; i>=0; i--) {
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
		INT_T m = P.size();

		std::pair<INT_T,INT_T> mp_l(0,0);
		std::pair<INT_T,INT_T> mp_r(n-1,r-1);
		std::pair<INT_T,INT_T> mp_sa_r(SA_sampl[r-1],SA_sampl_idx[r-1]);

		for (int64_t i=m-1; i>=0; i--) {
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
			for (int64_t i=1; i<occ.size(); i++) {
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

		out.write((char*)&n,sizeof(INT_T));
		w_bytes += sizeof(INT_T);

		out.write((char*)&r,sizeof(INT_T));
		w_bytes += sizeof(INT_T);

		w_bytes += bwt_rh.serialize(out);

		out.write((char*)&bwt_rh_s[0],r);
		w_bytes += r;

		out.write((char*)&SA_sampl[0],r*sizeof(INT_T));
		w_bytes += r*sizeof(INT_T);

		out.write((char*)&SA_sampl_idx[0],r*sizeof(INT_T));
		w_bytes += r*sizeof(INT_T);

		w_bytes += mds_LF.serialize(out);
		
		w_bytes += mds_phi.serialize(out);

		return w_bytes;
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {
		in.read((char*)&n,sizeof(INT_T));

		in.read((char*)&r,sizeof(INT_T));

		bwt_rh.load(in);

		bwt_rh_s.resize(r);
		in.read((char*)&bwt_rh_s[0],r);

		SA_sampl = std::vector<INT_T>(r);
		in.read((char*)&SA_sampl[0],r*sizeof(INT_T));

		SA_sampl_idx = std::vector<INT_T>(r);
		in.read((char*)&SA_sampl_idx[0],r*sizeof(INT_T));

		mds_LF = mds<INT_T>(in);

		mds_phi = mds<INT_T>(in);
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

private:

	huff_string bwt_rh;
	string bwt_rh_s;
	INT_T n;
	INT_T r;
	std::vector<INT_T> SA_sampl;
	std::vector<INT_T> SA_sampl_idx;
	mds<INT_T> mds_LF;
	mds<INT_T> mds_phi;
};

}

#endif /* R_INDEX_MDS_S_H_ */
