#pragma once

//requires:
/*
#include <fstream>
#include <string>
#include <vector>*/
//#include "stdafx.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>
using namespace std;

const double gd_eps = 0.0000000001; //1000*DBL_MIN;//0.00005;			//numerical tolerance
const double gd_bigM = 2000;			//big M value

class Params {
public:
	int n;			//number of movements
	int m;			//number of segments
	int q;			//number of transits
	double delta;	//following clearance
	double gamma;	//oncoming clearance
	double eps;		//berth clearance

	Params() {};
	Params(ifstream &af_in) {
		af_in >> n >> m >> q >> delta >> gamma >> eps;
	}//Ctor
};//Params

class Transit {
public:
	const int mi_id;
	const int mi_j;		//id of the ship this transit belongs to
	const int mi_i;		//id of the the channel segment for this transit (value of -1 indicates a berth allocation to be made)
	const int mi_k;		//k'th transit for this ship
	const int mi_dir;	//type of transit +1=in, 0=berth, -1=out
	const double md_p;	//transit time

	double md_cum_p;	//cumulative time from start of movement in this direction to the start of this transit

	Transit(int ai_id, int ai_j, int ai_i, int ai_k, int ai_dir, double ad_p, double ad_cum_p) :
		mi_id(ai_id), mi_j(ai_j), mi_i(ai_i), mi_k(ai_k), mi_dir(ai_dir), md_p(ad_p), md_cum_p(ad_cum_p) {}
};//Transit

class Segment {
public:
	const int mi_id;
	const bool mb_pass;

	Segment(int ai_id, bool ab_pass) : mi_id(ai_id), mb_pass(ab_pass) {}
};//Segment

class Ship {
public:
	const int mi_id;
	const double md_r;		//release time
	vector<Transit*> mv_Q;	//set of transits for the ship

	vector<Transit*> mv_B;	//set of transits describing options for berth occupations
							//mi_id should be the same as the Transit to be allocated
							//mi_i should be the id of the berth in the set of possible allocations

	double md_toa;			//this is used in a version of the MIP to sort vessels in turn-of-arrival order (default is md_r, but takes different values for bap-cr decomposition)

	Ship(int ai_id, double ad_r) : mi_id(ai_id), md_r(ad_r), md_toa(ad_r) {}
	~Ship() {
		for (int i = 0; i < mv_Q.size(); ++i) delete mv_Q[i];
		for (int i = 0; i < mv_B.size(); ++i) delete mv_B[i];
		mv_Q.clear();
		mv_Q.clear();
	}//dtor
	
	double md_tab() {
		if (md_berth_time >= 0)
			return md_berth_time;
		for (int i = 0; i < mv_Q.size(); ++i) {
			if (mv_Q[i]->mi_dir == 0) {
				md_berth_time = mv_Q[i]->md_p;
				return md_berth_time;
			}//if
		}//i
		return 0;
	}//md_tab

	int mi_berth() {
		if (mi_b == -2) {
			if (mv_B.size() > 0)
				mi_b = -1;		//no allocated berth for this vessel yet
			else {
				for (int a = 0; a < mv_Q.size(); ++a) {
					if (mv_Q[a]->mi_dir == 0) {
						mi_b = mv_Q[a]->mi_i;
						break;
					}//if
				}//a
			}//else
		}//if
		return mi_b;
	}//mi_berth

	void m_add_berth_option(int ai_k, int ai_berth, double ad_p) {
		Transit *p_tr = mv_Q[ai_k];
		assert(p_tr->mi_i < 0);		//this indicates that p_tr requires a berth allocation
		Transit *p_trb = new Transit(p_tr->mi_id, mi_id, ai_berth, ai_k, 0, ad_p, p_tr->md_cum_p);
		mv_B.push_back(p_trb);
	}//m_add_berth_option

private:
	double md_berth_time = -1;		//time at berth (or an estimate of earliest berth comletion time for BAP-CR)
	int mi_b = -2;
};//Ship

class Instance {
public:
	int mi_id = 0;				//general id for external use
	Params& mx_prm;
	vector<Ship*> mv_ships;
	vector<Segment*> mv_segs;
	vector<Transit*> mv_transits;

	double md_wt_penalty = 0;	//this is for decomposition methods which need to adjust their total waiting time values by a fixed amount for direct comparison to other methods

	Instance(Params& ax_prm) : mx_prm(ax_prm) {}

	Instance(Params& ax_prm, ifstream& af_in) : mx_prm(ax_prm) {
		string c_dummy;
		int i, i_ship, i_ch, i_dir, i_j, i_prev_dir;
		double d_r, d_p, d_cum_p;
		getline(af_in, c_dummy);		//skip header line for segment data
		int b_nopass;
		for (i = 0;i<ax_prm.m;++i) {
			af_in >> b_nopass;
			mv_segs.push_back(new Segment(i, !b_nopass));
		}//i
		c_dummy = af_in.get();			//need to read in the \n character
		getline(af_in, c_dummy);		//skip header line for ship data
		for (i = 0;i<ax_prm.n;++i) {
			af_in >> d_r;		//NOTE: assumes ships sorted in ascending order of r
			mv_ships.push_back(new Ship(i, d_r));
		}//i
		c_dummy = af_in.get();			//need to read in the \n character
		getline(af_in, c_dummy);		//skip header for transit data
		i_j = -1;
		for (i = 0;i<ax_prm.q;++i) {
			af_in >> i_ship >> i_ch >> i_dir >> d_p;
			if (i_ship != i_j || i_dir != i_prev_dir) {
				d_cum_p = 0;
				i_j = i_ship;
				i_prev_dir = i_dir;
			}//if
			mv_ships[i_ship]->mv_Q.push_back(new Transit(i, i_ship, i_ch, mv_ships[i_ship]->mv_Q.size(), i_dir, d_p, d_cum_p));
			d_cum_p += d_p;
		}//i
		m_initialise();
	}//Ctor	

	~Instance() {
		for (int i = 0; i < mv_ships.size(); ++i) delete mv_ships[i];
		for (int i = 0; i < mv_segs.size(); ++i) delete mv_segs[i];
		mv_ships.clear();
		mv_segs.clear();
		mv_transits.clear();
	}//dtor

	void m_initialise() {
		for (int j = 0; j < mv_ships.size(); ++j)
			for (int k = 0; k < mv_ships[j]->mv_Q.size(); ++k) {
				mv_transits.push_back(mv_ships[j]->mv_Q[k]);
				assert(mv_ships[j]->mv_Q[k]->mi_id + 1 == mv_transits.size());
			}//k
	}//m_initialise
};//Instance

  //this class holds minimal data to describe a solution to the ChSP
  //
class SolutionChSP {
public:
	string mc_tag;			//general tag for external use (e.g. method which generated the solution)
	Instance& mx_dat;
	vector<double> md_t;
	vector<int> mi_berth;	//berth allocated to each transit operation (-1 if transit in channel)
							//OR, could be used less generally on a per ship basis rather than per transit

	//TODO: provide structure for berth allocations to be recorded

	double md_ctot;
	double md_cmax;
	double md_wt_tot;
	double md_cpu_time = 0;

	double md_gap = -1;			//estimated optimality gap if known

	SolutionChSP(Instance& ax_dat) : mx_dat(ax_dat) {
		md_t.resize(ax_dat.mx_prm.q);
		mi_berth.resize(ax_dat.mx_prm.q);
		for (int a = 0; a < ax_dat.mx_prm.q; ++a)
			mi_berth[a] = -1;		//initialise to dummy value
	}//Ctor

	~SolutionChSP() {
		md_t.clear();
	}//dtor

	void m_metrics_calc() {
		md_ctot = 0;
		md_cmax = 0;
		md_wt_tot = mx_dat.md_wt_penalty;
		Transit *p_tr;
		for (int j = 0; j < mx_dat.mv_ships.size(); ++j) {
			p_tr = mx_dat.mv_ships[j]->mv_Q.back();
			md_ctot += md_t[p_tr->mi_id] + p_tr->md_p - mx_dat.mv_ships[j]->md_r;
			md_cmax = max(md_cmax, md_t[p_tr->mi_id] + p_tr->md_p);

			p_tr = mx_dat.mv_ships[j]->mv_Q.front();
			md_wt_tot += md_t[p_tr->mi_id] - mx_dat.mv_ships[j]->md_r;
			for (int k = 1; k < mx_dat.mv_ships[j]->mv_Q.size(); ++k) {
				p_tr = mx_dat.mv_ships[j]->mv_Q[k];
				if (p_tr->mi_dir == 0) {
					if (p_tr->mi_i >= 0)
						md_wt_tot += md_t[p_tr->mi_id + 1] - (md_t[p_tr->mi_id] + p_tr->md_p);
					else {
						double d_minp = mx_dat.mv_ships[j]->mv_B[0]->md_p;
						for (int h = 1; h < mx_dat.mv_ships[j]->mv_B.size(); ++h)
							d_minp = mx_dat.mv_ships[j]->mv_B[h]->md_p < d_minp ? mx_dat.mv_ships[j]->mv_B[h]->md_p : d_minp;
						md_wt_tot += md_t[p_tr->mi_id + 1] - (md_t[p_tr->mi_id] + d_minp);
					}//else
				}//if
			}//k
		}//j
	}//m_metrics_calc

	void m_dump(ostream& af_out) {
		m_metrics_calc();
		int a;
		for (a = 0; a < md_t.size(); ++a)
			af_out << md_t[a] << ",";
		af_out << endl << md_ctot << "," << md_cmax << "," << md_wt_tot << endl;
		for (a = 0; a < md_t.size(); ++a)
			af_out << "0,";
		af_out << endl;
	}//m_dump

	SolutionChSP* mp_clone(const string& ac_tag) {
		SolutionChSP* p_clone = new SolutionChSP(mx_dat);

		m_metrics_calc();
		p_clone->md_ctot = md_ctot;
		p_clone->md_cmax = md_cmax;
		p_clone->md_wt_tot = md_wt_tot;
		for (int q = 0; q < md_t.size(); ++q) {
			p_clone->md_t[q] = md_t[q];
			p_clone->mi_berth[q] = mi_berth[q];
		}//q
		p_clone->mc_tag = ac_tag;
		p_clone->md_cpu_time = md_cpu_time;
		p_clone->md_gap = md_gap;

		return p_clone;
	}//mp_clone
};//SolutionChSP
