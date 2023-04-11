//requires:
/*
#include <fstream>
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <cassert>
using namespace std;
*/

#pragma once

#include "ChSP_base.cpp"

class TimeWindow {
public:
	double md_start;
	double md_end;
	TimeWindow(double ad_start, double ad_end) : md_start(ad_start), md_end(ad_end) {}
	TimeWindow(const TimeWindow& ax_rhs) : md_start(ax_rhs.md_start), md_end(ax_rhs.md_end) {}

	static bool mb_lt(const TimeWindow& ax_l, const TimeWindow& ax_r) {
		if (ax_l.md_start < ax_r.md_start)
			return true;
		if (ax_l.md_start == ax_r.md_start && ax_l.md_end < ax_r.md_end)
			return true;
		return false;
	}//mb_lt
};//TimeWindow

  //this class holds working data for the FCFS ChSP heuristic
  //
class WorkingChSP {
public:
	Instance& mx_dat;		//reference to problem instance data
	SolutionChSP mx_sln;	//reference to in-progress solution

	bool mb_spt;			//this is set to true if SPT (for Ch-SP) or SCT (for BAP-CR) is to be applied between ships ready to move, 
							//otherwise just standard FCFS based on release times and arbitrary ordering is applied

	vector<Ship*> mv_queue;		//FCFS queue of ships to be scheduled
	vector<TimeWindow> mv_windows;		//forbidden start movement windows for current movement to be scheduled
	bool mb_bapcr;				//indicates if the problem is BAP-CR or ChSP
	int mi_qj;					//position in the queue currently being scheduled
	Ship* mp_sh;				//current ship being scheduled
	Transit* mp_tr;				//transit at the start of the next movement to be scheduled
	double md_earliest_t;		//earliest time at which the next movement can start
	bool mb_cycle_berths;		//if true, sequencing vessels will be done cyclically through berths (only relevant if mb_bapcr = true)
	int mi_berth;			//next berth to be allocated a ship (only relevant if mb_bapcr = true and mb_cycle_berth = true)

	vector<double> md_berth_ready;		//earliest time each berth is ready (indexed by channel segment id, so that some indices will be unused)

	WorkingChSP(Instance& ax_dat, bool ab_spt) : mx_dat(ax_dat), mx_sln(ax_dat), mb_spt(ab_spt), mb_cycle_berths(true) {
		md_berth_ready.resize(ax_dat.mx_prm.m);
	}//WorkingChSP

	double md_t(Transit& ax_tr) {
		return mx_sln.md_t[ax_tr.mi_id];
	}//md_t

	 //move on to the next movement/occupation to be scheduled, update working variables and determine earliest possible starting time for the movement
	 //
	void m_next() {
		if (mp_sh != NULL && mp_tr == NULL) {	//this condition is deliberately contrived to flag initialisation
			mp_tr = mp_sh->mv_Q[0];		//this should only ever happen at initialisation
			mi_qj = 0;
		}//if
		else {
			//figure out the transit to kick off the next movement/occupation to be scheduled
			//
			int i_tr = mp_tr->mi_k;
			mp_tr = NULL;
			for (int a = i_tr + 1; a < mp_sh->mv_Q.size(); ++a) {
				if (mp_sh->mv_Q[a]->mi_dir != mp_sh->mv_Q[i_tr]->mi_dir) {
					mp_tr = mp_sh->mv_Q[a];		//we have found the first transit which describes the start of the next movement or berth occupation
					break;
				}//if
			}//q
			if (mp_tr == NULL) {
				 //now just incrmenet to the first movement of the next ship
				 //
				mp_sh = mp_sh == mv_queue.back() ? NULL : mv_queue[mi_qj + 1];		//move onto the next ship
				mp_tr = mp_sh != NULL ? mp_sh->mv_Q[0] : NULL;
				++mi_qj;
			}//if
		}//else
		if (mp_tr != NULL) {
			//determine the earliest possible starting time for this movement/occupation
			//
			if (mp_tr->mi_dir == 0) {
				if (mp_tr->mi_k == 0)		//this is a ship that starts at berth
					md_earliest_t = mb_bapcr ? mp_sh->md_r : max(mp_sh->md_r, md_berth_ready[mp_tr->mi_i]);
				else	//this is a ship that has just arrived at berth
					md_earliest_t = mx_sln.md_t[mp_tr->mi_id - 1] + mp_sh->mv_Q[mp_tr->mi_k - 1]->md_p;	//berth arrival equal to completion time of previous transit
			}
			else if (mp_tr->mi_dir == 1) {			//this is an inbound movement, earliest start depends on when berth is available
				double d_sump = 0;
				md_earliest_t = mp_tr->mi_k == 0 ? mp_sh->md_r : 0;
				for (int a = mp_tr->mi_k; a < mp_sh->mv_Q.size(); ++a) {

					if (mp_sh->mv_Q[a]->mi_dir == 0) {
						if (mp_sh->mv_Q[a]->mi_i >= 0)
							md_earliest_t = max(md_earliest_t, md_berth_ready[mp_sh->mv_Q[a]->mi_i] - d_sump);	//ChSP (berth pre-allocated)
						else {
							//BAP-CR: choose the berth yielding the earliest likely completion time
							//
							double d_earliest_c, d_min_c = gd_bigM;
							int i_alloc_b;
							for (int b = 0; b < mp_sh->mv_B.size(); ++b) {
								if (mb_spt && mb_cycle_berths)
									b = mi_berth;
								d_earliest_c = max(md_earliest_t + d_sump, md_berth_ready[b]) + mp_sh->mv_B[b]->md_p;
								if (d_earliest_c < d_min_c - gd_eps) { //floating point '<'
									d_min_c = d_earliest_c;
									i_alloc_b = b;
								}//if
								if (mb_spt && mb_cycle_berths)
									break;
							}//b
							md_earliest_t = max(md_earliest_t, md_berth_ready[i_alloc_b] - d_sump);
							mx_sln.mi_berth[mp_sh->mv_Q[a]->mi_id] = i_alloc_b;
						}//else
						break;
					}//if
					else d_sump += mp_sh->mv_Q[a]->md_p;
				}//a
			}
			else {		//this is an outbound movement, can depart after berth occupation finished
				if (mp_tr->mi_k == 0)
					md_earliest_t = mp_sh->md_r;
				else {
					if (mp_sh->mv_Q[mp_tr->mi_k - 1]->mi_i >= 0)
						md_earliest_t = mx_sln.md_t[mp_tr->mi_id - 1] + mp_sh->mv_Q[mp_tr->mi_k - 1]->md_p;
					else
						md_earliest_t = mx_sln.md_t[mp_tr->mi_id - 1] + mp_sh->mv_B[mx_sln.mi_berth[mp_tr->mi_id - 1]]->md_p;
				}//else
			}//else
		}//if
	}//m_next

	 //This procedure determines the forbidden windows for the current ship (ax_wrk.mp_sh) about to start transit (ax_wrk.mp_tr)
	 //
	void m_calculate_forbidden_windows() {
		Transit* p_tr;
		Transit *p_sh_tr;
		Ship* p_sh;
		TimeWindow* p_tw;
		int i_mv_dir;
		double d_t;

		assert(mp_tr->mi_dir != 0);
		mv_windows.clear();
		p_tw = NULL;
		for (int j = 0; j < mi_qj; ++j) {	//iterate through ships already scheduled
			i_mv_dir = 0;
			p_sh = mv_queue[j];
			for (int k = 0; k < p_sh->mv_Q.size(); ++k) {	//iterate through the movements of the previously scheduled ship
				p_tr = p_sh->mv_Q[k];
				if (p_tr->mi_dir == 0) {
					i_mv_dir = 0;
					continue;		//berth occupation doesn't create conflict in the channel, so move on to the next transit
				}//if
				if (p_tr->mi_dir != i_mv_dir) {	//this is the start of a new movement, so initialise a new time window
					if (p_tw == NULL) {
						mv_windows.emplace_back(gd_bigM, -gd_bigM);
						p_tw = &(mv_windows.back());
					}//if
					else if (p_tw->md_start < p_tw->md_end) {		//only create a new window if the previous one had width >0
						mv_windows.emplace_back(gd_bigM, -gd_bigM);
						p_tw = &(mv_windows.back());
					}//if
					else {
						p_tw->md_start = gd_bigM;
						p_tw->md_end = -gd_bigM;
					}//else
					i_mv_dir = p_tr->mi_dir;
				}//b_new_window
				if (p_tr->mi_dir == mp_tr->mi_dir) {	//this is a following headway situation
														//find the transit of mp_ship's movment that occurs in this same segment
														//
					p_sh_tr = NULL;
					for (int kk = mp_tr->mi_k; kk < mp_sh->mv_Q.size(); ++kk) {
						if (mp_sh->mv_Q[kk]->mi_dir == i_mv_dir && mp_sh->mv_Q[kk]->mi_i == p_tr->mi_i) {
							p_sh_tr = mp_sh->mv_Q[kk];
							break;		//p_sh_tr is mp_ship's transit in the same direction and channel as p_tr
						}//if
					}//kk
					if (p_sh_tr != NULL) {
						d_t = md_t(*p_tr) - p_sh_tr->md_cum_p;									//check start of transit p_tr
						p_tw->md_start = min(p_tw->md_start, d_t - mx_dat.mx_prm.delta);		//use min because following time windows in different segments combine to form one window
						p_tw->md_end = max(p_tw->md_end, d_t + mx_dat.mx_prm.delta);
						d_t = md_t(*p_tr) + p_tr->md_p - p_sh_tr->md_cum_p - p_sh_tr->md_p;		//now check end of transit p_tr (typically a redundant step)
						p_tw->md_start = min(p_tw->md_start, d_t - mx_dat.mx_prm.delta);
						p_tw->md_end = max(p_tw->md_end, d_t + mx_dat.mx_prm.delta);
					}//if
				}//if
				else if (p_tr->mi_dir == -mp_tr->mi_dir && !mx_dat.mv_segs[p_tr->mi_i]->mb_pass) {		//this is an opposing conflict situation
																										//find the transit of mp_ship's movement that occurs in the same segment
																										//
					p_sh_tr = NULL;
					for (int kk = mp_tr->mi_k; kk < mp_sh->mv_Q.size(); ++kk) {
						if (mp_sh->mv_Q[kk]->mi_dir == -i_mv_dir && mp_sh->mv_Q[kk]->mi_i == p_tr->mi_i) {
							p_sh_tr = mp_sh->mv_Q[kk];
							break;		//p_sh_tr is mp_ship's transit in the opposite direction and same channel as p_tr
						}//if
					}//kk
					if (p_sh_tr != NULL) {
						d_t = md_t(*p_tr) - (p_sh_tr->md_cum_p + p_sh_tr->md_p);				//check start of transit p_tr
						p_tw->md_start = d_t - mx_dat.mx_prm.gamma;
						d_t = md_t(*p_tr) + p_tr->md_p - p_sh_tr->md_cum_p;					//now check end of transit p_tr
						p_tw->md_end = d_t + mx_dat.mx_prm.gamma;
						mv_windows.emplace_back(gd_bigM, -gd_bigM);
						p_tw = &(mv_windows.back());
					}//if
				}//else if
			}//k
		}//j
		if (p_tw != NULL && p_tw->md_start > p_tw->md_end)
			mv_windows.pop_back();		//clean up: get rid of the last timewindow as it is not valid

		sort(mv_windows.begin(), mv_windows.end(), TimeWindow::mb_lt);
		d_t = md_earliest_t;

		//now calculate earliest departure time
		for (int w = 0; w < mv_windows.size(); ++w) {
			if (mv_windows[w].md_end > d_t) {
				if (mv_windows[w].md_start < d_t)
					d_t = mv_windows[w].md_end;
				else
					break;
			}//if
		}//w
		md_earliest_t = d_t;
	}//m_calculate_forbidden_windows

	static bool mb_ship_order(Ship* ap_lhs, Ship* ap_rhs) { 
		if (ap_lhs->md_r == ap_rhs->md_r) {
			//arbitrary tie breaking
			//
			return ap_lhs->mi_id < ap_rhs->mi_id;
			//return ap_lhs->md_tab() < ap_rhs->md_tab();  //old SPT tie breaking
		}//if
		return ap_lhs->md_r < ap_rhs->md_r;
	}//mb_ship_order

	//select ship with next release date
	//for tie breaker, select ship with earliest possible berth completion time
	//
	void m_select_next_ship() {
		vector<Ship*> v_tmp;
		int j;

		assert(mb_spt);

		if (mi_qj == mv_queue.size() - 1)
			return;

		for (j = 0; j < mv_queue.size(); ++j)
			v_tmp.push_back(mv_queue[j]);
		mv_queue.clear();

		//copy the vessels already scheduled back into mv_queue
		//
		for (j = 0; j <= mi_qj; ++j)
			mv_queue.push_back(v_tmp[j]);

		//now find the unscheduled vessel to go next (earliest release time, tie breaker: SPT for ChSP, ECT for BAP-CR)
		//
		double d_min_r = gd_bigM;
		double d_min_c = gd_bigM;
		double d_earliest_c;
		int i_next_j, b;

		for (j = mi_qj + 1; j < v_tmp.size(); ++j) {
			if (mb_bapcr) {
				for (b = 0; b < v_tmp[j]->mv_B.size(); ++b) {
					if (mb_cycle_berths)
						b = mi_berth;		//override the loop if we are sequencing vessels cyclically through berths
					d_earliest_c = max(v_tmp[j]->md_r + v_tmp[j]->mv_B[b]->md_cum_p, md_berth_ready[b]) + v_tmp[j]->mv_B[b]->md_p;
					if (d_earliest_c < d_min_c - gd_eps) {
						i_next_j = j;
						d_min_c = d_earliest_c;
					}//if
					if (mb_cycle_berths) { //WARNING: this isn't bulletproof if different vessels have different sets of berth options
						break;
					}//if
				}//b
			}//if
			else {
				for (b = 0; b < v_tmp[j]->mv_Q.size(); ++b) {
					if (v_tmp[j]->mv_Q[b]->mi_dir == 0)
						break;
				}//b
				d_earliest_c = max(v_tmp[j]->md_r + v_tmp[j]->mv_Q[b]->md_cum_p, md_berth_ready[v_tmp[j]->mv_Q[b]->mi_i]) + v_tmp[j]->mv_Q[b]->md_p;
				if (d_earliest_c < d_min_c - gd_eps) {
					i_next_j = j;
					d_min_c = d_earliest_c;
				}//if
			}//else
		}//j

		//now reconstruct mv_queue with v_tmp[i_next_j] as the next ship to be scheduled
		//
		mv_queue.push_back(v_tmp[i_next_j]);
		for (j = mi_qj + 1; j < v_tmp.size(); ++j) {
			if (j != i_next_j)
				mv_queue.push_back(v_tmp[j]);
		}//j
	}//m_select_next_ship

	void m_ChSP_FCFS() {
		double t;

		//create ship queue
		//
		mv_queue.clear();
		mb_bapcr = false;
		for (int j = 0; j < mx_dat.mx_prm.n; ++j) {
			mv_queue.push_back(mx_dat.mv_ships[j]);
			mb_bapcr = mb_bapcr || mx_dat.mv_ships[j]->mv_B.size() > 0;
		}//j
		if (!mb_spt)
			sort(mv_queue.begin(), mv_queue.end(), WorkingChSP::mb_ship_order);	//fcfs

		//initialise
		//
		mi_berth = 0;		//next berth to be allocated a ship (only relevant for mb_bapcr)
		for (int i = 0; i < mx_dat.mx_prm.m; ++i)
			md_berth_ready[i] = 0;
		mi_qj = -1;
		if (mb_spt)
			m_select_next_ship();
		mp_sh = mv_queue[0];
		mp_tr = NULL;
		m_next();

		//main loop through each ship
		//
		bool b_do_next_ship;
		while (mp_sh != NULL) {
			if (mp_tr->mi_dir != 0) {
if (mp_tr->mi_dir == -1 && mp_tr->mi_j == 9)
int i=0;
				m_calculate_forbidden_windows();
				//
				//now go through successive transits in the movment and set the transit start times
				//
				mx_sln.md_t[mp_tr->mi_id] = md_earliest_t;
				t = md_earliest_t + mp_sh->mv_Q[mp_tr->mi_k]->md_p;
				b_do_next_ship = mp_tr->mi_k == mp_sh->mv_Q.size() - 1;
				for (int k = mp_tr->mi_k + 1; k < mp_sh->mv_Q.size(); ++k) {
					if (mp_sh->mv_Q[k]->mi_dir != mp_tr->mi_dir)
						break;
					mx_sln.md_t[mp_sh->mv_Q[k]->mi_id] = t;
					t += mp_sh->mv_Q[k]->md_p;
					b_do_next_ship = k == mp_sh->mv_Q.size() - 1;
				}//k
				if (b_do_next_ship) {
					//so the next movement to be scheduled belongs to the next ship, so first make required updates to md_berth_ready from scheduling the current ship
					//
					for (int a = 0; a < mp_sh->mv_Q.size(); ++a) {
						if (mp_sh->mv_Q[a]->mi_dir == 0) {
							int b = -1;
							if (mp_sh->mv_Q[a]->mi_i >= 0)
								b = mp_sh->mv_Q[a]->mi_i;		//ChSP
							else
								b = mx_sln.mi_berth[mp_sh->mv_Q[a]->mi_id];			//BAP-CR
							md_berth_ready[b] = mx_sln.md_t[mp_sh->mv_Q[a + 1]->mi_id] + mx_dat.mx_prm.eps;
							break;
						}//if
					}//a
					if (mb_spt) {
						if (mb_cycle_berths) { //WARNING: this isn't bulletproof if different vessels have different sets of berth options
							if (++mi_berth == mp_sh->mv_B.size())
								mi_berth = 0;
						}//if
						m_select_next_ship();
					}//if
				}//if
			}//if
			else
				mx_sln.md_t[mp_tr->mi_id] = md_earliest_t;		//this is just a berth occupation (NOTE: berth alloc previously written for BAP-CR in m_next)
			m_next();
		}//wend
	}//gp_ChSP_FCFS
};//WorkingChSP
