/***
 *  $Id: parabayl.cpp 565 2010-04-01 21:21:45Z olia $
 **
 *  File: parabayl.cpp
 *  Created: May 22, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *          Olga Nikolova <olga.nikolova@gmail.com>
 */

#include "parabayl.hpp"


template <typename T>
inline bool inset(T set, unsigned int i) { return (set & (static_cast<T>(1) << i)); }


// Method to generate elements which have one bit less set than val
// We are starting resetting from k-th bit of val are interesting
template <typename Int> class rm_bit {
public:
    rm_bit(Int val, unsigned int k) : val_(val), k_(k) {
	mask_ = (static_cast<Int>(1) << k);
    } // rm_bit

    bool next() {
	mask_ >>= 1;
	while ((mask_ != 0) && (val_ & (~mask_)) == val_) mask_ >>= 1;
	return (mask_ != 0);
    } // next

    Int operator()() { return (val_ & (~mask_)); }

private:
    Int mask_;
    Int val_;
    unsigned int k_;

}; // class rm_bit


// Method to generate elements which have one bit more set than val
// We are starting resetting from k-th bit of val are interesting
template <typename Int> class add_bit {
public:
    add_bit(Int val, unsigned int k) : val_(val), k_(k) {
	mask_ = (static_cast<Int>(1) << k);
    } // add_bit

    bool next() {
	mask_ >>= 1;
	while ((mask_ != 0) && (val_ | mask_) == val_) mask_ >>= 1;
	return (mask_ != 0);
    } // next

    Int operator()() { return (val_ | mask_); }

private:
    Int mask_;
    Int val_;
    unsigned int k_;

}; // class add_bit




bool parabayl(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
	      InputData& data) {
    std::cout << master << "* computing network..." << std::endl;
    TIMER_START;

    unsigned int rank = mpi_env.rank();

    // n is the number of variables
    unsigned int n = data.g_nrow;

    // k is s.t. number_of_processors = 2^k
    unsigned int k = static_cast<unsigned int>(log2(mpi_env.size()));

    // total number of parallel steps
    bit_t t_max = k + (static_cast<bit_t>(1) << (n - k));
    std::cout << master << "n: " << n << ", k: " << k << ", t_max: " << t_max << ", d: " 
	      << app_conf.deg_bound << std::endl;

    if (app_conf.sync_step > 0) {
      std::cout << master << " syncing every " << app_conf.sync_step << " iterations" << std::endl;
    }

    // progress bar
    unsigned int ps = static_cast<unsigned int>(ceil((1.0 * t_max - k) / 10));

    // each processor has fixed in- and out- degree
    // and communicates to predefined neighbors
    unsigned int in_degree = jaz::bit_count(rank);
    unsigned int out_degree = k - in_degree;

    // MPI tags used in this function
    const int PBL_F_TAG = 111;

#ifdef PARABAYL_BACKTRACK
    const int PBL_PI_TAG = 112;
    const int PBL_B_TAG = 113;
    const int PBL_NET_TAG = 114;
#endif // PARABAYL_BACKTRACK

    // MPI Request & Status arrays
#ifdef PARABAYL_BACKTRACK
    const unsigned int PBL_ST_RQ_SZ = 4*in_degree;
#else
    const unsigned int PBL_ST_RQ_SZ = in_degree;
#endif //PARABAYL_BACKTRACK

    MPI_Status* in_stats = new MPI_Status[PBL_ST_RQ_SZ];
    MPI_Request* in_reqs = new MPI_Request[PBL_ST_RQ_SZ];

    // each rank is stored on k bits
    // source ranks have exactly one bit less set to 1
    // destination ranks have exactly one bit more set to 1
    // destination and source ranks remain the same
    // at each time step for a given processor

    // create source ranks
    std::vector<unsigned int> in_ranks;
    in_ranks.reserve(in_degree);
    rm_bit<bit_t> rb(rank, k);
    while (rb.next() == true) in_ranks.push_back(rb());

    // create destination ranks
    std::vector<unsigned int> out_ranks;
    out_ranks.reserve(out_degree);
    add_bit<bit_t> ab(rank, k);
    while (ab.next() == true) out_ranks.push_back(ab());

    // counter on n - k bits
    // initial value is 0
    BitCounter<bit_t> cnt(n - k);

    // level iteration counter (same as \mu(cnt))
    // at each level we keep track of
    // iterations thanks to this we know
    // in which iteration which cnt occurred
    bit_t l_iter = 0;

    // here iteration 0 starts

    // set A which is currently being computed
    jaz::plain_array<unsigned int> A;

    // we compute all local scores; A is the empty set
    //MDL_score mdl(data.g_nrow, data.g_ncol, data.D);
    BDE_score bde(data.g_nrow, data.g_ncol, data.D);

    // structures to maintain the F-scores
    // (1) F-scores to be sent
    jaz::plain_array<double> F;
    // (2) F-scores received from source ranks
    jaz::plain_array<double> Fin;
    // (3) F-scores stored locally
    // to be used by the same processor
    jaz::plain_array<double> Fl;

#ifdef PARABAYL_BACKTRACK
    // structures to maintain the ordering array pi
    // (1) permutation to be sent
    jaz::plain_array<unsigned char> pi;
    // (2) permutations received from source ranks
    jaz::plain_array<unsigned char> piin;
    // (3) permutations stored locally
    jaz::plain_array<unsigned char> pil;

    // structures to maintain subsets B
    // which consist of the parents that gave
    // the corresponding optimal F-scores at the
    // previous computation step
    // (1) B's to be sent
    jaz::plain_array<bit_t> B;
    // (2) B's received from source ranks
    jaz::plain_array<bit_t> Bin;
    // (3) B's stored locally
    jaz::plain_array<bit_t> Bl;

    // structures for maintaining the optimal network so far
    // (1) net to be sent
    jaz::plain_array<bit_t> net;
    // (2) net's received from source ranks
    jaz::plain_array<bit_t> netin;
    // (3) net's stored locally
    jaz::plain_array<bit_t> netl;
#endif // PARABAYL_BACKTRACK

    // this map stores mapping from cnt to l_iter
    // thanks to this when processing Fl (pil, netl, resp.)
    // based on ncnt we can easily get the position
    // of the corresponding F chunk (pi, net, resp.)
    // because we only store the previous (cnt_idx_pos = 0)
    // and current level (cnt_idx_pos ^ 1) of the lattice
    // we use two such maps:
    unsigned int cnt_idx_pos = 0;
    std::map<bit_t, bit_t> cnt_idx[2];

    // Some commonly used constants:

    // constant to represent the empty set as a parent set
    const bit_t EMPTY_SET = std::numeric_limits<bit_t>::max();

    // VPARENT - dummy value used to simplify sending and
    // receiving of pi; we ignore it later
    const unsigned char VPARENT = std::numeric_limits<unsigned char>::max();

    // Currently optimal elements:

    // here we store index of the optimal element
    // (X_i^* in the paper) & the optimal score
    unsigned int g_opt = 0;
    double Q_opt = std::numeric_limits<double>::max();

#ifdef PARABAYL_BACKTRACK
    // optimal permutation so far  keep track of the
    // best permutation found so far
    unsigned char* pi_opt = 0;

    bit_t B_g_opt = EMPTY_SET;
    bit_t* net_opt = 0;
#endif // PARABAYL_BACKTRACK

    // temporary variable used for indexing arrays
    unsigned int pos = 0;

    // HERE WE GO LADIES!!!
    bit_t t = 0;

    for (; t < t_max; ++t) {
	// given process is active only for 2^(n - k),
	// which is how many times cnt will be incremented
	// note that cnt counts in the BitCounter way (see class)
	if ((in_degree <= t) && (t < (in_degree + t_max - k))) {
	    // lid (lattice-id) is what we call omega in the paper
	    // lower-order k bits represent the processor id,
	    // which is followed by (n - k) bits from cnt
	    // A consists of variables assigned to the processor
	    // |A| is equal to the number of 1's in lid (|A| = \mu(omega))
	    bit_t lid = ((cnt << k) | rank);
	    A.resize(jaz::bit_count(lid));

	    // STEP 1: start receiving

	    // this is how many element we get from each parent
	    // we must receive complement of parent_of_A and
	    // Q-score attached at the very end
	    unsigned int blck_in_sz = n - A.size() + 2;
	    Fin.resize(in_degree * blck_in_sz);
	    Fin.zero();

#ifdef PARABAYL_BACKTRACK
	    // pi consists of |A| - 1 elements and VPARENT at the beginning
	    unsigned int pi_in_sz = A.size();
	    piin.resize(in_degree * pi_in_sz);
	    piin.zero();

	    unsigned int B_in_sz = n - A.size() + 1;
	    Bin.resize(in_degree * B_in_sz);
	    Bin.zero();

	    unsigned int net_in_sz = pi_in_sz;
	    netin.resize(in_degree * net_in_sz);
	    netin.zero();
#endif // PARABAYL_BACKTRACK


	    for (unsigned int i = 0; i < in_degree; ++i) {
		MPI_Irecv(Fin.begin() + i * blck_in_sz, blck_in_sz, MPI_DOUBLE, in_ranks[i], PBL_F_TAG, mpi_env.comm(), &in_reqs[i]);

#ifdef PARABAYL_BACKTRACK
		MPI_Irecv(piin.begin() + i * pi_in_sz, pi_in_sz, MPI_UNSIGNED_CHAR, in_ranks[i], PBL_PI_TAG, mpi_env.comm(), &in_reqs[in_degree + i]);
		MPI_Irecv(Bin.begin() + i * B_in_sz, B_in_sz, MPI_BIT_T, in_ranks[i], PBL_B_TAG, mpi_env.comm(), &in_reqs[in_degree * 2 + i]);
		MPI_Irecv(netin.begin() + i * net_in_sz, net_in_sz, MPI_BIT_T, in_ranks[i], PBL_NET_TAG, mpi_env.comm(), &in_reqs[in_degree * 3 + i]);
#endif // PARABAYL_BACKTRACK

	    } // for i

	    // STEP 2: compute s and F

	    // generate set A
	    pos = 0;
	    for (unsigned int i = 0; i < n; ++i) if (inset(lid, i)) A[pos++] = i;

	    F.resize(n - A.size());
	    F.zero();

#ifdef PARABAYL_BACKTRACK
	    B.resize(n - A.size());
	    B.zero();
#endif // PARABAYL_BACKTRACK
	    
	    // for elements not in A, compute s(i, A)
	    // we hope that this part will take long enough
	    // to hide communication (if needed)
	    pos = 0;
	    
	    // compute s() only if in non-restricted in-degree case
	    if (A.size() <= app_conf.deg_bound) {
	      for (unsigned int i = 0; i < n; ++i) {
		if (!inset(lid, i)) {
		  //F[pos] = mdl(i, A.begin(), A.end());
		  F[pos] = bde(i, A.begin(), A.end());
		  
#ifdef PARABAYL_BACKTRACK
		  B[pos] = (lid == 0) ? EMPTY_SET : lid;
#endif // PARABAYL_BACKTRACK
		  
		  pos++;
		}
	      } // for i
	    } else{
	      for (unsigned int i = 0; i < n; ++i) {
		if (!inset(lid, i)) {
		  F[pos] = std::numeric_limits<double>::max();
		  
#ifdef PARABAYL_BACKTRACK
		  B[pos] = EMPTY_SET;
#endif // PARABAYL_BACKTRACK
		  
		  pos++;
		}
	      } // for i
	    }// if A.size()
	    
	    //}
	// reset optimal score and permutation pointer
	Q_opt = std::numeric_limits<double>::max();
	g_opt = 0;
	
#ifdef PARABAYL_BACKTRACK
	pi_opt = 0;
#endif // PARABAYL_BACKTRACK
	
	// (1) now, we consider F that happened to be
	// local to us :-( the previous level is always stored
	    // in the beginning of the Fl array
	    // everything is the same as above, except that
	    // instead of F we operate on Fl
	    rb = rm_bit<bit_t>(static_cast<bit_t>(cnt), n - k);
	    while (rb.next() == true) {
		// ncnt is the cnt from the previous level
		// because rank is the same (we are dealing with local data)
		// the only difference is between cnt and ncnt
		bit_t ncnt = rb();

		double* Fcur = Fl.begin() + cnt_idx[cnt_idx_pos][ncnt] * blck_in_sz;

#ifdef PARABAYL_BACKTRACK
		bit_t* Bcur = Bl.begin() + cnt_idx[cnt_idx_pos][ncnt] * B_in_sz;
		bit_t* netcur = netl.begin() + cnt_idx[cnt_idx_pos][ncnt] * net_in_sz;
#endif // PARABAYL_BACKTRACK

		// current position in Fcur & Bcur
		unsigned int cur_pos = 0;

		pos= 0;

		unsigned int var = static_cast<unsigned int>(log2((cnt - ncnt) << k));
		bit_t lid_in = ((ncnt << k) | rank);

		for (unsigned int j = 0; j < n; ++j) {
		    if (!inset(lid_in, j)) {

			if (j == var) {
			    double Q = Fcur[cur_pos] + Fcur[blck_in_sz - 1];
			    if (Q < Q_opt) {
				Q_opt = Q;
				g_opt = j;
#ifdef PARABAYL_BACKTRACK
				pi_opt = pil.begin() + cnt_idx[cnt_idx_pos][ncnt] * pi_in_sz;
				B_g_opt = Bcur[cur_pos];
				net_opt = netcur;
#endif // PARABAYL_BACKTRACK

			    }
			} else {
			    if (Fcur[cur_pos] < F[pos]) {
				F[pos] = Fcur[cur_pos];

#ifdef PARABAYL_BACKTRACK
				B[pos] = Bcur[cur_pos];
#endif // PARABAYL_BACKTRACK

			    }
			    pos++;
			} // if j
			cur_pos++;
		    } // if !inset
		} // for j

	    } // while rb

	    // here we finish communication
	    MPI_Waitall(PBL_ST_RQ_SZ, in_reqs, in_stats);


	    // (2) here we go over received data in order
	    // to find optimal score and compute part of F

	    // consider all received F, in processor order
	    for (unsigned int i = 0; i < in_degree; ++i) {
		// Fcur points to F scores received from
		// processor with rank in_ranks[i]
		double* Fcur = Fin.begin() + i * blck_in_sz;
		// cur position of pi will be computed on the fly

#ifdef PARABAYL_BACKTRACK
		bit_t* Bcur = Bin.begin() + i * B_in_sz;
		bit_t* netcur = netin.begin() + i * net_in_sz;
#endif // PARABAYL_BACKTRACK

		// current position in Fcur & Bcur
		unsigned int cur_pos = 0;
		pos = 0;

		// var is a little tricky bastard
		// it is the element that we received but it is already in A
		// hence it is not used in computing F but instead to find g_opt
		// since rank differs only by one bit from in_ranks[i] and
		// the difference is due to the tricky element we can
		// use log2 of difference to get it
		// we look into ranks and not lid because upper bits of
		// lid are the same
		unsigned int var = static_cast<unsigned int>(log2(rank - in_ranks[i]));

		// lid of the received set
		bit_t lid_in = ((cnt << k) | in_ranks[i]);

		// iterate over received elements
		// these elements are complement of (A \ {var})
		// so var is actually in this complement
		for (unsigned int j = 0; j < n; ++j) {
		    if (!inset(lid_in, j)) {
			if (j == var) {
			    // var (or j if you like) is in our A
			    // hence it contributes to the optimal
			    // score and not to F to be sent out
			    double Q = Fcur[cur_pos] + Fcur[blck_in_sz - 1];
			    if (Q < Q_opt) {
				g_opt = j;
				Q_opt = Q;

#ifdef PARABAYL_BACKTRACK
				pi_opt = piin.begin() + i * pi_in_sz;
				B_g_opt = Bcur[cur_pos];
				net_opt = netcur;
#endif // PARABAYL_BACKTRACK

			    }
			} else {
			    // here we compute F, which is min of
			    // s and all subsets of A
			    // keep in mind that F initially store s(i,A)
			    // because we recognize var we are sure
			    // that Fcur[cur_pos] and F[pos] are "aligned"
			    if (Fcur[cur_pos] < F[pos]) {
				F[pos] = Fcur[cur_pos];
#ifdef PARABAYL_BACKTRACK
				B[pos] = Bcur[cur_pos];
#endif // PARABAYL_BACKTRACK

			    }
			    pos++;
			} // if j

			cur_pos++;
		    } // if !inset
		} // for j

	    } // for i

	    // STEP 4: send F
	    // we first append Q_opt score to F
	    // if t == 0 we are at root and
	    // hence we have empty set with Q_opt = 0
	    if (t == 0) Q_opt = 0.0;
	    F.push_back(Q_opt);

//	    std::cout << "lid = " << lid << "\tQ-opt = " << Q_opt << "\tcnt = " << cnt << "\tt = " << t << std::endl;

#ifdef PARABAYL_BACKTRACK
	    // and then we create pi & net ;-)
	    pi.resize(pi_in_sz + 1);
	    memcpy(pi.begin(), pi_opt, pi_in_sz * sizeof(unsigned char));
	    pi[pi_in_sz] = g_opt;

	    net.resize(net_in_sz + 1);
	    memcpy(net.begin(), net_opt, net_in_sz * sizeof(bit_t));
	    net[net_in_sz] = B_g_opt;
#endif // PARABAYL_BACKTRACK

	    for (unsigned int i = 0; i < out_degree; ++i) {
		MPI_Send(F.begin(), F.size(), MPI_DOUBLE, out_ranks[i], PBL_F_TAG, mpi_env.comm());

#ifdef PARABAYL_BACKTRACK
		MPI_Send(pi.begin(), pi.size(), MPI_UNSIGNED_CHAR, out_ranks[i], PBL_PI_TAG, mpi_env.comm());
		MPI_Send(B.begin(), B.size(), MPI_BIT_T, out_ranks[i], PBL_B_TAG, mpi_env.comm());
		MPI_Send(net.begin(), net.size(), MPI_BIT_T, out_ranks[i], PBL_NET_TAG, mpi_env.comm());
#endif // PARABAYL_BACKTRACK

	    } // for i

	    // STEP 5: shoot yourself in the foot
	    // we have to store some Fs locally for the next level
	    // we keep track which cnt occurred in which iteration
	    // the same for pis and so on
	    Fl.push_back(F);

#ifdef PARABAYL_BACKTRACK
	    pil.push_back(pi);
	    Bl.push_back(B);
	    netl.push_back(net);
#endif // PARABAYL_BACKTRACK

	    cnt_idx[cnt_idx_pos ^ 1][cnt] = l_iter;

	    bit_t old_cnt = cnt;

	    ++l_iter;
	    ++cnt;

	    // when we hit the next level we can free some
	    // memory (i.e. remove previous level)
	    if (jaz::bit_count(old_cnt) < jaz::bit_count(static_cast<bit_t>(cnt))) {
		// remove from Fl
		unsigned int to_rm = cnt_idx[cnt_idx_pos].size() * blck_in_sz;

		for (unsigned int i = 0; i < Fl.size() - to_rm; ++i) {
		    Fl[i] = Fl[i + to_rm];
		}

		Fl.resize(Fl.size() - to_rm);

#ifdef PARABAYL_BACKTRACK
		// remove from pil
		to_rm = cnt_idx[cnt_idx_pos].size() * pi_in_sz;

		for (unsigned int i = 0; i < pil.size() - to_rm; ++i) {
		    pil[i] = pil[i + to_rm];
		}

		pil.resize(pil.size() - to_rm);

		// remove from Bl
		to_rm = cnt_idx[cnt_idx_pos].size() * B_in_sz;
		for (unsigned int i = 0; i < Bl.size() - to_rm; ++i) {
		    Bl[i] = Bl[i + to_rm];
		}

		Bl.resize(Bl.size() - to_rm);

		// remove from netl
		to_rm = cnt_idx[cnt_idx_pos].size() * net_in_sz;
		for (unsigned int i = 0; i < netl.size() - to_rm; ++i) {
		    netl[i] = netl[i + to_rm];
		}

		netl.resize(netl.size() - to_rm);
#endif // PARABAYL_BACKTRACK

		// update cnt_idx
		cnt_idx[cnt_idx_pos].clear();
		cnt_idx_pos ^= 1;
		l_iter = 0;
	    }
	} // if t

	if ((app_conf.sync_step > 0) && (t % app_conf.sync_step == 0)) {
	    MPI_Barrier(mpi_env.comm());
	}

	if ((rank == mpi_env.size() - 1) && ((t - k) % ps == 0)) {
	    std::cout << "." << std::flush;
	}
    } // for t

    delete[] in_reqs;
    delete[] in_stats;

    MPI_Barrier(mpi_env.comm());
    std::cout << master << "\ncomputing done: " << timer(TIMER_GET) << "\n" << std::endl;
    MPI_Barrier(mpi_env.comm());

    // generate output if we are the last CPU
    if (rank == mpi_env.size() - 1) {

	std::cout << "Q_opt: " << Q_opt << std::endl;

#ifdef PARABAYL_BACKTRACK
	std::cout << "\n";


	std::vector<std::string> names = data.probes;
	for (unsigned int i = 1; i < pi.size(); ++i) {
	    std::cout << i << " " << static_cast<int>(pi[i]) << ": " << names[static_cast<int>(pi[i])] << std::endl;
	}

	std::ofstream f(app_conf.output.c_str());
	if(!f) std::cout << error << "output file could not be created" << std::endl;

	for (unsigned int i = 1; i < pi.size(); ++i) {
	    if (net[i] == EMPTY_SET) {
		f << "*" << "\t" << "par" << "\t" << names[static_cast<int>(pi[i])]  << std::endl;
	    } else {
		for (unsigned int j = 0; j < n; ++j) {
		    if (inset(net[i], j)) {
			f << names[j] << "\tpar\t" << names[static_cast<int>(pi[i])]  << std::endl;
		    }
		}
	    }
	} // for i
	f.close();
#endif // PARABAYL_BACKTRACK

    } // if rank == mpi_env.size() - 1

    return true;
} // parabayl
