/*
 * @BEGIN LICENSE
 *
 * Forte: an open-source plugin to Psi4 (https://github.com/psi4/psi4)
 * that implements a variety of quantum chemistry methods for strongly
 * correlated electrons.
 *
 * Copyright (c) 2012-2021 by its authors (see COPYING, COPYING.LESSER, AUTHORS).
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 * @END LICENSE
 */
#include <map>
#include <numeric>
#include <regex>
#include <vector>
#include <cmath>

#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"

#include "base_classes/forte_options.h"

#include "fragment_projector.h"
#include "clustering.cc"

using namespace psi;

namespace forte {

std::pair<psi::SharedMatrix, int> make_fragment_projector(SharedWavefunction wfn,
                                                          std::shared_ptr<ForteOptions> options) {
    // Run this code only if user specified fragments
    std::shared_ptr<Molecule> molecule = wfn->molecule();
    int nfrag = molecule->nfragments();
    if (nfrag == 1) {
        throw PSIEXCEPTION("A input molecule with fragments (-- in atom list) is required "
                           "for embedding!");
    }
    outfile->Printf(
        "\n  The input molecule has %d fragments, treating the first fragment as the system.\n",
        nfrag);

    std::shared_ptr<BasisSet> prime_basis = wfn->basisset();

    // Create a fragmentprojector object
    FragmentProjector FP(molecule, prime_basis);

    // Create a fragmentprojector with the second constructor if we want to project to minAO or use
    // IAO procedure FragmentProjector FP(molecule, prime_basis, minao_basis);

    // Compute and return the projector matrix
    psi::SharedMatrix Pf = FP.build_f_projector(prime_basis, options);

    if (options->get_bool("AUTOFRAG")) {
        // Apply automatic fragmentation
        FP.build_auto_projector(wfn->Fa(), options);
        throw PSIEXCEPTION("The automatic fragmentation implementation is not done :( ...");
    }

    int nbfA = FP.get_nbf_A();
    std::pair<psi::SharedMatrix, int> Projector = std::make_pair(Pf, nbfA);

    return Projector;
}

FragmentProjector::FragmentProjector(std::shared_ptr<Molecule> molecule,
                                     std::shared_ptr<BasisSet> basis)
    : molecule_(molecule), basis_(basis) {
    startup();
}

void FragmentProjector::startup() {

    std::vector<int> sys_list = {0};  // the first fragment in the input is the system
    std::vector<int> ghost_list = {}; // leave empty to include no fragments with ghost atoms

    // extract the sys molecule objects
    std::shared_ptr<Molecule> mol_sys = molecule_->extract_subsets(sys_list, ghost_list);

    outfile->Printf("\n  System Fragment \n");
    mol_sys->print();

    nbf_ = basis_->nbf();
    outfile->Printf("\n  Number of basis on all atoms: %d", nbf_);

    natom_A_ = mol_sys->natom();
    int count_basis = 0;
    for (int mu = 0; mu < nbf_; mu++) {
        int A = basis_->function_to_center(mu);
        // outfile->Printf("\n  Function %d is on atom %d", mu, A);
        if (A < natom_A_) {
            count_basis += 1;
        }
    }
    outfile->Printf("\n  Number of basis in the system fragment: %d", count_basis);
    nbf_A_ = count_basis;
}

SharedMatrix FragmentProjector::build_f_projector(std::shared_ptr<psi::BasisSet> basis, std::shared_ptr<ForteOptions> options) {

    std::shared_ptr<IntegralFactory> integral_pp(new IntegralFactory(basis, basis, basis, basis));
    std::shared_ptr<OneBodyAOInt> S_int(integral_pp->ao_overlap());
    S_ = std::make_shared<psi::Matrix>("S_nn", nbf_, nbf_);
    S_int->compute(S_);

    std::vector<int> add_basis = options->get_int_list("EMBEDDING_CUSTOM_PARTITION");
    std::vector<int> basis_vec;

    if (options->get_bool("PURE_BASIS_PARTITION")) {
        outfile->Printf("\n Ignore the fragment and use pure basis set input as A. \n");
        nbf_A_ = 0;
    }

    for (int i = 0; i < nbf_A_; ++i) {
        basis_vec.push_back(i);
    }
    outfile->Printf("\n Will add %d basis to the fragment.", add_basis.size());
    for (int v : add_basis) {
        if (v < nbf_A_ || v >= nbf_) {
            outfile->Printf("\n Warning! Illegal basis number (%d) ! Check the input.", v);
        }
        else {
            outfile->Printf("\n Add basis %d to the fragment projector", v);
            basis_vec.push_back(v);
        }
    }
    nbf_A_ = basis_vec.size();

    if (nbf_A_ == 0) {
        throw PSIEXCEPTION("The fragment size is 0! Use ADD_FRAG_BASIS to assign fragment basis, or turn off PURE_BASIS_PARTITION.");
    }

    int nbf_new = basis_vec.size();
    // Construct the system portion of S (S_A)
    //SharedMatrix S_A = S_->get_block(fragA, fragA);
    SharedMatrix S_A(new Matrix("S_A block", nbf_new, nbf_new));
    for (int m = 0; m < nbf_new; ++m) {
        for (int n = 0; n < nbf_new; ++n) {
            S_A->set(m, n, S_->get(basis_vec[m], basis_vec[n]));
        }
    }

    // Construct S_A^-1 and store it in a matrix of size nbf x nbf
    S_A->general_invert();

    SharedMatrix S_A_nn(new Matrix("S system in fullsize", nbf_, nbf_));
    //S_A_nn->set_block(fragA, fragA, S_A);
    for (int p = 0; p < nbf_new; ++p) {
        for (int q = 0; q < nbf_new; ++q) {
            S_A_nn->set(basis_vec[p], basis_vec[q], S_A->get(p, q));
        }
    }

    // Evaluate AO basis projector  P = S^T (S_A)^{-1} S
    S_A_nn->transform(S_);

    return S_A_nn;
}

void FragmentProjector::build_auto_projector(SharedMatrix F_w, std::shared_ptr<ForteOptions> options) {
    auto dist_type = options->get_str("AUTOFRAG_DIST");
    double F_weight = options->get_double("FOCK_WEIGHT");
    double M_weight = options->get_double("GEOMETRY_WEIGHT");
    int num_clusters = options->get_int("N_FRAGMENT");
    bool enforce_fragment = options->get_bool("FRAGMENT_CONSTRAINED");
    bool use_custom = options->get_bool("USE_CUSTOM_WINDOW");

    outfile->Printf("\n  Computing full overlap distances from S: \n");

    std::vector<int> window_list;
    int atom_idx = 0;
    if (!use_custom) {
        outfile->Printf("\n  The partition will be based on basis located on every atoms. \n");
        int atom_window_count = 0;
        for (int ka = 0; ka < nbf_; ka++) {
            int A = basis_->function_to_center(ka);
            if (A == atom_idx) {
                atom_window_count += 1;
            }
            else {
                window_list.push_back(atom_window_count);
                atom_idx += 1;
                atom_window_count = 1;
            }
        }
        window_list.push_back(atom_window_count);
    }
    else {
        outfile->Printf("\n  The partition will be based on input basis indices from CUSTOM_FRAG_WINDOW. \n");
        window_list = options->get_int_list("CUSTOM_FRAG_WINDOW");
        outfile->Printf("\n  The window size is: %d \n", window_list.size());
        outfile->Printf("\n  The first window is: %d \n", window_list[0]);
    }

    outfile->Printf("\n  The atom/composite windows are: \n");
    atom_idx = 0;
    outfile->Printf("\n  The window count is: %d \n", window_list.size());
    for (auto val : window_list) {
        outfile->Printf("\n Atom/composite %d have %d AO on it.", atom_idx, val);
        atom_idx += 1;
    }
    outfile->Printf("\n");

    std::vector<int> zeropi(1, 0);
    Dimension S1_begin(zeropi);
    Dimension S1_end(zeropi);
    Dimension S2_begin(zeropi);
    Dimension S2_end(zeropi);

    int atom_idx_1 = 0;
    int atom_idx_2 = 0;

    int cum1 = 0;
    int cum2 = 0;

    int natom_full = window_list.size();
    //SharedMatrix dist_mat(new Matrix("Atomic distance matrix", natom_full, natom_full));

    Graph g_dist(natom_full, natom_full*natom_full);
    //SharedMatrix mol_dist_mat = std::make_shared<psi::Matrix>(molecule_->distance_matrix());

    for (auto val1 : window_list) {
        S1_begin[0] = cum1;
        S1_end[0] = cum1 + val1;
        Slice S1(S1_begin, S1_end);
        for (auto val2 : window_list) {
            S2_begin[0] = cum2;
            S2_end[0] = cum2 + val2;
            Slice S2(S2_begin, S2_end);

            SharedMatrix S_12 = S_->get_block(S1, S2);
            SharedMatrix F_12 = F_w->get_block(S1, S2);
            F_12->scale(F_weight);
            S_12->add(F_12);

            if (atom_idx_1 != atom_idx_2) {
                double Dist = S_12->trace();
                //double real_dist = mol_dist_mat->get(atom_idx_1, atom_idx_2) * M_weight;
                double real_dist = 0.0;
                if (dist_type == "TR") {
                    //outfile->Printf("\n  The trace distances between atom %d and %d is %8.8f: \n", atom_idx_1, atom_idx_2, Dist);
                    //dist_mat->set(atom_idx_1, atom_idx_2, Dist);
                    g_dist.addEdge(atom_idx_1, atom_idx_2, -1.0 * log(abs(Dist)) + real_dist);
                }
                double Dist_avg = Dist / std::min(val1, val2);
                if (dist_type == "TR_AVG") {
                    //outfile->Printf("\n  The average trace distances between atom %d and %d is %8.8f: \n", atom_idx_1, atom_idx_2, Dist_avg); 
                    //dist_mat->set(atom_idx_1, atom_idx_2, Dist_avg);
                    g_dist.addEdge(atom_idx_1, atom_idx_2, -1.0 * log(abs(Dist_avg)) + real_dist);
                }
                double Dist_sum_sq = S_12->sum_of_squares();
                if (dist_type == "SSQ") {
                    //outfile->Printf("\n  The SSQ distances between atom %d and %d is %8.8f: \n", atom_idx_1, atom_idx_2, Dist_sum_sq);
                    //dist_mat->set(atom_idx_1, atom_idx_2, Dist_sum_sq);
                    g_dist.addEdge(atom_idx_1, atom_idx_2, -1.0 * log(sqrt(abs(Dist_sum_sq))) + real_dist);
                }
		double Dist_sum_sq_avg = Dist_sum_sq/(val1 * val2);
		if (dist_type == "SSQ_AVG") {
                    //outfile->Printf("\n  The size-weighted SSQ distances between atom %d and %d is %8.8f: \n", atom_idx_1, atom_idx_2, Dist_sum_sq_avg);
                    //dist_mat->set(atom_idx_1, atom_idx_2, Dist_sum_sq_avg);
                    g_dist.addEdge(atom_idx_1, atom_idx_2, -1.0 * log(sqrt(abs(Dist_sum_sq_avg))) + real_dist);
                }
            }

            cum2 += val2;
            atom_idx_2 += 1;
        }
        cum1 += val1;
        cum2 = 0;
        atom_idx_2 = 0;
        atom_idx_1 += 1;
    }

    //outfile->Printf("\n The weighted distance matrix: \n");
    //dist_mat->print();

    outfile->Printf("\n Perform single-linkage HAC clustering ...");
    std::vector<int> enforce_list;
    if (enforce_fragment) {
        outfile->Printf("\n Linking %d atoms in the original fragment ...", natom_A_);
        for (int i = 0; i < natom_A_; ++i) { 
            enforce_list.push_back(i);
            outfile->Printf("\n Adding atom %d to fix list", i);
        }
    }
    std::unordered_map<int, std::vector<int>> res = g_dist.kruskalMST(num_clusters, enforce_list, window_list);
    outfile->Printf("\n Done ... \n");

    return;
}

} // namespace forte
