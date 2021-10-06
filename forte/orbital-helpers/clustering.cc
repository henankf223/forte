#include <unordered_map>
#include <vector>
#include <algorithm>

#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"

// Creating shortcut for an integer pair
// typedef  std::pair<int, int> iPair;

// Graph
struct Graph
{
    int V, E;
    std::vector<std::pair<double, std::pair<int, int>>> edges;

    // Constructor
    Graph(int V, int E)
    {
        this->V = V;
        this->E = E;
    }

    // Utility function to add an edge
    void addEdge(int u, int v, double w)
    {
        edges.push_back({w, {u, v}});
    }

    // Function to find MST using Kruskal's
    // MST algorithm
    std::unordered_map<int, std::vector<int>> kruskalMST(int n_c, std::vector<int> enforce, std::vector<int> window);
};

// Union-Find data structure for clustering
struct UnionFind
{
    int *parent, *rank;
    int n, count;

    // Constructor
    UnionFind(int n)
    {
        // Allocate memory
        this->n = n;
        this->count = n;
        parent = new int[n+1];
        rank = new int[n+1];

        // Initialize all entries
        for (int i = 0; i <= n; i++)
        {
            rank[i] = 0;
            parent[i] = i;
        }
    }

    // Find parent of u, with recurrsive path compression
    int find(int u)
    {
        if (u != parent[u])
            parent[u] = find(parent[u]);
        return parent[u];
    }

    // Union by rank
    void merge(int x, int y)
    {
        x = find(x), y = find(y);

        if (rank[x] > rank[y])
            parent[y] = x;
        else 
            parent[x] = y;

        if (rank[x] == rank[y])
            rank[y]++;
        count--;
    }

    // Check number of sets
    int count_set()
    {
        return count;
    }

    // Print all clusters
    void print_clusters()
    {
        for (int i = 0; i < n; ++i)
        {
            psi::outfile->Printf("\n Atom %d is in cluster %d. ", i, find(i));
        }
        psi::outfile->Printf("\n ");
    }
};

std::unordered_map<int, std::vector<int>> Graph::kruskalMST(int n_c, std::vector<int> enforce, std::vector<int> window)
{
    double mst_wt = 0.0;

    // Sort edges in increasing order on basis of cost
    sort(edges.begin(), edges.end());

    // Create UnionFind
    UnionFind uf(V);

    // Merge enforced elements
    if (enforce.size() > 0) {
        for(int k : enforce) {
            if (k != enforce[0]) {
                uf.merge(k, enforce[0]);
            }
        }
    }

    // Create a map to memorize whether an edge is in MST
    std::unordered_map<int, double> IsMST;

    // Iterate through all sorted edges
    std::vector< std::pair<double, std::pair<int, int>> >::iterator it;
    for (it=edges.begin(); it!=edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;

        int set_u = uf.find(u);
        int set_v = uf.find(v);

        // Check if the selected edge is creating
        // a cycle or not (Cycle is created if u
        // and v belong to same set)
        if (set_u != set_v)
        {
            // Current edge will be in the MST

            // Update MST weight
            mst_wt += it->first;

            // Merge two sets
            uf.merge(set_u, set_v);

            IsMST[u*V + v] = it->first;
        }
        int n_curr = uf.count_set();
        if (n_curr == n_c) { break; }
        // TODO: save the pointer here
    }
    //uf.print_clusters();
    std::unordered_map<int, int> index_list;
    std::unordered_map<int, std::vector<int>> cluster_list;

    for (int i = 0; i < V; ++i) {
        cluster_list[uf.find(i)].push_back(i);
    }

    if (enforce.size() > 0) { 
        index_list[0] = uf.find(enforce[0]);
        psi::outfile->Printf("\n The fragment atoms defined in input will be constrained in cluster 0. \n");
        int mask = 1;
        for (auto c : cluster_list) {
            if (c.first ==  uf.find(enforce[0])) { continue; }
            else {
                index_list[mask] = c.first;
                mask++;
            }
        }
    }
    else {
        int mask = 0;
        for (auto c : cluster_list) {
            index_list[mask] = c.first;
            mask++;
        }
    }

    std::unordered_map<int, std::vector<int>> final_cluster;

    psi::outfile->Printf("\n =============== Summary of clustering ===============");

    for (auto val : index_list) {
        int idx = val.first;
        int centroid = val.second;
        double mst_c = 0.0;
        int size_c = 0;
        int window_budget = 0;
        psi::outfile->Printf("\n Cluster %d (centroid %d):  ", idx, centroid);
        for (int at : cluster_list[centroid]) {
            psi::outfile->Printf(" %d", at);
            window_budget += window[at];
        }
        psi::outfile->Printf(" \n");
        final_cluster[idx] = cluster_list[centroid];
        for (int m : final_cluster[idx]) {
            for (int n : final_cluster[idx]) {
                int idx_c = m * V + n;
                if (IsMST.find(idx_c) != IsMST.end()) {
                    mst_c += IsMST[idx_c];
                    size_c++;
                }
            }
        }
        // TODO: if any combination of clusters have NOA smaller than max_budget, merge them
        psi::outfile->Printf("Size: %d; Single cluster inertia: %8.8f, cluster basis size: %d. \n", size_c + 1, mst_c / size_c, window_budget);
    }

    int n_cluster = uf.count_set();
    double inertia = mst_wt / (V - n_cluster);
    psi::outfile->Printf("\n The total number of clusters is %d, total inertia is %8.8f. \n ", n_cluster, inertia);
    psi::outfile->Printf("\n ======================================================");

    return final_cluster;
}


