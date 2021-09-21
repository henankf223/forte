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
typedef  std::pair<int, int> iPair;

// Structure to represent a graph
struct Graph
{
    int V, E;
    std::vector<std::pair<double, iPair>> edges;

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
    double kruskalMST(int n_c);
};

// To represent Disjoint Sets
struct DisjointSets
{
    int *parent, *rnk;
    int n;

    // Constructor.
    DisjointSets(int n)
    {
        // Allocate memory
        this->n = n;
        parent = new int[n+1];
        rnk = new int[n+1];

        // Initially, all vertices are in
        // different sets and have rank 0.
        for (int i = 0; i <= n; i++)
        {
            rnk[i] = 0;

            //every element is parent of itself
            parent[i] = i;
        }
    }

    // Find the parent of a node 'u'
    // Path Compression
    int find(int u)
    {
        /* Make the parent of the nodes in the path
           from u--> parent[u] point to parent[u] */
        if (u != parent[u])
            parent[u] = find(parent[u]);
        return parent[u];
    }

    // Union by rank
    void merge(int x, int y)
    {
        x = find(x), y = find(y);

        /* Make tree with smaller height
           a subtree of the other tree  */
        if (rnk[x] > rnk[y])
            parent[y] = x;
        else // If rnk[x] <= rnk[y]
            parent[x] = y;

        if (rnk[x] == rnk[y])
            rnk[y]++;
    }

    // Check number of sets
    int count_set()
    {
        std::unordered_map<int, bool> clusters;
        for (int i = 0; i < n; ++i)
        {
            clusters[find(i)] = true;
        }
        return clusters.size();
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

double Graph::kruskalMST(int n_c)
{
    double mst_wt = 0.0; // Initialize result

    // Sort edges in increasing order on basis of cost
    sort(edges.begin(), edges.end());

    // Create disjoint sets
    DisjointSets ds(V);

    // Iterate through all sorted edges
    std::vector< std::pair<double, iPair> >::iterator it;
    for (it=edges.begin(); it!=edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;

        int set_u = ds.find(u);
        int set_v = ds.find(v);

        // Check if the selected edge is creating
        // a cycle or not (Cycle is created if u
        // and v belong to same set)
        if (set_u != set_v)
        {
            // Current edge will be in the MST
            // so print it
            //cout << u << " - " << v << endl;

            // Update MST weight
            mst_wt += it->first;

            // Merge two sets
            ds.merge(set_u, set_v);
        }
        int n_curr = ds.count_set();
        if (n_curr == n_c) { break; }
    }
    ds.print_clusters();
    psi::outfile->Printf("\n The total inertia is %8.8f. \n ", mst_wt);

    return mst_wt;
}


