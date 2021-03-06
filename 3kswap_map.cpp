/* Compile by g++ -O3 3kswap_map.cpp -o 3kswap */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>

int main()
{
  string edgelist_file = "./Sets21o5/21o5sp30_c45"
  int step = 100000, dir_rho = 1;
  double rec_rho = -0.05;

  std::map<int, std::vector<int> > nn_map;
  std::string line;
  std::ifstream edgelist(edgelist_file);
  while (std::getline(edgelist,line)) {
    int i, j;
    std::istringstream iss(line);
    if (line[0] != '#' && (iss >> i >> j)) {
      nn_map[i].push_back(j);
      nn_map[j].push_back(i);
    }
  }
  edgelist.close();

  int n = nn_map.rbegin()->first + 1;
  std::vector<std::vector<int> > nn_seq(n);
  for (int i=0; i<n; i++) nn_seq[i].swap(nn_map[i]);
  std::vector<int> k(n);
  for (int i=0; i<n; i++) k[i] = nn_seq[i].size();;
  int e = std::accumulate(k.begin(), k.end(), 0)/2;
  int k_max = *std::max_element(k.begin(), k.end());
  std::vector<std::vector<int> > k_class(k_max+1);
  for (int i=0; i<n; i++) k_class[k[i]].push_back(i);
  std::cout << "No. nodes: " << n << ", No. edges: " << e << ", Max degree: " << k_max << std::endl;

  /* Assortativity */
  double c = 0.0, s = 0.0, avg_q = 0.0;
  std::vector<int> k_hist(k_max+1);
  std::for_each(k.begin(), k.end(), [&] (int i) { k_hist[i]++; });
  std::vector<double> Q(k_max+1);
  for (int i=0; i<k_max+1; i++) Q[i] = (double) i*k_hist[i]/e/2;
  for (int i=0; i<k_max+1; i++) avg_q += i*Q[i];
  for (int i=0; i<k_max+1; i++) s += (i-avg_q)*(i-avg_q)*Q[i];
  for (int i=0; i<n; i++) for (int j=0; j<k[i]; j++) c += (k[i]-avg_q)*(k[nn_seq[i][j]]-avg_q)/e/2;
  double r = c/s;
  std::cout << "Transsortativity: " << r << std::endl;

  std::vector<double> nn_avg(n);
  std::vector<double> nn_avg_sq(n);
  for (int i=0; i<n; i++) for (int j=0; j<k[i]; j++) {
    nn_avg[i] += (double) k[nn_seq[i][j]] / k[i];
    nn_avg_sq[i] += (double) k[nn_seq[i][j]] / k[i] * k[nn_seq[i][j]];
  }
  std::vector<std::vector<double> > P_ff(k_max+1, std::vector<double> (3, 0));
  for (int i=0; i<n; i++) {
    P_ff[k[i]][0] += (double) nn_avg[i] / k_hist[k[i]];
    P_ff[k[i]][1] += (double) nn_avg[i] * nn_avg[i] / k_hist[k[i]]; /* 3K structure */
    P_ff[k[i]][2] += (double) nn_avg_sq[i] / k_hist[k[i]]; /* 2K structure */
  }
  std::vector<double> rho_k(k_max+1);
  for (int i=2; i<k_max+1; i++) if (k_hist[i]) 
    rho_k[i] = (i * (P_ff[i][1] - P_ff[i][0] * P_ff[i][0]) / (P_ff[i][2] - P_ff[i][0] * P_ff[i][0]) - 1) / (i-1); /* Corr Coef */
  double rho_global = 0.0;
  for (int i=2; i<k_max+1; i++) rho_global += rho_k[i] * (double) k_hist[i] / n;
  std::cout << "Mean Transsortativity: " << rho_global << std::endl;

  double back = 0.5;
  int count = 0, swap = 0;
  bool status = false;
  while (true) {
    int v0, w0, v1, w1, kc, x;
    bool ext_edge;
    for (int i=0; i<step; i++) {

        while (true) {
            
            // a random node is chosen with degree > 0
            
            do v0 = rand()%n;
            while (k[v0] == 0);
            
            // the node's degree is recorded
            
            kc = k[v0];
            //a random node is selected
            v1 = nn_seq[v0][rand()%k[v0]];
            // the random edge has end nodes v0, v1
            
            // w0 is chosen as a random nodes with the same degree as v0
            w0 = k_class[k[v0]][rand()%k_hist[k[v0]]]; // let k[w0] == k[v0]
            // w1 is a randomly chosen node
            w1 = nn_seq[w0][rand()%k[w0]];
            // a new edge is w0-w1 where degree of w0 = w1 but degree of w1 is not necessarily v1's degree
            
            // check if external edge exists
            if (std::find(nn_seq[v0].begin(), nn_seq[v0].end(), w1) < nn_seq[v0].end()) ext_edge = true;
            if (std::find(nn_seq[v1].begin(), nn_seq[v1].end(), w0) < nn_seq[v1].end()) ext_edge = true;
            // stop searching when: edges are not the same, and they are not external edges
            if (v0 != w0 && v0 != w1 && v1 != w0 && v1 != w1 && !ext_edge) break;
        }
        
        // conditions for us to move transsortativity in the right direction:
        // if degree of v1 > w1, mean degree of neighbors of v0 < w0
        // OR degree of v1 < w1, mean degree of neighbors of v0 > w0 and we want to increase transsortativity
        // if degree of v1 > w1, mean degree of neighbors of v0 > w0
        // OR degree of v1 < w1, mean degree of neighbors of v0 < w0 and we want to decrease transsortativity
        if ((k[v1] > k[w1]) != (nn_avg[v0] > nn_avg[w0]) == (dir_rho > 0)) status = true;
        //else with p(back) = 0.5, back -> back -1/(sum of all edges);
        else if ((double) rand()/RAND_MAX < back) {
            status = true;
            if (back > 0) back -= (double) 1/e;
        }
        // if swapping moves transsortativity in the right direction...
        if (status) {
            // update statistics
            // P_ff (degree of v0) -= (<q(v0)>^2 + <q(w0)>^2)/# nodes with degree equal to v0
            P_ff[kc][1] -= (nn_avg[v0] * nn_avg[v0] + nn_avg[w0] * nn_avg[w0]) / k_hist[kc];
            // we record te new nn_av of v0, w0
            nn_avg[v0] += (double) (k[w1] - k[v1]) / kc;
            nn_avg[w0] += (double) (k[v1] - k[w1]) / kc;
            // P_ff (degree of v0) += (<q(v0)>^2 + <q(w0)>^2)/# nodes with degree equal to v0
            P_ff[kc][1] += (nn_avg[v0] * nn_avg[v0] + nn_avg[w0] * nn_avg[w0]) / k_hist[kc];
            // transsortativity for degree equal to v0
            rho_k[kc] = (kc * (P_ff[kc][1] - P_ff[kc][0] * P_ff[kc][0]) / (P_ff[kc][2] - P_ff[kc][0] * P_ff[kc][0]) - 1) / (kc-1);
            // SWAP
            *std::find(nn_seq[v0].begin(), nn_seq[v0].end(), v1) = w1;
            *std::find(nn_seq[v1].begin(), nn_seq[v1].end(), v0) = w0;
            *std::find(nn_seq[w0].begin(), nn_seq[w0].end(), w1) = v1;
            *std::find(nn_seq[w1].begin(), nn_seq[w1].end(), w0) = v0;
            swap++;
            status = false;
        }
        // update transsortativity
        rho_global = 0.0;
        for (int i=2; i<k_max+1; i++) rho_global += rho_k[i] * (double) k_hist[i] / n;
        // if transsortativity is decreasing when it should be, or increasing when it should be
        if ((rho_global - rec_rho) * dir_rho > 0) {
            // create file name
            std::string filename = "rdmswap_" + std::to_string(count++);
            std::ofstream outfile;
            // output statistics to this file name
            outfile.open(filename); /* Output file */
            // output:
            //    - # nodes
            //    - # edges
            //    - assortativity
            //    - transsortativity (called "neighbor-neighbor correlation" in the code)
            outfile << "# Nodes: " << n << ", # Edges: " << e << ", Assortativity: " << r << ", Nb-Nb corr: " << rho_global << std::endl;
            // over all nodes
            for (int i=0; i<n; i++) {
                // sort sequence of nearest neighbors
                std::sort (nn_seq[i].begin(), nn_seq[i].end());
                for (std::vector<int>::iterator it=nn_seq[i].begin(); it!=nn_seq[i].end(); ++it)
                    // if i < j print (i,j) pair
                    if (*it>i) outfile << i << '\t' << *it << std::endl;
            }
            outfile.close();
            // record direction transsortativity changed since last time: +0.05, 0, -0.05
            rec_rho += 0.05*dir_rho;
        }
        
    }
    std::cout << "Assort.: " << r << ", Transsortativity: " << rho_global << ", Swap: " << swap << ", Back: " << back << std::endl;
    swap = 0;
  }

  return 0;
}





// original, uncommented code, starting at line 80 - 163


/*
 while (true) {
 ext_edge = false;
 do v0 = rand()%n;
 while (k[v0] == 0);
 kc = k[v0];
 v1 = nn_seq[v0][rand()%k[v0]];
 w0 = k_class[k[v0]][rand()%k_hist[k[v0]]]; // let k[w0] == k[v0]
 w1 = nn_seq[w0][rand()%k[w0]];
 if (std::find(nn_seq[v0].begin(), nn_seq[v0].end(), w1) < nn_seq[v0].end()) ext_edge = true;
 if (std::find(nn_seq[v1].begin(), nn_seq[v1].end(), w0) < nn_seq[v1].end()) ext_edge = true;
 if (v0 != w0 && v0 != w1 && v1 != w0 && v1 != w1 && !ext_edge) break;
 }
 
 if ((k[v1] > k[w1]) != (nn_avg[v0] > nn_avg[w0]) == (dir_rho > 0)) status = true;
 else if ((double) rand()/RAND_MAX < back) {
 status = true;
 if (back > 0) back -= (double) 1/e;
 }
 
 if (status) {
 P_ff[kc][1] -= (nn_avg[v0] * nn_avg[v0] + nn_avg[w0] * nn_avg[w0]) / k_hist[kc];
 nn_avg[v0] += (double) (k[w1] - k[v1]) / kc;
 nn_avg[w0] += (double) (k[v1] - k[w1]) / kc;
 P_ff[kc][1] += (nn_avg[v0] * nn_avg[v0] + nn_avg[w0] * nn_avg[w0]) / k_hist[kc];
 rho_k[kc] = (kc * (P_ff[kc][1] - P_ff[kc][0] * P_ff[kc][0]) / (P_ff[kc][2] - P_ff[kc][0] * P_ff[kc][0]) - 1) / (kc-1);
 *std::find(nn_seq[v0].begin(), nn_seq[v0].end(), v1) = w1;
 *std::find(nn_seq[v1].begin(), nn_seq[v1].end(), v0) = w0;
 *std::find(nn_seq[w0].begin(), nn_seq[w0].end(), w1) = v1;
 *std::find(nn_seq[w1].begin(), nn_seq[w1].end(), w0) = v0;
 swap++;
 status = false;
 }
 
 rho_global = 0.0;
 for (int i=2; i<k_max+1; i++) rho_global += rho_k[i] * (double) k_hist[i] / n;
 if ((rho_global - rec_rho) * dir_rho > 0) {
 std::string filename = "rdmswap_" + std::to_string(count++);
 std::ofstream outfile;
 outfile.open(filename); //Output file
 outfile << "# Nodes: " << n << ", # Edges: " << e << ", Assortativity: " << r << ", Nb-Nb corr: " << rho_global << std::endl;
 for (int i=0; i<n; i++) {
 std::sort (nn_seq[i].begin(), nn_seq[i].end());
 for (std::vector<int>::iterator it=nn_seq[i].begin(); it!=nn_seq[i].end(); ++it)
 if (*it>i) outfile << i << '\t' << *it << std::endl;
 }
 outfile.close();
 rec_rho += 0.05*dir_rho;
 }
 */

