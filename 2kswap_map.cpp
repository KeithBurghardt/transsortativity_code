/* Compile by g++ -O3 2kswap_map.cpp -o 2kswap */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>

int main() {
  int dir = 1, step = 100000;
  double rec_ass = -0.20;

  std::map<int, std::vector<int> > nn_map;
  std::string line;
  std::ifstream edgelist("../Preprocess/Sets24/24sn25");
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
  std::cout << "Network assortativity: " << r << std::endl;

  double cond1, cond2, back = 0.5;
  int v0, w0, v1, w1, status = 0, count = 0, swap = 0;
  bool ext_edge1, ext_edge2;
  while (true) {
  // for (int z = 0; z < 10; ++z) {
    for (int i=0; i<step; i++) {
      while (true) {
        ext_edge1 = false, ext_edge2 = false;
        do v0 = rand()%n;
        while (k[v0] == 0);
        v1 = nn_seq[v0][rand()%k[v0]];
        do w0 = rand()%n;
        while (k[w0] == 0);
        w1 = nn_seq[w0][rand()%k[w0]];
        if (std::find(nn_seq[v0].begin(), nn_seq[v0].end(), w1) < nn_seq[v0].end()) ext_edge1 = true;
        if (std::find(nn_seq[v1].begin(), nn_seq[v1].end(), w0) < nn_seq[v1].end()) ext_edge1 = true;
        if (std::find(nn_seq[v0].begin(), nn_seq[v0].end(), w0) < nn_seq[v0].end()) ext_edge2 = true;
        if (std::find(nn_seq[v1].begin(), nn_seq[v1].end(), w1) < nn_seq[v1].end()) ext_edge2 = true;
        if (v0 != w0 && v0 != w1 && v1 != w0 && v1 != w1 && !(ext_edge1 && ext_edge2)) break;
      }

      cond1 = k[v0]*k[w1] + k[v1]*k[w0] - k[v0]*k[v1] - k[w0]*k[w1];
      cond2 = k[v0]*k[w0] + k[v1]*k[w1] - k[v0]*k[v1] - k[w0]*k[w1];
      if (!ext_edge1) 
        if (cond1*dir > 0) status = 1;
        else if ((double) rand()/RAND_MAX < back) {
          status = 1;
          if (back > 0) back -= (double) 1/e/100;
        }
      if (!ext_edge2)
        if (cond2*dir > 0) status = 2;
        else if ((double) rand()/RAND_MAX < back) {
          status = 2;
          if (back > 0) back -= (double) 1/e/100;
        }

      switch (status) {
        case 1:
          c += (double) cond1/e;
          std::iter_swap(std::find(nn_seq[v0].begin(), nn_seq[v0].end(), v1), std::find(nn_seq[w0].begin(), nn_seq[w0].end(), w1));
          std::iter_swap(std::find(nn_seq[v1].begin(), nn_seq[v1].end(), v0), std::find(nn_seq[w1].begin(), nn_seq[w1].end(), w0));
          swap++;
          status = 0;
          break;
        case 2:
          c += (double) cond2/e;
          std::iter_swap(std::find(nn_seq[v0].begin(), nn_seq[v0].end(), v1), std::find(nn_seq[w1].begin(), nn_seq[w1].end(), w0));
          std::iter_swap(std::find(nn_seq[v1].begin(), nn_seq[v1].end(), v0), std::find(nn_seq[w0].begin(), nn_seq[w0].end(), w1));
          swap++;
          status = 0;
          break;
        default:
          break;
      }

      r = c/s;
      if ((r-rec_ass)*dir > 0) {
        std::string filename = "rdmswap_" + std::to_string(count++);
        std::ofstream outfile;
        outfile.open(filename); /* Output file */
        outfile << "# Nodes: " << n << ", # Edges: " << e << ", Assortativity: " << r << std::endl;
        for (int i=0; i<n; i++) {
          std::sort (nn_seq[i].begin(), nn_seq[i].end());
          for (std::vector<int>::iterator it=nn_seq[i].begin(); it!=nn_seq[i].end(); ++it)
            if (*it>i) outfile << i << '\t' << *it << std::endl;
        }
        outfile.close();
        rec_ass += 0.05*dir;
      }
    } /* end for */
    std::cout << "Assort.: " << r << ", Swap: " << swap << ", Back: " << back << std::endl;
    swap = 0;
  } /* end while */

  return 0;
}
