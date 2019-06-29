/* Compile by g++ -O3 fspdxavg.cpp -o fspdx */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cmath>

int main() {
  auto tstart = std::chrono::system_clock::now();
  std::map<int, std::vector<int> > nn_map;
  std::string line;
  std::ifstream edgelist("./Google2");
  while (std::getline(edgelist,line)) {
    int i, j;
    std::istringstream iss(line);
    if (line[0] != '#' && (iss >> i >> j)) {
      nn_map[i].push_back(j);
      nn_map[j].push_back(i);
    }
  }
  edgelist.close();
  std::ofstream outfile;
  outfile.open("fspdx.out"); /* Output file */

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

  outfile << "No. nodes: " << n << ", No. edges: " << e << ", Assortativity: " << r << std::endl;

  /* Friendship Paradox */
  std::vector<std::vector<int> > nk_seq(n);
  std::vector<double> nn_median(n);
  std::vector<int> paradox_class(k_max+1);
  for (int i=0; i<n; i++) for (int j=0; j<k[i]; j++) nk_seq[i].push_back(k[nn_seq[i][j]]);
  for (int i=0; i<n; i++) if (k[i]) {
	std::sort(nk_seq[i].begin(), nk_seq[i].end());
    if (k[i]%2) nn_median[i] = (double) nk_seq[i][(k[i]-1)/2];
    else nn_median[i] = (double) (nk_seq[i][k[i]/2-1] + nk_seq[i][k[i]/2])/2.0;
  }
  for (int i=0; i<n; i++) if (nn_median[i]>k[i]) paradox_class[k[i]]++;

  outfile << "[";
  for (int i=0; i<k_max+1; i++) if (k_hist[i]) {
    outfile << "[" << i << ", " << k_hist[i] << ", " << paradox_class[i] << "]";
    if (i < k_max) outfile << ", ";
  }
  outfile << "]\n" << std::endl;

  /* Neighbor degree mean & variance */
  std::vector<std::vector<double> > P_ff(n,std::vector<double>(4));
  for (int i=0; i<n; i++) {
    double nn_avg = 0.0, nn_avg_sq = 0.0;
    for (int j=0; j<k[i]; j++) {
      nn_avg += (double) k[nn_seq[i][j]] / k[i];
      nn_avg_sq += (double) k[nn_seq[i][j]] / k[i] * k[nn_seq[i][j]];
    }
    P_ff[k[i]][0] += (double) nn_avg / k_hist[k[i]];
    P_ff[k[i]][1] += (double) nn_avg * nn_avg / k_hist[k[i]];
    P_ff[k[i]][2] += (double) nn_avg_sq / k_hist[k[i]];
  }
  for (int i=0; i<k_max+1; i++) P_ff[i][3] = std::sqrt(P_ff[i][1] - P_ff[i][0] * P_ff[i][0]);

/*
  outfile << "[";
  for (int i=0; i<k_max+1; i++) if (k_hist[i]) {
    outfile << "[" << i << ", " << k_hist[i] << "]";
    if (i < k_max) outfile << ", ";
  }
  outfile << "]\n" << std::endl;
*/

  /* Output Neighbor correlation */
  outfile << "[";
  for (int i=0; i<k_max+1; i++) if (k_hist[i]) {
    outfile << "[" << i << ", " << std::scientific << P_ff[i][0] << ", " << P_ff[i][1] << ", " << P_ff[i][2] << "]";
    if (i < k_max) outfile << ", ";
  }
  outfile << "]\n" << std::endl;

  outfile << "[";
  for (int i=0; i<k_max+1; i++) if (k_hist[i]) {
    outfile << "[" << i << ", " << P_ff[i][3]*P_ff[i][3]/2.0/(P_ff[i][0]-i)/(P_ff[i][0]-i) << ", " << (double) paradox_class[i]/k_hist[i] << "]";
    if (i < k_max) outfile << ", ";
  }
  outfile << "]"<< std::endl;
  outfile.close();

  std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - tstart;
  std::cout << "Time used: " <<  elapsed.count() << " sec" << std::endl;
  return 0;
}
