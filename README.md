## Code for measuring assortativity and transsortativity and algorithms to change assortativity and transsortativity in arbitrary networks
======

This is the code used to:
-  calculate transsortativity
-  increase and decrease transsortativity in any network
-  calculate the freindship paradox in a network
-  calculate the largest vulnerable component

All code developed by Xin-Zeng Wu (2019)
email: energiya@gmail.com

- 2k_swap_map.cpp  :  Performs 2K-rewiring on the input network (as edgelist) to achieve higher or lower assortativity based on Metropolis alg proposed by Mark Newman.

		      Input: 
		      
		      	- edgelist (change edgelist_file string)
			
		      Output: 
		      
		        - new edgelist files with correct transsortativity value
			
		        - output files: "rdmswap_" + count

- 3k_swap_map.cpp  :  Performs 3K-rewiring on the input network (as edgelist) to achieve higher or lower neoghbor assortativity based on greedy algorithm.

		      Input: 
		      
		      	- edgelist (change edgelist_file string)
			
		      Output: 
		      
		        - new edgelist files with correct transsortativity value
			
		        - output files: "rdmswap_" + count

- fspdxavg.cpp  :     Calculate the transsortativity of the input network (as edgelist)
		      Input:
		      
		  	- edgelist (change input_file string)
			
		      Output:
		      
		        - output as [k, r_{3k}(k), p(k)]. 
			
                        - output_file: "fspdx.out"
			
		        - mean transsortativity can be calculated using the output. 
			
	 	       (Computations of Figure 2 is performed on this code)

- vulnerable.py  :  Find the number of nodes and edges of the giant vulnerable component in the network using the Watts cascade model
		      Input:
		      
		  	- edgelist (change infile string)
			
		      Output:
		      
			- [1/phi, n_gvul, e_gvul]
			
			- (iterates sequentially through range of values of phi)
			
			- outfile: './vulnerable.out'
			
		        (Simulations of Figure 4 is performed on this code)

- zeroone_opt.c  :  Calculate the strength of the majority illusion on a series of assignments of attribute x=1 with increasing k-x correlation, and compare with theoretical results based on Eq (2.7). 

                    Input: 
		    
		        - edgelist (change input_file string)
			
		     Output:
		     
		        - Increase steps of corr, output in one file: [corr, paradox, P_maj]
			
			- corr: correlation between attribute and degree
			
			- P_maj: Theoretical majority illusion calculation
			
			- paradox: Number of agents who experience a majority illusion
			
			- output file ("of" in the code)
			
			(Simulations of Figure 9 is performed on this code)
