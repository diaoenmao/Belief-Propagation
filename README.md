# Belief-Propagation
Implementation of Generalized Belief Propagation and Convergence Rate Analysis


Please add folder ./topology analysis to Matlab path before run any test.

Files taken from Brown University CS 254:
add_facnode.m
add_varnode.m
belief_diff.m
get_beliefs.m
init_graph.m
marg_brute_force.m

topology analysis folder contains Matlab Tools for Network Analysis

Files I contributed:

Implementation:
adjGfacG.m %tranforms adjacency matrix into factor graphs
ave_scatter.m %average scatter plot vertically
cardi.m %retrieve cardinality of a factor graph
initialize.m %initialize incoming and outgoing message cell arrays
run_loopy_bp_parallel.m %the actual BP algorithm

Simulation:
test.m % test the rate of convergence wrt. different parameters
testdim.m % test the rate of convergence wrt. dimension of variable nodes
testvis.m % visualize convergence process

Data:
image.mat % Full workspace after run test.m
Stat.mat % Full output of test.m and we can just use this file to retrieve the plots in the report.