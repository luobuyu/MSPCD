The zip-file Bench.zip contains 180 completely random task graphs. 
The names of task graph files are t<n>_<r>_<i>.td, where <n> should be 
substituted with the number of task in that file, <r> with edge existence 
density, while $i$ is the index of graph with same $n$ and $r$ value (there
are 6 graphs for each ($n$, $r$) pair). The graphs in that
zip-file does not depend on multiprocessor architecture, they can
be scheduled to arbitrary multiprocessor system, but, since
optimal solution is not known, the solution quality can not be
estimated.


In each file, data are written in following format:

- first row contains the value for the number of tasks n;

- next n rows contain following data: task index i=1,...,n, task
duration t_i, i=1,...,n, number of successors n_{succ} and list of
n_{succ} pairs s_j  c_{ij} representing index of j-th successor
and corresponding communication.

