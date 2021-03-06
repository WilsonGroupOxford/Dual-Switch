# Overview 

Monte Carlo code to generate atomic networks which match given ring statistics and Aboav-Weaire parameter.

The algorithm used is "dual-switching" i.e. manipulating the dual graph.

## Compilation 

Compilation is easiest using CMake.
To compile the code execute the following in terminal:
```commandline
cd src/
cmake .
make
```
As the code is relatively computationally intensive, it is recommended that in the
```CMakeCache.txt``` you set ```CMAKE_BUILD_TYPE:STRING=Release``` before making.
You may also need to set ```CMAKE_CXX_FLAGS:STRING=-std=c++11```.

The generated executable is called ```dual_switch```.

## Input

The parameters for the calculation can be found in dual_switch.inpt which are read at runtime. These are explained below:

```text
I/O
2: ./NA       input prefix
3: ./path/to/output/test     output prefix
```

This is the path for the input (if using existing network) and output files and the prefix. 
In this case output file ```test*.dat``` would be placed in the directory ```/path/to/output/```.

```text
Network Properties
5: 1          periodic
6: 0           use existing network (P)
7: 0           override default periodic lattice dimensions (P+L)
8: 200.0 200.0  custom periodic lattice dimensions (P+L)
9: hexagonal   starting lattice type (P for non-hexagonal)
10: 10 10     starting bulk hexagonal lattice dimensions
11: 4 8     ring size limits
12: 0.20      target alpha
13: 0.05 0.15 0.6 0.15 0.05      ring statistics
```

These parameters set the fundamental properties of the target network. 
You probably want a hexagonal periodic lattice (not sure if the other options will work for you tbh), of the given n x m size.
The ring size limits set the allowable range of ring sizes, and the target alpha/ring statistics are the key part, 
these are the network properties the algorithm are aiming for.

```text
Simulation Parameters
17: 0 4         random seed range (inclusive)
18: 0.0050      temperature
19: 100000    maximum proposed monte carlo moves
20: 0.001        energy convergence criteria
21: 10.0          alpha energy scaling factor
```

These control the simulation, the number of seeds to try (i.e. the number of samples).
The others tune the MC process, determining how much the alpha term is weighted etc.

```text
Dual Potential Model and Geometry Optimisation
24: 3.031535        si-si distance
25: 25          harmonic k
26: 1          locally minimise every step
27: 100         maximum iterations
28: 0.01         convergence criteria
29: 0.001        line search step size
30: 100           globally minimise after simulation end
31: 1000       maximum iterations
32: 0.001         convergence criteria
```

The simulation happens in dual space with a simple harmonic potential acting between nodes.
These control that potential.

```text
Analysis Tools
1           convert to atomic network
1           periodic visualisation
0           partial spatial rdfs (P)
0.01        partial rdf bin width
50.0        partial rdf extent
0           partial topological rdfs (P)
12.0        number of topological shells to include
1           ring areas
1           cluster analysis (P)
1           assortative mixing (P)
```

Various analysis methods with various unhelpfully vague outputs.

```text
Atomic Potential Model and Geometry Optimisation
0           globally minimise atomic network
1.42        Keating a parameter
25.880      Keating alpha parameter
5.176       Keating beta parameter
1000         maximum iterations
0.00001        convergence criteria
0.00001        line search step size
```

You can also choose to optimise the atomic network with the Keating potential on completion.

### Runtime

Running ```dual_switch``` will read the input file and start the code. 
Log and output files will be written.

### Output

The code will produce output files in two categories: network description and analysis.
It will also continuously produce a log file.
These will have the names `_*.out`.

#### Network Description 

The network is described via the `<prefix>_<seed>_graph_*.out` and `<prefix>_<seed>_dual_*.out` files.

#### Analysis

Various analysis files can be generated `*_analysis_*.out`

### Visualisation

You will probably want to visualise the networks generated by the code at some point.
I have included a python script to help you get started with this, `visualise.py`.
They are really for my personal use, so good luck.

You probably want `python visualise.py <prefix>_<seed> -gp`
or `python visualise.py <prefix>_<seed> -gpd`