#Project Title

This software is a c++ implementation of the algorithm proposed by Jia Wang at [here](https://arxiv.org/pdf/1512.06021v1.pdf). To compile this project, the Eigen has to be put at the root directory, or you can create a link flag during the compilation. The software has been developed and tested in the following platform: Ubuntu 14.04 with gcc 4.8.4. 

## Getting Started

Build the project simply by typing

```
make
```

Then you could run the demo using a hospital dataset by

```
./bin/netcart
```

The resulting model parameters would be stored in four files: `Const.txt`, `R.txt`, `W.txt`, `X.txt`. You could find their meaning in the paper.

## Input Format

The input files are: `mygraph.edges`, `mygraph.nodefeat`, which describes your graph and node attributes.

In each file, we have rows of data. In `mygraph.edges`, a row `1 2` means there is an edge from 1 to 2. In `mygraph.nodefeat` it means node 1 has attribute 2.

In case that there is some isolated node, whose information is not appearing in the .edge file. One has to list all the nodes by indicating there is a link to the node itself. This is just for the sake of IO. The self referential edges would not be considered in the actual algorithm.

E.g. For the following graph, 1 is only connected to 2. And 3,4 are two isolated graph. Then in the .edge file, one should see
```
1 2
2 1
3 3
4 4
```
The same rule applies to .nodefeat file, if some attributes is not used by any node (which is highly impossible)
