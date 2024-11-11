# Core Maintenance
This repository contains C++ codes and datasets for the paper.

- [Core Maintenance](#core-maintenance)
  - [Introduction](#introduction)
  - [Environment](#environment)
  - [Datasets](#datasets)
  - [How to Run the Codes](#how-to-run-the-codes)
    - [A. Code Compilation](#a-code-compilation)
    - [B. Experimentation](#b-experimentation)
      - [1. Algorithms](#1-algorithms)
      - [2. Preparation](#2-preparation)
      - [3. Edge Insertions](#3-edge-insertions)
      - [4. Edge Deletion](#4-edge-deletion)
      - [5. Batch Update](#5-batch-update)
      - [6. Query](#6-query)
  - [Running Toy Graph](#running-toy-graph)
    - [1. Toy Graph](#1-toy-graph)
    - [2. Edge Insertion](#2-edge-insertion)
    - [3. Edge Deletion](#3-edge-deletion)


## Introduction
We present efficient algorithms for maintaining ($\alpha$, $\beta$)-core (bi-core) over bipartite graphs. Our approach introduces the concept of bi-core numbers and uses them to narrow down the scope of updates caused by edge insertions and deletions. Based on this concept, we propose efficient algorithms for maintaining ($\alpha$, $\beta$)-core, for both edge insertion and edge deletion. 



## Environment

The algorithms are implemented in C++ and compiled with the g++ compiler at -O3 optimization level. The development environment used for implementing and testing is:

- Linux version: Ubuntu 20.04.3 LTS
- Kernel version: 5.4.0-139-generic
- g++ version: 9.4.0



## Datasets

We use eight large bipartite graphs, including seven real graphs (IMDB, WC, AR, DBLP, DE, DTI, WT) obtained from [KONECT](http://konect.cc/networks/) and one synthetic graph (PL) generated according to the bipartite network model[6], following the power-law distribution. You can generate PL by running the following command in the terminal:

    1) cd random
    2) g++ bigraph.cpp gephi.cpp main.cpp random.cpp utility.cpp -O3 -o rand.o
    3) ./rand.o output_file_path 5000000 5000000 1000000000

Our paper contains detailed statistics and original sources of the above datasets.

**Dataset Format:** 

- `graph.meta` file containing the number of vertices on each side and the number of edges.
- `graph.e` file containing all the edges.

Example:

```
graph.meta:
3 //number of vertices in U-side
3 //number of vertices in V-side
6 //number of edges

graph.e:
//U V
1 1
1 3
2 1
2 2
2 3
3 1
```
**Generation of Insertion and Deletion Edges:** 
Generate the edge list file to be inserted/deleted randomly:

   ```bash
   python ../random_edges_generator.py graph_path edge_num edgelist_path
   ```    


## How to Run the Codes

### A. Code Compilation
#### 1. Prerequisites

1.1 C++ Compiler

1.2 Boost Library

1.2.1 Installing Boost on Ubuntu:

    sudo apt-get install libboost-all-dev
  
1.2.2 Installing Boost on MacOS using Homebrew:
     
    brew install boost

1.3 Configure the Makefile

Both the `insertion-maintenance` and `deletion-maintenance` folders contain a `Makefile`. Make sure to correctly configure the Boost include and library paths in the Makefile to ensure successful compilation.

`BOOST_ROOT=PATH_TO_BOOST_ROOT`: This variable defines the root path of the Boost installation. Ensure the BOOST_ROOT variable points to the correct Boost installation directory. For example, on macOS with Homebrew: `BOOST_ROOT=/opt/homebrew/opt/boost`; on Linux, it might look like: `BOOST_ROOT=/usr/local/boost`.

1.4 Use corresponding code to compile:

    1) cd insertion-maintenance/ or cd deletion-maintenance/
       make clean
       make
    2) cd bicore-index-maintenance
       make clean
       make

After compilation, executable files named `bicore` and `abcore` will be generated.



### B. Experimentation

The commands for different algorithms are stored in the folder "command". So, run ```cd command``` first.

#### 1. Algorithms

- Recompute: the state-of-the-art bi-core decomposition algorithm.
- BiCore-Index-Ins* (BII): the bi-core maintenance algorithm to handle an edge insertion [26].
- BiCore-Index-Rem* (BIR): the bi-core maintenance algorithm to handle an edge deletion [26].
- Edge-Insert (EI): our proposed bi-core maintenance algorithm for handling an edge insertion.
- Edge-Delete (ED): our proposed bi-core maintenance algorithm for handling an edge deletion.
- BI-batch: the batch algorithm in [26].
- Edge-batch: our batch algorithm.

|         Parameter         |                         Description                          |
| :-----------------------: | :----------------------------------------------------------: |
|       `graph_path`        |                    The path of the graph.                    |
|    `vertex_1/vertex_2`    |    (vertex_1, vertex_2) is the inserted or deleted edge.     |
| `insert_path/delete_path` |            The inserted/deleted edges list file.             |
|        `edge_num`         |       The total number of inserted and deleted edges.        |
|       `U_num/V_num`       |     The number of left/right side vertices in the graph.     |
|    `α/β`     |        ($\alpha$, $\beta$) is the given integer pair.        |
|         `vertex`          |                      The given vertex.                       |
|         `is_left`         | Indicates whether the given vertex belongs to U or V. 1 (True) for U, 0 (False) for V. |



#### 2. Preparation

| Parameter of Command Line |         Description         |          Command           |
| :--------------------: | :--------------------------: | :--------------------------: |
|           `BBI`           |     Build BiCore-Index      | ```sh BBI.sh graph_path``` |
|           `BBN`           | Compute all bi-core numbers | ```sh BBN.sh graph_path``` |



#### 3. Edge Insertions

3.1 Insert a specific edge (vertex_1, vertex_2)

| Parameter of Command Line |       Description       |                   Command                    |
| :-----------------------: | :---------------------: | :------------------------------------------: |
|           `RCI`           |        Recompute        | ```sh RCI.sh graph_path vertex_1 vertex_2``` |
|           `BII`           | BiCore-Index-Ins* (BII) | ```sh BII.sh graph_path vertex_1 vertex_2``` |
|           `EI`            |    Edge-Insert (EI)     | ```sh EI.sh graph_path vertex_1 vertex_2```  |



3.2 Insert edges from an edge list file

| Parameter of Command Line |      Description       |                 Command                 |
| :-----------------------: | :--------------------: | :-------------------------------------: |
|          `RCIS`           |       Recompute        | ```sh RCIS.sh graph_path insert_path``` |
|          `BIIS`           | BiCore-Index-Ins*(BII) | ```sh BIIS.sh graph_path insert_path``` |
|           `EIS`           |    Edge-Insert (EI)    | ```sh EIS.sh graph_path insert_path```  |



#### 4. Edge Deletion

4.1 Delete a specific edge (vertex_1, vertex_2)

| Parameter of Command Line |     Description     |                   Command                    |
| :---: | :------------------------: | :------------------------------------------: |
|           `RCR`           |       Recompute        | ```sh RCR.sh graph_path vertex_1 vertex_2``` |
|           `BIR`           | BiCore-Index-Rem* (BIR) | ```sh BIR.sh graph_path vertex_1 vertex_2``` |
|           `ED`            |    Edge-Delete (ED)    | ```sh ED.sh graph_path vertex_1 vertex_2```  |



4.2 Delete edges from an edge list file

| Parameter of Command Line |     Description       |                 Command                 |
| :-----------------------: | :--------------------: | :-------------------------------------: |
|          `RCRS`           |       Recompute        | ```sh RCRS.sh graph_path delete_path``` |
|          `BIRS`           | BiCore-Index-Rem* (BIR) | ```sh BIRS.sh graph_path delete_path``` |
|           `EDS`           |    Edge-Delete (ED)    | ```sh EDS.sh graph_path delete_path```  |



#### 5. Batch Update

| Parameter of Command Line | Description |                       Command                       |
| :-----------------------: | :---------: | :-------------------------------------------------: |
|         `BBatch`          |  BII-batch  | ``` sh BBatch.sh graph_path edge_num U_num V_num``` |
|         `EBatch`          |  EI-batch   | ``` sh EBatch.sh graph_path edge_num U_num V_num``` |



#### 6. Query

6.1 bi-core query

| Parameter of Command Line |  Description  |                    Command                    |
| :-----------------------: | :-----------: | :-------------------------------------------: |
|          `QCBI`           | BiCore-Index  | ``` sh QCBI.sh graph_path α β``` |
|          `QCBN`           | BiCore number | ```sh QCBN.sh graph_path α β```  |



6.2 community search

| Parameter of Command Line |  Description  |                           Command                            |
| :-----------------------: | :-----------: | :----------------------------------------------------------: |
|          `QSBI`           | BiCore-Index  | ``` sh QSBI.sh graph_path α β vertex is_left``` |
|          `QSBN`           | BiCore number | ```sh QSBN.sh graph_path α β vertex is_left```  |



6.3 $\alpha$ / $\beta$-offsets query

| Parameter of Command Line |  Description  |                   Command                   |
| :-----------------------: | :-----------: | :-----------------------------------------: |
|          `QOBI`           | BiCore-Index  | ``` sh QOBI.sh graph_path vertex is_left``` |
|          `QOBN`           | BiCore number | ```sh QOBN.sh graph_path vertex is_left```  |






## Running example on toy graph

### 1. Toy Graph

- Number of vertices: 6, number of edges: 6

- Edge list:
  - $u_1$ - $v_1$
  - $u_1$ - $v_3$
  - $u_2$ - $v_1$
  - $u_2$ - $v_2$
  - $u_2$ - $v_3$
  - $u_3$ - $v_1$

- We can obtain the original bi-core numbers of all vertices:
  
  | U                            | V                             |
  | :--------------------------- | ----------------------------- |
  | $u_1$: (1, 3), (2, 2)        | $v_1$: (1, 3), (2, 2), (3, 1) |
  | $u_2$: (1, 3), (2, 2), (3,1) | $v_2$: (3, 1)                 |
  | $u_3$: (1, 3)                | $v_3$: (2, 2), (3, 1)         |
  

### 2. Preparation

2.1 To build the BiCore-Index, we can run:

    sh BBI.sh ../data/toy-graph/

Then, we can get the BiCore-Index of the graph.


2.2 To compute all bi-core numbers, we can run:

    sh BBN.sh ../data/toy-graph/

Then, we can get the original bi-core numbers of all vertices.
  

### 3. Edge Insertion

3.1 Suppose we want to insert edge ($u_1$, $v_2$), we can run one of the three commands:

    sh RCI.sh ../data/toy-graph/ 1 2

    sh BII.sh ../data/toy-graph/ 1 2

    sh EI.sh ../data/toy-graph/ 1 2

The updated bi-core numbers of all vertices are stored in "bi-core-number-ins.txt". The content is listed as follows:

| U                     | V                    |
| :--------------------- | :-------------------- |
| $u_1$: (1, 3), (3, 2) | $v_1$:(1, 3), (3, 2) |
| $u_2$: (1, 3), (3, 2) | $v_2$: (3,2)         |
| $u_3$: (1, 3)         | $v_3$: (3, 2)        |


3.2 Suppose we want to insert edges from an edge list file, we can run one of the three commands:

    sh RCIS.sh ../data/toy-graph/ ../data/insert-edges.txt

    sh BIIS.sh ../data/toy-graph/ ../data/insert-edges.txt

    sh EIS.sh ../data/toy-graph/ ../data/insert-edges.txt

The updated bi-core numbers are listed as follows:

| U                     | V                    |
| :--------------------- | :--------------- |
| $u_1$: (3, 3) | $v_1$: (3, 3)  |
| $u_2$: (3, 3) | $v_2$: (3, 3)  |
| $u_3$: (3, 3) | $v_3$: (3, 3)  |


### 4. Edge Deletion

4.1 Suppose we want to delete edge ($u_1$, $v_2$) after the edge insertion, we can run one of the three commands:

    sh RCR.sh ../data/toy-graph-AfterIn/ 1 2
    
    sh BIR.sh ../data/toy-graph-AfterIn/ 1 2

    sh ED.sh ../data/toy-graph-AfterIn/ 1 2

The updated bi-core numbers of all vertices are stored in "bi-core-number-rem.txt". The content is listed as follows:

| U                            | V                           |
| :---------------------------- | :--------------------------- |
| $u_1$: (1, 3), (2, 2)        | $v_1$:(1, 3), (2, 2), (3,1) |
| $u_2$: (1, 3), (2, 2), (3,1) | $v_2$: (3,1)                |
| $u_3$: (1, 3)                | $v_3$: (2, 2), (3,1)        |


4.2. Suppose we want to delete edges from an edge list file, we can run one of the three commands:

    sh RCRS.sh ../data/toy-graph-AfterIn/ ../data/delete-edges.txt
  
    sh BIRS.sh ../data/toy-graph-AfterIn/ ../data/delete-edges.txt
  
    sh EDS.sh ../data/toy-graph-AfterIn/ ../data/delete-edges.txt

The updated bi-core numbers are listed as follows:

| U                     | V                    |
| :--------------------- | :--------------- |
| $u_1$: (1, 2) | $v_1$: (1, 1)  |
| $u_2$: (1, 2), (2, 1) | $v_2$: (1, 2), (2, 1)  |
| $u_3$: (1, 1) | $v_3$: (2, 1)  |


### 5. Batch Update

Suppose we want to insert or delete two edges randomly from the graph, we can run one of the two commands:

    sh BBatch.sh ../data/toy-graph/ 2 4 4

    sh EBatch.sh ../data/toy-graph/ 2 4 4
    
This process is random, so we don't list the result here.

### 6. Query

6.1 Suppose we want to query the (1, 1)-core in the graph, we can run one of the two commands:

    sh QCBI.sh ../data/toy-graph/ 1 1

    sh QCBN.sh ../data/toy-graph/ 1 1

The result is: $u_1$, $u_2$, $u_3$, $v_1$, $v_2$, $v_3$.
    

6.2 Suppose we want to query the community which is a (1, 1)-core and contains vertex $u_1$, we can run one of the two commands:

    sh QSBI.sh ../data/toy-graph/ 1 1 1 1

    sh QSBN.sh ../data/toy-graph/ 1 1 1 1

The result is: $u_1$, $u_2$, $u_3$, $v_1$, $v_2$, $v_3$.


6.3 Suppose we want to query the $\alpha$-offsets for vertex $u_1$, we can run one of the two commands:

    sh QOBI.sh ../data/toy-graph/ 1 1

    sh QOBN.sh ../data/toy-graph/ 1 1

The result is: 2.



## Citation
If you use the code, please cite our paper:

    @article{luo2023efficient,
      title={Efficient Core Maintenance in Large Bipartite Graphs},
      author={Luo, Wensheng and Yang, Qiaoyuan and Fang, Yixiang and Zhou, Xu},
      journal={Proceedings of the ACM on Management of Data},
      volume={1},
      number={3},
      pages={1--26},
      year={2023},
      publisher={ACM New York, NY, USA}
    }
