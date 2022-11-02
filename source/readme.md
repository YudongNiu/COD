# Discover Your Characteristic Community from Graph Hierarchical Communities

Build HIMOR-Index
-------
To build the HIMOR-Index for our COD problem, execute the following command on linux:

```sh
python3 himor.py dataset
```

The parameter "dataset" denotes the name of your network.

Running code
-------
To run the code for COD problem, execute the following command on linux:

```sh
python3 COD.py dataset labeled K cpr loc indexed
```

**There are 6 parameters:**
* dataset: the name of your network
* labeled: if node labels are considered
* K: the query node is top-K influential
* cpr: if compressed COD framework is used
* loc: if lore algorithm is used
* indexed: if the HIMOR-index is used

For example, following commands perform the corresponding algorithms in our experiments on dataset cora with the requirement that the query node is top-5 influential in the result charateristic community.

```sh
python3 COD.py cora 1 5 1 1 1
# CODL
python3 COD.py cora 1 5 1 1 0
# CODL-
python3 COD.py cora 1 5 1 0 0
# CODR
python3 COD.py cora 1 5 0 0 0
# Independent
python3 COD.py cora 0 5 0 0 0
# CODU
```


Input Files
-----------
**The program COD.py requires 3 input files:**
* dataset.txt: edge list of your network
* dataset.feat: node labels
* dataset.query_nodes: COD queries

All the input files should be placed under the folder labeled_graph/dataset/.