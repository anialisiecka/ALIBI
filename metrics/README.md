You can use [this script](./metrics.py) to calculate the following metrics in your GFA file:
- weighted feedback arc (WFA) – the sum of weights of all feedback arcs, i.e. backward pointing edges
- weighted reversing join (WRJ) – the sum of weights of all reversing joins, i.e. edges joining two in- or two out-sides.

It is assumed thet the order of the nodes is specified by the order of the "S" lines in the input file.

To run the script execute the following shell command:

```
python metrics.py <gfa_input_file>
```
