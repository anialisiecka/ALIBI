####Performance on *Escherichia coli* genomes
We compared both versions of ALIBI against two linearization algorithms implemented in [vg tool](https://github.com/vgteam/vg):
**Eades** is the implementation of a well-known heuristic for the feedback arc set problem of
[Eades et al. (1993)](https://scholar.google.com/scholar_lookup?journal=Inf.+Process.+Lett.&title=A+fast+and+effective+heuristic+for+the+feedback+arc+set+problem&author=P.+Eades&author=X.+Lin&author=W.F.+Smyth&volume=47&publication_year=1993&pages=319-323&)
and **FP** is the flow procedure proposed in [Haussler et al., 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6067104/).

The following figures present the weighted reversing join (WRJ), weighted feedback arc (WFA) and average cut width (ACW) results of the algorithms obtained
on *Escherichia coli* data set (the full list of genomes can be found [here](https://github.com/anialisiecka/ALIBI/blob/master/example_data/E.coli/acc_numbers.txt)).

ALIBI outperforms the competitors in both WRJ and WFA. However, the differences in WFA are surprisingly extreme here – in some cases,
Eades or FP algorithms have WFA over 500 times larger than ALIBI (note the logarithmic scale on the Y axis). ALIBI_v2 has the lowest ACW.

<img width="733" alt="acw" src="https://user-images.githubusercontent.com/44555028/173303875-806d7a9c-d070-4952-adb0-aaa99aabdbc7.png">
<img width="725" alt="wfa" src="https://user-images.githubusercontent.com/44555028/173303980-bb5ea9df-db2c-4dff-8da4-d679bb3ec6d1.png">
<img width="739" alt="wrj" src="https://user-images.githubusercontent.com/44555028/173304000-5a691d92-4e76-4831-bf5e-fdf178133ed1.png">

The following figures summarize computational efficiency of the algorithms. ALIBI performs best with respect to both time and memory resources.
All algorithms scale roughly linearly with respect to both number of edges and number of genomes. ALIBI_v2 performs slightly worse than ALIBI_v1
with respect to memory resources. However, the difference is negligible and nearly constant.

<img width="733" alt="time" src="https://user-images.githubusercontent.com/44555028/173303956-8f8e7d16-6071-4032-bec8-2a7f32e208fd.png">
<img width="733" alt="memory" src="https://user-images.githubusercontent.com/44555028/173303921-3581407c-122f-45e1-b93f-fd76e43d7936.png">
<img width="726" alt="memoryALIBI" src="https://user-images.githubusercontent.com/44555028/173303942-810c5bed-0870-44c6-aa70-ab2f10b175fc.png">
