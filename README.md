# mv-network-test-code
Code to reproduce simulations and data application from Gao, L.L., Witten, D., and Bien, J. Testing for Association in Multi-View Network Data.

## Required R libraries:
multiviewtest, simulator, randnet, doParallel, Matrix, igraph, dplyr, ggplot2, grid, gridExtra, scales, patchwork

## Code Dependencies: 
Code for some figures (e.g. Figures 2, 3 and 4 in the main text) depend on the results of auxiliary files. Before running code for the figures, first: 
* Run main.R in Sec_7_Subsec_1_2 with this.sim.id ranging from 1 to 12. 
* Run main.R in Sec_7_Subsec_3_Web_Appendix_C with this.sim.id ranging from 1 to 4.
* Run init.R in Web_Appendix_B, then run continue.R in Web_Appendix_B with batch ranging from 1 to 10. 
* Run Sec8.R.
* Run main.R in Web_Appendix_D with this.sim.id ranging from 1 to 3. 
* Run main.R in Web_Appendix_E with this.sim.id ranging from 1 to 3.

## Data Dependencies: 
The binary interaction and co-complex association human protein-protein interaction data sets can be downloaded from the public [HINT database](http://hint.yulab.org/download). Here are the direct links: 
1. [Binary interaction data set](http://hint.yulab.org/download/HomoSapiens/binary/hq)
2. [Co-complex association data set](http://hint.yulab.org/download/HomoSapiens/cocomp/hq)