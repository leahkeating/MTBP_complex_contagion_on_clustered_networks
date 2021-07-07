A Multi-Type Branching Process Method for Modelling Complex Contagion on
Clustered Networks
================

## Description of the files

1.  The
    [mtbp\_functions.R](https://github.com/leahkeating/MTBP_complex_contagion_on_clustered_networks/blob/main/mtbp_functions.R)
    file contains all of the functions required for the other files. The
    functions include:
      - A function to create the mean matrix (mean\_mat( ))
      - A function that calculates the expected cascade size
        analytically (expected\_size()). THe analytical calculation of
        the expected cascade size is discussed in Sec. IV A.
2.  The
    [cascade\_condition.R](https://github.com/leahkeating/MTBP_complex_contagion_on_clustered_networks/blob/main/cascade_condition.R)
    file contains code to generate the plots showing the critical
    boundaries for the 4 networks, similarly to Fig. 3 in the paper.
      - The code details how to generate the boundary lines shown in the
        plots below:
        <img src="20210624_newman_crit.png" width="45%" /><img src="20210624_clique_crit.png" width="45%" />
3.  The
    [MTBP\_simulations.R](https://github.com/leahkeating/MTBP_complex_contagion_on_clustered_networks/blob/main/MTBP_simulations.R)
    file contains the code for doing the branching process type
    simulations.
