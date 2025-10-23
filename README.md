# n5GDHC_sim

This repository provides the Modelica implementation of a fifth generation district heating and cooling (5GDHC) network, inspired by the network in Gutach-Bleibach, Germany.

The main components of the network are the ice storage, the pipe network (including the supply and return pipes and the surrounding soil) and the heat transfer stations located at consumers' premises. These can be found in the `components` folder.
To simulate a basic demo of a network model with an ice storage system that supplies heat and cold to buildings, run `n5GDHC.mo` in the `systems` folder.

-------------
The simulation of the model was tested with OpenModelica 1.25.0 and Modelica Standard Library 4.0.0.

## Citation

If you use this repository in your work, please cite the following publication:

> Kollmar, M., Bürger, A., Bohlayer, M., Altmann-Dieses, A., Braun, M., Diehl, M. *A detailed simulation model for fifth generation district heating and cooling networks with seasonal latent storage evaluated on field data*. Applied Energy, Volume 401, Part B, 2025, 126757, ISSN 0306-2619,[DOI: [10.48550/arXiv.2506.18528](https://doi.org/10.48550/arXiv.2506.18528)](https://doi.org/10.1016/j.apenergy.2025.126757.)
