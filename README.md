This is a framework that allows instantiating arbitrary combinations of existing Particle Swarm Optimization and Differential Evolution algorithms as well as hybrid PSO/DE algorithms in a very simple way. Algorithms are generated as combinations of serveral operator modules, such as mutation, crossover, velocity updates, ... etc. Many modules have been implemented, and it is simple to implement new variants. Thousands of different combinations are possible.

To compile and run your experiment:
```
$ make
$ ./experiment
```

To compile and run an MPI experiment (for parallelizing many algorithm instances):
```
$ make mpi
$ mpirun -np [# of instances in suite] mpi_experiment
```

See `experiment.cc` and `mpi_experiment.cc` for example experiments.

This framework uses IOHexperimenter for benchmarking. Please consult:
https://github.com/IOHprofiler/IOHexperimenter for installation instructions.

Please cite this work as follows:

```
@inproceedings{10.1145/3377929.3398123,
author = {Boks, Rick and Wang, Hao and B\"{a}ck, Thomas},
title = {A Modular Hybridization of Particle Swarm Optimization and Differential Evolution},
year = {2020},
isbn = {9781450371278},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3377929.3398123},
doi = {10.1145/3377929.3398123},
booktitle = {Proceedings of the 2020 Genetic and Evolutionary Computation Conference Companion},
pages = {1418–1425},
numpages = {8},
keywords = {metaheuristics, hybridization, differential evolution, swarm intelligence, particle swarm optimization},
location = {Canc\'{u}n, Mexico},
series = {GECCO ’20}
}
```





