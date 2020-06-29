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
https://github.com/IOHprofiler/IOHexperimenter for more information.

The framework uses a precompiled static library (`lib/libiohexperimenter.a`). If you get any linking errors when compiling, please try:
```
./get_ioh_lib
```
This will clone the IOHexperimenter repository, recompile the library, and place it in the `lib` folder.

Please cite this work as follows:

```
@misc{boks2020modular,
    title={A Modular Hybridization of Particle Swarm Optimization and Differential Evolution},
    author={Rick Boks and Hao Wang and Thomas Bäck},
    year={2020},
    eprint={2006.11886},
    archivePrefix={arXiv},
    primaryClass={cs.NE}
}
```





