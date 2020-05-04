This is a framework that allows instantiating arbitrary combinations of existing Particle Swarm Optimization and Differential Evolution algorithms as well as hybrid PSO/DE algorithms in a very simple way. Thousands of different combinations are possible.

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

See `experiment.cc` and `mpi\_experiment.cc` for example experiments.

This framework uses IOHexperimenter for benchmarking. Please consult:
```
https://github.com/IOHprofiler/IOHexperimenter
```

for more information.

The framework uses a precompiled static library (`lib/libiohexperimenter.a`). If you get any linking errors when compiling, please try:
```
./get_ioh_lib
```
This will clone the IOHexperimenter repository, recompile the library, and place it in the `lib` folder.



