To compile and run your experiment:
$ make
$ ./experiment

To compile and run an MPI experiment:
$ make mpi
$ mpirun -np [# of instances in suite] mpi_experiment

See experiment.c and mpi_experiment.c for example experiments.

To implement a new module, extend the base class of that module (e.g. for a new velocity update, extend the ParticleUpdateManager class) and implement the member functions. For example, when implementing a new ParticleUpdateManager, you need to extend the update() function. After implementing your ParticleUpdateManager, add a new entry to the UpdateManagerType enum in particleupdatemanager.h that represents your ParticleUpdateManager, and let ParticleUpdateManagerFactory return an instance of your ParticleUpdateManager when the enum entry is passed to createParticleUpdateManager(). You can now pass the integer value from the enum to the constructor of ParticleSwarm, and it will use your ParticleUpdateManager. You should also add your new module to the getIdString() function in the ParticleSwarm or DifferentialEvolution classes, so that the resulting algorithms get an appropriate name.

New modules will automatically be added to the PSO- and DE suites. If you only want to consider only a subset of the available modules, you can do so by using the member functions of ParticleSwarmSuite and DESuite. Please take a look at the header files particleswarmsuite.h and desuite.h for the functions.

This framework uses COCO (COmparing Continuous Optimizers) to measure the algorithms' performance in the experiments. Please consult:

https://github.com/numbbo/coco

for more information on COCO and how to post-process the data resulting from the experiments.



