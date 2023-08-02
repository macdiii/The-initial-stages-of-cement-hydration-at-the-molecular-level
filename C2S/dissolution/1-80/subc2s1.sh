
#!/bin/bash --login

# Replace "[your-project]" with the appropriate project code
# Maximum wall-clock time limit 24 hours (--time=24:00:00)
#SBATCH --job-name=c2s_py
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --account=pawsey0795
#SBATCH --time=23:33:55
#SBATCH --cpus-per-task=1


# Load the lammps module so we can find the "lmp_mpi"
# executable

module load python/3.9.7

# Launch with srun (essential) using 128 MPI tasks ("-n 128")

srun --export=all -N 1 -n 128 -c 1 lmp -in 1.py -log py.log
