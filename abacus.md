#Submitting jobs to ABACUS

[ABACUS](http://www.abacus.cinvestav.mx/) uses the Simple Linux Utility for Resource Management, aka **slurm**. Information on slurm options can be obtained at http://slurm.schedmd.com/sbatch.html


####CPU nodes:
The ABACUS supercomputer has 268 CPU nodes and each node has 2 [Intel Xeon E5-2697V3 @ 2.60 GHz CPUs](http://ark.intel.com/products/81059/Intel-Xeon-Processor-E5-2697-v3-35M-Cache-2_60-GHz). Each processor has 14 cores and 28 threads. This means that for each node you can request up to 56 (2 x 28) tasks (threads):
In slurm this is:
`-N [1..268]` or `--nodes=[1..268]`
and
`-n [1..56]` or `--ntasks=[1..56]`


####CPU time:
Abacus has a strict policy regarding CPU usage, and therefore it is imperative that requested CPU time is accurate. In slurm time is requested:
For jobs lasting up to 24 hours:
`-t HH:MM:SS` or `--time=HH:MM:SS`
For jobs lasting more than 24 hours:
`-t DD-HH:MM:SS` or `--time=DD-HH:MM:SS`


####Partitions (queues):
The number of hours or days requested will determine the partition (queue) to be used. Abacus has 3 CPU partitions: **1-28cpu-1d** (jobs lasting up to 24 hours), **1-28cpu-7d** (jobs lasting up to 7 days) and **1-28cpu-30d** (jobs lasting up to 1 month), which are requested using:
`-p partition_name` or `--partition=partition_name`

You should also indicate the name of your job, preferably the name of the program used:
`-J jobname` or `--job-name=jobname`

And standard output and standard error can be directed using:
`-e error_filename` or `--error=error_filename`
`-o output_filename` or `--output=output_filename`


####Examples:
A slurm script header for a blastn job that requests 1 node, 20 threads and 1 hour of CPU time should therefore look like this:
```bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p 1-28cpu-1d
#SBATCH -t 01:00:00
#SBATCH -J blastn
#SBATCH -e blastn.err
#SBATCH -o blastn.out
```
or in long form option:
```bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --partition=1cpu1-28
#SBATCH --time=01:00:00
#SBATCH --job-name=blastn
#SBATCH --error=blastn.err
#SBATCH --output=blastn.out
```
A slurm script header for a blastn job that requests 10 nodes, 56 threads per node and 6 days should look like this:
```bash
#SBATCH -N 10
#SBATCH -n 56
#SBATCH -p 1-28cpu-7d
#SBATCH -t 06-00:00:00
#SBATCH -J blastn
#SBATCH -e blastn.err
#SBATCH -o blastn.out
```
