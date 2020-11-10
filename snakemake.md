# Snakemake

Here I will demonstrate how we can
use [Snakemake](https://snakemake.readthedocs.io/en/stable/) to
process files and to run a simulation on the cluster.

On our cluster, Snakemake is pre-installed, all you need to do is
`module load python`, and you will have access to the `snakemake`
command. 

## Why use Snakemake?

* Helps others see exactly what steps are involved, helps with
  reproducibility. 
* Figures out job dependencies automatically - if a single job
  crashes, Snakemake will know what needs to be re-run.
* Works on a variety of clusters, can be customized as far as how the
  jobs are sent out.

## Snakefile

The key file that you need to run Snakemake is a `Snakefile`. This can
file will contain `rules` that tell Snakemake how to create the
necessary files. Here is a simple one:

```
rule all: 
  input: "final.txt"

rule make_txt:
  input: "{sample}.foo"
  output: "{sample}.txt"
  shell: "mv {input} {output}"

rule make_final:
  input: "a.txt", 
         "b.txt"
  output: "final.txt"
  shell:
      "cat a.txt b.txt > final.txt"
```

We can run Snakemake in a testing mode with `snakemake -np` which will
print out the lines it would run, but it won't actually run them. If
we do this, it will tell us that it can't find `a.foo` and `b.foo`. So
then create touch files with `touch a.foo` and `touch b.foo`. Now try
again with `snakemake -np`. 

## Running on the cluster

We can run these in an interactive session, or we can submit them to
the cluster scheduler. Nearly always for a real project you would want
to submit jobs to the cluster scheduler (try to aim for each job to
run for at least ~5 min to a few hours on the long end).

The command I use to run on the Longleaf cluster is, from an
interactive node:

```
snakemake -j 4 --latency-wait 30 --cluster "sbatch -N 1 -n 6"
```

This will send 4 jobs at a time, and each job will run on one node,
and require 6 cores (e.g. if I am using Salmon or multi-core R
commands, otherwise set `-n 1`).

This will run Snakemake on the command line, so it will last as long
as your interactive session lasts. You can also submit a job with the
`snakemake` call by putting the above into a bash script and
submitting that to the cluster with `sbatch`, see an example in the
[quantify transcripts](quantify.md) 
instructions.

## Running Snakemake across parameter values

The above works well for when we have input files and want to process
them uniformly, and there are various tricks using `glob_wildcards()`
to process e.g. all files that are in a directory with a certain file
ending (see quantification example).

What about if we don't have input files per se but want to run a
simulation across a grid of parameter values?
