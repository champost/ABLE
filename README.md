## ABLE

`ABLE` is a program for the joint inference of arbitrary population histories and the genome-wide recombination rate using data from multiple whole genome sequences or fragmented assemblies (e.g. UCE's, RADSeq, and targeted exomes). It makes use of the distribution of blockwise SFS (bSFS) patterns which retain information on the variation in genealogies spanning short-range linkage blocks across the genome. `ABLE` does not require phased data as the bSFS does not distinguish the sampled lineage in which a mutation has occurred. Like with the SFS, outgroup information can be also be ignored by folding the bSFS. `ABLE` takes advantage of `openmp` parallelization and is tailored for studying the population histories of model as well as non-model species.

`ABLE` stands for Approximate Blockwise Likelihood Estimation and is written in C/C++.

## InstallABLE

### Linux

It is easiest to build an `ABLE` binary under all flavours of Linux. `ABLE` requires the GNU Compiler Collection ([`gcc`](https://gcc.gnu.org/)) 
and GNU Make ([`make`](https://www.gnu.org/software/make/))
for a smooth installation and has been tested using `gcc 4.8.4` and `make 3.81`. It is also useful to have [`git`](https://git-scm.com/) installed for seamless updates of the latest version of `ABLE`. If you do not have `gcc`, `make` or `git`, you can use your OS specific package handling utility. 

Under Ubuntu this corresponds to the following in a terminal

    sudo apt-get install build-essential git


Other dependencies such as the GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)) and the Non-Linear Optimization ([NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt)) library are automatically installed by following the instructions outlined below.

1. Clone the `ABLE` repository. If `git` is not available see the next section (or [here](https://git-scm.com/downloads)).

        git clone https://github.com/champost/ABLE.git

2. Change directory

        cd ABLE

3. If you are installing `ABLE` for the **first time** you might have to compile the `GSL` and `NLopt` libraries from source. This can take some time as the command below performs a **static installation** of the libraries. You can skip this step if you already have these libraries installed system-wide or if you are simply updating `ABLE` to the latest version (see next section).

        make deps

4. Finally, build an `ABLE` binary

        make clean && make all

If you want `ABLE` to be accessible from everywhere, such as your data folder, you might want to

    cp ABLE ~/bin

This ensures that you can execute the program by specifying `ABLE ...` instead of `./ABLE ...` from the installation folder. This holds only if `~/bin` exists and is part of your `$PATH` environment variable.

## UpdatABLE

### With Git

As mentioned above, `git` considerably simplifies the process (or headache) of updating `ABLE`. You simply have to change your current working directory to that of `ABLE` and execute the following

    git pull
    make clean && make all

### Without Git

If `git` not available on your machine you can try looking for solutions [here](https://git-scm.com/downloads). Alternatively, you can download the `ABLE` repository somewhere on your computer as shown below.

    wget https://github.com/champost/ABLE/archive/master.tar.gz

Untar the archive using 

    tar -xzf master.tar.gz

and manually compare the contents of the `ABLE-master` directory with your `ABLE` installation and replace the files that have changed in size.

### Updating dependencies

With or without `git`, you will seldom need to recompile the libraries that `ABLE` depends upon and in which case you will have to run the `make deps` command before compiling an `ABLE` binary with `make clean && make all`. This information will be clearly stated in the [Release notes](https://github.com/champost/ABLE/releases) and/or the `ABLE` [documentation](https://github.com/champost/ABLE/blob/master/doc/helpABLE.pdf).


## EnjoyABLE

* You can read the full `ABLE` documentation [online](https://github.com/champost/ABLE/blob/master/doc/helpABLE.pdf) or downoad the [PDF](https://github.com/champost/ABLE/raw/master/doc/helpABLE.pdf).
* Please report bugs and feature requests at : [https://github.com/champost/ABLE/issues](https://github.com/champost/ABLE/issues).
* General questions/announcements regarding `ABLE` and for the benefit of all can be posted to [groupABLE](https://groups.google.com/forum/#!forum/groupable).
* An accompanying draft illustrating the performance of `ABLE` is available on [bioRxiv](http://dx.doi.org/10.1101/077958).
