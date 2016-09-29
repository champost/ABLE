##ABLE

`ABLE` is a program written in C/C++ for the joint inference of arbitrary population histories and the genome-wide recombination rate using data from multiple whole genome sequences or fragmented assemblies (e.g. UCE's, RADSeq, and targeted exomes). It makes use of the distribution of blockwise SFS (bSFS) patterns which retain information on the variation in genealogies spanning short-range linkage blocks across the genome. `ABLE` does not require phased data as the bSFS does not distinguish the sampled lineage in which a mutation has occurred. Like with the SFS, outgroup information can be also be ignored by folding the bSFS. `ABLE` takes advantage of `openmp` parallelization and is tailored for studying the population histories of model as well as non-model species.

`ABLE` stands for Approximate Blockwise Likelihood Estimation.

##InstallABLE

###Linux

It is easiest to build an `ABLE` binary under all flavours of Linux. `ABLE` requires the GNU Compiler Collection ([`gcc`](https://gcc.gnu.org/)) 
and GNU Make ([`make`](https://www.gnu.org/software/make/))
for a smooth installation and has been tested using `gcc 4.8.4` and `make 3.81`. If you don't have `gcc` or `make`, you can use your OS specific package handling utility. 

Under Ubuntu this corresponds to the following in a terminal

    sudo apt-get install build-essential


Other dependencies such as the GNU Scientific Library ([GSL](http://www.gnu.org/software/gsl/)) and the Non-Linear Optimization ([NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt)) library are automatically installed by following the instructions outlined below.

* Download the `ABLE` repository

        wget https://github.com/champost/ABLE/archive/master.tar.gz

* Untar the archive and change directory

        tar -xzf master.tar.gz && cd ABLE-master

* If you are installing `ABLE` for the first time you might have to install the `GSL` and `NLopt` libraries. This can take some time as the command below performs a **static installation** of the libraries. You can skip this step if you already have these libraries installed system-wide.

        make deps

* Finally, build an `ABLE` binary

        make clean && make all

If you want `ABLE` to be accessible from everywhere, such as your data folder, you might want to

    cp ABLE ~/bin

This ensures that you can execute the program by specifying `ABLE ...` instead of `./ABLE ...` from the installation folder. This holds only if `~/bin` exists and is part of your `$PATH` environment variable.

##EnjoyABLE

* `ABLE` documentation will be available shortly at : [https://github.com/champost/ABLE/tree/master/doc](https://github.com/champost/ABLE/tree/master/doc)
* Please report bugs and feature requests at : [https://github.com/champost/ABLE/issues](https://github.com/champost/ABLE/issues)
* General questions regarding `ABLE` and for the benefit of all can be posted to [groupABLE](https://groups.google.com/forum/#!forum/groupable).
* A [bioRxiv](http://biorxiv.org/) draft illustrating the performance of `ABLE` will be referenced here shortly.
