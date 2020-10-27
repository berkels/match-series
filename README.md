# Joint non-rigid registration and reconstruction of image series

Joint non-rigid registration aims to calculate non-linear deformation artifacts from a series of images taken of an object and reconstruct the true image. 
Applications include removal of scan noise in scanning transmission electron microscopy datasets.

This code implements the methods proposed in the following papers:

[1] Benjamin Berkels, Peter Binev, Douglas A. Blom, Wolfgang Dahmen, Robert C. Sharpley, and Thomas Vogt. Optimized imaging using non-rigid registration. *Ultramicroscopy*, 138:46--56, March 2014. [[DOI](https://dx.doi.org/10.1016/j.ultramic.2013.11.007) | [arXiv](https://arxiv.org/abs/1403.6774)]

[2] Benjamin Berkels and Christian H. Liebscher. Joint non-rigid image registration and reconstruction for quantitative atomic resolution scanning transmission electron microscopy. *Ultramicroscopy*, 198:49--57, March 2019. [[DOI](https://dx.doi.org/10.1016/j.ultramic.2018.12.016) | [arXiv](https://arxiv.org/abs/1901.01709)]

[3] Sönke Reiche and Benjamin Berkels. Automated stacking of seismic reflection data based on non-rigid image matching. Geophysics, 83(3):V171--V183, May 2018. [[DOI](http://dx.doi.org/10.1190/geo2017-0189.1)]

[4] Sönke Reiche, Benjamin Berkels, and Benedikt Weiß. Automated static and moveout corrections of high-resolution seismic data from the Baltic Sea. Near Surface Geophysics, 2019. in press. [[DOI](http://dx.doi.org/10.1002/nsg.12068)]

`matchSeries` implements [1-2], `matchSeismicSeries` implements [3-4]. Note that the latter needs the options `USE_BLAS`, `USE_LAPACK` and `DUSE_SUITESPARSE` to be activated.

It is based on the [QuocMesh software library](https://archive.ins.uni-bonn.de/numod.ins.uni-bonn.de/software/quocmesh/index.html), AG Rumpf, Institute for Numerical Simulation, University of Bonn, and distributed under the terms of the [Common Development and Distribution License](LICENSE.txt). Particularly, you can **download and use** this software **free of cost**.

Some subdirectories contain code from other software projects, which is redistributed under the respective licenses:

* `quocmesh/external/bz2/`

  Code from the [bzip2 library by Julian R Seward](http://www.bzip.org/), License under `external/bz2/LICENSE`

* `quocmesh/external/dirent/dirent.h`

  Code from the [Dirent API by Toni Ronkko](https://github.com/tronkko/dirent), License at the beginning of external/dirent/dirent.h

* `quocmesh/cmake/FindSUITESPARSE.cmake`

  CMake module from [OpenFlipper](https://www.openflipper.org/), GNU Lesser General Public License

We appreciate any feedback on your experience with our methods. We could also appreciate if you cite the above mentioned papers when you use the software in your work. In case you encounter any problems when using this software, please do not hesitate to contact us: <berkels@aices.rwth-aachen.de>

## Installation
Binaries for matchSeries are available for Linux, Windows and Mac OSX through [conda](https://docs.conda.io/en/latest/) on the [conda-forge](https://anaconda.org/conda-forge/match-series) channel. Install with:

```
$ conda install -c conda-forge match-series
```
An experimental python wrapper meant for use in Jupyter notebooks is available on [pypi](https://pypi.org/project/pyMatchSeries/):
```
$ pip install --user pyMatchSeries
```
You may find [an example notebook](https://github.com/din14970/pyMatchSeries/blob/master/examples/example.ipynb) on the respective [Github repo](https://github.com/din14970/pyMatchSeries).

## Using the software

If you installed match-series with conda, `matchSeries` should be on the path. If you open up a terminal (Linux/OSX) or Anaconda prompt (Windows) you should be able to run:

    $ matchSeries <path-to-parameter-file>

If you do not provide a path to a parameter fill, `matchSeries` will give an error.

If you compiled from source (see [below](#compiling)), the executable `matchSeries` is in the `projects/electronMicroscopy` subdirectory of the `quocGCC` directory.

The program is controlled over the command line and a parameter file, there is no GUI.

An example parameter file can be found at `examples/example_parameters.param`. An explanation of each parameter is provided in a comment string above each parameter. Be sure to set the `templateNamePattern` and `levels` (`precisionLevel`, `stopLevel`, `refineStopLevel`, `startLevel`, `refineStartLevel`) for each calculation. If you use the python wrapper, most of these will be filled out automatically.


## <a name="compiling">Compiling</a>

If you wish to compile the software yourself please follow the following instructions.

### Prerequisites
To build the software, cmake and a compiler like GCC or clang need to be installed. Furthermore, a few external libraries are needed. Under Linux, usually all of the dependencies can be installed with the package manager that comes with the Linux distribution. Under OS X/macOS, the installation of a third party package manager like [MacPorts](https://www.macports.org/) or [Homebrew](https://brew.sh/) is recommended. For instance, if MacPorts is installed, one can install cmake with:

    sudo port install cmake

### Procedure

First create a clone of our git repository e.g. by

    git clone https://github.com/berkels/match-series

In the clone directory, there should be two directories, `quocmesh` and `quocGCC`. The source code is in `quocmesh`, the compiled binaries will go into `quocGCC`.

To compile the code, go into the directory `quocGCC` and call the script called `goLinux.sh` (or `goWindows.bat` if compiling under Windows with MinGW). This invokes cmake with all necessary options and will create makefiles.

After the go script has been run, the source can be compiled with `make` like any other project with makefiles. There is also a test target to check whether the basics of the code actually work.

In short, to get started the following steps should be enough: Decompress the archive and change into the directory `quocGCC`. In there, call

    ./goLinux.sh
    make
    make test

Under Windows with MinGW call

    ./goWindows.bat
    mingw32-make
    mingw32-make test

instead.

If all dependencies are installed correctly, this should work without any errors (although there may be warnings during the first two steps).

To speed up the compilation, one can have make use as many threads as desired, by calling `make -jN` (where `N` is the number of threads) instead of just `make`.

