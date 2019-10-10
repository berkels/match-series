# Joint non-rigid registration and reconstruction of image series

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

We appreciate any feedback on your experience with our methods. In case you encounter any problems when using this software, please do not hesitate to contact us: <berkels@aices.rwth-aachen.de>

## Prerequisites
To build the software, cmake and a compiler like GCC or clang need to be installed. Furthermore, a few external libraries are needed. Under Linux, usually all of the dependencies can be installed with the package manager that comes with the Linux distribution. Under OS X/macOS, the installation of a third party package manager like [MacPorts](https://www.macports.org/) or [Homebrew](https://brew.sh/) is recommended. For instance, if MacPorts is installed, one can install cmake with:

    sudo port install cmake

## Compiling

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

## Running
After the code is successfully compiled, there is an executable called `matchSeries` in the `projects/electronMicroscopy` subdirectory of the `quocGCC` directory.

The program is controlled over the command line and a parameter file, there is no GUI.

To start the program, change into the subdirectory with the executable and start the program with

    ./matchSeries parameterfilename

An example parameter file can be found at `quocmesh/projects/electronMicroscopy/matchSeries.par`. The first important parameter one needs to adjust in the parameter file is `templateNamePattern`. It specifies where the input images are located and assumes that all files have the same name except for a frame counter.

For instance, a valid name pattern is

    "/Volumes/Data/Cpp/input/StEM/ADF%04d.dm3"

The frame counter is encoded in "printf" style, e.g. `%04d` means here is the counter, it has four digits and unused digits are filled with zeros. The actual values are configured with the parameters `templateNumOffset`, `templateNumStep` and `numTemplates`. For instance,

    templateNumOffset 0
    templateNumStep 1
    numTemplates 20

leads to the values 0000, 0001, 0002, ..., 0019.

The second important set of parameters to set are the so-called levels, which specify the resolution of the input images and the parameters of the multilevel minimization. The multilevel algorithm expects the input images to be quadratic and to have a dyadic resolution, i.e. the number of pixels in one coordinate direction is a power of two. The exponent is called level. For instance, level 8 means that the input images have $2^8=256$ pixels in each direction. Thus, 9 and 10 encode 512 and 1024 pixels respectively.

The parameter `precisionLevel` specifies the level used to store and output the data. For instance, if the input images are of size 512x512, precisionLevel has to be set to 9. There are numerous other levels that also need to be set properly. Usually, `stopLevel`, `refineStartLevel` and `refineStopLevel` can be chosen to be the same as `precisionLevel`. The last level you need to change is `startLevel`. This controls how coarse the multilevel hierarchy is stated. So, `startLevel` should be smaller than `stopLevel`/`precisionLevel`. Just small enough that you still see the individual atoms if the input images are downsampled to the `startLevel`. Usually `startLevel` is one or two smaller than `stopLevel`.

To be more precise: In terms of Algorithm 1 of the Ultramicroscopy paper, we have $m_0=$ `startLevel` and $m_1=$ `stopLevel` in the first for loop of Algorithm 2, and $m_0=$ `refineStartLevel` and $m_1=$ `refineStopLevel` for the second for loop in Algorithm 2.

Thus, a good starting guess for input images of size 512x512 would be

    startLevel 7
    stopLevel 9
    precisionLevel 9
    refineStartLevel 9
    refineStopLevel 9
