
# NuDEX (a Nuclear DE-eXcitation code)

Authors: Emilio Mendoza Cembranos and Daniel Cano Ott

Contact: [emilio.mendoza@ciemat.es](mailto:emilio.mendoza@ciemat.es) and [daniel.cano@ciemat.es](mailto:daniel.cano@ciemat.es)

Date: 10 July 2024

GitHub: https://github.com/UIN-CIEMAT/NuDEX

## TABLE OF CONTENTS

- [About](#About)
- [Directories](#Directories)
- [Compilation](#Compilation)
- [Data library](#Data-library)
- [Usage](#Usage)
- [How to reference](#How-to-reference)
- [License](#License)
- [References](#References)

## About

NuDEX is a computer toolkit which allows to generate de-excitation cascades of atomic nuclei. It is written in C++ programming language and can be used to create different applications, adapted to specific needs. At the moment, NuDEX is distributed together with two specific applications: one to generate cascades emitted after neutron capture (`NuDEX_NCaptureCascadeGenerator01`), and the other to generate cascades starting from a certain set of excited levels (`NuDEX_DecayCascadeGenerator01`).

NuDEX operates in a similar way as DICEBOX [1] and DEGEN [2]. To perform the γ-ray decay of a nucleus from a certain excited level, NuDEX generates the full level scheme and associated branching ratios, together with the internal conversion coefficients.  To do this, it takes the known values from a database based mainly in RIPL-3 [3] and ENSDF [4], and fills in the missing information with statistical models. The database (or data library) is distributed together with the NuDEX code and allows to generate automatically the de-excitation cascades of ∼3000 nuclei.

More information concerning the performance of NuDEX can be found in the user's manual.

## Directories

- `include`: NuDEX source code - header files
- `src`: NuDEX source code
- `docs`: contains the user's manual
- `applications`: includes the `NuDEX_NCaptureCascadeGenerator01` and `NuDEX_DecayCascadeGenerator01` applications

## Compilation

NuDEX applications have to be compiled along with ROOT (https://root.cern/).

After downloading the source code, the NuDEX applications have to be compiled with a command similar to:

```sh
cd applications/
g++ -std=c++11 ../src/*.cc NuDEX_NCaptureCascadeGenerator01.cc `root-config --libs --cflags` -I../include/ -o NuDEX_NCaptureCascadeGenerator01
g++ -std=c++11 ../src/*.cc NuDEX_DecayCascadeGenerator01.cc `root-config --libs --cflags` -I../include/ -o NuDEX_DecayCascadeGenerator01
```

## Data library

NuDEX data library is available for download from https://github.com/UIN-CIEMAT/NuDEXlib

## Usage

`NuDEX_NCaptureCascadeGenerator01` is used to generate nuclear de-excitation cascades emitted after neutron capture cascades. It can be executed first with:

```sh
./NuDEX_NCaptureCascadeGenerator01 [outputfilebase]  [LIBDIR]  [ZA]
```

where `[outputfilebase]` is the path+base name for the output files; `[LIBDIR]` is the path to the NuDEX database; and `[ZA]` defines the target nucleus (not the compound nucleus) as ZA=1000·Z+A, with Z the number of protons and A the number of nucleons.

For example, to model Pu-239(n,g) cascades:

```sh
NuDEX_NCaptureCascadeGenerator01 output01  [...]/NuDEXlib-1.0 94239
```

Then we will get three different output files: output01.out, with general information; output01.cas, with the cascades; and output01.inp, an input file with all the parameters used by NuDEX to generate these cascades (the default ones, in this case). This input file can be then modified by hand, if needed, to generate new cascades via:

```sh
./NuDEX_NCaptureCascadeGenerator01  [outputfilebase]  [inputfile]
```

In this case:

```sh
./NuDEX_NCaptureCascadeGenerator01  output02 output01.inp
```

And different output files will be generated: output02.out, output02.cas and output02.inp.


Alternatively, some of the parameters can also be changed through the command line. For example:

```sh
./NuDEX_NCaptureCascadeGenerator01 output01  [...]/NuDEXlib-1.0 94239 NCASCADES 1000000
```

will produce 1000000 cascades instead of the default number of cascades (100).

Another example:

```sh
./NuDEX_NCaptureCascadeGenerator01  output02 output01.inp NCASCADES 1000000 SEED4 123147
```

Will use the parameters from output01.inp, but will replace the values there for `NCASCADES` and `SEED4` with the values provided in the command line. Here `SEED4` is the seed to compute the branching ratios of the primary transitions.


The application `NuDEX_DecayCascadeGenerator01` works in a very similar way.

## How to reference

The user can reference NuDEX with the following publication:

E. Mendoza et al., Neutron capture measurements with high efficiency detectors and the Pulse Height Weighting Technique, Nucl. Instr. Meth. A 1047, 167894 (2023)
https://doi.org/10.1016/j.nima.2022.167894


## License

See the LICENSE.txt file for license rights and limitations (GNU GPLv3).


## References

[1] F. Bečvář, Simulation of γ cascades in complex nuclei with emphasis on assessment of uncertainties of cascade-related quantities, Nucl. Instrum. Methods A 417 (2) (1998) 434 – 449. https://doi.org/10.1016/S0168-9002(98)00787-6

[2] D. Jordan et al., An event generator for simulations of complex β-decay experiments, Nucl. Instrum. Methods A 828 (2016) 52 – 57. https://doi.org/10.1016/j.nima.2016.05.034

[3] R. Capote et al., RIPL – reference input parameter library for calculation of nuclear reactions and nuclear data evaluations, Nuclear Data Sheets 110 (12) (2009) 3107 – 3214. https://doi.org/10.1016/j.nds.2009.10.004

[4] https://www.nndc.bnl.gov/ensdf/
