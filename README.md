# Progenitors of low-mass binary black-hole mergers in the isolated binary evolution scenario

This repository contains the MESA template used to obtain the population of BBH mergers in the paper from García, Simaz Bunzel, Chaty et al. 2021.

## Pre-requisites

For computing detailed binary evolutionary models, [MESA](http://mesa.sourceforge.net/) must
be installed (supported version: **r10398**), and [MESASDK](http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk)
(version: **20180822**).

## Repository organization

```
.
├── README.md
└── template_r10398
    ├── column_lists
    ├── make
    └── src
        ├── ce
        │   ├── make
        │   ├── private
        │   └── public
        └── core_collapse
            ├── make
            ├── private
            └── public
```

Inside the `src` folder two modules, `ce` and `core_collapse`, can be found which were developed in this work.

The `ce` module handles the evolution during a common-envelope phase (using the energy formalism, [Ivanova et al. 2003](https://ui.adsabs.harvard.edu/abs/2013A%26ARv..21...59I/abstract)) while the
`core_collapse` simply computes relevant masses at the core-collapse stage of a massive star ([Fryer et al. 2012](https://ui.adsabs.harvard.edu/abs/2012ApJ...749...91F/abstract)).

`ce` is controlled via the file `inlist_ce` while `core_collapse` with controls found in `inlist_cc`, both using
namelist as every MESA user control.

## How to run

To compile and run this template, set up the environment variables needed by MESA: `MESA_DIR`, `MESASDK_ROOT`, source
the correspoding `mesasdk_init.*` and run `./mk_mods && ./mk`. This will create the proper libraries for the `ce` and
`core_collapse` modules and link everything, creating a binary file called `bin2dco`. Finally, run a binary model with `./rn`
