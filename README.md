# rt

A multi-threaded Monte-Carlo geometric ray tracer to calculate the spectrally resolved flux density distribution on targets. This is __not__ about nicely rendering scenes like does POV-Ray. Details see the supplied (sparse) documentation and the examples given.


## Features


## Installation

### Requirements

```
gcc
libconfig
gsl
BLAS library
Geomview (or similar to visualize scene from OFF files)
```
Download zip and extract or clone repository. From the resulting folder run

```bash
$ xmkmf -a
$ make
```

## Usage

```bash
$ rt -h

rt Version: v2.1-dirty  (2017-01-25) released by <Ivo Alxneit> ivo.alxneit@psi.ch

Usage: rt
       --append, -a      append to output files. new seed must be given.
       --keep_closed, -k keep output files closed. [keep open].
       --Log, -L         log path of first n rays per source [0].
                         Raw format used.
       --log, -l         log path of first l rays per source [0].
                         OFF format used.
       --mode, -m        select run mode [0].
                         0: check and print input.
                         1: output geometry (OFF files).
                         2: run simulation.
       --threads, -t     number of threads to use [1]
       --Version, -V     Print version number
       --help, -h        Print this help message

```

or with one of the provided tests

```
$ rt < test_screen.cfg
seed = 1234;
P_factor = 0.0001;
sources = ( 
  {
    name = "s_1";
    type = "uniform point source";
    origin = [ 0.0, 0.0, 1.0 ];
    power = 10.0;
    spectrum = "lamp_spec_1.dat";
  }, 
  {
    name = "s_2";
    type = "uniform point source";
    origin = [ 0.0, 1.0, 1.0 ];
    power = 50.0;
    spectrum = "lamp_spec_1.dat";
  } );
targets = ( 
  {
    name = "screen1";
    type = "two-sided plane screen";
    point = [ 0.0, 0.5, 1.0 ];
    normal = [ 0.0, 1.0, 0.0 ];
    x = [ 1.0, 0.0, 0.0 ];
  }, 
  {
    name = "screen2";
    type = "one-sided plane screen";
    point = [ 0.0, 2.0, 1.0 ];
    normal = [ 0.0, -1.0, 0.0 ];
    x = [ 1.0, 0.0, 0.0 ];
  }, 
  {
    name = "screen3";
    type = "one-sided plane screen";
    point = [ 0.0, -1.0, 1.0 ];
    normal = [ 0.0, 1.0, 0.0 ];
    x = [ 1.0, 0.0, 0.0 ];
  }, 
  {
    name = "screen4";
    type = "one-sided plane screen";
    point = [ 0.0, 0.5, 0.0 ];
    normal = [ 0.0, 0.0, 1.0 ];
    x = [ 1.0, 0.0, 0.0 ];
  } );
```
$ rt < input_file
```

