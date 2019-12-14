# gaussdecomp
Gaussian Decomposition program mainly for HI spectra

This is the Gaussian Decomposition software described in Nidever et al. (2008) following the algorithm from Haud (2000).  While it was designed to be used for HI spectra, it can be used for other types of data like Halpha spectra.

The software was originally run on the Leiden-Argentine-Bonn all-sky HI survey and there are hard-coded settings that need to be modified for each data-set.  My plan is to modify the code at some point in the future to be more general and to allow these values to be more configurable.  Currently, I've been making data-set specific copies of five of the programs and modifying the settings.  The programs are:
- gdriver.pro: Main driver program.  There are a few hard-coded defaults that you might want to change.
- gloadspec.pro: Loads the cube.  You might want to make some modifications to how it loads the data.
- gincrement.pro: Increments the current position.  Depending on the data the step-size will change.
- parcheck.pro: Checks if Gaussians have "bad" parameters.  The thresholds for "bad" might need to be tweaked.
- setlimits.pro: Set limits on all of the Gaussian parameters (height, velocity, width). 
- hinoise.pro: The program that calculates the noise in each spectrum.  You should set which velocity range to use.


## Example

The main program is gdriver.pro.  I normally create a small IDL batch script to run segments of a cube.  An example one is provided using a small GASS cube downloaded from https://www.astro.uni-bonn.de/hisurvey/gass/ using these parameters.

```
l = 295.0 deg
b = -41.0 deg
width in l = 1 deg
width in b = 1 deg
```

`gass.in`

spawn,'echo $HOST',host
print,'RUNNING THIS PROGRAM ON ',host
@compile_all
; The cube is [551,1251,1137]
gdriver,lonr=[0,1250],latr=[0,550],file='gass3.dat',/noplot,$
        btrack=btrack,gstruc=gstruc,/backret,savestep=100,subcube=2

## Plotting the Results

The repository includes a plotting routine called `ghess.pro`.