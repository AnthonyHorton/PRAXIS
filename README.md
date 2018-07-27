# PRAXIS Viewer

This programme is a quick look viewer for PRAXIS data. If produces a reconstructed IFU image to help
with target acquistion.

## Requirements

### Python libraries (essential)

In addition to Python this programme requires the `numpy`. `astropy` & `matplotlib` Python
libraries. These can be installed with the operating system's package manager, `conda`, HomeBrew or
```
pip install -r requirements.txt
```

### Environment variable (recommended)

It is possible to reduce the amount of typing required by setting the `$PRAXIS` environment variable
to the path of the directory that contains the PRAXIS exposures, e.g.
```
export PRAXIS=/run/user/1234/gvfs/smb-share:server=127.0.0.1,share=hxrg/H2RG-C001-ASIC-SC2-B-010-JWST/FSRamp
```
Add this to your `~.bashrc` file to automatically set `$PRAXIS` in every new shell.

## Usage

To run the PRAXIS viewer on the most recent image:
```
./PRAXIS_viewer.py
```
A window will appear with the reconstructed IFU image and the spectrum from the brightest fibre. To
continue close the window. When run without sky subtraction the IFU image is normalised by
subtracting the flux value from the faintest fibre and dividing by the brightest. **Note, automatic location of the latest image requires the `$PRAXIS` environment variable to have been set.**

To run on the most recent image with sky subtraction:
```
./PRAXIS_viewer.py --subtract 20180727103250
```
This will subtract the image taken at a datetime of `20180727103250` from the more recent image
before doing the IFU image reconstruction. When run with sky subtraction the IFU image is normalised
by dividing by the brightest fibre, no zero point correction is done.

If `$PRAXIS` is not set, or if you want to use a sky image that is not in the PRAXIS raw data
directory then the full path must be given:
```
./PRAXIS_viewer.py --subtract ~praxis/stuff/sky_10_seconds.fits
```

To run on an earlier image specify its datetime:
```
./PRAXIS_viewer.py --filename 20180727121602 --subtract 20180727103250
```
