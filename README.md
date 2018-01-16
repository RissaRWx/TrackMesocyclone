# TrackMesocyclone
NCL code to detect and track midlevel mesocyclones in model outpuit

## Main Function (track_meso.ncl)
This is a callable NCL function that reads in a time-varying updraft helicty (UH) array and uses it to detect and track mesocyclones in the model output. The user provides a minimum UH value required to be considered part of a mesocyclone, as well as a minimim number of contiguous gridpoints that must exceed this UH value to be counted as a mesocyclone. The algorithm also requires an estimation of storm motion, typically computed as the Bunkers right-deviant vector (Bunkers et al. 2000) from a pre-stom model sounding.

## Version Notes
This code was written using NCL v6.3.0, but with a beta code of some added functions for v6.4.0 (contributed.ncl_beta_640). If v6.4.0 or later is used, this code and the line calling it (Line 1) are unnecessary.

### Maintainer and Contributer
  - [Larissa Reames] (larissa.reames@noaa.gov)
