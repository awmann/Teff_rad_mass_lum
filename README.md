Teff_rad_mass_lum
=================

Note: the original version of this code was based on the Mann et al. (2013c) calibration. 
While this calibration is fine, it's somewhat obsolete as it failed on extremely metal-rich (>+0.3), metal-poor (<-0.4), or very cool (~<M5) stars. 
So I've replaced it with some more accurate relations. Additionally, this code can use simple photometry to get out the effective temperatures.  
If you make use of this code cite Mann et al. 2015 (http://adsabs.harvard.edu/abs/2015ApJ...804...64M)

Code is split into different components:
- Phot_Teff.pro, which gives out Teff values from input photometry.   
- Rad_teff.pro, which gives out a radius given a Teff (and [Fe/H] if available)  
- Rad_MK, which gives a radius from an absolute K-band magnitude (and [Fe/H] if available)  
See the individual code for syntax/instructions. Bolometric relation formulae are coming soon!

All programs take errors on the given values as well, if no error is provided it assumes the error is 0 and only returns errors based on the calibration error. 


Important: this assumes V magnitude is on the Johnson system, I is on Cousins, rz are on SDSS, and JHK are on Twomass system. Keep in mind that many sources of photometry claim to be on these systems but are not, or at least have larger errors than purported.\
For example, here are APASS gri magnitudes (which are presumably the same as SDSS) compared to the synthetic magnitudes used for our calibrations as a function of V-J color. As you can see, APASS r does appear to match our system, 
but the errors are underestimated (APASS notes this in their documentation as many times they have too few measurements for an accurate error estimate). 
A bigger concern is the fact that g and i show 'color' terms, suggesting the system throughput is not exactly the same as SDSS. This is part of the reason we don't have such relations.

![](images/apass_comp2.png?raw=true)
