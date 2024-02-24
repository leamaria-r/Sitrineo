# Our work

In this project, we tried to improve the existing code by added a few features that will be described below.

## Source files

### Mraw.cxx

* Sitrineo Align :  
We added the gaussian fit for the fourth plane and corrected a few typos (l.4218-4341 and l.18290-18400)  
We then added the uncertainty on the gaussian fit.  
We also automatized the alignment procedure by updating the config file with the new position values

* SitrineoCumul and SitrineoAnalysisFromHits
For the cuts (around l.17650), we defined the positions on z-axis of planes 1 and 4 using the positions in the config,
choose a Chi2 using the plots obtained in the analysis with SitrineoCumul(), and defined a geometrical cut based on the hits spot taking into account the diffusion and multiple scattering on the planes in the back.

More generally we linked the planes positions and magnetic field magnitude (that were hardcoded before) to the config file

### DSetup.cxx :
We added BFieldMagnitude as a tracker parameter

### DGlobalTools.cxx : 
We corrected a typo in the muon mass

## Include files

MRaw.h and DSetup.h were modified to match the modifications of the source files
