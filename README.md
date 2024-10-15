# LATTICE_PLANE_ANGLES

This script can be used to help with indexing of diffraction patters. 
Input parameters are controlled in the area marked by the ######'s at the top of the script.

----------------------
## Inputs:
 - lattice parameters of candidate structure (obtained e.g. from powder diffraction refinement)
 - type of Bravais lattice (currently only primitive, body-centred, and face-centred are implemented)
 - measured angle between two diffraction spots and the ratio between the lengths of the corresponding vectors
 - some other control parameters, which can be used to fine-tune the output (default values should be okay for most cases)
 
 ## Output:
 - miller indices of candidate lattice planes
 - angle and ratio between the two planes (note that d1 is the spacing for (h1,k1,l1) and vice versa)
 - calculated viewing orthogonal direction
 - candidate combinations of planes are ranked based on the squared deviation from the input angle and ratio, weighted by the respective tolerance (the 'score') 
----------------------

## Notes:
 - this code currently only accounts for systematic absences due to Bravais lattice, so there are additional space group-dependent reflection conditions (see ITC vol A, table 3.2)
 - for electron diffraction you minght have to consider double diffraction, so there might be spots that appear to break the reflection conditions.
   If a spot can't be indexed with I- or F-centering, it might be wirth trying P-centering (since that shows all forbidden reflections).
 - the code can calculate the lattice spacings for a candidate lattice plane, but these values might differ significantly from measured values (e.g. due to measurement calibration or deviations in lattice parameters)
   so it's usually more reliable to just compare the ratios and angles.
