# README: Interactables
### Author: T. Hankins, Date Last Modified: 240718

## Quick Access

- [energy_calibration](#energy_calibration)
- [perform_thickness_measurement](#perform_thickness_measurement)

## energy_calibration

- Interaction with this executable is only needed if the fit to the energy spectrum returns incorrect results, such as misidentifying or missing peaks.
- `energy_calibration` makes use of the ROOT `TSpectrum` class to automatically identify peaks. The algorithm relies on a detection threshold, defined as a value between 0 and 1, which determines which features in a spectrum can be identified as peaks. A value of 0 means "anything goes", while 1 limits to the maximum bin in the spectrum.
- A common value for the detection threshold is anywhere between 0.15-0.45; at the time of writing, the value is set to 0.35 (line 110). This value is set based on the structure of the true source spectrum as well as statistical effects and may need to be modified. If a peak is missing, the threshold needs to be lowered, and if too many peaks are identified, the threshold needs to be increased.
- If the threshold is changed, don't forget to `make install` to compile the changes.

## perform_thickness_measurement

### General Usage

- Because the material that is being measured likely changes between targets, the source code for `perform_thickness_measurement` has to be modified and recompiled. The only thing that needs to be changed is the `CycSrim` target material, and the line for this is the very first one in the body of the executable (line ~20). *Only change the* `CycSrim` *material; nothing else about the target in the code needs to be changed.*
- `perform_thickness_measurement` also has an additional option that the other executables do not have, and that is option `c`. This defines the coarseness of the thickness map, or the position discretization.
    - The option is given in millimeters and defines how wide each histogram XY bin relates to physical space. As an example, the DADL physical dimensions are 20 mm x 20 mm, so a `c` value of 1 would divide the face into 1 mm x 1 mm segments.
    - Individual energy histograms are filled based on the position, and all spectra are written to the output file, in addition to the position map. The fits are also present on these, so any issues that arise can be reasonably diagnosed.
        - A frequent issue that arises is missing the highest energy peak in <sup>228</sup>Th in position bins with low statistics.
    - The position resolution of a ~8 MeV $\alpha$ particle is about 0.5 mm, and thus a warning was added to the program. If the value of `c` is less than 0.5, the program will print a warning.
    - The naming and consolidation scheme of the program also relies on the discretization of the DADL. A requirement of `c` is that it evenly divides into 20. Therefore, values of 0.1, 0.2, 0.25, 0.4, 0.5, etc. are all numerically valid, but 0.3 is not.
    - The runtime of the program scales as N<sup>2</sup>. Additionally, the fitting routine struggles with lesser statistics, so the two go hand-in-hand. It is recommended that the program is run coarsely to begin with, only increasing the discretization once confirming that the more coarse analog worked properly.

### Modifications
- As mentioned in the previous subsection, the program struggles with lesser statistics. As with `energy_calibration`, `perform_thickness_measurement` also uses `TSpectrum` for identifying peaks. A threshold value exists for this executable (at the time of writing, line 186). If problems are encountered, this value may be varied through trial and error to obtain a better result. However, generally, the problem arises for spectra with limited statistics.
- Explicit treatment of low statistic spectra has been included in this program (lines 181 and 182). These values/treatments can also be changed as desired.
- Two versions of the thickness map exist - one that considers all valid fit bins, and another that filters bin content in comparison to the rough average of the target to make quick viewing easier. When this second histogram is filled, the individual contents are added only if they are within a certain window relative to the rough average. The thresholds for this acceptance can also be modified (at the time of writing, line 273) and may be necessary for a highly non-uniform target.
- If these are changed, don't forget to `make install` to compile the changes.
- [statement about various modifications - 228th number of events cutoff due to thoron emanation, threshold, etc.]