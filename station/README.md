# Target Testing Station Run Guide
### Author: T. Hankins, Date Last Modified: N/A

- This guide serves as an introduction and reference for using the ANTLER target testing station. The publication discussing this station and measurement technique is here: https://doi.org/10.1016/j.nimb.2025.165788.

## Quick Reference Checklist
- The detector **must be unbiased** before opening the chamber; do not apply bias unless the chamber is light-tight and the pressure is less than 10 mT.
- The full-size roughing pump and the pumping station should **never** be open to the chamber at the same time. 
- Do not open the pumping station to the chamber if the foreline pressure exceeds 500 mT, preferably 200 mT.
- Keep vent cycles as short as possible to reduce subsequent pumpdown time.

## Hardware
- An image of the full target testing setup is given below. The station is comprised of a shoebox vacuum chamber attached to a full-size roughing pump and a HiCube 80 Eco turbomolecular pumping station. The data acquisition hardware, detector bias module, and preamplifier power supply are stored in the blue rack next to the setup.

<!-- <center>
<img src="source/full_station.jpg" width="400">
</center> -->



### Detector
- The detector used in the testing station is a dual-axis duo-lateral (DADL) position-sensitive silicon detector. See the following publications for supplementary information:
    - https://doi.org/10.1016/j.nima.2009.11.053
    - https://doi.org/10.1016/j.nima.2020.164674
    - https://doi.org/10.1016/j.nima.2023.168130
- The DADL has six output signals: four face contacts (front face 1 and 2 (F1, F2), back face 1 and 2 (B1, B2)) and two guard rings (front (FGR) and back (BGR)). The four face contacts are used in the determination of energy and position.
    - Front contacts produce negative signals, while back contacts produce positive signals.
- The detector signals are passed through the rear feedthrough and directly into a preamplifier box containing four ~45 mV/MeV charge-sensitive preamplifiers. LEMO inputs and outputs on this box are appropriately labeled. Images of the feedthrough and signal identification diagram are given below.

<!-- <center>
<img src="source/feedthrough.jpg" height="200"><img src="source/feedthrough_label.jpg" height="200">
</center> -->

- For the detector in use at the time of writing, the front face is biased to -65 V and the front guard ring to -58.5 V. Both the back face and back guard ring are held at ground. The measured leakage currents for both are consistently measured to be less than 0.1 uA each. If this is not the case, **immediately bias down and inform Alan McIntosh.** A bias curve for this detector is provided below for reference.
    - If the detector is ever exchanged, the bias curve will need to be reproduced and any statements of bias voltages throughout this guide will need to be replaced with values from the new detector.

<!-- <center>
<img src="source/250709_BiasCurve.png" height="300">
</center> -->

- Inside the chamber, the DADL is mounted to an aluminum mount that can be moved. In general, the mount shouldn't need to be moved unless the measurement demands it. If moved, *always* check the four face signals before and throughout the pumpdown procedure by connecting the LEMO cables fed into the 3316 into an oscilloscope (see [Rack](#rack) for discussion). DADL signals are known to disappear if the Molex connector (directly off the DADL) is moved to a poor position; this requires an additional vent cycle if it's not noticed until the data acquisition stage.
    - With all four signals fed into the oscilloscope, it is very obvious when a signal is missing. To recover the signal, it is enough to adjust the Molex connector until the signal returns. Exercise caution when performing adjustments.



### Chamber
- TODO: detector mount, again; dimensions, statement of position accuracy depending on source-detector distance.
- The chamber has a single actuator arm that moves an aluminum blocker in front of the alpha source. An image of the arm is given below; note the two faint marks. The mark toward the outside of the arm is the position for *blocking the alpha source*, while the inner mark is the position for *unblocking* (as an example, the arm as shown below is in the unblocked position). 

<!-- <center>
<img src="source/actuator_arm.jpg" width="400">
</center> -->

- When moving the arm to a new position, hold the chamber in place, as reasonable effort is required. Move the arm slowly so that vacuum will not be seriously affected (see [Vacuum](#vacuum)).
- The arm is in the correct position when the corresponding line is flush with the knurled brass cylinder attached to the flange.
    - When moving the arm to a new position, get the position correct to within a few millimeters.
- When not taking data, keep the actuator arm in the blocking position to prevent unnecessary detector exposure to source.



### Vacuum
- Performing vacuum cycles and alternating between the pumps used to achieve high vacuum is the most common physical operation performed on the testing station short of opening the chamber. **Read this section carefully, as a lot of damage can be done very quickly if not careful.**
- As mentioned previously, the vacuum system for the testing station is comprised of a full-size rotary vane roughing pump and a turbomolecular pumping station. The full-size roughing pump is used to quickly achieve rough vacuum, after which the pumping station takes over.
    - The pumping station is comprised of a diaphragm backing pump and a small turbo. Conventionally, a station such as this would be used independently to achieve high vacuum, but due to the low pumping speed and high base pressure of diaphragm pumps (further, the size of the pump itself), an independent, full-size roughing pump was used.
    - In addition to accelerating the pumpdown procedure, this choice significantly reduces the number of spin-up and spin-down cycles for the turbo.
- Vacuum for the chamber is monitored using a thermocouple gauge. An image of the controller/readout is given below.

<!-- <center>
<img src="source/gauge_readout.jpg" width="400">
</center> -->

- The chamber connections to the main roughing pump and the pumping station are controlled separately using two manual isolation valves, shown in the image below; the roughing pump is controlled using the tan valve, while the pumping station is controlled using the black valve.

<!-- <center>
<img src="source/valves.jpg" width="200">
</center> -->

- An image of the turbomolecular pumping station is given below. The turbo is the small cylinder on top of the red chassis; the diaphragm pump is located inside the chassis. The main power is controlled using the green switch on the front. The tan panel is the digital control unit (DCU) which controls and monitors the status of the turbo.
    - The nominal rotational speed of the turbo is 1500 Hz with an base current draw of around 1.15-1.3 A. If the current is much higher than this value when it shouldn't be (e.g., pulling on the vacuum chamber when it's fully evacuated), **immediately contact Alan McIntosh.** Note that the current draw does depend on the load; spinning up, the turbo can exceed 3.0 A, and spinning down, the current may be close to zero, but otherwise should be near above.
    - The instrument is designed to automatically shutdown if an error occurs to reduce chance of damage. If the turbo is not running, contact Alan McIntosh.

<!-- <center>
<img src="source/turbo.jpg" width="200">
</center> -->

- **Under no circumstances** should the main roughing pump and the pumping station be open to the chamber at the same time; the high vacuum pulled by the pumping station can cause backstreaming of oil from the roughing pump, which can deposit inside the chamber and on the detector.
- While the general rule-of-thumb for using turbomolecular pumps is to open the valve once the foreline pressure is less than 100 mT, this turbo can be opened at a higher pressure. The key difference is the use of a manual valve, which is finely adjustable (in constrast to a beamline valve, which is either fully open or closed and quickly alternates between these states), as well as the small chamber volume.
    - To do this, monitor the current draw on the DCU and **slowly** open the valve to the turbo once a reasonable foreline pressure is attained and the valve to the full-size roughing pump is closed. The current draw will increase as the turbo further evacuates the chamber; use this to determine whether to further open or re-close the valve. Do not open the valve further than what corresponds to a current draw of ~1.4 A.
    - **Absolutely** do not attempt to open the turbo valve until the foreline pressure is less than 500 mT; if this cannot be achieved using a roughing pump, the chamber may be leaking. Further, if the roughing pump is able to continue evacuating the chamber in excess of 1 mT/second, let the roughing pump do the work until it is incapable.
- Getting a feel for pumping down the station will take a couple of cycles. In general, never pump down faster than ~50 torr/s when beginning at atmosphere; this limit may need to be lowered for fragile targets. **Breaking a target is possible and has occurred when pumping down too quickly.** Similarly, be careful when venting.
- Try and keep vent cycles as short as possible. Time spent at air will gradually lengthen the time it takes to pump back down due to an increasing need to outgas.
- Pulling the actuator arm out of the chamber creates a mild, temporary leak; the severity of the leak is related to how quickly the arm is pulled out of the chamber. When pulling the arm, go no faster than a centimeter per second so that the turbo isn't strained.
- Finally, several manuals for the pumping station components are in a black binder next to the safe.



### Rack
- The rack (image below) holds the NIM bin, VME crate, and Topward power supply for the detector bias, data acquisition hardware, and preamplifier power, respectively.

<!-- <center>
<img src="source/rack.jpg" width="400">
</center> -->

- Power for the preamplifiers is supplied using a Topward 6306D dual-tracking DC power supply. The voltage for each output is set to 12.0 ± 0.1 V; because of the orientation of the dual banana to BNC adapter, the left channel supplies +12 V, while the right channel supplies -12 V. A typical current draw is < 0.1 A per channel.
    - If the current is higher than this (often while being incapable of supplying the full 12 V on one or both channels), contact Alan McIntosh.
- DADL bias is supplied using a singlewide Tennelec TC-953 dual-channel high-voltage power supply, which is housed in a standard NIM bin. The NIM bin power switch is located on the front right panel.
- DADL signals are recorded using a Struck Innovative Systeme (SIS) 3316 14-bit VME digitizer. This module is housed in a W-IE-NE-R VME crate in the top of the rack. Each of the face signals from the DADL are fed into a single channel of the 3316. The crate is connected to the acquisition computer using a SIS3104-2 crate controller and an orange fiber optic cable. **Do not bend this cable.**
    - The crate power switch is located on the front panel. When turning the crate off, the cooling fans are programmed to run for an additional minute before turning off.
    - Channels 1-4 are, in order: F1, F2, B1, and B2. **Do not change the signal order.** For completeness, the 3316 trigger out is fed into channel 8. Channels 9-16 do not appear to work on the 3316 currently in use.



### Hardware: Operation Steps
- In general, most concerns when working with the testing setup are addressed through the three V's: **visual, vacuum, and voltage**. When working with one of these, remain keenly aware of the status of the other two.
    - **Visual**: detector exposure to light.
    - **Voltage**: application of bias to the detector.
#### Evacuating
- Assuming that the chamber is at atmosphere and nothing is powered on in the rack:
    1. If the DADL mount has been moved since the previous signal check, unplug the four DADL signals from the 3316 and plug them into the oscilloscope. Turn the Topward power supply on, then slowly turn both current knobs clockwise until 12 V is applied to both outputs. Verify that signals are present on the scope.
    2. Load the vacuum chamber with a target to measure (if not collecting reference data). 
    3. Ensure that the valve to the roughing pump is fully closed. The valve to the turbo pump should already be closed. If not attached, attach the KF-10 vacuum tube from the full-size roughing pump to the tan valve using an O-ring with a plastic centering ring. If the chamber isn't closed, close it. Use a couple of screws to help align the lid; tightening is not necessary.
    4. Plug in the full-size roughing pump to begin pulling vacuum on the foreline. Then, **slowly** open the tan valve to begin evacuating the chamber, using the thermocouple readout as a guide for how far to open the valve. As the chamber pumps down, further open the valve. If the chamber was open for longer than a couple of minutes, do not be surprised if the evacuation stalls several times at various intermediate pressures.
    5. Once the chamber is sufficiently evacuated for the turbo, close the valve to the roughing pump (do not overtighten). **Slowly** open the valve to the turbo, using the *DCU current* as a guide rather than the thermocouple readout. Once the current rises slightly, check the thermocouple readout to verify that the vacuum is improving. Again, as the chamber pumps down, further open the valve while closely monitoring the current draw.
    6. At full vacuum, unplug the full-size roughing pump. If the DADL signals were plugged into the oscilloscope, again verify that they are present, then return the outputs to the 3316 without altering the signal order; otherwise, turn the Topward on as described in step 1. Power on the fan underneath the NIM bin, the NIM bin itself, and the VME crate.
    7. In the NIM bin, check that the "NEG" lights on both channels of the Tennelec TC-953 are lit (if not, contact Alan McIntosh). Turn on both channels if they are not already on. Switch the LCD panel viewing option underneath each panel to "V" if they aren't already set to this. Both channels should read 0.0.
    8. If the chamber is light-tight (it should be when at vacuum and the viewport cover is taped to the chamber), begin turning both voltage knobs clockwise to apply bias. Do not increase the voltage faster than 5 V/s. At 30 V, switch the LCD panels to "uA" to check the leakage currents; if they're okay (as compared to the bias curve provided in [Detector](#detector)), switch back to "V" and continue biasing up.
    9. Stop increasing the bias on the lower channel (front guard ring) at -58.5 V. Continue increasing the upper channel (front face) until at -65 V. Check the leakage currents once more, and if they're still good, proceed with data acquisition.

#### Venting
- Venting the chamber from high vacuum is largely the evacuation procedure performed in reverse. Assuming that everything is powered on and the chamber is at high vacuum:
    1. Ramp down the bias on the detector by turning the voltage knobs on the Tennelec counter-clockwise until both channels are zeroed. The Topward current will fluctuate during the ramp down. Do not decrease the voltage faster than 5 V/s. Checking the leakage currents is not necessary.
    2. If done collecting data for a prolonged period, turn off the NIM bin, then the fan underneath the bin. Close out any open acquisition software on the computer, then turn off the VME crate. On the Topward, turn both current knobs counter-clockwise until zeroed, then turn the unit off. Otherwise, leave all units on and continue with step 3.
    3. Close the valve to the turbo pump; the valve to the roughing pump should already be closed and the roughing pump itself should be unplugged.
    4. Disconnect the KF-10 vacuum tube from the roughing pump valve; this will vent the roughing pump foreline. Leave the tube disconnected.
    5. **Slowly** open the tan valve to begin venting the chamber. As the chamber gets closer to atmosphere, re-close the valve to slow venting; use the thermocouple readout as a guide. At atmosphere, the chamber can be opened.



## Software
- The target testing station uses standard Cyclotron Applications ROOT-based GUIs to control the acquisition hardware, collect data, monitor this collection and its status, and perform preliminary analysis. The software used to control the hardware and collect data is called the transport manager, abbreviated "Tmr". The software used to monitor collection by means of diagnostic plots is called the analysis manager, abbreviated "Anl".
- Versions of Tmr and Anl are created for each experiment, which are identified by their start date. For versions that do not correspond to an experiment, a fake date is generally chosen. At the time of writing, the control and acquisition computer used is `cycfe6` and the "experiment" is `063123`. The user for the computer is `sjygroup` and the password is `sjysjy`.
- An image of both Tmr (left) and Anl (right) is given below.

<!-- <center>
<img src="source/software.png" width="600">
</center> -->

### Tmr
- Tmr is controlled from `063123_frontend/`. It can be activated by `cd`-ing into the directory and running `Tmr sis` once the VME crate is turned on. An example process would look like:
    ```
    cd 063123_frontend/
    Tmr sis
    ```
    Note that Tmr will crash if the VME crate is not powered on.
- New runs are started by clicking "Start ACQ". When this is clicked, a "Run Information" dialog box will appear (image below). The next run number, if data is being written to disk, is shown in the run number area. The type of data that is being collected is listed below that; this will always be "Source" for 060925. Finally, the option to write data to disk is given by the checkbox labeled "Disk Output".
    - If "Disk Output" is checked, additional areas to fill out will appear. Do not change the "Directory" or "Filename" fields, as these are automatically generated. The "Supervisor" field is for the user's initials. Finally, the "Comment" field provides an area to put notes for the given run. Include the target type, alpha source, and DADL bias voltages for completeness.
        - When typing in the box, keep the mouse hovering within the box; otherwise, text will not be written.
    - It is considered good practice to start a run without writing to disk to check the summaries in Anl, then begin writing after verifying everything is working. Anl will collect data from Tmr independent of Tmr's write state.

<center>
<img src="source/dialog.png" width="350">
</center>

- When finished collecting data, clicking "Stop ACQ" will stop the acquisition and, if writing to disk, close the current file. After a few seconds, Tmr will be ready to run again.
- Use the "Exit ROOT" button to close Tmr rather than clicking "X". If applicable, close Anl before closing Tmr to prevent a segmentation fault (see [Known Issues](#known-issues))
- Advanced settings (including, but not limited to the "SIS3316 Control" tab) should not be modified without consulting Alan McIntosh or Kris Hagel first.



### Anl
- Anl is controlled from `063123/`. It can be activated by `cd`-ing into the directory and running `Anl` once an instance of Tmr is running. An example process would look like:
    ```
    cd 063123/
    Anl
    ```
    Once the program is running, Anl will begin receiving data when the "Begin" button is clicked.
- A number of plots are already built into the `063123`'s Anl; the ones that are most important are:
    - `h_DADL_X_ch` - energy spectrum, in channels, for the front (F), back (B), or sum (S) of DADL signals.
    - `h_FrontVsBack` - plot of front versus back face energies. This should be a tight correlation with a slope close to one.
    - `h_DADL_X/Y` - 1-D histograms of particle position for the front (Y) or back (X) face.
    - `h_DADL_XY` - 2-D scatter of particle position. Ideally, this is rectangular for reference data and ovular for target measurements.
    - `h_DADL_ENvsM` - 2-D scatter of face energy versus position (front-Y, back-X). Ideally, this is flat as a function of position, but is often not.
- Clicking "Clear Arrays" will clear any previously collected data, if needed; this has no effect on the raw event information written to disk from Tmr.
- Use the "Exit ROOT" button to close Anl rather than clicking "X".

### Known Issues
- This version of Anl is known to be susceptible to crashing, especially when viewing plots while data is being sent from Tmr. Pausing the data transfer between Tmr and Anl by clicking the "pause" button in Anl helps significantly - just don't forget to resume afterward, if desired.
- If either Anl or Tmr crashes while running and does not automatically close, open a new terminal and run either `pkill Anl` or `pkill Tmr`. Further, if Tmr is closed, either via a crash or manually, while Anl is open, closing Anl will lead to a segmentation fault that often requires `pkill`.
    - In some instances, `pkill` will not be sufficient to stop the process. In this case, run the following:
    ```
    ps aux | grep "X"
    ```
    where X is either Tmr or Anl. This will return a list of active processes with the X specifier. In the following terminal output, find the corresponding process. It will look something like:
    ```
    sjygroup 3455482 6.4 2.3 699136 372752 pts/5 Sl+ 09:55 0:01 Anl
    ```
    What matters most is the user (leftmost) and the process ID (second to leftmost). With the process ID, run `kill -9`. For the ongoing example, this would look like:
    ```
    kill -9 3455482
    ```
    This should handle the faulting program.

### Accessing Data for Further Analysis
- To be added...