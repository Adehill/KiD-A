# KiD-A intecomparison testcases

Phase 1 of the KiD-A intercomparison is be based on the following:

warm1 (1D) - simple 1D updraft case described in Shipway and Hill (2012)
2D stratocumulus (Sc 2D) - based on the case 1 from the 8th International Cloud modelling workshop (2012)

In both cases aerosol are assumed to be soluble ammonium sulfate particles. The initial aerosol distribution is assumed to be a single mode lognormal distribution with a lognormal geometric mean diameter = 0.08 * 10-6m and a log standard deviation = 1.4. The initial aerosol number concentrations (Na) are defined in the test case descriptions below. These parameters are defined in the test case namelists. If a participants model does not include aerosol, please set the initial cloud drop number concentration (Nd) to the initial Na defined below.

### 1D case (warm1) - no aerosol processing

This is the initial cloud microphysics test, which will be run with a fixed aerosol or Nd depending on scheme. The vertical velocity setup for this case is based on the warm1 case, where maximum vertical velocity is set in the namelist to 2 (W2) and 3 (W3) m s-1. For both W2 and W3, we require a simulation for Na or Nd (depending on the scheme) = 50, 150 and 300 cm-3. For each updraft and Nd, you are asked to simulate the following:

cond-evap case - only aersol-activation (dependent on scheme), condensation and evaporation are switched on, i.e. sedimentation, collision-coalescence, breakup, etc. are switched off. Hopefully, this case will show that all schemes produce the same amount of liquid water at about the same rate.
Precipitating case - all cloud microphysics processes are switched on; however, no removal or replenishment of aerosol and no in-cloud aerosol processing (in schemes that have this functionality) should be simulated.
Templates for a bin and bulk namelist required to run the 1D case, i.e. namelists/kida_1D_bin_template and namelists/kida_1D_bulk_template.nml, are included in KiD Source code download .

### Sc 2D case - with and without aerosol processing

For the 2D SC case, we require a simulation to be run with maximum vertical velocity (wctrl) = 0.25 and 1 m s-1. Once again wctrl is set in the namelist for the case. Na (or Nd depending on scheme) = 50, 150 and 300 cm-3 for the following configurations:

cond-evap case only aersol-activation (dependent on scheme), condensation and evaporation are switched on, i.e. sedimentation, collision-coalescence, breakup, etc. are switched off. Hopefully, this case will show that all schemes produce the same amount of liquid water at about the same rate.
Precipitating case - no removal or replenishment of aerosol and no in-cloud aerosol processing.
Precip case with aerosol processing - this will include the removal of aerosol by activation scavenging, replenishment of aerosol by evaporation and in-cloud aerosol processing.
In all cases, we apply the following spin-up
- 0 - 10 mins - During this period maximum supersaturation is restricted to 0.5%, all velocities are 0 and precipitation processes are switched off.
- 10 - 20 mins - supersaturation is no longer restricted and velocities gradually increase to the wctrl but precipitation processes are switched off
- 20 - 30 mins - run for 10 minutes with standard wctrl (0.25 or 1 m s-1) with precipitation switched off
- 30 - 3 hours 30 mins - run with precipitation switched on
The reason for this 30 minute initialisation is because the case starts with a supersaturation, which can cause some models some problems. Also, the 30 muinute initialisation permits the development of the cloud drop distribution in the bin and super-droplet schemes, prior to precipitation processes.
To set-up this initialisation, we have added the following switches to the KiD model and the namelists:

no_precip_time = 1800.0 ! The time from beginning of simulation with no precipitation processes or sedimentation
smax_limit_time = 600.0 ! The time from the beginning of the simulation with in which a maximum supersaturation is applied
smax = 0.5 ! the maximum percentage supersaturation during smax_limit_time
Participants need to add these switches to the appropriate parts of the microphysics code so that the recommended initialisation is used. The switches are stored in the namelists module.

Templates for a bin and bulk namelist required to run the 2D Sc case, i.e. namelists/kida_wmoSC_2d_bin_template.nml and namelists/kida_wmoSC_2d_bulk_template.nml, are included in KiD Source code download . Finally, to build the executable for the WMO case use CASE=WMO_CASE1 in the make statement (see the documentation).

## Diagnostics

The example diagnostic outputs (see links above) present standard diagnostic output for the KiD model. Once a microphysics scheme is coupled with the KiD framework, we believe that all cloud microphysics fields accept the process rates and precipitation rates will be automatically output. The participant will need to add a save_dg call to their scheme to store the process and precipitation rates that are output from their scheme. Details about adding diagnostics to a scheme are presented in the documentation, with examples available in the Morrison, Thompson and TAU schemes that are available with the download.

A diagnostic requirement for all schemes (bulk and detailed bin microphysics) is the provision of autoconversion, accretion and rain evaporation process rates. This requires the definition of a boundary between cloud and rain drops. For these diagnostics, we are suggesting that rates be provided for a 20 and 32 micron (radius) cut-off, where all bins below are described as cloud and all bins above are described as rain. We accept that this is artificial; however, such a cut-off permits comparison with the bulk schemes. If you are unable to provide such partitioning and process rates, or if you want to provide a best estimate partition, please do so but also provide a definition of the partition.

In addition to cloud microphysics fields, we also require the following aersol diagnostics from schemes that include aerosol and aersol processing.

- Aerosol mass activation rate
- Aerosol converted from cloud droplets to rain drops via autoconversion
- Aerosol converted from cloud droplets to rain drops via accretion
- Aerosol lost to the surface via precipitation removal
- Aerosol regeneration rate via droplet evaporation
Finally, for the bin microphysics schemes, we have requested that bin diagnostics be included so that we are able to analyse the particle size distributions. This can lead to very large diagnostic files; thereofre, we kindly ask that zipped files be submitted.

Once diagnostics files are available and participants are ready to submit, participants should contact Adrian Hill (adrian.hill@metoffice.gov.uk) and Zach Lebo ( zlebo@uwyo.edu ).