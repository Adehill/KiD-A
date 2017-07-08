## Introduction

A major uncertainty in numerical weather prediction (NWP) and climate prediction is the response of clouds and precipitation to changes in aerosol concentrations. Previous intercomparison work has shown that the simulation of precipitation for a given cloud drop number concentration (Nd) and the response of precipitation to changes in Nd in cloud resolving models (CRMs) and idealised kinematic frameworks is very sensitive to both the physical representation of microphysics and numerical complexity of the scheme cloud microphysics (e.g. ADD REFERENCE). The purpose of this microphysics intercomparison project is to compare detailed size resolved and bulk parametrised microphysics schemes to understand how they simulate aerosol-cloud-precipitation interactions. The project consists of two phases:

- Phase 1 - Employ the Kinematic Driver model (KiD) to compare microphysics schemes when considering the a) the simulation of warm rain b) the response of warm rain to changes in Nd and c) in-cloud processing of aerosol 
- Phase 2 - investigate a, b and c from phase 1 with a dynamic case in CRMs, where the dynamic feedbacks are considered 

This site presents the details for the Phase 1 of the project, KiD-A 

## Overview of the KiD-A project
The overarching aim of KiD-A is to use the Kinematic Driver model (KiD) to compare detailed and bulk microphysics schemes in a dynamically consistent framework without the complication of dynamic feedbacks and numerical issues, that have been experienced in previous intercomparison projects. The main aims are:
1. Undertake the first kinematic intercomparison of detailed microphysics schemes, i.e. size resolved bin microphysics schemes, superdroplet schemes and 2D aerosol-cloud schemes 
This will be a sanity check to make sure that the schemes, which are consistently used to develop simpler bulk schemes produce similar results, when forced with the same dynamics. 
Tests will be performed for a range of initial aerosol concentrations and 1D and 2D kinematic cases 
Tests will exclude in-cloud aersol processing and focus solely on the precipitation processes, timing and amount 
2. Examine and compare in-cloud aerosol processing from both detailed microphysics and bulk microphysical representations. 
This stage of the project will compare the detailed and bulk microphysics schemes that can include in-cloud aersol processing. 
This will be based on only the 2D kinematic cases. 
The results from this stage will be used to benchmark present modeling capability when considering aerosol-cloud-precipitation interactions. 

## Kinematic Driver Model (KiD)

The KiD model is the basis of Phase 1 of the KiD-A project. The source code is available on this site and documentation contains all of the details about how to build and run the KiD model as well as details about haow to add a microphysics scheme and include diagnostics. This version of the code includes old versions of the Morrison microphysics and Thompson microphysics schemes (both single and double moment rain schemes) and a dummy subroutine that contains details about how to obtain and include the Tel-Aviv University (TAU) bin microphysics schemes.

It is expected that participants will add their microphysics scheme to the KiD model. This will involve writing an interface that couples the generic KiD prognostics to the microphysics variables. The documentation describes how to add a scheme to the KiD framework, while the Morrison, Thompson and TAU provide schemes are provided as examples of adding a microphysics code to the KiD model. In general, there is no need to change the microphysics code so that it will work with the KiD model. Further, diagnostic calls may be required in the microphysics code so that process rates and precipitation rates are output. Otherwise, all required diagnostics will be automatically output once a scheme is successfully coupled. The required diagnostics are discussed below.

Please do not submit results from the provided schemes (i.e., Morrison, Thompson or TAU schemes) unless you have modified the codes to enhance/change the aerosol-cloud or precipitation functionality.

## Intercomparison Test-Cases

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
0 - 10 mins - During this period maximum supersaturation is restricted to 0.5%, all velocities are 0 and precipitation processes are switched off.
10 - 20 mins - supersaturation is no longer restricted and velocities gradually increase to the wctrl but precipitation processes are switched off
20 - 30 mins - run for 10 minutes with standard wctrl (0.25 or 1 m s-1) with precipitation switched off
30 - 3 hours 30 mins - run with precipitation switched on
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

Aerosol mass activation rate
Aerosol converted from cloud droplets to rain drops via autoconversion
Aerosol converted from cloud droplets to rain drops via accretion
Aerosol lost to the surface via precipitation removal
Aerosol regeneration rate via droplet evaporation
Finally, for the bin microphysics schemes, we have requested that bin diagnostics be included so that we are able to analyse the particle size distributions. This can lead to very large diagnostic files; thereofre, we kindly ask that zipped files be submitted.

Once diagnostics files are available and participants are ready to submit, participants should contact Adrian Hill (adrian.hill@metoffice.gov.uk) and Zach Lebo ( zlebo@uwyo.edu ).
