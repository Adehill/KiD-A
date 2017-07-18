- 
{:toc}

# Introduction

A major uncertainty in numerical weather prediction (NWP) and climate prediction is the response of clouds and precipitation to changes in aerosol concentrations. Previous intercomparison work has shown that the simulation of precipitation for a given cloud drop number concentration (Nd) and the response of precipitation to changes in Nd in cloud resolving models (CRMs) and idealised kinematic frameworks is very sensitive to both the physical representation of microphysics and numerical complexity of the scheme cloud microphysics (e.g. ). The purpose of this microphysics intercomparison project is to compare detailed size resolved and bulk parametrised microphysics schemes to understand how they simulate aerosol-cloud-precipitation interactions. The first phase of this project employs the Kinematic Driver model (KiD) to compare microphysics schemes when considering the simulations of and the response of warm rain to changes in N_d. Further this intercomparison will investigate the role of in-cloud processing of aerosol  

# Overview of the KiD-A project
The overarching aim of KiD-A is to use the Kinematic Driver model (KiD) to compare detailed and bulk microphysics schemes in a dynamically consistent framework without dynamic feedbacks, which have complicated the interpretation of previous intercomparison projects. The main aims are:
1. Undertake the first kinematic intercomparison of detailed microphysics schemes, i.e. size resolved bin microphysics schemes, superdroplet schemes and 2D aerosol-cloud schemes 
- This will be a sanity check to make sure that the schemes, which are consistently used to develop simpler bulk schemes produce similar results, when forced with the same dynamics. 
- Tests will be performed for a range of initial aerosol concentrations and 1D and 2D kinematic cases 
- Tests will exclude in-cloud aersol processing and focus solely on the precipitation processes, timing and amount 
2. Examine and compare in-cloud aerosol processing from both detailed microphysics and bulk microphysical representations. 
- This stage of the project will compare the detailed and bulk microphysics schemes that can include in-cloud aersol processing. 
- This will be based on only the 2D kinematic cases. 
- The results from this stage will be used to benchmark present modeling capability when considering aerosol-cloud-precipitation interactions. 

# Kinematic Driver Model (KiD)

The KiD model is the computational and numerical basis the KiD-A project. The source code is available from the [main KiD-A github repository]((https://github.com/Adehill/KiD-A) and can be obtained by either of the following methods
- downloading the zip file from [main KiD-A github repository](https://github.com/Adehill/KiD-A)
- clone the repository to your system using `
```
git clone https://github.com/Adehill/KiD-A.git
```
The ```docs``` directory in the download or cloned repository contains the [KiD documentation](https://github.com/Adehill/KiD-A/blob/master/docs/KiD_2.3.2625.pdf),  which details how to build and run the KiD model as well as details about haow to add a microphysics scheme and include diagnostics. This version of the code includes old versions of the Morrison microphysics and Thompson microphysics schemes (both single and double moment rain schemes) and a dummy subroutine that contains details about how to obtain and include the Tel-Aviv University (TAU) bin microphysics schemes.

It is expected that participants will add their microphysics scheme to the KiD model. This will involve writing an interface that couples the generic KiD prognostics to the microphysics variables. The documentation describes how to add a scheme to the KiD framework, while the Morrison, Thompson and TAU provide schemes are provided as examples of adding a microphysics code to the KiD model. In general, there is no need to change the microphysics code so that it will work with the KiD model. Further, diagnostic calls may be required in the microphysics code so that process rates and precipitation rates are output. Otherwise, all required diagnostics will be automatically output once a scheme is successfully coupled. The required diagnostics are discussed below.

Please do not submit results from the provided schemes (i.e., Morrison, Thompson or TAU schemes) unless you have modified the codes to enhance/change the aerosol-cloud or precipitation functionality.

# KiD-A intecomparison testcases

The KiD-A intercomparison uses following testcases, which can all be found in [KiD testcase directory](https://github.com/Adehill/KiD-A/blob/master/namelists)

## 1D and 2D kinematic cases

### Aerosol specifications for 1D and 2D case
In both cases aerosol are assumed to be soluble ammonium sulfate particles. The initial aerosol distribution is assumed to be a single mode lognormal distribution with a lognormal geometric mean diameter = 0.08 * 10-6m and a log standard deviation = 1.4. The initial aerosol number concentrations (Na) are defined in the test case descriptions below. These parameters are defined in the test case namelists. If a participants model does not include aerosol, please set the initial cloud drop number concentration (Nd) to the initial Na defined below.


### 1D case
- This is the initial cloud microphysics case with simple 1D updraft, which is based on Shipway and Hill (2012). The case employs a fixed aerosol or Nd depending on scheme (no aerosol processing)
- Two vertical velocity set-ups are requested for this case, i.e.
   - W1p25 where ```wctrl = 1.25 m s^{-1}```
   - W2 where ```wctrl = 2 m s^{-1}```
      - NOTE: The updraft velocities are lower than previous iterations of this project and Shipway and Hill (2012) because the divergence term has been switched. This is required for a fair comparison between the bin and lagrangian model. 
- For both W1p25 and W2, we require a simulation for Na or Nd (depending on the scheme) = 50, 150 and 300 cm-3.
- For each updraft and Nd, you are asked to simulate the following:
   - cond-evap case - only aersol-activation (dependent on scheme), condensation and evaporation are switched on, i.e. sedimentation, collision-coalescence, breakup, etc. are switched off.
   - Precipitating case - all cloud microphysics processes are switched on; however, no removal or replenishment of aerosol and no in-cloud aerosol processing (in schemes that have this functionality) should be simulated.
- Templates for a bin and bulk namelist required to run the 1D case, i.e. namelists/kida_icmw1D_bin_template and namelists/kida_icmw1D_bulk_template.nml are provided in the repository (see [1D bulk namelist template](https://github.com/Adehill/KiD-A/blob/master/namelists/kida_icmw1D_bulk_template.nml) and a [1D bin namelist template](https://github.com/Adehill/KiD-A/blob/master/namelists/kida_icmw1D_bin_template.nml) ).
- To build and run this case with gfortran type the following (assuming you are in the KiD-A directory) 
```
make COMPILER=gfortran CASE=1D all
./bin/KiD_1D.exe namelists/kida_icmw1D_bulk_template.nml
```
   
### 2D stratocumulus (Sc 2D)
- based on the case 4 from the 9th International Cloud modelling workshop (2016).
- we require a simulation to be run with maximum vertical velocity (wctrl) = 0.25 and 1 m s-1. Once again wctrl is set in the namelist for the case. Na (or Nd depending on scheme) = 50, 150 and 300 cm-3 for the following configurations:
   - cond-evap case only aersol-activation (dependent on scheme), condensation and evaporation are switched on, i.e. sedimentation, collision-coalescence, breakup, etc. are switched off.
   - Precipitating case - no removal or replenishment of aerosol and no in-cloud aerosol processing (where this is possible with your scheme)
   - Precip case with aerosol processing - this will include the removal of aerosol by activation scavenging, replenishment of aerosol by evaporation and in-cloud aerosol processing.
- In all cases, we apply the following start-up to prevent numerical problems the initial supersaturation 
   - 0 - 10 mins - During this period maximum supersaturation is restricted to 0.5%, all velocities are 0.
   - 10 - 20 mins - supersaturation is no longer restricted and velocities gradually increase to the wctrl
   - 20 - 3 hours - run simulation
   - To set-up this initialisation, we have added the following switches to the KiD model and the namelists:
```
smax_limit_time = 600.0 ! The time from the beginning of the simulation with in which a maximum supersaturation is applied
smax = 0.5 ! the maximum percentage supersaturation during smax_limit_time
```
   - Participants need to add these switches to the appropriate parts of the microphysics code so that the recommended initialisation is used. The switches are stored in the namelists module.
- Templates for a bin and bulk namelist required to run the 2D Sc case, i.e. namelists/kida_icmwSC_2d_bin_template.nml and namelists/kida_icmwSC_2d_bulk_template.nml, are provided in the repository (see [Sc 2D bin template](https://github.com/Adehill/KiD-A/blob/master/namelists/kida_icmwSC_2d_bin_template.nml) or [Sc 2d bulk template](https://github.com/Adehill/KiD-A/blob/master/namelists/kida_icmwSC_2d_bulk_template.nml) ).
- To build the 2D case with gfortran type the following (assuming you are in the KiD-A directory) 
```
make COMPILER=gfortran CASE=ICMW_SC all
```

## Box model tests with KiD
In order to understand potential differences between the detailed schemes, we propose running some simple box cases which will test the evolution of the drop size distrbution resulting from 1) collision-coalescence and 2) condensational growth. The aim of both of these tests is to prescribe the initial cloud drop size distribution using bulk variables of total mass, number and shape of the distribution, so that all models start from the same point. Once initialised the cases either test collision-coalescence or condensational growth, to understand how different schemes evolve given the same starting point. To implement these test each scheme will need to be able to prescribe gamma drop size distribution and discretize this onto their respective bins or into lagrangian parcels. The routine ```src/init_bins.F90``` has been developed to do this for the 2-moment TAU scheme, and this routine has been adapted and applied to the 2d-bin scheme and a lagrangian model.

To build the box model case type
```
make COMPILER=gfortran CASE=ICMW_BOX all
```

### Box - Condensational growth
- A template for this case is provided [namelists/kida_icmwBOX_CE_bin_template.nml](https://github.com/Adehill/KiD-A/blob/master/namelists/kida_icmwBOX_CE_bin_template.nml). In this namelists, the initial bulk distribution is specified as follows
```
pctrl_v(1) = 1. ! hydrometeor to initialise
pctrl_v(2) = 0.0001 ! bulk mass kg/kg (integrate over the distribution)
pctrl_v(3) = 300.e6 ! bulk Nc /kg
pctrl_v(4) = 0.0    ! mu, shape parameter
```
- Then the following switches force only condensation to occur
```
! TAU switches
l_act=.false.       ! Allow activation
l_cond_evap=.true. ! Allow condensation and evaporatio (bin model)
l_coll_coal=.false. ! Allow collision-coalescence (bin model)
l_break=.False.     ! Allow collisional breakup (bin model)
l_fix_supersat = .True. ! Allow user to prescribe supersaturation. This is
                         ! only applicable to box condensational growth case. 
                         ! Be careful setting this for other cases!
```
- Condensational growth is then forced by ```l_fix_supersat = .True.``` and a prescribed supersaturation of 0.1% for 10 minutes, using ```smax = 0.1```
- In order to test condensational growth alone, it is important that sedimentation, collision-coalescence and activation are switched off in your scheme (as demonstrated in the above code block).
- As with previous tests, we would like participants to provide a tests with
   - Nd = 50, 150 and 300 cm^{-1}
   - mu = 0.0 and 2.5

### Box - Collision-coalescence growth
- A template for this case is provided [namelists/kida_icmwBOX_COLL_bin_template.nml](https://github.com/Adehill/KiD-A/blob/master/namelists/kida_icmwBOX_COLL_bin_template.nml). In this namelists, the initial bulk distribution is specified as follows
```
pctrl_v(1) = 1. ! hydrometeor to initialise
pctrl_v(2) = 0.001 ! bulk mass kg/kg (integrate over the distribution)
pctrl_v(3) = 300.e6 ! bulk Nc /kg
pctrl_v(4) = 0.0    ! mu, shape parameter
```
- The above code block is same as the condensation case except that the bulk mass is increased by an order of magnitude, i.e. ```pctrl_v(2) = 0.001 ! bulk mass kg/kg ```
- In order to test collision-coalescence growth alone, it is important that activation, condensation, sedimentation are switched off, e.g.
```
! TAU switches
l_act=.false.       ! Allow activation
l_cond_evap=.false. ! Allow condensation and evaporatio (bin model)
l_coll_coal=.true. ! Allow collision-coalescence (bin model)
l_break=.False.     ! Allow collisional breakup (bin model)
l_fix_supersat = .false. ! Allow user to prescribe supersaturation. This is
                         ! only applicable to box condensational growth case. 
                         ! Be careful setting this for other cases!
```

- As with previous tests, we would like participants to provide a tests with
   - Nd = 50, 150 and 300 cm^{-1}
   - mu = 0.0 and 2.5

# Diagnostics

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


