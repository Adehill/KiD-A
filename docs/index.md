- 
{:toc}

# Introduction

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



