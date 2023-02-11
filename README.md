# THESEUS
Energy Storage Downselection

This repository is the public accessible repository for the THESEUS software.

## Getting started

THESEUS.exe is the executable file to run the THESEUS software. It provides a user interface to provide the input data, such as power plant type and conditions, ambient temperature, and upload a discretized load profile as a .csv file. The outputs include the optimal selection of energy storage technologies, their design and time-varying operation. 

The backend files constituting the energy storage models are in [backend](backend).

### Overview

1. The backend files for the THESEUS framework are divided into battery and non-battery models. 

2. The operational models for the battery technologies have a similar form, such as the dependence of voltage on the internal battery state. Thus, they are represented in a similar way, and are given in [battery](backend/battery). This folder contains the files for Li-ion, NaS and VRFB batteries with representative inputs. Here, battery_test.gms is the main model file and the rest are files for input specifications.

3. The models for the rest of the non-battery technologies can be found in [model](backend/model). Here, the specific storage models are in: [f1.gms](backend/model/pars/f1.gms) for power functions, [f2.gms](backend/model/pars/f2.gms) for energy functions, and [stor_pars.gms](backend/model/pars/stor_pars.gms) for the storage parameters.
