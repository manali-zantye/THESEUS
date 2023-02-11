# THESEUS
Energy Storage Downselection

This repository is the public accessible repository for the THESEUS software.

## Getting started

THESEUS.exe is the executable file to run the THESEUS software. It provides a user interface to provide the input data, such as power plant type and conditions, ambient temperature, and upload a discretized load profile as a .csv file. The outputs include the optimal selection of energy storage technologies, their design and time-varying operation. 

The backend files constituting the energy storage models are in [backend](backend).

### Overview

The backend files for the THESEUS framework are divided into battery and non-battery models. The operational models for the battery technologies have a similar form, such as the dependence of voltage on the internal battery state. Thus, they are represented in a similar way, and are given in [battery](battery).
