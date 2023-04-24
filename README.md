[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# Runing Susi Peatland forest productivity model for Swedish sites

# What to expect:
This scripts provides some functionalities to help running the susi model by using a spreadsheet, it also automates the weather data collection using the SMHI public data.  Note that it doesn't modify any functionality inside the susi.

Thake in consideration that this is an exercise dveloped for specific sites in Sweden and it has some gaps of funcionality that can be improved. It is possible to replicate it and be used for different sites, but it should be used carefully.  

# About SUSI:
Here is some extract from the original description:

---
># Peatland simulator SUSI, thoroughly revised version 2022
>## - Mechanistic model for peatland hydrology, biogeochemistry and forest growth
><h3><center>Annamari (Ari) Laurén, University of Eastern Finland </center></h3>
>Original version based on:
>Laurén, A.; Palviainen, M.; Launiainen, S.; Leppä, K.; Stenberg, L.; Urzainki, I.; Nieminen, M.; Laiho, R.; Hökkä, H. Drainage and Stand Growth Response in Peatland Forests—Description, Testing, and Application of Mechanistic Peatland Simulator SUSI. Forests 2021, 12, 293. https://doi.org/10.3390/f12030293
>
>## Short model description
>SUSI describes hydrology, biogeochemical processes and stand growth along a 2D cross-section between two parallel ditches (**Figure 1a**). The structure of SUSI is modular (**Figure 1b**, modules in green boxes). Hydrology modules simulate above-ground and below-ground water fluxes and storages and compute daily water table ($WT$). The peat temperature module simulates temperatures at different depths in the peat ($T_{peat}$). $WT$ and $T_{peat}$ affect the organic matter decomposition rate ($OM_{decroot}$) and the supply of nutrients. Furthermore, $WT$ directly controls net primary production ($NPP$) by scaling it down in case $WT$ is too high through air-filled posrosity function ($f(ϵ_A)$) or if too low with drought function ($f_w$). The nutrient balance module allocates the nutrient supply to stand, litter and ground vegetation, and together with the $NPP$ module they provide constraints for the stand growth module. The new stand volume follows Liebig’s law so that the new stand volume, at the end of the annual time step, is set to the minimum supported by $NPP$ ($V_{NPP}$) or by the supply of N, P or K ($V_{N,P,K}$). The new stand dimensions resulting from the growth are determined by allometric functions and have a feedback to the hydrology, decomposition, ground vegetation and NPP modules. The main computational loop uses a daily time step while the stand growth, ground vegetation and nutrient balance modules are updated annually. SUSI outputs include $WT$, stand growth rates, above-ground biomass and nutrient and water balance components. The model variables and their units are presented in Abbreviations of Laurén et al 2021.

[Read the article](<https://doi.org/10.3390/f12030293>)\
[Download the susi source code from GitHub](https://github.com/annamarilauren/susi_2022)

---

# About this script
Broadly, it can be said that SUSI model requires three types of input data at the site level:
- Site parameters and description
- Daily weather data
- Yearly modeled stand data and biomass 

This script aims to provide SUSI users with a simple tool to fulfill some of the model requirements, adding two basic functionalities. 

The first functionality is provided by a parser that reads one spreadsheet allowing the user to modify the parameters without the need to modify the susi_para.py file, it also offers a brief description of the parameters. The user can modify all the parametrized including the default, such as species, ditch depth, peat type, organic layers, fertilization, nutrient deposition, and others that describe predefined sites, which applies mainly to Finland. 

The second functionality is the automatization of the weather data collection.  This separate module uses the site’s location of interest and the simulation period defined in the parameter file to collect the weather data from the nearest SMHI stations. The script compiles all the data required by the model, sets a continuous timeline, fills the data gaps, and delivers the weather file in the format required by the SUSI model.

>Note that the SUSI requires projected stand data, which is not considered in this script.







---
///  notes from here

To run all sites at once:
    
    swe_susi.py

To run sites individually

    swe_susi.ipynb


