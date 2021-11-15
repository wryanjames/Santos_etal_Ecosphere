# Santos_etal_Ecosphere
[![DOI](https://zenodo.org/badge/428325942.svg)](https://zenodo.org/badge/latestdoi/428325942)

This repository contains the data and code that was used in Santos et al. (2021). 

Santos RO, James WR, Nelson JA, Rehage JS, Serafy J, Pittman SJ, & D Lirman (2021) Influence of seascape spatial pattern on the trophic niche of an omnivorous fish. Ecosphere

BBpinfish.csv - File containing raw data of pinfish
  - Site - site where pinfish collected
  - Length(mm) - Total length of pinfish in mm
  - Weight(g) - Weight of pinfish in g
  - Seascape - Type of seascape
  - Zone - Salinity zone (zone 1 = High/Stable, zone 2 = Low/Variable)
  - ID - unique ID of pinfish
  - d13C - d13C stable isotope value
  - d15N - d15N stable isotope value
  - Benthic algae - benthic algae mixing model estimated source contribution
  - Drift - drift algae mixing model estimated source contribution
  - Epiphytes - epiphytic algae mixing model estimated source contribution
  - Seagrass - seagrass mixing model estimated source contribution
  - TL - calculated trophic level based on mixing model results 

sz1.csv - source isotope values for zone 1 (high and stable salinity)
  - Species - ID of source
  - Zone - Salinity zone (zone 1 = High/Stable, zone 2 = Low/Variable) 
  - Meand13C - mean d13C source stable isotope value
  - Meand15N - mean d15N source stable isotope value
  - Concd13C - mean percentage of carbon concentation in source
  - Concd15N - mean percentage of nitrogen concentation in source
  - SDd13C - standard deviation of d13C values
  - SDd15N - standard deviation of d15N values
  
sz2.csv - source isotope values for zone 2 (low and variable salinity)
  - Species - ID of source
  - Zone - Salinity zone (zone 1 = High/Stable, zone 2 = Low/Variable) 
  - Meand13C - mean d13C source stable isotope value
  - Meand15N - mean d15N source stable isotope value
  - Concd13C - mean percentage of carbon concentation in source
  - Concd15N - mean percentage of nitrogen concentation in source
  - SDd13C - standard deviation of d13C values
  - SDd15N - standard deviation of d15N values
  
  
pinBB_hypervolumes.R - script used for hypervolume analysis of pinfish