# This is the repository of the article *General Laws of Biodiversity: Climatic Niches Predict Plant Range Size and Ecological Dominance Globally*, accepted for publication in PNAS

- Citation: Moulatlet, G. M., Merow, C., Maitner, B., Boyle, B., Feng, X., Frazier, A. E., Hinojo-Hinojo, C., Newman, E. A., Roehrdanz, P. R., Song, L., Villalobos, F., Marquet, P. A., Svenning, J.-C., & Enquist, B. J. (2025). General laws of biodiversity: Climatic niches predict plant range size and ecological dominance globally. Proceedings of the National Academy of Sciences.
﻿
- Contact: Gabriel M. Moulatlet (mandaprogabriel@gmail.com)
﻿
- Abstract: A longstanding question in ecology asks whether or not species that achieve large geographic ranges also have large climatic niche breadths. Using a dataset of ~250,000 terrestrial plant species spanning diverse clades (bryophytes, ferns, gymnosperms, and flowering plants), we demonstrate a consistent positive relationship between geographic range size and climatic niche breadth across latitudinal and elevational gradients. This relationship holds across major phylogenetic groups, suggesting a general biogeographical rule for range size variation. Our findings indicate that latitudinal and elevational gradients in range size arise from selective pressures and species sorting based on climatic tolerance. Additionally, we show that species with larger range sizes tend to be ecologically dominant, supporting a long-suspected connection between range size, niche breadth, and local and regional abundance. Our results suggest a spectrum of dominance, where species with extensive geographic ranges and broader climatic tolerances tend to be more abundant. We posit that the relationship between range size, niche breadth, and ecological dominance is an emergent macroecological pattern that can be used for understanding and predicting the impacts of climate change on species distributions.
﻿
- Responsible for writing code: Gabriel M. Moulatlet (mandaprogabriel@gmail.com)
﻿
- Folders and files: There are two folders, data (which contain the data necessary to run the R codes) and code (which contain the R scripts).
﻿
* Data folder:
  * "data.csv": contains the columns
    
			- "species" which has species names
			- "convex" which has the niche breadth calculations using the convex hull method
			- "convex2" sqrt-root transformed convex
			- "midX" and "midX" which have the centroid coordinates of species geographic ranges, in Mollweide projection
			- "area" which has the geographic range area in km2
			- "area10" log10-transformed area
			- "region" defines the latitutinal zones, if the species range is located in tropical, temperate or in both tropical/temperate zones
			- "band" which defines the longitudinal bands, if the species ranges is mostly in the America, Europe/Africa or Asia/Oceania
			- "elevation" is the elevation at the centroid of the species ranges
			- "higher_plant_group" is the taxonomic group, ferns and lycophytes, flowering plants, gymnosperms or bryophytes
			- "family" is the botanical family
			- "midX84" and "midX94" which have the centroid coordinates of species geographic ranges, in WGS84 projection
			- "ele"  defines if the species is upland or lowland species, as set by the 1000 m.a.s.l. threshold
			- "cLat"  defines the latitudinal bin every 5 degrees
			- "source" only for tropical trees, indicates the dominance status, according to Cooper et al.2024 Nature Ecology and Evolution

	* "dataNULL.csv": contains the columns
     
			- "species which has the species names
			- "cnv" which is the convex hull simulation number
			- "cv" is the convex hull calculation value
			- "type" depicts if the cv value is the observed or the simulated one.
﻿
* Code folder:
    * Analysis_figures_stats_v2025_10_09_toRepo.R: contains the codes to generate the figures and the results of the manuscript

 
* Software version:
   - R version 4.5.1 (2025-06-13 ucrt)
   - Platform: x86_64-w64-mingw32/x64
   - Running under: Windows 11 x64 (build 26200)
   - Rstudio 2025.9.1.401 (Cucumberleaf Sunflower)

* Funding sources: NSF (grants 2225078 and 2225076). CONAHCyT Mexico, the U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research, Earth and Environmental Systems Sciences Division and Data Management Program, under Award Number DEAC02- 05CH11231, as part of the Watershed Function Scientific Focus Area and the ExaShed project.

* Please also check biendata.org
