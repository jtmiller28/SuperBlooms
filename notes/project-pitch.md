### Short Project Proposal 
#### Research Motivation
Superblooms are widely recognized as major tourist attractions, drawing significant
public interest. However, their ecological importance extends far beyond their visual appeal.
These rare events play a vital role in desert ecosystems by supporting plant reproduction,
sustaining wildlife populations, and preventing soil erosion. Additionally, some groups have
strong cultural ties to superblooms (Winkler & Brooks 2020). To monitor and preserve these
events we propose the utilization of community science data from iNaturalist to identify what
species appear during superblooms and how the blooms are changing over time. While some of
this data may be skewed by observer bias we intend to implement strategies to mitigate this
inaccuracy. By identifying the core taxa in superblooms we will obtain greater insights into how
environmental factors may be impacting their occurrence, how we can protect these events
throughout global warming, and how these blooms may affect local ecosystems. 
#### Target Audience 
Our findings would be most useful for scientists and resource managers working to better
understand superbloom events and how best to protect them to ensure the health of desert
ecosystems at large. The audience would also include citizen scientists and outdoor enthusiasts
that may have contributed to the data or just have a vested interest in the conservation of
superblooms. The phenomenon of the superbloom is highly charismatic and holds cultural
significance for desert communities. It draws in many tourists and helps shape national
perspectives regarding desert ecosystems (Winkler & Brooks 2020). For this reason, the
audience of our research extends outside of the scientific community, as our work may be helpful
in educating the public about what species they are observing, how they can sustainably interact
with superblooms, and how to avoid biases in data collection when contributing to citizen
science.
#### Logistics & Timeline
_Our Definition of a Superbloom_: A short temporal window where an abnormally high density of flowering occurs within a spatial area in arid regions. <br> 
Hypothesis: Superblooms in the North American SW will contain consistent core taxa across
independent events, and community science data from iNaturalist can be used to identify these
taxa despite observer bias. <br> 
All methods for this project will be collaborated on collectively, and code openly shared via github: https://github.com/jtmiller28/SuperBlooms . For all analyses we will be utilizing the programming language R (v.4.3). There will be 4 major phases of data wrangling, analysis, and interpretation: Data acquisition, Spatial & Temporal Data Subsetting, Species Extractions, Assessments of Biases, and Data Visualizations. 
Data collection - We will use phenovision’s latest data download to acquire computer vision phenology annotated records for plant iNaturalist records worldwide (Dinnage et al. 2025). This download was acquired and shared by Rob Guralnick 01-2025. 
Spatial & Temporal Subsetting -  using the programming language R, we will load the csv file using the package data.table (v.1.16.4) and subset the data to the Mojave Desert & Sonoran Desert level 3 ecoregion shapefiles (Omernik & Griffith 2014) with the package sf (v.1.0-19). We will then overlay 5x5km hex-cells onto the ecoregions’ shapefiles, and assign occurrence’s spatial IDs according to these hex-cells using sf & dplyr (v.1.1.4). Within each hexcell, we will bin the data by year. These spatial & hex style bins may be adjusted if warranted during downstream analysis steps. <br> 
Species Extractions - Using dplyr we will group data by cell-time bin and sort by high density observation cells. Per cell, we will summarize the number of observations per species for both all time and for each time bin. To determine what qualifies as superbloom events, we will use a statistical test that identifies what time bins have significantly higher flowering observations as compared to the mean of all time bins in the dataset. This statistical test will be checked for assumptions,  and further evaluated after the bias assessment phase. <br>
Assessments of Biases: Temporal, spatial, and taxonomic biases will be assessed using the above tables of data noted above in the species extractions phase. Methods for mitigating/correcting these biases are still in development, but as a goal we want to assess whether the taxonomic representation of a cell is weighted by observation effort. <br>
Data Visualizations - Using ggplot2 (v.3.5.1), we will construct visuals including a time series of the number of observations per species per cell, a 5x5 km hex gridded map of the Mojave & Sonoran regions with observation density, and a figure (TBD) illustrating the presumed overlap between cells. <br>
Extra: If we have extra time we will consider using environmental data sourced from PRISM to model if the core superbloom taxa triggers flowering in response to similar environmental variables. <br>
Timeline: <br>
Week 0 (03/03/2025 - 03/07/2025) - Data Collection - Finished, data provided by Phenovision <br>
Week 1 (03/10/2025 - 03/14/2025) - Project Coordination - In Progress (JT, Kenna, & Lou) <br>
Week 2 (03/17/2025 - 03/21/2025) - Spring Break - Spatial & Temporal Data Subsetting - JT  <br>
Week 3 (03/24/2025 - 03/28/2025) - Species Extractions - Kenna & Lou (JT will provide help) <br>
Week 4 (03/31/2025 - 04/04/2025) - Assessments of Biases - Kenna, Lou, & JT <br>
Week 5 (04/07/2025 - 04/11/2025) - Data Visualizations - Kenna, Lou, JT <br>
Week 6 (04/13/2025 - 04/17/2025) - Preliminary Project Findings for Class - Kenna, Lou, JT <br>
Week 7 (04/20/2025 - 04/23/2025) - Wrap-up & Present - Kenna, Lou, JT <br>
Bibliography <br>
1. Dinnage, R., Grady, E., Neal, N., Deck, J., Denny, E., Walls, R., Seltzer, C., Guralnick, R., & Li, D. (2024). PhenoVision: A framework for automating and delivering research-ready plant phenology data from field images. bioRxiv. https://doi.org/10.1101/2024.10.10.617505 <br>
2. Omernik, J. M., & Griffith, G. E. (2014). Ecoregions of the conterminous United States: Evolution of a hierarchical spatial framework. Environmental Management, 54(6), 1249–1266. https://doi.org/10.1007/s00267-014-0364-1 <br>
3. Winkler, E, D., Brooks, E. (2020). Tracing Extremes across Iconic Desert Landscapes: Socio-Ecological and Cultural Responses to Climate Change, Water Scarcity, and Wildflower Superblooms. Human Ecology, 48, 211-223. https://doi.org/10.1007/s10745-020-00145-5 <br>
