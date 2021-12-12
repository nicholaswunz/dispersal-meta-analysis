# Dataset for dispersal meta-analysis 

## data/ind_disp_raw_data.csv  
### Variable description
*study_ID*: Unique identifier for each paper;  
*author_first*: Surname of the first author of the paper;  
*year*: Year of publication;  
*extracted*: Where the data was extracted;  
*thermal strategy*: Study species identified as an ectotherm or endotherm;  
*taxa*: Generalised grouping of the species based on phylogenetic relatedness;  
*species*: Binomial nomenclature for species;  
*species_OTL:* Scientific name to match the Open tree of life database (https://tree.opentreeoflife.org);  
*origin*: Study species sourced either from the wild or raised in captive settings (past F1) such as farms, aquaculture, laboratories;  
*setting:* Experimental design of the study either in indoor laboratory conditions, outdoor laboratory conditions (e.g., mesocosms, aviaries), or in the wild (no movement restriction);  
*age*: Study species categorised by life stage, either adult (reproductive) or juvenile (non-reproductive);  
*sex*: Sex of the study species (male, female, mixed);  
*mean_mass_g*: Mean body mass of the study species (in g);  
*lon*: Longitude coordinate of the study site in decimal degree (applicable for animals collected from the wild);  
*lat*: Latitude coordinate of the study site in decimal degree (applicable for animals collected from the wild);  
*trt_temp*: Mean temperature of the experiments (in °C);  
*bio_hier*: Biological hierarchy of organisation for the specific response variable (organism, cellular, or macromolecule);  
*trait*: Generalised functional traits of the response variable (see **Table S1** for full description of traits);  
*response*: Specific response measured in the study (see **Table S1** for full description of responses);  
*response_unit*: Unit for the response presented in the study; 
*mass_corrected*: Whether the response was corrected for mass or independent of body mass (y = yes, n = no);  
*dis_trait*: Generalised category of animal movement based on response description in the methods section (activity, exploration, dispersal);  
*dis_unit*: Unit for the animal movement response presented in the study;  
*measure_time_h*: Time duration of the response measurements (in hours);  
*time_lag_days*: Time difference between the response measurements and animal movement measurements (in days);  
*sample_size*: Sample size of the measurements;  
*chi*: Chi-square statistics value;  
*F*: F statistics value;  
*t*: t statistics value;  
*df*: Degrees of freedom value;  
*corr_coeff*: Correlation coefficient value;  
*corr_corrected*: Binary whether the correlation coefficient value was corrected for direction (y = yes, n = no);  
*dispersal_def*: Study authors description of dispersal;  
*notes*: General comments;  
*title*: The title of the paper extracted;  
*link*: Link to the online article.

## data/pop_disp_raw_data.csv
### Variable description
*study_ID*: Unique identifier for each paper;  
*author_first*: Surname of the first author of the paper;  
*year*: Year of publication;  
*extracted*: Where the data was extracted;  
*thermal strategy*: Study species identified as an ectotherm or endotherm;  
*taxa*: Generalised grouping of the species based on phylogenetic relatedness;  
*species*: Binomial nomenclature for species;  
*species_OTL:* Scientific name to match the Open tree of life database (https://tree.opentreeoflife.org);  
*disp_mode*: Mode of dispersal for the study species (aerial, aquatic, or terrestrial);  
*type*: Dispersal type (invasive, where a species is introduced to a new area and dispersing from the invasion centre, or native, where a species is expanding from a known range);   
*core*: Location of the original population (relative to study);  
*front*: Location of the population at the dispersal front;  
*age*: Study species categorised by life stage, either adult (reproductive) or juvenile (non-reproductive);  
*sex*: Sex of the study species (male, female, mixed);  
*mean_mass_g*: Mean body mass of the study species (in g);  
*lon_core*: Longitude coordinate of the original population site in decimal degree;  
*lat_core*: Latitude coordinate of the original population site in decimal degree;  
*lon_front*: Longitude coordinate of the range expanding population site in decimal degree;  
*lat_front*: Latitude coordinate of the range expanding population site in decimal degree;  
*dist_km*: Distance between the core population and dispersal front;  
*time_core*: Year the core population established (relative to study);  
*time_front*: Year the dispersal population established;  
*time_diff*: Time between the core population and dispersal front (in years);  
*dispersal_rate_km_y*: Rate of dispersal from the core population to the dispersal front as dist_km over time_diff (in km/year);  
*rainfall_core_mm_y*: Mean yearly precipitation of the core population site (in mm/year);  
*rainfall_front_mm_y*: Mean yearly precipitation of the dispersal front site (in mm/year);  
*rainfall_diff*: Difference in mean yearly precipitation between the core population and dispersal front (in mm/year);  
*temp_core*: Mean yearly temperature of the core population site (in °C);  
*temp_front*: Mean yearly temperature of the dispersal front site (in °C);  
*temp_diff*: Difference in mean yearly temperature between the core population and dispersal front (in °C);  
*exp_temp*: Study temperature when response was measured (in °C);  
*bio_hier*: Biological hierarchy of organisation for the specific response variable (organism, cellular, or macromolecule);  
*trait*: Generalised functional traits of the response variable (see **Table S2** for full description of traits);  
*response*: Specific response measured in the study (see **Table S2** for full description of responses);  
*response_unit*: Unit for the response presented in the study;  
*mean_core*: Mean response of the core population;  
*sd_core*: Response standard deviation of the core population;  
*se_core*: Response standard error of the core population;  
*n_core*: Sample size of the core population;  
*mean_front*: Mean response of the dispersal front population;  
*sd_front*: Response standard deviation of the dispersal front population;  
*se_front*: Response standard error of the dispersal front population;  
*n_front*: Sample size of the dispersal front population;  
*lnRR*: Natural log of response ratio (mean_front / mean_core);  
*notes*: General comments;  
*title*: Title of the paper extracted;  
*link*: Link to the online article.
