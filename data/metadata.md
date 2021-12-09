# Dataset for dispersal meta-analysis 

## data/ind_disp_raw_data.csv  
### Variables descriptions
*bisphenol:* Bisphenol treatments (BPA, BPF or BPS);  
*flume ID:* Unique flume identifier;  
*ID:* Unique fish identifier;  
*ID_2:* Unique fish identifier (excluding bisphenols);  
*treatment:* Temperature-exposure grouping;  
*temp:* Acclimation temperature (°C);  
*exposure:* Control or exposed to bisphenol (30 µg l⁻¹);  
*temp.acute:* Acute test temperature exposure (°C);  
*body.mass:* Wet weight (g) prior to temperature-bisphenol exposure;  
*length:* Standard lenght (tip of the snout to the posterior end of the last vertebra; cm) prior to temperature-bisphenol exposure;  
*length.2:* Standard lenght (m) prior to temperature-bisphenol exposure;  
*condition:* Fish condition (100 x (mass_g / (length_m * 100) ^ 3) prior to temperature-bisphenol exposure;  
*ucrit_pre:* Critical sustained swimming speed (m s⁻¹) prior to temperature-bisphenol exposure;
*ucrit2_pre:* Critical sustained swimming speed corrected for body length (BL s⁻¹) prior to temperature-bisphenol exposure;
*body.mass_post:* Wet weight (g) 21 days post-temperature-bisphenol exposure;  
*length_post:* Standard lenght (tip of the snout to the posterior end of the last vertebra; cm) 21 days post-temperature-bisphenol exposure;  
*length2_post:* Standard lenght (m) 21 days post-temperature-bisphenol exposure;  
*condition_post:* Fish condition (100 x (mass_g / (length_m * 100) ^ 3) 21 days post-temperature-bisphenol exposure;  
*ucrit_post:* Critical sustained swimming speed (m s⁻¹) 21 days post-temperature-bisphenol exposure;  
*ucrit2_post:* Critical sustained swimming speed corrected for body length (BL s⁻¹) 21 days post-temperature-bisphenol exposure; 
*ucrit_diff:* ucrit2_post - ucrit2_pre; 
*ucrit_change:* (ucrit2_post-ucrit2_pre) / ucrit2_pre x 100.

## data/pop_disp_raw_data.csv
### Variables descriptions 
*bisphenol:* Bisphenol treatments (BPA, BPF or BPS); 
*ID:* Unique fish identifier;  
*ID.2:* Unique fish identifier (excluding bisphenols);  
*treatment:* Temperature-exposure grouping;  
*temp:* Acclimation temperature (°C);  
*exposure:* Control or exposed to bisphenol (30 µg l⁻¹);  
*temp.acute:* Acute test temperature exposure (°C);  
*CS.slope.1:* CS activity replicate 1 (µmol min⁻¹ g⁻¹);  
*CS.slope.2:* CS activity replicate 2 (µmol min⁻¹ g⁻¹);  
*CS.dilution:* Dilution factor for CS;  
*LDH.activity:* Mean LDH activity from CS.slope.1 and CS.slope.2 (µmol min⁻¹ g⁻¹);  
*LDH.slope.1:* LDH activity replicate 1 (µmol min⁻¹ g⁻¹);  
*LDH.slope.2:* LDH activity replicate 2 (µmol min⁻¹ g⁻¹);  
*LDH.dilution:* Dilution factor for LDH;  
*LDH.activity:* Mean LDH activity from LDH.slope.1 and LDH.slope.2 (µmol min⁻¹ g⁻¹);  

