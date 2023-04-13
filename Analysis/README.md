## Code used to run and analyse soybean FBA model
1. model_with_databases_subs.mat = Soybean FBA model used in Holland et al., (2023) 
2. combined_model_with_databases_subs.mat = Soybean FBA combined model (with ureide & amide export) used in Holland et al., (2023)
3. NvsMass.m = matlab file to run analysis for varying soil ammonium uptake against plant RGR
4. HighRGR_amidevureideFVA.m = analysis to run flux variability analysis with combined and ureide soybean FBA model
5. PlimitedAmidevUreideFVA.m similar analysis to HighRGR_amidevureideFVA.m except when simulating a phosphorus limitation
6. NvsNfix.m = analysis for varying soil ammonium uptake against nitrogen fixation rate 
7. VaryingSR_massbalance.m supplementary analysis to test the effect of varied shoot:root ratio on carbon cost of N fixation and RGR
8. photo_test_mass_balanced.m supplementary analysis to test the effect of CO2 uptake rate against carbon cost of N fixation 

