load('combined_with_formulas.mat');
load('USDA110_model.mat');
brady = USDA110;
fixing_mass_balance;
save('Combined_mass_charge_balanced_model.mat', 'model');
