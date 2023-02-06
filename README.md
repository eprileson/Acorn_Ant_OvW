# Overwintering Project Scripts
Is there evolutionary divergence in Acorn Ant Overwintering physiology?
To find out, we ran a common garden rearing design for field caught urban and rural ants under overwintering conditions and then tested their heat tolerance, cold resistance, metabolic rate and survival

## Scripts
R code for each physiological measure under overwintering conditions

- Worker count / Survival
- Metabolic Rate
- Chill Coma Recovery Time (CCRT)
- Critical Thermal Maximum (CTmax)

## Metadata for each file
description of variables within each dataset

### CCRT
- CCRT csv file in minutes; code changes to seconds for modeling as an integer
- Treatment (source population) and Collection season are factor variables
- Colony ID is randomly generated value

### CTmax
- CTmax csv file: ctmax in degrees C
- Treatment (source population) and Collection season are factor variables
- Colony ID is randomly generated value

### Worker Count / Survival
- Worker start and end are the number of workers at the beginning (0) and end (4) months of the experiment
- Source population is a factor variable
- Colony ID is randomly generated value

### metabolic rate
- mean MR = average metabolic rate in ppm CO2
- mean MR RA = avg metabolic rate in ppm CO2 rolling average (not used in analyses)
- Chamber index = chamber where colony was; chamber 8 = control or blank container
- Q10 = thermal sensitivity quotient of the two metabolic rates at two temperatures raised to the 10th power divided by the temperature difference
- Colletion date = similar to collection season factor from ccrt and ctmax files above
- Source.pop = factor variable
- Colony ID is randomly generated value
