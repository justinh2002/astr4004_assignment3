# Assignment 3
### Messier 67 Observation & Radial Metallicity Analysis
## Overview
This project contains two main scripts:

```astro_query.py```: Queries Gaia DR3 for stars within 1 degree of Messier 67, applies quality cuts (2MASS photometry and parallax), and generates color-magnitude diagrams.

```metallicity.py```: Analyses radial metallicity in simulated Milky Way-like galaxy data, fits a linear model, and generates residual plots.

## Requirements
Make sure you have the following Python packages installed:


```pip install numpy matplotlib pandas astropy scipy```

## Usage
# Messier 67 Observation Analysis:

Run the script to query Gaia DR3, crossmatch with 2MASS, and plot the results:

```python astro_query.py```

# Radial Metallicity Analysis:
Ensure the FITS file is placed in the data/ directory and run the analysis:

```python metallicity.py```

## Output
The generated figures are saved in the figures/ directory.
