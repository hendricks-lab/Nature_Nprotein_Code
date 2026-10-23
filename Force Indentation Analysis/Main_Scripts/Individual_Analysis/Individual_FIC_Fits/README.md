# N Protein Condensate Force–Indentation Data

This repository contains MATLAB (`.mat`) structures with force–indentation curve (FIC) fitting results for individual N protein condensate droplets under different RNA conditions.

## Experimental conditions

The dates used to identify the datasets correspond to the following RNA conditions:

| Experimental date | RNA condition |
| --- | --- |
| March 5 | No RNA |
| March 6 and April 9 | vgRNA |
| March 12 and April 12 | F3RNA |
| March 18 | FCRNA |

## Contents of the MATLAB structures

Each droplet-level structure contains:

- **Force–indentation fits:** fitting results for all FICs performed on a single droplet.
- **Droplet radius:** the radius of the droplet associated with those measurements.
- **Estimated prestress:** an estimate derived from the optical trap signal offset recorded before the measurement.

The structures retain the association between the individual FIC fits, droplet size, and estimated prestress, allowing measurements to be compared across droplets and RNA conditions. Multiple FICs within a structure are repeated measurements on the same droplet.

## Loading the data

Load a file in MATLAB using:

```matlab
data = load('filename.mat');
```

Replace `filename.mat` with the desired filename. To inspect the variables stored in a file before loading it, use:

```matlab
whos('-file', 'filename.mat')
```

The loaded structures can then be explored in the MATLAB Workspace to access the fitting results, droplet radius, and estimated prestress.
