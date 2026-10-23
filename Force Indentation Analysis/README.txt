N PROTEIN CONDENSATE FORCE-INDENTATION ANALYSIS
================================================

This repository contains the main MATLAB scripts, required functions, and a subset of the raw data used to extract, construct, and analyze force-indentation curves from N protein condensates.

The analysis consists of two complementary workflows, followed by a statistical analysis of condensate Young's moduli:

1. Individual force-indentation curve analysis

   Hertzian_Hookian_Model_Automated_Analysis detects individual indentation events, extracts the corresponding force-indentation curves, and fits each curve individually using a Hertzian contact model.

2. Ensemble-average force-indentation analysis

   Extracting_Saving_Indentation_Curves extracts and saves individual force-indentation curves without fitting them. These curves are subsequently analyzed as ensembles using loglog_Ensemble_Average_Analysis and avgBinFitting.

3. Bootstrap comparison of Young's moduli

   BootStrappedDistributionDifferences uses experimentally determined Young's moduli for individual condensates to construct bootstrap distributions of the mean and compare selected protein/RNA conditions.


REPOSITORY ORGANIZATION
=======================

Main Scripts
------------

Individual curve analysis:
1. Hertzian_Hookian_Model_Automated_Analysis

Ensemble-average analysis, in execution order:
2. Extracting_Saving_Indentation_Curves
3. loglog_Ensemble_Average_Analysis
4. avgBinFitting

Statistical analysis:
5. BootStrappedDistributionDifferences


Required Functions
------------------

Force-indentation analysis:
1. fitHertzModel
2. binIndentationData
3. loglogfitting_ensemble
4. Interpolation
5. residuals

Bootstrap analysis:
6. boostrap
7. DifferenceOfDistributions

The function name boostrap is spelled as shown to match the call in BootStrappedDistributionDifferences.


Sample Data
-----------

1. feb_12_analysis

   Contains:
   - A .csv file containing condensate geometries.
   - N .h5 files containing the raw force and trap-position signals exported directly from the LUMICKS C-Trap system.

2. feb_14_analysis

   Contains:
   - A .csv file containing condensate geometries.
   - N .h5 files containing the raw force and trap-position signals exported directly from the LUMICKS C-Trap system.

Here, N denotes the number of .h5 files in each sample-data folder.


SOFTWARE REQUIREMENTS AND SETUP
===============================

The force-indentation analysis code was written and tested using MATLAB 2023/2024b under an Academic license. The MATLAB Signal Processing Toolbox is required.

All required functions must be accessible to MATLAB when running the main scripts. They can either be placed in the same directory as the main scripts or added to the MATLAB path. The bootstrap analysis additionally requires boostrap and DifferenceOfDistributions; any toolbox dependencies of these helper functions must also be satisfied.

File paths to the raw and processed data are system-dependent and can be edited directly within the scripts. Example file paths are provided in the code.

The required scripts, functions, and sample data should take substantially less than 2 minutes to download from the repository. Once downloaded, add them to a MATLAB working directory or the MATLAB path and run the relevant workflow below.

For BootStrappedDistributionDifferences, experimental Young's moduli are entered directly as vectors within the script rather than loaded from external files.


WORKFLOW 1: INDIVIDUAL FORCE-INDENTATION CURVE ANALYSIS
=======================================================

Hertzian_Hookian_Model_Automated_Analysis
-----------------------------------------

This script performs automated analysis of individual force-indentation curves. It detects indentation events from the raw C-Trap data, extracts the corresponding force-indentation curves, and fits a Hertzian contact model to each curve individually.


DETECTION OF INDENTATION EVENTS

The time-series force and trap-position data are divided into subsections of length N. Indentation events are initiated by displacement of the optical trap, resulting in an increase in the recorded trap position.

Each trap-position subsection is evaluated for a significant increase in position. If a subsection contains a significant trap displacement, the time point corresponding to the beginning of that subsection is identified as the start of an indentation event. The same time point is then used to identify the corresponding force data.


FORCE-INDENTATION CURVE EXTRACTION

Polystyrene probes were optically trapped and calibrated near the condensates. Instrument calibration assumed passive thermal motion near a boundary (the coverslip).

Using the trap calibration, measured force, and trap position, condensate indentation is calculated as:

Xind = Xtrap - deltaXforce

where Xind is the condensate indentation, Xtrap is the trap displacement, and deltaXforce is the probe displacement determined from the measured force and trap calibration.

This calculation is performed for every detected indentation event.


HERTZIAN CONTACT MODEL FITTING

A Hertzian contact model is fitted to every force-indentation curve identified from an indentation event.

Sequential measurements are performed on individual condensates. Consequently, each condensate can contain multiple indentation events and therefore multiple independently fitted force-indentation curves.


INPUTS

1. .h5 files containing C-Trap force data and calibration information.
2. .csv files containing N protein condensate geometries (radius).


OUTPUTS

1. An experimental data structure containing the fitted parameters for all individual indentation events, grouped by condensate:

   Experiment
       -> Condensate ID
           -> Individual indentation event
               -> Fitted parameters

The extracted fit parameters can subsequently be used for statistical analysis. Statistical tests and analyses applied to these parameters are described in the main manuscript.

The bootstrap script described below expects one Young's modulus value per condensate. It does not automatically import this structure or combine multiple indentation-event fits into a single condensate-level value.


WORKFLOW 2: ENSEMBLE-AVERAGE FORCE-INDENTATION ANALYSIS
=======================================================

The ensemble-average analysis consists of three stages:

1. Extracting_Saving_Indentation_Curves
2. loglog_Ensemble_Average_Analysis
3. avgBinFitting

The first script extracts and saves the individual curves required to construct ensemble averages. The subsequent scripts operate on ensemble-average force-indentation curves rather than fitting individual indentation events.


STEP 1: Extracting_Saving_Indentation_Curves
--------------------------------------------

This script detects individual indentation events and extracts the corresponding force-indentation curves for subsequent ensemble-average analysis. It generates the collection of curves required for downstream ensemble analysis without fitting each event individually.


DETECTION OF INDENTATION EVENTS

The time-series force and trap-position data are divided into subsections of length N. Indentation events are initiated by displacement of the optical trap, resulting in an increase in the recorded trap position.

Each trap-position subsection is evaluated for a significant increase in position. If a significant trap displacement is detected, the time point corresponding to the beginning of that subsection is identified as the start of an indentation event. The same time point is then applied to the force data.


FORCE-INDENTATION CURVE EXTRACTION

Polystyrene probes were optically trapped and calibrated near the condensates. Instrument calibration assumed passive thermal motion near a boundary (the coverslip).

Condensate indentation is calculated from the trap calibration, measured force, and trap position according to:

Xind = Xtrap - deltaXforce

This analysis is performed for every detected indentation event.


INPUTS

1. .h5 files containing C-Trap force data and calibration information.
2. .csv files containing N protein condensate geometries (radius).


OUTPUTS

1. A 1 x 1 MATLAB data structure containing N condensate-specific 1 x M data structures, where:

   - N is the number of condensates measured for a given condition.
   - M is the number of indentation events detected for a given condensate and may differ between condensates.

Each indentation event contains 1 x S force and trap-position vectors, where S is the number of data points associated with that event.

The resulting data structure is used as an input for loglog_Ensemble_Average_Analysis.


STEP 2: loglog_Ensemble_Average_Analysis
----------------------------------------

This script uses the force-indentation curves generated by Extracting_Saving_Indentation_Curves to construct and analyze ensemble-average force-indentation curves.

The analysis determines the residual between the experimental ensemble-average force-indentation response and the initial linear region of the force-indentation relationship in log-log space.


INPUTS

1. .mat file containing the extracted force-indentation curves.
2. .mat file containing fitted Young's moduli and condensate geometries for a given N protein condensate condition.


OUTPUTS

1. A data structure containing the ensemble-average force-indentation curves for a given protein/RNA condition. Save this data structure for subsequent analysis using avgBinFitting.

2. The residual between the ensemble-average force-indentation curve and the Hertzian contact model in log-log space as a function of indentation.


STEP 3: avgBinFitting
---------------------

This script fits a Hertzian contact model to the ensemble-average force-indentation curve generated during the ensemble analysis.


INPUTS

1. .mat file containing the extracted force-indentation curves.
2. .mat file containing condensate geometries for a given N protein condensate condition.
3. Ensemble-average force-indentation data generated during the preceding analysis stage.


OUTPUTS

1. Fitted Hertzian contact model parameters, including Young's modulus, for the ensemble-average force-indentation curve.

2. A theoretical force-indentation curve constructed using the fitted Hertzian contact model parameters.


COMMENTS

The workflow provided in this script can be copied and repeated for any set of ensemble-average force-indentation curves.


STATISTICAL ANALYSIS: BOOTSTRAP COMPARISON OF YOUNG'S MODULI
=========================================================

BootStrappedDistributionDifferences
----------------------------------

This script performs bootstrap analysis of experimentally determined Young's moduli for individual condensates under four conditions: noRNA, FCRNA, F3RNA, and vgRNA. It constructs bootstrap distributions of the mean Young's modulus and compares selected condition pairs.


INPUTS

1. Four vectors containing experimental Young's moduli for individual condensates:

   - EXP_noRNA
   - EXP_FCRNA
   - EXP_F3RNA
   - EXP_vgRNA

Each entry represents one condensate. The vectors are defined directly within the script and can be edited to analyze another dataset. All values must use the same units across conditions.

The script does not automatically load fitted-parameter structures or calculate condensate-level values from individual indentation events.


BOOTSTRAP ANALYSIS

The experimental Young's moduli are first rounded to the nearest multiple of 10 using round(values, -1).

The number of bootstrap iterations is set by nBootstrap, which defaults to 10,000. The script calls boostrap separately for each condition, returning a bootstrap mean estimate, a confidence interval, and a distribution of bootstrap means.

The sampling procedure and confidence-interval definition are determined by the boostrap helper function.


COMPARISON BETWEEN CONDITIONS

The script passes the bootstrap distributions to DifferenceOfDistributions for five comparisons:

1. FCRNA versus noRNA.
2. F3RNA versus noRNA.
3. vgRNA versus noRNA.
4. F3RNA versus FCRNA.
5. vgRNA versus FCRNA.

For each comparison, the script counts the number of values in the returned difference distribution that are less than or equal to zero. These results are stored in the P_count variables.

The P_count variables contain raw counts, not normalized p-values. The script does not convert these counts into p-values or apply a multiple-comparison correction. The definition and direction of each difference are determined by DifferenceOfDistributions.


OUTPUTS

1. Bootstrap mean estimates, confidence intervals, and distributions of bootstrap means for each condition.

2. Difference distributions for the five condition pairs.

3. Counts of non-positive differences, stored in the P_count variables.

4. A figure showing the four bootstrap distributions of the mean Young's modulus in vertically arranged panels with shared x-axis limits and common y-axis limits. Histogram bin counts are divided by nBootstrap to display relative frequencies.

Results remain in the MATLAB workspace. The script does not automatically save data or figures. Additional plotting code for the difference distributions is included but commented out.


REQUIRED FUNCTIONS

1. boostrap
2. DifferenceOfDistributions

Both functions must be accessible on the MATLAB path. Their implementations determine the bootstrap sampling method, confidence-interval calculation, and difference-distribution calculation.


RUNTIME
=======

The force-indentation analysis scripts run in under 1 minute for datasets comparable in size to the provided sample data.

The bootstrap script uses 10,000 iterations per condition by default. Its runtime depends on the dataset sizes, helper-function implementations, and computer hardware.


REPRODUCING THE MANUSCRIPT ANALYSIS
===================================

Applying the force-indentation workflows to the complete C-Trap dataset reproduces the quantitative force-indentation results presented in Figure 2 of the manuscript.

BootStrappedDistributionDifferences provides the associated bootstrap comparison of condensate Young's moduli across the conditions defined in the script. To repeat this analysis, use the corresponding condensate-level Young's modulus values and ensure that both bootstrap helper functions are available.