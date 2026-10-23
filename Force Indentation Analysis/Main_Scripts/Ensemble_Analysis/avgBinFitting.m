% Ben G. 2024

% avgBinFitting:
% avgBinFitting will fit the a Herzitan contact model to the ensemble
% average force indentaion curve.
%
% INPUTS:
% 1- .mat file of the extracted force indentation curves
% 2- .mat file droplet geometry for a given N
%    protein condensate.
%
% OUTPUTS:
% 1- Fitted Herztian contact model parameters for the ensemble average
%    force indentation curves.
% 2- A theoretical force indentation curve contructed using the fitted
%    Hertzian contact model parameters.
%
% COMMENTS:
% The workflow provide in this script can be copied and repeated for
% any set of ensemble average force indentation curves.
%

clear all

% load the required data structures for fitting
loadedData_noRNA = load("C:\Users\Ben_G\OneDrive\Documents\MatLab\Indentation Analysis\M_5_binData.mat"); % INPUT 1
loadedData_noRNA_radius = load("C:\Users\Ben_G\OneDrive\Documents\MatLab\Indentation Analysis\M_5_radiusData.mat"); % INPUT 2
avgRadius_noRNA = mean(loadedData_noRNA_radius.radiusData);

% Extract the indentation array - effetively the bin centers from the
% ensemble average curve
indentation = loadedData_noRNA.binCenters;

[E_data_hertz,dropletYoungModulusFit_noRNA, x, y, eff_R_noRNA] = fitHertzModel(indentation,loadedData_noRNA.avgForces,0.5,avgRadius_noRNA,0.45); % OUTPUT 1
theoreticalFIC_noRNA = (4/3)*(dropletYoungModulusFit_noRNA/10^5)*sqrt(eff_R_noRNA)*((indentation/1000).^(3/2)); % OUTPUT 2

