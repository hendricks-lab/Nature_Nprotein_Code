% Ben G. 2024

% INDENTATION ANALYSIS:
% Using the MATLAB outputs from "Extracting_Saving_Indentation_Curves" - to
% determine the residual between the experimental ensemble average and the
% initial linear region of the loglog force indenation curve. 
% 
% INPUTS:
% 1- .mat file of the extracted force indentation curves
% 2- .mat file of fitted Young's Moduli and droplet geometry for a given N
%    protein condensate.
%
% OUTPUTS:
% 1- A vector containing the ensemble average force indentation for a given
%    protein\RNA condition - used in avgBinFitting
% 2- The residual between the ensemble average force indentation curve and
%    the Herztian contact model in loglog space.

clear all

loadedData_curves = load("C:\Users\Ben_G\OneDrive\Documents\MATLAB\MatLab\Indentation Analysis\IndentationDataSets\February_12_indentation_curves.mat"); % INPUT 1
loadedData_DropletParameters = load("C:\Users\Ben_G\OneDrive\Documents\MATLAB\MatLab\Indentation Analysis\AllDropletData\F12_AllDropletData.mat"); % INPUT 2
dataStruct_curves = loadedData_curves.DataStruct;
dataStruct_DropletParameters = loadedData_DropletParameters.dropletResultDataStruct;
DataStructSaveName = "F_12";

% saving commands
residualDataSave = 1;
binDataSave = 1;
radiusDataSave = 0;

% Initiliaze a data structure with all of the condensate geometrie for a
% given condition protein/RNA condition. 
radiusData = [];
dropletFieldNames = fieldnames(loadedData_DropletParameters.dropletResultDataStruct);

for i = 1:length(fieldnames(loadedData_DropletParameters.dropletResultDataStruct))

    radiusData = [radiusData,loadedData_DropletParameters.dropletResultDataStruct.(dropletFieldNames{i}).Var5];

end

if radiusDataSave == 1
     % Save the data into a .mat file
    save(DataStructSaveName + "_radiusData.mat", "radiusData");
    fprintf("Radius data saved to %s\n", DataStructSaveName + "_radiusData.mat");
end


% get the names of the different groups. Fieldname{i} = droplets{i}. these
% correspond to the same droplet.
fieldnames_list = fieldnames(dataStruct_curves);
droplets_list = fieldnames(dataStruct_DropletParameters);

allIndentation_validID = [];
allForce_validID = [];
allIndentation = [];
allForce = [];
allDropletPreStress = []; % this array will be populated with prestress values for each [indentation, force]. So all [indentation, force] for a given droplet will be assigned the same PreStress value
allLogIndentation_validID = [];
allLogForce_validID = [];

% The following section constructs an ensemble average force indentation
% curve for a given droplet by aliging and averaging all of the sequential
% indentation curves. 
for i = 1:length(fieldnames_list)
    experimentName = fieldnames_list{i}; % experimentName --> protein/RNA condition
    dropletName = droplets_list{i}; % dropletName --> droplet ID from experiment
    curves = dataStruct_curves.(experimentName); % curves --> set of force indentation curves extracted from "Extracting_Saving_Indentation_Curves" that are associated with that droplet.
    
    if length(curves) == 1
        % only consider condensates that have more than one sequential
        % indentation curves.
        continue
    else

        droplet_data = dataStruct_DropletParameters.(dropletName); % get the fitted parameters for that droplets and those curves
        droplet_radius = droplet_data.Var5; 
        droplet_prestress = droplet_data.Var6;
    
        for j = 1:length(curves)
    
            force = smooth(curves(j).force);
            indentation = smooth(curves(j).indendation); % typo in the curves group names "indendation" == "indentation"
    
            % Creating a dummy array that will be used to populate the
            % allDropletPreStress array
            preStress_dummy = zeros(length(force),1);
            preStress_dummy(:) = droplet_prestress;
    
            % Creates a vector containing all of the data points for a given
            % experimental conditions
            allForce = [allForce; force];
            allIndentation = [allIndentation; indentation];
            allDropletPreStress = [allDropletPreStress;preStress_dummy];
    
            % aligns the force curves with zero and only consider the indexes
            % that are above zero for the log log plot. These vectors will be
            % used to determine the fit of the linear region of FIC
            validIdx = force > 0 & indentation > 30 & indentation < 100 ;
            force = force(validIdx); % pN 
            indentation = indentation(validIdx); % nm
            
            % creating an vector thatcombines all the force indentation curves
            % for later fitting
            allIndentation_validID = [allIndentation_validID; indentation]; % pN
            allForce_validID = [allForce_validID; force]; % nm (should be on the order of 100s of nm)
            allLogIndentation_validID = [allLogIndentation_validID; log10(indentation)];
            allLogForce_validID = [allLogForce_validID; log10(force)];
        end
    end
end

% Sorting the vectors so that they all match up when plotting
[sortedIndentation,sortIdx] = sort(allIndentation);
sortedExperimentalForce = allForce(sortIdx);
sortedDropletPreStress = allDropletPreStress(sortIdx);

% Calculate the avg force value over a range - will be used for fitting
[binCenters,avgForces,stdForces,semForces] = binIndentationData(sortedIndentation,sortedExperimentalForce,10,7);
[slope,intercept]  = loglogfitting_ensemble(binCenters,avgForces);

if binDataSave == 1
     % Save the data into a .mat file
    save(DataStructSaveName + "_binData.mat", "binCenters", "avgForces", "stdForces"); % OUTPUT 1
    fprintf("Binning data saved to %s\n", DataStructSaveName + "_binData.mat");
end

interpolateLogForce = Interpolation(binCenters,slope,intercept);

residuals = residuals(avgForces,interpolateLogForce);

if residualDataSave == 1
    save(DataStructSaveName + "_residuals.mat", "residuals"); % OUTPUT 2
end