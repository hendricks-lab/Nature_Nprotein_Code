% Ben G. 2024

% DETECTION OF INDENTATION EVENTS:
% The time series force and trap position data are sectioned into
% sub-sections of lenght N. Indentation events are initiated by the Trap displacement -
% causing an increase in the Trap position data. To detect the begining of
% an indentation event, each trap position sub-section is evaluated as
% having or not having a significant increase in position. If a sub-section
% is determined to have a significant Trap dispclament - the time label
% corresponding to the start of the sub-section is marked as being the
% begining of an indentation event. That same time label is then applied to
% the force data. 

% FORCE INDENTATION EXTRACTION:
% Polystyrene probes were trapped and calibrated near the condensates - the
% instrument calibration assumed passive thermal motion close to a boundary
% (the coverslip). Using the trap calibration, force measurement, and trap
% position - we determine the condesate indentation as follows --> 
% Xind = Xtrap - deltaXforce. 
% This analysis is performed for every detected indentation event

% HERTZIAN CONTACT MODEL FITTING OF INDIVIDUAL FORCE INDENTATION CURVES:
% A Hertzian contact model if fit gor every force indentation curve
% identified at an indentation event. Sequential measurement are performed
% on a single condensate - consequently, each condensate has multiple
% indentation events and fitted force indentation curves. 

% INPUTS:
% 1- H5 file of CTrap force data and calibration
% 2- CSV files of N protein condensate geometries (radius)
%
% OUTPUTS:
% 1- An experiment data sructure containing the fitted parameters of all
% indentation events grouped by their respective condensate:
% data structure --> condensate ID --> fitted parameters of individual
% force indentation curves

clear all

directoryPath = "C:\Users\Ben_G\OneDrive\Documents\MatLab\DATA\jan_19_analysis"; % INPUT 1
dropletRadiusDirectory = "C:\Users\Ben_G\OneDrive\Documents\MatLab\DATA\jan_19_analysis\DropletRadius\January_19_droplet_radii"; % INPUT 2
keyword = "measurement";

radiusMatrix = readmatrix(dropletRadiusDirectory);
files = dir(directoryPath);
files = files(~[files.isdir]); 
probeDiameter = 0.5; %

%% Create the data structure that will hold all of parameter values associated with each droplet
numGroups = length(files);
numParameters = 6; % this can be changed based on the number of the indentation curves that we opt for

dropletResultDataStruct = struct();
allDropletYoungModulusFit = [];
allDropletLogLogFit = [];
allDropletLogLogTransition = [];
allDropletRadius = [];
averageModulus = [];
trapStiffness = [];
dropletPreStress = [];
maxDropletStiffness = [];

for i = 1:numGroups
    groupName = sprintf('Group%d', i);
    dropletResultDataStruct.(groupName) = struct();

    for j = 1:numParameters
        varName = sprintf('Var%d', j);
        dropletResultDataStruct.(groupName).(varName) = NaN;
    end
end    

for k = 1:length(files)
    
    file = files(k).name;
   
    if contains(file,keyword)
        
        % information pertaining to the CTrap calibration associated with the force measurements. 
        fullFilePath = fullfile(directoryPath, file);
        h5disp(fullFilePath)
        infoStructure = h5info(fullFilePath);
        calibrationName = infoStructure.Groups(2).Groups.Name;
        group_name = sprintf('Group%d', k); % assign the group name for the data structure
    
        % extracting relevant force, position, and calibration data from
        % the H5 file.
        force_1x = h5read(fullFilePath, '/Force HF/Force 1x');
        num_points = length(force_1x);
        sampling_rate = h5readatt(fullFilePath,'/Force HF/Force 2x','Sample rate (Hz)');
        time_step = 1/sampling_rate;
        t = [0:time_step:0+(num_points-1)*time_step];
    
        xAOD = h5read(fullFilePath, '/Trap position/1X'); 
        yAOD = h5read(fullFilePath, '/Trap position/1Y');
        offset = h5readatt(fullFilePath, string(calibrationName) + '/Force 1x','Offset (pN)');
       
        RF_1x = h5readatt(fullFilePath, string(calibrationName) + '/Force 1x','Response (pN/V)');
        RF_1y = h5readatt(fullFilePath, string(calibrationName) + '/Force 1y','Response (pN/V)');
        TK_calibrated = h5readatt(fullFilePath, string(calibrationName) + '/Force 1x','kappa (pN/nm)');
        trapStiffness = [trapStiffness, TK_calibrated];

        if TK_calibrated > 0.08
            disp("test");
            continue
        end
        
        % set to "true" for visual confirmation that the force and trap
        % positon data is correctly being selected
        if true == false 
            figure(1)
            plot(t, movmean(force_1x,5000));
            ylabel('Force (pN)');
            plot(t, xAOD);
            ylabel('Trap position (um)');
            xlabel("Time (s)");
            title("Smooth Zeroed Force and Trap Position");
        end
    
        % Automated detection of sequential indentation events collected within a single H5 file.
        trapDisplacement = 0.38; % um
        displacementTime = 0.48; % seconds
        numIndices = round(displacementTime * sampling_rate);
   
        possibleStartTimes = [];
        possibleEndTimes = [];
        chuncksize = 3500; % can be fine tuned. 
    
        % Identification of possible start and end times for a force
        % indentation curve
        for i = 1:chuncksize:(length(t) - numIndices)
            increase = xAOD(i + numIndices) - xAOD(i);
           
            % Evaluates if the increase over a sub-section is sufficient to
            % be considered an indentation event
            if increase >= trapDisplacement
                possibleStartTimes = [possibleStartTimes, t(i)];
                possibleEndTimes = [possibleEndTimes, t(i + numIndices)];
            end
        end
    
        % select the actual start and end times
        trueStartTimes = [];
        trueEndTimes = [];
        for i = 2:3:length(possibleEndTimes)
            trueStartTimes = [trueStartTimes, possibleStartTimes(i)];
            trueEndTimes = [trueEndTimes, possibleEndTimes(i)];
        end 
    
        % Fit a Hertzian contact model to every detected indentation event
        % for a given condensate.
        dropletRadius = radiusMatrix(k); % condensate geometry --> radius
        DropletYoungModulusFit_IndividualCurves = []; % initialize the Young's Modulus array that will store the moduli for the individual force indentation curves

        for i = 1:length(trueEndTimes)
           
            var_name = sprintf('Var%d', i);

            window_size = 5000;
            start = trueStartTimes(i);
            ends = trueEndTimes(i);
            
            zeroing_position = xAOD(round(start*sampling_rate) - 100);
            zeroing_force_average = mean(force_1x(round(start*sampling_rate)-5000:round(start*sampling_rate)-500));
            zeroed_position = xAOD(round(start*sampling_rate)-5000:round(ends*sampling_rate)+5000) - zeroing_position;
            zeroed_force = force_1x(round(start*sampling_rate)-5000:round(ends*sampling_rate)+5000) - zeroing_force_average;
            time = ([0:(1/sampling_rate):(length(zeroed_force) -1)*(1/sampling_rate)]);
            smoothed_data = movmean(zeroed_force, window_size);
           
            indent = 1;
            if indent == 1
                center_disp = smoothed_data/TK_calibrated; % nm
                droplet_indentation = (zeroed_position*1000) - center_disp;
            end

            [E_data_hertz, dropletYoungModulusFit, x, y, eff_R] = fitHertzModel(droplet_indentation, smoothed_data, probeDiameter, dropletRadius, 0.45); % Fit the Hertzian model to every indentation event
            DropletYoungModulusFit_IndividualCurves = [DropletYoungModulusFit_IndividualCurves, dropletYoungModulusFit]; % Add the fitted parameters to an array - each entry corresponds to an indentation event.
     
            dropletResultDataStruct.(group_name).(var_name) = dropletYoungModulusFit; % add the fitted parrameters to a variable in a data structure - will be used in subsequent scripts

        end

        averageModulus = mean(DropletYoungModulusFit_IndividualCurves); % Determine the mean Young's Modulus across all indentation events for a given condensate.
        
        allDropletYoungModulusFit = [allDropletYoungModulusFit,averageModulus];
        
        % Determine the prestress on a droplet
        if true == true
            averageEndForce = mean(force_1x((length(force_1x)-5000):(length(force_1x))));
            dropletPreStress = [dropletPreStress, (zeroing_force_average - averageEndForce)];
            preStress = zeroing_force_average - averageEndForce;
            dropletResultDataStruct.(group_name).(sprintf('Var%d',numParameters)) = preStress;
        end
    end

    if true == false
        saveFileName = "DropletResult.mat";
        saveFilePath = fullfile(directoryPath, saveFileName);
        save(saveFilePath,'dropletResultDataStruct'); % OUTPUT 1
    end

end
