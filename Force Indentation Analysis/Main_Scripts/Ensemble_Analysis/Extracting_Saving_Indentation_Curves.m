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

% INPUTS:
% 1- H5 file of CTrap force data and calibration
% 2- CSV files of N protein condensate geometries (radius)
%
% OUPUTS:
% 1- A 1x1 data structure. Containing N 1xM data structures - where N is the
% number of condensates measured for a given condition and M is the number of
% indentation events detected. Each indentation event is composed of a 1xS
% force and trap position vector - where S is the number of data point
% collected by the LUMICKS system.
%

clear all

% Data directories 
directoryPath = "C:\Users\Ben_G\OneDrive\Documents\MATLAB\MatLab\DATA\mar_5_analysis"; % INPUT 1
dropletRadiusDirectory = "C:\Users\Ben_G\OneDrive\Documents\MATLAB\MatLab\DATA\mar_5_analysis\DropletRadius\March_05_droplet_radii"; % INPUT 2
keyword = "measurement"; 

radiusMatrix = readmatrix(dropletRadiusDirectory);
files = dir(directoryPath);
files = files(~[files.isdir]); 

% Indentation event data structure
DataStruct = struct(); % OUTPUT 1

for k = 1:length(files)
    
     file = files(k).name;
     experimentName = sprintf("Experiment_%d", k);
     DataStruct.(experimentName) = struct();

     if contains(file,keyword)
        
        % information pertaining to the CTrap calibration associated with the force measurements. 
        fullFilePath = fullfile(directoryPath, file);
        infoStructure = h5info(fullFilePath);
        calibrationName = infoStructure.Groups(2).Groups.Name;
        group_name = sprintf('Group%d', k); 
    
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
        
        if TK_calibrated > 0.08
            disp("test");
            continue
        end

        % Automated detection of sequential indentation events collected within a single H5 file.
        trapDisplacement = 0.38; % um
        displacementTime = 0.98; % seconds
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
    
        trueStartTimes = [];
        trueEndTimes = [];
        for i = 2:3:length(possibleEndTimes)
            trueStartTimes = [trueStartTimes, possibleStartTimes(i)];
            trueEndTimes = [trueEndTimes, possibleEndTimes(i)];
        end 
        
        % Extraction of the condensate indentation
        for i = 1:length(trueEndTimes)
           
            var_name = sprintf('Var%d', i);
    
            window_size = 5000;
            start = trueStartTimes(i);
            ends = trueEndTimes(i);
            
            zeroing_position = xAOD(round(start*sampling_rate) - 100); % the Trap position prior to the indentation event
            zeroing_force_average = mean(force_1x(round(start*sampling_rate)-5000:round(start*sampling_rate)-500)); % the force signal is zeroed relative to its position prior to the indentation event
            zeroed_position = xAOD(round(start*sampling_rate)-5000:round(ends*sampling_rate)+5000) - zeroing_position; % the Trap position is zeroed relative to its position prior to the indentation event
            zeroed_force = force_1x(round(start*sampling_rate)-5000:round(ends*sampling_rate)+5000) - zeroing_force_average;
            time = ([0:(1/sampling_rate):(length(zeroed_force) -1)*(1/sampling_rate)]);
            smoothed_data = movmean(zeroed_force, window_size); % data was smoothed using an averaging window to minimize high frequency noise. 
           
            indent = 1;
            if indent == 1
                center_disp = smoothed_data/TK_calibrated; % nm
                droplet_indentation = (zeroed_position*1000) - center_disp;
            end

            DataStruct.(experimentName)(i).force = smoothed_data;
            DataStruct.(experimentName)(i).indendation = droplet_indentation;
              
        end
    end
end

display("DONE - DataStruct has been saved");