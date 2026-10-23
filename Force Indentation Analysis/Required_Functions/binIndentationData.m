function [binCenters,avgForces,stdForces,semForces] = binIndentationData(indentation,force,binwidth,binextension)
    % BININDENTATIONDATA Bins indentation-force pairs and computes the average force in each bin
    %
    % Inputs:
    %   indentation - Array of indentation values (1D numeric array)
    %   force       - Array of corresponding force values (1D numeric array)
    %
    % Outputs:
    %   binCenters  - Centers of the bins (1D numeric array)
    %   avgForces   - Average force in each bin (1D numeric array)
    
    % Ensure inputs are column vectors
    indentation = indentation(:);
    force = force(:);

    % Define the bin centers (10 to 400 in nanometer steps)
    binCenters = 10:binwidth:350; % Bin centers from 10 to 400 nm
    
    % Preallocate array for average forces
    avgForces = nan(size(binCenters)); % Use NaN to mark bins with no data
    stdForces = nan(size(binCenters)); % Standard deviation
    semForces = nan(size(binCenters)); % Standard error of the mean
    
    % Loop through each bin
    for i = 1:length(binCenters)
        % Define bin range: [binCenter - 0.5, binCenter + 0.5]
        binStart = binCenters(i) - binextension;
        binEnd = binCenters(i) + binextension;
        
        % Find indices of indentation values within the current bin
        binIndices = (indentation >= binStart) & (indentation < binEnd);
        
        % Compute the average force for the current bin
        if any(binIndices) % Check if there are points in the bin
            avgForces(i) = mean(force(binIndices));
            stdForces(i) = std(force(binIndices)); % Standard deviation
            semForces(i) = stdForces(i) / sqrt(sum(binIndices)); % Standard error of the mean
        end
    end
end

