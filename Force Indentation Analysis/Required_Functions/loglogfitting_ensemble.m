function [slope,intercept] = loglogfittingensemble(indentation,force)
    % LOGLOGFITTINGENSEMBLE Fits log-log data for indentations under 100 nm
    %
    % Inputs:
    %   indentation - Array of indentation values (1D numeric array)
    %   force       - Array of corresponding force values (1D numeric array)
    %
    % Outputs:
    %   slope       - Slope of the log-log fit
    %   intercept   - Intercept of the log-log fit

    validIndices = 25 < indentation & indentation < 100;
    filteredIndentation = indentation(validIndices);
    filteredForce = force(validIndices);

     % Take the log of the filtered data
    log_indentation = log10(filteredIndentation);
    log_force = log10(filteredForce);
    
    p = polyfit(log_indentation,log_force,1);
    slope = p(1);
    intercept = p(2);
end

