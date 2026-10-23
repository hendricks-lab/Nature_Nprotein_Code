function [residuals] = residuals(experimental,interpolated)

    residuals = abs(log10(experimental) - interpolated);
    
end

