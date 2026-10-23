function [interpolatedLogForce] = Interpolation(bincenters,slope,intercept)
    
    interpolatedLogForce = slope.*(log10(bincenters)) + intercept;
    
end

