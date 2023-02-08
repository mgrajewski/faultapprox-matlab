% This function creates a set of points using the Halton sequence for
% each coordinate in [0,1]^dim. For details, see e.g.
% https://en.wikipedia.org/wiki/Halton_sequence
%
% Input:
% - numSamples: desired number of sampling points
% - dim: dimension of hypercube [0,1]^dim to sample
% - seed: seed for Halton sequence (useful if one wnts to continue with an
%   existing Halton sequence)
%
% Output:
% - PointSet: (numSamples x dim)-array containing the cartesian
%   coordinates of the sampling points

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function PointSet = CreateHaltonSet(numSamples, dim, seed)

    % primes as starting point for each coordinate (here: implicitly
    % hard-coded to 8, as I do not see realistic scenarios with dim > 8)
    Primes = [2,3,5,7,11,13,17,19];

    if (nargin == 2)
        seed = 0;
    end
    
    % allocate the point set
    PointSet = zeros(numSamples, dim);
    
    for idim = 1: dim
        basis = Primes(idim);
    
        for isample = 1:numSamples    
            n0 = isample + seed;
            HaltonNumber = 0;
            f = 1/basis;
            while (n0>0)
                n1 = floor(n0/basis);
                r = n0-n1*basis;
                HaltonNumber = HaltonNumber + f*r;
                f = f/basis;
                n0 = n1;
            end
            PointSet(isample, idim) = HaltonNumber;
        end
    end
end