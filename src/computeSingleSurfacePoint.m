% This function computes two points PointLeftFinal, PointRightFinal
% near the fault between classes classOfLeft and classOfRight by
% bisection. The new points belong to to class classOfLeft and
% classOfRight, resp. and have a mutual distance of at most
% 2*FaultApproxParams.abstolBisection.
% This function assumes two initial points PointLeftStart and
% PointRightStart which belong to class classOfLeft and classOfRight,
% resp.
% If the search was finished sucessfully, this is indicated by
% finished = true, otherwise false. In the latter case, PointLeftFinal
% and PointRightFinal are dummies only, so do not use them.
%
% Input:
% - PointLeftStart: first of the two starting points
% - PointRightStart: second of the two stating points
% - classOfLeft: class value of PointLeftStart
% - classOfRight: class value of PointRightStart
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% - PointLeftFinal: final point which belongs to classOfLeft
% - PointRightFinal: final point which belongs to classOfRight
% - finished: true, if a pair (PointLeftFinal, PointRightFinal) could be
%   constructed; false if not.

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [PointLeftFinal, PointRightFinal, finished] = ...
    computeSingleSurfacePoint(PointLeftStart, PointRightStart, ...
                              classOfLeft, classOfRight, ProblemDescr, ...
                              FaultApproxParams)

    % stopping criterion for bisection algorithm
    abstolBisection = FaultApproxParams.abstolBisection;
    
    % maximum number of bisection iterations
    maxiterBisection = FaultApproxParams.maxiterBisection;

    finished = false;
    PointLeftFinal = 0;
    PointRightFinal = 0;
    normError = norm(PointLeftStart -PointRightStart);
    iiter = 1;
    
    while (iiter <= maxiterBisection)
        
        % The desired accuracy FaultApproxParams.abstolBisection has been
        % met for the mean: we are done.
        if (normError < 2.0*abstolBisection)

            % add points to the list of surface points and return
            % note that both points are in different classes
            PointLeftFinal = PointLeftStart;
            PointRightFinal = PointRightStart;
            finished = true;
            
            % In this case, we do not need to think about issuing a
            % warning, but can go straight back to the calling function.
            return
        else
        
            % new midpoint
            PointMid = 0.5*(PointLeftStart + PointRightStart);

            % the interval length is halved
            normError = 0.5*normError;

            classOfMid = computeClassification(PointMid, ProblemDescr);

            % Same class as left point: The fault crosses the right part
            % of the interval.
            if(classOfMid == classOfLeft)
                PointLeftStart = PointMid;
                classOfLeft = classOfMid;

            % Same class as right point: The fault crosses the left part
            % of the interval.
            elseif(classOfMid == classOfRight)
                PointRightStart = PointMid;
                classOfRight = classOfMid;

            % This case indicates that there are at least two sections with
            % faults on the line from Left to Right. We skip search
            % then and return back to the calling function.
            else
                warnMessage = ['Bisection in computeSingleSurfacePoint failed, res = ' ...
                               num2str(normError)...
                               ', due to a third class ' ...
                               int2str(classOfMid)];
                warning(warnMessage)
                return
            end
        end
        
        iiter = iiter+1;
    end
    
    warnMessage = ['Bisection in computeSingleSurfacePoint failed, res = ' ...
                   num2str(normError) ...
                   ', consider enlarging maxiterBisection.'];
    warning(warnMessage)
end
