% This function computes two arrays of points PointsLeft, PointsRight near
% the fault between the classes stored in ClassesOfLeft and ClassesOfRight
% by bisection. The points in PointsLeft, PointsRight belong to the classes
% stored in ClassesOfLeft and ClassesOfRight, respectively. Hereby, the
% mutual distance of a point in PointsLeft and its counterpart in
% PointsRight is at most 2*FaultApproxParams.abstolBisection.
% This function assumes two arrays of initial points PointsLeftStart and
% PointsRightStart with points belonging to the classes stored in
% ClassesOfLeft and ClassesOfRight, resp.
% If the search for a certain point i was finished sucessfully, we set
% Finished(i) = true, otherwise false. In the latter case, the
% corresponding entries in PointsLeft and PointsRight are dummies only, so
% do not use them.
% If during bisection a third class occurs, we stop the computation of that
% point pair and proceed. We set Finished(i) = false in this case.
%
% Input:
% - PointsLeftStart: first of the two arrays of starting points
% - PointsRightStart: second of the two arrays of starting points
% - ClassesOfLeft: class values of PointsLeftStart
% - ClassesOfRight: class values of PointsRightStart
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.
% - FaultApproxParams: structure containing all parameters relevant for
%   the algorithm. We refer to FaultApproxParameters.m for details.
%
% Output:
% - PointsLeft: final points which belong to classOfLeft
% - PointsRight: final points which belong to classOfRight
% - Finished: Finished(i) is true, if a pair (PointsLeft(i,:),
%   PointsRight(i,:)) could be constructed; false if not.

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [PointsLeft, PointsRight, Finished] = ...
    computeSurfacePoints(PointsLeftStart, PointsRightStart, ...
                         ClassesOfLeft, ClassesOfRight, ...
                         ProblemDescr, FaultApproxParams)

    global ncalls
    
    [numPoints, ndim] = size(PointsLeftStart);        
    
    % stopping criterion for bisection algorithm
    abstolBisection = FaultApproxParams.abstolBisection;
    
    % maximum number of bisection iterations
    maxiterBisection = FaultApproxParams.maxiterBisection;

    Finished = zeros(numPoints,1);
    FinishedAux = zeros(numPoints,1);

    PointsMid = zeros(numPoints,ndim);
    ClassesOfMid = zeros(numPoints,1);
    
    anotherclass = false;
    PointsLeft = zeros(numPoints, ndim);
    PointsRight = zeros(numPoints, ndim);
    iiter = 1;
    
    while (iiter <= maxiterBisection)

        Finished = (vecnorm((PointsLeftStart -PointsRightStart)') < ...
                    2*abstolBisection)';
        FinishedAux = Finished | FinishedAux;
        
        if (all(FinishedAux))

            % add points to the list of surface points and break
            % note that both points are in different classes
            PointsLeft = PointsLeftStart;
            PointsRight = PointsRightStart;
            
            if anotherclass
                warnMessage = ['Bisection for at least one point' ...
                               ' in computeSurfacePoints failed' ...
                               ' due to a third class.'];
                warning(warnMessage)
            end
                        
            return
        else
        
            PointsMid(~FinishedAux,:) = 0.5*(PointsLeftStart(~FinishedAux,:) + ...
                                             PointsRightStart(~FinishedAux,:));

        
            ClassesOfMid(~FinishedAux) = ...
                computeClassification(PointsMid(~FinishedAux,:), ...
                                      ProblemDescr);
 
            % same class as left point: fault line is on the right
            PointsLeftStart((~FinishedAux & (ClassesOfMid == ClassesOfLeft)),:) = ...
                PointsMid((~FinishedAux & (ClassesOfMid == ClassesOfLeft)),:);
            PointsRightStart((~FinishedAux & (ClassesOfMid == ClassesOfRight)),:) = ...
                PointsMid((~FinishedAux & (ClassesOfMid == ClassesOfRight)),:);
 
            % sort out all points where a third class is involved
            IidxAnotherClass = ~FinishedAux & ...
                               (ClassesOfMid ~= ClassesOfLeft) & ...
                               (ClassesOfMid ~= ClassesOfRight);
            if (any(IidxAnotherClass))
                anotherclass = true;
            end
            FinishedAux(IidxAnotherClass) = 1;
         
        end
        iiter = iiter + 1;
    end
    
    if ~all(Finished)
        warnMessage = 'Bisection for at least one point in computeSurfacePoints failed ';
        warnMessage = [warnMessage ', consider enlarging maxIterBisection.'];
        warning(warnMessage)
    end
end