% This function creates pairs of points near a fault line, one of
% which is in class iclass and the other in class jclass. These
% point pairs and their respective classes are stored in
% PointPairs and ClassPointPairs. It serves for finding
% starting values for bisection. The function requires an initial guess
% for the point pairs given in PointsLeft and PointsRight.
% If a point in PointsLeft is in class iclass and its counterpart (aka the
% point at the same position in PointsRight) is in class jclass,
% fine. If not, we consider two cases: Either both points belong to
% the same class iclass or jclass. Then, we mirror these points in
% PointsRight at the respective points in PointsLeft and continue. If
% this does not lead to a suitable point pair, we enlarge the
% perturbation for PointsRight and continue up to maxNumTrials
% times.
% As second case, it may happen that one point of the pair neither
% belongs to iclass or jclass but to another class. This may happen if
% the fault line is about to leave the domain or if three or more
% classes meet on in one point. Then, we scatter along the
% normal to the point on the curve in both directions. If one of these
% points leads to a valid pair of points, fine. If not, we give up for
% that point pair.
%
% Input:
% PointsLeft: Points near a curve approximating some part of the
% fault line (usually generated by subdividing a line segment or
% polynomial extrapolation)
% PointsRight: Points near a curve, ideally on the other side of the
% curve as PointsLeft.
% ClassLeft, ClassRight: classes of the respective points
% iclass, jclass: indices of the classes the fault line is in between
% - ProblemDescr: structure containing all problem-relevant parameters.
%   We refer to its documentation in ProblemDescr.m for details.

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [PointPairs, Succeeded, ClassPointPairs] = startPairs(PointsLeft, ...
        PointsRight, ClassLeft, ClassRight, iclass, jclass, ProblemDescr)

    
    % maximum number of trials in case of failure mode 1 (both points
    % belong to the same class) 
    maxNumTrials = 3;

    epsLoc = 1e-10;
        
    % dimension of data
    ndim = size(PointsLeft, 2);

    % array of vaild point pairs
    PointPairs = zeros(size(PointsLeft,1),2*ndim);
        
    % corresponding classes
    ClassPointPairs = zeros(size(PointsLeft,1),2);

    % true, if for the corresponding point pair a valid starting pair could
    % be found (initially false)
    PointsOldEven = PointsLeft;
    
    % finding valid pairs my fail due to presence of a third class very
    % near the fault line. However, this may turn out only during the first
    % phase and cannot by derived from the initial choice of points.
    % We mark all points where a third class occurred during the first
    % phase. 
    ThirdClassInvolved = zeros(size(PointsLeft,1),1);
    
    PointsRightOrig = PointsRight;

    % first pass: consider case both classes are the same

    % If a point is in class iclass and its counterpart (aka the
    % point at the same position in PointsRight) is in class
    % jclass or vice versa: search of point pair succeeded
    Succeeded = or( ClassLeft == iclass & ClassRight == jclass, ...
                    ClassLeft == jclass & ClassRight == iclass ...
                  );

    % we mark all points where at some point in the algorithm a third
    % class has been involved: It may happen that initially, a third
    % class has been involved, but at the current stage, the points of
    % the pair belong to the same class. This case requires a different
    % treatment than the case where the both points of the starting
    % pair just belong to iclass or jclass. If initially a third class
    % has been involved, but by coincide, a suitable point pair has
    % been found, then this pair is marked as succeeded and is excluded
    % from any further treatment anyway.
    ThirdClassInvolved = ThirdClassInvolved | (ClassLeft ~= iclass & ClassLeft ~= jclass) | (ClassRight ~= iclass & ClassRight ~= jclass);

    % if StillTheSameClass, both points belong to the same class iclass
    % or jclass. Only in this case, the expansive treatment in the first
    % part of this algorithm makes sense.
    StillTheSameClass = (ClassLeft == ClassRight) & (ClassLeft == iclass | ClassLeft == jclass);

    PointPairs(Succeeded, :) = [PointsLeft(Succeeded,:), PointsRight(Succeeded,:)];
    ClassPointPairs(Succeeded, :) = [ClassLeft(Succeeded), ClassRight(Succeeded)];

    % if for all points class change takes place:
    % finish, if not: try opposite direction for remaining points
    if any(StillTheSameClass)
        
        for j = 1:maxNumTrials

            % try opposite direction: mirror on line, unless one leaves the
            % domain
            
            % This is to choose optimal point pairs only
            PointsOldOdd = PointsRight;

            PointsRight(StillTheSameClass,:) = 2*PointsLeft(StillTheSameClass,:) - PointsRight(StillTheSameClass,:);

            % restrict the points to the domain
            PointsRight(StillTheSameClass,1) = min(max(ProblemDescr.Xmin(1) + epsLoc, PointsRight(StillTheSameClass,1)), ProblemDescr.Xmax(1)-epsLoc);
            PointsRight(StillTheSameClass,2) = min(max(ProblemDescr.Xmin(2) + epsLoc, PointsRight(StillTheSameClass,2)), ProblemDescr.Xmax(2)-epsLoc);
            ClassRight(StillTheSameClass) = computeClassification(PointsRight(StillTheSameClass,:), ProblemDescr);

            SucceededNow = or(ClassLeft == iclass & ClassRight == jclass, ...
                              ClassLeft == jclass & ClassRight == iclass) & ~Succeeded;

            StillTheSameClass = (ClassLeft == ClassRight) & (ClassLeft == iclass | ClassLeft == jclass);
            PointPairs(SucceededNow, :) = [PointsOldEven(SucceededNow,:), PointsRight(SucceededNow,:)];
            ClassPointPairs(SucceededNow, :) = [ClassLeft(SucceededNow), ClassRight(SucceededNow)];
            Succeeded = Succeeded | SucceededNow;
            ThirdClassInvolved = ThirdClassInvolved | (ClassLeft ~= iclass & ClassLeft ~= jclass) | (ClassRight ~= iclass & ClassRight ~= jclass);

            % if for all points class change between iclass and jclass takes place:
            % finish, if not: double distance to straight line
            if ~any(StillTheSameClass)
                break
            end
                
            % PointsOdd
            PointsOldEven = PointsRight;

            % if still points are unfinished: go into the opposite
            % direction
            PointsRight(StillTheSameClass,:) = 3*PointsLeft(StillTheSameClass,:) - 2*PointsRight(StillTheSameClass,:);
            PointsRight(StillTheSameClass,1) = min(max(ProblemDescr.Xmin(1) + epsLoc, PointsRight(StillTheSameClass,1)), ProblemDescr.Xmax(1)-epsLoc);
            PointsRight(StillTheSameClass,2) = min(max(ProblemDescr.Xmin(2) + epsLoc, PointsRight(StillTheSameClass,2)), ProblemDescr.Xmax(2)-epsLoc);
            ClassRight(StillTheSameClass) = computeClassification(PointsRight(StillTheSameClass,:), ProblemDescr);

            SucceededNow = or(ClassLeft == iclass & ClassRight == jclass, ...
                              ClassLeft == jclass & ClassRight == iclass) & ~Succeeded;

            StillTheSameClass = (ClassLeft == ClassRight) & (ClassLeft == iclass | ClassLeft == jclass);
            PointPairs(SucceededNow, :) = [PointsOldOdd(SucceededNow,:), PointsRight(SucceededNow,:)];
            ClassPointPairs(SucceededNow, :) = [ClassLeft(SucceededNow), ClassRight(SucceededNow)];
            Succeeded = Succeeded | SucceededNow;
            ThirdClassInvolved = ThirdClassInvolved | (ClassLeft ~= iclass & ClassLeft ~= jclass) | (ClassRight ~= iclass & ClassRight ~= jclass);        

            if ~any(StillTheSameClass)
                break
            end
        end
    end
    
    % If there are any points still unfinished due to a third class
    % involved, we try scattering around the midpoint of the original point
    % pair along the approximate normal.
    IdxFailed2 = ~Succeeded & ThirdClassInvolved;
    if any(IdxFailed2)
    
        % go back to the original points ()PointsLeft has never been
        % touched)
        PointsRight(IdxFailed2,:) = PointsRightOrig(IdxFailed2,:);
        IidxVec = find(IdxFailed2 & ~Succeeded ==1);

        % scatter around the normal to the points
        pars = [-0.8, -0.4, -0.2 0, 0.2, 0.4, 0.8];
        parsize = size(pars, 2);
        AuxLeft = PointsLeft(IdxFailed2,:);
        AuxRight = PointsRight(IdxFailed2,:);
        Normals = AuxLeft - AuxRight;
        
        AuxPoints = zeros(parsize*size(AuxLeft, 1), ndim);
        for iloc = 1: size(AuxLeft, 1)
            AuxPoints((iloc-1)*parsize+1:iloc*parsize, :)  = 0.5*(AuxLeft(iloc, :) + AuxRight(iloc,:)) + pars'*Normals(iloc,:);
        end

        % indices of the auxiliry points inside the domain
        idxInside = all(AuxPoints >= ProblemDescr.Xmin & AuxPoints <= ProblemDescr.Xmax, 2);

        % we classify only the points inside the domain
        ClassAux = -ones(1, size(AuxPoints,1));
        ClassAux(idxInside) = computeClassification(AuxPoints(idxInside,:), ProblemDescr);

        for iloc = 1: size(AuxLeft, 1)

            ClassAuxLoc = ClassAux((iloc-1)*parsize+1:iloc*parsize);
            AuxPointsLoc = AuxPoints((iloc-1)*parsize+1:iloc*parsize, :);
            
            % there is a suitable point pair
            if (any(ClassAuxLoc == iclass) && any(ClassAuxLoc == jclass))

                % This is not optimal: We just take one suitable pair and
                % not the one with the nearest points if there are more
                % than one
                idxIclassMin = find(ClassAuxLoc == iclass, 1, 'first');
                idxJclassMin = find(ClassAuxLoc == jclass, 1, 'first');

                idxIclassMax = find(ClassAuxLoc == iclass, 1, 'last');
                idxJclassMax = find(ClassAuxLoc == jclass, 1, 'last');

                % find points in the Aux vector with the correct
                % classes
                if (idxIclassMax < idxJclassMin)
                    idxIclass = idxIclassMax;
                    idxJclass = idxJclassMin;
                else
                    idxIclass = idxIclassMin;
                    idxJclass = idxJclassMax;
                end
                PointPairs(IidxVec(iloc),:) = [AuxPointsLoc(idxIclass,:), AuxPointsLoc(idxJclass,:)];
                ClassPointPairs(IidxVec(iloc),:) = [iclass, jclass];
                Succeeded(IidxVec(iloc)) = 1;
            end
        end
    end
end
