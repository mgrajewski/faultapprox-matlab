% This function computes the centre and the radius of a circle defined by
% three points and returns the reciprocal of the radius as curvature. 
% The corresponding algorithm can be derived as follows:
% Let be a = (x1, y1), O = (x2, y2), b = (x3,y3) three distinct points
% and
%
%        / x^2 +y^2   x  y  1\
%        | x1^2+y1^2  x1 y1 1|
%    A = | x2^2+y2^2  x1 y2 1|  .
%        \ x3^2+y3^2  x3 y3 1/
%
% Then, direct computation by a Laplace expansion of A along the
% first row reveals that {det(A) = 0} is a circle through these points
% around (d2/2d1, -d3/2d1) with radius
%    r = sqrt((d2^2+d3^2)/(4d1^2) + d4/d1), where
% 
%         |x1 y1 1|        |x1^2+y1^2 y1 1|         |x1^2+y1^2 x1 1| 
%    d1 = |x2 y2 1|,  d2 = |x2^2+y2^2 y2 1|,   d3 = |x2^2+y2^2 x2 1|,
%         |x3 y3 1|        |x3^2+y3^2 y3 1|         |x3^2+y3^2 x3 1| 
% 
%         |x1^2+y1^2 x1 y1|
%    d4 = |x2^2+y2^2 x2 y2|
%         |x3^2+y3^2 x3 y3|
% Its reciprocal is the radius we are looking for.
% However, we apply this function to shifted coordinates such that
%(x2,y2) = 0. Therefore, all this simplifies to
%    d1 = <a', b> with a' = (y1, -x1),
%    d2 = y1 ||b||^2 - y3 ||a||^2,
%    d3 = x1 ||b||^2 - x3 ||a||^2,
%    d4 = 0,
% and ultimately
%
%    curvature = 1/r = (2d1)/(||a||^2 ||b||^2 ||a-b||)              (*)
%
% For the angle alpha between two vectors a and b, we have
%    cos alpha = <a,b>/(||a|| ||b||) or
%    cos(90- alpha) = <a', b>/(||a|| ||b||) =  0.5 d1
%
% Apart of scaling, computing the curvature with (*) becomes unstable,
% if both the nominator and the denominator tend to zero. Our analysis
% reveals that this happens, if approx. a = b (note that a' is
% perpendicular to a).
%
% Input:
% - threePoints: The columns contain the cordinates of the three points
%   The second point MUST be 0! It would have been sufficient to define the
%   first and third point only, but it is more convenient to pss all
%   three points.
%
% Output:
% - curvature: the reciprocal of the radius of the circle
% - center: the center of the circle

% Author: Matthias Grajewski (grajewski@fh-aachen.de)
% This file is part of faultapprox-matlab
% (https://github.com/mgrajewski/faultapprox-matlab)
function [curvature, center] = curvCenterFromThreePoints(threePoints)

    a = threePoints(1,:);
    b = threePoints(3,:);

    d1 = a(2)*b(1) - a(1)*b(2);
    d2 = a(2)*(b*b') - b(2)*(a*a');
    d3 = a(1)*(b*b') - b(1)*(a*a');

    % The points are not (approximately) on a straight line: It makes sense
    % to compute a center of the circle.
    if (d1 > 1e-8)
        center = [d2/(2*d1), -d3/(2*d1)];

    % The notion of a circle does not make too much sense. As we have to
    % return anything, we return a vector with zeroes.
    else
        % This makes the center have the right size as array.
        center = 0*a;
    end
    curvature = 2.0*abs(d1)/(norm(a-b, 2) * norm(a, 2) * norm(b, 2));
end