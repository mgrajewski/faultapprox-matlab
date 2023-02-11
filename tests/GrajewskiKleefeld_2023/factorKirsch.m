% Implementation of Andreas Kirsch's factorization method.

% Author: Andreas Kleefeld (a.kleefeld@fz-juelich.de)
% This file is part of faultapprox-matlab with Andreas Kleefeld's
% permission.
% (https://github.com/mgrajewski/faultapprox-matlab)
function val=factorKirsch(x, ProblemDescr)
    global M
    global d
    global V
    global sigma
    global wavenumber

    numPoints = size(x,1);
    val = zeros(numPoints,1);
    
    for i = 1: numPoints
        xaux = x(i,:);
        rz=zeros(M,1);
        for j=1:M
            rz(j)=exp(-1i*wavenumber*xaux*d((j-1)*M+1,:)');
        end
        rhoz=V'*rz;

        val(i)=1/sum(abs(rhoz).^2./abs(sigma));
        if val(i)<=0.026
            val(i)=1;
        else
            val(i)=2;
        end
    end
end