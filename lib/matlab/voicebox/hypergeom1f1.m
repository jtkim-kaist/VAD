function [h,l,alg]=hypergeom1f1(a,b,z,tol,maxj)
% HYPERGEOM1F1 Confluent hypergeometric function, 1F1 a.k.a Kummer's M function [h,l]=(a,b,z,tol,maxj)
%
%  Inputs: a,b,z  input arguments
%                 a and b must be real scalars, z can be a real matrix
%          tol    Optional tolerance [default 1e-10]
%          maxj   Optional iteration limit [default 500]
%
% Outputs: h      output result (size of z): 1F1(a;b;z) = M(a;b;z)
%          l      actual number of iterations (size of z)
%
% This routine is functionally equivalent to the symbolic toolbox function hypergeom(a,b,z) but much faster.
% The function calculated is the solution to z M'' + (b-z) M' - a M = 0 where
% M' and M'' are 1st and 2nd derivatives of M with respect to z.
% This routine is closely based on taylorb1f1.m from [4] which is explained in Sec 3.2 of [3].
%
% Special cases [2]:
%        M(a;b;z) = Inf for integer b<=0
%        M(a;b;0) = 1
%        M(0;b;z) = 1
%        M(a;a;z) = exp(z)
%        M(1;2;2z) = exp(z)*sinh(z)/z
%        M(a;b;z) = exp(z)*M(b-a;b;-z)
%        M(a;a+1;-z) = exp(z)*M(1;a+1;z) = a*gamma(z,a)/z^a
%        M(a;b;z) = M(a-1;b;z) + M(a;b+1;z)*z/b
%
% References:
% [1]	M. Abramowitz and I. A. Stegun, editors.
%       Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables.
%       Dover, New York, 9 edition, 1964. ISBN 0-486-61272-4.
% [2]	F. W. J. Olver, D. W. Lozier, R. F. Boisvert, and C. W. Clark, editors.
%       NIST Handbook of Mathematical Functions. CUP, 2010. ISBN 978-0-521-14063-8.
% [3]	J. Pearson. Computation of hypergeometric functions. Master’s thesis, Oxford University, 2009.
%       URL http://¬people.maths.ox.ac.uk/¬porterm/¬research/-pearson_final.pdf.
% [4]	J. Pearson. Computation of hypergeometric functions. Matlab code, Oxford University, 2009.
%       URL http://¬people.maths.ox.ac.uk/¬porterm/¬research/-hypergeometricpackage.zip.
%

%      Copyright (C) Mike Brookes 2016
%      Version: $Id: hypergeom1f1.m 8911 2016-10-28 14:58:42Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5
    maxj=500;
    if nargin<4
        tol=1e-10;
    end
end

h=zeros(size(z));
l=zeros(size(z));
nz=numel(z);
for i=1:nz
    y=z(i);
    q=0;  % break criterion initially false
    if abs(y)<10            % small |y|
        alg=1;
        d=y^2*(a*(a+1)/(b*(2*b+2)));
        g=1+y*a/b+d;
        for j=3:maxj
            d=d*y*(a+j-1)/(j*(b+j-1));  % d = A(j)-A(j-1) = (A(j-1)-A(j-2))*r(j)*z from [4]
            g=g+d;                      % g = A(j) from [4]
            p=abs(d)<tol*abs(g);
            if q && p
                break
            end
            q=p;
        end
    elseif y>0 % big positive y
        alg=2;
        r=(b-a)*(1-a);
        d=r*(b-a+1)*(1-0.5*a)/y^2;
        g=1+r/y+d;
        for j=3:maxj
            d=d*(b-a+j-1)*(j-a)/(j*y);
            g=g+d;
            p=abs(d)<tol*abs(g);
            if q && p
                break
            end
            q=p;
        end
        g=exp(y+gammaln(b)-gammaln(a)+(a-b)*log(y))*g;
    else   % big negative y
        alg=3;
        r=a*(a-b+1);
        d=r*(a+1)*(a-b+2)/(2*z^2);
        g=1-r/z+d;
        for j=3:maxj
            r=(a+j-1)*(a-b+j)/j;
            d=-d*r/z;
            g=g+d;
            p=abs(d)<tol*abs(g);
            if q && p
                break
            end
            q=p;
        end
        g=exp(y+gammaln(b)-gammaln(b-a)-a*log(-y))*g;
    end
h(i)=g;
l(i)=j;
end