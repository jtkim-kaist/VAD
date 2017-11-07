function [eta,etaf]=finishat(frac,tol,fmt)
%FINISHAT print estimated finish time of a long computation (FRAC,TOL,FMT)
% Usage: (1)  finishat(0);
%              for i=1:many
%                  ... computation ...
%                  finishat(i/many);
%              end
%
%         (2)  finishat(0);
%              for i=1:NI
%                  for j=1:NJ
%                      for k=1:NK
%                          ... computation ...
%                          finishat((((i-1)*NJ+j-1)*NK+k)/(NI*NJ*NK));
%                      end
%                  end
%              end
%          
% Inputs: FRAC = fraction of total comutation that has been completed
%                As a special case, 0 is used to initialize the routine
%         TOL  = Tolerance in minutes. If the estimated time has changed by less
%                than this, then nothing will be printed. [default 10% of remaining time]
%         FMT  = Format string which should include %s for estimated finish time, %d for remaining minutes and %f for fraction complete
%
% Output: ETA  = string containing the expected finish time
%                specifying this will suppress printing message to std err (fid=2)
%         ETAF = expected finish time as a daynumber
%
% Example:       finishat(0);
%                for i=1:many
%                    long computation;
%                    finishat(i/many);
%                end

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: finishat.m 8893 2016-10-27 16:52:30Z dmb $
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

persistent oldt oldnw
if nargin<3
    fmt='Estimated finish at %s (%.2f done, %d min remaining)\n';

end
if frac<=0
    oldt=0;
    eta='Unknown';
    tic;
else
    nw=now;                             % current time as serial number
    sectogo=(1/frac-1)*toc;      % seconds to go
    newt=nw+sectogo/86400;       % add estimated time in days
    if nargin<2 || ~numel(tol)
        tol=max(0.1*(newt-nw)*1440,1);
    end
    if ~exist('oldt','var') || oldt==0 || (abs(newt-oldt)>tol/1440 && (nw-oldnw)>10/86400) || (nw-oldnw)>10/1440 || nargout>0
        oldt=newt;
        if floor(oldt)==floor(nw)
            df='HH:MM';
        else
            df='HH:MM dd-mmm-yyyy';
        end
        eta=datestr(oldt,df);
        if ~nargout
            ix=find(fmt=='%',1);
            while ~isempty(ix)
                fprintf(2,fmt(1:ix-1));
                fmt=fmt(ix:end);
                ix=find(fmt>='a' & fmt<='z',1); % find letter
                switch fmt(ix)
                    case 's'
                        fprintf(2,fmt(1:ix),eta);
                    case 'd'
                        fprintf(2,fmt(1:ix),round(sectogo/60));
                    case 'f'
                        fprintf(2,fmt(1:ix),frac);
                end
                fmt=fmt(ix+1:end);
                ix=find(fmt=='%',1);
            end
            fprintf(2,fmt);
        end
        oldnw=nw;                           %
    end
end
etaf=oldt;