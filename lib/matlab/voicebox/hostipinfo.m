function [ipinfo,tt]=hostipinfo(m,t)
%HOSTIPINFO get host name and internet connection information
%  Usage:    hostname   = hostipinfo;        % computer name
%            ip_address = hostipinfo('i');   % IP address
%
%  Inputs:   m  mode input selects the information to output (see below)
%            t  optional output from an "ipconfig -all" command
%
%  Outputs:  ipinfo   cell array with one row per IP connection. The number
%                     of columns is detemined by the m input (see below).
%                     If only a single output value is selected by the m input, then
%                     ipinfo will be a string rather than a cell array.
%            tt       output from the "ipconfig -all" command
%
%  The "m" input selects whch IP connections to supply information about
%  (rows of IPINFO) and what information to supply (columns of IPINFO).
%  
%  Row selection: A   include all connections
%                 L   Local Internet Connection [default if m non-empty]
%                 B   Bluetooth Connection
%                 P   PPP connection (e.g. a VPn connection)
%
%  Col selection: h   hostname (same for all connections) [default if m is empty]
%                 c   Connection type code: L=local, B=bluetooth, P=PPP or ?=other
%                 t   Connection type (e.g. "Ethernet" or "PPP")
%                 n   Connection name (e.g. "Local Area Connection")
%                 m   MAC address (e.g. "00-53-45-00-00-00"
%                 i   IP address (e.g. "192.168.1.1") [default if m non-empty]

%      Copyright (C) Mike Brookes 2012
%      Version: $Id: hostipinfo.m 1443 2012-02-22 09:10:03Z dmb $
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

if nargin<2
    t=evalc('!ipconfig -all');
    if nargin<1
        m='';
    end
end
tt=t;
if ~numel(m)
    m='Lh';
end
host=regexp(t,'Host Name[ .]*:\s*([\w-]*)','tokens','once'); % get host name
if ~numel(m)
    ipinfo=host{1};
else
    [iad,ad]=regexp(t,'\n(\w*) adapter ([ \w]*):','start','tokens'); % get adapters
    [iadx,adx]=regexp(t,'\n(\w* adapter [ \w]*):','start','tokens'); % get entire adapter names
    [imac,mac]=regexp(t,'Physical Address[ .]*:\s*(\S*)','start','tokens'); % get MAC addresses
    [iip,ip]=regexp(t,'IP Address[ .]*:\s*(\S*)','start','tokens'); % get IP addresses
    nad=numel(iad);  % number of adapters listed in ipconfig
    iad(nad+1)=length(t); % add an extra fake adapter at the end of the list
    fty='hctnmi'; % possible output fields
    nfty=numel(fty);    % number of fields
    ipinfo=cell(nad,nfty); % cell array to store the fields for each adapter
    nr=nad;
    aty='LBP';   % types of adapter
    naty=numel(aty);
    assx={'Ethernet adapter Local Area Connection'; 'Ethernet adapter Bluetooth Network Connection' ; 'PPP adapter '};
    for i=1:nad
        ipinfo{i,1}=host{1};
        ipinfo(i,3:4)=ad{i};
        ix=find((imac>iad(i)) & (imac<iad(i+1)),1);
        if numel(ix)
            ipinfo{i,5}=mac{ix}{1};
        else
            ipinfo{i,5}='';
        end
        ix=find((iip>iad(i)) & (iip<iad(i+1)),1);
        if numel(ix)
            ipinfo{i,6}=ip{ix}{1};
        else
            ipinfo{i,6}='';
        end
        jx=cell(length(assx),1);
        for j=1:length(assx)
            jx{j}=strfind(adx{i}{1},assx{j});    % compare adapter strings with adapter types
        end
        jx=find(~cellfun('isempty',jx));
        if numel(jx)
            ipinfo{i,2}=aty(jx);
        else
            ipinfo{i,2}='?';
        end
    end
    % Select the requested adapters from the rows of ipinfo{}
    if ~any(m=='A')
        ix=(1:naty)*(aty(ones(numel(m),1),:)'==m(ones(naty,1),:));
        ix=ix(ix>0);
        nr=numel(ix);
        if ~nr
            ix=1;  % default to 'L' = Local Internet
        end
        jx=ix;      % jx will hold first match
        for i=1:numel(ix)
            a=strmatch(aty(ix(i)),ipinfo(:,2));
            if numel(a)
                jx(i)=a(1);
            else
                jx(i)=0;
            end
        end
        jx=jx(jx>0);
        ipinfo=ipinfo(jx,:);
    end
    % Select the requested fields from the columns of ipinfo{}
    ix=(1:nfty)*(fty(ones(numel(m),1),:)'==m(ones(nfty,1),:));
    ix=ix(ix>0);
    nc=numel(ix);
    if ~nc
        ix=6;  % default to 'i' = IP address
    end
    ipinfo=ipinfo(:,ix);
    if nr<2 && nc<2
        ipinfo=ipinfo{1};  % just output a string if only one element requested
    end
end
