function [vlist,edgeq,flist,info]=upolyhedron(w,md)
% UPOLYHEDRON calculate uniform polyhedron characteristics
%
% Inputs:  W    Specifies the desired polyhedron in one of three forms:
%                 (a) name e.g. W='cube'; precede by 'dual' for dual or
%                 'laevo' or 'dextro' [default] for chiral polyhedra
%                 (b) index in the list given below e.g. W=6 is the cube; negative for dual
%                     n.{1,2,3,4,5} gives n-sided prism, antiprism, grammic prism, grammic antiprism or grammic crossed antiprism
%                 (c) Wythoff symbol (see below) e.g. W=[3 0 2 4] is the cube
%                          |p q r, p|q r, p q|r or p q r| respectively with 0 for |
%                         using -1 instead of 0 gives the dual of the polyhedron
%                         using -2 (or -3 for dual) gives a reflected version of a snub polyhedron
%          MD   specifies a mode string
%                 'g' plot image
%                 'w' plot wireframe
%                 'f' plot coloured faces
%                 'v' number vertices
%                 't' plot vertex figure (a slice centred on a vertex)
%
%                 ? homogeneous output coordinates
%                 ? create net
%                 ? segment faces to remove internal portions of the surfaces
%                 ? size: [max] diameter=1, [longest] edge=1, include anisotropic normalization
%                   to minimize variance of vertex radius or edge vector length
%                 ? orientation: vertex at the top, largest stable base at bottom
%
% Outputs:
%          VLIST(:,7)  gives the [x y z d n e t] for each vertex
%                       x,y,z = position, d=distance from origin, n=valency, e=edge index, t=type (-ve for reflected)
%          EDGEQ(:,9) has one row for each direction of each edge:
%                         1  v1   first vertex (normally start)
%                         2  v2   second vertex
%                         3  f1   first face (normally on left)
%                         4  f2   second face
%                         5  ev1  next edge around vertex 1 (normally anticlockwise)
%                         6  ef1  next edge around f1 (normally anticlockwise)
%                         7  er   reverse edge
%                         8  z    twisted edge: clockwise neighbours around v1 and v2 are on the same face
%                         9  sf   swap face order: ???
%                        10  sv   swap vertex order: v2 preceeds v1 around f1
%          FLIST(:,7)  gives the [x y z d n e t] for each face
%                       x,y,z = unit normal, d=distance from origin, n=valency, e=edge index, t=type (-ve for reflected)
%          INFO        structure containing the following fields:
%                         hemi      true if faces are hemispherical (i.e. pass through the origin)
%                         onesided  true is one-sided (like a moebius strip)
%                         snub      true if a snub polyhedron

% This software is based closely on a Mathmatica program described in [1] which, in turn, was based on a
% C program described in [2].
%
% Wythoff Symbol
%    p,q,r define a spherical triangle whose angles at the corners are pi/p, pi/q and pi/r;
%    this triangle tiles the sphere if repeatedly reflected in its sides. The
%    polyhedron vertices are at the reflections of a seed vertex as follows
%    where the Vertex configuration gives the polygon orders of the faces around each vertex:
%       |p q r  : Vertex at a point such that when rotated around any of p,q,r by twice the angle at that
%                 corner is displaced by the same distance for each corner. This is a snub polyhedron and
%                 only even numbers of reflections are used to generate vertices. Configuration {3 p 3 q 3 r}
%        p|q r  : Vertex at p. Configuration = {q r q r ... q r} with 2n terms where n is the numerator of p
%        p q|r  : Vertex on pq and on the bisector of angle r. Configuration {p 2r q 2r}
%        p q r| : Vertex at incentre: meeting point of all the angle bisectors. Configuration {2p 2q 2r}
%        |3/2 5/3 3 5/2  This special case is the great dirhombicosidodecahedron. It is a bit weird because
%                 many of the edges are shared by four faces instead of the usual two.
%    If two of p,q,r = 2 then the third is arbitrary (prisms and antiprisms), otherwise only the numerators
%    1:5 can occur, and 4 and 5 cannot occur together. If all are integers then the poyhedron is convex.
%
% References:
%   [1] R. E. Maeder. Uniform polyhedra. The Mathematica Journal, 3 (4): 48–57, 1993.
%   [2] Z. Har’El. Uniform solution for uniform polyhedra. Geometriae Dedicata, 47: 57–110, 1993.
%   [3] H. S. M. Coxeter, M. S. Longuet-Higgins, and J. C. P. Miller. Uniform polyhedra.
%       Philosophical Transactions of the Royal Society A, 246 (916): 401–450, May 1954.
%   [4] P. W. Messer. Closed-form expressions for uniform polyhedra and their duals.
%       Discrete and Computational Geometry, 27 (3): 353–375, Jan. 2002.

%%%% BUGS and SUGGESTIONS %%%%%%
% (1) we should ensure the "first" edges of the vertices and faces are consistent
% (2) need to sort faces and vertices into a type order
% (3) should ensure that for a non-chiral polyhedron, the vertex polarity alternates
% (4) w=75 does not work
% (5) dual of henispherical poyhedron
% (6) flist not calculated correctly for duals
% (7) add additional stuff into info.*
% (8) vertex figures seem to include additional (duplicated) lines
% (9) sort out when reflected vertices are really rotationally congruent
% (10) could optionally colour by face type
% (11) order alphabetically by noun
% (12) allow abbreviated names + prisms with preceding decimal number
% (13) calculate correct names
% (14) include names for duals
% (15) make "names" etc persistent

% Example slide show:
% for i=1:74, disp(num2str(i)); upolyhedron(i); pause(2); end
% for i=5:10, disp(num2str(i)); for j=1:5, upolyhedron(i+j/10); pause(2); end, end

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: upolyhedron.m 713 2011-10-16 14:45:43Z dmb $
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

% [line nn]  referes to the kaleido Mathematica program on which this is based
% see http://www.mathconsult.ch/showroom/unipoly/unipoly.html#Images

% Variables used
%   adjacent    adjacency matrix
%   chi         characteristic
%   cosa        cos of the angle subtended by a half edge
%   d           density
%   e           number of edges
%   even        number of even faces to remove
%   f           number of faces
%   fi(n)       number of faces of type i
%   g           order of group
%   gamma(n)    included spherical triangle angle between centre of face,
%               vertex and edge for face type i
%   hemiQ       =1 for hemispherical faces (includes polyhedron centre)
%   incidence   faces incident at vertex
%   k           type of group (2=dihedral, 3=tetrahedral, 4=octahedral, 5=icosahedral)
%   m           number of edges meeting at each vertex (valence)
%   mi(n)       Number of faces of type i
%   n           number of face types
%   ni(n)       number of edges on i'th face type
%   p           1st Wythoff number
%   q           2nd Wythoff number
%   r           3rd Wythoff number
%   rot(m)      vertex configuration: anti-clockwise list of face types
%               around a vertex
%   snub(m)     =1 for snub triangles
%   snubQ       =1 for snub polyhedron
%   v           number of vertices
%   vcoord      vertex coordinates
%   wy(4)       Wythoff Symbol
%   wyd(4)      Wythoff Symbol Denominators
%   wyn(4)      Wythoff Symbol Numerators

% names: [name abbreviation synonym dual dual-synonym]
persistent names wys prefs
if ~numel(names)
    prefs={'dual '; 'dextro '; 'laevo '};
    names={
        'Tetrahedron' 'Tet' '' '' '';                              %  1
        'truncated tetrahedron' 'Tut' '' '' '';                    %  2
        'octahemioctahedron' 'Oho' '' '' '';                       %  3
        'tetrahemihexahedron' 'Thah' '' '' '';                      %  4
        'Octahedron' 'Oct' '' 'Cube' '';                               %  5
        'cube' 'Cube' '' 'Octahedron' '';                                     %  6
        'cuboctahedron' 'Co' '' 'rhombic dodecahedron' '';                            %  7
        'truncated octahedron' 'Toe' '' 'tetrakishexahedron' '';                     %  8
        'truncated cube' 'Tic' '' 'triakisoctahedron' '';                           %  9
        'rhombicuboctahedron' 'Sirco' '' 'deltoidal icositetrahedron' '';                      % 10
        'truncated cuboctahedron' '' '' 'disdyakisdodecahedron' '';                  % 11
        'snub cube' 'Snic' '' 'pentagonal icositetrahedron' '';                                % 12
        'small cubicuboctahedron' 'Socco' '' 'small hexacronic icositetrahedron' '';                  % 13
        'great cubicuboctahedron' 'Gocco' '' 'great hexacronic icositetrahedron' '';                  % 14
        'cubohemioctahedron' 'Cho' '' 'hexahemioctacron' '';                       % 15
        'cubitruncated cuboctahedron' 'Cotco' '' 'tetradyakishexahedron' '';              % 16
        'great rhombicuboctahedron' 'Girco' '' 'great deltoidal icositetrahedron' '';                % 17 or Querco
        'small rhombihexahedron' 'Sroh' '' 'small rhombihexacron' '';                   % 18
        'stellated truncated hexahedron' 'Quith' '' 'great triakisoctahedron' '';          	% 19
        'great truncated cuboctahedron' 'Quitco' '' 'great disdyakisdodecahedron' '';         	  % 20
        'great rhombihexahedron' 'Groh' '' 'great rhombihexacron' '';                   % 21
        'icosahedron' 'Ike' '' 'dodecahedron' '';                              % 22
        'dodecahedron' 'Doe' '' 'icosahedron' '';                             % 23
        'icosidodecahedron' 'Id' '' 'rhombic triacontahedron' '';                        % 24
        'truncated icosahedron' 'Ti' '' 'pentakisdodecahedron' '';                    % 25
        'truncated dodecahedron' 'Tid' '' 'triakisicosahedron' '';                   % 26
        'rhombicosidodecahedron' 'Srid' '' 'deltoidal hexecontahedron' '';                   % 27
        'truncated icosidodechedon' 'Grid' 'great rhombiicosidodecahedron' 'disdyakistriacontahedron' '';                % 28
        'snub dodecahedron' 'Snid' '' 'pentagonal hexecontahedron' '';                        % 29
        'small ditrigonal icosidodecahedron' 'Sidtid' '' 'small triambic icosahedron' '';     	% 30
        'small icosicosidodecahedron' 'Siid' '' 'small icosacronic hexecontahedron' '';              % 31
        'small snub icosicosidodecahedron' 'Seside' '' 'small hexagonal hexecontahedron' '';        	% 32
        'small dodecicosidodecahedron' 'Saddid' '' 'small dodecacronic hexecontahedron' '';            	% 33
        'small stellated dodecahedron' 'Sissid' '' 'great dodecahedron' '';            	% 34
        'great dodecahedron' 'Gad' '' 'small stellated dodecahedron' '';                       % 35
        'dodecadodecahedron' 'Did' '' 'medial rhombic triacontahedron' '';                       % 36
        'truncated great dodecahedron' 'Tigid' '' 'small stellapentakisdodecahedron' '';            	% 37
        'rhombidodecadodecahedron' 'Raded' '' 'medial deltoidal hexecontahedron' '';                 % 38
        'small rhombidodecahedron' 'Sird' '' 'small rhombidodecacron' '';                 % 39
        'snub dodecadodecahedron' 'Siddid' '' 'medial pentagonal hexecontahedron' '';                  % 40
        'ditrigonal dodecadodecahedron' 'Ditdid' '' 'medial triambic icosahedron' '';          	% 41
        'great ditrigonal dodecicosidodecahedron' 'Gidditdid' '' 'great ditrigonal dodecacronic hexecontahedron' '';	% 42
        'small ditrigonal dodecicosidodecahedron' 'Sidditdid' '' 'small ditrigonal dodecacronic hexecontahedron' '';	% 43
        'icosidodecadodecahedron' 'Ided' '' 'medial icosacronic hexecontahedron' '';                  % 44
        'icositruncated dodecadodecahedron' 'Idtid' '' 'tridyakisicosahedron' '';        % 45
        'snub icosidodecadodecahedron' 'Sided' '' 'medial hexagonal hexecontahedron' '';             % 46
        'great ditrigonal icosidodecahedron' 'Gidtid' '' 'great triambic icosahedron' '';       % 47
        'great icosicosidodecahedron' 'Giid' '' 'great icosacronic hexecontahedron' '';              % 48
        'small icosihemidodecahedron' 'Seihid' '' 'small icosihemidodecacron' '';              % 49
        'small dodecicosahedron' 'Siddy' '' 'small dodecicosacron' '';                   % 50
        'small dodecahemidodecahedron' 'Sidhid' '' 'small dodecahemidodecacron' '';             % 51
        'great stellated dodecahedron' 'Gissid' '' 'great icosahedron' '';           	% 52
        'great icosahedron' 'Gike' '' 'great stellated dodecahedron' '';                        % 53
        'great icosidodecahedron' 'Gid' '' 'great rhombic triacontahedron' '';                  % 54
        'great truncated icosahedron' 'Tiggy' '' 'great stellapentakisdodecahedron' '';              % 55
        'rhombicosahedron' 'Ri' '' 'rhombicosacron' '';                         % 56
        'great snub icosidodecahedron' 'Gosid' '' 'great pentagonal hexecontahedron' '';             % 57
        'small stellated truncated dodecahedron' 'Quitsissid' '' 'great pentakisdodekahedron' ''; 	% 58
        'truncated dodecadodecahedron' 'Quitdid' '' 'medial disdyakistriacontahedron' '';           	% 59
        'inverted snub dodecadodecahedron' 'Isdid' '' 'medial inverted pentagonal hexecontahedron' '';       	% 60
        'great dodecicosidodecahedron' 'Gaddid' '' 'great dodecacronic hexecontahedron' '';           	% 61
        'small dodecahemicosahedron' 'Sidhei' '' 'small dodecahemicosacron' '';               % 62
        'great dodecicosahedron' 'Giddy' '' 'great dodecicosacron' '';                   % 63
        'great snub dodecicosidodecahedron' 'Gisdid' '' 'great hexagonal hexecontahedron' '';      	% 64
        'great dodecahemicosahedron' 'Gidhei' '' 'great dodecahemicosacron' '';               % 65
        'great stellated truncated dodecahedron' 'Quitgissid' '' 'great triakisicosahedron' '';  	% 66
        'great rhombicosidodecahedron' 'Qrid' '' 'great deltoidal hexecontahedron' '';            	% 67 or Qrid
        'great truncated icosidodecahedron' 'Gaquatid' '' 'great disdyakistriacontahedron' '';      	% 68
        'great inverted snub icosidodecahedron' 'Gisid' '' 'great inverted pentagonal hexecontahedron' '';   	% 69
        'great dodecahemidodecahedron' 'Gidhid' '' 'great dodecahemidodecacron' '';          	  % 70
        'great icosihemidodecahedron' 'Geihid' '' 'great icosihemidodecacron' '';              % 71
        'small retrosnub icosicosidodecahedron' 'Sirsid' '' 'small hexagrammic hexecontahedron' '';   	% 72
        'great rhombidodecahedron' 'Gird' '' 'great rhombidodecacron' '';                 % 73
        'great retrosnub icosidodecahedron' 'Girsid' '' 'great pentagrammic hexecontahedron' '';        % 74
        'great dirhombicosidodecahedron' 'Gidrid' '' 'great dirhombicosidodecacron' '';         	% 75
        'pentagonal prism' 'Pip' '' '' '';                         % 76
        'pentagonal antiprism' 'Pap' '' '' '';                     % 77
        'pentagrammic prism' 'Stip' '' '' '';                       % 78
        'pentagrammic antiprism' 'Stap' '' '' '';                   % 79
        'pentagrammic crossed antiprism'  'Starp'   '' '' '';      	% 80
        'triangular prism' '' '' '' '';      % 81
        'triangular antiprism' '' '' '' '';     % 82
        'square prism' '' '' '' '';     % 83
        'square antiprism' '' '' '' ''};    % 84

    % could multiply by 12 to make integer
    wys=[
        2, 3, 2, 3;        %  1   'tetrahedron' 'Tet';
        3, 2, 3, 3;        %  2   'truncated tetrahedron' 'Tut';
        3, 3/2, 3, 3;      %  3   'octahemioctahedron' 'Oho';
        3, 3/2, 3, 2;      %  4   'tetrahemihexahedron' 'Thah';
        2, 4, 2, 3;        %  5   'octahedron' 'Oct';
        2, 3, 2, 4;        %  6   'cube' 'Cube';
        2, 2, 3, 4;        %  7   'cuboctahedron' 'Co';
        3, 2, 4, 3;        %  8   'truncated octahedron' 'Toe';
        3, 2, 3, 4;        %  9   'truncated cube' 'Tic';
        3, 3, 4, 2;        % 10   'rhombicuboctahedron' 'Sirco';
        4, 2, 3, 4;        % 11   'truncated cuboctahedron' '?';
        1, 2, 3, 4;        % 12   'snub cube' 'Snic';
        3, 3/2, 4, 4;      % 13   'small cubicuboctahedron' 'Socco';
        3, 3, 4, 4/3;      % 14   'great cubicuboctahedron' 'Gocco';
        3, 4/3, 4, 3;      % 15   'cubohemioctahedron' 'Cho';
        4, 4/3, 3, 4;      % 16   'cubitruncated cuboctahedron' 'Cotco';
        3, 3/2, 4, 2;      % 17   'great rhombicuboctahedron' 'Girco?';
        4, 3/2, 2, 4;      % 18   'small rhombihexahedron' 'Sroh';
        3, 2, 3, 4/3;      % 19   'stellated truncated hexahedron' 'Quith';
        4, 4/3, 2, 3;      % 20   'great truncated cuboctahedron' 'Quitco';
        4, 4/3, 3/2, 2;    % 21   'great rhombihexahedron' 'Groh';
        2, 5, 2, 3;        % 22   'icosahedron' 'Ike';
        2, 3, 2, 5;        % 23   'dodecahedron' 'Doe';
        2, 2, 3, 5;        % 24   'icosidodecahedron' 'Id';
        3, 2, 5, 3;        % 25   'truncated icosahedron' 'Ti';
        3, 2, 3, 5;        % 26   'truncated dodecahedron' 'Tid';
        3, 3, 5, 2;        % 27   'rhombicosidodecahedron' 'Srid';
        4, 2, 3, 5;        % 28   'truncated icosidodecahedron' 'Grid'; or 'great rhombiicosidodecahedron' Grid
        1, 2, 3, 5;        % 29   'snub dodecahedron' 'Snid';
        2, 3, 5/2, 3;      % 30   'small ditrigonal icosidodecahedron' 'Sidtid';
        3, 5/2, 3, 3;      % 31   'small icosicosidodecahedron' 'Siid';
        1, 5/2, 3, 3;      % 32   'small snub icosicosidodecahedron' 'Seside';
        3, 3/2, 5, 5;      % 33   'small dodecicosidodecahedron' 'Saddid';
        2, 5, 2, 5/2;      % 34   'small stellated dodecahedron' 'Sissid';
        2, 5/2, 2, 5;      % 35   'great dodecahedron' 'Gad';
        2, 2, 5/2, 5;      % 36   'dodecadodecahedron' 'Did';
        3, 2, 5/2, 5;      % 37   'truncated great dodecahedron' 'Tigid';
        3, 5/2, 5, 2;      % 38   'rhombidodecadodecahedron' 'Raded';
        4, 2, 5/2, 5;      % 39   'small rhombidodecahedron' 'Sird';
        1, 2, 5/2, 5;      % 40   'snub dodecadodecahedron' 'Siddid';
        2, 3, 5/3, 5;      % 41   'ditrigonal dodecadodecahedron' 'Ditdid';
        3, 3, 5, 5/3;      % 42   'great ditrigonal dodecicosidodecahedron' 'Gidditdid'
        3, 5/3, 3, 5;      % 43   'small ditrigonal dodecicosidodecahedron' 'Sidditdid'
        3, 5/3, 5, 3;      % 44   'icosidodecadodecahedron' 'Ided';
        4, 5/3, 3, 5;      % 45   'icositruncated dodecadodecahedron' 'Idtid';
        1, 5/3, 3, 5;      % 46   'snub icosidodecadodecahedron' 'Sided';
        2, 3/2, 3, 5;      % 47   'great ditrigonal icosidodecahedron' 'Gidtid';
        3, 3/2, 5, 3;      % 48   'great icosicosidodecahedron' 'Giid';
        3, 3/2, 3, 5;      % 49   'small icosihemidodecahedron' 'Seihid';
        4, 3/2, 3, 5;      % 50   'small dodecicosahedron' 'Siddy';
        3, 5/4, 5, 5;      % 51   'small dodecahemidodecahedron' 'Sidhid';
        2, 3, 2, 5/2;      % 52   'great stellated dodecahedron' 'Gissid';
        2, 5/2, 2, 3;      % 53   'great icosahedron' 'Gike';
        2, 2, 5/2, 3;      % 54   'great icosidodecahedron' 'Gid';
        3, 2, 5/2, 3;      % 55   'great truncated icosahedron' 'Tiggy';
        4, 2, 5/2, 3;      % 56   'rhombicosahedron' 'Ri';
        1, 2, 5/2, 3;      % 57   'great snub icosidodecahedron' 'Gosid';
        3, 2, 5, 5/3;      % 58   'small stellated truncated dodecahedron' 'Quitsissid'
        4, 5/3, 2, 5;      % 59   'truncated dodecadodecahedron' 'Quitdid';
        1, 5/3, 2, 5;      % 60   'inverted snub dodecadodecahedron' 'Isdid';
        3, 5/2, 3, 5/3;    % 61   'great dodecicosidodecahedron' 'Gaddid';
        3, 5/3, 5/2, 3;    % 62   'small dodecahemicosahedron' 'Sidhei';
        4, 5/3, 5/2, 3;    % 63   'great dodecicosahedron' 'Giddy';
        1, 5/3, 5/2, 3;    % 64   'great snub dodecicosidodecahedron' 'Gisdid';
        3, 5/4, 5, 3;      % 65   'great dodecahemicosahedron' 'Gidhei';
        3, 2, 3, 5/3;      % 66   'great stellated truncated dodecahedron' 'Quitgissid'
        3, 5/3, 3, 2;      % 67   'quasirhombicosidodecahedron' 'Qrid'; 'great rhombicosidodecahedron' 'Nonconvex great rhombicosidodecahedron
        4, 5/3, 2, 3;      % 68   'great truncated icosidodecahedron' 'Gaquatid';
        1, 5/3, 2, 3;      % 69   'great inverted snub icosidodecahedron' 'Gisid';
        3, 5/3, 5/2, 5/3;  % 70   'great dodecahemidodecahedron' 'Gidhid';
        3, 3/2, 3, 5/3;    % 71   'great icosihemidodecahedron' 'Geihid';
        1, 3/2, 3/2, 5/2;  % 72   'small retrosnub icosicosidodecahedron' 'Sirsid';
        4, 3/2, 5/3, 2;    % 73   'great rhombidodecahedron' 'Gird';
        1, 3/2, 5/3, 2;    % 74   'great retrosnub icosidodecahedron' 'Girsid';
        5, 3/2, 5/3, 3;    % 75   'great dirhombicosidodecahedron' 'Gidrid';
        3, 2, 5, 2;        % 76   'pentagonal prism' 'Pip';
        1, 2, 2, 5;        % 77   'pentagonal antiprism' 'Pap';
        3, 2, 5/2, 2;      % 78   'pentagrammic prism' 'Stip';
        1, 2, 2, 5/2;      % 79   'pentagrammic antiprism' 'Stap';
        1, 2, 2, 5/3;     % 80   'pentagrammic crossed antiprism'  'Starp'
        3, 2, 3, 2;        % 81   'triangular prism' ;
        1, 2, 2, 3;        % 82   'triangular antiprism' ;
        3, 2, 4, 2;        % 83   'square prism' ;
        1, 2, 2, 4];        % 84   'square antiprism' ;
end

if nargin<2
    md='';
end
dual=0;  % dual polyhedron
dextro=0; % reflected (for snub polyhedra only)
if ischar(w)
    i=0;
    while i<length(prefs)
        i=i+1;
        if length(w)>length(prefs{i}) && strcmpi(w(1:length(prefs{i})),prefs{i})
            switch i
                case 1
                    dual=1;
                case 2
                    dextro=1;
            end
            w=w(length(prefs{i})+1:end);
            i=0;
        end
    end

    wyidx=find(strcmpi(w,names(:)),1);
    if isempty(wyidx)
        % check for prism names here
        error('Cannot find %s',w);
    end
    if wyidx>3*size(names,1)
        dual=1-dual;
    end
    wyidx=1+rem(wyidx-1,size(names,1));
    wy=wys(wyidx,:);
else
    if length(w)==1
        dual=w<0;
        wyidx=abs(w);
        wyjdx=floor(wyidx);
        wykdx=round(10*(wyidx-wyjdx));
        switch wykdx
            case 1  % prism
                wy=[3 2 wyjdx 2];
            case 2  % antiprism
                wy=[1 2 2 wyjdx];
            case 3  % grammic prism
                wy=[3 2 wyjdx/2 2];
            case 4  % grammic antiprism
                wy=[1 2 2 wyjdx/2];
            case 5  % grammic crossed antiprism
                wy=[1 2 2 wyjdx/3];
            otherwise
                if wyjdx>=1 && wyjdx<=size(wys,1)
                    wy=wys(wyjdx,:);
                else
                    error('Polyhedron number out of range');
                end
        end
    elseif length(w)==4
        vbar=find(w<=0,1);
        if ~numel(vbar)
            error('Invalid Wythoff symbol: %g %g %g %g %g',w);
        end
        dual=mod(w(vbar),2);	% least significant bit indicates dual
        dextro=mod(1+w(vbar),4)>=2; % second bit indicates reflected version
        wy=[vbar w(1:vbar-1) w(vbar+1:4)];
        [xx,wyidx]=min(sum((wys-repmat(wy,size(wys,1),1)).^2,2));
        if abs(xx)>1e-3
            if sum(wy==2)<2
                error('Invalid Wythoff symbol: %g %g %g %g %g',w);
            else
                % need to sort out prism names here
            end
        end
    elseif length(w)==5
        if w(1)<=0 && all(round(12*w(2:end))==[18 40 36 30])   % [3/2 5/3 3 5/2]
            dual=w(1)<0;
            wyidx=75;
            wy=wys(wyidx,:);
        else
            error('Invalid Wythoff symbol: %g %g %g %g %g',w);
        end
    else
        error('Invalid polyhedron specification');
    end
end


[wyn,wyd]=rat(wy(2:4));    % convert to rational numbers
wy(2:4)=wyn./wyd;           % force exact rational values
p=wy(2);
q=wy(3);
r=wy(4);
hemiQ = 0;          % includes hemispherical faces that go through polyhedron centre [ p q | r] with pq=p+q
onesidedQ = 0;      % one-sided polyhedron (aka: non-orientable) [ p q r |] with exactly one of p,q,r having an even denominator
evenQ = 0;           % identifies which of p, q, r has an even denominator
even=1;
snubQ = 0;          % snub polyhedron: [ | p q r ]
allrot = 1;         % all vertices a congruent with rotations (no reflections needed)

% call to AnalyseWythoff[line 105]

k= max(wyn);
if sum(wy(2:4)==2)>=2      % check if it is a prism
    % this is a prism - treat specially [note mathematica line 86 has redundant chack for 2 >5]
    g=4*k;      % order of group
    k=2;
else                       % not a prism: only numerators 1,2,3,4,5 are allowed
    if (any(wy(2:4)<=1) || (k>5) || any(wyn==4) && any(wyn==5))
        error('Invalid Wythoff numbers [%d/%d %d/%d %d/%d]',[wyn;wyd]);
    end
    g=24*k/(6-k);       % order of group
end
if abs(wy(1))==5 % special case
    chi = 0;
    d = 0;
else
    chi = (sum(wyn.^(-1))-1)*g/2;
    d = (sum(wy(2:4).^(-1))-1)*g/4;
    if d<0
        error('density < 0');
    end
end
if abs(wy(1))==5
    error('great dirhombicosidodecahedron not implemented');
end
switch abs(wy(1))
    case {1,5}    % [ | p q r ] snub polyhedron
        n=4;            % number of face types
        m=6;            % valence of vertices
        v=g/2;          % number of vertices
%         ni=[3 wy(2:4)]; % type of each face
%         mi=[3 1 1 1];   % number of faces of each type
        nimi=[3 wy(2:4); 3 1 1 1]';
%         rot = [1 2+dextro 1 3-dextro 1 4];
        snubQ=1;
%         snub=[1 0 1 0 1 0];
        rotsnub=[1 2+dextro 1 3-dextro 1 4; repmat([1 0],1,3)]';
    case 2  % [ p | q r ]
        n=2;
        m=2*wyn(1);
        v=g/m;
%         ni=wy(3:4);
%         mi=[p p];         % these face counts may be fractional: m counts their numerators
% we could ill in the numerator denominator columns, nimi(:,3:4) here instead of later
        nimi=[wy(3:4); p p]';
%         rot = repmat(1:2,1,wyn(1));
        rotsnub = repmat([1 0; 2 0],wyn(1),1);
    case 3  % [ p q | r ]
        n=3;
        m=4;
        v=g/2;
%         ni=[2*r wy(2:3)];
%         mi = [2 1 1];
        nimi=[2*r wy(2:3); 2 1 1]';
%         rot = [1 2 1 3];
        rotsnub = [1 2 1 3; zeros(1,4)]';
        if abs(p-q/(q-1))<1e-6
            hemiQ=1;
            d=0;
            if (p~=2 && ~(wy(4)==3 && any(wy(2:3)==3)))
                onesidedQ=1;
                v=v/2;
                chi=chi/2;
            end
        end
    case 4  % [ p q r | ]
        n=3;
        m=3;
        v=g;
%         ni=2*wy(2:4);
%         mi=ones(1,3);
        nimi=[2*wy(2:4); ones(1,3)]';
%         rot=1:3;
        rotsnub=[1:3; zeros(1,3)]';
        allrot=(wy(2)==wy(3)) || (wy(2)==wy(4)) || (wy(3)==wy(4)); % isosceles so all
        evenden=find(~rem(wyd,2)); % check for even denominators
        if (length(evenden)>1)
            error('Multiple even denominators are not allowed');
        end
        if ~isempty(evenden)
            even = evenden;
            evenQ=1;
            onesidedQ = 1;
            d=0;
            v=v/2;
            chi = chi - g/wyn(even)/2;  % check this *****
        end
    otherwise
        error('Invalid polyhedron type %d',wy(1));
end
% [line 155] call sortAndMerge [sortAndMerge: line 220]
% save stuff to prevent overwriting
n_1=n;
m_1=m;
% now do sort merge based on rotsnub and nimi
% rotsnub(:,2) = [rot snub]
% nimi(:,4) = [ni mi mi-num mi-den]

nimi(nimi(:,1)==2,1)=-2;            % make equal to -2 to force to bottom of list
[nimi,inim]=sortrows(nimi,-1);
jnim=zeros(n,1);
jnim(inim)=1:n;
msk=[~0; nimi(2:end,1)<nimi(1:end-1,1)]; % identify distinct values
cmsk=cumsum(msk); % destination of each value
nimi(msk,2)=sparse(1,cmsk,nimi(:,2));  % add up repreated counts
nimi=nimi(msk,:);       % select just the distinct values
% cmsk(jnim) maps original to new positions
[nimi(:,3),nimi(:,4)]=rat(nimi(:,2));
rotsnub(:,1)=cmsk(jnim(rotsnub(:,1)));
even=cmsk(jnim(even));
if nimi(end,1)<0        % remove digons
    if size(nimi,1)==1
        error('Degenerate polyhedron (digons only)');
    end
    m=m-nimi(end,3);  % note that m=sum(nimi(:,3)) not sum(nimi(:,2)) as you might expect
    nimi(end,:)=[];
    n=size(nimi,1);
    rotsnub(rotsnub(:,1)>n,:)=[];         % abolish references to digons
end

% original sort and merge

% [ss,tr]=sort(-ni);          % sort into descending order
% itr=zeros(1,n);
% itr(tr)=1:n;
% ni=ni(tr);
% mi=mi(tr);
% rot = itr(rot);
% if evenQ
%     even = itr(even);
% end
% % now merge equal faces
% i=1;
% while i<n
%     if ni(i)==ni(i+1)
%         mi(i)=mi(i)+mi(i+1);
%         mi(i+1)=[];
%         ni(i+1)=[];
%         n=n-1;
%         even=even-(even>i);
%         rot = rot - (rot>i);
%     else
%         i=i+1;
%     end
% end
% [mii,mij]=rat(mi);    % find numerator
% i=find(ni==2);  % find digons
% if ~isempty(i)
%     if n==1
%         error('Degenerate (digons only)');
%     end
%     i=i(1);
%     m=m-mii(i);      % reduce total valance
%     ni(i)=[];
%     mi(i)=[];
%     mii(i)=[];
%     even=even-(even>i);
%     if snubQ
%         snub(rot==i)=[];
%     end
%     rot(rot==i)=[];
%     rot = rot - (rot>i);
%     n=n-1;
% end
% 
% % check that new version is correct
% 
% 
% if any(nimi(:,1:3)~=[ni' mi' mii'])
%     error('nimi has errors');
% end
% if snubQ
%     if any(rotsnub~=[rot' snub'])
%         error('rotsnub has errors');
%     end
% else
%     if any(rotsnub(:,1)~=rot')
%         error('rotsnub has errors');
%     end
% end
% if evenQ && even~=even_1
%     error('even has errors');
% end

% [line 159] Solve fundamental equations [line 266]
% doesn't converge well for a dodecahedron
% we could avoid the iteration when n=1
%
% We want to solve the following set of equations for cosa and gammai:
%  (a) sum(mi.*gammai)=pi
%  (b) cos(alphai) = cosa*sin(gammai)
%      where ni, mi, alphai and gammai are vectors of length m, the number of edges at a vertex
%      ni gives the number of sides in each face type in decreasing order (possibly fractional)
%      mi is the repretition count of each face type (possibly fractional)
%      alphai=pi * ni.^(-1) is the angle subtended by a half edge at the centre of a face
%      gammai (> pi/2 - alphai) is half the spherical angle subtended by a face at the vertex
%      cosa (< 1) is the distance of the edge centre from the origin if the radius of the vertex is unity
% we solve this iteratively using gammai(1) as the independent variable


alphai = pi*nimi(:,1).^(-1); % half of external angle of polygon
calphai = cos(alphai);
ca1 = calphai(1);
calphai(1)=[];
% initial values
gammai = pi/2 - alphai; % half of internal angle of polygon
delta = pi - nimi(:,2)'*gammai;
% we could initialize it to be exact for a single face type
% g1n=gammai(1)*pi/(mi*gammai')
% cosa = ca1/sin(g1n);
% gammai = [g1n asin(calphai/cosa)];
% delta = pi - mi*gammai';
% should check here to see if delta is better than before and g1n is valid
iter=0;
while (abs(delta) > 1e-10)
    g1n = gammai(1) + delta*tan(gammai(1))/(nimi(:,2)'*tan(gammai));
    if (g1n<0 || g1n>pi)
        error ('Gamma out of range');
    end
    cosa = ca1/sin(g1n);
    gammai = [g1n; asin(calphai/cosa)];
    delta = pi - nimi(:,2)'*gammai;
    iter=iter+1;
end
gamma = gammai;
cosa = ca1/sin(gammai(1));

% [line 165] postprocess special cases

if evenQ
    nimi(even,:)=[];
    gamma(even)=[];
    nh=n-1;
    n=2*nh;
    m=4;
%     ni=[ni 1+fliplr(ni-1).^(-1)];
nimi=repmat(nimi,2,1);
nimi(nh+1:end,1)=1+(nimi(nh:-1:1,1)-1).^(-1);
    gamma=[gamma; -flipud(gamma)];
%     mi=repmat(3-nh,1,n);
%     mii=mi;
nimi(:,2:3)=repmat(3-nh,n,2);
    rotsnub=[1 nh n 1+nh; zeros(1,4)]';
end
if wy(1)==5
    error('not yet implemented');
    % needs to use rotsnub and nimi
    n=5;
    m=8;
    hemiQ = 1;
    d=0;
    mi=[4 1 1 1 1];
    mii=mi;
    ni=[4 ni(1:3) ni(1)/(ni(1)-1)];
    gamma=[pi/2 gamma(1:3) -gamma(1)];
    rot=[rot+[0 1 0 1 0 1] 1 5]; % [1 4 1 3 1 2 1 5]
    snub=[snub 1 0];
end



% [line 197] count vertices and faces [count: line 286]

e = m*v/2;          % number of edges = valency * vertices /2
[nii,nij]=rat(nimi(:,1));  % find numerators
fi = v*nimi(:,3)./nii;
f = sum(fi);
if (d>0) && (gamma(1)> pi/2) % [might be better to use cycles rather than radians for gamma]
    d = fi(1)-d;
end
if wy(1)==5
    chi = v - e + f; % Euler characteristic
end

% [line 199] generate vertex coordinates and adjacency matrix [vertices: line 297]
% we currently scale it so the edge mid-points have unit radius

adj=zeros(v,m);         % j'th vertex adjacent to vertex i anti-clockwise;
frot = zeros(1,v);
rev = zeros(v,1);
vlist = zeros(v+1,7); % vertex information: [x y z d n e t] for each vertex
%                       x,y,z = position, d=distance from origin, n=valency, e=edge index, t=type (-ve for reflected)
vlist(:,5)=m;           % all vertices have the same valency
v1 = [0 0 1]/cosa;        % first vertex
frot(1) = 1;
adj(1,1) = 2;
v2 = [2*cosa*sqrt(1-cosa^2), 0, 2*cosa^2-1]/cosa;   % second vertex
if snubQ
    frot(2) = m*(1-rotsnub(m,2))+rotsnub(m,2); % inefficient
    adj(2,1)=1;
else                    % reflexible
    frot(2) = 1;
    rev(2) = 1;
    adj(2,m)=1;
end
vlist(1,1:3)=v1;
vlist(2,1:3)=v2;
nv = 2;
i = 1;
skew=zeros(3,3);
skewp=[6 7 2];  % index of positive entries
skewn=[8 3 4];  % index of negative entries
cosg=cos(2*gamma);
sing=sin(2*gamma);
eye3=eye(3);
% veqth=cos(acos(cosa)/(d+1));   % threshold for vertex equality test = cosa * |vlist(1,1:3)|^2
veqth=0.999999/(cosa^2);  % threshold for vertex equality test = cosa * |vlist(1,1:3)|^2
while i<=nv  % loop for each vertex
    if rev(i)
        one = -1;
        start = m-1;
        limit = 1;
    else
        one = 1;
        start = 2;
        limit = m;
    end
    k = frot(i);    % rotation to use first
    v1 = vlist(i,1:3);    % the centre of rotation
    v1 = v1/sqrt(v1*v1');   % normalize to unit length [not clear why we don't make unit length in the first place]
    v2 = vlist(adj(i,start-one),1:3);
    for j=start:one:limit
        %             R=wwT+cos(x)(I-wwT)+sin(x)SKEW(w)
        %             SKEW(a) =  [0 -a3 a2; a3 0 -a1; -a2 a1 0]
        skew(skewp)=v1;
        skew(skewn)=-v1;
        rsym=v1'*v1;
        rotk=rotsnub(k,1);
        rotmat=rsym+cosg(rotk)*(eye3-rsym)-one*sing(rotk)*skew;
        v2=v2*rotmat;  % rotate v2 poition around v1
        vlist(nv+1,1:3)=v2;   % add into list in case it is good
        nvi=find(vlist(:,1:3)*v2'>veqth,1); % take the first matching vertex
        %         cmatch=1-vlist(nvi,:)*v2'*cosa^2      % diagnostic printout for vertex match test (0 for perfect match)
        adj(i,j)=nvi; % save as next vertex
        lastk = k;
        k=1+mod(k,m); % increment k circularly in the range 1:m
        if nvi>nv       % we have a new vertex
            if snubQ % if a snub polyhedron
                % Mathematica: frot[[nvi]] = If[ !sn[[lastk]], lastk, If[ !sn[[k]], next[k, m], k ] ];
                % A snub triangle plays the role of the next available snub triangle for the next vertex
                frot(nvi)=lastk+rotsnub(lastk,2)*(rotsnub(k,2)*k-lastk+(1-rotsnub(k,2))*(1+mod(k,m)));
                rev(nvi)=0;   % snub polyhedra always have rev=0
                adj(nvi,1)=i;
            else
                frot(nvi)=k;
                rev(nvi)=(1+one)/2; % = 1-rev(i)
                adj(nvi,1+(m-1)*rev(nvi))=i;
            end
            nv=nvi;
            if nv>v
                error('Too many vertices found');
            end
        end
    end
    i=i+1;
end
if (nv~=v)
    error('Not all vertices found');
end
vlist(v+1,:)=[];    % remove the dummy extra vertex
vlist(:,7)=1-2*rev;     % vertex types all all +1 or -1

% construct edge map
% edgeq: 1=v1 2=v2 3=f1 4=f2 5=ev1 6=ef1 7=er 8=z 9=sf 10=sv]
%    1  v1   first vertex (normally start)
%    2  v2   second vertex
%    3  f1   first face (normally on left)
%    4  f2   second face
%    5  ev1  next edge around vertex 1 (normally anticlockwise)
%    6  ef1  next edge around f1 (normally anticlockwise)
%    7  er   reverse edge
%    8  z    twisted edge: clockwise neighbours around v1 and v2 are on the same face
%    9  sf   swap face order: ???
%   10  sv   swap vertex order: v2 preceeds v1 around f1

edgeq=zeros(v*m,10); % OLD *** ===[reverse_edge start_vertex left_face next_edge_at_vertex next_edge_on_face]
edgeq(:,1:2)=[repmat((1:v)',m,1) adj(:)];  % 2=v_start 3=v_end
edgeq(:,3:4)=edgeq(:,1:2)*[v 1; 1 v];      % encode start and end as a single integer
[xx,ia]=sort(edgeq(:,3));
[xx,ib]=sort(edgeq(:,4));
edgeq(:,3:4)=0;
ib(ib)=1:v*m;               % make reverse index
edgeq(:,7)=ia(ib);          % index to reverse direction edge: 1=e_reverse
edgeq(:,5)=[v+1:v*m 1:v];   % adjacent edge anti-clockwise around vertex: 6=e_next_at_start
ia(edgeq(:,5))=1:v*m;       % adjacent edge clockwise around vertex; alternatively ia = mod((1:v*m)-v,v*m)'
edgeq(:,6)=ia(edgeq(:,7));  % adjacent edge anti-clockwise around face [but not necessarily correct
fpt=zeros(f,2);             % edge count and pointer to an edge on a face
iff=0;                      % face index
if onesidedQ
    edgeq(:,8)=reshape(mod(repmat(frot',1,m)+repmat((1:m),v,1),2),v*m,1); %face type modulo 2
    edgeq(:,8)=abs(edgeq(:,8)-edgeq(edgeq(:,6),8)); % 1 = twisted edge
    while iff<f
        iee=find(edgeq(:,3)==0,1); % find an edge without a left face
        if ~numel(iee)
            error('Not enough faces found');
        end
        iff=iff+1;
        fpt(iff,2)=iee;
        nee=0;
        flip=abs((edgeq(iee,6)==edgeq(edgeq(iee,7),5))-edgeq(iee,8)); % check if this initial edge is flipped on output
        while edgeq(iee,3)~=iff
            if edgeq(iee,3)  % use the reverse edge if this one is already taken
                iee=edgeq(iee,7);
                edgeq(pee,6)=iee;           % correct the previous exit edge taken
                edgeq(iee,10)=1;     % use it in the reverse direction
                edgeq(iee,6)=ia(iee); % normal exit for this direction
                jee=edgeq(iee,5);  % alternate exit edge if flipping
            else
                if flip && ~edgeq(iee,4) % we are entering on the usual reverse exit
                    edgeq(edgeq(iee,7),6)=edgeq(iee,5); % So give the reverse exit a valid path
                end
                jee=edgeq(edgeq(iee,7),5); % alternate exit edge if flipping
            end
            if edgeq(iee,3)
                error('Edge already in use');
            end
            edgeq(iee,3)=iff;
            edgeq(edgeq(iee,7),4)=iff;      % mark right face in reverse edge
            nee=nee+1;
            flip=abs(flip-edgeq(iee,8));
            if flip   % use alternate edge
                edgeq(iee,6)=jee;           % mark actual edge used
            end
            pee=iee;            % save previous iee value
            iee=edgeq(iee,6);   % step around the face
        end
        fpt(iff,1)=nee;
    end
else  % a normal two-sided polygon
    while iff<f
        iee=find(edgeq(:,3)==0,1); % find an edge without a left face
        iff=iff+1;
        fpt(iff,2)=iee;
        nee=0;
        while edgeq(iee,3)~=iff
            if edgeq(iee,3)
                error('Edge already in use');
            end
            edgeq(iee,3)=iff;
            edgeq(edgeq(iee,7),4)=iff;      % mark right face in reverse edge
            nee=nee+1;
            iee=edgeq(iee,6);   % step around the face
        end
        fpt(iff,1)=nee;
    end
end
% [(1:v*m)' edgeq]

vlist(:,4)=sqrt(sum(vlist(:,1:3).^2,2));        % fill in the radii
vlist(:,6)=(1:v)';
% now fill in the flist array
% [x y z d n e t] for each face
flist=zeros(f,7);
[ia,ib]=sort(edgeq(:,3));           % finst all the edges belonging to each face
ix=[1; 1+find(ia(2:end)>ia(1:end-1)); 2*e+1];
flist(:,6)=ib(ix(1:f));         % pointer to first edge for each face
flist(:,5)=ix(2:end)-ix(1:end-1);       % size of each face
% edgeq: 1=v1 2=v2 3=f1 4=f2 5=ev1 6=ef1 7=er 8=z 9=sf 10=sv]
tedge=flist(:,6);           % this edge
nedge=edgeq(tedge,6);       % next edge
vmid=vlist(edgeq(tedge+2*e*(1-edgeq(tedge,10))),1:3);
flist(:,1:3)=cross(vlist(edgeq(nedge+2*e*(1-edgeq(nedge,10))),1:3)-vmid,vlist(edgeq(tedge+2*e*edgeq(tedge,10)),1:3)-vmid,2);
flist(:,1:3)=flist(:,1:3)./repmat(sqrt(sum(flist(:,1:3).^2,2)),1,3);
flist(:,4)=sum(flist(:,1:3).*vmid,2);

% now check for dual

if dual
    vlisto=vlist;        % save vlist
    vlist=flist;         % old faces are new vertices
    flist=vlisto;
    % edgeq: 1=v1 2=v2 3=f1 4=f2 5=ev1 6=ef1 7=er 8=z 9=sf 10=sv]
    edgeq=edgeq(:,[3 4 2 1 6 5 7 8 10 9]);  % swap vertices and faces in the edge list
    erev=edgeq(:,7);        % reverse edge
    edgeq(:,6)=erev(edgeq(erev,6));
    flist(:,6)=erev(flist(:,6));
    if hemiQ             % cannot invert through centre
        error('Cannot take dual');
    else
        vlist(:,4)=vlist(:,4).^(-1);
    end
    vlist(:,1:3)=vlist(:,1:3).*repmat(vlist(:,4),1,3);
    v=size(vlist,1);
    vlist(:,1:3)=vlist(:,1:3)-repmat(mean(vlist(:,1:3),1),v,1);
    f=size(flist,1);            % recalculate face positions from scratch
    tedge=flist(:,6);           % this edge
    nedge=edgeq(tedge,6);       % next edge
    vmid=vlist(edgeq(tedge+2*e*(1-edgeq(tedge,10))),1:3);
    flist(:,1:3)=cross(vlist(edgeq(nedge+2*e*(1-edgeq(nedge,10))),1:3)-vmid,vlist(edgeq(tedge+2*e*edgeq(tedge,10)),1:3)-vmid,2);
    flist(:,1:3)=flist(:,1:3)./repmat(sqrt(sum(flist(:,1:3).^2,2)),1,3);
    flist(:,4)=sum(flist(:,1:3).*vmid,2);
end

info.snub=snubQ>0;
info.onesided=onesidedQ>0;
info.hemi=hemiQ>0;

% create Wythoff string
wystr='';
j=0;
for i=1:4
    if i==abs(wy(1))
        wystr=[wystr '| '];
        j=1;
    else
        if (wyd(i-j)==1)
            wystr=[wystr num2str(wyn(i-j)) ' '];
        else
            wystr=[wystr num2str(wyn(i-j)) '/' num2str(wyd(i-j)) ' '];
        end
    end
end
info.wythoff=wystr(1:end-1);
info.vef=[v e f];
info.chi=chi;   % Euler Characteristic

% now draw an image if requested

if ~nargout || any(md=='g')
    clf;
    if any(md=='t')
        [nii,nij]=rat(nimi(:,1));      % face sizes as rational numbers
        veca=0:1;
        plot(0,0,'ok');
        hold on
        for j=1:m
            i=rotsnub(j);
            veca=veca(2)*[1 exp(2i*gamma(i))];
            th=(0:nii(i))*2*pi/nimi(i,1);
            sc=[sin(th); 1-cos(th)];
            scv=sc(:,[2 end-1]);  % first and last vertices
            veci=veca/scv;
            plot(real(veci*sc), imag(veci*sc),'-k',real(veci*scv), imag(veci*scv),':b');
        end
        hold off
    else
        if (all(wyd==1) && ~any(md=='w')) || any(md=='f')

            % now render the polyhedron

            for i=1:f
                %     patch('Faces',fa{i},'Vertices',vlist,'FaceVertexCData',rgb(i,:),'FaceColor','flat')
                fa=zeros(flist(i,5),1);
                ix=flist(i,6);
                for j=1:numel(fa)
                    fa(j)=edgeq(ix,1+edgeq(ix,10));   % first vertex of this edge
                    ix=edgeq(ix,6); % next edge around the face
                end
                patch(vlist(fa,1),vlist(fa,2),vlist(fa,3),1-((flist(i,1:3)+1)/2).^2,'facealpha',0.95);
            end
        else
            edges=edgeq(edgeq(:,1)<edgeq(:,2),1:2);
            plot3(reshape(vlist(edges,1),size(edges))',reshape(vlist(edges,2),size(edges))',reshape(vlist(edges,3),size(edges))','k');
        end
        if any(md=='v')
            for i=1:v
                text(vlist(i,1),vlist(i,2),vlist(i,3),sprintf('%d',i));
            end
        end
    end
    axis equal
    if dual
        tdual='Dual ';
    else
        tdual='';
    end
    if wyidx==round(wyidx) && wyidx>=1 && wyidx<=size(wys,1)
        title(sprintf('%s%s: %s (%s)',tdual,info.wythoff,names{wyidx,1:2}));
    else
        title(sprintf('%s%s',tdual,info.wythoff));
    end
end
