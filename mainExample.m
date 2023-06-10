close all; clear all; clc;
obj = CoordConv(false); % Instantiate Class Object
%% Configure cell arrays for parameter input to Class CoordConv Method
A = {[3,5,8],[5,3*pi/4,4],[7,pi/3,pi/4]};   % radians input
A_ = {[3,5,8],[5,135,4],[7,60,45]};         % degrees input
% Constants to loop through A
Coords = {["xyz", "cyl", "sph"]};
outCKEY = {[2, 3], [1, 3], [1, 2]};
% Loop 1 convert to/from base system in radians and degrees per vector in A
for i=1:length(A)
    Vrad = A{i}; 
    Vdeg = A_{i};  
    inC = Coords{1}(i);
    currKEY = outCKEY{i};
    for k=1:2        
        outC = Coords{1}(currKEY(k));
        [B(1),B(2),B(3)] = obj.convCoord_r(Vrad,inC,outC);
        [C(1),C(2),C(3)] = obj.convCoord_r([B(1),B(2),B(3)],outC,inC);    
        [D(1),D(2),D(3)] = obj.convCoord_d(Vdeg,inC,outC);
        [E(1),E(2),E(3)] = obj.convCoord_d([D(1),D(2),D(3)],outC,inC);
    end
end

