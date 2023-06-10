classdef CoordConv    
%% Use of Class that Converts a point/vector's Coordinate representations
%   1) obj = CoordConv(print);    % instantiate the class
%       print = true = disp("results") in CWindow
%`      print = false = no display, only output converted point/vector
  %% Parameter Limitations
%   V input vector must be a 1x3
%   U output vector will be 3 vector components
%       and can be assigned to a 1x3 vector
%   V and U vectors component order...
%       Cartesian format "xyz" : V(x,y,z)
%       Cylindrical format "cyl" : V(s,theta,z)
%       Spherical format "sph" : V(r,phi,theta)
%   inC is the input vectors coordinate system...
%   outC is the output vectors coordinate system...
%       Cartesian format "xyz"
%       Cylindrical format "cyl"
%       Spherical format "sph"  
  %% Class Methods
%   [u1,u2,u3] = convCoord_d(V, inC, outC)
%       V = [x,y,z] with angles in degrees
%       if inC = "xyz" then outC = "cyl" or outC = "sph"
%       if inC = "cyl" then outC = "xyz" or outC = "sph"
%       if inC = "sph" then outC = "xyz" or outC = "cyl"
%       [u1,u2,u3] is output in degrees
%   [u1,u2,u3] = convCoord_r(V, inC, outC)
%       V = [x,y,z] with angles in radians
%       if inC = "xyz" then outC = "cyl" or outC = "sph"
%       if inC = "cyl" then outC = "xyz" or outC = "sph"
%       if inC = "sph" then outC = "xyz" or outC = "cyl"
%       [u1,u2,u3] is output in radians
  %% Class Properties
%   boolean to turn printing off
%       only the vector is returned w/o disp("results")
  %% Example Application
% %     close all; clear all; clc;
% %     obj = CoordConv(); % Instantiate Class Object
% %     %% Configure cell arrays for parameter input to Class CoordConv Method
% %     A = {[3,5,8],[5,3*pi/4,4],[7,pi/3,pi/4]};   % radians input
% %     A_ = {[3,5,8],[5,135,4],[7,60,45]};         % degrees input
% %     % Constants to loop through A
% %     Coords = {["xyz", "cyl", "sph"]};
% %     outCKEY = {[2, 3], [1, 3], [1, 2]};
% %     % Loop 1 convert to/from base system in radians and degrees per vector in A
% %     for i=1:length(A)
% %         Vrad = A{i}; 
% %         Vdeg = A_{i};  
% %         inC = Coords{1}(i);
% %         currKEY = outCKEY{i};
% %         for k=1:2        
% %             outC = Coords{1}(currKEY(k));
% %             [B(1),B(2),B(3)] = obj.convCoord_r(Vrad,inC,outC);
% %             [C(1),C(2),C(3)] = obj.convCoord_r([B(1),B(2),B(3)],outC,inC);    
% %             [D(1),D(2),D(3)] = obj.convCoord_d(Vdeg,inC,outC);
% %             [E(1),E(2),E(3)] = obj.convCoord_d([D(1),D(2),D(3)],outC,inC);
% %         end
% %     end
  %% PUBLIC------------------------------------------------------------        
properties
    configPrint;    
end
methods
    %% Constructor (bool print)  disp ON/OFF
    function obj = CoordConv(print)
        if nargin < 1
            obj.configPrint = false;
            disp("""Pass True on CoordConv Class Instance to Print Results""")
            disp("    obj = CoordConv(true) instead of obj = CoordConv()")
            disp("    ...default is false, results are not printed")
            disp(" ")
        else
            obj.configPrint = print;
        end
    end
    %% Degree Coordinate Conversion (numerical 1x3 vector,"string","string")
    function [u1, u2, u3] = convCoord_d(obj, V, inC, outC)        
        % Validate Input Parameters
        flag = obj.validateInput(V, inC, outC); 
        if ~flag
            return; % invalid input parameters
        end
        % Convert Point/Vector Coordinates (bool) deg/rad I/O
        [u1, u2, u3] = obj.convCoord(V, inC, outC, true); 
    end
    function [u1, u2, u3] = convCoord_r(obj, V, inC, outC)        
        % Validate Input Parameters
        flag = obj.validateInput(V, inC, outC); 
        if ~flag
            return; % invalid input parameters
        end
        % Convert Point/Vector Coordinates (bool) deg/rad I/O
        [u1, u2, u3] = obj.convCoord(V, inC, outC, false);            
    end
end
%% PRIVATE-----------------------------------------------------------  
methods (Access = private)
    %% Convert a Point (Vector) Coordinate System (bool unit)
    function [u1, u2, u3] = convCoord(obj, V, inC, outC, unit)                           
        % Assign Conversion Function an index per input parameters
        funcIndex = obj.configOutput(inC, outC);
        % Display Conversion Function and Input Parameters in Print ON
        if obj.configPrint
            disp(obj.Prompt(funcIndex));
            disp(V)
        end
        % Convert Degree input into Radians
        if unit    % unit=true=degrees I/O
            if strcmp(inC,"cyl")
                V(2)= V(2) * pi/180;
            elseif strcmp(inC,"sph")
                V(2)= V(2) * pi/180;
                V(3)= V(3) * pi/180;
            end
        end
        % Calculate the Converted Point/Vector in radians        
        [u1, u2, u3] = obj.selectConversion(funcIndex,V);            
        % Convert Radian output into Degrees
        if unit    % unit=false=radians I/O
            if strcmp(outC,"cyl")
                u2 = u2 * 180/pi;
            elseif strcmp(outC,"sph")
                u2 = u2 * 180/pi;
                u3 = u3 * 180/pi;
            end
        end
        % Display Conversion Output
        if obj.configPrint
            disp([u1, u2, u3]);                       
        end
    end
    %% Condition Input Parameters
    function done = validateInput(obj,V, inC, outC)
        % check that V is 1x3 vector             
        s1x3 = [1, 3];
        if ~isvector(V) || ~isequal(size(V),[1 3])
           disp("Input Vector does not identify as a 1x3 Vector")
           done = false;
           return; 
        end            
        % check Input Coordinate System inC is valid
        if ~strcmp(inC, "xyz") && ~strcmp(inC, "cyl") && ~strcmp(inC, "sph")
            disp("Input Coordinate System does not identify as ""xyx"" or ""cyl"" or ""sph""")
            done = false;
            return;
        end
        % check Output Coordinate System outC is valid
        if ~strcmp(outC, "xyz") && ~strcmp(outC, "cyl") && ~strcmp(outC, "sph")
            disp("Input Coordinate System does not identify as ""xyx"" or ""cyl"" or ""sph""")
            done = false;
            return;
        elseif strcmp(inC,outC)
            disp("The Input and Output Coordinate Systems identify as one, no conversion")
            done = false;
            return;
        end
        done = true;
    end           
    %% Set the requested conversion function index
    function funcIndex = configOutput(obj, inC, outC)
        switch inC
            case "xyz"
                switch outC
                    case "cyl"
                        funcIndex = 1;
                    case "sph"
                        funcIndex = 2;
                end
            case "cyl"
                switch outC
                    case "xyz"
                        funcIndex = 3;
                    case "sph"
                        funcIndex = 4;
                end
            case "sph"
                switch outC
                    case "xyz"
                        funcIndex = 5;
                    case "cyl"
                        funcIndex = 6;
                end
        end
    end        
    %% Type of Conversion with output angles in degrees
    function [u1,u2,u3] = selectConversion(obj,i, V)
        switch i
            case 1
                [u1, u2, u3] = obj.xyz2cyl(V);
                %u2 = u2*180/pi;
            case 2
                [u1, u2, u3] = obj.xyz2sph(V);
                %u2 = u2*180/pi;
                %u3 = u3*180/pi;
            case 3
                [u1, u2, u3] = obj.cyl2xyz(V);
            case 4
                [u1, u2, u3] = obj.cyl2sph(V);
                %u2 = u2*180/pi;
                %u3 = u3*180/pi;
            case 5
                [u1, u2, u3] = obj.sph2xyz(V);
            case 6
                [u1, u2, u3] = obj.sph2cyl(V);
                %u2 = u2*180/pi;            
            otherwise
                u3 = 0; u2 = u3; u1 = u2;
        end
    end
    %% Cartesian to Cylindrical
    function [s, theta, z] = xyz2cyl(obj, V) 
        x = V(1); y = V(2); z_ = V(3);
        s = sqrt(x^2+y^2);
        theta = atan2(y,x);
        z = z_;
    end
    %% Cartesian to Spherical
    function [r, phi, theta] = xyz2sph(obj, V) 
        x = V(1); y = V(2); z = V(3);
        r = sqrt(x^2+y^2+z^2);
        phi = atan2(sqrt(x^2+y^2),z);
        theta = atan2(y,x);    
    end
    %% Cylindrical to Cartesian
    function [x, y, z] = cyl2xyz(obj, V) 
        s = V(1); theta = V(2); z_ = V(3);
        x = s*cos(theta);
        y = s*sin(theta);
        z = z_;
    end
    %% Cylindrical to Spherical
    function [r, phi, theta] = cyl2sph(obj, V) 
        s = V(1); theta_ = V(2); z_ = V(3);
        r = sqrt(s^2+z_^2);
        phi = atan2(s,z_);
        theta = theta_;
    end
    %% Spherical to Cartesian
    function [x, y, z] = sph2xyz(obj, V) 
        r = V(1); phi = V(2); theta = V(3);
        x = r*sin(phi)*cos(theta);
        y = r*sin(phi)*sin(theta);
        z = r*cos(phi);
    end
    %% Spherical to Cylindrical
    function [s, theta, z] = sph2cyl(obj, V) 
        r = V(1); phi = V(2); theta_ = V(3);
        s = r*sin(phi);
        theta = theta_;
        z = r*cos(phi);
    end
    %% Output Prompt List
    function prompt = Prompt(obj, i)
    charPrompt = cell(7,1);
        charPrompt{1} = 'Cartesian(x,y,z) to Cylindrical(s,theta,z)';
        charPrompt{2} = 'Cartesian(x,y,z) to Spherical(r,phi,theta)';
        charPrompt{3} = 'Cylindrical(s,theta,z) to Cartesian(x,y,z)';
        charPrompt{4} = 'Cylindrical(s,theta,z) to Spherical(r,phi,theta)';
        charPrompt{5} = 'Spherical(r,phi,theta) to Cartesian(x,y,z)';
        charPrompt{6} = 'Spherical(r,phi,theta) to Cylindrical(s,theta,z)';
        charPrompt{7} = 'Printing is turned off, make ''print'' = true to Print';
        
        strPrompt = cellstr(charPrompt);    
        prompt = char(strPrompt{i});
    end
end
end
