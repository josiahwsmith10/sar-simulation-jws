%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

%% 1. Add the Necessary Folders to Path (Run First)
%-------------------------------------------------------------------------%
addpath(genpath("../"))

%% 2. Load iParams, fParams, and p
%-------------------------------------------------------------------------%
load fParamsAll; load iParamsAll; load pAll
fParams = fParamsAll.v0;                    % Frequency Parameters
iParams = iParamsAll.SISO_CSAR;             % Image and Scanning Parameters (360deg of rotation)
p = pAll.CSAR_Grid3D;                       % Reflectivity p(x,y,z) parameters
clear fParamsAll iParamsAll pAll

%% 3. Get the SISO-CSAR Echo Signal s(theta,k,y): csarData
%-------------------------------------------------------------------------%
% p.pxyz(end/2,end/2,end/2) = 1;
p.pxyz(end/4,end/2,end/4) = 1;
iParams.showP = true;
csarData = CSAR_2D_createEcho_SISO(iParams,fParams,p);

%% 4. 3D Image Reconstruction using Polar Formatting Algorithm (PFA)
%-------------------------------------------------------------------------%
iParams.nFFT = 512;
iParams.xyzSizeT_m = 0.4;
iParams.PFA = 'linear';
iParams.xU = 2;
iParams.yU = 1;
iParams.zU = 2;
csarImage3D_PFA = CSAR_2D_reconstructImage_3D_PFA_JWS(csarData,iParams,fParams);