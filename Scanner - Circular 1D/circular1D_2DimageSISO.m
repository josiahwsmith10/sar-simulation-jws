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

%% 3. Get the SISO-CSAR Echo Signal s(theta,k): csarData
%-------------------------------------------------------------------------%
% p.pxyz(end/2,end/2,end/2) = 1;
p.pxyz(end/4,end/2,end/4) = 1;
iParams.showP = true;

csarData = CSAR_1D_createEcho_SISO(iParams,fParams,p);

%% 4. 2D Image Reconstruction by Back Projection Algorithm (BPA)
%-------------------------------------------------------------------------%
iParams.isAmplitudeFactor = true;
csarImage2D_BPA = CSAR_1D_reconstructImage_2D_BPA(s_theta_k,iParams,fParams,p);

%% 5. 2D Image Reconstruction Polar Formatting Algorithm (PFA)
%-------------------------------------------------------------------------%
iParams.nFFT = 512;
iParams.xzSizeT_m = 0.4;
iParams.xU = 4;
iParams.zU = 4;
csarImage2D_PFA = CSAR_1D_reconstructImage_2D_PFA(csarData,iParams,fParams);

%% 6. 2D Image Reconstruction Segmented Polar Formatting Algorithm (PFA)
%-------------------------------------------------------------------------%
iParams.nFFT = 512;
iParams.xzSizeT_m = 0.4;
iParams.xU = 4;
iParams.zU = 4;
csarImage2D_SegmentedPFA = CSAR_1D_reconstructImage_2D_SegmentedPFA(csarData,iParams,fParams);
