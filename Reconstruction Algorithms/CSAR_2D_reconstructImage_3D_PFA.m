%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function [csarImage3D,xRangeT_m,yRangeT_m,zRangeT_m] =  CSAR_2D_reconstructImage_3D_PFA(csarData3D,iParams,fParams)
% Reconstructs a 3D image from CSAR echo data using the polar reformatting
% algorithm (PFA), performed by built-in MATLAB interpolation.

%% Input Arguments
% csarData: iParams.nAngMeasurement x fParams.adcSample x iParams.nVerMeasurement
%
% iParams: struct with fields
%   nAngMeasurements    :   Number of angular measurements
%   tStepM_deg          :   Step size of angular measurements in degrees
%   nVerMeasurements    :   Number of vertical measurements
%   yStepM_mm           :   Step size of vertical axis in mm
%   R0_mm               :   Scanning radius
%   xyzSizeT_m          :   Expected size of the target in the xz domain in m
%   scanName            :   Name of the scan
%   PFA                 :   Interpolation method for PFA
%   nFFT                :   Number of FFT bins to use (before upsampling)
%   xU                  :   Upsampling factor in x-axis
%   yU                  :   Upsampling factor in y-axis
%   zU                  :   Upsampling factor in z-axis
%
% fParams: struct with fields
%   K                   :   Chrip Slope (Hz/s)
%   fS                  :   Sampling frequency (1/s)
%   adcSample           :   Number of samples
%   ADCStartTime        :   Time it takes to start sampling (s)
%   f0                  :   Chirp start frequency (Hz)

%% Declare Optional Parameters
%-------------------------------------------------------------------------%
if ~isfield(iParams,'scanName')
    iParams.scanName = 'CSAR-PFA';
end
if ~isfield(iParams,'PFA')
    iParams.PFA = 'nearest';
end
if ~isfield(iParams,'xU')
    iParams.xU = 4;
end
if ~isfield(iParams,'yU')
    iParams.yU = 1;
end
if ~isfield(iParams,'zU')
    iParams.zU = 4;
end
if ~isfield(iParams,'xyzSizeT_m')
    iParams.xyzSizeT_m = 0.4;
end

%% Declare Wavenumber Vector
%-------------------------------------------------------------------------%
f0 = fParams.f0 + fParams.ADCStartTime*fParams.K; % This is for ADC sampling offset
f = f0 + (0:fParams.adcSample-1)*fParams.K/fParams.fS; % wideband frequency

c = 299792458; % physconst('lightspeed'); in m/s
k = 2*pi*f/c;
k = reshape(k,1,[]);
clear f f0 c

%% Declare theta_rad Synthetic Aperture Vector
%-------------------------------------------------------------------------%
theta_rad = (0:iParams.tStepM_deg:((iParams.nAngMeasurement-1)*iParams.tStepM_deg))*pi/180;
theta_rad = reshape(theta_rad,[],1);

%% Declare kY and kR
%-------------------------------------------------------------------------%
%%% MEY Correction Start
% kY = reshape(linspace(-2*max(k),2*max(k),iParams.yU*iParams.nFFT),1,1,[]);
kSy = 2*pi/(iParams.yStepM_mm*1e-3);
kY = linspace(-kSy/2,kSy/2,iParams.nFFT);
kY = reshape(kY,1,1,[]);
%%% MEY Correction End
kR = sqrt(4*k.^2 - kY.^2);

%% Declare kX and kZ
%-------------------------------------------------------------------------%
kX = 2*kR.*cos(theta_rad);
kZ = 2*kR.*sin(theta_rad);
clear kR

kXmax = max(kX,[],'all');
kZmax = max(kZ,[],'all');
kXmin = min(kX,[],'all');
kZmin = min(kZ,[],'all');
clear kX kZ

%% Declare kXU and kZU
%-------------------------------------------------------------------------%
% According to Ozdemir
kXU = reshape(linspace(iParams.xU*kXmin,iParams.xU*kXmax,iParams.xU*iParams.nFFT),[],1);
kZU = reshape(linspace(iParams.zU*kZmin,iParams.zU*kZmax,iParams.zU*iParams.nFFT),1,[]);
clear kXmax kZmax kXmin kZmin

kXUsteps = mean(diff(kXU));
kZUsteps = mean(diff(kZU));

%% Declare Uniform theta_radU and kU for Interpolation
%-------------------------------------------------------------------------%
theta_radU = atan2(kZU,kXU);
kRU = sqrt(kXU.^2 + kZU.^2);
clear kXU kZU

%% Upsample csarData using interp2
%-------------------------------------------------------------------------%
kUpsample = reshape(linspace(min(k),max(k),iParams.zU*iParams.nFFT),1,[]);
theta_radUpsample = reshape(linspace(min(theta_rad),max(theta_rad),iParams.xU*iParams.nFFT),[],1);
kRUpsample = sqrt(4*kUpsample.^2 - kY.^2);
clear kY

csarDataUpsampled = zeros(length(theta_radUpsample),length(kUpsample),iParams.nVerMeasurement*iParams.yU);

parfor iY = 1:size(csarDataUpsampled,3)
    csarDataUpsampled(:,:,iY) = interp2(k,theta_rad,csarData3D(:,:,iY),kUpsample,theta_radUpsample,iParams.PFA,0);
end
clear k theta_rad csarData3D

%% Compute Azimuth Filter: h(Theta,k,kY)
%-------------------------------------------------------------------------%
azimuthFilterFFT = fft(exp(1j*2*kRUpsample*(iParams.R0_mm*1e-3).*cos(theta_radUpsample)),[],1);
clear kRUpsample

%% Compute Azimuth Filtered Data: p(theta,k,kY) = IFT[ S(Theta,k,kY) * H*(Theta,k,kY) ]
%-------------------------------------------------------------------------%
%%% MEY Correction Start
% azimuthFiltered = ifft(fft(fft(csarDataUpsampled,[],1),[],3) .* conj(azimuthFilterFFT));
azimuthFiltered = ifft(fftshift(fft(fft(csarDataUpsampled,[],1),iParams.nFFT,3),3) .* conj(azimuthFilterFFT));
%%% MEY Correction End
clear csarDataUpsampled azimuthFilterFFT

%% Interpolate Azimuth Filtered Data to CSAR Image FFT: p(kX,kZ)
%-------------------------------------------------------------------------%

csarImageFFT = zeros(length(theta_radUpsample),length(kUpsample),iParams.yU*iParams.nFFT);
parfor iY = 1:size(csarImageFFT,3)
    csarImageFFT(:,:,iY) = interp2(kUpsample,theta_radUpsample,azimuthFiltered(:,:,iY),kRU,theta_radU,iParams.PFA,0);
end
clear theta_radU kRU theta_radUpsample kUpsample

%% Recover CSAR Image: p(x,z)
%-------------------------------------------------------------------------%
csarImage3D = ifftshift(ifftn(csarImageFFT));
clear csarImageFFT

%% Reorient CSAR Image: p(x,z,y) -> p(x,y,z)
%-------------------------------------------------------------------------%
csarImage3D = permute(csarImage3D,[1,3,2]);

%% Crop the Image for Related Region
%-------------------------------------------------------------------------%
xRangeT_m = (-size(csarImage3D,1)/2:(size(csarImage3D,1)/2-1))*(2*pi/(kXUsteps*iParams.xU*iParams.nFFT));
yRangeT_m = (1:iParams.nFFT)*iParams.lambda/(2*iParams.nFFT);
zRangeT_m = (-size(csarImage3D,2)/2:(size(csarImage3D,2)/2-1))*(2*pi/(kZUsteps*iParams.zU*iParams.nFFT));
clear kXUSteps kZUSteps

if (iParams.xzSizeT_m ~= -1)
    indX = xRangeT_m>(-iParams.xyzSizeT_m/2) & xRangeT_m<(iParams.xyzSizeT_m/2);
    indY = yRangeT_m>(-iParams.xyzSizeT_m/2) & yRangeT_m<(iParams.xyzSizeT_m/2);
    indZ = zRangeT_m>(-iParams.xyzSizeT_m/2) & zRangeT_m<(iParams.xyzSizeT_m/2);
    xRangeT_m = xRangeT_m(indX);
    yRangeT_m = yRangeT_m(indY);
    zRangeT_m = zRangeT_m(indZ);
    csarImage3D = csarImage3D(indX,indY,indZ);
    clear indX indZ
end

%% Display the Result
%-------------------------------------------------------------------------%
volumeViewer(abs(csarImage3D));