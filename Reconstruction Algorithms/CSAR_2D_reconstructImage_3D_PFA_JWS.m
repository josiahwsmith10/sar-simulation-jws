%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function [csarImage3D,xRangeT_m,yRangeT_m,zRangeT_m] =  CSAR_2D_reconstructImage_3D_PFA_JWS(csarData3D,iParams,fParams)
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

%% Declare y Synthetic Aperture Vector

y = (0:iParams.nVerMeasurement-1)*iParams.yStepM_mm*1e-3;
if mod(iParams.nVerMeasurement,2) == 0
    y = y - y(end/2);
else
    y = y - y( (end+1)/2 );
end
y = reshape(y,1,1,[]);

%% Declare kY and kR
%-------------------------------------------------------------------------%
kSy = 2*pi/(iParams.yStepM_mm*1e-3);
kY = linspace(-kSy/2,kSy/2,iParams.yU*iParams.nFFT);
kY = reshape(kY,1,1,[]);
kR = sqrt(4*k.^2 - kY.^2);
clear kSy

%% Declare kX and kZ
%-------------------------------------------------------------------------%
kX = 2*k.*cos(theta_rad);
kZ = 2*k.*sin(theta_rad);
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
% kRU = sqrt(kXU.^2 + kZU.^2);
kU = (1/2)*sqrt(kXU.^2 + kZU.^2 + kY.^2);
clear kXU kZU

%% Upsample csarData using interpn
%-------------------------------------------------------------------------%
theta_radUpsample = reshape(linspace(min(theta_rad),max(theta_rad),iParams.xU*iParams.nFFT),[],1);
kUpsample = reshape(linspace(min(k),max(k),iParams.zU*iParams.nFFT),1,[]);
yUpsample = reshape(linspace(min(y),max(y),iParams.yU*iParams.nFFT),1,1,[]);
kRUpsample = sqrt(4*kUpsample.^2 - kY.^2);

% csarDataUpsampled = zeros(length(theta_radUpsample),length(kUpsample),iParams.nVerMeasurement*iParams.yU);

[theta_radG,kG,yG] = ndgrid(theta_rad,k,y);
clear k theta_rad y
[theta_radUpsampleG,kUpsampleG,yUpsampleG] = ndgrid(theta_radUpsample,kUpsample,yUpsample);
clear yUpsample

csarDataUpsampled = interpn(theta_radG,kG,yG,csarData3D,theta_radUpsampleG,kUpsampleG,yUpsampleG,iParams.PFA,0);
clear theta_radG kG yG csarData3D theta_radUpsampleG kUpsampleG yUpsampleG

%% Compute Azimuth Filter: h(Theta,k,kY)
%-------------------------------------------------------------------------%
azimuthFilterFFT = fft(exp(1j*2*kUpsample*(iParams.R0_mm*1e-3).*cos(theta_radUpsample)),[],1);
clear kRUpsample

%% Compute Azimuth Filtered Data: p(theta,k,kY) = IFT[ S(Theta,k,kY) * H*(Theta,k,kY) ]
%-------------------------------------------------------------------------%
azimuthFiltered = ifft(fftshift(fft(fft(csarDataUpsampled,[],1),[],3),3) .* conj(azimuthFilterFFT),[],1);
clear csarDataUpsampled azimuthFilterFFT

%% Interpolate Azimuth Filtered Data to CSAR Image FFT: p(kX,kZ,kY)
%-------------------------------------------------------------------------%

[theta_radUpsampleG,kUpsampleG,~] = ndgrid(theta_radUpsample,kUpsample,kY);
clear theta_radUpsample kUpsample
% Make my own ndgrid because theta_radU and kU are multidimensional already
theta_radUG = repmat(theta_radU,[1,1,size(kU,3)]);
clear theta_radU
kYG = repmat(kY,[size(kU,1),size(kU,2),1]);
clear kY

csarImageFFT = interpn(theta_radUpsampleG,kUpsampleG,kYG,azimuthFiltered,theta_radUG,kU,kYG,iParams.PFA,0);
clear theta_radUpsampleG kUpsampleG kYG azimuthFiltered theta_radUG kU kYG

%% Recover CSAR Image: p(x,z)
%-------------------------------------------------------------------------%
csarImage3D = ifftshift(ifftshift(ifftn(csarImageFFT),1),2);
clear csarImageFFT

%% Reorient CSAR Image: p(x,z,y) -> p(x,y,z)
%-------------------------------------------------------------------------%
csarImage3D = permute(csarImage3D,[1,3,2]);

%% Crop the Image for Related Region
%-------------------------------------------------------------------------%
xRangeT_m = (-(size(csarImage3D,1))/2:(size(csarImage3D,1)/2-1))*(2*pi/(kXUsteps*iParams.xU*iParams.nFFT));
yRangeT_m = (-(size(csarImage3D,2))/2:(size(csarImage3D,2)/2-1))*iParams.yStepM_mm*1e-3;
zRangeT_m = (-(size(csarImage3D,3))/2:(size(csarImage3D,3)/2-1))*(2*pi/(kZUsteps*iParams.zU*iParams.nFFT));

if (iParams.xyzSizeT_m ~= -1)
    indX = xRangeT_m>(-iParams.xyzSizeT_m/2) & xRangeT_m<(iParams.xyzSizeT_m/2);
    indY = yRangeT_m>(-iParams.xyzSizeT_m/2) & yRangeT_m<(iParams.xyzSizeT_m/2);
    indZ = zRangeT_m>(-iParams.xyzSizeT_m/2) & zRangeT_m<(iParams.xyzSizeT_m/2);
    xRangeT_m = xRangeT_m(indX);
    yRangeT_m = yRangeT_m(indY);
    zRangeT_m = zRangeT_m(indZ);
    csarImage3D = csarImage3D(indX,indY,indZ);
    clear indX indZ
end
clear kXUsteps kZUsteps

%% Display the Result
%-------------------------------------------------------------------------%
volumeViewer(abs(csarImage3D));