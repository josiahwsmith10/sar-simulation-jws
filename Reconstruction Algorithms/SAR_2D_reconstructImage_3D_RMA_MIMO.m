%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function sarImage = SAR_2D_reconstructImage_3D_RMA_MIMO(sarData,iParams,fParams)
% Reconstructs a 3D image from SAR echo data using the range migration
% algorithm (RMA). Expects a 2D rectilinear scan across the x-axis with a
% SISO transceiver and the y-axis with a MIMO transceiver

%% Input Arguments
%-------------------------------------------------------------------------%
% sarData: iParams.nVerMeasurement*iParams.nTx*iParams.nRx x iParams.nUsefulHorMeasurement x fParams.adcSample single, double, or gpuArray
% sarData: s(y,x,k)
%
% iParams: struct with fields
%   z0_mm               :   Distance of target
%   xStepM_mm           :   Step size of x-axis in mm
%   yStepM_mm           :   Step size of y-axis in mm
%   scanName            :   Name of the scan
%   Stolt               :   Interpolation method for RMA
%
% fParams: struct with fields
%   K                   :   Chrip Slope (Hz/s)
%   fS                  :   Sampling frequency (1/s)
%   adcSample           :   Number of samples
%   ADCStartTime        :   Time it takes to start sampling (s)
%   f0                  :   Chirp start frequency (Hz)

%% Declare Optional Parameters
%-------------------------------------------------------------------------%
if ~isfield(iParams,"Stolt")
    iParams.Stolt = "v5cubic";
end
if ~isfield(iParams,"mex")
    iParams.mex = true;
end
if 2*fParams.adcSample > iParams.nFFT
    iParams.nFFT = 2*fParams.adcSample;
end

%% Zeropad sarData
%-------------------------------------------------------------------------%
sarDataPadded = sarData;
if (iParams.nFFT > size(sarData,1))
    sarDataPadded = padarray(sarDataPadded,[floor((iParams.nFFT-size(sarData,1))/2) 0],0,'pre');
    sarDataPadded = padarray(sarDataPadded,[ceil((iParams.nFFT-size(sarData,1))/2) 0],0,'post');
else
    iParams.nFFT = size(sarData,1);
end
if (iParams.nFFT > size(sarData,2))
    sarDataPadded = padarray(sarDataPadded,[0 floor((iParams.nFFT-size(sarData,2))/2)],0,'pre');
    sarDataPadded = padarray(sarDataPadded,[0 ceil((iParams.nFFT-size(sarData,2))/2)],0,'post');
else
    iParams.nFFT = size(sarData,2);
end
clear sarData

%% Obtain S(kY,kX,k): sarDataFFT
%-------------------------------------------------------------------------%
sarDataFFT = fftshift(fftshift(fft(fft(conj(sarDataPadded),iParams.nFFT,1),iParams.nFFT,2),1),2);
clear sarDataPadded

%% Define Interpolation Parameters
%-------------------------------------------------------------------------%
f0 = fParams.f0 + fParams.ADCStartTime*fParams.K; % This is for ADC sampling offset
f = f0 + (0:fParams.adcSample-1)*fParams.K/fParams.fS; % wideband frequency
c = 299792458; % m/s
k = reshape(2*pi*f/c,1,1,[]); % Wavenumber Vector

kSx = 2*pi/(iParams.xStepM_mm*1e-3);
kX = reshape(linspace(-kSx/2,kSx/2,iParams.nFFT),1,[]);

kSy = 2*pi/(iParams.yStepM_mm*1e-3);
kY = reshape(linspace(-kSy/2,kSy/2,iParams.nFFT),[],1);

% Consider only visible spectrum of kZ
kZU = reshape(linspace(0,2*max(k),iParams.nFFT),1,1,[]);
KU = 1/2 * sqrt(kY.^2 + kX.^2 + kZU.^2);
clear c f f0 kSx kX kSy kY kZU

%% Attempt Stolt Interpolation S(kY,kX,k) -> S(kY,kX,kZ)
%-------------------------------------------------------------------------%
disp("Attempting 3D Stolt interpolation using " + iParams.Stolt + " method")
tic
if iParams.mex
    sarImageFFT = stoltInterp3D_mex(sarDataFFT,k,KU);
else
    sarImageFFT = stoltInterp3D(sarDataFFT,k,KU);
end
disp("Stolt interpolation completed in " + toc + " seconds")
% Works with: linear,nearest,next,previous,v5cubic
sarImageFFT(isnan(sarImageFFT)) = 0;
clear sarData KU k fParams sizeKU2

%% Recover the Reflectivity Function
%-------------------------------------------------------------------------%
sarImage = ifftn(sarImageFFT,[iParams.nFFT,iParams.nFFT,iParams.nFFT]);

%% Display the Result
%-------------------------------------------------------------------------%
volumeViewer(abs(sarImage));