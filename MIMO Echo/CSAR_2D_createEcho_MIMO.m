%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function csarData = CSAR_2D_createEcho_MIMO(iParams,fParams,p)
%% Declare s(theta,k,y)
%-------------------------------------------------------------------------%
csarData = zeros(iParams.nAngMeasurement, fParams.adcSample, iParams.nRx*iParams.nTx*iParams.nVerMeasurement);

%% Declare Wavenumber Vector
%-------------------------------------------------------------------------%
f0 = fParams.f0 + fParams.ADCStartTime*fParams.K; % This is for ADC sampling offset
f = f0 + (0:fParams.adcSample-1)*fParams.K/fParams.fS; % wideband frequency

c = 299792458; % physconst('lightspeed'); in m/s
k = 2*pi*f/c;
k = reshape(k,1,[]);
clear f f0

%% Declare theta Synthetic Aperture Vector
%-------------------------------------------------------------------------%
thetaM_rad = 0:iParams.tStepM_deg:((iParams.nAngMeasurement-1)*iParams.tStepM_deg);
thetaM_rad = reshape(thetaM_rad,[],1);

%% Get Array Dimensions
%-------------------------------------------------------------------------%
[yr_orig_m,yt_orig_m,yv_orig_m] = awr1443ArrayDimensions(false);

yr_orig_m = reshape(yr_orig_m,1,1,[]);
yt_orig_m = reshape(yt_orig_m,1,1,[]);
yv_orig_m = reshape(yv_orig_m,1,1,[]);

%% Create Echo Signal
%-------------------------------------------------------------------------%
pxyz = p.pxyz;
R0 = iParams.R0_mm*1e-3; % m
cthetaM = cos(thetaM_rad);
sthetaM = sin(thetaM_rad);

sizeSarData = size(csarData);

parfor ixP = 1:p.xLim
    for iyP = 1:p.yLim
        for izP = 1:p.zLim
            if pxyz(ixP,iyP,izP) > 1e-8
                stot = zeros(sizeSarData);
                for delz = (-(iParams.nVerMeasurement-1)/2:1:(iParams.nVerMeasurement-1)/2) * 2*iParams.lambda_mm * 1e-3
                    
                    yr = yr_orig_m + delz - mean(yv_orig_m);
                    yt = yt_orig_m + delz - mean(yv_orig_m);
                    
                    Rr = sqrt( (R0*cthetaM - p.xT(ixP)).^2 + (R0*sthetaM - p.zT(izP)).^2 + (yr - p.yT(iyP)).^2);
                    Rt = sqrt( (R0*cthetaM - p.xT(ixP)).^2 + (R0*sthetaM - p.zT(izP)).^2 + (yt - p.yT(iyP)).^2);
                    
                    s = zeros(length(thetaM_rad),length(k),iParams.nRx*iParams.nTx);
                    s(:,:,1:iParams.nRx*iParams.nTx/2) = pxyz(ixP,iyP,izP) .* (Rt(:,:,1).*Rr).^(-1) .* exp(1j*k.*(Rt(:,:,1) + Rr) - 1j*pi*fParams.K/(c^2).*(Rt(:,:,1).*Rr).^2);
                    s(:,:,iParams.nRx*iParams.nTx/2+1:end) = pxyz(ixP,iyP,izP) .* (Rt(:,:,2).*Rr).^(-1) .* exp(1j*k.*(Rt(:,:,2) + Rr) - 1j*pi*fParams.K/(c^2).*(Rt(:,:,2).*Rr).^2);
                    
                    idxStart = (delz/(2*iParams.lambda_mm*1e-3) + (iParams.nVerMeasurement-1)/2)*iParams.nTx*iParams.nRx + 1;
                    idxEnd = (delz/(2*iParams.lambda_mm*1e-3) + (iParams.nVerMeasurement-1)/2 + 1)*iParams.nTx*iParams.nRx;
                    
                    stot(:,:,idxStart:idxEnd) = s;
                end
                
                csarData = csarData + stot;
            end
        end
    end
end