%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function isarData = CSAR_1D_createEcho_SISO(iParams,fParams,p)

%% Define Optional Parameters
%-------------------------------------------------------------------------%
if ~isfield(iParams,"showP")
    iParams.showP = true;
end

%% Maintain everything in the size of:
% nAngMeasurement x adcSample
%-------------------------------------------------------------------------%
isarData = complex(zeros(iParams.nAngMeasurement, fParams.adcSample));

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
thetaM_rad = (0:iParams.tStepM_deg:((iParams.nAngMeasurement-1)*iParams.tStepM_deg))*pi/180;
thetaM_rad = reshape(thetaM_rad,[],1);

%% Look at the Reflectivity of the Target Scene
%-------------------------------------------------------------------------%
if ~isfield(p,"pxz")
    if mod(p.yLim,2) == 0
        pxz = squeeze(p.pxyz(:,end/2,:));
    else
        pxz = squeeze(p.pxyz(:,(end+1)/2,:));
    end
else
    pxz = p.pxz;
end

if iParams.showP
    figure;
    mesh(p.zT,p.xT*1e3,pxz);
    title("Reflectivity Function of Target Scene")
    xlabel("z (m)")
    ylabel("x (mm)")
end

%% Create Echo Signal
%-------------------------------------------------------------------------%
R0 = iParams.R0_mm*1e-3; % m
if ~isempty(gcp('nocreate')) % if parallel pool is open
    parfor ixP = 1:p.xLim
        for izP = 1:p.zLim
            if pxz(ixP,izP) > 1e-8
                R = sqrt( (R0*cos(thetaM_rad) - p.xT(ixP)).^2 + (R0*sin(thetaM_rad) - p.zT(izP)).^2 );
                isarData = isarData + pxz(ixP,izP) .* (R.^(-2)) .* exp(1j*2*k.*R);
            end
        end
    end
else
    for ixP = 1:p.xLim
        for izP = 1:p.zLim
            if pxz(ixP,izP) > 1e-8
                R = sqrt( (R0*cos(thetaM_rad) - p.xT(ixP)).^2 + (R0*sin(thetaM_rad) - p.zT(izP)).^2 );
                isarData = isarData + pxz(ixP,izP) .* (R.^(-2)) .* exp(1j*2*k.*R);
            end
        end
    end
end

% Alternate Method (Does not take advantage of the fact that p(x,z) is
% known and we can ignore all the zeros in p(x,z)
% p.xT = reshape(p.xT,[],1);
% p.zT = reshape(p.zT,1,[]);
% 
% lengthk = length(k);
% parfor iT = 1:length(thetaM)
%     for iK = 1:lengthk
%         R = sqrt( (R0*cos(thetaM(iT)) - p.xT).^2 + (R0*sin(thetaM(iT)) - p.zT).^2 );
%         isarData2(iT,iK) = sum( pxz .* (R.^(-2)) .* exp(1j*2*k(iK).*R), 'all');
%     end
% end
