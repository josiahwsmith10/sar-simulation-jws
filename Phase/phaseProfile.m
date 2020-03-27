%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function phaseProfile(sarDataIn)

% Can handle both 2D and 3D sarDataIn arrays

if ndims(sarDataIn) == 3
    rangeData3D = fft(sarDataIn,[],3);
    rangeData3D = permute(rangeData3D,[3,1,2]);
    tempRange = mag2db(squeeze(abs(rangeData3D(1:50,:))));
    maxIdx = zeros(1,size(tempRange,2));
    
    for ii = 1:size(tempRange,2)
        [~,maxIdx(ii)] = max(tempRange(:,ii));
    end
    
    rangeBinIdx = input("What is the range bin index? (Recommended: " + mode(maxIdx) + ") ");
    
    rangeDataPeak = squeeze(rangeData3D(rangeBinIdx,:,:));
    
    phaseData2D = unwrap(unwrap(angle(rangeDataPeak),[],1),[],2);
    
    figure
    subplot(121)
    plot(phaseData2D(:,end/2))
    title("Phase Profile at Midpoint")
    
    subplot(122)
    mesh(phaseData2D)
    title("2D Phase Profile")
elseif ismatrix(sarDataIn)
    rangeData2D = fft(sarDataIn,[],2);
    rangeData2D = permute(rangeData2D,[2,1]);
    tempRange = mag2db(squeeze(abs(rangeData2D(1:50,:))));
    maxIdx = zeros(1,size(tempRange,2));
    for ii = 1:size(tempRange,2)
        [~,maxIdx(ii)] = max(tempRange(:,ii));
    end
    
    rangeBinIdx = input("What is the range bin index? (Recommended: " + mode(maxIdx) + ") ");
    
    rangeDataPeak = squeeze(rangeData2D(rangeBinIdx,:));
    
    figure
    plot(unwrap(unwrap(angle(rangeDataPeak))))
    title("Phase Profile")
    xlabel("Capture Index")
end