%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

%% Load pAll
%-------------------------------------------------------------------------%
load("../../Algorithms/pAll")

%% Save new pAll
%-------------------------------------------------------------------------%
save("../../Algorithms/pAll","pAll")

%% Create 2D Airplane Reflectivity Function
%-------------------------------------------------------------------------%
p.pxy = imread("airplane.png");
p.pxy = imresize(p.pxy,[32,32]);
p.pxy = double(mean(p.pxy,3));
p.pxy = p.pxy - min(p.pxy(:));

[p.xLim,p.yLim] = size(p.pxy);

p.pxy(p.pxy>1) = 1;
p.pxy(p.pxy<1) = 0;

p.xT = linspace(-0.1,0.1,p.xLim);

p.yT = linspace(-0.1,0.1,p.yLim);

% Look at the Reflectivity Function
figure
imagesc(squeeze(p.yT),squeeze(p.xT),squeeze(p.pxy))
title("Reflectivity Function")
xlabel("x")
ylabel("y")

%% Create 3D Airplane Reflectivity Function
%-------------------------------------------------------------------------%
p.pxy = imread("airplane.png");
p.pxy = imresize(p.pxy,[32,32]);
p.pxy = double(mean(p.pxy,3));
p.pxy = imrotate(p.pxy,180);
p.pxy = p.pxy - min(p.pxy(:));

[p.xLim,p.yLim] = size(p.pxy);

p.thr = 1;

p.pxy(p.pxy>p.thr) = 1;
p.pxy(p.pxy<p.thr) = 0;

p.xT = linspace(-0.1,0.1,p.xLim);
p.yT = linspace(-0.1,0.1,p.yLim);

p.zLim = 100;
p.zT = linspace(0,1,p.zLim);

%% Put the Airplane different places in the scene
%-------------------------------------------------------------------------%
p.pxyz = zeros([size(p.pxy) p.zLim]);
iParams.z0 = 500;
[~,iy0] = min(abs(p.zT - iParams.z0*1e-3));
p.pxyz(:,:,iy0) = p.pxy;

y1 = 750;
[~,iy1] = min(abs(p.zT - y1*1e-3));
p.pxyz(:,:,iy1) = imrotate(p.pxy,90);

y2 = 250;
[~,iy2] = min(abs(p.zT - y2*1e-3));
p.pxyz(:,:,iy2) = imrotate(p.pxy,-90);

clear iz0 z1 iz1 z2 iz2
pAll.Airplanes = p

%% Create 1D/2D PSF Reflectivity Function
%-------------------------------------------------------------------------%
p.xLim = 225;
p.yLim = 225;

p.pxy = zeros(p.xLim,p.yLim);
p.pxy(round(end/2),round(end/2)) = 1;
% p.pxy(round(end/2),round(end/2)+5) = 1;

p.xT = linspace(-0.1,0.1,p.xLim);
p.yT = linspace(-0.1,0.1,p.yLim);

p.py = zeros(1,p.yLim);
p.py( (end+1)/2 ) = 1;

figure
plot(p.yT,p.py)

% Look at the PSF
figure
mesh(squeeze(p.yT),squeeze(p.xT),squeeze(p.pxy))
title("Reflectivity Function")
xlabel("x")
ylabel("y")
p
%% Create 3D PSF
%-------------------------------------------------------------------------%
p.xLim = 225;
p.yLim = 225;

p.pxy = zeros(p.xLim,p.yLim);
p.pxy(round(end/2),round(end/2)) = 1;
% p.pxy(round(end/2),round(end/2)+5) = 1;

p.xT = linspace(-0.1,0.1,p.xLim);
p.yT = linspace(-0.1,0.1,p.yLim);

p.zLim = 45;
p.zT = linspace(0,2*iParams.z0*1e-3,p.zLim);

p.pxyz = cat(3, zeros( [size(p.pxy) (p.zLim-1)/2] ), p.pxy , zeros( [size(p.pxy) (p.zLim-1)/2 ]));

pAll.PSF = p

%% Create 3D Grid of Point Reflectors
%-------------------------------------------------------------------------%
addpath("../../Algorithms");
load fParamsAll; load iParamsAll;
fParams = fParamsAll.Muhammet;
iParams = iParamsAll.MIMO;
clear fParamsAll iParamsAll

p.zLim = 256;
p.xLim = 256;
p.yLim = 256;
p.zT = linspace(1/p.zLim,1,p.zLim);
p.xT = linspace(-1+2/p.xLim,1,p.xLim);
p.yT = linspace(-1+2/p.yLim,1,p.yLim);
iParams.z0 = 500; % mm

p.pxyz = zeros([length(p.xT),length(p.yT),length(p.zT)]);
p.pxyz((end/2-16):16:(end/2+16),(end/2-16):16:(end/2+16),(end/2-32):32:(end/2+32)) = 1;

pAll.Grid3D = p
% volumeViewer(permute(p.pxyz,[2 1 3]))

%% Look at the non-zero z-ranges
%-------------------------------------------------------------------------%
for izP = 1:p.zLim
    if squeeze(mean(p.pxyz(:,:,izP),[1,2])) > 0
        figure
        mesh(flip(squeeze(p.xT)),flip(squeeze(p.yT)),flipud(squeeze(p.pxyz(:,:,izP)).'))
        title("Reflectivity Function")
        view(2)
        xlabel("x")
        ylabel("y")
        title("2D Image at z = " + p.zT(izP))
    end
end

clear izP