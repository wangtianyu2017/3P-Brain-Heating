function [cors,dirs]=photoninitGauss(num,r,zcenter,fl)
% This function initalizes the position and propagation direction of each
% photon at the brain surface, according to the Gaussian beam intensity
% profile at the objective lens back aperture.

% num: number of photons, e.g. number of photons in one photon packet which
% is set to 10000 by nphotons.
% th_out: half focusing angle, defined as arcsin(NA/n).
% zcenter: the z location where the light ray simulation starts
% fl: the distance between the focus and zcenter
% r: the radius of the light cone at objective back aperture, use 1/e^2 radius of Gaussian beam. 

global params;

rinit = sqrt(-0.5*log(rand(1,3*num)))*r; %radial location for 2D Gaussian with 1/e^2 radius of r. Generate redundant so that there is still enough after clipping.
rinit = rinit(rinit <= params.geo.f*params.geo.NA); % clipping by the back aperature
rinit = rinit(1:num); 
thinit = asind(rinit/params.geo.f/params.opt.nWater); %sine condition for inifinitely conjugated system, used to calculate theta angle for each photon.
thinitpos = -180+2*rand(1,num)*180; %initial position (angular coordinate)

%directions are encoded in spherical coordinates:
% thinit = cosd(th_out)+rand(1,num)*(1-cosd(th_out)); %length of the Z-vector
% phiinit=-180+2*rand(1,num)*180; %inital angle in the z plane
%initialize x,y, and z direction matrices such that [x y z] is a unit
%vector
zdir=cosd(thinit);
temp = sqrt(1-zdir.^2);
xdir=-temp.*cosd(thinitpos);
ydir=-temp.*sind(thinitpos);

%simulating scanning by translating the fiber in xy @ zcenter, NOT angular
%scanning. Give each photon a random translation in xy.
FOV = params.geo.FOV;
dx = (FOV*rand(1,num)-0.5*FOV);
dy = (FOV*rand(1,num)-0.5*FOV);

%initializing x, y, and z coordinates at zcenter
xcor=fl*tand(thinit).*cosd(thinitpos);
ycor=fl*tand(thinit).*sind(thinitpos);
zcor=ones(1,num)*zcenter;
cors=[xcor+dx;ycor+dy;zcor];
dirs=[xdir;ydir;zdir];

end