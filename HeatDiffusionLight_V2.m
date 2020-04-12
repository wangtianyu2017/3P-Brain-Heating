function [u_time,u_space,t,r,depth, uRAW]=HeatDiffusionLight_V2(frac_abs,u_step,t_on,t_off,t_max,Power,t_save,r_avg)

%This is a modified version based on Kaspar Podgorski's version (Podgorski
%& Ranganathan, 2016), with the following updates: 
%Organize parameter definitions in a global struct
%in the script calls this function, instead of passing parameters as a
%long list to this function.

%Parameters required for execuation of this function:
% Discretization step for the light model: params.geo.dstep (unit:mm) 
% The starting depth for simulation: zmin=params.geo.zrange(1) (unit:mm);
% The end depth for simulation: zmax=params.geo.zrange(3) (unit:mm);
% The radius in the x-y plane to calculate heat: params.geo.rmax (unit:mm); 
% Thickness of the cover glass: params.geo.d_glass (unit:mm)
% Thickness of the immersion water: params.geo.d_water (unit:mm)
% Thickness of the skull: params.geo.d_skull (unit:mm)
% Density of the brain: params.mech.pTissue (unit:kg/mm^3) 
% Specific heat of brains: params.therm.cTissue (unit:mJ/kg/C) 
% Whole brain thermal diffusivity constant: params.therm.kTissue (unit:mW/mm/C)
% Coverslip thermal diffusivity: params.therm.kGlass (unit:mW/mm/C) 
% Water thermal diffusivity: params.therm.kWater (unit:mW/mm/C)  
% Blood perfusion rate of brain: params.mech.wBlood (unit:/s)
% Specific heat of blood: params.therm.cBlood (unit:mJ/kg/C) 
% Density of blood: params.mech.pBlood (unit:kg/mm^3) 

% Tianyu Wang 11/13/2019

% COMMENTS FROM THE PREVIOUS VERSIONS ARE LISTED BELOW:

%A slightly changed version of Stujenske's code, by Kaspar Podgorski 2015
%adds a cranial window and immersion fluid above the brain volume simulated
%in the earlier code. The water is held at a fixed temperature at the surface, mimicking
%the effect of the objective. Convection and evaporation are not simulated.

%[u_time,t,u_space]=HeatDiffusionLight(frac_abs,params,u_step,t_on,t_off,t_max,Power,t_save,r_avg)
%This function models temperature change (in degrees Celsius) due to light
%from an optical fiber in a block of brain tissue based on the Pennes
%Bio-heat Equation (Pennes, 1948). Must be run after MonteCarloLight.
%
%Inputs:
%1. FRAC_ABS - Fraction absorbed (output of MonteCarloLight)
%2. PARAMS - Model parameters (output of MonteCarloLight)
%3. U_STEP - Discretization step in mm for model (does not need to match
%light model, and in fact, it is better for the light model to have a
%smaller discretization factor. This discretization step should be larger
%because the heat equation solution is more computationally taxing. Default
%is .03 mm.
%4. T_ON - Times when the laser turns on (e.g. [0, 1, 2]) in s. Default is 0.
%5. T_OFF - Times when the laser turns off (e.g. [0.05 1.05 2.05]) in s.
%Default is t_max.
%6. T_MAX - Time to end simulation. Default is 120 or the largest value of
%t_on/t_off (whichever is larger).
%7. POWER - Laser power in mW. Default is 1 mW.
%8. T_SAVE - Times to save heat in u_space in s. Default is t_max.
%9. R_AVG - Radial distance in mm over which to average for u_time. Default
%is .25 mm.
%
%Outputs:
%1. U_TIME - matrix with temperature change (degrees Celsius) over time (s)
%at each depth (depth x time) relative to baseline. 0 indicates no change.
%2. U_SPACE - matrix with temperature values in the plane of the fiber
% [depth x radial distance x length(t_save)]
%3. T - times (s) corresponding to second dimension of u_time.
%4. R - radial distance (mm) corresponding to second dimension of u_space.
%5. DEPTH - depths (mm) corresponding to first dimension of u_time and
%u_space.
%
%Written by Joseph M. Stujenske
%Stujenske et al., Cell Reports

% [u_time,u_space,t,r,depth]=HeatDiffusionLight(frac_abs,params,.03,0,20,20,200,20,.2);

global params;

if nargin<2
    error('You must specify absorbed fraction and params.');
end
if ~exist('u_step','var') || isempty(u_step)
    u_step=.03; %discretization step in mm
end
if ~exist('Power','var') || isempty(Power)
    Power=1; %output power in mW
end
if ~exist('t_on','var') || isempty(t_on)
    t_on=0; %times when the laser is on
end
if ~exist('t_max','var') || isempty(t_max)
    t_max=max([120 max(t_on) max(t_off)]); %time to end simulation
end
if ~exist('t_off','var') || isempty(t_off)
    t_off=t_max; %time when the laser is off; default is on and then off
end
if ~exist('t_save','var') || isempty(t_save)
    t_save=t_max; %output power in mW
end
if ~exist('r_avg','var') || isempty(r_avg)
    r_avg=.25; %output power in mW
end
if r_avg<u_step/2
    r_avg=u_step/2;
    disp('Warning: Changing r_avg value to include central two voxels (minimum r)')
end

%Initialize Variables
dr=params.geo.dstep; %discretization step for the light model
zmin=params.geo.zrange(1);
zmax=params.geo.zrange(3);
rmax=params.geo.rmax; %mm; how far in the x-y plane to calculate heat
corrector=rem(rmax,u_step);
rmax=rmax+corrector;

glass_pixels = floor(params.geo.d_glass/u_step); %number of pixels at surface that are the coverslip
water_pixels = floor(params.geo.d_water/u_step); %number of pixels at surface that are water
skull_pixels = floor(params.geo.d_skull/u_step);
crust_pixels = floor(params.geo.d_boneCortical/u_step);

%Constants from Literature (Elwassif et al., 2006; Janssen et al., 2005)
density=params.mech.pTissue; %density of brain
spheat=params.therm.cTissue; %specific heat of brain
kbrain=params.therm.kTissue; %whole brain thermal diffusivity constant;

kglass = params.therm.kGlass; %coverslip thermal diffusivity
kwater = params.therm.kWater; %water thermal diffusivity ; convection is not simulated

wb=params.mech.wBlood; %blood perfusion rate of brain (units converted)
pblood=params.therm.cBlood; %specific heat of blood (units converted)
densblood=params.mech.pBlood; %density of blood (units converted)
%(units converted)
Tinit=37; %Temperature of brain at rest
Ta=36.7; %Arterial temperature
Tsurface = 25; %Fixed temperature at surface of water; convection is not simulated
qm=(Tinit-Ta)*wb*pblood*densblood; %value chosen to balance equation

%Time Step
deltat= 0.15* (u_step^2)/(6*(kbrain/(spheat*density))); %s; this is the
%minimum value needed for numerical stability

%Initialize vectors and matrices
depth=(zmin:u_step:zmax)';
r=0:u_step:rmax;
I=zeros(length(r),length(depth));

%CALCULATING LIGHT INTENSITY FOR EVERY X AND Z:
lightmodelrs=0:dr:rmax;
for rep=1:length(r)
    [~,in]=min(abs(r(rep)-lightmodelrs));
    I(rep,:)=frac_abs(1:u_step/dr:end,in)*Power;
end

%Initialize constants for later calculations:
uchange=(I/(density*spheat))*deltat; %change in temperature due to light
u=ones(length(r),length(depth))*Tinit; %temperature at t = 0
rg = params.geo.r_glass; % Radius of cover glass: mm

%u has dimensions radius, depth
if isfield(params, 'u_start')
    u = params.u_start;
else
    u(:,1:water_pixels) = 30; %water temperature is roughly 30c
%     u(:,1) = Tsurface;
    center = ones(ceil(rg/u_step),1)*Tsurface; %area under glass
    outside = ones(ceil(rg/u_step),1)*Tinit; %boundary zone
    edge = linspace(Tinit,Tsurface,(size(u,1)-length(center)-length(outside)))'; %remaining area
    u(:,1) = [center ; flipud(edge); outside]; %Tsurface; %IMMERSION FLUID CLAMPS SURFACE TEMP OF GLASS
end

u_init = u;

pc=1/(spheat*density); %constant for later calculations
pc_glass = 1/(params.therm.cGlass*params.mech.pGlass);
pc_boneCortical = 1/(params.therm.cBoneCortical*params.mech.pBoneCortical);
pc_boneCancellous = 1/(params.therm.cBoneCancellous*params.mech.pBoneCancellous);

pc_matrix= ones(length(r),length(depth))*pc;
pc_matrix(1:ceil(rg/u_step),water_pixels+(1:glass_pixels)) = pc_glass; 
pc_matrix((ceil(rg/u_step)+1):length(r),water_pixels+glass_pixels-(0:(skull_pixels-1))) = pc_boneCancellous;
%pc_matrix((ceil(rg/u_step)+1):length(r),[water_pixels+glass_pixels-(0:(crust_pixels-1)),water_pixels+glass_pixels-skull_pixels+(1:crust_pixels)]) = pc_boneCortical;


k=ones(length(r),length(depth))*kbrain; %This matrix can be made non-homogenous
k(1:ceil(rg/u_step),water_pixels+(1:glass_pixels)) = kglass; % Heat conductivity beyond the edge of the glass is assumed to be the same as brain.
k((ceil(rg/u_step)+1):length(r),water_pixels+glass_pixels-(0:(skull_pixels-1))) = params.therm.kBoneCancellous;
%k((ceil(rg/u_step)+1):length(r),[water_pixels+glass_pixels-(0:(crust_pixels-1)),water_pixels+glass_pixels-skull_pixels+(1:crust_pixels)]) = params.therm.kBoneCortical;


%for your applications, but numerical stability will not be assured. For
%this reason, this functionality is not built-in.

%Initalize Variable to Save Data
u_time=zeros(length(depth),length(0:deltat:t_max));
u_time(1,:)=Tinit;
t=0:deltat:t_max;
stepper=1;
ri=repmat(r',1,length(depth))+u_step/2;

%Initialize Progress Output
u_space=zeros(length(depth),length(r),length(t_save));
fprintf('\nProgress:\n')
stringout=['t=',num2str(0),' out of ',num2str(t_max),' s'];
fprintf(stringout)
tplot=0;
tsavecount=1;

%%RUN SIMULTATION%%
for t_loop=deltat:deltat:t_max;
    stepper=stepper+1;
    [m,n] = size(u);
    
    %Finite Difference Method for Heat Equation
    vrr = ((k(1:m,1:n).*(ri+u_step/2)).*u([2:m m],1:n) -(k(1:m,1:n).*(ri+u_step/2)+k([1 1:m-1],1:n).*(ri-u_step/2)).*u(1:m,1:n) +...
        (k([1 1:m-1],1:n).*(ri-u_step/2)).*u([1 1:m-1],1:n))./(u_step.^2*ri);
    vzz = ((k(1:m,1:n)).*u(1:m,[2:n n]) -(k(1:m,1:n)+k(1:m,[1 1:n-1])).*u(1:m,1:n) +...
        (k(1:m,[1 1:n-1])).*u(1:m,[1 1:n-1]))./(u_step.^2);
    
    %Discretization of Laplacian: Deltav = v_xx + v_yy + v_zz
    deltau = (vrr+vzz)*deltat.*pc_matrix; %*pc;  %.*pc_matrix;
    
    %Heat Change Due to Perfusion by Blood and Metabolic Heat
    deltaperfusion=((Ta-u)*pblood*densblood*wb+qm)*pc*deltat;
    
    %no perfusion in water or glass
    deltaperfusion(:,1:water_pixels) = 0;
    deltaperfusion(1:ceil(rg/u_step),water_pixels+(1:glass_pixels)) = 0;
    
    %Check if light is on or off
    %Total Heat Change
    
    t_ondiff=t_on-t_loop;
    t_offdiff=t_off-t_loop;
    if (isempty(t_offdiff(t_offdiff<=0)) && ~isempty(max(t_ondiff(t_ondiff<=0)))) || (~isempty(max(t_ondiff(t_ondiff<=0))) && max(t_ondiff(t_ondiff<=0))>max(t_offdiff(t_offdiff<=0)))
        uchange1=(deltau+uchange+deltaperfusion);
    else
        uchange1=(deltau+deltaperfusion);
    end
    
    %Absorbing ends
    %uchange1(end,:)=0; % comment out for insulating lateral edge
    uchange1(:,1)=0; %clamp the temperature at the surface
    uchange1(:,end)=0; %clamp the tempeature deep
    
    %Increase Temperature
    u=u+uchange1;
    
    %Save Data
    u_time(:,stepper)=(nansum(u((r+u_step/2)<=r_avg,:).*repmat(r(r+u_step/2<=r_avg)',1,size(u,2)))./(nansum(r(r+u_step/2<=r_avg),2)));  %-Tinit
    if tsavecount<=length(t_save) && t_save(tsavecount)-t_loop<=deltat && t_save(tsavecount)-t_loop>0
        u_space(:,:,tsavecount)= u';     %(u-u_init)';        %u'-Tinit;
        uRAW = u;
        tsavecount=tsavecount+1;
    end
    
    %Progress Update
    if t_loop>tplot
        removeout=repmat('\b',1,length(stringout));
        stringout=['t=',num2str(floor(t_loop)),' out of ',num2str(t_max),' s'];
        fprintf([removeout stringout])
        tplot=tplot+1; %Update Every Second
    end
    
end

%keyboard
%%END SIMULATION%%