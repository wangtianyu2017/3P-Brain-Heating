% Calculation of the temperature profile in the mouse neocortex beneath a
% cranial window and immersion water under 920-nm 2-photon imaging with
% different scanning depth and power. This code was used to produce part of
% Figure 3 figure supplement 3.
% The value of 'depthArray' indicates the array of scanning depth;
% 'powerPercent' indicates the array of power after the objective lens. The
% actual power after the objective lens is calculated as
% (12mW)*exp(depthArray(i)*utot). utot is the attenuation coefficient in
% the unit of mm.

% Tianyu Wang 11/13/2019

% Set up all model parameters in a global struture
clear all;
global params;

% Optical parameters and properties of media
params.opt.wavelength = 920; % nm
params.opt.plotWavelength = []; % nm, plot wavelength may be different from the actual wavelength, if it is invisible.
params.opt.tissueModelType = 3; %Data from Johansson, 2010 (Model 1) and Yaroslavsky, 2002 (Model 2)
params.opt.absTissue = 0.039;  % 1/mm. If 'nan', will be calculated and assigned with a model in MonteCarloLight_V2
params.opt.scatterTissue = 6.66; % 1/mm, our own measurement
params.opt.gTissue = 0.9; % anistropy factor
params.opt.nTissue = 1.36; % Refractive index of brain
params.opt.absGlass = [];
params.opt.nGlass = [];
params.opt.absWater = [];
params.opt.nWater = 1.328;
params.opt.scatterFrontCort = @(x) 10.9*((x/500).^(-0.334)); % Human data, x normalized to 500nm, Bevilacqua et al 2000
params.opt.scatterTempCort = @(x) 11.6*((x/500).^(-0.601)); % Human data, x normalized to 500nm, Bevilacqua et al 2000

% Mechanical properties of media
params.mech.pTissue = 1.04E-6; % density: kg/mm^3
params.mech.pBlood = 1.06E-6; % density: kg/mm^3
params.mech.wBlood = 8.5E-3; % blood perfusion rate: /s 
params.mech.pGlass = 2.23E-6; % density: kg/mm^3
params.mech.pWater = 1E-6; % density: kg/mm^3
params.mech.pBoneCancellous = 1.178E-6; % density: kg/mm^3
params.mech.pBoneCortical = 1.908E-6; % density: kg/mm^3

% Thermal properties of media
params.therm.cTissue = 3.65E6; % specific heat: mJ/kg/C
params.therm.kTissue = 0.527; % heat conductivity: mW/m/C
params.therm.cBlood = 3.6E6; % specific heat: mJ/kg/C
params.therm.qBrain = 9.7E-3; % Brain metabolic heat mW/mm^3
params.therm.cGlass = 0.647E6; % specific heat: mJ/kg/C
params.therm.kGlass = 0.8; % heat conductivity: mW/m/C
params.therm.cWater = 4.184E6; % specific heat: mJ/kg/C
params.therm.kWater = 0.6; % heat conductivity: mW/m/C
params.therm.cBoneCancellous = 2.274E6; % specific heat: mJ/kg/C
params.therm.kBoneCancellous = 0.31; % heat conductivity: mW/m/C
params.therm.cBoneCortical = 1.313E6; % specific heat: mJ/kg/C
params.therm.kBoneCortical = 0.32; % heat conductivity: mW/m/C
params.therm.uStep = 0.03; % minimum length unit for heat conduction simulation: mm

% Parameters for sample geometry
params.geo.dstep=[];
params.geo.zrange=[];
params.geo.rmax=[];
params.geo.d_glass = 0.16; % Thinkness of cover glass: mm
params.geo.r_glass = 2; % Radius of cover glass: mm
params.geo.d_water = []; % Thickness of immersion water: mm
params.geo.d_skull = 0.14; % Thinkness of skull: mm
params.geo.d_boneCortical = 0.01; % Thinkness of skull crust: mm

% Parameters for focusing and scanning geometry
params.geo.NA = 1.05;
params.geo.f = 7.2; % Focal length of the objective: mm
params.geo.wd = 2; % Work distance of the objective: mm
params.geo.w0 = 5.3; % 1/e^2 radius of the beam at objective back aperture, usually matched to 70% of back aperture diameter for deep imaging. unit: mm
params.geo.FOV = 0.3; % Linear FOV of focal scan: mm
params.geo.focalDepthTissue = 0; % Focal depth in the tissue: mm 
params.geo.surfIllumRadius = []; % The raius of the instersection circle of light cone and brain tissue

% Simulate light propagation for different focal depth
utot = params.opt.scatterTissue+params.opt.absTissue;
depthArray = [0:7]./utot;
depthArray(1) = 0.01; % avoid computation error at zero depth
powerPercent = [0.02 0.04 0.06 0.08 0.1 0.15 0.2 0.25 0.3 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2.0];
maxTemp = zeros(length(powerPercent),length(depthArray));
absPortion = zeros(size(depthArray));
lostPortions = zeros(4,length(depthArray));

mkdir ./920GaussDepthScan;
mkdir ./920GaussDepthScan/TempMaps;
mkdir ./920GaussDepthScan/DeltaTempMaps;

for i = 1:length(depthArray)    
    params.geo.focalDepthTissue = depthArray(i);
    [frac_abs,frac_trans,r1,d1,catcher, nlaunched, lostphotons] = MonteCarloLight_V3; 
    params.opt.plotWavelength = [];
    fig = figure; 
    LightHeatPlotter_V2(r1,d1,frac_trans./max(frac_trans(:)),[0.01 0.1 0.5]); %wavelength for plotting set to 630 so it looks red, not relevant to simulation
    drawnow
    disp([num2str(i),' scattering length=', num2str(1000*i/utot),'um; abs_frac=',num2str(sum(sum(catcher))/5000000)]);
    absPortion(i) = sum(sum(catcher))/nlaunched;
    lostPortions(:,i) = lostphotons'/nlaunched;
    saveas(fig,['./920GaussDepthScan/light_920nm_',num2str(i-1),'attleng.fig']);
    close(fig);
    save(['./920GaussDepthScan/light_920nm_',num2str(i-1),'attleng.mat'],'frac_abs','frac_trans','r1','d1','catcher');
    
    %simulation parameters for HeatDiffusionLight_V2
    params.geo.d_water = ceil((params.geo.wd - params.geo.d_glass - depthArray(i))/params.geo.dstep)*params.geo.dstep;
    added_Z = (params.geo.d_water+params.geo.d_glass);
    params.geo.zrange(1) = -added_Z;
    times = 60; % s
    maxPower = 12*exp(depthArray(i)*utot); % 12 means 12mW at focus, which is 0.15nJ at 80MHz rep rate
    powers = powerPercent*maxPower; % mW

    %extend the simulation volume to the glass and water above the brain (not including this volume in the call to MonteCarloLight_V2 above models them as transparent)
    frac_abs_padded = [zeros(ceil(added_Z/params.geo.dstep), size(frac_abs,2)) ; frac_abs];

    %extize
    u_time = cell(1,length(powers));
    u_space = cell(1,length(powers));
    t = cell(1,length(powers));
    r2 = cell(1,length(powers));
    depth = cell(1,length(powers));

    %Initialization: run the model long enough to allow it to reach steady state at the clamped initial conditions
    if isfield(params,'u_start')
        params = rmfield(params,'u_start');
    end
    [O1,O2,O3,O4,O5, u_start] = HeatDiffusionLight_V2(frac_abs_padded,.03,0,60,60,0,60,.1); 
    params.u_start = u_start;


    %Simulate heating at the various power levels
    for power = 1:length(powers) 
        if powers(power) > 500
            continue;
        end
        [u_time{power},u_space{power},t{power},r2{power},depth{power}]=HeatDiffusionLight_V2(frac_abs_padded,.03,0,max(times),max(times),powers(power),times,.1);
        u_time{power} = u_time{power}(:,2:((size(u_time{power},2)-1)/2000):end); %cut down this gigantic matrix
        maxTemp(power,i) = prctile(reshape(u_space{power}(:,:,1),1,numel(u_space{power})),99); %99 percentil as maximum temperature
    end
    
    save(['./920GaussDepthScan/heat_920nm_',num2str(i-1),'attleng.mat'],'u_time','u_space','t','r2','depth');
     
    %plot temperature at different powers and times
    params.opt.plotWavelength = 'hot';
    for power = 1:length(powers)
        if powers(power) > 500
            continue;
        end
        for time = 1:length(times)
            figName = ['920nm_',num2str(i-1),'le_',int2str(powers(power)) 'mW_' int2str(times(time)) 's'];
            fig = figure('Name', figName);
            LightHeatPlotter_V2(r2{power},depth{power},u_space{power}(:,:,time),[25:1:100]); caxis([30 42]);
            set(gca, 'ylim', [-added_Z 6], 'xlim', [-6 6])
            saveas(fig,['./920GaussDepthScan/TempMaps/',figName]);
            close(fig);
        end
    end

    %plot deltaT; i.e. difference in temperature between 0 and specified power
    for power = 1:length(powers)
        if powers(power) > 500
            continue;
        end
        for time = 1:length(times)
            figName = ['deltaT_920nm_',num2str(i-1),'le_',int2str(powers(power)) 'mW_' int2str(times(time)) 's'];
            fig = figure('Name', figName);
            LightHeatPlotter_V2(r2{power},depth{power},u_space{power}(:,:,time) - u_start',[0:1:70]);
            set(gca, 'ylim', [-added_Z 6], 'xlim', [-6 6])
            saveas(fig,['./920GaussDepthScan/DeltaTempMaps/',figName]);
            close(fig);
        end
    end
end

save('./920GaussDepthScan/params_920nm.mat','params','maxTemp','absPortion','lostPortions');