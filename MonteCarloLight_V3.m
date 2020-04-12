function [frac_abs,frac_trans,r,depth,catcher,nlaunched, lostphotons] = MonteCarloLight_V2(varargin)

%This function outputs a matrix of light intensity values for light
%emitted out of a fiber optic, simulated in cylindrical coordinates.
%This model is based on the work of Jacques and colleagues (see Jacques,
%1998; Jacques, 2003; Jacques, 2010)

%This is a modified version based on Kaspar Podgorski's version (Podgorski
%& Ranganathan, 2016), with the following major updates: 
%1. Updated tissue optical parameters for 1280 nm and 1320 nm
%for 3PM simulation based experimental atttenuation length measurement. 
%2. Organize parameter definitions in a global struct
%in the script calls this function, instead of passing parameters as a
%long list to this function.

%Parameters required for execuation of this function:
%1. simulation wavelength: params.opt.wavelength (unit:nm)
%2. tissue absorption coefficient: params.opt.absTissue (unit:1/mm)
%3. tissue scattering coefficient: params.opt.scatterTissue (unit:1/mm)
%4. tissue scattering anisotropy coefficient: params.opt.gTissue (unit:dimensionless)
%5. tissue refractive index: params.opt.nTissue (unit:dimensionless)
%6. 1/e^2 beam radius at the objective lens back aperture: params.geo.w0 (unit:mm)
%7. tissue focal depth: params.geo.focalDepthTissue (unit:mm)
%8. objective lens numerical aperture: params.geo.NA (unit:unitless)
%9. objective lens focal length: params.geo.f (unit:mm)
%10. objective lens work distance: params.geo.wd (unit:mm)
%11. scanning field-of-view diameter: params.geo.FOV (unit:mm)

%Outputs:
%1. frac_abs - absorbed fraction, A [1/mm^3] (depth x radial distance)
%2. frac_trans - fractional transport, T [1/mm^2] (depth x radial distance)
%3. r - radial distances in mm corresponding to first dimension of
%   frac_abs and frac_trans
%4. depth - depth in mm, corresponding to second dimension of these
%matrices.
%5. catcher - the matrix recording the absorbed photons at each point in the
%simulation volume.
%6. nlaunched - the total number of photons used for the simulation
%7. lostphotons - book-keeping on the number of photons hitting the
%simulation volume boundaries or backscattered to the window.

%by Tianyu Wang 2019

% COMMENTS FROM THE PREVIOUS VERSIONS ARE LISTED BELOW:

%A slightly modified version of Joseph Stujenske's simulation code. 
%Changes by Kaspar Podgorski 2015

%This function outputs a matrix of light intensity values for light
%emitted out of a fiber optic, simulated in cylindrical coordinates.
%This model is based on the work of Jacques and colleagues (see Jacques,
%1998; Jacques, 2003; Jacques, 2010)

% The function requires the following parameters as inputs, which are now
% defined as global parameters in the scripts executing the function for
% the purpose of clarity and tractability (instead of being passed as input
% parameters):

%1. MODELNUMBER - 1 = Johansson, 2010 (default); 2 = Yaroslavsky et al.,
%2002; Data from one of these papers is used to linearly interpolate the
%scattering and absorption coefficient and the anisotropy parameter for the
%wavelength of choice modelnumber = 3 for custom input
%parameters, which is specified by the global struct params as params
%where a, s, and g are the absorption coefficient, scattering coefficient,
%and anisotropy parameter, respectively
%2. F - distance of focus from the top of the volume, unit: mm.
%3. NA - numerical aperture of fiber (between 0 and 1).
%4. WAVELENGTH - wavelength (nm) of light to use (model 1,2) -OR-
%   OPTICAL_PARAMETERS - [a s g] (model 3)
%
%Optional Inputs:
%[]=montecarlolight_function(modelnumber,...,nphotonpackets,zrange,rmax,
%   step)
%5. NPHOTONPACKETS - 100,000 * nphotonpackets are launched (requires
%integer value; rounds up if non-integer value is given). Packets are
%launched in multiples of 100,000 for computational efficiency. 100 is
%default.
%6. ZRANGE - [A B C], where A is the minimum depth value, B is the location
%of the fiber tip, and C is the maximum depth value. [0 0 6] is default. In this modified code, B
%should always be 0.
%Out of bounds photon packets are absorbed and removed from the simulation.
%7. RMAX - maximum lateral distance away from the fiber. 6 mm is default.
%8. STEP - step to use for discretization of space (mm). .01 is default.
%
%Outputs:
%1. FRAC_ABS - absorbed fraction, A [1/mm^3] (depth x radial distance)
%2. FRAC_TRANS - fractional transport, T [1/mm^2] (depth x radial distance)
%3. R - radial distances in mm corresponding to first dimension of
%frac_abs and frac_trans
%4. DEPTH - depth in mm, corresponding to second dimension of these matrices
%
%Equations:
%T = FRAC_ABS/absorption coefficient
%Fluence Rate=Power * FRAC_TRANS [mW/mm^2]
%[Note: script takes ~20 minutes to run 100*100,000 (10^7) photon packets on
%a computer with Intel Core i7-2600 CPU @ 3.4 GHz with 32 GB RAM.]
%
%Written by Joseph M. Stujenske
%Stujenske et al., Cell Reports

tic %Initialize timing
%Error checking on inputs and sorting out the variable arguments:
global params;
wavelength = params.opt.wavelength;
modelnumber = params.opt.tissueModelType;

if length(varargin)>=1
    nphotonpackets=varargin{1};
    if nphotonpackets<1
        nphotonpackets=1;
        disp('Warning: Number of Photon Packets Increased to 1*100000 (minimum)')
    end
end
if length(varargin)>=2
    zrange=varargin{2};
end
if length(varargin)>=3
    rmax=varargin{3};
end
if length(varargin)>=4
    dstep=varargin{4};
end
if ~exist('nphotonpackets','var') || isempty(nphotonpackets)
    nphotonpackets=500;
end
if ~exist('zrange','var') || isempty(zrange)
    zrange= [0 0 6];
end
if ~exist('rmax','var') || isempty(rmax)
    rmax=6;
end
if ~exist('dstep','var') || isempty(dstep)
    dr=.01; %discretization parameter for space in mm
else
    dr=dstep;
end

%Input their custom input parameters or interpolate from the models
if modelnumber==3
    absorption=params.opt.absTissue;
    scattering=params.opt.scatterTissue;
    g=params.opt.gTissue;
else
    if ((wavelength<480 || wavelength >900) && modelnumber==1) || ((wavelength<450 || wavelength>1064) && modelnumber==2)
        disp('Warning: Extrapolating Outside the Range of Experimentally Recorded Data')
    end
    %Data from Johansson, 2010 (Model 1) and Yaroslavsky, 2002 (Model 2)
    absorptionmat{1}=[480 0.37;... % unit: 1/mm; solved from scattering coefficient for in vivo human cortex
        560	0.26;...
        580	0.19;...
        640	0.05;...
        780	0.02;...
        900	0.02];
    absorptionmat{2}=[450 0.07;... unit: 1/mm; human ex vivo no blood 
        510 0.04;...
        630 0.02;...
        670 0.02;...
        1064 0.05];
    scatteringmat{1}=[480 11;... % unit:1/mm human ex vivo cortex
        580 9.7;...
        640 9.0;...
        700 8.2;...
        780 7.8;...
        900 6.6];
    scatteringmat{2}=[450 11.7;... % unit:1/mm human ex vivo no blood
        510 10.6;...
        630 9.0;...
        670 8.4;...
        1064 5.7];
    gmat{1}=[480 .89;...
        580 .89;...
        640 .89;...
        700 .90 ;...
        780 .90;...
        900 .90];
    gmat{2}=[450 .88;...
        510 .88;...
        630 .89;...
        670 .91;...
        1064 .9];
    
    %Interpolate from models
    absorption=interp1(absorptionmat{modelnumber}(:,1),...
        absorptionmat{modelnumber}(:,2),wavelength,'linear','extrap'); % absorption coefficient in mm^-1
    scattering=interp1(scatteringmat{modelnumber}(:,1),...
        scatteringmat{modelnumber}(:,2),wavelength,'linear','extrap'); %scattering coefficient in mm^-1
    g=interp1(gmat{modelnumber}(:,1),gmat{modelnumber}(:,2),wavelength,...
        'linear','extrap'); %anisotropy factor (between 0 and 1)
    g(g<0)=0;
    g(g>1)=1;
end

%Initialize variables
nphotonstotal_touse= 10000*nphotonpackets;   %100000*nphotonpackets; %Total number of photon packets to launch
nphotons=  10000; %100000; %The number of photons to have "in play" at any particular time.
lostphotons = [0 0 0 0]; % The count of photons lost to different boundaries.

% %for a lightweight simulation:
% nphotonstotal_touse= 1000*nphotonpackets;
% nphotons=  1000;

nlaunched=nphotons;
dz=dr; % discretization parameter for space in mm

n=params.opt.nTissue; % refractive index of the brain
zmin=zrange(1); zcenter=zrange(2); zmax=zrange(3); % define bounds in z

% The following two lines are only used for tophat beam profile:
% th_out=asind(NA/n); % half maximum angle of objective focus 
% r = f*tand(th_out); 

catcher=zeros(length(0:dr:rmax),length(zmin:dz:zmax)); %matrix to accept energy from photons

newplotthreshold=nphotons; %How often to give command line progress update

%Initialize Progress Update
fprintf('\nProgress:\n')
stringout=[num2str(nphotons./...
    nphotonstotal_touse),'%% of photons launched\nEstimated Time Remaining: Estimating'];
fprintf(stringout)

%Launch first batch of photons
%[cors,dirs]=photoninit(nphotons,th_out,r,zcenter);
fll = params.geo.focalDepthTissue; % focal depth inside the tissue
[cors,dirs]=photoninitGauss(nphotons,params.geo.w0,zcenter,fll);

w=ones(1,nphotons); %weight matrix which indicates the percentage of
%photon energy left in each packet

% Pseudo-code for the loop below:
% Each photon has 3 attributes: weight, position, and direction.
% Randomly generate free path length for each photon before the next
% colision
% Update photon position by moving along their dir for the randomly generated
% free path length
% Set weights of out-of-bound photons to zero, and return them to the
% center of fiber.
% Reduce in-bound photon weights by absorption, and then add lost
% energy to catcher matrix.
% Randomly generating scattered direction for each photon, and update
% directions.
% 
while 1
    %choose step size for photon packets
    s=-log(rand(1,nphotons))/(absorption+scattering);
    
    %move photon packets in the x,y, and z direction by the step size s
    cors=cors+repmat(s,3,1).*dirs;
    
    % Count out of bounds photons:
    lostphotons(1) = lostphotons(1)+sum(w(sqrt(sum(cors(1:2,:).^2))>=rmax+dr)); % Out of radial boundaries
    lostphotons(2) = lostphotons(2)+sum(w(cors(3,:)>=zmax+dz)); % too deep
    lostphotons(3) = lostphotons(3)+sum(w(cors(3,:)<zmin & sqrt(sum(cors(1:2,:).^2))<params.geo.r_glass)); %Back reflected to coverglass
    lostphotons(4) = lostphotons(4)+sum(w(cors(3,:)<zmin & sqrt(sum(cors(1:2,:).^2))>=params.geo.r_glass)); %Back reflected to skull
    
    %out of bounds photon packets are killed:
    outofbounds=sqrt(sum(cors(1:2,:).^2))>=rmax+dr | cors(3,:)>=zmax+dz | cors(3,:)<zmin;    
    w(outofbounds)=0;
    
    %put out of bounds photon packets "inbounds" to avoid errors later on
    cors(:,w==0)=[zeros(2,sum(w==0));ones(1,sum(w==0))*zcenter]; % Returned to the center of the fiber, but their weights remain zero.
    rcors=sqrt(sum(cors(1:2,:).^2));
    
    %photon packet weight is stored in the catcher matrix:
    zs=floor((cors(3,:)-zmin)/dz)+1; % How many dz from the top for each photon
    rs=floor(rcors/dr)+1; % How many dr from the center
    ins=sub2ind(size(catcher),rs,zs); % Returns the position each photon corresponds to in the catcher matrix, calculated by radial and vertical position.
    wacumm = accumarray(ins',w); % wacumm has the same size as catcher, but linearized. Each element of wacumm is the sum of the weights of the photons in the corresponding position in catcher.
    ns=(wacumm~=0);
    catcher(ns)=catcher(ns)+wacumm(ns)*(absorption/(absorption+scattering));
    
    %reduce weight of the photon packets by what was absorbed:
    w=w*(scattering/(absorption+scattering));
    
    %kill packets that are within rmax and rmax+dr and zmax and zmax+dz
    lostphotons(1) = lostphotons(1)+sum(w(rcors>=rmax)); % Out of radial boundaries
    lostphotons(2) = lostphotons(2)+sum(w(cors(3,:)>=zmax)); % too deep
    outofbounds2=(rcors)>=rmax | cors(3,:)>=zmax;
    w(outofbounds2)=0;
    
    %change direction of the photon packets:
    if g>0 && g<=1
        costh=((1+g.^2-((1-g.^2)./(1-g+2*g.*rand(1,nphotons))).^2)./(2.*g)); %Heyey Greenstein scattering
    elseif g==0
        costh=2*rand(1,nphotons)-1;
    else
        error('g must be between 0 and 1');
    end
    phi=2*pi*rand(1,nphotons);
    sinth=sqrt(1-costh.^2);
    temp=sqrt(1-dirs(3,:).^2);
    
    uxyz=repmat(sinth,3,1).*[(dirs(1:2,:).*repmat(dirs(3,:).*cos(phi),2,1)+...
        [-dirs(2,:);dirs(1,:)].*repmat(sin(phi),2,1))./repmat(temp,2,1);-cos(phi).*temp]+...
        dirs.*repmat(costh,3,1);
    %Above line is equivalent to:
    %             uxx=sinth.*(xdir.*zdir.*cos(phi)-ydir.*sin(phi))./temp+xdir.*costh;
    %             uyy=sinth.*(ydir.*zdir.*cos(phi)+xdir.*sin(phi))./temp+ydir.*costh;
    %             uzz=-sinth.*cos(phi).*temp+zdir.*costh;
    
    %photon packets very close to being vertical are treated as being
    %vertical:
    tofix=abs(dirs(1,:))<.0001 & abs(dirs(2,:))<.0001;
    if any(tofix)
        uxyz(:,tofix)=[repmat(sinth(tofix),2,1).*[cos(phi(tofix));...
            sin(phi(tofix))];costh(tofix).*sign(dirs(3,tofix))];
        %Above line is equilavent to:
        %             uxx(tofix)=sinth(tofix).*cos(phi(tofix));
        %             uyy(tofix)=sinth(tofix).*sin(phi(tofix));
        %             uzz(tofix)=costh(tofix).*sign(zdir(tofix));
    end
    dirs=uxyz;
    
    %correct dirs to be a unit vector if it has deviated
    %due to rounding errors:
    mag=sqrt(sum(dirs.^2));
    dirs=dirs./repmat(mag,3,1);
    
    %roulette to see if packets with low weight die while preserving
    %conservation of energy:
    chance=rand(1,nphotons);
    w(w<1e-4 &chance<=.1& w>0)=w(w<1e-4 &chance<=.1& w>0)./.1; % 90% will be destoried, and the other 10% gain the energy lost.
    todestroy=(w<1e-4 & chance>=.1 & w>0) | outofbounds | outofbounds2;
    ntodestroy=sum(todestroy);
    
    %replace destroyed photon packets with new photon packets
    if ntodestroy>0
        if ntodestroy+nlaunched<=nphotonstotal_touse 
        %    [cors(:,todestroy),dirs(:,todestroy)]=photoninit(ntodestroy,th_out,r,zcenter);
            [cors(:,todestroy),dirs(:,todestroy)]=photoninitGauss(ntodestroy,params.geo.w0,zcenter,fll);
            w(todestroy)=1;
            nlaunched=nlaunched+ntodestroy;
        elseif nlaunched<nphotonstotal_touse 
            which=find(todestroy);
            replaceins=(which(1:nphotonstotal_touse-nlaunched));
            %[cors(:,replaceins),dirs(:,replaceins)]=photoninit(length(replaceins),th_out,r,zcenter);
            [cors(:,replaceins),dirs(:,replaceins)]=photoninitGauss(length(replaceins),params.geo.w0,zcenter,fll);
            w(replaceins)=1;
            nlaunched=nlaunched+length(replaceins);
            w(which(nphotonstotal_touse-nlaunched+1:end))=0;
        else
            w(w<1e-4 &chance>=.1 & w>0)=0;
        end
    end
    
    %terminate the simulation if all of the packets have no weight:
    if all(w==0)
        break
    end
    
    %Progress output
    if nlaunched>=max(newplotthreshold,nphotons+1)
        tout=toc;
        
        %Estimate Time to Destroy the Rest of the Photon Packets
        timeremaining=tout.*((nphotonstotal_touse/(nlaunched-nphotons))-1); 
        
        newplotthreshold=newplotthreshold+nphotons;
        removeout=repmat('\b',1,length(stringout)-2);
        stringout=[num2str(nlaunched./nphotonstotal_touse*100),...
            '%% of photons launched\nEstimated Time Remaining: ',...
            num2str(timeremaining),' s'];
        fprintf([removeout stringout])
    end
end

%Clean Up - Generate Outputs
totaltime=toc;
fprintf(['\nSimulation Time: ',num2str(totaltime),' s','\n'])

%Save data!
frac_abs=(catcher./repmat(2*pi*((1:size(catcher,1))-.5)',1,size(catcher,2)))';
frac_abs=frac_abs./(dr.^2*dz)./(nphotonstotal_touse); %Correct for cylindrical geometry and make fraction
frac_trans=frac_abs/absorption;

params.geo.dstep=dr;
params.geo.zrange=zrange;
params.geo.rmax=rmax;
params.opt.absTissue=absorption; % Save parameters if they were obtained by interpretation.
params.opt.scatterTissue=scattering;
params.opt.gTissue=g;

r=0:dr:rmax;
depth=zmin:dz:zmax;


% This function is used for objective lens with an over-filling beam.
% photon initialization with an under-filling beam is implemented in photoninitGauss.m
% Tianyu Wang 11/13/2019

function [cors,dirs]=photoninit(num,th_out,r,zcenter)
% num: number of photons, e.g. number of photons in one photon packet which
% is set to 10000 by nphotons.
% th_out: half focusing angle, defined as arcsin(NA/n).
% zcenter: the z location where the light ray simulation starts
% r: the radius of the light cone at zcenter, 

fl = r/tand(th_out);

rinit = sqrt(rand(1,num))*r; %radial location %sqrt for uniform in 2D
thinitpos = -180+2*rand(1,num)*180; %initial position (angular coordinate)

thinit = fl./sqrt(fl^2 + rinit.^2);
phiinit = -thinitpos;

%directions are encoded in spherical coordinates:
% thinit = cosd(th_out)+rand(1,num)*(1-cosd(th_out)); %length of the Z-vector
% phiinit=-180+2*rand(1,num)*180; %inital angle in the z plane
%initialize x,y, and z direction matrices such that [x y z] is a unit
%vector
xdir=sqrt(1-thinit.^2).*cosd(phiinit);
ydir=sqrt(1-thinit.^2).*sind(phiinit);
zdir=thinit;

%simulating scanning by translating the fiber in xy @ zcenter, NOT angular
%scanning. Give each photon a random translation in xy.
FOV = params.geo.FOV;
dx = (FOV*rand(1,num)-0.5*FOV);
dy = (FOV*rand(1,num)-0.5*FOV);

%initializing x, y, and z coordinates at zcenter
xcor=rinit.*cosd(thinitpos);
ycor=rinit.*sind(thinitpos);
zcor=ones(1,num)*zcenter;
cors=[xcor+dx;ycor+dy;zcor];
dirs=[xdir;ydir;zdir];

end

end