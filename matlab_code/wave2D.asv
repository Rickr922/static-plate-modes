clear all; 
close all;

boundaryCond = 0;   %0 = Dirichlet; 1 = Neumann
init = 0;           %sets if to initialize the states or not

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Excitation Signal
%'impulse', 'sin', 'square', 'sweep', 'audiofile'
inputType = 'impulse'; 
%in case inputType is 'audiofile', specify file name and path
audiofileName = 'dry_samples/DrumReference.wav';
%amplification factor
osFac = 1;
%if durSec is set to 0 when audiofileName is 'audiofile', the entire file
%is played
durSec = 2;
[excit,SR,k,timeSamples,timeVec,durSec] = ExcitSignal(1,osFac,durSec,inputType);

% Define variables of the system
Lx = 1;
Ly = 2;
Lz = 9e-4;                   %[m] Thickness

rho = 1400; %Mylar density in kg/m^3

rhoH = rho*Lz;

T = 3000; %1000;                      %[N] Tension
c = T/(rho*Lz);

% Calculate grid spacing from variables
h = c * k*sqrt(2);
Nx = floor(Lx/h);
Ny = floor(Ly/h);

% Nx = 70;
% Ny = 70;

lambda = c*k/h;

inPoint = [0.52*Lx,0.53*Ly];
outPoint = [0.47*Lx,0.62*Ly];

lo = floor(inPoint(1)/h);
mo = floor(inPoint(2)/h);
alphax = inPoint(1)/h - lo;
alphay = inPoint(2)/h - mo;

Jcoeff = zeros(Nx-1,Ny-1);
Jcoeff(lo,mo) = (1-alphax) * (1-alphay) / h^2;
Jcoeff(lo, mo+1) = (1-alphax)*alphay / h^2;
Jcoeff(lo+1, mo) = alphax*(1-alphay) / h^2;
Jcoeff(lo+1,mo+1) = alphax*alphay / h^2;

loO = floor(outPoint(1)/h);
moO = floor(outPoint(2)/h);
alphaxO = outPoint(1)/h - lo;
alphayO = outPoint(2)/h - mo;

JcoeffO = zeros(Nx-1,Ny-1);
JcoeffO(loO,moO) = (1-alphaxO) * (1-alphayO);
JcoeffO(loO, moO+1) = (1-alphaxO)*alphayO;
JcoeffO(loO+1, moO) = alphaxO*(1-alphayO);
JcoeffO(loO+1,moO+1) = alphaxO*alphayO;

uNext = zeros(Nx,Ny);
u = zeros(Nx,Ny);
if init
    width = 10;
    u((1:width) + floor(Nx/2)-(width/2),(1:width) + floor(Ny/2)-(width/2)) = window2(width,width,'hann');
    %u(floor(Nx/2),floor(Ny/2)) = 1;
    excit = zeros(timeSamples,1);
end
uPrev = u;
out = zeros(durSec,1);

for n=1:timeSamples    
    exc = excit(n);
    if boundaryCond == 0
%     %Full Dirichlet
        for l=2:Nx-1
            for m = 2:Ny-1
                uNext(l,m) = 2*(1-2*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1)) + Jcoeff(l,m)*exc;
            end
        end
    elseif boundaryCond == 1
        %Full Neumann
        for l=1:Nx
            if l==1
                uNext(l,1) = 2*(1-lambda^2)*u(l,1) - uPrev(l,1) + lambda^2*(u(l+1,1) + u(l,2)) + Jcoeff(l,m)*exc;
                for m = 2:(Ny-1) %taking corners into account
                    uNext(l,m) = (2-3*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l+1,m) + u(l,m+1) + u(l,m-1))  + Jcoeff(l,m)*exc;
                end
                uNext(l,Ny) = 2*(1-lambda^2)*u(l,Ny) - uPrev(l,Ny) + lambda^2*(u(l+1,Ny) + u(l,Ny-1)) + Jcoeff(l,m)*exc;
            elseif l==Nx
                uNext(l,1) = 2*(1-lambda^2)*u(l,1) - uPrev(l,1) + lambda^2*(u(l-1,1) + u(l,2)) + Jcoeff(l,m)*exc;
                for m = 2:(Ny-1) %taking corners into account
                    uNext(l,m) = (2-3*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l-1,m) + u(l,m+1) + u(l,m-1)) + Jcoeff(l,m)*exc;
                end
                uNext(l,Ny) = 2*(1-lambda^2)*u(l,Ny) - uPrev(l,Ny) + lambda^2*(u(l-1,Ny) + u(l,Ny-1)) + Jcoeff(l,m)*exc;
            else
                for m = 1:Ny
                    if m==1
                        uNext(l,m) = (2-3*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l+1,m) + u(l-1,m) + u(l,m+1)) + Jcoeff(l,m)*exc;
                    elseif m==Ny
                        uNext(l,m) = (2-3*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l+1,m) + u(l-1,m) + u(l,m-1)) + Jcoeff(l,m)*exc;
                    else
                        uNext(l,m) = 2*(1-2*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1)) + Jcoeff(l,m)*exc;
                    end
                end
            end
        end  
    end
    
%     if n==5
%         u((1:width) + floor(Nx/2)-(width/2),(1:width) + floor(Ny/2)-(width/2)) = window2(width,width,'hann');
%     end
    
    %osservo solo UN punto della "corda"
    out(n)= uNext(floor(outPoint(1)/h),floor(outPoint(2)/h));
    
    surf(uNext)
    %plot(u);
    zlim([-10,10]);
    view([45 45]);
    %xlim([50,150]);
    drawnow;
    pause(50/1000)
    
    uPrev = u;
    u = uNext;
end

soundsc(out,SR)


figure(1)
plot(out)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%Functions
function [excit,SR,k,timeSamples,timeVec,durSec]=ExcitSignal(amp,OSFac,durSec,excitType,filename)
    switch nargin
        case 4
            if ~durSec
                disp('Zero input duration');
                return
            end
            SR = OSFac*44100;
            timeSamples = durSec*SR;
            k = 1/SR;
            timeVec = (1:timeSamples)*k;
            if excitType == "impulse"
                excit = zeros(timeSamples,1);
                excit(5)=amp*1;
            elseif excitType == "sin"
                excit(1:timeSamples) = amp*sin(200*2*pi*timeVec);
            elseif excitType == "square"
                excit(1:timeSamples) = amp*square(200*2*pi*timeVec);
            elseif excitType == "sweep"
                startFreq = 200;
                endFreq = 5000;
                chirpLengthSec = durSec;%floor(2*durSec/3);
                chirpLength = timeVec(1:SR*chirpLengthSec);
                
                envelope = [linspace(0,amp,SR*chirpLengthSec/10), amp*ones(1,SR*8*chirpLengthSec/10), linspace(amp,0,SR*chirpLengthSec/10)];
                
                forcing = envelope.*chirp(chirpLength,startFreq,chirpLengthSec,endFreq,'quadratic',90);
                excit = zeros(1,timeSamples);
                excit(1:length(forcing)) = forcing;
            else 
                disp('Wrong input type');
                return;
            end
        case 5
            if excitType == "audiofile"
                [excit,SR] = audioread(filename);
                if OSFac > 1
                    excit = resample(excit,OSFac,1);
                    SR = SR*OSFac;
                end
                k = 1/SR;
                if durSec
                    excit = amp*excit(1:floor(SR*durSec),1);
                else
                    excit = amp*excit;
                end
                timeSamples = length(excit);
                timeVec = (1:timeSamples)*k;
                durSec = timeSamples/SR;
            end
        otherwise
            disp("wrong input type");
            return
    end
end