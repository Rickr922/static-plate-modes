%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   Modal plate simulation with automatic c++ header output
%                    Riccardo Russo
%                 University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Excitation Signal
%'impulse', 'sin', 'square', 'sweep', 'audiofile'
inputType = 'impulse'; 
%in case inputType is 'audiofile', specify file name and path
audiofileName = 'DrumReference.wav';
%amplification factor
osFac = 1;
%if durSec is set to 0 when audiofileName is 'audiofile', the entire file
%is played
durSec = 2;
%[excit,SR,k,timeSamples,timeVec,durSec] = ExcitSignal(1,osFac,durSec,inputType,audiofileName);
[excit,SR,k,timeSamples,timeVec,durSec] = ExcitSignal(1,osFac,durSec,inputType);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
play = true;
playDry = false;
saveAudio = false;

exactOsc = true;

%Normalize output (before visualization)
normalizeOut = true;

smSolver = false;

% Lx = 0.4;                       %[m] Hor length
% Ly = 0.8;                       %[m] Ver lentgh
Lx = 3.5;                       %[m] Hor length
Ly = 3;                       %[m] Ver lentgh
Lz = 5e-4;                    %[m] Thickness

materialData = GetMaterialData("gold");

E = materialData(1);
ni = materialData(2);
rho = materialData(3);

T = 600;                      %[N] Tension
D = E*Lz^3/(12*(1-ni^2));     %Flexural Rigidity

omegaLim = 20000*2*pi;

inPoint = [0.52*Lx,0.53*Ly];
outPoint1 = [0.47*Lx,0.62*Ly];

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Computing eigenfrequencies and eigenvectors

% Computing max number of modes
Nx = 0; Ny = 0; omegaCurr = 0;
while omegaCurr<omegaLim
    Nx = Nx + 1; 
    omegaCurr = ComputeOmega(Nx,0,T,rho,D,Lx,Ly,Lz);
end
omegaCurr = 0;
while omegaCurr<omegaLim
    Ny = Ny + 1;
    omegaCurr = ComputeOmega(0,Ny,T,rho,D,Lx,Ly,Lz);
end

% Computing eigenfrequencies
eigenFreqs = zeros(Nx*Ny,3);

for i=1:Nx
    for j=1:Ny
        omegaCurr = ComputeOmega(i,j,T,rho,D,Lx,Ly,Lz);
        if omegaCurr<omegaLim && omegaCurr>(20*2*pi)
            
            eigenFreqs(j+(i-1)*Ny,1)=omegaCurr;
        
            eigenFreqs(j+(i-1)*Ny,2) = i;
            eigenFreqs(j+(i-1)*Ny,3) = j;
        end
    end
end
%removing invalid frequencies
eigenFreqs(~any(eigenFreqs,2),:) = [];
eigenFreqs = sortrows(eigenFreqs,1);
modesNumber = size(eigenFreqs,1);

%Computing Modes for in and out points
modesIn = zeros(modesNumber,1);
modesOut = zeros(modesNumber,1);

for i=1:modesNumber
    modesIn(i) = ComputeMode(inPoint(1),inPoint(2),eigenFreqs(i,2),...
        eigenFreqs(i,3),Lx,Ly);
    modesOut(i) = ComputeMode(outPoint1(1),outPoint1(2),eigenFreqs(i,2),...
        eigenFreqs(i,3),Lx,Ly);
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Damping Coefficients
decaySec = [8,7,8,6,5,6,3,2];
decayCoeffs = 6*log(10)./decaySec;
freqBands = [62,125,250,500,1000,2000,4000,omegaLim/2/pi];
decayCoeffsFull = zeros(modesNumber,1);

for i=1:modesNumber
    currFreqHz = eigenFreqs(i,1)/2/pi;
    count = 1;
    while count <= length(freqBands)
        if currFreqHz <= freqBands(count)
            decayCoeffsFull(i) = decayCoeffs(count);
            break;
        end
        count = count + 1;
    end 
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%% C++ Header output
fileID = fopen('../xcodeproj/plate-modes/plateModalData.h', 'wt');

if fileID < 3
    disp('error opening file');
end

fprintf(fileID, '#ifndef plateModalData\n#define plateModalData\n\n');

fprintf(fileID, '#define modesNumberFull %i\n\n', modesNumber);

fprintf(fileID, 'static float eigenFreqs[%i] = {', modesNumber);
for i = 1:modesNumber
    if i < modesNumber
        fprintf(fileID,'%f,\n',eigenFreqs(i,1));
    elseif i == modesNumber
        fprintf(fileID,'%f',eigenFreqs(i,1));
    end
end
fprintf(fileID, '};\n\n');

fprintf(fileID, 'static float modesIn[%i] = {', modesNumber);
for i = 1:modesNumber
    if i < modesNumber
        fprintf(fileID,'%f,\n',modesIn(i)/rho/Lz);
    elseif i == modesNumber
        fprintf(fileID,'%f',modesIn(i)/rho/Lz);
    end
end
fprintf(fileID, '};\n\n');

fprintf(fileID, 'static float modesOut[%i] = {', modesNumber);
for i = 1:modesNumber
    if i < modesNumber
        fprintf(fileID,'%f,\n',modesOut(i));
    elseif i == modesNumber
        fprintf(fileID,'%f',modesOut(i));
    end
end
fprintf(fileID, '};\n\n');

fprintf(fileID, 'static float dampCoeffs[%i] = {', modesNumber);
for i = 1:modesNumber
    if i < modesNumber
        fprintf(fileID,'%f,\n',decayCoeffsFull(i));
    elseif i == modesNumber
        fprintf(fileID,'%f',decayCoeffsFull(i));
    end
end
fprintf(fileID, '};\n\n');

fprintf(fileID, '#endif');

status = fclose(fileID);

%% Simulation
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
if ~exactOsc
    for i=1:modesNumber
        if eigenFreqs(i)>(2/k)
            modesNumber = i-1;
            break;
        end
    end
end

%%%%% Initializing vectors
sNext = zeros(modesNumber,1);

x = zeros(modesNumber,1);
xNext = zeros(modesNumber,1);
xPrev = zeros(modesNumber,1);

if exactOsc
    c1 = 2*exp(-decayCoeffsFull*k).*cos(sqrt(eigenFreqs(:,1).^2 - decayCoeffsFull.^2)*k);
    c2 = -exp(-2*decayCoeffsFull*k);

    c1exc = (1+k*decayCoeffsFull)/12;
    c2exc =  - 5/6 - decayCoeffsFull*k + (2/3)*k^2*decayCoeffsFull.^2 - k^2*eigenFreqs(:,1).^2/12;
    c3exc = (1-k*decayCoeffsFull)/12;

    c3 = k^2*modesIn/(rho*Lz);
    

    excit = [0;0;excit];
else
    c1 = (2-eigenFreqs(:,1).^2*k^2)./(decayCoeffsFull*k+1);
    c2 = (decayCoeffsFull*k-1)./(decayCoeffsFull*k+1); 
    c3 = k^2*modesIn./((decayCoeffsFull*k+1)*rho*Lz);
end

output = zeros(1,timeSamples);
tic
for n = 1:timeSamples

    if exactOsc
        for i = 1:modesNumber
            exc = excit(n+2)*c1exc(i) + excit(n+1)*c2exc(i) + excit(n)*c3exc(i);
            xNext(i) = c1(i)*x(i) + c2(i)*xPrev(i) + c3(i)*exc;
    
            output(n) = output(n) + xNext(i)*modesOut(i);
        end
    else
        exc = excit(n);
        for i = 1:modesNumber
            xNext(i) = c1(i)*x(i) + c2(i)*xPrev(i) + c3(i)*exc;
    
            output(n) = output(n) + xNext(i)*modesOut(i);
        end
    end

    xPrev = x;
    x = xNext;
end
realTimeFrac = toc/durSec

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Plotting & Saving Audio 
figure(1)
plot(output);

if normalizeOut
    output = output/max(abs(output));
end

if ~strcmp(inputType,'audiofile')
    maxFreq = omegaLim/2/pi; %Hz
    fontSize = 15;
    
    figure(2)
    windowLength = 1024;
    [s,freqs,t,p] = spectrogram(output,blackmanharris(windowLength),floor(windowLength/2),windowLength,SR);
    colormap hot
    mesh([0,t(1:end-1)],freqs,20*log10(p));
    view(2)
    ylim([20,20000]);
    xlim([0,t(end)]);
    yticks([20 10000 20000])
    xticks([0 1 t(end)])
    set(gca,'FontSize',fontSize);
    ylabel("Frequency (Hz)");
    xlabel("Time (s)");
end

if play
    if playDry
        soundsc(excit,SR);
        pause(durSec)
    end
    diffOut = diff(output);
    soundsc(diffOut,SR);
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%Functions
function o = ComputeOmega(m1,m2,T,rho,D,Lx,Ly,Lz)
    C1 = T/rho/Lz;
    C2 = D/rho/Lz;
    C3 = (m1^2*pi^2)/Lx^2 + (m2^2*pi^2)/Ly^2;
    o = sqrt(C1*C3+C2*C3^2);
end
function m = ComputeMode(xp,yp,m1,m2,Lx,Ly)
    m = sqrt(4/Lx/Ly)*sin(m1*pi*xp/Lx)*sin(m2*pi*yp/Ly); 
end
function d = ComputeDamp(t60) 
    d = 6*log(10)./t60;
end
function g = GetMaterialData(material)
%Lz = 0.5mm
    switch material
        case "steel1"
            E = 2e11;
            ni = 0.3;   
            rho = 8.05e3;
            R1 = 9.416e-3; 
            C1 = 0.14965e-3;
        case "steel2"
            E = 2e11;
            ni = 0.3;   
            rho = 7.872e3;
            R1 = 9.664e-3; 
            C1 = 0.1855e-3;
        case "gold"
            E = 79e9;
            ni = 0.42;   
            rho = 7.872e3;
            R1 = 4.727e-3; 
            C1 = 1.270e-3;
        case "silver"
            E = 1e11;
            ni = 0.38;   
            rho = 10.49e3;
            R1 = 8.403e-3; 
            C1 = 1.679e-3;
        case "copper"
            E = 1.4e11;
            ni = 0.35;   
            rho = 8.96e3;
            R1 = 5.691e-3; 
            C1 = 1.148e-3;
        case "aluminium"
            E = 8e10;
            ni = 0.34;   
            rho = 2.7e3;
            R1 = 9.975e-3; 
            C1 = 0.976e-3;
        otherwise
            disp("wrong material");
            return;
    end
    g = [E,ni,rho,R1,C1];
end

function g=GetFreqBands(omegaVect,bandsType)
    if bandsType == "bark"
        g = [0;100;200;300;400;510;630;770;920;1080;1270;1480;1720;...
            2000;2320;2700;3150;3700;4400;5300;6400;7700;...
            9500;12000;omegaVect(end)/2/pi];
    elseif bandsType == "linear"
        g = (0:100:omegaVect(end)/2/pi)';
    end
end

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
                excit = excit(:,1);
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