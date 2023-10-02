%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%               FDTD plate simulation
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
audiofileName = 'dry_samples/DrumReference.wav';
%amplification factor
osFac = 1;
%if durSec is set to 0 when audiofileName is 'audiofile', the entire file
%is played
durSec = 2;
% [excit,SR,k,timeSamples,timeVec,durSec] = ExcitSignal(amp,osFac,durSec,inputType,audiofileName);
[excit,SR,k,timeSamples,timeVec,durSec] = ExcitSignal(1,osFac,durSec,inputType);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
play = true;
playDry = false;
saveAudio = false;

%Normalize output (before visualization)
normalizeOut = false;

smSolver = false;

Lx = 0.05;                       %[m] Hor length
Ly = 0.1;                       %[m] Ver lentgh
Lz = 5e-4;                    %[m] Thickness

sigma0 = 2;
sigma1 = 0.5;

materialData = GetMaterialData("steel2");

E = materialData(1);
ni = materialData(2);
rho = materialData(3);

T = 600;                      %[N] Tension
D = E*Lz^3/(12*(1-ni^2));     %Flexural Rigidity

kappa = sqrt(D/(rho*Lz));
c = sqrt(T/(rho*Lz));

inPoint = [0.52*Lx,0.53*Ly];
outPoint1 = [0.47*Lx,0.62*Ly];

stab1 = 2*sqrt(kappa*k);
stab2 = c*k/sqrt(2);
h = max([stab1,stab2]) + 0.0005;

Nx = floor(Lx/h);
Ny = floor(Ly/h);

lo = floor(inPoint(1)/h);
mo = floor(inPoint(2)/h);
alphax = inPoint(1)/h - lo;
alphay = inPoint(2)/h - mo;

Jcoeff = zeros(Nx-1,Ny-1);
Jcoeff(lo,mo) = (1-alphax) * (1-alphay) / h^2;
Jcoeff(lo, mo+1) = (1-alphax)*alphay / h^2;
Jcoeff(lo+1, mo) = alphax*(1-alphay) / h^2;
Jcoeff(lo+1,mo+1) = alphax*alphay / h^2;

Jvec = zeros((Nx-1)*(Ny-1),1);

for i = 1:Nx-1
    for j = 1:Ny-1
        Jvec(j + (Ny-1)*(i-1)) = Jcoeff(i,j);
    end
end

totN = (Nx-1)*(Ny-1);

Dxx = sparse(toeplitz([-2/h^2;1/h^2;zeros(Nx-3,1)]));
Dyy = sparse(toeplitz([-2/h^2;1/h^2;zeros(Ny-3,1)]));
D = kron(eye(Nx-1), Dyy)+kron(Dxx, eye(Ny-1)); DD = D*D;
B = sparse((2*eye(totN)-kappa^2*k^2*DD)/(1+sigma0*k));
I = sparse(eye(totN));
C = ((1-sigma0*k)/(1+sigma0*k))*I;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%% C++ Header output
fileID = fopen('../xcodeproj/plate-fdtd-stiff/plateFDTDStiffData.h', 'wt');

if fileID < 3
    disp('error opening file');
end

fprintf(fileID, '#ifndef plateFDTDStiffData\n#define plateStiffFDTDData\n\n');

fprintf(fileID, '#define SR %i\n', SR);
fprintf(fileID, '#define h %i\n\n', h);

fprintf(fileID, '#define Nx %i\n', Nx);
fprintf(fileID, '#define Ny %i\n\n', Ny);
fprintf(fileID, '#define totN %i\n\n', totN);

fprintf(fileID, '#define Lx %i\n', Lx);
fprintf(fileID, '#define Ly %i\n\n', Ly);

fprintf(fileID, 'static float B[%i][%i] = {', size(B,1), size(B,2));
for i = 1:size(B,1)
    for j = 1:size(B,2)
        fprintf(fileID,'%f, ',full(B(i,j)));
    end
   fprintf(fileID,'\n');
end
fprintf(fileID, '};\n\n');

fprintf(fileID, 'static float D[%i][%i] = {', size(D,1), size(D,2));
for i = 1:size(D,1)
    for j = 1:size(D,2)
        fprintf(fileID,'%f, ',full(D(i,j)));
    end
   fprintf(fileID,'\n');
end
fprintf(fileID, '};\n\n');

fprintf(fileID, '#endif');
status = fclose(fileID);

%% Simulation
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing vectors
% u = zeros(Nx, Ny);
% uPrev = zeros(Nx, Ny);
% uNext = zeros(Nx, Ny);
coeff = 1+sigma0*k;
muSq = kappa^2*k^2/h^4;
u = zeros(totN, 1);
uPrev = zeros(totN, 1);
uNext = zeros(totN, 1);

output = zeros(1,timeSamples);
tic
for n = 1:timeSamples
    exc = excit(n);

    % for l = 2:Nx-1
    %     for m = 2:Ny-1
    %         A = u(l,m+1) + u(l,m-1) + u(l+1,m) + u(l-1,m);
    %         Aprev = uPrev(l,m+1) + uPrev(l,m-1) + uPrev(l+1,m) + uPrev(l-1,m);
    % 
    %         B = u(l+1,m+1) + u(l+1,m-1) + u(l-1,m+1) + u(l-1,m-1);
    % 
    %         if l == 2
    %             if m == 2
    %                 C = u(l,m+2) - u(l,2) + u(l+2,m) - u(2,m);
    %             elseif m == Ny-1
    %                 C = - u(l,Ny-1) + u(l,m-2) + u(l+2,m) - u(2,m);
    %             else
    %                 C = u(l,m+2) + u(l,m-2) + u(l+2,m) - u(2,m);
    %             end
    %         elseif m == 2
    %             if l == Nx-1
    %                 C = u(l,m+2) - u(l,2) - u(Nx-1,m) + u(l-2,m);
    %             else
    %                 C = u(l,m+2) - u(l,2) + u(l+2,m) + u(l-2,m);
    %             end
    %         elseif l == Nx-1
    %             if m == 2
    %                 C = u(l,m+2) - u(l,2) - u(Nx-1,m) + u(l-2,m);
    %             elseif m == Ny-1
    %                 C = - u(l,Ny-1) + u(l,m-2) - u(Nx-1,m) + u(l-2,m);
    %             else
    %                 C = u(l,m+2) + u(l,m-2) - u(Nx-1,m) + u(l-2,m); 
    %             end
    %         elseif m == Ny-1
    %             if l == Nx-1
    %                 C = - u(l,Ny-1) + u(l,m-2) - u(Nx-1,m) + u(l-2,m);
    %             else
    %                 C = - u(l,Ny-1) + u(l,m-2) + u(l+2,m) + u(l-2,m); 
    %             end
    %         else
    %             C = u(l,m+2) + u(l,m-2) + u(l+2,m) + u(l-2,m);
    %         end
    % 
    %         uNext(l,m) = (2 - 20*muSq)*u(l,m)/coeff + 8*muSq*A/coeff - 2*muSq*B/coeff - muSq*C/coeff + c^2*k^2*A/h^2/coeff - (sigma0*k-1)*uPrev(l,m)/coeff + k^2*Jcoeff(l,m)*exc/coeff;% + 2*sigma1*(A-Aprev)/(h^2*coeff*k); 
    %     end
    % end

    % uNext = B*u-C*uPrev + k^2*Jvec*exc/coeff;
    % uNext = B*u + (k^2*c^2 + 2*sigma1*k)/(1 + sigma0*k)*D*u + (sigma0*k - 1)/(1 + sigma0*k)*I*uPrev - (2*k*sigma1)/(1 + sigma0*k)*D*uPrev + k^2*Jvec*exc/coeff;
    uNext = B*u + k^2*c^2/(1 + sigma0*k)*D*u + (sigma0*k - 1)/(1 + sigma0*k)*I*uPrev + k^2*Jvec*exc/coeff;

    uPrev = u;
    u = uNext;
    
    % u2D = reshape(u, Nx-1,Ny-1);
    % surf(u2D)
    % drawnow
    % pause
    output(n) = u(floor(totN/2));
end
realTimeFrac = toc/durSec

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Plotting & Saving Audio 
if normalizeOut
    output = output/max(abs(output));
end

figure(1)
plot(output);

if ~strcmp(inputType,'audiofile')
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
            E = 1.6e11;
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