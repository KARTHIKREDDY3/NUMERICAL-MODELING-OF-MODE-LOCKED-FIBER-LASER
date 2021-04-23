function []=Simulator()
% Input pulse parameters
shape = 1; % 1 for Gaussian, 2 for sech^2, 3 for noise
duration = 20.8; % Pulse duration in [fs], discarded for noise
energy = 0.8; % Pulse energy in [nJ]
lambda_c = 1030; % Central wavelength of the input pulse in [nm]
frep = 75; % Pulse repetition frequency in [MHz]
% Simulation parameters
timewidth = 10; % Width of the time window in [ps]
tres =10; % Resolution of the time window in [fs]
roundtrips = 100;

duration = duration*1e-3; % Conversion to [ps]
frep = frep*1e-6; % Conversion to [THz]
tres = tres*1e-3; % Conversion to [ps]
c = 299792.458; % [nm/ps]

t = -timewidth*0.5:tres:timewidth*0.5; % Time in [ps]
timesteps=length(t);

if shape==1
 E=sqrt(energy/duration)*0.969215.*exp(-0.5*(t./duration).^2*4*log(2)); % [kW^0.5]
elseif shape==2
 E=sqrt(energy/duration)*0.9417*sech(t*2*sech(0.5)/duration); % [kW^0.5]
else
 E=zeros(1,timesteps);
 
 for i=1:(timesteps+1)/2
 E(i)=rand(1,1)*2-1;
 end
 
 for i=1:(timesteps-1)/2
 E(i+(timesteps+1)/2)=E((timesteps+1)/2-i);
 end
 E = E.*sqrt(energy/trapz(t,abs(E).^2));
end
% figure;
% plot(t,abs(E).^2);
% grid on;
% axis tight;

fprintf('Energy: %e\n',trapz(t,abs(E).^2));

fs=1/timewidth;
freq=c/lambda_c+fs*linspace(-(timesteps-1)/2,(timesteps-1)/2,timesteps); % [THz]
wave=c./freq; % [nm]
w=2*pi*c/lambda_c; % [THz]
omegas = 2*pi*freq; % [THz]

%Fibers definitions
fiberlengths = [30,20,6.6,350].*1e7; %Fiber lengths in [nm]. Units are [cm] in the paranthesis
nsteps = [60,40,1,700];
gammas = [1.147973417107121e-08,5.904095002854001e-09,0,5.904095002854001e-09];
betas = [26.16e-12,13.04e-15;24.8e-12,23.3e-15;-1509.2e-12,2841.0e-15;24.8e-12,23.3e-15];
%betas = [-26.16e-12,13.04e-15;-24.8e-12,23.3e-15;-24.8e-12,23.3e-15];
Bs = [B_calc(betas(1,:),omegas,w);B_calc(betas(2,:),omegas,w);B_calc(betas(3,:),omegas,w);B_calc(betas(4,:),omegas,w)];
dopings = [8e-2,0,0,0];
Pps = [150e-6,0,00,];
%Pps = [15e-6,0,0];
lambda_ps = [976,976,976,976];
dl_ps = [0.3,0.3,0.3,0.3];
frs = [0.18,0.18,0.18,0.18];

comp_pos=[56.25,100.00,231.25,362.50,463.00,497.00,611.50,691.67,793.67,835.50,930.00,1080.00,1215.00].*1e-7; %1/nm
peak_int=[1.00,11.40,36.67,67.67,74.00,4.50,6.80,4.60,4.20,4.50,2.70,3.10,3.00]; %Unitless
gau_FWHM=[52.10,110.42,175.00,162.50,135.33,24.50,41.50,155.00,59.50,64.30,150.00,91.00,160.00].*1e-7; %1/nm
lor_FWHM=[17.37,38.81,58.33,54.17,45.11,8.17,13.83,51.67,19.83,21.43,50.00,30.33,53.33].*1e-7; %1/nm
hr=zeros(1,timesteps);
for a=ceil(timesteps/2):timesteps
    for b=1:13
        hr(a)=hr(a)+peak_int(b)*exp(-pi*c*t(a)*lor_FWHM(b))*exp(-((pi*c*gau_FWHM(b))^2)*(t(a)^2)/4)*sin(2*pi*c*comp_pos(b)*t(a));
    end
end
hr_integral=trapz(t,hr);
hr=hr./hr_integral;
hr=fft(hr);

hrs=[hr,hr,hr,hr];
OCs=[1,0.5,0.1,1];
%OCs=[1,1,1];

h1 = figure;
h2 = figure;
h3 = figure;
width = 944;
height = 428;
set(h1, 'Position', [8 569 width height]);
set(h2, 'Position', [width+24 height+141 width height]);
set(h3, 'Position', [width+24 49 width height]);

for lap=1:roundtrips
    fprintf('Roundtrip %d\n',lap);
    for segment=1:length(fiberlengths)
        if dopings(segment)==0
            E=Passive(E,w,gammas(segment),Bs(segment,:),frs(segment),hrs(segment),fiberlengths(segment),nsteps(segment),tres);
       
        end
        fprintf('Energy Out Segment %d: %e\n',segment,trapz(t,abs(E).^2));
        E=E.*sqrt(OCs(segment));
        fprintf('Energy Out Segment %d OC: %e\n',segment,trapz(t,abs(E).^2));
        if(segment==2)
            Esave = E;
            spect=fftshift(fft(Esave));
            spect=spect./max(spect);            
            
            plot(wave,abs(spect).^2);
            grid on;
            xlabel('Wavelength (nm)');
            ylabel('Amplitude (a.u.)');
            xlim([lambda_c-100 lambda_c+100]);
            title(['Roundtrip ',num2str(lap)]);
            drawnow
            
            figure(h2);
            plot(t,abs(Esave).^2);
            grid on;
            xlabel('Time Delay (ps)');
            ylabel('Power (kW)');
            axis tight;
            title(['Roundtrip ',num2str(lap)]);
            drawnow
            %if(lap>5)
                %Energies(lap-5) = trapz(t,abs(E).^2);
                Energies(lap) = trapz(t,abs(E).^2);
                figure(h3);
                plot(Energies);
                grid on;
                axis tight;
                xlabel('Roundtrip Number');
                ylabel('Energy (nJ)');
                drawnow
            %end
            
            %E=SESAM(E,t,1,0,2,1);
            E=NPE(E,0.05,0.95,3.0,0);
            
            fprintf('Energy out NPE: %e\n',trapz(t,abs(E).^2));
            E=Filter(E,wave,1,1030,10);
            fprintf('Energy out Filter: %e\n',trapz(t,abs(E).^2));
        end
%         for j=1:16
%             E(j)=E(j)*exp(-17+j);
%             E(timesteps-j+1)=E(timesteps-j+1)*exp(-17+j);
%         end
    end
end

% spect=fftshift(fft(Esave));
% spect=spect./max(spect);
% figure;
% plot(wave,abs(spect).^2);
% grid on;
% xlabel('Wavelength (nm)');
% ylabel('Amplitude (a.u.)');
% xlim([lambda_c-100 lambda_c+100]);
% 
% figure;
% plot(t,abs(Esave).^2);
% grid on;
% xlabel('Time Delay (ps)');
% ylabel('Power (kW)');
% axis tight;

%fprintf('Energy out: %e\n',trapz(t,abs(E).^2));

end
