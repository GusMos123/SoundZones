clc
clear

%Generate a very simple impulse response
fs=44100; % sample rate
mic=[19 18 1.6]; % mic position (x,y,z)
n=12; %
r=0.3; %reflection coefficient for the walls, in general -1<R<1.
rm=[20 19 21]; %reflection coefficient for the walls, in general -1<R<1.
src=[5 2 1]; %row vector giving the x,y,z coordinates of the sound source.
h=rir(fs, mic, n, r, rm, src);


%%




[x,fs]=audioread("dangerzone.mp3");
x=x(:,1);
bin_duration=fs/50; %antal sample points som mostvarar 20 ms (ish stationär process)
nfft=bin_duration;

%% Testar 1 kHz

% % Duration of the signal (seconds)
% duration = 1;  % Example duration, you can change it as needed
% 
% % Time vector
% t = 0:1/fs:duration-1/fs;
% 
% % Frequency of the tone (Hz)
% f_tone = 1000;  % 1 kHz tone
% 
% % Generate the sinusoidal signal
% x = sin(2*pi*f_tone*t);

H=fft(h,2^16);

bin_length=floor(length(H)/nfft); %antal binns vi behöver
% hitta gränserna för H
limits=round(linspace(-fs/2,fs/2,882));

H_temp=zeros(1,nfft);
for i=1:length(limits)-2
    H_temp(i)=sum(H((i-1)*bin_length+1:i*bin_length));
end

% ett känt fel - sista binnarna blir åt helvete


no_bins=ceil(length(x)/bin_duration); %antal binns vi behöver
X=zeros(nfft,no_bins); %här lagrar vi alla sekvenser i frekvensplanet

for bin=1:no_bins-1
    X(:,bin)=fft(x((bin-1)*bin_duration+1:bin*bin_duration,1),nfft).*H_temp'; %Om nfft=no_bins blir det bra
end

%%
soundout = zeros(length(X)*nfft,1); %allokerar en vektor att spara resultatet i
for bin = 1:no_bins-1
    soundout((bin-1)*nfft+1:bin*nfft) = ifft(X(:,bin),nfft);
end

sound(real(soundout),round(fs*(nfft/bin_duration))) %lyssna på topgun

%%
indexes=30000:300100;
norm(x(indexes,1)-soundout(indexes))



plot(indexes,[x(indexes,1)/norm(x(indexes,1)) soundout(indexes)/norm(soundout(indexes))])