%Top gun test
clc
clear
close all;

[x,fs]=audioread("dangerzone.mp3");
x=x(:,1);

nfft=2^14;
time_per_segment=1000*nfft/fs
bin_duration=nfft; %om detta är lika med bin_duration blir det bra
no_bins=ceil(length(x)/bin_duration); %antal binns vi behöver
X=zeros(nfft,no_bins); %här lagrar vi alla sekvenser i frekvensplanet

for bin=1:no_bins-1
    X(:,bin)=fft(x((bin-1)*bin_duration+1:bin*bin_duration,1),nfft); %taking fourier transform
end

%%
soundout = zeros(length(X)*nfft,1); %allokerar en vektor att spara resultatet i
for bin = 1:no_bins-1
    soundout((bin-1)*nfft+1:bin*nfft) = ifft(X(:,bin),nfft); %inverstransform
end

sound(real(soundout),fs*nfft/bin_duration) %lyssna på topgun

%%
indexes=30000:300100;
norm(x(indexes,1)-soundout(indexes))



plot(indexes,[x(indexes,1) soundout(indexes)])