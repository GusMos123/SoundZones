%Top gun test
clc
clear
close all;

[x,fs]=audioread("dangerzone.mp3");



bin_duration=fs/50; %antal sample points som mostvarar 20 ms (ish stationär process)
no_bins=ceil(length(x)/bin_duration); %antal binns vi behöver
X=zeros(bin_duration,no_bins); %här lagrar vi alla sekvenser i frekvensplanet

for bin=1:no_bins-1
    X(:,bin)=fft(x((bin-1)*bin_duration+1:bin*bin_duration,1),bin_duration); %Om nfft=no_bins blir det bra
end

%%
soundout = zeros(length(X)*bin_duration,1); %allokerar en vektor att spara resultatet i
for bin = 1:no_bins-1
    soundout((bin-1)*bin_duration+1:bin*bin_duration) = ifft(X(:,bin),bin_duration);
end

sound(real(soundout),fs) %lyssna på topgun

%%
indexes=30000:300100;
norm(x(indexes,1)-soundout(indexes))



plot(indexes,[x(indexes,1) soundout(indexes)])