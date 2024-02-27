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


%% Laddar in top gun


[x,fs]=audioread("dangerzone.mp3");
x=x(:,1); % vill bara ha ena spåret

Ly=length(x)+length(h)-1;  % 
Ly2=pow2(nextpow2(Ly));    % Find smallest power of 2 that is > Ly
X=fft(x, Ly2);		   % Fast Fourier transform
H=fft(h, Ly2);	           % Fast Fourier transform
Y=X.*H;        	           % 
y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

sound(y,fs);

%%
indices=1:40000;
plot([x(indices),y(indices)]);
legend(["original", "with RIR"])