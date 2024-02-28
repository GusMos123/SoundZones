%% Testing one source and one microphone
clc
clear

%Generate a very simple impulse response
fs=44100; % sample rate
mic=[19 19 1.8]; % mic position (x,y,z)
n=5; % hur nogrant den ska simulera studs i väggar etc
r=0.8; %reflection coefficient for the walls, in general -1<R<1. nära 0 ger ingen reflektion
rm=[20 20 3]; %row vector giving the dimensions of the room.
src=[1 1 1]; %row vector giving the x,y,z coordinates of the sound source.


%räknar ut impulssvar
h=rir(fs, mic, n, r, rm, src); 

[x,fs]=audioread("dangerzone.mp3");
x=x(:,1); % vill bara ha ena spåret

x=x(find(abs(x)>1e-3,1):end); %tar bort lite onödigt ljud i början
x=x(1:fs*5); %tar första 5 sekunderna

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

plot([y(indices)/max(abs(y(indices))),x(indices)/max(abs(x(indices)))]);
legend(["with RIR","original"])

%% Trying out a filter to just compensate for the issues of the room FUNKAR EJ!

%find when RIR reaches speaker
time_delay=find(abs(h)>1e-2,1);

desired_sound=[zeros(time_delay,1);x];
desired_sound_fft=fft(desired_sound,Ly2);
desired_sound_fft=desired_sound_fft/max(abs(desired_sound_fft));

Q=zeros(Ly2,1);

XH=X.*H;
XH=XH/max(abs(XH));

figure
plot(abs(desired_sound_fft-XH))
title("Difference between desired sound and X*H before opt")


%% Optimization
cvx_begin
    variable Q(Ly2)
    minimize(norm(desired_sound_fft-XH.*Q))
cvx_end

figure
plot(abs(desired_sound_fft-XH.*Q))
title("Difference between desired sound and X*H after opt")

filtered_sound=real(ifft(XH.*Q,Ly2));
filtered_sound=filtered_sound(1:1:Ly);               % Take just the first N elements
filtered_sound=filtered_sound/max(abs(filtered_sound));           % Normalize the output

sound(filtered_sound,fs)

figure
plot([y, filtered_sound])
legend(["Output with RIR", "Output with RIR and Filter"])

figure
plot([desired_sound/max(abs(desired_sound)), filtered_sound(1:length(desired_sound))])
legend(["Desired Sound", "Output with RIR and filtering"])

%% Testing 16 speakers and 2 microphones. From now on code is sort of good
mic=[6 19 1.8; 14 19 1.8];
n=5; %
r=0.8; %reflection coefficient for the walls, in general -1<R<1. nära 0 ger ingen reflektion
rm=[20 20 3]; %row vector giving the dimensions of the room.
src=zeros(16,3); %row vector giving the x,y,z coordinates of the sound source.
src(:,2:3)=ones;
src(:,1)=linspace(5,15,16);

%plotting room
figure
hold on
plot([0 0],[0 20],'black','LineWidth',2) % room walls
plot([0 20],[0 0],'black','LineWidth',2)
plot([0 20],[20 20],'black','LineWidth',2)
plot([20 20],[0 20],'black','LineWidth',2)
plot(mic(:,1),mic(:,2),'*','LineWidth',5) % microphones
plot(src(:,1),src(:,2),'o','LineWidth',2) %Speakers

hold off

%% Generate RIRs

H=zeros(length(src(:,1)),length(mic(:,1)),Ly2);
for speaker=1:length(src)
    for microphone=1:length(mic(:,1))
        simulated_rir=rir(fs,mic(microphone,:),n,r,rm,src(speaker,:));
        H(speaker,microphone,:)=fft(simulated_rir,Ly2);
    end
end
%% Listening time! Ok funkar rätt dåligt atm
mic=1; %här lyssnar vi

X=fft(x, Ly2);		   % Fast Fourier transform
Y=zeros(size(X));

for speaker=1:16
    transformed_rir=squeeze(H(speaker,mic,:));
    Y=Y+X.*transformed_rir;
end

y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output




sound(y,fs);
%%
indices=1:40000;
plot([x(indices)/max(abs(x(indices))),y(indices)/max(abs(y(indices)))]);
legend(["original", "with RIR"])
