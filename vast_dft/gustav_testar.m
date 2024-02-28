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

%% Trying out a filter to just compensate for the issues of the room

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



% ------------------------------------------------------------
%% Testing 16 speakers and 2 microphones. From now on code is sort of good
% ------------------------------------------------------------


mic=[6 19 1.8; 14 19 1.8];
n=5; %
r=1; %reflection coefficient for the walls, in general -1<R<1. nära 0 ger ingen reflektion
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
h=zeros(length(src(:,1)),length(mic(:,1)),Ly2);
for speaker=1:length(src)
    for microphone=1:length(mic(:,1))
        simulated_rir=rir(fs,mic(microphone,:),n,r,rm,src(speaker,:));
        h(speaker,microphone,1:length(simulated_rir))=simulated_rir;
        H(speaker,microphone,:)=fft(h(speaker,microphone,:),Ly2);
    end
end
%% Listening time! Ok funkar rätt dåligt atm
mic=2; %här lyssnar vi

X=fft(x, Ly2);		   % Fast Fourier transform
Y=zeros(size(X));

for speaker=1:16
    transformed_rir=squeeze(H(speaker,mic,:));
    Y=Y+X.*transformed_rir;
end

y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

sound(y, fs)
%%
indices=1:40000;
plot([x(indices)/max(abs(x(indices))),y(indices)/max(abs(y(indices)))]);
legend(["original", "with RIR"])

%% Gustav

%params
mu = 0.3;
Q_vec = zeros(16, Ly2);

for speaker=1:16
    Q=ones(Ly2,1);
    transformed_rir=squeeze(H(speaker,1,:));
    XH1=X.*transformed_rir;
    XH1=XH1/max(abs(XH1));
    
    transformed_rir=squeeze(H(speaker,2,:));
    XH2=X.*transformed_rir;
    XH2=XH2/max(abs(XH2));
    
    time_delay=find(abs(h(speaker,1,:))>1e-2,1); %hitta time delay till mic 1
    desired_sound=[zeros(time_delay,1);x];
    desired_sound_fft=fft(desired_sound,Ly2);
    desired_sound_fft=desired_sound_fft/max(abs(desired_sound_fft));
    fprintf("Before opt: %f \n", norm(desired_sound_fft-XH1.*Q) + mu * norm(-XH2.*Q));
    
    cvx_begin
        variable Q(Ly2)
        minimize(norm(desired_sound_fft-XH1.*Q) + mu * norm(-XH2.*Q))
    cvx_end
    Q_vec(speaker,:) = Q;
    speaker
end

%% Lyssna

mic=1; %här lyssnar vi

X=fft(x, Ly2);		   % Fast Fourier transform
Y=zeros(size(X));

for speaker=1:16
    transformed_rir=squeeze(H(speaker,mic,:));
    Y=Y+X.*transformed_rir.*Q(speaker,:)';
end

y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output
%y=y*10^9;
sound(y, fs)

%% Lyssna på vad en högtalare spelar
speaker=11;
transformed_rir=squeeze(H(speaker,mic,:));
Y=X.*transformed_rir;
y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

sound(y, fs)



%%
hold on
for i = 1:16
    plot(Q(i,:))
end


%% Köra allt i ett svep

%params
mu = 5;

XH1 = zeros(16,Ly2);
XH2 = zeros(16,Ly2);
desired_sound_fft = zeros(16,Ly2);

for speaker=1:16
    transformed_rir=squeeze(H(speaker,1,:));
    XH1(speaker, :)=X.*transformed_rir;
    XH1(speaker, :)=XH1(speaker, :)/max(abs(XH1(speaker, :)));
    
    transformed_rir=squeeze(H(speaker,2,:));
    XH2(speaker, :)=X.*transformed_rir;
    XH2(speaker, :)=XH2(speaker, :)/max(abs(XH2(speaker, :)));
    
    time_delay=find(abs(h(speaker,1,:))>1e-2,1); %hitta time delay till mic 1
    desired_sound=[zeros(time_delay,1);x];
    desired_sound_fft(speaker, :)=fft(desired_sound,Ly2);
    desired_sound_fft(speaker, :)=desired_sound_fft(speaker, :)/max(abs(desired_sound_fft(speaker, :)));
    
end

cvx_begin
    variable Q(16, Ly2)
    minimize(sum(sum(abs(desired_sound_fft-XH1.*Q) + mu * abs(-XH2.*Q))))
cvx_end

%% Köra allt i ett svep 2

%params
mu = 1;

XH1 = zeros(16,Ly2);
XH2 = zeros(16,Ly2);
desired_sound_fft = zeros(16,Ly2);

for speaker=1:16
    transformed_rir=squeeze(H(speaker,1,:));
    XH1(speaker, :)=X.*transformed_rir;
    XH1(speaker, :)=XH1(speaker, :)/max(abs(XH1(speaker, :)));
    
    transformed_rir=squeeze(H(speaker,2,:));
    XH2(speaker, :)=X.*transformed_rir;
    XH2(speaker, :)=XH2(speaker, :)/max(abs(XH2(speaker, :)));
    
    time_delay=find(abs(h(speaker,1,:))>1e-2,1); %hitta time delay till mic 1
    desired_sound=[zeros(time_delay,1);x];
    desired_sound_fft(speaker, :)=fft(desired_sound,Ly2);
    desired_sound_fft(speaker, :)=desired_sound_fft(speaker, :)/max(abs(desired_sound_fft(speaker, :)));
    
end

cvx_begin
    variable Q2(16, Ly2)
    minimize((desired_sound_fft-XH1.*Q2) + mu * (-XH2.*Q2))
cvx_end
