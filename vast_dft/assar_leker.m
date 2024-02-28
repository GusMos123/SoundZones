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
h=zeros(length(src(:,1)),length(mic(:,1)),Ly2);
for speaker=1:length(src)
    for microphone=1:length(mic(:,1))
        simulated_rir=rir(fs,mic(microphone,:),n,r,rm,src(speaker,:));
        h(speaker,microphone,1:length(simulated_rir))=simulated_rir;
        H(speaker,microphone,:)=fft(simulated_rir,Ly2);
    end
end
%% Listening time! Ok funkar rätt dåligt atm
selected_mic=1; %här lyssnar vi

X=fft(x, Ly2);		   % Fast Fourier transform
Y=zeros(size(X));

for speaker=1:16
    transformed_rir=squeeze(H(speaker,selected_mic,:));
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

%% Filter design time!
ideal_response=rir(fs, mic(1,:), n, 0, rm, src(8,:)); 

figure
plot(ideal_response);
title("Ideal room impulse reponse")

ideal_response_fft_bright=fft(ideal_response,Ly2);

H_bright=squeeze(H(:,1,:)); %choosing bright zone for top gun
H_dark=squeeze(H(:,2,:));

%% Optimization
mu=0.3;

cvx_precision('low')
cvx_begin
  
 variable Q_start(16,Ly2)
 minimize(norm(sum(H_bright.*Q_start,1)' -ideal_response_fft_bright) + mu*norm(sum(H_dark.*Q_start,1)))

cvx_end

%%

Q_topgun=Q_start;

%% Listening time!
mic=1; %här lyssnar vi

X=fft(x, Ly2);		   % Fast Fourier transform
Y=zeros(size(X));

for speaker=1:16
    transformed_rir=squeeze(H(speaker,mic,:));
    Y=Y+X.*transformed_rir.*Q_end(speaker,:)';
end

y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

sound(y,fs)

%% Now we optimize for Taylor Swift
[xt,fs]=audioread("blankspace.mp3");
xt=xt(:,1);
xt=xt(find(abs(xt)>1e-3):end);
xt=xt(fs*5:fs*10);
%% Filter design time!
ideal_response=rir(fs, mic(1,:), n, 0, rm, src(8,:)); 

figure
plot(ideal_response);
title("Ideal room impulse reponse")

ideal_response_fft_bright=fft(ideal_response,Ly2);

H_bright=squeeze(H(:,2,:)); %bright zone for swift is zone 2
H_dark=squeeze(H(:,1,:));

%% Optimization
mu=0.3;
cvx_begin
 variable Q_taylor(16,Ly2)
 minimize(norm(sum(H_bright.*Q_taylor,1)' -ideal_response_fft_bright) + mu*norm(sum(H_dark.*Q_taylor,1)))

cvx_end

%%
mic=1; %här lyssnar vi

X_topgun=fft(x, Ly2);		   % Fast Fourier transform
X_taylor=fft(xt,Ly2);
Y=zeros(size(X));

Effective_filter=zeros(Ly2,1);

for speaker=1:16
   
    transformed_rir=squeeze(H(speaker,mic,:));
    Y=Y+X_topgun.*transformed_rir.*Q_topgun(speaker,:)' + X_taylor.*transformed_rir.*Q_taylor(speaker,:)';

     Effective_filter=Effective_filter+transformed_rir.*Q_topgun(speaker,:)';
end
effective_filter_time_domain=real(ifft(Effective_filter,Ly2));
effective_filter_time_domain=effective_filter_time_domain(1:1:Ly);
effective_filter_time_domain=effective_filter_time_domain(1:1:Ly)/max(abs(effective_filter_time_domain));

y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

sound(y,fs)

%% Lyssna på output från en högtalare
speaker=4;
Y=X_topgun.*Q_topgun(speaker,:)' + X_taylor.*Q_taylor(speaker,:)';
y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

sound(y,fs)


