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

y_distorted=y;
sound(y,fs);

%%
indices=1:40000;

plot([y(indices)/max(abs(y(indices))),x(indices)/max(abs(x(indices)))]);
legend(["with RIR","original"])

%% Testing super simple filter
ideal_response=rir(fs, mic, n, 0, rm, src); 
ideal_response_fft=fft(ideal_response,Ly2);

cvx_solver_settings( 'max_iter', 2 )
cvx_begin
  variable Q_start(Ly2,1)
 minimize(norm(H.*Q_start-ideal_response_fft))
cvx_end

%% Testing filter in both domains
Y=X.*Q_start;
y=real(ifft(Y,Ly2));
y=y(1:Ly);
y=conv(y,h);
y=y(1:Ly);
y=y/max(abs(y));

q=ifft(Q_start,Ly2);
q=q(1:Ly/8);
y2=conv(conv(x,q),h);
y2=y2(1:Ly);
y2=y2/max(abs(y2));


sound(y,fs);

diff=y-y2;
%% Testing 16 speakers and 2 microphones. From now on code is sort of good
mic=[6 19 1.8; 14 19 1.8;6.2,19,1.8];
n=5; %
r=0.15; %reflection coefficient for the walls, in general -1<R<1. nära 0 ger ingen reflektion
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

H=zeros(length(mic(:,1)),length(src(:,1)),Ly2);
h=zeros(length(mic(:,1)),length(src(:,1)),Ly2);
for microphone=1:length(mic(:,1))
    for speaker=1:length(src)
        simulated_rir=rir(fs,mic(microphone,:),n,r,rm,src(speaker,:))';
        h(speaker,microphone,1:length(simulated_rir))=simulated_rir;
        H(speaker,microphone,:)=fft(simulated_rir,Ly2);
    end
end
%% Listening time!
selected_mic=3; %här lyssnar vi

X=fft(x, Ly2);		   % Fast Fourier transform
Y=zeros(size(X));
y_time_domain=zeros(Ly,1);
for speaker=1:16
    transformed_rir=squeeze(H(speaker,selected_mic,:));
    Y=Y+X.*transformed_rir;
    soundout=conv(x,squeeze(h(speaker,selected_mic,:)));
    y_time_domain=y_time_domain+soundout(1:Ly);
end

y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

y_time_domain=y_time_domain/max(abs(y_time_domain));

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

cvx_solver_settings( 'max_iter', 2 )
cvx_begin
  
 variable Q_start(16,Ly2)
 minimize(norm(sum(H_bright.*Q_start,1)' -ideal_response_fft_bright) + mu*norm(sum(H_dark.*Q_start,1)))

cvx_end

%%

load("filters.mat")
%Q_topgun=Q_start;

%% Listening time!
desired_mic=1; %här lyssnar vi

X=fft(x, Ly2);		   % Fast Fourier transform
Y=zeros(size(X));


for speaker=1:16
    transformed_rir=squeeze(H(speaker,desired_mic,:));
    filter_component=reshape(Q_topgun(speaker,:),[Ly2 1]);
    Y=Y+X.*filter_component.*transformed_rir;
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
mic=3; %här lyssnar vi

X_topgun=fft(x, Ly2);		   % Fast Fourier transform
X_taylor=fft(xt,Ly2);
Y=zeros(size(X));

Effective_filter=zeros(Ly2,1);

for speaker=1:16
   
    transformed_rir=squeeze(H(speaker,mic,:));
    topgun_filter=reshape(Q_topgun(speaker,:),[Ly2 1]);
    taylor_filter=reshape(Q_taylor(speaker,:),[Ly2 1]);
    Y=Y+X_topgun.*transformed_rir.*topgun_filter + X_taylor.*transformed_rir.*taylor_filter;

    Effective_filter=Effective_filter+topgun_filter;
end
effective_filter_time_domain=real(ifft(Effective_filter,Ly2));
effective_filter_time_domain=effective_filter_time_domain(1:1:Ly);
effective_filter_time_domain=effective_filter_time_domain(1:1:Ly)/max(abs(effective_filter_time_domain));

y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

sound(y,fs)

%% Test i tidsdomänen
q_taylor=zeros(16,round(Ly/2)); % de behöver ej vara så långa
q_topgun=zeros(16,round(Ly/2));
y=zeros(482644,1);
for speaker=1:16
    q_taylor_temp=real(ifft(Q_taylor(speaker,:),Ly2));
    q_topgun_temp=real(ifft(Q_topgun(speaker,:),Ly2));
    q_taylor(speaker,:)=q_taylor_temp(1:1:round(Ly/2));
    q_topgun(speaker,:)=q_topgun_temp(1:1:round(Ly/2));
    output_topgun=conv(x,1);
    output_taylor=conv(xt,1);
    sound_topgun=conv(output_topgun,squeeze(h(speaker,mic,:)));
    sound_taylor=conv(output_taylor,squeeze(h(speaker,mic,:)));
    y=y+[sound_topgun;0]+sound_taylor;
end
y=y/max(abs(y));

sound(y,fs)

%% Testar att mixa domäner
mic=1; %här lyssnar vi

X_topgun=fft(x, Ly2);		   % Fast Fourier transform
X_taylor=fft(xt,Ly2);
Y=zeros(size(X));
output=zeros(16,Ly);
for speaker=1:16
    topgun_filter=reshape(Q_topgun(speaker,:),[Ly2 1]);
    taylor_filter=reshape(Q_taylor(speaker,:),[Ly2 1]);

    transformed_rir=squeeze(H(speaker,mic,:));
    Y=Y+X_topgun.*topgun_filter + X_taylor.*taylor_filter;
    output_temp=real(ifft(Y,Ly2));
    output_temp=output_temp(1:Ly);
    output_temp=output_temp/max(abs(output_temp));
    output(speaker,:)=output_temp;
end

y=zeros(Ly,1);

for speaker=1:16
    output=conv(output(speaker,:),squeeze(h(speaker,mic,:)));
    y=y+output(1:Ly);
end

sound(y,fs)

%% Lyssna på output från en högtalare
speaker=2;
Y=X_topgun.*topgun_filter + X_taylor.*taylor_filter;
y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

sound(y,fs)

%% Possibility to play entire song
mic=2;
x_topgun=audioread("dangerzone.mp3");
x_topgun=x_topgun(:,1);
x_topgun=x_topgun(find(abs(x_topgun)>1e-2):end,1);

x_taylor=audioread("blankspace.mp3");
x_taylor=x_taylor(:,1);
x_taylor=x_taylor(find(abs(x_taylor)>1e-2):end,1);

%dela upp problemet
bin_length=Ly;
no_bins=floor(min(length(x_taylor),length(x_topgun))/bin_length);
%no_bins=5;

y=zeros(bin_length*no_bins,1);
index=find(abs(h)>1e-3,1);

for bin=1:no_bins
    indices=(bin-1)*bin_length+1:bin*bin_length;
    % Find smallest power of 2 that is > Ly
    X_topgun=fft(x_topgun(indices), Ly2);		   % Fast Fourier transform
    X_taylor=fft(x_taylor(indices), Ly2);	
    
    Y=zeros(size(X));

    for speaker=1:16
        topgun_filter=reshape(Q_topgun(speaker,:),[Ly2 1]);
        taylor_filter=reshape(Q_taylor(speaker,:),[Ly2 1]);
        transformed_rir=squeeze(H(speaker,mic,:));
        Y=Y+X_topgun.*transformed_rir.*topgun_filter + X_taylor.*transformed_rir.*taylor_filter;

    end

    y_temp=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
    y_temp=y_temp(1:1:Ly);               % Take just the first N elements
    y_temp=y_temp/max(abs(y_temp));           % Normalize the output
  
    y(indices)=y_temp(1:Ly);
end

sound(y,fs)