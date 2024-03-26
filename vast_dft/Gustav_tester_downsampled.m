%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Testing 16 speakers and 2 microphones. From now on code is sort of good
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

[x,fs]=audioread("dangerzone.mp3");
x=x(:,1); % vill bara ha ena spåret

x=resample(x,16000,fs); 
fs=16000;

x=x(find(abs(x)>1e-3,1):end); %tar bort lite onödigt ljud i början
x=x(1:fs*5); %tar första 5 sekunderna. Det är detta vi faktiskt jobbar med

%mic=[6 19 1.8; 14 19 1.8;6.2,19,1.8]; %mic positions
mics_per_zone=2;

bright_zone=generate_zound_zone([6,15,1.8],mics_per_zone,0.2);
dark_zone=generate_zound_zone([14,15,1.8],mics_per_zone,0.2);
mic=[bright_zone;dark_zone];

n=5; %accuracy of RIR
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

h=rir(fs, mic, n, r, rm, src); 
Ly=length(x)+length(h)-1;  %Hitta längden på faltningen som vi kommer få
Ly2=pow2(nextpow2(Ly));    % Find smallest power of 2 that is > Ly

%% Generate RIRs

H=zeros(length(mic(:,1)),length(src(:,1)),Ly2);
h=zeros(length(mic(:,1)),length(src(:,1)),Ly2);

for microphone=1:length(mic(:,1))
    for speaker=1:length(src(:,1))
        simulated_rir=rir(fs,mic(microphone,:),n,r,rm,src(speaker,:))';
        h(microphone,speaker,1:length(simulated_rir))=simulated_rir;
        H(microphone,speaker,:)=fft(simulated_rir,Ly2);
    end
end

%%
plot(abs(squeeze(H(1,1,:))'))
%% Listening time!
selected_mic=4; %här lyssnar vi

X=fft(x, Ly2);		   % Fast Fourier transform
Y=zeros(size(X));
y_time_domain=zeros(80000,1);

for speaker=1:16
    rir_fft_domain=squeeze(H(selected_mic,speaker,:)); %finding the RIR
    Y=Y+X.*rir_fft_domain;
    y_time_domain=y_time_domain+filter(squeeze(h(selected_mic,speaker,:)),1,x);
end

y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly);               % Take just the first N elements
y=y/max(abs(y));           % Normalize the output

y_time_domain=y_time_domain/max(abs(y_time_domain));

sound(y_time_domain,fs);
y=y(1:80000);
diff=y-y_time_domain;
%%
indices=1:40000;
plot([x(indices)/max(abs(x(indices))),y(indices)/max(abs(y(indices)))]);
legend(["original", "with RIR"])

%% Filter design time!
ideal_response=rir(fs, mic(1,:), n, 0, rm, src(16,:)); 

figure
plot(ideal_response);
title("Ideal room impulse reponse")

ideal_response_fft_bright=fft(ideal_response,Ly2);

H_bright=squeeze(H(1:mics_per_zone,:,:)); %choosing bright zone for top gun
H_dark=squeeze(H(mics_per_zone+1:end,:,:));

%% Optimization
mu=0.3; %how much dark zone should matter

cvx_solver_settings( 'max_iter', 2 )
cvx_begin
  
 variable Q_start(16,Ly2) %ett filter per högtalare (16 st)
 %minimize( norm( mean( sum( H.*reshape( repmat( Q_start, [2*mics_per_zone,1,1]), size(H))-reshape(repmat(ideal_response_fft, [2*mics_per_zone, speaker, 1]), size(H)) ) ) ) )
 minimize(goal_function(Q_start, H_bright, H_dark, ideal_response_fft_bright, mu)) %lite oklart med vilka inputs...
cvx_end

%%

Q_topgun=full(Q_start);

%% Listening time!
desired_mic=3; %här lyssnar vi

X=fft(x, Ly2);		   % Fast Fourier transform
Y=zeros(size(X));
y_time=zeros(80000,1);

for speaker=1:16
    rir_fft_domain=squeeze(H(desired_mic,speaker,:));
    filter_component=reshape(Q_topgun(speaker,:),[Ly2 1]); %picking out correct filter
    Y=Y+X.*rir_fft_domain.*filter_component;
    
    rir_time_domain=real(ifft(rir_fft_domain,Ly2));
    filter_time=real(ifft(filter_component,Ly2));
    filter_time=filter_time(1:Ly);
    temp=filter(filter_time,1,x);
    y_time=y_time + filter(rir_time_domain,1,temp);

end

y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform           % Take just the first N elements
y=y(1:1:Ly);
y=y*40;           % Normalize the output

y_time=real(y_time);

y_time=y_time*500;

sound(y_time,fs)



%% INTE KOMMIT LÄNGRE ÄN HIT
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
   
    rir_fft_domain=squeeze(H(speaker,mic,:));
    topgun_filter=reshape(Q_topgun(speaker,:),[Ly2 1]);
    taylor_filter=reshape(Q_taylor(speaker,:),[Ly2 1]);
    Y=Y+X_topgun.*rir_fft_domain.*topgun_filter + X_taylor.*rir_fft_domain.*taylor_filter;

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
mic=2; %här lyssnar vi

X_topgun=fft(x, Ly2);		   % Fast Fourier transform
X_taylor=fft(xt,Ly2);
Y=zeros(size(X));
output=zeros(16,Ly);
for speaker=1:16
    topgun_filter=reshape(Q_topgun(speaker,:),[Ly2 1]);
    taylor_filter=reshape(Q_taylor(speaker,:),[Ly2 1]);

    rir_fft_domain=squeeze(H(speaker,mic,:));
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
mic=1;
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
        rir_fft_domain=squeeze(H(speaker,mic,:));
        Y=Y+X_topgun.*rir_fft_domain.*topgun_filter + X_taylor.*rir_fft_domain.*taylor_filter;

    end

    y_temp=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
    y_temp=y_temp(1:1:Ly);               % Take just the first N elements
    y_temp=y_temp/max(abs(y_temp));           % Normalize the output
  
    y(indices)=y_temp(1:Ly);
end

sound(y,fs)

%% Lyssna på hela låten
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