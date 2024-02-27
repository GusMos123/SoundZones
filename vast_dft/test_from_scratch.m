clc
clear

%TODO
%1. Fixa så vi får med den relevanta tidsförskjutningen av impulssvaret,
%behövs egentligen bara för första binnen
%2. Fixa ett fungerande filterfan

%% Testing one source and one microphone
%Generate a very simple impulse response
fs=44100; % sample rate
mic=[19 19 1.8]; % mic position (x,y,z)
n=5; % hur nogrant den ska simulera studs i väggar etc
r=0.15; %reflection coefficient for the walls, in general -1<R<1. nära 0 ger ingen reflektion
rm=[20 20 3]; %row vector giving the dimensions of the room.
src=[1 1 1]; %row vector giving the x,y,z coordinates of the sound source.


%räknar ut impulssvar
h=rir(fs, mic, n, r, rm, src); 

[x,fs]=audioread("dangerzone.mp3");
x=x(:,1); % vill bara ha ena spåret

x=x(find(abs(x)>1e-3,1):end); %tar bort lite onödigt ljud i början

%dela upp problemet
bin_length=length(h);
no_bins=floor(length(x)/bin_length);

Ly=bin_length+length(h)-1;
Ly2=pow2(nextpow2(Ly));

y=zeros(bin_length*no_bins,1);
index=find(abs(h)>1e-3,1);

for bin=1:no_bins
    indices=(bin-1)*bin_length+1:bin*bin_length;
    % Find smallest power of 2 that is > Ly
    X=fft(x(indices), Ly2);		   % Fast Fourier transform
    H=fft(h, Ly2);	           % Fast Fourier transform
    Y=X.*H;        	           %
    y_temp=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
    y_temp=y_temp(1:1:Ly);               % Take just the first N elements
    y_temp=y_temp/max(abs(y_temp));           % Normalize the output

    result_indices=(bin-1)*bin_length+1:bin*bin_length; %picking out correct indeces
    result_indices=result_indices+index-1;
    y(result_indices)=y_temp(index:index+bin_length-1);
end

sound(y,fs);

%%
indices=1:40000;

plot([y(indices)/max(abs(y(indices))),x(indices)/max(abs(x(indices)))]);
legend(["with RIR","original"])

%% Trying out a filter to just compensate for the issues of the room FUNKAR EJ!

%choose a random index:
bin=2;
indeces=(bin-1)*bin_length+1:bin*bin_length;


index=find(abs(h)>1e-3,1); %finding when the sound reaches the mic
X=fft(x(indeces),Ly2);

desired_sound=conv(h,x(indeces));

desired_indeces=index:bin_length+index-1;
desired_sound=desired_sound(desired_indeces);

figure
hold on
plot(desired_sound./norm(desired_sound))
plot(x(indeces)./norm(x(indeces)))
hold off

desired_sound_fft=fft(desired_sound/max(abs(desired_sound)),Ly2);
desired_sound_fft=desired_sound_fft/max(abs(desired_sound_fft));

Q=zeros(Ly2,1);


%%
%do transform
%find when to start

XH=X.*H;
XH=XH/max(abs(desired_sound_fft));
figure
plot(abs(desired_sound_fft-XH))


%%
cvx_begin
    variable Q(Ly2)
    minimize(norm(desired_sound_fft-X.*H.*Q))
cvx_end

figure
plot(abs(desired_sound_fft-XH.*Q))

filtered_sound=real(ifft(XH.*Q,Ly2));
filtered_sound=filtered_sound(1:1:Ly);               % Take just the first N elements
filtered_sound=filtered_sound/max(abs(filtered_sound));           % Normalize the output
filtered_sound=filtered_sound(1:length(indeces));


figure
plot([filtered_sound,x(desired_indeces)/max(abs(x(desired_indeces)))])

%% Testing using filter

Ly=bin_length+length(h)-1;
Ly2=pow2(nextpow2(Ly));

y=zeros(bin_length*no_bins,1);

for bin=1:no_bins
    indeces=(bin-1)*bin_length+1:bin*bin_length;
    % Find smallest power of 2 that is > Ly
    X=fft(x(indeces), Ly2);		   % Fast Fourier transform
    H=fft(h, Ly2);	           % Fast Fourier transform
    Y=X.*H.*Q;        	           %
    y_temp=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
    y_temp=y_temp(1:1:Ly);               % Take just the first N elements
    y_temp=y_temp/max(abs(y_temp));           % Normalize the output
    index=find(abs(y_temp)>1e-2,1); %finding when the sound reaches the mic
    result_indices=(bin-1)*bin_length+1:bin*bin_length;
    result_indices=result_indices+index-1;
    y(result_indices)=y_temp(index:index+bin_length-1);
end


sound(y,fs);
%% Testing 16 speakers and 2 microphones. From now on code is sort of good
mic=[6 19 1.8; 14 19 1.8];
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

H=zeros(length(src(:,1)),length(mic(:,1)),Ly2);
for speaker=1:length(src)
    for microphone=1:length(mic(:,1))
        simulated_rir=rir(fs,mic(microphone,:),n,r,rm,src(speaker,:));
        H(speaker,microphone,:)=fft(simulated_rir,Ly2);
    end
end
%% Listening time! Ok funkar rätt dåligt atm
mic=1;
y=zeros(size(x));
for bin=1:no_bins
    indeces=(bin-1)*bin_length+1:bin*bin_length;
    X=fft(x(indeces), Ly2);		   % Fast Fourier transform
    y_temp=zeros(size(y_temp));
    for speaker=1:16
        transformed_rir=squeeze(H(speaker,mic,:));
        Y=X.*transformed_rir;
        y_temp_2=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
        y_temp_2=y_temp_2(1:1:Ly);               % Take just the first N elements
        y_temp_2=y_temp_2/max(abs(y_temp_2));           % Normalize the output
        
        y_temp=y_temp+y_temp_2;
    end
    index=find(abs(y_temp_2)>1e-2,1); %finding when the sound reaches the mic
    result_indices=(bin-1)*bin_length+1:bin*bin_length;
    result_indices=result_indices+index-1;
    y(result_indices)=y_temp(index:index+bin_length-1);

end

sound(y,fs);
%%
indices=1:40000;
plot([x(indices)/max(abs(x(indices))),y(indices)/max(abs(y(indices)))]);
legend(["original", "with RIR"])
