
clc
clear
close all

%generera 1 kHz
fs=16000;
% Duration of the signal (seconds)
duration = 1;  % Example duration, you can change it as needed

% Time vector
t = 0:1/fs:duration-1/fs;

% Frequency of the tone (Hz)
f_tone = 1000;  % 1 kHz tone

% Generate the sinusoidal signal
x = sin(2*pi*f_tone*t);

nfft=4096;


X=fft(x,nfft);

% Frequency axis
f_axis = (0:length(X)-1)*(fs/length(X));

std_periodogram=[abs(X(2049:end)) abs(X(1:2048))];

%korrekt periodogram
plot(linspace(-8000,8000,4096),std_periodogram);


%%
y=ifft(X,nfft);
sound(repmat(y,[1 500]),fs);

%% Testa nu med bins
% Add paths for relevant files
addpath(fullfile(pwd,'subcodes'))

% Initialization
datafilename = 'srir_fs16kHz_vast_dft__rir_generator.mat';
simulator_room_impulse_response = load(fullfile(pwd,'rirdata',datafilename));

frequency_bin_edges=simulator_room_impulse_response.varout.freq;
clear simulator_room_impulse_response;
clear datafilename;

X_selected=X(1:2048);

selected_frequencies=linspace(0,8000,2048);

X_binned=zeros(1,121);

for i=1:120
    X_binned(i)=sum( X_selected(selected_frequencies>frequency_bin_edges(i) & selected_frequencies<frequency_bin_edges(i+1)));
end




%vika tillbaka skiten
%X_binned=[X_binned flip(X_binned)];

%borde ändå inte funka
y=ifft(X_binned,4096);


y=repmat(y,[1 10000]);
sound(real(y),fs);
