%% Kopia av kod, mest kladd
close all;clear;clc

% hej hej

% Add paths for relevant files
addpath(fullfile(pwd,'subcodes'))

% Initialization
datafilename = 'srir_fs16kHz_vast_dft__rir_generator.mat';

simulator_room_impulse_response = load(fullfile(pwd,'rirdata',datafilename));


% Processing details
general = simulator_room_impulse_response.general;

% Array geometry
loudspeaker_array = simulator_room_impulse_response.array;

% Room
room = simulator_room_impulse_response.room;

% Zone
zones = simulator_room_impulse_response.zone;

% Simulated room impulse responses (RIRs) by the rir_generator
% Note that rir_generator can be downloaded from the following link:
% https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
% https://github.com/ehabets/RIR-Generator
impulse_response_measured = simulator_room_impulse_response.irMeasured;

% The impulse response of the virtual source
% This is the RIR of the 8th loudspeaker in irMeasured
impulse_response_virtual_score = simulator_room_impulse_response.irVirsrc;

% Some variables
output_variables = simulator_room_impulse_response.varout;

clear simulator_room_impulse_response


% lägga till frequency-axeln:
frequencies=output_variables.dF;

%ok så det verkar som att varje kolumn i Impulse-response motsvarar en
%mikrofons impulssvar från samtliga 16 högtalare.'


%följande bit är otroligt oklart
control_filter = get_control_filter(general, output_variables);
targetdB = -10;
for jj = general.idx.vast_nf:general.idx.vast_t
    control_filter{jj}.cvxopt_properties.findopt = false;
    control_filter{jj}.cvxopt_properties.opttype = 'min_sd';
    control_filter{jj}.cvxopt_properties.const = 'nsb';
    control_filter{jj}.cvxopt_properties.tarval = 10^(targetdB/10);
    control_filter{jj}.mu = 1;
    control_filter{jj}.V = control_filter{jj}.Vmax;
end

% Initialization for VAST-NF
% 1 kHz
target_options={};
target_options.target_index = 16; %target frequency finns i bin no 16
target_options.target_frequency = (target_options.target_index-1)*output_variables.dF;
target_options.journal_exp_1 = false;


experiment_1_control_filter = control_filter{general.idx.vast_nf};
experiment_1_target_options = target_options;
experiment_1_target_options.journal_exp_1 = false;
experiment_1_control_filter.include_dc_and_nyqvist_frequencies = true;
experiment_1_target_options.taridx=16;

desired_room_impulse_responses = impulse_response_virtual_score; %otroligt oklart

%[ctrfilter, perform, rankcheck,q]=calculatefVAST(general, loudspeaker_array, zones, experiment_1_control_filter, impulse_response_measured, desired_room_impulse_responses, [], 'narrow', experiment_1_target_options);

%% Här börjar nåt viktigt

%generera 1 kHz

fs=general.fs;

% Duration of the signal (seconds)
duration = 1;  % Example duration, you can change it as needed

% Time vector
t = 0:1/fs:duration-1/fs;

% Frequency of the tone (Hz)
f_tone = 1000;  % 1 kHz tone

% Generate the sinusoidal signal
x = sin(2*pi*f_tone*t);

nfft=4096;

% Take the Fourier transform
tone_freq_domain = fft(x,nfft)';

% Frequency axis
frequencyAxis = linspace(-fs/2, fs/2, nfft);


 % Plot the magnitude spectrum
figure;
plot(frequencyAxis, abs(tone_freq_domain));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum of 1 kHz Tone, not shifted');
grid on;

%Steg 1: separera impulssvaren
all_impulse_responses_time_domain_zone_1=impulse_response_measured{1,1};
all_impulse_responses_time_domain_zone_2=impulse_response_measured{2,1};

impulse_responses_time_domain_zone_1=zeros(16,37,2967);
impulse_responses_time_domain_zone_2=zeros(16,37,2967);

impulse_responses_freq_domain_zone_1=zeros(16,37,nfft);
impulse_responses_freq_domain_zone_2=zeros(16,37,nfft);


%svårt att fouriera till frequency domain, eftersom vi vet inte riktigt
%vilken frekvensaxel vi får i frekvensdomänen.
control_filter=control_filter{1,1};
[Hml, Dm, hml, dm] = getTransferFunction(general,loudspeaker_array,zones,impulse_response_measured,desired_room_impulse_responses,false);
[q control_filter]=getfVASTnarrow(control_filter,Hml,Dm,experiment_1_target_options);

%% Hämta ut impulssvar och konvertera till användbara h:n
for speaker=1:16
    for microphone=1:37
         impulse_responses_time_domain_zone_1(speaker,microphone,:)=all_impulse_responses_time_domain_zone_1((speaker-1)*2967 +1:speaker*2967,microphone);
         impulse_responses_time_domain_zone_2(speaker,microphone,:)=all_impulse_responses_time_domain_zone_2((speaker-1)*2967 +1:speaker*2967,microphone);

         impulse_responses_freq_domain_zone_1(speaker,microphone,:)=fft(impulse_responses_time_domain_zone_1(speaker,microphone,:),nfft);
         impulse_responses_freq_domain_zone_2(speaker,microphone,:)=fft(impulse_responses_time_domain_zone_2(speaker,microphone,:),nfft);
    end
end

clear all_impulse_responses_time_domain_zone_2;
clear all_impulse_responses_time_domain_zone_1;

filter_sound_1=q{1,1};
filter_sound_2=q{2,1};


%% summering för att få sound out (pm) i en mikrofon

%testa fucka upp signalen först
%och det borde typ funka ändå??
bright_mic=5;
dark_mic=7;

output_bright_zone=multiply_with_bins(tone_freq_domain,squeeze(impulse_responses_freq_domain_zone_1(:,bright_mic,:)),filter_sound_1,nfft);
output_dark_zone=multiply_with_bins(tone_freq_domain,squeeze(impulse_responses_freq_domain_zone_2(:,dark_mic,:)),filter_sound_1,nfft);


%%
sound_out=ifft(output_bright_zone, nfft);
sound_out_dark=ifft(output_dark_zone, nfft);

sound(repmat(real(sound_out),[1 4]),fs);
%%
sound(repmat(real(sound_out_dark),[1 4]),fs);

%% Plotta impulssvaren och jämföra
hold on
for i = 1:16
    plot(squeeze(impulse_responses_time_domain_zone_1(i,37,:))) % ser bra ut (lika ut för alla högtalare på varje mikronfon)
end

%%
figure
hold on
plot(0:66.66:8000,abs(filter_sound_1(10,:)),'green')
plot(0:66.66:8000,abs(filter_sound_2(1,:)),'black')
