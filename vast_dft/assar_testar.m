%% Initialization
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

desired_room_impulse_responses = impulse_response_virtual_score;

[ctrfilter, perform, rankcheck,q]=calculatefVAST(general, loudspeaker_array, zones, experiment_1_control_filter, impulse_response_measured, desired_room_impulse_responses, [], 'narrow', experiment_1_target_options);

%%

loudspeaker_filter_zone_1=q{1,1}(1,:);
loudspeaker_filter_zone_2=q{2,1}(1,:);

desired_microphone_1_from_zone_1=desired_room_impulse_responses{1,1}(:,1);
desired_microphone_1_from_zone_2=desired_room_impulse_responses{2,1}(:,1);

figure
plot([desired_microphone_1_from_zone_1, desired_microphone_1_from_zone_2])