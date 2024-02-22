function [ctrfilt, general] = get_control_filter(general, output_variables)

% Set convex optimization properties
convex_optimizer_properties.findopt = true;
convex_optimizer_properties.opttype = 'min_sb';  % 'min_sb' , 'min_sd'
if strcmpi(convex_optimizer_properties.opttype, 'min_sb')
    convex_optimizer_properties.const = 'nsd';        % 'nsd', 'nsb', 'sd'
    convex_optimizer_properties.tarval = 10^(-3.5);   % an norm energy, hence, -35 dB
else
    % min_sd
    convex_optimizer_properties.const = 'nsb';
    convex_optimizer_properties.tarval = 10^(-1.5); % an energy, hence, -15 dB
end
convex_optimizer_properties.initpara = zeros(1,output_variables.nzones);
convex_optimizer_properties.optpara = zeros(1,output_variables.nzones);

domaintype = {'frequency', 'time', 'stft'};

conFilter = cellfun(@(x) zeros(output_variables.LJ, 1), ...
    cell(output_variables.nzones,1), 'UniformOutput', false);

% Initialization for VAST V, mu
vast_structure = struct('Vmax', output_variables.LJ, 'mu', 0, 'V', 1, ...
    'conFilter', {conFilter}, 'cvxopt_properties', convex_optimizer_properties, ...
    'system', 'LTI', 'domain', 'time', 'method', '', 'incl_dcnyq', true);

% Make a set of control filters for various methods
% For plotting the performance metrics and exporting figures

jn_no = output_variables.jn_no;

general = rmfield(general, 'idx');
switch jn_no
    case 1
    case 2
        % Index for each method
        general.idx.vast_nf = 1;
        general.idx.vast_bf = 2;
        general.idx.vast_t = 3;
        general.idx.pm = 4;
        general.idx.accpm = 5;
        general.idx.acc = 6;
        general.idx.nc = 7;
        general.groupnames = {'VAST_NF', 'VAST_BF', 'VAST_T', ...
            'PM', 'ACCPM', 'ACC', 'NC'}';
        idxdomain = [1 1 2 1 1 1 1];
        broadband = [false true true false false false false];
    otherwise
end
idxnames = fieldnames(general.idx);

general.nmethods = length(general.groupnames);

general.filenames_ctr = append(string(general.groupnames), ...
    repmat("_c_",length(general.groupnames),1));

ctrfilt = cellfun(@(x) vast_structure, cell(general.nmethods,1), 'UniformOutput', false);

for nameidx = 1:general.nmethods
    ctrfilt{nameidx}.method = idxnames{nameidx};
    ctrfilt{nameidx}.domain = domaintype{idxdomain(nameidx)};
    if strcmpi(ctrfilt{nameidx}.domain, 'time')
        ctrfilt{nameidx}.Vmax = output_variables.LJ;
        % Target Sd is -20 dB
        ctrfilt{nameidx}.cvxopt_properties.tarval = 10^(-2.0);
        ctrfilt{nameidx}.incl_dcnyq = false;
    elseif strcmpi(ctrfilt{nameidx}.domain,'frequency') && broadband(nameidx)
        ctrfilt{nameidx}.Vmax = output_variables.nloudspks*output_variables.Kbins;
    else % Narrowband
        ctrfilt{nameidx}.Vmax = output_variables.nloudspks;
    end
    ctrfilt{nameidx}.V = ctrfilt{nameidx}.Vmax;
    ctrfilt{nameidx}.broadband = broadband(nameidx);
    if ~broadband(nameidx)
        ctrfilt{nameidx}.cvxopt_properties.optpara = ...
            zeros(output_variables.Kbins,output_variables.nzones);
    end
    ctrfilt{nameidx} = orderfields(ctrfilt{nameidx});
end

ctrfilt{general.idx.accpm}.mu = 0.8;
ctrfilt{general.idx.pm}.mu = 1;

end