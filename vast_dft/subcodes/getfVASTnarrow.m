function [q, ctrfilter] = getfVASTnarrow(ctrfilter, Hml, Dm, taroption)
% Independent Kbins frequency bins
taridx = taroption.taridx;

try
    journal_exp_1 = taroption.journal_exp_1;
catch
    journal_exp_1 = false;
end

nzones = size(Hml,1);
[Kbins, number_of_microphones, number_of_loudspeakers] = size(Hml{nzones});

allones = ones(number_of_loudspeakers,1);

q = cellfun(@(x) zeros(number_of_loudspeakers, Kbins), ...
    cell(nzones,1), 'UniformOutput', false);

if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sb')
    minSb = true;
else
    minSb = false;
end

% Since this is the narrowband approach
taroption.broadband = false;
taroption.tightfigure = true;
taroption.nloudspks = number_of_loudspeakers;
if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sd')
    if strcmpi(ctrfilter.cvxopt_properties.const, 'sb')
        taroption.nsb = false;
    else
        taroption.nsb = true;
    end
end
if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sb')
    if strcmpi(ctrfilter.cvxopt_properties.const, 'sd')
        taroption.nsd = false;
    else
        taroption.nsd = true;
    end
end
taroption.frequency = true;

if journal_exp_1
    ff = taridx;
    H = cell(nzones,1);
    ssidx = flipud(perms(1:nzones));
    for sound_region_index = 1:nzones
        H{sound_region_index} = squeeze(Hml{sound_region_index}(ff,:,:));
    end
    
    for sound_region_index = 1:nzones
        Hb = H{ssidx(sound_region_index,1)};
        Hd = H{ssidx(sound_region_index,2)};
        
        d = Dm{sound_region_index}(ff,:).';
        
        Rb = Hb'*Hb;
        Rd = Hd'*Hd;
        rb = Hb'*d;
        
        % Joint diagonalization
        [U, D] = jdiag(Rb, Rd, 'vector', true);
        
        calinfo.mucan = logspace(-15,15,101);
        calinfo.Vcan = 1:16;
        
        calinfo.wRd = Rd;
        calinfo.whz = d;
        calinfo.wrb = rb;
        calinfo.U = U;
        calinfo.D = D;
        calinfo.sridx = sound_region_index;
        calinfo.tarfreq = taroption.tarfreq;
        
        pfmmtx = get_mse_pfm(calinfo);
        q = pfmmtx;
    end

else
    
    for ff = 1:Kbins
        H = cell(nzones,1);
        ssidx = flipud(perms(1:nzones));
        for sound_region_index = 1:nzones
            H{sound_region_index} = squeeze(Hml{sound_region_index}(ff,:,:));
        end

        for sound_region_index = 1:nzones
            Hb = H{ssidx(sound_region_index,1)};
            Hd = H{ssidx(sound_region_index,2)};

            d = Dm{sound_region_index}(ff,:).';

            Rb = Hb'*Hb;
            Rd = Hd'*Hd;
            rb = Hb'*d;

            % Joint diagonalization
            [U, D] = jdiag(Rb, Rd, 'vector');

            V = ctrfilter.V;

            d_part = real(D(1:V));
            U_part = U(:,1:V);
            uvrb = U_part'*rb;

            % Calculate the optimal coefficient a
            if ctrfilter.cvxopt_properties.findopt
                if strcmpi(ctrfilter.cvxopt_properties.opttype, 'min_sd')
                    ctrfilter = getOptPara(uvrb,d_part,d,ctrfilter,[ff, sound_region_index]);
                else
                    iRdi = abs(allones'*Rd*allones);
                    ctrfilter = getOptPara(uvrb,d_part,iRdi,ctrfilter,[ff, sound_region_index]);
                end

                mu_new = ctrfilter.cvxopt_properties.optpara(ff, sound_region_index);
            else
                mu_new = ctrfilter.mu;
            end

            ctrfilter.mu = mu_new;

            qf = U_part*(uvrb./(mu_new + d_part));
            q{sound_region_index}(:,ff) = qf;

%             % Plot the cost functions
%             if sridx == 1 && ff == taridx
%                 costfcnplot(U,D,wrb,abs(whz'*whz),minSb,true,taroption)
%             end

        end

    end
end

end