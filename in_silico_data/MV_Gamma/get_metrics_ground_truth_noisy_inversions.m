function structure = get_metrics_ground_truth_noisy_inversions(system_name, method_names, param_names, xps, SNR, Nbs, run_time_path, path)

if nargin < 8
   path = ''; 
end

Nmethod = numel(method_names);
Nparam = numel(param_names);

% Ground truth
[dpar, dperp, theta, phi, w] = choose_simulated_system(system_name, path);
opt = get_opt;
opt.dtd_gamma.do_weight = 1; %%%%%%%%%%%%%%%%%%%%%%

dtd = dtd_par2dist(dpar,dperp,theta,phi,w);
m_true = dtd_dtd2m(dtd,opt);
struct_true = dtd_1d_fit2param(m_true);

for nparam = 1:Nparam
    if strcmp(param_names{nparam},'chisqn')
        structure.chisqn = (1/SNR)^2;
    else
        eval(['structure.' param_names{nparam} '_true = struct_true.' param_names{nparam} ';'])
    end
end

% Store additional ground truth metrics
missing_metrics = {'msddelta' ; 'vsddelta' ; 'vsddeltan' ; 'vsdaniso' ; 'vsdanison' ; 'cvdisosdaniso' ; 'cvdisosddelta'};
for c = 1:length(missing_metrics)
    metric_char = missing_metrics{c};
    if ~any(strcmp(param_names,metric_char))
        eval(['structure.' metric_char '_true = struct_true.' metric_char ';'])
    end   
end

s_true = dtd_1d_fit2data(m_true, xps);
[s_pa_true, xps_pa] = mdm_powder_average_1d(s_true,xps);

% Noisy inversions
run_time = zeros(Nmethod,1);
zeroarray = zeros(Nbs,1);
for nmethod = 1:Nmethod
    eval(['fitted_signal_' method_names{nmethod} ' = zeros(xps_pa.n,Nbs);'])
    for nparam = 1:Nparam
        eval([param_names{nparam} '_' method_names{nmethod} ' = zeroarray;'])
    end
    if strcmp(method_names{nmethod},'dtd')
        if ~any(strcmp(param_names,'vsdanison'))
            vsdanison_dtd = zeroarray;
        end
        if ~any(strcmp(param_names,'vsdaniso'))
            vsdaniso_dtd = zeroarray;
        end
    end
    if strcmp(method_names{nmethod},'dtd_pa')
        if ~any(strcmp(param_names,'vsdanison'))
            vsdanison_dtd_pa = zeroarray;
        end
        if ~any(strcmp(param_names,'vsdaniso'))
            vsdaniso_dtd_pa = zeroarray;
        end
    end
end

for nbs = 1:Nbs
    if strcmp(method_names{nmethod},'dtd_mv_gamma')
       nbs 
    end
    s = sqrt((s_true + 1/SNR*randn([xps.n 1])).^2 + (1/SNR*randn([xps.n 1])).^2); % Adding Rician noise
    
    for nmethod = 1:Nmethod
        method_name = method_names{nmethod};
        
        if strcmp(method_name,'dtd') || strcmp(method_name,'dtd_covariance') || strcmp(method_name,'dtd_mv_gamma')
            [s_pa, ~] = mdm_powder_average_1d(s,xps);
            tic
            m = feval([method_name '_1d_data2fit'],s,xps,opt);
            run_t = toc;
            dps = mio_1d_fit2param(method_name,m);
            s_fit = feval([method_name '_1d_fit2data'],m,xps);
            [s_pa_fit, ~] = mdm_powder_average_1d(s_fit,xps);
            chisqn = msf_chisqn(s,s_fit);
            
            eval(['fitted_signal_' method_names{nmethod} '(:,nbs) = s_pa_fit;'])
            for nparam = 1:Nparam
                if strcmp(param_names{nparam},'chisqn')
                    chisqn_dtd(nbs) = chisqn;
                else
                    eval([param_names{nparam} '_' method_names{nmethod} '(nbs) = dps.' param_names{nparam} ';'])
                end
            end
            
            if strcmp(method_name,'dtd')
                if ~any(strcmp(param_names,'vsdanison'))
                    vsdanison_dtd(nbs) = dps.vsdanison;
                end
                if ~any(strcmp(param_names,'vsdaniso'))
                    vsdaniso_dtd(nbs) = dps.vsdaniso;
                end
            end
        
        elseif strcmp(method_name,'dtd_gamma') || strcmp(method_name,'dtd_pa') 
            [s_pa, ~] = mdm_powder_average_1d(s,xps);
            tic
            m = feval([method_name '_1d_data2fit'],s_pa,xps_pa,opt);
            run_t = toc;
            dps = mio_1d_fit2param(method_name,m);
            s_fit = feval([method_name '_1d_fit2data'],m,xps_pa);
            chisqn = msf_chisqn(s_pa,s_fit);
            
            eval(['fitted_signal_' method_names{nmethod} '(:,nbs) = s_fit;'])
            for nparam = 1:Nparam
                if strcmp(param_names{nparam},'chisqn')
                    chisqn_dtd(nbs) = chisqn;
                else
                    eval([param_names{nparam} '_' method_names{nmethod} '(nbs) = dps.' param_names{nparam} ';'])
                end
            end 
            
            if strcmp(method_name,'dtd_pa')
                if ~any(strcmp(param_names,'vsdanison'))
                    vsdanison_dtd_pa(nbs) = dps.vsdanison;
                end
                if ~any(strcmp(param_names,'vsdaniso'))
                    vsdaniso_dtd_pa(nbs) = dps.vsdaniso;
                end
            end
        end
        
        run_time(nmethod) = run_time(nmethod) + run_t;
    end
end

structure.signal_pa_true = s_pa_true;
for nmethod = 1:Nmethod
    eval(['structure.fitted_signal_pa_' method_names{nmethod} ' = fitted_signal_' method_names{nmethod} ';'])
    for nparam = 1:Nparam    
        eval(['structure.' param_names{nparam} '_' method_names{nmethod} ' = ' param_names{nparam} '_' method_names{nmethod} ';'])
    end
    if strcmp(method_names{nmethod},'dtd')
        if ~any(strcmp(param_names,'vsdanison'))
            structure.vsdanison_dtd = vsdanison_dtd;
        end
        if ~any(strcmp(param_names,'vsdaniso'))
            structure.vsdaniso_dtd = vsdaniso_dtd;
        end
    end
    if strcmp(method_names{nmethod},'dtd_pa')
        if ~any(strcmp(param_names,'vsdanison'))
            structure.vsdanison_dtd_pa = vsdanison_dtd_pa;
        end
        if ~any(strcmp(param_names,'vsdaniso'))
            structure.vsdaniso_dtd_pa = vsdaniso_dtd_pa;
        end
    end
end

run_time = run_time/Nbs*1e3;
txt_to_write = strings(Nmethod,2);
for nmethod = 1:Nmethod
    txt_to_write(nmethod,:) = [string(method_names{nmethod}) , string(run_time(nmethod))];
end

fileID = fopen(strcat(run_time_path, '/run_time_ms_', system_name, '.txt'),'w');
% fprintf(fileID,'%15s %15s %15s %15s\n',method_names{1},method_names{2},method_names{3},method_names{4});
fprintf(fileID,'%15s %12.4f\n',txt_to_write');
fclose(fileID);

