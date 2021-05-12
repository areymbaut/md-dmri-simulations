function structure = perform_inversions(Nbs,method_names,param_names,m_true,xps,opt,SNR)

Nmethod = numel(method_names);
Nparam = numel(param_names);

zeroarray = zeros(Nbs,1);
for nmethod = 1:Nmethod
    for nparam = 1:Nparam
        eval([param_names{nparam} '_' method_names{nmethod} ' = zeroarray;'])
    end
end

for nbs = 1:Nbs
    s_true = dtd_1d_fit2data(m_true, xps);
    s = sqrt((s_true + 1/SNR*randn([xps.n 1])).^2 + (1/SNR*randn([xps.n 1])).^2); % Adding Rician noise
    
    for nmethod = 1:Nmethod
        method_name = method_names{nmethod};
        
        if strcmp(method_name,'dtd') || strcmp(method_name,'dtd_covariance') || strcmp(method_name,'dtd_mv_gamma')
            m = feval([method_name '_1d_data2fit'],s,xps,opt);
            dps = mio_1d_fit2param(method_name,m);
            s_fit = feval([method_name '_1d_fit2data'],m,xps);
            chisqn = msf_chisqn(s,s_fit);
            
            for nparam = 1:Nparam
                if strcmp(param_names{nparam},'chisqn')
                    chisqn_dtd(nbs) = chisqn;
                else
                    eval([param_names{nparam} '_' method_names{nmethod} '(nbs) = dps.' param_names{nparam} ';'])
                end
            end
        
        elseif strcmp(method_name,'dtd_gamma') || strcmp(method_name,'dtd_pa') 
            [s_pa, xps_pa] = mdm_powder_average_1d(s,xps);
            m = feval([method_name '_1d_data2fit'],s_pa,xps_pa,opt);
            dps = mio_1d_fit2param(method_name,m);
            s_fit = feval([method_name '_1d_fit2data'],m,xps_pa);
            chisqn = msf_chisqn(s_pa,s_fit);
            
            for nparam = 1:Nparam
                if strcmp(param_names{nparam},'chisqn')
                    chisqn_dtd(nbs) = chisqn;
                else
                    eval([param_names{nparam} '_' method_names{nmethod} '(nbs) = dps.' param_names{nparam} ';'])
                end
            end 
        end
    end
end

for nmethod = 1:Nmethod
    for nparam = 1:Nparam    
        eval(['structure.' param_names{nparam} '_' method_names{nmethod} ' = ' param_names{nparam} '_' method_names{nmethod} ';'])
    end
end