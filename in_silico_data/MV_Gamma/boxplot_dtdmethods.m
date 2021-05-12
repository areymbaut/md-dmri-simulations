function boxplot_dtdmethods(system_name_list, method_names, param_names, method_names_plot, param_names_plot, xps, SNR, Nbs, path, labels_system, output_path, path_pdf_name, method_color, char_labels_system)

path_to_pdf = strcat(output_path, path_pdf_name);
saved_metrics_name = strrep(path_pdf_name,"boxplot_","");
saved_metrics_name = strrep(saved_metrics_name,".pdf","");
xps_pa = mdm_xps_pa(xps);

set(0, 'defaultLegendInterpreter','latex');

if nargin < 14
   char_labels_system = ''; 
end

Nparam = numel(param_names);
Nmethod = numel(method_names);
nb_chunks = length(system_name_list);

zeroarray_true = zeros(1,nb_chunks);
zeroarray_inversion = zeros(Nbs,nb_chunks);
zeroarray_true_signal = zeros(xps_pa.n,nb_chunks);
zeroarray_signal_inversion = zeros(xps_pa.n,Nbs,nb_chunks);
for nparam = 1:Nparam
    eval([param_names{nparam} '_true = zeroarray_true;'])
    if nparam == 1
        signal_true = zeroarray_true_signal;
    end
    for nmethod = 1:Nmethod
        eval([param_names{nparam} '_' method_names{nmethod} ' = zeroarray_inversion;'])
        if nparam == 1
            eval(['signal_' method_names{nmethod} ' = zeroarray_signal_inversion;'])
        end
    end
end

if strcmp(char_labels_system,'MeanDiso') && strcmp(char_labels_system,'MeanDaniso') && strcmp(char_labels_system,'Viso') && strcmp(char_labels_system,'VarDaniso')
    labels_system = {};
end

param_names_0 = param_names;
Nparam_0 = Nparam;

for it_chunks = 1:nb_chunks
    param_names = param_names_0;
    Nparam = Nparam_0;
    system_name = system_name_list(it_chunks);
    structure_metrics = get_metrics_ground_truth_noisy_inversions(system_name, method_names, param_names, xps, SNR, Nbs, output_path, path);
    save(strcat(output_path, '/', saved_metrics_name, '_system', num2str(it_chunks), '.mat'),'-struct','structure_metrics');
    
    if strcmp(char_labels_system,'MeanDiso')
        new_label_system = structure_metrics.mdiso_true;
        labels_system{it_chunks} = ['$\mathrm{E}_{\mathrm{I}}=' num2str(round(new_label_system*1e9,2)) '$'];
    end
    
    if strcmp(char_labels_system,'MeanDaniso')
        new_label_system = structure_metrics.msdanison_true;
        labels_system{it_chunks} = ['$\tilde{\mathrm{E}}_{\mathrm{A}} =' num2str(round(new_label_system,2)) '$'];
    end
    
    if strcmp(char_labels_system,'Viso')
        new_label_system = structure_metrics.vdiso_true;
        labels_system{it_chunks} = ['$\mathrm{V}_{\mathrm{I}} =' num2str(round(new_label_system*1e18,2)) '$'];
    end
    
    if strcmp(char_labels_system,'VarDaniso')
        new_label_system = structure_metrics.vsdaniso_true;
        labels_system{it_chunks} = ['$\mathrm{V}_{\mathrm{A}} =' num2str(round(new_label_system*1e18,2)) '$'];
    end
    
    for nparam = 1:Nparam
        eval([param_names{nparam} '_true(it_chunks) = structure_metrics.' param_names{nparam} '_true;'])
        if nparam == 1
            signal_true(:,it_chunks) = structure_metrics.signal_pa_true;
        end
        if strcmp(param_names{nparam},'vsdaniso') || strcmp(param_names{nparam},'vsdanison')
            for nmethod = 1:Nmethod_Vaniso
                eval([param_names{nparam} '_' method_names_Vaniso{nmethod} '(:,it_chunks) = structure_metrics.' param_names{nparam} '_' method_names_Vaniso{nmethod} ';'])
            end
        else
            for nmethod = 1:Nmethod
                eval([param_names{nparam} '_' method_names{nmethod} '(:,it_chunks) = structure_metrics.' param_names{nparam} '_' method_names{nmethod} ';'])
                if nparam == 1
                    eval(['signal_' method_names{nmethod} '(:,:,it_chunks) = structure_metrics.fitted_signal_pa_' method_names{nmethod} ';'])
                end
            end
        end
    end
end

for nparam = 1:Nparam
    if strcmp(param_names{nparam},'vsdaniso') || strcmp(param_names{nparam},'vsdanison')
        eval(['grouped_' param_names{nparam} ' = cell([1,Nmethod_Vaniso]);'])
        for nmethod = 1:Nmethod_Vaniso
            eval(['grouped_' param_names{nparam} '{nmethod} = ' param_names{nparam} '_' method_names_Vaniso{nmethod} ';'])
        end
    else
        eval(['grouped_' param_names{nparam} ' = cell([1,Nmethod]);'])
        for nmethod = 1:Nmethod
            eval(['grouped_' param_names{nparam} '{nmethod} = ' param_names{nparam} '_' method_names{nmethod} ';'])
        end
    end
end

%% BOXPLOTS
lw = 2;
width = 0.5;
height = 0.2;
delta = linspace(-.3,.3,Nmethod); % Define offsets to distinguish plots
box_width = .2; % Small width to avoid overlap
x = (1-width)/2;
inter_v = 0.03;
inter_v_top_bottom = (1 - 2*inter_v - 3*height)/2;
y1 = inter_v_top_bottom;
y2 = y1 + inter_v + height;
y3 = y2 + inter_v + height;

f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off'); 
f.PaperOrientation = 'landscape';
f.PaperUnits = 'normalized';
f.PaperPosition = [0 0 1 1];

axh_mdiso = axes('position',[x y3 width height]);
axh_msdanison = axes('position',[x y2 width height]);
axh_vdiso = axes('position',[x y1 width height]);

for nparam = 1:Nparam
    
    eval(['axh = axh_' param_names{nparam} ';'])
    hold(axh,'on')
    
    if strcmp(param_names{nparam},'mdiso')
        correct_unit = 1e9;
    elseif strcmp(param_names{nparam},'vdiso')
        correct_unit = 1e18;
    elseif strcmp(param_names{nparam},'vsdaniso')
        correct_unit = 1e36;
    else
        correct_unit = 1;
    end
    
    for it_chunks = 1:nb_chunks
        eval(['xdat = ' param_names{nparam} '_true(it_chunks);'])
        plot(axh,[it_chunks+delta(1) it_chunks+delta(Nmethod)],correct_unit*xdat*[1 1],'-k','LineWidth',1)
    end
    
    eval(['GroupedData = grouped_' param_names{nparam} ';'])
    potential_max = [];
    potential_min = [];
    
    if strcmp(param_names{nparam},'vsdaniso')
        
        for nmethod=1:Nmethod_Vaniso
            if nparam == Nparam
                boxplot(axh, correct_unit*GroupedData{nmethod},'Color', method_color_Vaniso{nmethod}, 'Whisker', 1.5, 'boxstyle','filled', 'Symbol', '', 'OutlierSize', 3, 'position',(1:nb_chunks)+delta(nmethod), 'widths',box_width, 'labels',labels_system);
                %bp = gca;
                axh.XAxis.TickLabelInterpreter = 'latex';
                axh.YAxis.TickLabelInterpreter = 'latex';
                set(findobj(axh,'tag','Median'),'linewidth',1);
                set(findobj(axh,'tag','Whisker'),'linewidth',1);
            else
                boxplot(axh, correct_unit*GroupedData{nmethod},'Color', method_color_Vaniso{nmethod}, 'Whisker', 1.5, 'boxstyle','filled', 'Symbol', '', 'OutlierSize', 3, 'position',(1:nb_chunks)+delta(nmethod), 'widths',box_width);
                axh.YAxis.TickLabelInterpreter = 'latex';
                set(findobj(axh,'tag','Median'),'linewidth',1);
                set(findobj(axh,'tag','Whisker'),'linewidth',1);
                set(axh,'XTickLabel',{' '})
            end
            h = findobj(axh,'tag','Whisker');
            h_cell = num2cell(h);
            for it = 1:length(h_cell)
                y = h_cell{it}.YData;
                potential_max = [potential_max max(y)];
                potential_min = [potential_min min(y)];
            end
        end
    else
        for nmethod=1:Nmethod
            if nparam == Nparam
                boxplot(axh, correct_unit*GroupedData{nmethod},'Color', method_color{nmethod}, 'Whisker', 1.5, 'boxstyle','filled', 'Symbol', '', 'OutlierSize', 3, 'position',(1:nb_chunks)+delta(nmethod), 'widths',box_width, 'labels',labels_system);
                %bp = gca;
                axh.XAxis.TickLabelInterpreter = 'latex';
                axh.YAxis.TickLabelInterpreter = 'latex';
                set(findobj(axh,'tag','Median'),'linewidth',1);
                set(findobj(axh,'tag','Whisker'),'linewidth',1);
            else
                boxplot(axh, correct_unit*GroupedData{nmethod},'Color', method_color{nmethod}, 'Whisker', 1.5, 'boxstyle','filled', 'Symbol', '', 'OutlierSize', 3, 'position',(1:nb_chunks)+delta(nmethod), 'widths',box_width);
                axh.YAxis.TickLabelInterpreter = 'latex';
                set(findobj(axh,'tag','Median'),'linewidth',1);
                set(findobj(axh,'tag','Whisker'),'linewidth',1);
                set(axh,'XTickLabel',{' '})
            end
            h = findobj(axh,'tag','Whisker');
            h_cell = num2cell(h);
            for it = 1:length(h_cell)
                y = h_cell{it}.YData;
                potential_max = [potential_max max(y)];
                potential_min = [potential_min min(y)];
            end
            L(nmethod) = plot(axh, NaN,1,'color',method_color{nmethod}); %// dummy plot for legend
        end
    end
    box(axh,'off')
    grid(axh,'on')
    
    set(axh, 'Fontsize', 14)
    
    ylab = ylabel(axh,param_names_plot{nparam}, 'Fontsize', 16, 'Interpreter', 'latex');
    set(ylab, 'Units', 'Normalized', 'Position', [-0.07, 0.5, 0]);
    xlim(axh, [1+2*delta(1) nb_chunks+2*delta(Nmethod)])

    y_min = min(potential_min);
    y_max = max(potential_max);
    ylim(axh, [y_min y_max])
   
    set(axh, 'LineWidth', lw)
    set(axh, 'TickDir','out');
    if nparam==1
        [h_legend, hObj] = legend(axh, L, method_names_plot);
        pos = get(h_legend,'Position');
        set(h_legend,'Position', [pos(1)+1.2*pos(3), pos(2:4)]);
        hl = findobj(hObj,'type','line');
        set(hl,'LineWidth',1.5);
    end
end
saveas(gcf,path_to_pdf)
clf(f)
clear L














%% HISTOGRAMS
lw = 2;
width = 0.4;
height = 0.15;
delta = linspace(-.3,.3,Nmethod); % Define offsets to distinguish plots
box_width = .2; % Small width to avoid overlap
x = (1-width)/2;
inter_v = 0.06;
inter_v_top_bottom = (1 - 2*inter_v - 3*height)/2;
y1 = inter_v_top_bottom;
y2 = y1 + inter_v + height;
y3 = y2 + inter_v + height;

nbins = 50;

for nmethod=1:Nmethod
    for it_chunks = 1:nb_chunks
        
        path_to_pdf = strcat(output_path, '/', saved_metrics_name, '_hist_system', num2str(it_chunks), '_', method_names{nmethod}, '.pdf');
        
        f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off');
        f.PaperOrientation = 'landscape';
        f.PaperUnits = 'normalized';
        f.PaperPosition = [0 0 1 1];
        
        axh_mdiso = axes('position',[x y3 width height]);
        axh_msdanison = axes('position',[x y2 width height]);
        axh_vdiso = axes('position',[x y1 width height]);
        
        for nparam = 1:Nparam
            eval(['axh = axh_' param_names{nparam} ';'])
            hold(axh,'on')
            
            if strcmp(param_names{nparam},'mdiso')
                correct_unit = 1e9;
            elseif strcmp(param_names{nparam},'vdiso')
                correct_unit = 1e18;
            elseif strcmp(param_names{nparam},'vsdaniso')
                correct_unit = 1e36;
            else
                correct_unit = 1;
            end
            
            eval(['xdat = ' param_names{nparam} '_true(it_chunks);'])
            plot(axh,correct_unit*xdat*[1 1], [0 1],'-','LineWidth',2, 'Color', [0.75, 0, 0.75])
            
            eval(['Data = grouped_' param_names{nparam} '{nmethod};'])
            
            histogram(axh, correct_unit*Data(Nbs*(it_chunks-1)+1:Nbs*it_chunks), nbins, 'Normalization', 'probability', 'FaceColor', method_color{nmethod}, 'EdgeColor', method_color{nmethod});
            plot(axh,median(correct_unit*Data(Nbs*(it_chunks-1)+1:Nbs*it_chunks))*[1 1], [0 1],'--','LineWidth',2, 'Color', method_color{nmethod})
            axh.XAxis.TickLabelInterpreter = 'latex';
            axh.YAxis.TickLabelInterpreter = 'latex';
            
            box(axh,'off')
            grid(axh,'off')
            
            set(axh, 'Fontsize', 14)
            
            ylab = ylabel(axh,param_names_plot{nparam}, 'Fontsize', 16, 'Interpreter', 'latex');
            set(ylab, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
            %             xlim(axh, [1+2*delta(1) nb_chunks+2*delta(Nmethod)])
            
            %             y_min = min(potential_min);
            %             y_max = max(potential_max);
            %             ylim(axh, [y_min y_max])
            
            set(axh, 'LineWidth', lw)
            set(axh, 'TickDir','out');
            %             if nparam==1
            %                 [h_legend, hObj] = legend(axh, L, method_names_plot);
            %                 pos = get(h_legend,'Position');
            %                 set(h_legend,'Position', [pos(1)+1.2*pos(3), pos(2:4)]);
            %                 hl = findobj(hObj,'type','line');
            %                 set(hl,'LineWidth',1.5);
            %             end
        end
        
        saveas(gcf,path_to_pdf)
        clf(f)
        
    end
end















%% SIGNAL PLOTS
path_to_pdf = strcat(output_path, strrep(path_pdf_name,"boxplot","signal"));

% Detect the relevant encodings
potential_encoding_types = {'linear' ; 'spherical' ; 'planar'};
potential_encoding_subscripts = {'lin.' ; 'sph.' ; 'pla.'};
ind_linear = round(xps_pa.b_delta) == 1;
ind_spherical = round(xps_pa.b_delta) == 0;
ind_planar = round(xps_pa.b_delta,1) == -0.5;
is_encoding = [nnz(ind_linear),nnz(ind_spherical),nnz(ind_planar)];
nb_encodings = nnz(is_encoding);

encoding_types = cell([nb_encodings,1]);
encoding_subscripts = cell([nb_encodings,1]);
count = 0;
for it_encoding = 1:length(potential_encoding_types)
    if is_encoding(it_encoding) ~= 0
        count = count + 1;
        encoding_types{count} = potential_encoding_types{it_encoding};
        encoding_subscripts{count} = potential_encoding_subscripts{it_encoding};
    end
end

% Correct the powder-averaged xps
max_b = [];
encoding_added_b0 = [];
pa_weights = xps_pa.pa_w;
for it_encoding = 1:nb_encodings
    encoding = encoding_types{it_encoding};

    eval(['b_' encoding ' = xps_pa.b(ind_' encoding ');'])
    eval(['pa_weights_' encoding ' = pa_weights(ind_' encoding ');'])
    eval(['pa_weights_' encoding ' = pa_weights_' encoding '/sum(pa_weights_' encoding ');'])
    eval(['[b_' encoding ', ind_b_' encoding '_sorted] = sort(b_' encoding ');'])
    eval(['pa_weights_' encoding ' = pa_weights_' encoding '(ind_b_' encoding '_sorted);'])    
    
    eval(['unique_b_' encoding '_ind = {};'])
    eval(['unique_b_' encoding ' = [];'])
    
    count = 0;
    it = 1;
    eval(['b_value_list = b_' encoding ';'])
    eval(['weight_list = pa_weights_' encoding ';'])
    
    while it <= length(b_value_list)-1
        if abs(b_value_list(it+1)-b_value_list(it)) < 0.01*1e9
            count = count + 1;
            temporary_ind_list = [];
            temporary_b_list = [];
            while abs(b_value_list(it+1)-b_value_list(it)) < 0.01*1e9
                if ismember(it, temporary_ind_list)
                    temporary_ind_list = [temporary_ind_list it+1];
                    temporary_b_list = [temporary_b_list b_value_list(it+1)];
                else
                    temporary_ind_list = [temporary_ind_list it it+1];
                    temporary_b_list = [temporary_b_list b_value_list(it) b_value_list(it+1)];
                end
                it = it + 1;
                if it > length(b_value_list)-1
                    break
                end
            end
            unique_b_list = temporary_b_list';
            unique_ind_list = temporary_ind_list;
            unique_w_list = weight_list(unique_ind_list);
            new_b = sum(unique_w_list.*unique_b_list)/sum(unique_w_list);
            eval(['unique_b_' encoding '_ind{count} = unique_ind_list;'])
            eval(['unique_b_' encoding ' = [unique_b_' encoding ' ; new_b];'])   
            it = it + 1;
        else
            count = count + 1;
            eval(['unique_b_' encoding '_ind{count} = it;'])
            eval(['unique_b_' encoding ' = [unique_b_' encoding ' ; b_value_list(it)];']) 
            if it == length(b_value_list)-1
                count = count + 1;
                eval(['unique_b_' encoding '_ind{count} = it+1;'])
                eval(['unique_b_' encoding ' = [unique_b_' encoding ' ; b_value_list(it+1)];'])
            end
            it = it + 1;
        end
    end  
    
    eval(['max_b = [max_b max(unique_b_' encoding ')];'])
    
    eval(['b_vals = unique_b_' encoding ';'])
    if ~ismember(0,b_vals)
        eval(['unique_b_' encoding ' = [0 ; unique_b_' encoding '];'])
        encoding_added_b0 = [encoding_added_b0 string(encoding)];
    end
    
end
max_b = max(max_b);

lw = 2;  
height = 0.15;
width = 0.75*height;
inter_h = 0.0375;
inter_v = 0.0375;
width_median = 0.17/2;

lwm = 4;
width_method = [lwm, 0.77*lwm, 0.62*lwm, 0.45*lwm];

y_spherical = (1 - nb_encodings*height - (nb_encodings-1)*inter_v)/2;
y_planar = y_spherical + (nb_encodings-2)*(height + inter_v);
y_linear = y_spherical + (nb_encodings-1)*(height + inter_v);

f = figure('Position', [0,0,900,706], 'Units', 'pixels', 'visible', 'off'); 
f.PaperOrientation = 'landscape';
f.PaperUnits = 'normalized';
f.PaperPosition = [0 0 1 1];

min_S = 1;
min_MD_b = 1;

FLAG_removal_true = 1;
FLAG_removal_inversion = 1;

for it_chunks = 1:nb_chunks
    x = (1-nb_chunks*width-(nb_chunks-1)*inter_h)/2 + (it_chunks-1)*(width+inter_h);
    MD = mdiso_true(it_chunks);
    
    for it_encoding = 1:nb_encodings
        encoding = encoding_types{it_encoding};
        
        eval(strcat('axh_', num2str(it_chunks), '_', encoding, ' = axes(', '''position''', ',[x y_', encoding, ' width height]);'))
        eval(['axh = axh_' num2str(it_chunks) '_' encoding ';'])
           
        % Ground truth signal and inversion signals
        eval(['signal_' encoding '_true = signal_true(ind_' encoding ',it_chunks);']) 
        eval(['signal_' encoding '_true = signal_' encoding '_true(ind_b_' encoding '_sorted);'])
        for nmethod = 1:Nmethod        
            eval(['signal_' method_names{nmethod} '_' encoding '_inversion = squeeze(signal_' method_names{nmethod} '(ind_' encoding ',:,it_chunks));'])
            eval(['signal_' method_names{nmethod} '_' encoding '_inversion = signal_' method_names{nmethod} '_' encoding '_inversion(ind_b_' encoding '_sorted,:);'])
        end
        
        eval(['N = length(unique_b_' encoding '_ind);'])
        unique_signal_true = zeros([1,N]);
        for nmethod = 1:Nmethod 
            eval(['unique_signal_' method_names{nmethod} '_inversion = zeros([Nbs,N]);'])
            eval(['unique_median_signal_' method_names{nmethod} '_inversion = zeros([1,N]);'])
            eval(['unique_q1_signal_' method_names{nmethod} '_inversion = zeros([1,N]);'])
            eval(['unique_q3_signal_' method_names{nmethod} '_inversion = zeros([1,N]);'])
        end
        
        for c = 1:N
            eval(['weights = pa_weights_' encoding '(unique_b_' encoding '_ind{c});'])
            eval(['unique_signal_true(c) = sum(weights.*signal_' encoding '_true(unique_b_' encoding '_ind{c}))/sum(weights);'])        
            for nmethod = 1:Nmethod 
                for m =1:Nbs
                    eval(['unique_signal_' method_names{nmethod} '_inversion(m,c) = sum(weights.*signal_' method_names{nmethod} '_' encoding '_inversion(unique_b_' encoding '_ind{c},m))/sum(weights);'])
                end
                eval(['quartiles = quantile(unique_signal_' method_names{nmethod} '_inversion(:,c),3);'])
                eval(['unique_median_signal_' method_names{nmethod} '_inversion(c) = quartiles(2);'])
                eval(['unique_q1_signal_' method_names{nmethod} '_inversion(c) = quartiles(1);'])
                eval(['unique_q3_signal_' method_names{nmethod} '_inversion(c) = quartiles(3);'])
            end
        end
        
        % Update the .mat file with the good powder-averaged signals
        structure_metrics = load(strcat(output_path, '/', saved_metrics_name, '_system', num2str(it_chunks), '.mat'));
        if FLAG_removal_true
            structure_metrics = rmfield(structure_metrics, 'signal_pa_true');
            FLAG_removal_true = 0;
        end
        eval(['structure_metrics.b_' encoding ' = unique_b_' encoding ';'])
        eval(strcat('structure_metrics.signal_pa_true_', encoding, ' = unique_signal_true;'))
        for nmethod = 1:Nmethod 
            if FLAG_removal_inversion
                eval(strcat('structure_metrics = rmfield(structure_metrics,', '''fitted_signal_pa_', method_names{nmethod}, '''', ');'))
            end
            eval(strcat('structure_metrics.signal_pa_', method_names{nmethod}, '_', encoding, ' = unique_signal_', method_names{nmethod}, '_inversion;'))
        end
        FLAG_removal_inversion = 0;
        save(strcat(output_path, '/', saved_metrics_name, '_system', num2str(it_chunks), '.mat'),'-struct','structure_metrics');
    
        eval(['b_vals = unique_b_' encoding ';'])
        if ismember(string(encoding),encoding_added_b0)
            unique_signal_true = [1 unique_signal_true];
            for nmethod = 1:Nmethod
                eval(['unique_median_signal_' method_names{nmethod} '_inversion = [1 unique_median_signal_' method_names{nmethod} '_inversion];'])
                eval(['unique_q1_signal_' method_names{nmethod} '_inversion = [1 unique_q1_signal_' method_names{nmethod} '_inversion];'])
                eval(['unique_q3_signal_' method_names{nmethod} '_inversion = [1 unique_q3_signal_' method_names{nmethod} '_inversion];'])
            end
        end
        
        unique_signal_true = unique_signal_true'; 
        for nmethod = 1:Nmethod
            eval(strcat('unique_median_signal_', method_names{nmethod}, '_inversion = unique_median_signal_', method_names{nmethod}, '_inversion', '''', ';'))
            eval(strcat('unique_q1_signal_', method_names{nmethod}, '_inversion = unique_q1_signal_', method_names{nmethod}, '_inversion', '''', ';'))
            eval(strcat('unique_q3_signal_', method_names{nmethod}, '_inversion = unique_q3_signal_', method_names{nmethod}, '_inversion', '''', ';'))
        end
        
        if min(unique_signal_true) < min_S
            min_S = min(unique_signal_true);
        end
        
        x_fit = linspace(0, 1e-9*max_b, 150);
        eval(strcat('signal_true_fit = interp1(1e-9*unique_b_', encoding, ', log(unique_signal_true), x_fit,', '''spline''', ');'))
        
        eval(strcat('plot(axh, 1e-9*unique_b_', encoding, ', exp(-MD*unique_b_', encoding, '), ', '''--''', ', ', '''Color''', ', ', string('[0.75,0,0.75]'), ', ', '''Linewidth''', ', 3*lw/4)')) 
        hold on
        plot(x_fit, exp(signal_true_fit), 'k-', 'Linewidth', 3*lw/4)
        eval(strcat('scatter(axh, 1e-9*unique_b_', encoding, ',unique_signal_true, 25,', '''k''', ')'))
%         eval(strcat('scatter(axh, 1e-9*unique_b_', encoding, ',unique_signal_true, 30,', '''k''', ', ', '''filled''', ')'))
        
        eval(['min_MD_b_maybe = exp(-MD*max(unique_b_' encoding '));'])
        
        if min_MD_b_maybe < min_MD_b
            min_MD_b = min_MD_b_maybe;
        end

        eval(['bvals = unique_b_' encoding ';'])
        for nmethod = 1:Nmethod
            %eval(strcat('scatter(axh, 1e-9*unique_b_', encoding, ', unique_median_signal_', method_names{nmethod}, '_inversion, 20,', string('method_color{nmethod}'), ', ', '''filled''', ')'))
            
            if it_chunks == nb_chunks && strcmp(encoding,'linear')
                L(nmethod) = plot(axh, NaN,1,'color',method_color{nmethod}); %// dummy plot for legend
            end
            
            for i = 1:length(bvals)
                if bvals(i) ~= 0
                    eval(strcat('plot(axh, [1e-9*unique_b_', encoding, '(i)-width_median,1e-9*unique_b_', encoding, '(i)+width_median],[unique_median_signal_', method_names{nmethod}, '_inversion(i),unique_median_signal_', method_names{nmethod}, '_inversion(i)],', '''-''', ', ', '''Linewidth''', ', ', num2str(lw/2), ', ', '''Color''', ', ', string('method_color{nmethod}'), ')'))
                    eval(strcat('plot(axh, [1e-9*unique_b_', encoding, '(i),1e-9*unique_b_', encoding, '(i)],[unique_q1_signal_', method_names{nmethod}, '_inversion(i),unique_q3_signal_', method_names{nmethod}, '_inversion(i)],', '''-''', ', ', '''Linewidth''', ', ', num2str(width_method(nmethod)), ', ', '''Color''', ', ', string('method_color{nmethod}'), ')'))
                end
            end
        end
        
        set(axh, 'Fontsize', 14)
        
        if it_chunks == 1
            ylab = ylabel(axh,['$A_{\mathrm{' encoding_subscripts{it_encoding} '}}(b)$'], 'Fontsize', 16, 'Interpreter', 'latex');
            set(axh,'YTick',[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
            set(axh,'YTickLabel',{'0.1', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '1'})
        else
            set(axh,'YTickLabel',{' '})
        end
        
        if it_encoding == nb_encodings
            if it_chunks == 1
                xlab = xlabel(axh,'$b$ / ms$\cdot\mu$m$^{-2}$', 'Fontsize', 16, 'Interpreter', 'latex');
            end
        else
            tl = title(axh, labels_system{it_chunks}, 'Fontsize', 16, 'Interpreter', 'latex');
            set(tl, 'Units', 'Normalized', 'Position', [0.5, 1.05, 0]);
            set(axh,'XTickLabel',{' '})            
        end
        
       
        box(axh,'off')
        set(axh,'yscale','log')
        set(axh, 'LineWidth', lw)
        set(axh, 'TickLength',[0.03 0.025]);
        set(axh, 'TickDir','out');
        
        xlim(axh, [0 round(1e-9*max_b)+width_median])
        
        if it_chunks == nb_chunks && strcmp(encoding,'linear')
            [h_legend, hObj] = legend(axh, L, method_names_plot);
            pos = get(h_legend,'Position');
            set(h_legend,'Position', [pos(1)+1.25*pos(3), pos(2:4)]);
            hl = findobj(hObj,'type','line');
            set(hl,'LineWidth',1.5);
        end
    end
end

for it_chunks = 1:nb_chunks
    for it_encoding = 1:nb_encodings
        encoding = encoding_types{it_encoding};
        eval(strcat('ylim(axh_', num2str(it_chunks), '_', encoding, ', [min_MD_b/1.5 1])'))
    end
end

saveas(gcf,path_to_pdf)
clf(f)
