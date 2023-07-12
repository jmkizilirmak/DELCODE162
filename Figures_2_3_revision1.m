% DELCODE-SVM: Figures 2 & 3
% _
% written by Joram Soch <Joram.Soch@DZNE.de>, 13/07/2022, 19:49
% adaptations by Jasmin Kizilirmak, 2022-07-20; last edit 2023-06-20 to add
% GDS-GAI SVM

clear
close all

%%% Step 1: specify plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root_dir = 'C:\Users\kizilirmakj\ownCloud\DELCODE_JK\DELCODE_162_Personality\SVM_PerAF\revision1_additionalSVMs';

% specify results directories
res_dirs = {fullfile(root_dir,'SVMs_same_subjs_noADrel\'), ...
            fullfile(root_dir,'SVMs_same_subjs_CSFxtra_noADrel_no_DMN\'), ...
            fullfile(root_dir,'SVMs_same_subjs_CSFxtra_riskclass_noADrel_no_DMN\')};
%            fullfile(root_dir,'SVMs_same_subjs_riskclass_noADrel\'), ...
        
% specify analysis folders
ana_folds = {'008_SVM_counfounds_only', ...
             '002_SVM_NI_only', ...
             '003_SVM_BFI_only', ...
             '004_SVM_GDS_and_GAI', ...
             '005_SVM_BFIext_only', ...
             '006_SVM_ApoE_only', ...
             '007_SVM_CSF_only', ...
             '001_SVM_allinone'};

% specify analysis names
long_titles  = {'base model [age, gender, site]', ...
                'voxel-wise resting-state mPerAF (DMN)', ...
                'personality [BFI-10]', ...
                'depression, anxiety [GDS, GAI-SF]', ...
                'personality extended [BFI-10, GDS, GAI-SF]', ...
                'ApoE [E4 presence]', ...
                'CSF biomarkers [tTau, pTau-187, Abeta-42/40]', ...
                'all-in-one'};
short_titles = {'base model', 'mPerAF', 'personality', 'depr. & anx.', 'personality ext.', 'ApoE', 'CSF', 'all-in-one'};
class_names  = {{'HC', 'SCD', 'aMCI', 'AD'},{'HC', 'at-risk for AD', 'AD'}};

% create figures for all analysis variants (-> res_dirs)
for k=1:numel(res_dirs)
    
    if k > 1, ana_folds{8} = '001_SVM_all_wo_CSF'; long_titles{8} = 'all-in-one without CSF'; short_titles{8} = 'all w/o CSF';  end
    %if mod(k,2)==1, cnames = class_names{1}; else cnames = class_names{2}; end
    if k~=3, cnames = class_names{1}; else cnames = class_names{2}; end
    
    %%% Step 2: load results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get analysis dimensions
    load(strcat(res_dirs{k},ana_folds{1},'/','SVM.mat'));
    num_anas  = numel(ana_folds);
    num_subs  = SVM.pars.subs;
    num_class = SVM.data.m;
    clear SVM

    % preallocate results
    DAs   = zeros(num_subs,num_anas);
    DA_CI = zeros(2,num_anas);
    CAs   = cell(1,num_anas);
    CA_CI = cell(1,num_anas);

    % load all SVM analyses
    for i = 1:num_anas

        % load SVM.mat
        load(strcat(res_dirs{k},ana_folds{i},'/','SVM.mat'));
        CAs{i}   = zeros(num_subs,num_class);
        CA_CI{i} = zeros(2,num_class);

        % get accuracies
        DAs(:,i)   = reshape(SVM.perf.DA(1,1,:), [num_subs 1]);
        DA_CI(:,i) = mean(SVM.perf.DA_CI,3)';
        for j = 1:num_class
            CAs{i}(:,j)   = reshape(SVM.perf.CA(j,1,:), [num_subs 1]);
            CA_CI{i}(:,j) = mean(SVM.perf.CA_CI(j,:,:),3)';
        end

    end


    %%% Figure 2: DA vs. feature set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % open figure
    figure('Name', 'DELCODE-SVM: Figure 2', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
    fc = [0.42,0.76,0.98];  % face color bars
    ec = [0.05,0.50,0.79];  % edge color bars

    % plot DAs
    hold on;
    DA = mean(DAs);
    bar([1:num_anas], DA, 'FaceColor', fc, 'EdgeColor', ec, 'LineWidth', 1.5);
    errorbar([1:num_anas], DA, DA-DA_CI(1,:), DA_CI(2,:)-DA, '.k', 'LineWidth', 2, 'CapSize', 16, 'Color', ec);
    plot([0, num_anas+1], [1/num_class, 1/num_class], ':k', 'LineWidth', 2);
    if num_subs > 0
        rng(size(DAs,1));
        for i = 1:size(DAs,2)
            plot(i+3/4*(rand(size(DAs,1),1)-1/2), DAs(:,i), '.k', 'LineWidth', 2, 'MarkerSize', 10);
        end
    end
    axis([0, num_anas+1, 0, 2/3]);
    set(gca,'Box','On');
    set(gca,'FontSize',16);
    set(gca,'XTick',[1:num_anas],'XTickLabel',short_titles,'XTickLabelRotation',0);
    legend('decoding accuracy', '90% confidence intervals', 'chance level', 'Location', 'NorthWest', 'EdgeColor', 'none');
    xlabel('feature set', 'FontSize', 20);
    ylabel('predictive accuracy', 'FontSize', 20);
    %title('all feature sets', 'FontSize', 22);

    % save figure
    saveas (gcf, [res_dirs{k},'/Decoding_Accuracies.fig']); 
    F = getframe(gcf);
    orient landscape;
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [42.33 23.81]);
    print ('-dpdf', [res_dirs{k},'/Decoding_Accuracies.pdf'],'-bestfit');
    close; % close figure

    %%% Figure 3: CAs vs. feature set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fc = [0.30,0.65,0.93];  % face color bars
    %ec = [0.05,0.50,0.79];  % edge color bars
    
    % open figure
    figure('Name', 'DELCODE-SVM: Figure 3', 'Color', [1 1 1], 'Position', [50 0 1000 1200]);
    sp  = [4,2];
    sps = [1:8];

    % plot CAs
    for i = 1:num_anas
        subplot(sp(1), sp(2), sps(i));
        hold on;
        CA = mean(CAs{i});
        bar([1:num_class], CA, 'FaceColor', fc, 'EdgeColor', ec, 'LineWidth', 1);
        errorbar([1:num_class], CA, CA-CA_CI{i}(1,:), CA_CI{i}(2,:)-CA, '.k', 'LineWidth', 1, 'CapSize', 10, 'Color', ec);
        plot([0, num_class+1], [1/num_class, 1/num_class], ':k', 'LineWidth', 0.5);
        if num_subs > 0
            rng(size(CAs{i},1));
            for j = 1:size(CAs{i},2)
                plot(j+3/4*(rand(size(CAs{i},1),1)-1/2), CAs{i}(:,j), '.k', 'LineWidth', 1.5, 'MarkerSize', 6);
            end
        end
        axis([0, num_class+1, 0, 0.9]);
        set(gca,'Box','On');
        set(gca,'FontSize',10);
        set(gca,'XTick',[1:num_class],'XTickLabel', cnames,'XTickLabelRotation',0);
        %if i == 1, legend('class accuracies', '90% confidence intervals', 'chance level', 'Location', 'NorthWest', 'EdgeColor', 'none', 'Color', 'none'); end
        xlabel('participant group', 'FontSize', 12);
        ylabel('predictive accuracy', 'FontSize', 12);
        title(long_titles{i}, 'FontSize', 12);
    end
    saveas (gcf, [res_dirs{k},'/Class_Accuracies.fig']); 
    F = getframe(gcf);
    orient portrait;
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [23.81 31]);
    print ('-dpdf', [res_dirs{k},'/Class_Accuracies.pdf'],'-bestfit');
    close; % close figure
end