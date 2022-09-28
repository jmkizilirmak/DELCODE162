% Wrapper for comparing feature sets of DELCODE-162-personality via
% t tests or rather GLMs that correct for correlations between subsamples,
% which are naturally higher for small classes.
% Uses Joram Soch's ML4ML toolbox (https://github.com/JoramSoch/ML4ML) and
% its function ML_SVC_ttest.
%
% Written by Jasmin.Kizilirmak@dzne.de, 2022-07-26, 13:51
% last edit: 2022-07-27 10:58, added Bonferroni-Holm correction

%% set context variables
% specify results directories
res_dirs = {'C:\Users\kizilirmakj\ownCloud\DELCODE_JK\DELCODE_162_Personality\SVM_PerAF\SVMs_same_subjs_noADrel\', ...
            'C:\Users\kizilirmakj\ownCloud\DELCODE_JK\DELCODE_162_Personality\SVM_PerAF\SVMs_same_subjs_riskclass_noADrel\', ...
            'C:\Users\kizilirmakj\ownCloud\DELCODE_JK\DELCODE_162_Personality\SVM_PerAF\SVMs_same_subjs_CSFxtra_noADrel_no_DMN\', ...
            'C:\Users\kizilirmakj\ownCloud\DELCODE_JK\DELCODE_162_Personality\SVM_PerAF\SVMs_same_subjs_CSFxtra_riskclass_noADrel_no_DMN\'};

% specify subfolders
ana_folds = {'007_SVM_counfounds_only', ...
             '002_SVM_NI_only', ...
             '003_SVM_BFI_only', ...
             '004_SVM_BFIext_only', ...
             '005_SVM_ApoE_only', ...
             '006_SVM_CSF_only', ...
             '001_SVM_allinone'};        

%% run statistical comparisons (Mc Nemar's tests for two SVCs each)
alpha = 0.05; % significance threshold
tab_head = {strcat(ana_folds{7},' vs'),'df','t','p','lower_CI','upper_CI'};

for i=1:numel(res_dirs)
    %SVC1 = importdata(fullfile(res_dirs{i},ana_folds{7},'SVM.mat')); % base model for comparison
    % Note to myself: 'importdata' loads all contents of the matfile
    % directly into the variable, in this case, the fields of the SVM
    % struct, whereas the 'load' command would create a SVC1.SVM struct
    if i>2  % i.e. if the CSF feature set has a different sample size than the rest
        ana_folds{7} = '001_SVM_all_wo_CSF';
        SVC1 = importdata(fullfile(res_dirs{i},ana_folds{7},'SVM.mat')); % base model for comparison
        temp = cell(numel(ana_folds)-1,6); % table to be filled with stats
        temp(1,:) = tab_head;
        temp(2:end,1) = ana_folds(1:5);
        for j=1:5
            SVC2 = importdata(fullfile(res_dirs{i},ana_folds{j},'SVM.mat'));
            [h, p, ci, stats] = ML_SVC_ttest(SVC2, SVC1, alpha);
            temp(j+1,2:end) = {stats.df, stats.tstat, p, ci(1), ci(2)};
        end
    else    % i.e. if all feature sets, incl. CSF, have the same (small) sample size
        SVC1 = importdata(fullfile(res_dirs{i},ana_folds{7},'SVM.mat')); % base model for comparison
        temp = cell(numel(ana_folds),6); % table to be filled with stats
        temp(1,:) = tab_head;
        temp(2:end,1) = ana_folds(1:6);
        for j=1:6
            SVC2 = importdata(fullfile(res_dirs{i},ana_folds{j},'SVM.mat'));
            [h, p, ci, stats] = ML_SVC_ttest(SVC2, SVC1, alpha);
            temp(j+1,2:end) = {stats.df, stats.tstat, p, ci(1), ci(2)};
        end 
    end
    % perform Bonferroni-Holm correction and add to cell array
    [isSignificant,adjusted_pvals,alpha]= bonferroni_holm(cell2mat(temp(2:end,4)), alpha);
    temp(1,7)     = {'p_adjusted'};
    temp(2:end,7) = num2cell(adjusted_pvals);
    % save results cell array 'temp' as mat and excel
    save(fullfile(res_dirs{i},'infer_stats_ttests_all_vs_allinone.mat'),'temp');
    xlswrite(fullfile(res_dirs{i},'infer_stats_ttests_all_vs_allinone.xlsx'),temp);
end
        