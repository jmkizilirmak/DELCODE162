% Design for Konrad Waschkies's dissertation: 162_DELCODE_Personality SVC
% using resting-state fMRI mPerAF, psychometric scores from personality
% assessments, genetic and biomarkers associated with Alzheimer's disease.

% In this version, the exact same number of subjects are used for all
% variants. NaN exclusion accoring to SVM variant with the most features.
% However, CSF features are excluded from all variants except for the last,
% which only includes CSF. An extra NaN exclusion is used for the latter.
% This means that sample size is the same for variants 1-5, but different
% for 6. 4 classes (CN, SCD, MCI, AD) are included, but SCD and MCI are
% collapsed to form 1 class "AD risk". AD relatives are excluded.

% written by Jasmin Kizilirmak and Konrad Waschkies
% last edit: 2022-07-12 by JK -> corrected subject exclusion for base model

clear

%% specify directories
SVM_dir = '\\dzne-go-csv\Schott-Lab\commons\162_DELCODE_Personality\SVM_PerAF\SVMs_same_subjs_CSFxtra_riskclass_noADrel_noDMN\'; % output directory
imgs_dir = '\\dzne-go-csv\Schott-Lab\commons\162_DELCODE_Personality\PerAF_Results\mPerAF'; % where the imgs are
scripts_dir = '\\dzne-go-csv\Schott-Lab\commons\162_DELCODE_Personality\scripts';
subj_info_dir = '\\dzne-go-csv\Schott-Lab\commons\162_DELCODE_Personality';

%% specify SVM analysis parameters
options.SVM_type    = 'SVC';
options.C           = 1;        % soft margin
options.CV_mode     = 'kfc';    % kfc = k-fold cross-validation
options.k           = 10;       % number of CV folds
options.subs        = 30;       % number of subsamples
options.perm        = 1000;      % number of permutations

%% load subject file
subj_info = readtable(fullfile(subj_info_dir,'Antrag_162_Schott_Personality_20191115_SubjInfo_incl2.xlsx'));
num_subj = height(subj_info);

mask_img = fullfile('\\dzne-go-csv\Schott-Lab\commons\published_fROIs', 'DMN_mask_cut.nii');

%% get class names and non-imaging predictors (incl. confounding variables like sex, site)
x = str2num(str2mat(subj_info.prmdiag)); % CN = 0, SCD = 1, MCI = 2, AD = 5, AD relatives = 100
x(x==5) = 3; % AD
%x(x==2) = 2; % MCI, high risk group
x(x==1) = 2; % SCD, high risk group
x(x==0) = 1; % HC

% All categorical variables have to be dummycoded !
X = [subj_info.gender_male, subj_info.gender_female, subj_info.age,...
    subj_info.site10, subj_info.site11, subj_info.site13, subj_info.site14, ...
    subj_info.site16, subj_info.site17, subj_info.site18, subj_info.site2, ...
    subj_info.site5, subj_info.site8, ...
    subj_info.BFI_E, subj_info.BFI_A, subj_info.BFI_O, subj_info.BFI_N, subj_info.BFI_C,...
    subj_info.gds_total, subj_info.gaisf_total,...
    subj_info.ApoE_rec, subj_info.totaltau, ...
    subj_info.phosphotau181, subj_info.ratio_Abeta42_40];

%% prepare data images
data_imgs = cell(num_subj,1);
for i = 1:num_subj
    img_path = dir(strcat(imgs_dir,'/mPerAF_', subj_info.Repseudonym{i},'-M0_T1*.nii'));
    if numel(img_path) > 1, img_path = img_path(2); end
    data_imgs{i} = fullfile(img_path.folder,img_path.name);
end

%% exclude AD relatives
data_imgs = data_imgs(x~=100);
X         = X(x~=100,:);
x         = x(x~=100);

%% specify model variants
variants = {'_all_wo_CSF','_NI_only','_BFI_only', '_BFIext_only','_ApoE_only','_CSF_only','_counfounds_only'};

% We need the following variants of X:
Xselect{1} = deal(X(:,1:21));              % all included, except CSF
Xselect{2} = deal(X(:,1:13));              % confounding + neuroimaging only
Xselect{3} = deal(X(:,1:18));              % confounding, BFI; exclude neuroimaging !!!
Xselect{4} = deal(X(:,1:20));               % confounding, BFI, GDS, GAISF; exclude neuroimaging !!!
Xselect{5} = deal([X(:,1:13),X(:,21)]);     % confounding, APOE; exclude neuroimaging !!!
Xselect{6} = deal([X(:,1:13),X(:,22:24)]);  % confounding, CSF; exclude neuroimaging !!! 
Xselect{7} = deal(X(:,1:13));               % confounding only (base model)

%% exclude subjects with NaNs in Xselect from c, Xselect, and data_imgs_select

% loop through variants of X and create variants of c and data_imgs
for i = 1:length(Xselect)
    if (i~=6)
        if (i==1)
            empty = zeros(height(x),1);
            for j=1:height(x)
                if (sum(isnan(Xselect{i}(j,:)))>0)
                    empty(j) = 1;
                end
            end
            data_imgs    = data_imgs(empty==0);
        end
        Xselect{i} = Xselect{i}(empty==0,:);
        c{i}       = x(empty==0);
    else % CSF 'gold standard' variant
        empty2 = zeros(height(x),1);
        for j=1:height(x)
            if sum(isnan(Xselect{i}(j,:)))>0
                empty2(j) = 1;
            end
        end
        Xselect{i} = Xselect{i}(empty2==0,:);
        c{i}       = x(empty2==0);
    end
    if (i~=2)
        % standardize specific features for each sample size separately
        standardize_X = [3,14:width(Xselect{i})];
    else
        standardize_X = [3];
    end
    for k = 1:length(standardize_X)
        % mean center specific non-imaging features of Xselect
        Xselect{i}(:,standardize_X(k)) = Xselect{i}(:,standardize_X(k)) - repmat(mean(Xselect{i}(:,standardize_X(k)),'omitnan'),[height(Xselect{i}) 1]);
        % standardize SD of specified columns
        Xselect{i}(:,standardize_X(k)) = Xselect{i}(:,standardize_X(k))./ repmat(std(Xselect{i}(:,standardize_X(k)),'omitnan'),[height(Xselect{i}) 1]);
    end
end

%% perform SVM analysis

for i = 1:length(c)
    index = sprintf('%03d',i);
    curr_dir = [SVM_dir,index,'_SVM',variants{i}];
    if ~exist(curr_dir,'dir'), mkdir(curr_dir); end
    % if one of the two variants with neuroimaging...
    if (i < 3)
        % preprocessing steps for covariates
        preproc(1).op  = 'add_X';            % add non-imaging features X to imaging features Y
        preproc(1).cov = [1:size(Xselect{i},2)];
        preproc(2).op  = 'mc_Y';             % normalize neuroimaging data
        preproc(3).op  = 'std_Y';
        % perform SVM for NI
        ML4NI_SVM(curr_dir, mask_img, data_imgs, c{i}, [], Xselect{i}, preproc, options);
    % if no neuroimaging variant...
    else
        % Perform CV on points per class
        CV = ML_CV(c{i}, options.k, options.CV_mode);
        % Perform SVC  ML_SVC(x, Y, CV, C, perm, subs)
        SVM = ML_SVC(c{i}, Xselect{i}, CV,  options.C, options.perm, options.subs);
        % Save SVM results
        save(strcat(curr_dir,'/','SVM.mat'), 'SVM');
        % Save SVM design
        design.SVM_dir   = curr_dir;
        design.c         = c{i};
        design.X         = Xselect{i};
        design.preproc   = preproc;
        design.options   = options;
        save(strcat(curr_dir,'/','design.mat'), 'design');
    end
end