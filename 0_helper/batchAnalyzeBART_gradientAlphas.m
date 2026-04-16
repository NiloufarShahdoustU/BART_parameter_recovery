% % Script for looking at batch TDlearn

% modelRS_Metric is the learning rates for positive and negative alphas (ratio)
% taskRS_Metric is the difference between active and passive yellow balloons (-log)

%% IF TDLEARN RAN AS GRADIENT ALPHAS

% Set Up
close all;
sigBHFEEGcorr = struct();
unitIDs_Pts = {};
clusterTbl = [];
sigContacts_allPts = false(0,4);
moreImpulsive = false(0,1);
MNItrodeLocs_allPTs = [];
trodeLabels_allPts = {''};
anatomicalLocs_allPts = {''};
ptIDs_all = {''};
leftHemiContacts_allPts = [];
noWMpts = {};
[ptArray,bhvStruct,hazEEG] = BARTnumbers; 
nPts = length(ptArray);
set(0,'defaultfigurerenderer','painters')
pts = 1:nPts;
ptMap = cbrewer('qual','Set3',nPts);

% Loading in TDdataGradients for each patient
mainPath = '\\155.100.91.44\D\Data\preProcessed\BART_preprocessed\';
neuralData_All = {};

% Loading NeuralData All
for p = 1:length(ptArray)
    ptID = ptArray{p};
    dataPath = fullfile(mainPath, ptID, 'Data');
    fileInfo = dir(fullfile(dataPath, '*TDdataGradients.mat'));

    if ~isempty(fileInfo)
        fileToLoad = fullfile(dataPath, fileInfo(1).name);
        learningrate_data = load(fileToLoad);
        neuralData_All{p} = learningrate_data;
    else
        warning('No TDdataGradients.mat file found for %s', ptID);
        neuralData_All{p} = [];
    end
end

% Getting model learning rates and risk sensitivity 
for k = 1:nPts
    % Getting positive an negative learning rates
    rstd_positivePE(k) = neuralData_All{k}.TDdataGradients.positiveBetaMaxIdx/100; % positivePE dividing by 100 because it needs to go from 18 to 0.18
    rstd_negativePE(k) = neuralData_All{k}.TDdataGradients.negativeBetaMaxIdx/100; % negativePE
    % rstd metric: ratio of pos/neg
    modelRS_Metric(k) = (rstd_negativePE(k) - rstd_positivePE(k)) / (rstd_negativePE(k) + rstd_positivePE(k)); % calculating ratio.
end

% Getting task-derived risk sensitivity 
load('riskSensitivity_TD_bhv.mat'); % to get taskRS_Metric from behavior
impulsivityKLD_yellows = [bhvStruct.impulsivityKLD_yellows]; % apZVal_Yellow which is our impulsivity/risk aversion metric
taskRS_Metric = -impulsivityKLD_yellows;

%% 2: Getting data from all patients

% Getting data from models, RSTD R2 and pValues
for k = 1:nPts % patients
    for j = 1:length(neuralData_All{:,k}.TDdataGradients.neuralFit) % electrodes

        % ~~~~~~~~~~~~~~ P-VALUES ~~~~~~~~~~~~~~~~~~~~~~~ %
        % RSTD CUE MODEL ANOVA
        rstd_VE_Cue_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Cue.pValue(2);
        rstd_trialColor_Cue_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Cue.pValue(4);
        rstd_trialType_Cue_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Cue.pValue(5);
        % RSTD CUE SUCCESS MODEL ANOVA
        rstd_VE_CueSuccess_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_CueSuccess.pValue(2);
        rstd_trialColor_CueSuccess_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_CueSuccess.pValue(4);
        rstd_trialType_CueSuccess_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_CueSuccess.pValue(5);
        % RSTD OUTCOME MODEL ANOVA
        rstd_PE_Outcome_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Outcome.pValue(3);
        rstd_trialColor_Outcome_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Outcome.pValue(4);
        rstd_trialType_Outcome_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Outcome.pValue(5);
        rstd_outcome_Outcome_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Outcome.pValue(6);
        % RSTD OUTCOME SUCCESS MODEL ANOVA
        rstd_PE_OutcomeSuccess_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_OutcomeSuccess.pValue(3);
        rstd_trialColor_OutcomeSuccess_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_OutcomeSuccess.pValue(4);
        rstd_trialType_OutcomeSuccess_pValue(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_OutcomeSuccess.pValue(5);

        % ~~~~~~~~~~~~~~ R2 ADJUSTED ~~~~~~~~~~~~~~~~~~~~~~~ %
        % RSTD CUE MODEL
        rstd_Cue_R2adj(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueModel.Rsquared.Adjusted;
        rstd_CueSuccess_R2adj(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueSuccessModel.Rsquared.Adjusted;
        rstd_Outcome_R2adj(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeModel.Rsquared.Adjusted;
        rstd_OutcomeSuccess_R2adj(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeSuccessModel.Rsquared.Adjusted;

        % ~~~~~~~~~~~~~~ T-STATISTIC ~~~~~~~~~~~~~~~~~~~~~~~ %
        % RSTD CUE MODEL
        rstd_VE_Cue_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueModel.Coefficients(2,4);
        rstd_trialColorO_Cue_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueModel.Coefficients(3,4);
        rstd_trialColorR_Cue_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueModel.Coefficients(4,4);
        rstd_trialColorY_Cue_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueModel.Coefficients(5,4);
        rstd_trialType_Cue_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueModel.Coefficients(6,4);
        % RSTD CUE SUCCESS MODEL
        rstd_VE_CueSuccess_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueSuccessModel.Coefficients(2,4);
        rstd_trialColorO_CueSuccess_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueSuccessModel.Coefficients(3,4);
        rstd_trialColorR_CueSuccess_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueSuccessModel.Coefficients(4,4);
        rstd_trialColorY_CueSuccess_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueSuccessModel.Coefficients(5,4);
        rstd_trialType_CueSuccess_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_CueSuccessModel.Coefficients(6,4);
        % RSTD OUTCOME MODEL
        rstd_PE_Outcome_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeModel.Coefficients(3,4);
        rstd_trialColorO_Outcome_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeModel.Coefficients(4,4);
        rstd_trialColorR_Outcome_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeModel.Coefficients(5,4);
        rstd_trialColorY_Outcome_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeModel.Coefficients(6,4);
        rstd_trialType_Outcome_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeModel.Coefficients(7,4);
        rstd_outcome_Outcome_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeModel.Coefficients(8,4);
        % RSTD OUTCOME SUCCESS MODEL
        rstd_PE_OutcomeSuccess_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeSuccessModel.Coefficients(3,4);
        rstd_trialColorO_OutcomeSuccess_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeSuccessModel.Coefficients(4,4);
        rstd_trialColorR_OutcomeSuccess_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeSuccessModel.Coefficients(5,4);
        rstd_trialColorY_OutcomeSuccess_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeSuccessModel.Coefficients(6,4);
        rstd_trialType_OutcomeSuccess_tStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).rstd_OutcomeSuccessModel.Coefficients(7,4);

        % ~~~~~~~~~~~~~~ F-STATISTIC ~~~~~~~~~~~~~~~~~~~~~~~ %
        % RSTD CUE MODEL
        rstd_VE_Cue_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Cue(2,2);
        rstd_trialColor_Cue_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Cue(4,2);
        rstd_trialType_Cue_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Cue(5,2);
        % RSTD CUE SUCCESS MODEL
        rstd_VE_CueSuccess_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_CueSuccess(2,2);
        rstd_trialColor_CueSuccess_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_CueSuccess(4,2);
        rstd_trialType_CueSuccess_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_CueSuccess(5,2);
        % RSTD OUTCOME MODEL
        rstd_PE_Outcome_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Outcome(3,2);
        rstd_trialColor_Outcome_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Outcome(4,2);
        rstd_trialType_Outcome_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Outcome(5,2);
        rstd_outcome_Outcome_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_Outcome(6,2);
        % RSTD OUTCOME SUCCESS MODEL
        rstd_PE_OutcomeSuccess_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_OutcomeSuccess(3,2);
        rstd_trialColor_OutcomeSuccess_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_OutcomeSuccess(4,2);
        rstd_trialType_OutcomeSuccess_fStat(k,j) = neuralData_All{k}.TDdataGradients.neuralFit(j).ANOVA_rstd_OutcomeSuccess(5,2);
    end % neuralFit loop
end % pt loop

%% How many significant variables are there

alpha = 0.05; % Significance threshold

% Create logical matrices for significance
sig_rstd_VE_Cue              = (rstd_VE_Cue_pValue > 0)              & (rstd_VE_Cue_pValue < alpha);
sig_rstd_trialColor_Cue      = (rstd_trialColor_Cue_pValue > 0)      & (rstd_trialColor_Cue_pValue < alpha);
sig_rstd_trialType_Cue       = (rstd_trialType_Cue_pValue > 0)       & (rstd_trialType_Cue_pValue < alpha);
sig_rstd_VE_CueSuccess       = (rstd_VE_CueSuccess_pValue > 0)       & (rstd_VE_CueSuccess_pValue < alpha);
sig_rstd_trialColor_CueSuccess = (rstd_trialColor_CueSuccess_pValue > 0) & (rstd_trialColor_CueSuccess_pValue < alpha);
sig_rstd_trialType_CueSuccess  = (rstd_trialType_CueSuccess_pValue > 0)  & (rstd_trialType_CueSuccess_pValue < alpha);
sig_rstd_PE_Outcome          = (rstd_PE_Outcome_pValue > 0)          & (rstd_PE_Outcome_pValue < alpha);
sig_rstd_trialColor_Outcome  = (rstd_trialColor_Outcome_pValue > 0)  & (rstd_trialColor_Outcome_pValue < alpha);
sig_rstd_trialType_Outcome   = (rstd_trialType_Outcome_pValue > 0)   & (rstd_trialType_Outcome_pValue < alpha);
sig_rstd_outcome_Outcome     = (rstd_outcome_Outcome_pValue > 0)     & (rstd_outcome_Outcome_pValue < alpha);
sig_rstd_PE_OutcomeSuccess   = (rstd_PE_OutcomeSuccess_pValue > 0)   & (rstd_PE_OutcomeSuccess_pValue < alpha);
sig_rstd_trialColor_OutcomeSuccess = (rstd_trialColor_OutcomeSuccess_pValue > 0) & (rstd_trialColor_OutcomeSuccess_pValue < alpha);
sig_rstd_trialType_OutcomeSuccess  = (rstd_trialType_OutcomeSuccess_pValue > 0)  & (rstd_trialType_OutcomeSuccess_pValue < alpha);

%% Getting data of patient-wise significant contacts
% set up
pt_sig_rstd_VE_Cue               = sum(sig_rstd_VE_Cue,2);
pt_sig_rstd_trialColor_Cue       = sum(sig_rstd_trialColor_Cue,2);
pt_sig_rstd_trialType_Cue        = sum(sig_rstd_trialType_Cue,2);
pt_sig_rstd_PE_Outcome           = sum(sig_rstd_PE_Outcome,2);
pt_sig_rstd_trialColor_Outcome   = sum(sig_rstd_trialColor_Outcome,2);
pt_sig_rstd_trialType_Outcome    = sum(sig_rstd_trialType_Outcome,2);
pt_sig_rstd_outcome_Outcome      = sum(sig_rstd_outcome_Outcome,2);
pt_sig_rstd_VE_CueSuccess               = sum(sig_rstd_VE_CueSuccess,2);
pt_sig_rstd_trialColor_CueSuccess       = sum(sig_rstd_trialColor_CueSuccess,2);
pt_sig_rstd_trialType_CueSuccess        = sum(sig_rstd_trialType_CueSuccess,2);
pt_sig_rstd_PE_OutcomeSuccess           = sum(sig_rstd_PE_OutcomeSuccess,2);
pt_sig_rstd_trialColor_OutcomeSuccess   = sum(sig_rstd_trialColor_OutcomeSuccess,2);
pt_sig_rstd_trialType_OutcomeSuccess    = sum(sig_rstd_trialType_OutcomeSuccess,2);

dataLength = rstd_VE_Cue_pValue > 0; % everything that is not a zero represents an electrode
lengthElectrodes_PerPatient = sum(dataLength, 2);  % 71x1 vector of electrode lengths

% Proportion of significant p-values per patient
%(1) Cue aligned Model
sigVECue_ProportionPerPatient = pt_sig_rstd_VE_Cue ./ lengthElectrodes_PerPatient;
sigColorCue_ProportionPerPatient = pt_sig_rstd_trialColor_Cue ./ lengthElectrodes_PerPatient;
sigTypeCue_ProportionPerPatient = pt_sig_rstd_trialType_Cue ./ lengthElectrodes_PerPatient;
%(2) Outcome aligned Model
sigPEOutcome_ProportionPerPatient = pt_sig_rstd_PE_Outcome ./ lengthElectrodes_PerPatient;
sigColorOutcome_ProportionPerPatient = pt_sig_rstd_trialColor_Outcome ./ lengthElectrodes_PerPatient;
sigTypeOutcome_ProportionPerPatient = pt_sig_rstd_trialType_Outcome ./ lengthElectrodes_PerPatient;
sigOutcomeOutcome_ProportionPerPatient = pt_sig_rstd_outcome_Outcome ./ lengthElectrodes_PerPatient;

% Exclude patients with 100% VE cue proportion
include_idx = sigVECue_ProportionPerPatient < 1;

% Filtered data
filtered_propVE = sigVECue_ProportionPerPatient(include_idx);
filtered_propColor = sigColorCue_ProportionPerPatient(include_idx);
filtered_propType = sigTypeCue_ProportionPerPatient(include_idx);
filtered_rstd = modelRS_Metric(include_idx);
filtered_rstdPos = rstd_positivePE(include_idx);
filtered_rstdNeg = rstd_negativePE(include_idx);
filtered_taskRS_Metric = taskRS_Metric(include_idx);
filtered_propPE = sigPEOutcome_ProportionPerPatient(include_idx);
filtered_propColorOutcome = sigColorOutcome_ProportionPerPatient(include_idx);
filtered_propTypeOutcome = sigTypeOutcome_ProportionPerPatient(include_idx);
filtered_propOutcomeOutcome = sigOutcomeOutcome_ProportionPerPatient(include_idx);

% Summing over these variables:
n_sig_rstd_VE_Cue               = sum(pt_sig_rstd_VE_Cue(include_idx));
n_sig_rstd_trialColor_Cue       = sum(pt_sig_rstd_trialColor_Cue(include_idx));
n_sig_rstd_trialType_Cue        = sum(pt_sig_rstd_trialType_Cue(include_idx));
n_sig_rstd_VE_CueSuccess        = sum(pt_sig_rstd_VE_CueSuccess(include_idx));
n_sig_rstd_trialColor_CueSuccess = sum(pt_sig_rstd_trialColor_CueSuccess(include_idx));
n_sig_rstd_trialType_CueSuccess  = sum(pt_sig_rstd_trialType_CueSuccess(include_idx));
n_sig_rstd_PE_Outcome           = sum(pt_sig_rstd_PE_Outcome(include_idx));
n_sig_rstd_trialColor_Outcome   = sum(pt_sig_rstd_trialColor_Outcome(include_idx));
n_sig_rstd_trialType_Outcome    = sum(pt_sig_rstd_trialType_Outcome(include_idx));
n_sig_rstd_outcome_Outcome      = sum(pt_sig_rstd_outcome_Outcome(include_idx));
n_sig_rstd_PE_OutcomeSuccess    = sum(pt_sig_rstd_PE_OutcomeSuccess(include_idx));
n_sig_rstd_trialColor_OutcomeSuccess = sum(pt_sig_rstd_trialColor_OutcomeSuccess(include_idx));
n_sig_rstd_trialType_OutcomeSuccess  = sum(pt_sig_rstd_trialType_OutcomeSuccess(include_idx));

sigTotals = struct('VE_Cue', n_sig_rstd_VE_Cue, ...
    'TrialColor_Cue', n_sig_rstd_trialColor_Cue, ...
    'TrialType_Cue', n_sig_rstd_trialType_Cue, ...
    'VE_CueSuccess', n_sig_rstd_VE_CueSuccess, ...
    'TrialColor_CueSuccess', n_sig_rstd_trialColor_CueSuccess, ...
    'TrialType_CueSuccess', n_sig_rstd_trialType_CueSuccess, ...
    'PE_Outcome', n_sig_rstd_PE_Outcome, ...
    'TrialColor_Outcome', n_sig_rstd_trialColor_Outcome, ...
    'TrialType_Outcome', n_sig_rstd_trialType_Outcome, ...
    'Outcome_Outcome', n_sig_rstd_outcome_Outcome, ...
    'PE_OutcomeSuccess', n_sig_rstd_PE_OutcomeSuccess, ...
    'TrialColor_OutcomeSuccess', n_sig_rstd_trialColor_OutcomeSuccess, ...
    'TrialType_OutcomeSuccess', n_sig_rstd_trialType_OutcomeSuccess);

% ~~~~~~~~~~~ Supplementary Figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% plot significant numbers:
sig_counts = [n_sig_rstd_VE_Cue; n_sig_rstd_trialColor_Cue; n_sig_rstd_trialType_Cue;
    n_sig_rstd_VE_CueSuccess; n_sig_rstd_trialColor_CueSuccess; n_sig_rstd_trialType_CueSuccess;
    n_sig_rstd_PE_Outcome; n_sig_rstd_trialColor_Outcome; n_sig_rstd_trialType_Outcome; n_sig_rstd_outcome_Outcome;
    n_sig_rstd_PE_OutcomeSuccess; n_sig_rstd_trialColor_OutcomeSuccess; n_sig_rstd_trialType_OutcomeSuccess];
labels = {
    'VE Cue', 'Color Cue', 'Type Cue', ...
    'VE CueSucc', 'Color CueSucc', 'Type CueSucc', ...
    'PE Out', 'Color Out', 'Type Out', 'Outcome Out', ...
    'PE OutSucc', 'Color OutSucc', 'Type OutSucc'};
figure(1);
subplot(2,2,1)
bar(sig_counts(1:3));
set(gca, 'XTick', 1:3, 'XTickLabel', labels);
xtickangle(45);
ylabel('Number of Significant Units');
ylim([0 2000])
title('Cue Model');
grid on;
for i = 1:3
    text(i, sig_counts(i) + 1.5, num2str(sig_counts(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
axis square
subplot(2,2,2)
bar(sig_counts(4:6));
set(gca, 'XTick', 1:3, 'XTickLabel', labels(4:6));
xtickangle(45); % Rotate x-axis labels for readability
ylabel('Number of Significant Units');
ylim([0 2000])
title('Cue Success Model');
grid on;
for i = 1:3
    text(i, sig_counts(3+i) + 1.5, num2str(sig_counts(3+i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
axis square
subplot(2,2,3)
bar(sig_counts(7:10));
set(gca, 'XTick', 1:4, 'XTickLabel', labels(7:10));
xtickangle(45); % Rotate x-axis labels for readability
ylabel('Number of Significant Units');
ylim([0 2000])
title('Outcome Model');
grid on;
for i = 1:4
    text(i, sig_counts(6+i) + 1.5, num2str(sig_counts(6+i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
axis square
subplot(2,2,4)
bar(sig_counts(11:13));
set(gca, 'XTick', 1:3, 'XTickLabel', labels(11:13));
xtickangle(45);
ylabel('Number of Significant Units');
ylim([0 2000])
title('Outcome Success Model');
grid on;
for i = 1:3
    text(i, sig_counts(10+i) + 1.5, num2str(sig_counts(10+i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
axis square

saveas(1,'Z:\Data\Rhiannon\RSTD_BART\sFig_numberContacts_encoding_RSTDmodels.pdf') % save supplementary Fig
close(1)


%% Plot Alpha Values and RSTD metric

learningrate_data = {modelRS_Metric, rstd_positivePE, rstd_negativePE};
group = [];
values = [];
for i = 1:length(learningrate_data)
    values = [values; learningrate_data{i}(:)];
    group = [group; i * ones(size(learningrate_data{i}(:)))];
end

% plot (part of Manuscript Figure 1)
figure(2);
ax1 = axes;
hold(ax1, 'on');
box(ax1, 'on');
boxplot(values(group == 1), group(group == 1), ...
    'Colors', 'k', 'Symbol', '', 'Widths', 0.5, 'Positions', 1);
colors = lines(3);
jitterAmount = 0.1;
x1 = 1 + (rand(size(learningrate_data{1})) - 0.5) * jitterAmount;
scatter(ax1, x1, learningrate_data{1}, 30, 'filled', ...
    'MarkerFaceColor', colors(1, :), ...
    'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6);
set(ax1, 'XLim', [0.5 3.5], 'YColor', 'k');
xticks(ax1, 1:3);
xticklabels(ax1, {'RSTD Neural', 'Positive Alpha', 'Negative Alpha'});
ylabel(ax1, 'RSTD Value (Neural)');
title(ax1, 'RSTD Metrics with Split Y-Axis');
ax1_pos = get(ax1, 'Position');
axis square
ax2 = axes('Position', ax1_pos, ...
    'Color', 'none', 'YAxisLocation', 'right', ...
    'XAxisLocation', 'bottom', 'XColor', 'none', 'Box', 'off');
hold(ax2, 'on');
boxplot(ax2, values(ismember(group, [2 3])), group(ismember(group, [2 3])), ...
    'Colors', 'k', 'Symbol', '', 'Widths', 0.5, 'Positions', [2 3]);
for i = 2:3
    x = i + (rand(size(learningrate_data{i})) - 0.5) * jitterAmount;
    scatter(ax2, x, learningrate_data{i}, 40, 'filled', ...
        'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6);
end
set(ax2, 'XLim', [0.5 3.5], 'YColor', 'r');
ylabel(ax2, 'RSTD Value (Alpha)', 'Color', 'r');
linkaxes([ax1, ax2], 'x');
axis square

saveas(2,'Z:\Data\Rhiannon\RSTD_BART\mFig_alphas_ratio_RSTDmodels.pdf') % save supplementary Fig
close(2)

%% GOT TO HERE> REMOVE BAD PATIENTS: includeIdx

sig_rstd_VE_Cue = sig_rstd_VE_Cue(include_idx,:);
sig_rstd_trialColor_Cue = sig_rstd_trialColor_Cue(include_idx,:);
sig_rstd_trialType_Cue = sig_rstd_trialType_Cue(include_idx,:);
sig_rstd_VE_CueSuccess = sig_rstd_VE_CueSuccess(include_idx,:);
sig_rstd_trialColor_CueSuccess = sig_rstd_trialColor_CueSuccess(include_idx,:);
sig_rstd_trialType_CueSuccess = sig_rstd_trialType_CueSuccess(include_idx,:);

%% Set up data 
padded_sig_rstd_VE_Cue = NaN(size(sig_rstd_VE_Cue));  % size: 69x119
padded_sig_rstd_trialColor_Cue = NaN(size(sig_rstd_trialColor_Cue));  % size: 69x119
padded_sig_rstd_trialType_Cue = NaN(size(sig_rstd_trialType_Cue));  % size: 69x119
padded_sig_rstd_VE_CueSuccess = NaN(size(sig_rstd_VE_CueSuccess));  % size: 69x119
padded_sig_rstd_trialColor_CueSuccess = NaN(size(sig_rstd_trialColor_CueSuccess));  % size: 69x119
padded_sig_rstd_trialType_CueSuccess = NaN(size(sig_rstd_trialType_CueSuccess));  % size: 69x119
% Treat 0s as invalid/missing electrode entries

lengthElectrodes_PerPatient = lengthElectrodes_PerPatient(include_idx) % 69 patients 

% Loop over each row to add in NaNs
for i = 1:length(lengthElectrodes_PerPatient)
    nCols = lengthElectrodes_PerPatient(i);
    padded_sig_rstd_VE_Cue(i, 1:nCols) = sig_rstd_VE_Cue(i, 1:nCols);
    padded_sig_rstd_trialColor_Cue(i, 1:nCols) = sig_rstd_trialColor_Cue(i, 1:nCols);
    padded_sig_rstd_trialType_Cue(i, 1:nCols) = sig_rstd_trialType_Cue(i, 1:nCols);
    padded_sig_rstd_VE_CueSuccess(i, 1:nCols) = sig_rstd_VE_CueSuccess(i, 1:nCols);
    padded_sig_rstd_trialColor_CueSuccess(i, 1:nCols) = sig_rstd_trialColor_CueSuccess(i, 1:nCols);
    padded_sig_rstd_trialType_CueSuccess(i, 1:nCols) = sig_rstd_trialType_CueSuccess(i, 1:nCols);
end
% Flatten each matrix row-wise (patient-wise)
flat_sig_rstd_VE_Cue              = reshape(padded_sig_rstd_VE_Cue', [], 1);
flat_sig_rstd_trialColor_Cue      = reshape(padded_sig_rstd_trialColor_Cue', [], 1);
flat_sig_rstd_trialType_Cue       = reshape(padded_sig_rstd_trialType_Cue', [], 1);
flat_sig_rstd_VE_CueSuccess              = reshape(padded_sig_rstd_VE_CueSuccess', [], 1);
flat_sig_rstd_trialColor_CueSuccess      = reshape(padded_sig_rstd_trialColor_CueSuccess', [], 1);
flat_sig_rstd_trialType_CueSuccess       = reshape(padded_sig_rstd_trialType_CueSuccess', [], 1);
% Drop NaN rows
flat_sig_rstd_VE_Cue_clean = flat_sig_rstd_VE_Cue(~any(isnan(flat_sig_rstd_VE_Cue), 2), :);
flat_sig_rstd_trialColor_Cue_clean = flat_sig_rstd_trialColor_Cue(~any(isnan(flat_sig_rstd_trialColor_Cue), 2), :);
flat_sig_rstd_trialType_Cue_clean = flat_sig_rstd_trialType_Cue(~any(isnan(flat_sig_rstd_trialType_Cue), 2), :);
flat_sig_rstd_VE_CueSuccess_clean = flat_sig_rstd_VE_CueSuccess(~any(isnan(flat_sig_rstd_VE_CueSuccess), 2), :);
flat_sig_rstd_trialColor_CueSuccess_clean = flat_sig_rstd_trialColor_CueSuccess(~any(isnan(flat_sig_rstd_trialColor_CueSuccess), 2), :);
flat_sig_rstd_trialType_CueSuccess_clean = flat_sig_rstd_trialType_CueSuccess(~any(isnan(flat_sig_rstd_trialType_CueSuccess), 2), :);
% Change to logical
sig_rstd_VE_Cue_logical = logical(flat_sig_rstd_VE_Cue_clean);
sig_rstd_trialColor_Cue_logical = logical(flat_sig_rstd_trialColor_Cue_clean);
sig_rstd_trialType_Cue_logical = logical(flat_sig_rstd_trialType_Cue_clean);
sig_rstd_VE_CueSuccess_logical = logical(flat_sig_rstd_VE_CueSuccess_clean);
sig_rstd_trialColor_CueSuccess_logical = logical(flat_sig_rstd_trialColor_CueSuccess_clean);
sig_rstd_trialType_CueSuccess_logical = logical(flat_sig_rstd_trialType_CueSuccess_clean);

%% Plotting R2 for all models. Norm and Success

rstd_Cue_R2adj = rstd_Cue_R2adj(include_idx,:);
rstd_Outcome_R2adj = rstd_Outcome_R2adj(include_idx,:);
rstd_CueSuccess_R2adj = rstd_CueSuccess_R2adj(include_idx,:);
rstd_OutcomeSuccess_R2adj = rstd_OutcomeSuccess_R2adj(include_idx,:);

% Initialize output matrix with NaNs
padded_rstd_Cue_R2adj = NaN(size(rstd_Cue_R2adj));  % size: 69x119
padded_rstd_Outcome_R2adj = NaN(size(rstd_Outcome_R2adj));  % size: 69x119
padded_rstd_CueSuccess_R2adj = NaN(size(rstd_CueSuccess_R2adj));  % size: 69x119
padded_rstd_OutcomeSuccess_R2adj = NaN(size(rstd_OutcomeSuccess_R2adj));  % size: 69x119
% Loop over each row to add in NaNs
for i = 1:length(lengthElectrodes_PerPatient)
    nCols = lengthElectrodes_PerPatient(i);
    padded_rstd_Cue_R2adj(i, 1:nCols) = rstd_Cue_R2adj(i, 1:nCols);
    padded_rstd_Outcome_R2adj(i, 1:nCols) = rstd_Outcome_R2adj(i, 1:nCols);
    padded_rstd_CueSuccess_R2adj(i, 1:nCols) = rstd_CueSuccess_R2adj(i, 1:nCols);
    padded_rstd_OutcomeSuccess_R2adj(i, 1:nCols) = rstd_OutcomeSuccess_R2adj(i, 1:nCols);
end
% Flatten each matrix row-wise (patient-wise)
flat_rstd_Cue_R2adj              = reshape(padded_rstd_Cue_R2adj', [], 1);
flat_rstd_Outcome_R2adj      = reshape(padded_rstd_Outcome_R2adj', [], 1);
flat_rstd_CueSuccess_R2adj       = reshape(padded_rstd_CueSuccess_R2adj', [], 1);
flat_rstd_OutcomeSuccess_R2adj       = reshape(padded_rstd_OutcomeSuccess_R2adj', [], 1);
% Drop NaN rows
flat_rstd_Cue_R2adj_clean = flat_rstd_Cue_R2adj(~any(isnan(flat_rstd_Cue_R2adj), 2), :);
flat_rstd_CueSuccess_R2adj_clean = flat_rstd_CueSuccess_R2adj(~any(isnan(flat_rstd_CueSuccess_R2adj), 2), :);
flat_rstd_Outcome_R2adj_clean = flat_rstd_Outcome_R2adj(~any(isnan(flat_rstd_Outcome_R2adj), 2), :);
flat_rstd_OutcomeSuccess_R2adj_clean = flat_rstd_OutcomeSuccess_R2adj(~any(isnan(flat_rstd_OutcomeSuccess_R2adj), 2), :);

% Plot R2s for CUE models
figure(3)
subplot(2,3,1);
histogram(flat_rstd_Cue_R2adj_clean(~sig_rstd_VE_Cue_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_Cue_R2adj_clean(sig_rstd_VE_Cue_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('VE R^2');
xlabel('Cue Model R^2 Adj');
ylabel('Frequency');
axis square
subplot(2,3,2);
histogram(flat_rstd_Cue_R2adj_clean(~sig_rstd_trialColor_Cue_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_Cue_R2adj_clean(sig_rstd_trialColor_Cue_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('Color R^2');
xlabel('Cue Model R^2 Adj');
ylabel('Frequency');
axis square
subplot(2,3,3);
histogram(flat_rstd_Cue_R2adj_clean(~sig_rstd_trialType_Cue_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_Cue_R2adj_clean(sig_rstd_trialType_Cue_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('Type R^2');
xlabel('Cue Model R^2 Adj');
ylabel('Frequency');
axis square
% CUE SUCCESS MODEL
subplot(2,3,4);
axis square
histogram(flat_rstd_CueSuccess_R2adj_clean(~sig_rstd_VE_CueSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_CueSuccess_R2adj_clean(sig_rstd_VE_CueSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('VE R^2');
xlabel('Cue Sucess Model R^2 Adj');
ylabel('Frequency');
axis square
subplot(2,3,5);
histogram(flat_rstd_CueSuccess_R2adj_clean(~sig_rstd_trialColor_CueSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_CueSuccess_R2adj_clean(sig_rstd_trialColor_CueSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('Color R^2');
xlabel('Cue Sucess Model R^2 Adj');
ylabel('Frequency');
axis square
subplot(2,3,6);
histogram(flat_rstd_CueSuccess_R2adj_clean(~sig_rstd_trialType_CueSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_CueSuccess_R2adj_clean(sig_rstd_trialType_CueSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('Type R^2');
xlabel('Cue Success Model R^2 Adj');
ylabel('Frequency');
axis square

saveas(3,'Z:\Data\Rhiannon\RSTD_BART\mFig_R2s_CueModels.pdf') % save supplementary Fig
close(3)

%% OUTCOME R2s

sig_rstd_PE_Outcome = sig_rstd_PE_Outcome(include_idx,:);
sig_rstd_trialColor_Outcome = sig_rstd_trialColor_Outcome(include_idx,:);
sig_rstd_trialType_Outcome = sig_rstd_trialType_Outcome(include_idx,:);
sig_rstd_outcome_Outcome = sig_rstd_outcome_Outcome(include_idx,:);
sig_rstd_PE_OutcomeSuccess = sig_rstd_PE_OutcomeSuccess(include_idx,:);
sig_rstd_trialColor_OutcomeSuccess = sig_rstd_trialColor_OutcomeSuccess(include_idx,:);
sig_rstd_trialType_OutcomeSuccess = sig_rstd_trialType_OutcomeSuccess(include_idx,:);

% Initialize output matrix with NaNs
padded_sig_rstd_PE_Outcome = NaN(size(sig_rstd_PE_Outcome));  % size: 69x119
padded_sig_rstd_trialColor_Outcome = NaN(size(sig_rstd_trialColor_Outcome));  % size: 69x119
padded_sig_rstd_trialType_Outcome = NaN(size(sig_rstd_trialType_Outcome));  % size: 69x119
padded_sig_rstd_Outcome_Outcome = NaN(size(sig_rstd_outcome_Outcome));  % size: 69x119
padded_sig_rstd_PE_OutcomeSuccess = NaN(size(sig_rstd_PE_OutcomeSuccess));  % size: 69x119
padded_sig_rstd_trialColor_OutcomeSuccess = NaN(size(sig_rstd_trialColor_OutcomeSuccess));  % size: 69x119
padded_sig_rstd_trialType_OutcomeSuccess = NaN(size(sig_rstd_trialType_OutcomeSuccess));  % size: 69x119
% Loop over each row to add in NaNs
for i = 1:length(lengthElectrodes_PerPatient)
    nCols = lengthElectrodes_PerPatient(i);
    padded_sig_rstd_PE_Outcome(i, 1:nCols) = sig_rstd_PE_Outcome(i, 1:nCols);
    padded_sig_rstd_trialColor_Outcome(i, 1:nCols) = sig_rstd_trialColor_Outcome(i, 1:nCols);
    padded_sig_rstd_trialType_Outcome(i, 1:nCols) = sig_rstd_trialType_Outcome(i, 1:nCols);
    padded_sig_rstd_Outcome_Outcome(i, 1:nCols) = sig_rstd_outcome_Outcome(i, 1:nCols);
    padded_sig_rstd_PE_OutcomeSuccess(i, 1:nCols) = sig_rstd_PE_OutcomeSuccess(i, 1:nCols);
    padded_sig_rstd_trialColor_OutcomeSuccess(i, 1:nCols) = sig_rstd_trialColor_OutcomeSuccess(i, 1:nCols);
    padded_sig_rstd_trialType_OutcomeSuccess(i, 1:nCols) = sig_rstd_trialType_OutcomeSuccess(i, 1:nCols);
end
% Flatten each matrix row-wise (patient-wise)
flat_sig_rstd_PE_Outcome              = reshape(padded_sig_rstd_PE_Outcome', [], 1);
flat_sig_rstd_trialColor_Outcome      = reshape(padded_sig_rstd_trialColor_Outcome', [], 1);
flat_sig_rstd_trialType_Outcome       = reshape(padded_sig_rstd_trialType_Outcome', [], 1);
flat_sig_rstd_Outcome_Outcome       = reshape(padded_sig_rstd_Outcome_Outcome', [], 1);
flat_sig_rstd_PE_OutcomeSuccess              = reshape(padded_sig_rstd_PE_OutcomeSuccess', [], 1);
flat_sig_rstd_trialColor_OutcomeSuccess      = reshape(padded_sig_rstd_trialColor_OutcomeSuccess', [], 1);
flat_sig_rstd_trialType_OutcomeSuccess       = reshape(padded_sig_rstd_trialType_OutcomeSuccess', [], 1);
% Drop NaN rows
flat_sig_rstd_PE_Outcome_clean = flat_sig_rstd_PE_Outcome(~any(isnan(flat_sig_rstd_PE_Outcome), 2), :);
flat_sig_rstd_trialColor_Outcome_clean = flat_sig_rstd_trialColor_Outcome(~any(isnan(flat_sig_rstd_trialColor_Outcome), 2), :);
flat_sig_rstd_trialType_Outcome_clean = flat_sig_rstd_trialType_Outcome(~any(isnan(flat_sig_rstd_trialType_Outcome), 2), :);
flat_sig_rstd_Outcome_Outcome_clean = flat_sig_rstd_Outcome_Outcome(~any(isnan(flat_sig_rstd_Outcome_Outcome), 2), :);
flat_sig_rstd_PE_OutcomeSuccess_clean = flat_sig_rstd_PE_OutcomeSuccess(~any(isnan(flat_sig_rstd_PE_OutcomeSuccess), 2), :);
flat_sig_rstd_trialColor_OutcomeSuccess_clean = flat_sig_rstd_trialColor_OutcomeSuccess(~any(isnan(flat_sig_rstd_trialColor_OutcomeSuccess), 2), :);
flat_sig_rstd_trialType_OutcomeSuccess_clean = flat_sig_rstd_trialType_OutcomeSuccess(~any(isnan(flat_sig_rstd_trialType_OutcomeSuccess), 2), :);
% Change to logical
sig_rstd_PE_Outcome_logical = logical(flat_sig_rstd_PE_Outcome_clean);
sig_rstd_trialColor_Outcome_logical = logical(flat_sig_rstd_trialColor_Outcome_clean);
sig_rstd_trialType_Outcome_logical = logical(flat_sig_rstd_trialType_Outcome_clean);
sig_rstd_Outcome_Outcome_logical = logical(flat_sig_rstd_Outcome_Outcome_clean);
sig_rstd_PE_OutcomeSuccess_logical = logical(flat_sig_rstd_PE_OutcomeSuccess_clean);
sig_rstd_trialColor_OutcomeSuccess_logical = logical(flat_sig_rstd_trialColor_OutcomeSuccess_clean);
sig_rstd_trialType_OutcomeSuccess_logical = logical(flat_sig_rstd_trialType_OutcomeSuccess_clean);

% Plot R2s for OUTCOME models
figure(4)
subplot(2,4,1);
histogram(flat_rstd_Outcome_R2adj_clean(~sig_rstd_PE_Outcome_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_Outcome_R2adj_clean(sig_rstd_PE_Outcome_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('PE R^2');
xlabel('Outcome Model R^2 Adj');
ylabel('Frequency');
axis square
subplot(2,4,2);
histogram(flat_rstd_Outcome_R2adj_clean(~sig_rstd_trialColor_Outcome_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_Outcome_R2adj_clean(sig_rstd_trialColor_Outcome_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('Color R^2');
xlabel('Outcome Model R^2 Adj');
ylabel('Frequency');
axis square
subplot(2,4,3);
histogram(flat_rstd_Outcome_R2adj_clean(~sig_rstd_trialType_Outcome_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_Outcome_R2adj_clean(sig_rstd_trialType_Outcome_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('Type R^2');
xlabel('Outcome Model R^2 Adj');
ylabel('Frequency');
axis square
subplot(2,4,4);
histogram(flat_rstd_Outcome_R2adj_clean(~sig_rstd_Outcome_Outcome_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_Outcome_R2adj_clean(sig_rstd_Outcome_Outcome_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('Outcome R^2');
xlabel('Outcome Model R^2 Adj');
ylabel('Frequency');
axis square
% OUTCOME SUCCESS
subplot(2,4,5);
histogram(flat_rstd_OutcomeSuccess_R2adj_clean(~sig_rstd_PE_OutcomeSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_OutcomeSuccess_R2adj_clean(sig_rstd_PE_OutcomeSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('PE R^2');
xlabel('Outcome Success Model R^2 Adj');
ylabel('Frequency');
axis square
subplot(2,4,6);
histogram(flat_rstd_OutcomeSuccess_R2adj_clean(~sig_rstd_trialColor_OutcomeSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_OutcomeSuccess_R2adj_clean(sig_rstd_trialColor_OutcomeSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('Color R^2');
xlabel('Outcome Success Model R^2 Adj');
ylabel('Frequency');
axis square
subplot(2,4,7);
histogram(flat_rstd_OutcomeSuccess_R2adj_clean(~sig_rstd_trialType_OutcomeSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'k');
hold on
histogram(flat_rstd_OutcomeSuccess_R2adj_clean(sig_rstd_trialType_OutcomeSuccess_logical), 'BinMethod', 'auto', 'FaceColor', 'b');
title('Type R^2');
xlabel('Outcome Success Model R^2 Adj');
ylabel('Frequency');
axis square

saveas(4,'Z:\Data\Rhiannon\RSTD_BART\mFig_R2s_OutcomeModels.pdf') % save mFig
close(4)

%% ~~~~~~~~~` mean of top 100 R2's ~~~~~~~~~~~~ %%
% take top R2 values from all contacts (should be significant?)
top100_rstd_Cue_R2adj = maxk(flat_rstd_Cue_R2adj_clean(:), 100);
top100_rstd_Outcome_R2adj = maxk(flat_rstd_Outcome_R2adj_clean(:), 100);

mean(top100_rstd_Cue_R2adj)
mean(top100_rstd_Outcome_R2adj)

% subset by Value and PE?
top100_rstd_CueSuccess_VE_R2adj = maxk(flat_rstd_CueSuccess_R2adj_clean(sig_rstd_VE_CueSuccess_logical), 100);
mean(top100_rstd_CueSuccess_VE_R2adj) % as an example - most significant contacts have the highest R2s.
[~, idx] = max(flat_rstd_Cue_R2adj_clean);
maxR2_index = false(size(flat_rstd_Cue_R2adj_clean));
maxR2_index(idx) = true;
[sorted_vals, sorted_idx] = sort(flat_rstd_Cue_R2adj_clean, 'descend');

%% Plot t-statistics
rstd_VE_Cue_tStat = rstd_VE_Cue_tStat(include_idx,:);
rstd_trialColorO_Cue_tStat = rstd_trialColorO_Cue_tStat(include_idx,:);
rstd_trialColorR_Cue_tStat = rstd_trialColorR_Cue_tStat(include_idx,:);
rstd_trialColorY_Cue_tStat = rstd_trialColorY_Cue_tStat(include_idx,:);
rstd_trialType_Cue_tStat = rstd_trialType_Cue_tStat(include_idx,:);
rstd_PE_Outcome_tStat = rstd_PE_Outcome_tStat(include_idx,:);
rstd_trialColorO_Outcome_tStat = rstd_trialColorO_Outcome_tStat(include_idx,:);
rstd_trialColorR_Outcome_tStat = rstd_trialColorR_Outcome_tStat(include_idx,:);
rstd_trialColorY_Outcome_tStat = rstd_trialColorY_Outcome_tStat(include_idx,:);
rstd_trialType_Outcome_tStat = rstd_trialType_Outcome_tStat(include_idx,:);
rstd_outcome_Outcome_tStat = rstd_outcome_Outcome_tStat(include_idx,:);

% Initialize output matrix with NaNs
padded_rstd_VE_Cue_tStat = NaN(size(rstd_VE_Cue_tStat));  % size: 69x119
padded_rstd_trialColorO_Cue_tStat = NaN(size(rstd_trialColorO_Cue_tStat));  % size: 69x119
padded_rstd_trialColorR_Cue_tStat = NaN(size(rstd_trialColorR_Cue_tStat));  % size: 69x119
padded_rstd_trialColorY_Cue_tStat = NaN(size(rstd_trialColorY_Cue_tStat));  % size:69x119
padded_rstd_trialType_Cue_tStat = NaN(size(rstd_trialType_Cue_tStat));  % size: 69x119
padded_rstd_PE_Outcome_tStat = NaN(size(rstd_PE_Outcome_tStat));  % size:69x119
padded_rstd_trialColorO_Outcome_tStat = NaN(size(rstd_trialColorO_Outcome_tStat));  % size: 69x119
padded_rstd_trialColorR_Outcome_tStat = NaN(size(rstd_trialColorR_Outcome_tStat));  % size: 69x119
padded_rstd_trialColorY_Outcome_tStat = NaN(size(rstd_trialColorY_Outcome_tStat));  % size:69x119
padded_rstd_trialType_Outcome_tStat = NaN(size(rstd_trialType_Outcome_tStat));  % size:69x119
padded_rstd_outcome_Outcome_tStat = NaN(size(rstd_outcome_Outcome_tStat));  % size: 69x119
% Loop over each row to add in NaNs
for i = 1:length(lengthElectrodes_PerPatient)
    nCols = lengthElectrodes_PerPatient(i);
    padded_rstd_VE_Cue_tStat(i, 1:nCols) = rstd_VE_Cue_tStat(i, 1:nCols);
    padded_rstd_trialColorO_Cue_tStat(i, 1:nCols) = rstd_trialColorO_Cue_tStat(i, 1:nCols);
    padded_rstd_trialColorR_Cue_tStat(i, 1:nCols) = rstd_trialColorR_Cue_tStat(i, 1:nCols);
    padded_rstd_trialColorY_Cue_tStat(i, 1:nCols) = rstd_trialColorY_Cue_tStat(i, 1:nCols);
    padded_rstd_trialType_Cue_tStat(i, 1:nCols) = rstd_trialType_Cue_tStat(i, 1:nCols);
    padded_rstd_PE_Outcome_tStat(i, 1:nCols) = rstd_PE_Outcome_tStat(i, 1:nCols);
    padded_rstd_trialColorO_Outcome_tStat(i, 1:nCols) = rstd_trialColorO_Outcome_tStat(i, 1:nCols);
    padded_rstd_trialColorR_Outcome_tStat(i, 1:nCols) = rstd_trialColorR_Outcome_tStat(i, 1:nCols);
    padded_rstd_trialColorY_Outcome_tStat(i, 1:nCols) = rstd_trialColorY_Outcome_tStat(i, 1:nCols);
    padded_rstd_trialType_Outcome_tStat(i, 1:nCols) = rstd_trialType_Outcome_tStat(i, 1:nCols);
    padded_rstd_outcome_Outcome_tStat(i, 1:nCols) = rstd_outcome_Outcome_tStat(i, 1:nCols);
end
% Flatten each matrix row-wise (patient-wise)
flat_rstd_VE_Cue_tStat              = reshape(padded_rstd_VE_Cue_tStat', [], 1);
flat_rstd_trialColorO_Cue_tStat      = reshape(padded_rstd_trialColorO_Cue_tStat', [], 1);
flat_rstd_trialColorR_Cue_tStat       = reshape(padded_rstd_trialColorR_Cue_tStat', [], 1);
flat_rstd_trialColorY_Cue_tStat       = reshape(padded_rstd_trialColorY_Cue_tStat', [], 1);
flat_rstd_trialType_Cue_tStat              = reshape(padded_rstd_trialType_Cue_tStat', [], 1);
flat_rstd_PE_Outcome_tStat      = reshape(padded_rstd_PE_Outcome_tStat', [], 1);
flat_rstd_trialColorO_Outcome_tStat       = reshape(padded_rstd_trialColorO_Outcome_tStat', [], 1);
flat_rstd_trialColorR_Outcome_tStat              = reshape(padded_rstd_trialColorR_Outcome_tStat', [], 1);
flat_rstd_trialColorY_Outcome_tStat      = reshape(padded_rstd_trialColorY_Outcome_tStat', [], 1);
flat_rstd_trialType_Outcome_tStat       = reshape(padded_rstd_trialType_Outcome_tStat', [], 1);
flat_rstd_outcome_Outcome_tStat      = reshape(padded_rstd_outcome_Outcome_tStat', [], 1);
% Drop NaN rows
flat_rstd_VE_Cue_tStat_clean = flat_rstd_VE_Cue_tStat(~any(isnan(flat_rstd_VE_Cue_tStat), 2), :);
flat_rstd_trialColorO_Cue_tStat_clean = flat_rstd_trialColorO_Cue_tStat(~any(isnan(flat_rstd_trialColorO_Cue_tStat), 2), :);
flat_rstd_trialColorR_Cue_tStat_clean = flat_rstd_trialColorR_Cue_tStat(~any(isnan(flat_rstd_trialColorR_Cue_tStat), 2), :);
flat_rstd_trialColorY_Cue_tStat_clean = flat_rstd_trialColorY_Cue_tStat(~any(isnan(flat_rstd_trialColorY_Cue_tStat), 2), :);
flat_rstd_trialType_Cue_tStat_clean = flat_rstd_trialType_Cue_tStat(~any(isnan(flat_rstd_trialType_Cue_tStat), 2), :);
flat_rstd_PE_Outcome_tStat_clean = flat_rstd_PE_Outcome_tStat(~any(isnan(flat_rstd_PE_Outcome_tStat), 2), :);
flat_rstd_trialColorO_Outcome_tStat_clean = flat_rstd_trialColorO_Outcome_tStat(~any(isnan(flat_rstd_trialColorO_Outcome_tStat), 2), :);
flat_rstd_trialColorR_Outcome_tStat_clean = flat_rstd_trialColorR_Outcome_tStat(~any(isnan(flat_rstd_trialColorR_Outcome_tStat), 2), :);
flat_rstd_trialColorY_Outcome_tStat_clean = flat_rstd_trialColorY_Outcome_tStat(~any(isnan(flat_rstd_trialColorY_Outcome_tStat), 2), :);
flat_rstd_trialType_Outcome_tStat_clean = flat_rstd_trialType_Outcome_tStat(~any(isnan(flat_rstd_trialType_Outcome_tStat), 2), :);
flat_rstd_outcome_Outcome_tStat_clean = flat_rstd_outcome_Outcome_tStat(~any(isnan(flat_rstd_outcome_Outcome_tStat), 2), :);

%% plotting tstats against each other
% Combinations
sigC_both = sig_rstd_VE_Cue_logical & sig_rstd_trialType_Cue_logical;
sigC_ve_only = sig_rstd_VE_Cue_logical & ~sig_rstd_trialType_Cue_logical;
sigC_trialType_only = sig_rstd_trialType_Cue_logical & ~sig_rstd_VE_Cue_logical;

model_colorMap = [
    0.0000, 0.2471, 0.3608;  % Dark Blue: VE
    0.4761, 0.4510, 0.6423;  % Purple: PE/VE + Type
    0.7373, 0.3137, 0.5647;  % Pink: COLOR
    1.0000, 0.3882, 0.3804;  % Coral: PE/VE + Color
    1.0000, 0.6510, 0.0000;  % Orange Yellow: TYPE
    0.0000, 0.6235, 0.4784;   % Teal Green: PE
    0.678, 0.847, 0.902;   % Light Blue: Outcome
    0.5, 0, 0   % Maroon: PE & Outcome
];

% CUE: VE vs TRIAL TYPE
figure(5); 
subplot(2,2,1)
hold on;
scatter(flat_rstd_VE_Cue_tStat_clean, flat_rstd_trialType_Cue_tStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_VE_Cue_tStat_clean(sigC_ve_only), flat_rstd_trialType_Cue_tStat_clean(sigC_ve_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(1,:), 'MarkerEdgeColor','k');       % VE only - blue 
scatter(flat_rstd_VE_Cue_tStat_clean(sigC_trialType_only), flat_rstd_trialType_Cue_tStat_clean(sigC_trialType_only), 30, 'b', 'filled','o',  'MarkerFaceColor',  model_colorMap(5,:), 'MarkerEdgeColor','k');  % TrialType only - orange
scatter(flat_rstd_VE_Cue_tStat_clean(sigC_both), flat_rstd_trialType_Cue_tStat_clean(sigC_both), 30, 'y', 'filled','o', 'MarkerFaceColor',  model_colorMap(2,:), 'MarkerEdgeColor','k');             % Both - purple
% Labels and formatting
xlabel('VE Cue t-Statistic');
ylabel('TrialType Cue t-Statistic');
title('Scatter: VE vs TrialType Cue t-Statistics');
xlim([-20 10]);
ylim([-5 5]);
grid on;
legend({'NS','Sig VE','Sig TrialType','Sig Both'}, 'Location', 'best');
hold off;
axis square

% CUE: VE vs COLOR Orange
subplot(2,2,2)
% Combinations
sigC_both = sig_rstd_VE_Cue_logical & sig_rstd_trialColor_Cue_logical;
sigC_ve_only = sig_rstd_VE_Cue_logical & ~sig_rstd_trialColor_Cue_logical;
sig_trialColor_only = sig_rstd_trialColor_Cue_logical & ~sig_rstd_VE_Cue_logical;
hold on;
scatter(flat_rstd_VE_Cue_tStat_clean, flat_rstd_trialColorO_Cue_tStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_VE_Cue_tStat_clean(sigC_ve_only), flat_rstd_trialColorO_Cue_tStat_clean(sigC_ve_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(1,:), 'MarkerEdgeColor','k');       % VE only - blue 
scatter(flat_rstd_VE_Cue_tStat_clean(sig_trialColor_only), flat_rstd_trialColorO_Cue_tStat_clean(sig_trialColor_only), 30, 'b', 'filled','o',  'MarkerFaceColor', model_colorMap(3,:), 'MarkerEdgeColor','k');  % TrialColor only - pink
scatter(flat_rstd_VE_Cue_tStat_clean(sigC_both), flat_rstd_trialColorO_Cue_tStat_clean(sigC_both), 30, 'y', 'filled','o', 'MarkerFaceColor', model_colorMap(4,:), 'MarkerEdgeColor','k');             % Both - cyan
% Labels and formatting
xlabel('VE Cue t-Statistic');
ylabel('TrialType Cue t-Statistic');
title('Scatter: VE vs TrialColorO Cue t-Statistics');
xlim([-20 10]);
ylim([-5 5]);
grid on;
legend({'NS','Sig VE','Sig TrialColor','Sig Both'}, 'Location', 'best');
hold off;
axis square

% CUE: VE vs COLOR Red
subplot(2,2,3)
% Combinations
sigC_both = sig_rstd_VE_Cue_logical & sig_rstd_trialColor_Cue_logical;
sigC_ve_only = sig_rstd_VE_Cue_logical & ~sig_rstd_trialColor_Cue_logical;
sig_trialColor_only = sig_rstd_trialColor_Cue_logical & ~sig_rstd_VE_Cue_logical;
hold on;
scatter(flat_rstd_VE_Cue_tStat_clean, flat_rstd_trialColorO_Cue_tStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_VE_Cue_tStat_clean(sigC_ve_only), flat_rstd_trialColorR_Cue_tStat_clean(sigC_ve_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor',  model_colorMap(1,:), 'MarkerEdgeColor','k');       % VE only - green 
scatter(flat_rstd_VE_Cue_tStat_clean(sig_trialColor_only), flat_rstd_trialColorR_Cue_tStat_clean(sig_trialColor_only), 30, 'b', 'filled','o',  'MarkerFaceColor', model_colorMap(3,:), 'MarkerEdgeColor','k');  % TrialColor only - 
scatter(flat_rstd_VE_Cue_tStat_clean(sigC_both), flat_rstd_trialColorR_Cue_tStat_clean(sigC_both), 30, 'y', 'filled','o', 'MarkerFaceColor', model_colorMap(4,:), 'MarkerEdgeColor','k');             % Both - yellow
% Labels and formatting
xlabel('VE Cue t-Statistic');
ylabel('TrialType Cue t-Statistic');
title('Scatter: VE vs TrialColorR Cue t-Statistics');
xlim([-20 10]);
ylim([-5 5]);
grid on;
legend({'NS','Sig VE','Sig TrialColor','Sig Both'}, 'Location', 'best');
hold off;
axis square

% CUE: VE vs COLOR Yellow
subplot(2,2,4)
% Combinations
sigC_both = sig_rstd_VE_Cue_logical & sig_rstd_trialColor_Cue_logical;
sigC_ve_only = sig_rstd_VE_Cue_logical & ~sig_rstd_trialColor_Cue_logical;
sig_trialColor_only = sig_rstd_trialColor_Cue_logical & ~sig_rstd_VE_Cue_logical;
hold on;
scatter(flat_rstd_VE_Cue_tStat_clean, flat_rstd_trialColorO_Cue_tStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_VE_Cue_tStat_clean(sigC_ve_only), flat_rstd_trialColorY_Cue_tStat_clean(sigC_ve_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(1,:), 'MarkerEdgeColor','k');       % VE only - green 
scatter(flat_rstd_VE_Cue_tStat_clean(sig_trialColor_only), flat_rstd_trialColorY_Cue_tStat_clean(sig_trialColor_only), 30, 'b', 'filled','o',  'MarkerFaceColor', model_colorMap(3,:), 'MarkerEdgeColor','k');  % TrialColor only - 
scatter(flat_rstd_VE_Cue_tStat_clean(sigC_both), flat_rstd_trialColorY_Cue_tStat_clean(sigC_both), 30, 'y', 'filled','o', 'MarkerFaceColor', model_colorMap(4,:), 'MarkerEdgeColor','k');             % Both - yellow
% Labels and formatting
xlabel('VE Cue t-Statistic');
ylabel('TrialType Cue t-Statistic');
title('Scatter: VE vs TrialColorY Cue t-Statistics');
xlim([-20 10]);
ylim([-5 5]);
grid on;
legend({'NS','Sig VE','Sig TrialColor','Sig Both'}, 'Location', 'best');
hold off;
axis square

saveas(5,'Z:\Data\Rhiannon\RSTD_BART\mFig_tstats_CueModels.pdf') % save mFig
close(5)

%% PLOT TSTATS FOR OUTCOME.

% Combinations
sigO_both = sig_rstd_PE_Outcome_logical & sig_rstd_trialType_Outcome_logical;
sigO_pe_only = sig_rstd_PE_Outcome_logical & ~sig_rstd_trialType_Outcome_logical;
sigO_trialType_only = sig_rstd_trialType_Outcome_logical & ~sig_rstd_PE_Outcome_logical;

% OUTCOME: PE vs TRIAL TYPE
figure(6); 
subplot(2,2,1)
hold on;
scatter(flat_rstd_PE_Outcome_tStat_clean, flat_rstd_trialType_Outcome_tStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_PE_Outcome_tStat_clean(sigO_pe_only), flat_rstd_trialType_Outcome_tStat_clean(sigO_pe_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(1,:), 'MarkerEdgeColor','k');       % PE only - blue 
scatter(flat_rstd_PE_Outcome_tStat_clean(sigO_trialType_only), flat_rstd_trialType_Outcome_tStat_clean(sigO_trialType_only), 30, 'b', 'filled','o',  'MarkerFaceColor',  model_colorMap(5,:), 'MarkerEdgeColor','k');  % TrialType only - orange
scatter(flat_rstd_PE_Outcome_tStat_clean(sigO_both), flat_rstd_trialType_Outcome_tStat_clean(sigO_both), 30, 'y', 'filled','o', 'MarkerFaceColor',  model_colorMap(2,:), 'MarkerEdgeColor','k');             % Both - purple
% Labels and formatting
xlabel('PE Outcome t-Statistic');
ylabel('TrialType Outcome t-Statistic');
title('Scatter: PE vs TrialType Outcome t-Statistics');
xlim([-5 6]);
ylim([-5 5]);
grid on;
legend({'NS','Sig PE','Sig TrialType','Sig Both'}, 'Location', 'best');
hold off;
axis square

% Outcome: PE vs COLOR Orange
subplot(2,2,2)
% Combinations
sigO_both = sig_rstd_PE_Outcome_logical & sig_rstd_trialColor_Outcome_logical;
sigO_pe_only = sig_rstd_PE_Outcome_logical & ~sig_rstd_trialColor_Outcome_logical;
sig_trialColor_only = sig_rstd_trialColor_Outcome_logical & ~sig_rstd_PE_Outcome_logical;
hold on;
scatter(flat_rstd_PE_Outcome_tStat_clean, flat_rstd_trialColorO_Outcome_tStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_PE_Outcome_tStat_clean(sigO_pe_only), flat_rstd_trialColorO_Outcome_tStat_clean(sigO_pe_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(1,:), 'MarkerEdgeColor','k');       % PE only - blue 
scatter(flat_rstd_PE_Outcome_tStat_clean(sig_trialColor_only), flat_rstd_trialColorO_Outcome_tStat_clean(sig_trialColor_only), 30, 'b', 'filled','o',  'MarkerFaceColor', model_colorMap(3,:), 'MarkerEdgeColor','k');  % TrialColor only - pink
scatter(flat_rstd_PE_Outcome_tStat_clean(sigO_both), flat_rstd_trialColorO_Outcome_tStat_clean(sigO_both), 30, 'y', 'filled','o', 'MarkerFaceColor', model_colorMap(4,:), 'MarkerEdgeColor','k');             % Both - cyan
% Labels and formatting
xlabel('PE Outcome t-Statistic');
ylabel('TrialType Outcome t-Statistic');
title('Scatter: PE vs TrialColorO Outcome t-Statistics');
xlim([-5 6]);
ylim([-5 5]);
grid on;
legend({'NS','Sig PE','Sig TrialColor','Sig Both'}, 'Location', 'best');
hold off;
axis square

% Outcome: VE vs COLOR Red
subplot(2,2,3)
% Combinations
sigO_both = sig_rstd_PE_Outcome_logical & sig_rstd_trialColor_Outcome_logical;
sigO_pe_only = sig_rstd_PE_Outcome_logical & ~sig_rstd_trialColor_Outcome_logical;
sig_trialColor_only = sig_rstd_trialColor_Outcome_logical & ~sig_rstd_PE_Outcome_logical;
hold on;
scatter(flat_rstd_PE_Outcome_tStat_clean, flat_rstd_trialColorO_Outcome_tStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_PE_Outcome_tStat_clean(sigO_pe_only), flat_rstd_trialColorR_Outcome_tStat_clean(sigO_pe_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor',  model_colorMap(1,:), 'MarkerEdgeColor','k');       % PE only - green 
scatter(flat_rstd_PE_Outcome_tStat_clean(sig_trialColor_only), flat_rstd_trialColorR_Outcome_tStat_clean(sig_trialColor_only), 30, 'b', 'filled','o',  'MarkerFaceColor', model_colorMap(3,:), 'MarkerEdgeColor','k');  % TrialColor only - 
scatter(flat_rstd_PE_Outcome_tStat_clean(sigO_both), flat_rstd_trialColorR_Outcome_tStat_clean(sigO_both), 30, 'y', 'filled','o', 'MarkerFaceColor', model_colorMap(4,:), 'MarkerEdgeColor','k');             % Both - yellow
% Labels and formatting
xlabel('PE Outcome t-Statistic');
ylabel('TrialType Outcome t-Statistic');
title('Scatter: PE vs TrialColorR Outcome t-Statistics');
xlim([-5 6]);
ylim([-5 5]);
grid on;
legend({'NS','Sig PE','Sig TrialColor','Sig Both'}, 'Location', 'best');
hold off;
axis square

% Outcome: VE vs COLOR Yellow
subplot(2,2,4)
% Combinations
sigO_both = sig_rstd_PE_Outcome_logical & sig_rstd_trialColor_Outcome_logical;
sigO_pe_only = sig_rstd_PE_Outcome_logical & ~sig_rstd_trialColor_Outcome_logical;
sig_trialColor_only = sig_rstd_trialColor_Outcome_logical & ~sig_rstd_PE_Outcome_logical;
hold on;
scatter(flat_rstd_PE_Outcome_tStat_clean, flat_rstd_trialColorO_Outcome_tStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_PE_Outcome_tStat_clean(sigO_pe_only), flat_rstd_trialColorY_Outcome_tStat_clean(sigO_pe_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(1,:), 'MarkerEdgeColor','k');       % PE only - green 
scatter(flat_rstd_PE_Outcome_tStat_clean(sig_trialColor_only), flat_rstd_trialColorY_Outcome_tStat_clean(sig_trialColor_only), 30, 'b', 'filled','o',  'MarkerFaceColor', model_colorMap(3,:), 'MarkerEdgeColor','k');  % TrialColor only - 
scatter(flat_rstd_PE_Outcome_tStat_clean(sigO_both), flat_rstd_trialColorY_Outcome_tStat_clean(sigO_both), 30, 'y', 'filled','o', 'MarkerFaceColor', model_colorMap(4,:), 'MarkerEdgeColor','k');             % Both - yellow
% Labels and formatting
xlabel('PE Outcome t-Statistic');
ylabel('TrialType Outcome t-Statistic');
title('Scatter: PE vs TrialColorY Outcome t-Statistics');
xlim([-5 6]);
ylim([-5 5]);
grid on;
legend({'NS','Sig PE','Sig TrialColor','Sig Both'}, 'Location', 'best');
hold off;
axis square

% maybe add outcome_outcome if needed.
% check color is correct.. tstats seem low.

saveas(6,'Z:\Data\Rhiannon\RSTD_BART\mFig_tstats_OutcomeModels.pdf') % save mFig
close(6)


%% F-STATISTICS
%% Plot f-statistics

rstd_VE_Cue_fStat = rstd_VE_Cue_fStat(include_idx,:);
rstd_trialColor_Cue_fStat = rstd_trialColor_Cue_fStat(include_idx,:);
rstd_trialType_Cue_fStat = rstd_trialType_Cue_fStat(include_idx,:);
rstd_PE_Outcome_fStat = rstd_PE_Outcome_fStat(include_idx,:);
rstd_trialColor_Outcome_fStat = rstd_trialColor_Outcome_fStat(include_idx,:);
rstd_trialType_Outcome_fStat = rstd_trialType_Outcome_fStat(include_idx,:);
rstd_outcome_Outcome_fStat = rstd_outcome_Outcome_fStat(include_idx,:);

% Initialize output matrix with NaNs
padded_rstd_VE_Cue_fStat = NaN(size(rstd_VE_Cue_fStat));  % size: 69x119
padded_rstd_trialColor_Cue_fStat = NaN(size(rstd_trialColor_Cue_fStat));  % size: 69x119
padded_rstd_trialType_Cue_fStat = NaN(size(rstd_trialType_Cue_fStat));  % size: 69x119
padded_rstd_PE_Outcome_fStat = NaN(size(rstd_PE_Outcome_fStat));  % size: 69x119
padded_rstd_trialColor_Outcome_fStat = NaN(size(rstd_trialColor_Outcome_fStat));  % size: 69x119
padded_rstd_trialType_Outcome_fStat = NaN(size(rstd_trialType_Outcome_fStat));  % size: 69x119
padded_rstd_outcome_Outcome_fStat = NaN(size(rstd_outcome_Outcome_fStat));  % size: 69x119
% Loop over each row to add in NaNs
for i = 1:length(lengthElectrodes_PerPatient)
    nCols = lengthElectrodes_PerPatient(i);
    padded_rstd_VE_Cue_fStat(i, 1:nCols) = rstd_VE_Cue_fStat(i, 1:nCols);
    padded_rstd_trialColor_Cue_fStat(i, 1:nCols) = rstd_trialColor_Cue_fStat(i, 1:nCols);
    padded_rstd_trialType_Cue_fStat(i, 1:nCols) = rstd_trialType_Cue_fStat(i, 1:nCols);
    padded_rstd_PE_Outcome_fStat(i, 1:nCols) = rstd_PE_Outcome_fStat(i, 1:nCols);
    padded_rstd_trialColor_Outcome_fStat(i, 1:nCols) = rstd_trialColor_Outcome_fStat(i, 1:nCols);
    padded_rstd_trialType_Outcome_fStat(i, 1:nCols) = rstd_trialType_Outcome_fStat(i, 1:nCols);
    padded_rstd_outcome_Outcome_fStat(i, 1:nCols) = rstd_outcome_Outcome_fStat(i, 1:nCols);
end
% Flatten each matrix row-wise (patient-wise)
flat_rstd_VE_Cue_fStat              = reshape(padded_rstd_VE_Cue_fStat', [], 1);
flat_rstd_trialColor_Cue_fStat      = reshape(padded_rstd_trialColor_Cue_fStat', [], 1);
flat_rstd_trialType_Cue_fStat              = reshape(padded_rstd_trialType_Cue_fStat', [], 1);
flat_rstd_PE_Outcome_fStat      = reshape(padded_rstd_PE_Outcome_fStat', [], 1);
flat_rstd_trialColor_Outcome_fStat       = reshape(padded_rstd_trialColor_Outcome_fStat', [], 1);
flat_rstd_trialType_Outcome_fStat       = reshape(padded_rstd_trialType_Outcome_fStat', [], 1);
flat_rstd_outcome_Outcome_fStat      = reshape(padded_rstd_outcome_Outcome_fStat', [], 1);
% Drop NaN rows
flat_rstd_VE_Cue_fStat_clean = flat_rstd_VE_Cue_fStat(~any(isnan(flat_rstd_VE_Cue_fStat), 2), :);
flat_rstd_trialColor_Cue_fStat_clean = flat_rstd_trialColor_Cue_fStat(~any(isnan(flat_rstd_trialColor_Cue_fStat), 2), :);
flat_rstd_trialType_Cue_fStat_clean = flat_rstd_trialType_Cue_fStat(~any(isnan(flat_rstd_trialType_Cue_fStat), 2), :);
flat_rstd_PE_Outcome_fStat_clean = flat_rstd_PE_Outcome_fStat(~any(isnan(flat_rstd_PE_Outcome_fStat), 2), :);
flat_rstd_trialColor_Outcome_fStat_clean = flat_rstd_trialColor_Outcome_fStat(~any(isnan(flat_rstd_trialColor_Outcome_fStat), 2), :);
flat_rstd_trialType_Outcome_fStat_clean = flat_rstd_trialType_Outcome_fStat(~any(isnan(flat_rstd_trialType_Outcome_fStat), 2), :);
flat_rstd_outcome_Outcome_fStat_clean = flat_rstd_outcome_Outcome_fStat(~any(isnan(flat_rstd_outcome_Outcome_fStat), 2), :);

% CUE: VE vs TRIAL TYPE (F-STATISTIC)
figure(7); 
subplot(1,2,1)
hold on;
scatter(flat_rstd_VE_Cue_fStat_clean, flat_rstd_trialType_Cue_fStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_VE_Cue_fStat_clean(sigC_ve_only), flat_rstd_trialType_Cue_fStat_clean(sigC_ve_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(1,:), 'MarkerEdgeColor','k'); % VE only - dark blue 
scatter(flat_rstd_VE_Cue_fStat_clean(sigC_trialType_only), flat_rstd_trialType_Cue_fStat_clean(sigC_trialType_only), 30, 'b', 'filled','o',  'MarkerFaceColor', model_colorMap(5,:), 'MarkerEdgeColor','k'); % TrialType only - orange
scatter(flat_rstd_VE_Cue_fStat_clean(sigC_both), flat_rstd_trialType_Cue_fStat_clean(sigC_both), 30, 'y', 'filled','o', 'MarkerFaceColor', model_colorMap(2,:), 'MarkerEdgeColor','k'); % Both - purple
% Labels and formatting
xlabel('VE Cue f-Statistic');
ylabel('TrialType Cue f-Statistic');
title('Scatter: VE vs TrialType Cue f-Statistics');
xlim([0 100]);
ylim([0 100]);
grid on;
legend({'NS','Sig VE','Sig TrialType','Sig Both'}, 'Location', 'best');
hold off;
axis square
% CUE: VE vs COLOR (F-STATISTIC)
subplot(1,2,2)
% Combinations
sigC_both = sig_rstd_VE_Cue_logical & sig_rstd_trialColor_Cue_logical;
sigC_ve_only = sig_rstd_VE_Cue_logical & ~sig_rstd_trialColor_Cue_logical;
sig_trialColor_only = sig_rstd_trialColor_Cue_logical & ~sig_rstd_VE_Cue_logical;
hold on;
scatter(flat_rstd_VE_Cue_fStat_clean, flat_rstd_trialColor_Cue_fStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_VE_Cue_fStat_clean(sigC_ve_only), flat_rstd_trialColor_Cue_fStat_clean(sigC_ve_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(1,:), 'MarkerEdgeColor','k');       % VE only - blue 
scatter(flat_rstd_VE_Cue_fStat_clean(sig_trialColor_only), flat_rstd_trialColor_Cue_fStat_clean(sig_trialColor_only), 30, 'b', 'filled','o',  'MarkerFaceColor',  model_colorMap(3,:), 'MarkerEdgeColor','k');  % TrialColor only - pink
scatter(flat_rstd_VE_Cue_fStat_clean(sigC_both), flat_rstd_trialColor_Cue_fStat_clean(sigC_both), 30, 'y', 'filled','o', 'MarkerFaceColor',  model_colorMap(4,:), 'MarkerEdgeColor','k');             % Both - coral
% Labels and formatting
xlabel('VE Cue f-Statistic');
ylabel('TrialType Cue f-Statistic');
title('Scatter: VE vs TrialColor Cue f-Statistics');
xlim([0 100]);
ylim([0 100]);
grid on;
legend({'NS','Sig VE','Sig TrialColor','Sig Both'}, 'Location', 'best');
hold off;
axis square

saveas(7,'Z:\Data\Rhiannon\RSTD_BART\mFig_fstats_CueModels.pdf') % save mFig
close(7)

% OUTCOME: PE vs TRIAL TYPE (F-STATISTIC)
% need to create with outcome too^
sigO_both = sig_rstd_PE_Outcome_logical & sig_rstd_trialType_Outcome_logical;
sigO_pe_only = sig_rstd_PE_Outcome_logical & ~sig_rstd_trialType_Outcome_logical;
sigO_trialType_only = sig_rstd_trialType_Outcome_logical & ~sig_rstd_PE_Outcome_logical;
sigO_trialColor_only = sig_rstd_trialColor_Outcome_logical & ~sig_rstd_PE_Outcome_logical;

figure(8); 
subplot(1,3,1)
hold on;
scatter(flat_rstd_PE_Outcome_fStat_clean, flat_rstd_trialType_Outcome_fStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_PE_Outcome_fStat_clean(sigO_pe_only), flat_rstd_trialType_Outcome_fStat_clean(sigO_pe_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(6,:), 'MarkerEdgeColor','k'); % PE only - green 
scatter(flat_rstd_PE_Outcome_fStat_clean(sigO_trialType_only), flat_rstd_trialType_Outcome_fStat_clean(sigO_trialType_only), 30, 'b', 'filled','o',  'MarkerFaceColor', model_colorMap(5,:), 'MarkerEdgeColor','k'); % TrialType only - Pink
scatter(flat_rstd_PE_Outcome_fStat_clean(sigO_both), flat_rstd_trialType_Outcome_fStat_clean(sigO_both), 30, 'y', 'filled','o', 'MarkerFaceColor', model_colorMap(2,:), 'MarkerEdgeColor','k'); % Both - purple
% Labels and formatting
xlabel('PE Outcome f-Statistic');
ylabel('TrialType Outcome f-Statistic');
title('Scatter: PE vs TrialType Outcome f-Statistics');
xlim([0 100]);
ylim([0 100]);
grid on;
legend({'NS','Sig PE','Sig TrialType','Sig Both'}, 'Location', 'best');
hold off;
axis square
% OUTCOME: PE vs COLOR (F-STATISTIC)
subplot(1,3,2)
% Combinations
sigO_both = sig_rstd_PE_Outcome_logical & sig_rstd_trialColor_Outcome_logical;
sigO_pe_only = sig_rstd_PE_Outcome_logical & ~sig_rstd_trialColor_Outcome_logical;
sigO_trialColor_only = sig_rstd_trialColor_Outcome_logical & ~sig_rstd_PE_Outcome_logical;
hold on;
scatter(flat_rstd_PE_Outcome_fStat_clean, flat_rstd_trialColor_Outcome_fStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_PE_Outcome_fStat_clean(sigO_pe_only), flat_rstd_trialColor_Outcome_fStat_clean(sigO_pe_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(6,:), 'MarkerEdgeColor','k');       % PE only -  green 
scatter(flat_rstd_PE_Outcome_fStat_clean(sigO_trialColor_only), flat_rstd_trialColor_Outcome_fStat_clean(sigO_trialColor_only), 30, 'b', 'filled','o',  'MarkerFaceColor',  model_colorMap(3,:), 'MarkerEdgeColor','k');  % TrialColor only - 
scatter(flat_rstd_PE_Outcome_fStat_clean(sigO_both), flat_rstd_trialColor_Outcome_fStat_clean(sigO_both), 30, 'y', 'filled','o', 'MarkerFaceColor',  model_colorMap(4,:), 'MarkerEdgeColor','k');             % Both - 
% Labels and formatting
xlabel('PE Outcome f-Statistic');
ylabel('TrialColor Outcome f-Statistic');
title('Scatter: PE vs TrialColor Outcome f-Statistics');
xlim([0 100]);
ylim([0 100]);
grid on;
legend({'NS','Sig PE','Sig TrialColor','Sig Both'}, 'Location', 'best');
hold off;
axis square
% OUTCOME: PE vs OUTCOME (F-STATISTIC)
subplot(1,3,3)
% Combinations
sigO_both = sig_rstd_PE_Outcome_logical & sig_rstd_Outcome_Outcome_logical;
sigO_pe_only = sig_rstd_PE_Outcome_logical & ~sig_rstd_Outcome_Outcome_logical;
sigO_outcome_only = sig_rstd_Outcome_Outcome_logical & ~sig_rstd_PE_Outcome_logical;
hold on;
scatter(flat_rstd_PE_Outcome_fStat_clean, flat_rstd_outcome_Outcome_fStat_clean, 30, [0.8 0.8 0.8], 'filled');
% Plot significant groups
scatter(flat_rstd_PE_Outcome_fStat_clean(sigO_pe_only), flat_rstd_outcome_Outcome_fStat_clean(sigO_pe_only), 30, 'r', 'filled', 'o', 'MarkerFaceColor', model_colorMap(6,:), 'MarkerEdgeColor','k');       % PE only -  green 
scatter(flat_rstd_PE_Outcome_fStat_clean(sigO_outcome_only), flat_rstd_outcome_Outcome_fStat_clean(sigO_outcome_only), 30, 'b', 'filled','o',  'MarkerFaceColor',  model_colorMap(7,:), 'MarkerEdgeColor','k');  % Outcome only - light blue
scatter(flat_rstd_PE_Outcome_fStat_clean(sigO_both), flat_rstd_outcome_Outcome_fStat_clean(sigO_both), 30, 'y', 'filled','o', 'MarkerFaceColor',  model_colorMap(8,:), 'MarkerEdgeColor','k');             % Both - marroon
% Labels and formatting
xlabel('PE Outcome f-Statistic');
ylabel('Outcome Outcome f-Statistic');
title('Scatter: PE vs Outcome Outcome f-Statistics');
xlim([0 200]);
ylim([0 200]);
grid on;
legend({'NS','Sig PE','Sig Outcome','Sig Both'}, 'Location', 'best');
hold off;
axis square

saveas(8,'Z:\Data\Rhiannon\RSTD_BART\mFig_fstats_OutcomeModels.pdf') % save mFig
close(8)

%% ~~~~~~~~~ model-derived risk sensitivity metric by # of contacts encoded ~~~~~~~~~~~~~~~~ %%
% regress rstdMetric by Cue Model # of significant contacts
% plot significant cue-aligned linear models:
figure(9)
% VE CUE
subplot(1,3,1)
lm = fitlm(sigVECue_ProportionPerPatient(include_idx),modelRS_Metric(include_idx)); % Significant
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Cue: VE prop')
ylabel('modelRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
ylim([-1 1])
% COLOR CUE
subplot(1,3,2)
lm = fitlm(sigColorCue_ProportionPerPatient(include_idx),modelRS_Metric(include_idx)); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Cue: Color prop')
ylabel('modelRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
ylim([-1 1])
% CUE TYPE
subplot(1,3,3)
lm = fitlm(sigTypeCue_ProportionPerPatient(include_idx),modelRS_Metric(include_idx)); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Cue: Type prop')
ylabel('modelRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
ylim([-1 1])

% test to remove the one outlier in color cue prop (NS)
% test_cuecolor= sigColorCue_ProportionPerPatient(include_idx)
% test_cuecolor_logical = test_cuecolor>0.6
% modelrstest =modelRS_Metric(include_idx)
% lm = fitlm(test_cuecolor(~test_cuecolor_logical),modelrstest(~test_cuecolor_logical)); % NS
% h = lm.plot;
% h(1).Marker = '.';
% h(1).MarkerEdgeColor = 'k';
% h(1).MarkerSize = 25;
% xlabel('Cue: Color prop')
% ylabel('modelRS_Metric')

saveas(9,'Z:\Data\Rhiannon\RSTD_BART\mFig_regressions_CueModels_modelRSMetric.pdf') % save mFig
close(9)

% VE ENCODING BY positive and negative LEARNING RATE
figure(10)
% positive PE
subplot(1,2,1)
lm = fitlm(sigVECue_ProportionPerPatient(include_idx),rstd_positivePE(include_idx)); % sig
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Cue: VE prop')
ylabel('positive alpha')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
ylim([-1 1])
% negative PE
subplot(1,2,2)
lm = fitlm(sigVECue_ProportionPerPatient(include_idx),rstd_negativePE(include_idx)); % sig
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Cue: VE prop')
ylabel('positive alpha')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
ylim([-1 1])

saveas(10,'Z:\Data\Rhiannon\RSTD_BART\mFig_regressions_CueModels_pos_neg_alphas.pdf') % save mFig
close(10)

% other regressions non-sig
lm = fitlm(sigColorCue_ProportionPerPatient(include_idx),rstd_positivePE(include_idx)); % NS
lm = fitlm(sigColorCue_ProportionPerPatient(include_idx),rstd_negativePE(include_idx)); % NS
lm = fitlm(sigTypeCue_ProportionPerPatient(include_idx),rstd_positivePE(include_idx)); % NS
lm = fitlm(sigTypeCue_ProportionPerPatient(include_idx),rstd_negativePE(include_idx)); % NS

%% ~~~~~~~~~~~~~~~~~~~~~~ (3) Outcome aligned Model regressions ~~~~~~~~~~~~~~~~~~~~ %%

% plot significant cue-aligned linear models:
figure(11)
% PE OUTCOME
subplot(1,4,1)
lm = fitlm(sigPEOutcome_ProportionPerPatient(include_idx),modelRS_Metric(include_idx)); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Outcome: PE prop')
ylabel('modelRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
ylim([-1 1])
% COLOR OUTCOME
subplot(1,4,2)
lm = fitlm(sigColorOutcome_ProportionPerPatient,modelRS_Metric); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Outcome: Color prop')
ylabel('modelRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
ylim([-1 1])
% OUTCOME TYPE
subplot(1,4,3)
lm = fitlm(sigTypeOutcome_ProportionPerPatient,modelRS_Metric); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Outcome: Type prop')
ylabel('modelRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
ylim([-1 1])
% OUTCOME OUTCOME
subplot(1,4,4)
lm = fitlm(sigOutcomeOutcome_ProportionPerPatient,modelRS_Metric); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Outcome: Outcome prop')
ylabel('modelRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
ylim([-1 1])

saveas(11,'Z:\Data\Rhiannon\RSTD_BART\mFig_regressions_OutcomeModels_Props.pdf') % save mFig
close(11)


%% ~~~~~~~~~ task-derived risk sensitivity metric by # of contacts encoded ~~~~~~~~~~~~~~~~ %%

% Does the encoding correlate with the behavioral RSTD metric?
% task-RS vs. model-RS (already in manuscript figure 1)
figure;
subplot(1,3,1)
lm = fitlm(taskRS_Metric,modelRS_Metric); % Significant
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 15;
xlabel('taskRS_Metric')
ylabel('modelRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
subplot(1,3,2)
lm = fitlm(taskRS_Metric,rstd_positivePE); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 15;
xlabel('taskRS_Metric')
ylabel('positive Alpha')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
subplot(1,3,3)
lm = fitlm(taskRS_Metric,rstd_negativePE); % Significant
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 15;
xlabel('taskRS_Metric')
ylabel('negative Alpha')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square

% dont save. 

% Cue Model significant contacts vs.taskRS_Metric
figure(12)
% VE CUE
subplot(1,3,1)
lm = fitlm(sigVECue_ProportionPerPatient(include_idx),taskRS_Metric(include_idx)); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 15;
xlabel('Cue: VE prop')
ylabel('taskRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
% COLOR CUE
subplot(1,3,2)
lm = fitlm(sigColorCue_ProportionPerPatient,taskRS_Metric); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 15;
xlabel('Cue: Color prop')
ylabel('taskRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
% CUE TYPE
subplot(1,3,3)
lm = fitlm(sigTypeCue_ProportionPerPatient,taskRS_Metric); % Significant
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 15;
xlabel('Cue: Type prop')
ylabel('taskRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square

saveas(12,'Z:\Data\Rhiannon\RSTD_BART\mFig_regressions_CueModels_taskRS.pdf') % save mFig
close(12)

% Outcome Model significant contacts vs. taskRS_Metric
figure(13)
% PE Outcome
subplot(1,4,1)
lm = fitlm(sigPEOutcome_ProportionPerPatient(include_idx),taskRS_Metric(include_idx)); % 
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Outcome: PE prop')
ylabel('taskRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
% COLOR Outcome
subplot(1,4,2)
lm = fitlm(sigColorOutcome_ProportionPerPatient,taskRS_Metric); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Outcome: Color prop')
ylabel('taskRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
% OUTCOME TYPE
subplot(1,4,3)
lm = fitlm(sigTypeOutcome_ProportionPerPatient,taskRS_Metric); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Outcome: Type prop')
ylabel('taskRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square
% OUTCOME OUTCOME
subplot(1,4,4)
lm = fitlm(sigOutcomeOutcome_ProportionPerPatient,taskRS_Metric); % NS
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'k';
h(1).MarkerSize = 25;
xlabel('Outcome: Outcome prop')
ylabel('taskRS_Metric')
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue))
legend off
axis tight square

saveas(13,'Z:\Data\Rhiannon\RSTD_BART\mFig_regressions_OutcomeModels_taskRS.pdf') % save mFig
close(13)


%% ~~~~~~~~~~~~~~ BRAIN REGIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

% Figure out where in the brain these significant contacts are:

%  4a) incorporate new contact labels
for k = 1:nPts
    for j = 1:length(neuralData_All{k}.TDdataGradients.neuralFit)
        trodeLabels_All{k,j} = neuralData_All{k}.TDdataGradients.neuralFit(j).new_trodeLabel;
    end
end

% (a) How many contacts in each brain region?
% number of contacts
% (A) Amygdala
amygdala_Contacts = strcmp(trodeLabels_All, 'Amygdala');
sum(sum(amygdala_Contacts)) % = 221
% (B) Cingulate
cingulate_Contacts = strcmp(trodeLabels_All, 'Cingulate');
sum(sum(cingulate_Contacts)) % = 122
% (C) Hippocampus
hippocampus_Contacts = strcmp(trodeLabels_All, 'Hippocampus');
sum(sum(hippocampus_Contacts)) % = 616
% (D) Inf.Frontal
InfFrontal_Contacts = strcmp(trodeLabels_All, 'InfFrontal');
sum(sum(InfFrontal_Contacts)) % = 188
% (E) Insula
insula_Contacts = strcmp(trodeLabels_All, 'Insula');
sum(sum(insula_Contacts)) % = 180
% (F) Med.Frontal
MedFrontal_Contacts = strcmp(trodeLabels_All, 'MedFrontal');
sum(sum(MedFrontal_Contacts)) % = 318
% (G) N.Accumbens
NAccumbens_Contacts = strcmp(trodeLabels_All, 'N.Accumbens');
sum(sum(NAccumbens_Contacts)) % = 5
% (H) Orb.Frontal
OrbFrontal_Contacts = strcmp(trodeLabels_All, 'OrbFrontal');
sum(sum(OrbFrontal_Contacts)) % = 314
% (I) Striatum
striatum_Contacts = strcmp(trodeLabels_All, 'Striatum');
sum(sum(striatum_Contacts)) % = 165
% (J) Sup.Frontal
SupFrontal_Contacts = strcmp(trodeLabels_All, 'SupFrontal');
sum(sum(SupFrontal_Contacts)) % = 149
% (K) Thalamus
thalamus_Contacts = strcmp(trodeLabels_All, 'Thalamus');
sum(sum(thalamus_Contacts)) % = 123
% (L) White Matter
whiteMatter_Contacts = strcmp(trodeLabels_All, 'WhiteMatter');
sum(sum(whiteMatter_Contacts)) % = 940
% (K) Unknown
unknown_Contacts = strcmp(trodeLabels_All, 'Unknown');
sum(sum(unknown_Contacts)) % = 64

% how many contacts total?
sum_allContacts = sum(sum(amygdala_Contacts)) +  sum(sum(cingulate_Contacts)) + sum(sum(hippocampus_Contacts)) + sum(sum(InfFrontal_Contacts)) + sum(sum(insula_Contacts))...
    + sum(sum(MedFrontal_Contacts)) + sum(sum(NAccumbens_Contacts)) + sum(sum(OrbFrontal_Contacts)) + sum(sum(striatum_Contacts)) + sum(sum(SupFrontal_Contacts))...
    + sum(sum(thalamus_Contacts)) + sum(sum(whiteMatter_Contacts)) + sum(sum(unknown_Contacts)); % 3405 (only the 11 ROIs)

% plotting total contacts by ROIs:
regions = {'Amygdala', 'Cingulate', 'Hippocampus', 'Inf Frontal', 'Insula', 'Med Frontal', 'NAccumbens', 'Orb Frontal', 'Striatum', 'Sup Frontal', 'Thalamus', 'White Matter', 'Unknown'};

% Sum the contacts for each region
totals = [sum(amygdala_Contacts(:)), sum(cingulate_Contacts(:)), sum(hippocampus_Contacts(:)), ...
    sum(InfFrontal_Contacts(:)), sum(insula_Contacts(:)), sum(MedFrontal_Contacts(:)), sum(NAccumbens_Contacts(:)), sum(OrbFrontal_Contacts(:)), ...
    sum(striatum_Contacts(:)), sum(SupFrontal_Contacts(:)), sum(thalamus_Contacts(:)), sum(whiteMatter_Contacts(:)), sum(unknown_Contacts(:))];

% Create the bar plot
figure;
bar(totals);
set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 45);
ylabel('Total Contacts');
title('Total Contacts per Brain Region');
for i = 1:length(totals)
    text(i, totals(i) + max(totals)*0.02, num2str(totals(i)), 'HorizontalAlignment', 'center');
end

%% ~~~~~~~~~~~~~~~~~   PROPORTION ENCODING IN THE BRAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

propBrainEncoding = true;

if propBrainEncoding
    % working through ROIs

    % plot details
    cue_labels = {'VE', 'Color', 'Type'};
    outcome_labels = {'PE', 'Color', 'Type', 'Outcome'};
    cueSuccess_labels = {'VE', 'Color', 'Type'};
    outcomeSuccess_labels = {'PE', 'Color', 'Type'};

    % (1) Amygdala Encoding:
    AMY_VE_Cue = sum(sum(amygdala_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(amygdala_Contacts))*100;
    AMY_Color_Cue = sum(sum(amygdala_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(amygdala_Contacts))*100;
    AMY_Type_Cue = sum(sum(amygdala_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(amygdala_Contacts))*100;

    AMY_PE_Outcome = sum(sum(amygdala_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(amygdala_Contacts))*100;
    AMY_Color_Outcome = sum(sum(amygdala_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(amygdala_Contacts))*100;
    AMY_Type_Outcome = sum(sum(amygdala_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(amygdala_Contacts))*100;
    AMY_Outcome_Outcome = sum(sum(amygdala_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(amygdala_Contacts))*100;

    AMY_VE_CueSuccess = sum(sum(amygdala_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(amygdala_Contacts))*100;
    AMY_Color_CueSuccess = sum(sum(amygdala_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(amygdala_Contacts))*100;
    AMY_Type_CueSuccess = sum(sum(amygdala_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(amygdala_Contacts))*100;

    AMY_PE_OutcomeSuccess = sum(sum(amygdala_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(amygdala_Contacts))*100;
    AMY_Color_OutcomeSuccess = sum(sum(amygdala_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(amygdala_Contacts))*100;
    AMY_Type_OutcomeSuccess = sum(sum(amygdala_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(amygdala_Contacts))*100;

    % Plotting amygdala encoding of task variables
    % Labels and data
    cue_data = [AMY_VE_Cue, AMY_Color_Cue, AMY_Type_Cue];
    outcome_data = [AMY_PE_Outcome, AMY_Color_Outcome, AMY_Type_Outcome, AMY_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% Amygdala Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% Amygdala Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Amygdala Encoding of Task Variables');
    % SUCCESS MODELS ENCODING
    cueSuccess_data = [AMY_VE_CueSuccess, AMY_Color_CueSuccess, AMY_Type_CueSuccess];
    outcomeSuccess_data = [AMY_PE_OutcomeSuccess, AMY_Color_OutcomeSuccess, AMY_Type_OutcomeSuccess];
    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% Amygdala Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% Amygdala Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Amygdala Encoding of Task Variables');

    % (2) Cingulate Encoding:
    CING_VE_Cue = sum(sum(cingulate_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(cingulate_Contacts)) * 100;
    CING_Color_Cue = sum(sum(cingulate_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(cingulate_Contacts)) * 100;
    CING_Type_Cue = sum(sum(cingulate_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(cingulate_Contacts)) * 100;

    CING_PE_Outcome = sum(sum(cingulate_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(cingulate_Contacts)) * 100;
    CING_Color_Outcome = sum(sum(cingulate_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(cingulate_Contacts)) * 100;
    CING_Type_Outcome = sum(sum(cingulate_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(cingulate_Contacts)) * 100;
    CING_Outcome_Outcome = sum(sum(cingulate_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(cingulate_Contacts)) * 100;

    CING_VE_CueSuccess = sum(sum(cingulate_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(cingulate_Contacts)) * 100;
    CING_Color_CueSuccess = sum(sum(cingulate_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(cingulate_Contacts)) * 100;
    CING_Type_CueSuccess = sum(sum(cingulate_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(cingulate_Contacts)) * 100;

    CING_PE_OutcomeSuccess = sum(sum(cingulate_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(cingulate_Contacts)) * 100;
    CING_Color_OutcomeSuccess = sum(sum(cingulate_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(cingulate_Contacts)) * 100;
    CING_Type_OutcomeSuccess = sum(sum(cingulate_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(cingulate_Contacts)) * 100;

    % Plotting cingulate encoding of task variables
    % Labels and data
    cue_data = [CING_VE_Cue, CING_Color_Cue, CING_Type_Cue];
    outcome_data = [CING_PE_Outcome, CING_Color_Outcome, CING_Type_Outcome, CING_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% Cingulate Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% Cingulate Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Cingulate Encoding of Task Variables');

    % SUCCESS MODELS ENCODING
    cueSuccess_data = [CING_VE_CueSuccess, CING_Color_CueSuccess, CING_Type_CueSuccess];
    outcomeSuccess_data = [CING_PE_OutcomeSuccess, CING_Color_OutcomeSuccess, CING_Type_OutcomeSuccess];

    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% Cingulate Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% Cingulate Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Cingulate Encoding of Task Variables');

    % (3) Hippocampus Encoding:
    HIPP_VE_Cue = sum(sum(hippocampus_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(hippocampus_Contacts)) * 100;
    HIPP_Color_Cue = sum(sum(hippocampus_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(hippocampus_Contacts)) * 100;
    HIPP_Type_Cue = sum(sum(hippocampus_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(hippocampus_Contacts)) * 100;

    HIPP_PE_Outcome = sum(sum(hippocampus_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(hippocampus_Contacts)) * 100;
    HIPP_Color_Outcome = sum(sum(hippocampus_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(hippocampus_Contacts)) * 100;
    HIPP_Type_Outcome = sum(sum(hippocampus_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(hippocampus_Contacts)) * 100;
    HIPP_Outcome_Outcome = sum(sum(hippocampus_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(hippocampus_Contacts)) * 100;

    HIPP_VE_CueSuccess = sum(sum(hippocampus_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(hippocampus_Contacts)) * 100;
    HIPP_Color_CueSuccess = sum(sum(hippocampus_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(hippocampus_Contacts)) * 100;
    HIPP_Type_CueSuccess = sum(sum(hippocampus_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(hippocampus_Contacts)) * 100;

    HIPP_PE_OutcomeSuccess = sum(sum(hippocampus_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(hippocampus_Contacts)) * 100;
    HIPP_Color_OutcomeSuccess = sum(sum(hippocampus_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(hippocampus_Contacts)) * 100;
    HIPP_Type_OutcomeSuccess = sum(sum(hippocampus_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(hippocampus_Contacts)) * 100;

    % Plotting hippocampus encoding of task variables
    % Labels and data
    cue_data = [HIPP_VE_Cue, HIPP_Color_Cue, HIPP_Type_Cue];
    outcome_data = [HIPP_PE_Outcome, HIPP_Color_Outcome, HIPP_Type_Outcome, HIPP_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% Hippocampus Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% Hippocampus Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Hippocampus Encoding of Task Variables');

    % SUCCESS MODELS ENCODING
    cueSuccess_data = [HIPP_VE_CueSuccess, HIPP_Color_CueSuccess, HIPP_Type_CueSuccess];
    outcomeSuccess_data = [HIPP_PE_OutcomeSuccess, HIPP_Color_OutcomeSuccess, HIPP_Type_OutcomeSuccess];

    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% Hippocampus Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% Hippocampus Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Hippocampus Encoding of Task Variables');

    % (4) Inferior Frontal Encoding:
    INF_VE_Cue = sum(sum(InfFrontal_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(InfFrontal_Contacts)) * 100;
    INF_Color_Cue = sum(sum(InfFrontal_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(InfFrontal_Contacts)) * 100;
    INF_Type_Cue = sum(sum(InfFrontal_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(InfFrontal_Contacts)) * 100;

    INF_PE_Outcome = sum(sum(InfFrontal_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(InfFrontal_Contacts)) * 100;
    INF_Color_Outcome = sum(sum(InfFrontal_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(InfFrontal_Contacts)) * 100;
    INF_Type_Outcome = sum(sum(InfFrontal_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(InfFrontal_Contacts)) * 100;
    INF_Outcome_Outcome = sum(sum(InfFrontal_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(InfFrontal_Contacts)) * 100;

    INF_VE_CueSuccess = sum(sum(InfFrontal_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(InfFrontal_Contacts)) * 100;
    INF_Color_CueSuccess = sum(sum(InfFrontal_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(InfFrontal_Contacts)) * 100;
    INF_Type_CueSuccess = sum(sum(InfFrontal_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(InfFrontal_Contacts)) * 100;

    INF_PE_OutcomeSuccess = sum(sum(InfFrontal_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(InfFrontal_Contacts)) * 100;
    INF_Color_OutcomeSuccess = sum(sum(InfFrontal_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(InfFrontal_Contacts)) * 100;
    INF_Type_OutcomeSuccess = sum(sum(InfFrontal_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(InfFrontal_Contacts)) * 100;

    % Plotting Inferior Frontal encoding of task variables
    % Labels and data
    cue_data = [INF_VE_Cue, INF_Color_Cue, INF_Type_Cue];
    outcome_data = [INF_PE_Outcome, INF_Color_Outcome, INF_Type_Outcome, INF_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% Inferior Frontal Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% Inferior Frontal Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Inferior Frontal Encoding of Task Variables');

    % SUCCESS MODELS ENCODING
    cueSuccess_data = [INF_VE_CueSuccess, INF_Color_CueSuccess, INF_Type_CueSuccess];
    outcomeSuccess_data = [INF_PE_OutcomeSuccess, INF_Color_OutcomeSuccess, INF_Type_OutcomeSuccess];

    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% Inferior Frontal Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% Inferior Frontal Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Inferior Frontal Encoding of Task Variables');

    % (5) Insula Encoding:
    INS_VE_Cue = sum(sum(insula_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(insula_Contacts)) * 100;
    INS_Color_Cue = sum(sum(insula_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(insula_Contacts)) * 100;
    INS_Type_Cue = sum(sum(insula_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(insula_Contacts)) * 100;

    INS_PE_Outcome = sum(sum(insula_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(insula_Contacts)) * 100;
    INS_Color_Outcome = sum(sum(insula_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(insula_Contacts)) * 100;
    INS_Type_Outcome = sum(sum(insula_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(insula_Contacts)) * 100;
    INS_Outcome_Outcome = sum(sum(insula_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(insula_Contacts)) * 100;

    INS_VE_CueSuccess = sum(sum(insula_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(insula_Contacts)) * 100;
    INS_Color_CueSuccess = sum(sum(insula_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(insula_Contacts)) * 100;
    INS_Type_CueSuccess = sum(sum(insula_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(insula_Contacts)) * 100;

    INS_PE_OutcomeSuccess = sum(sum(insula_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(insula_Contacts)) * 100;
    INS_Color_OutcomeSuccess = sum(sum(insula_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(insula_Contacts)) * 100;
    INS_Type_OutcomeSuccess = sum(sum(insula_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(insula_Contacts)) * 100;

    % Plotting Insula encoding of task variables
    % Labels and data
    cue_data = [INS_VE_Cue, INS_Color_Cue, INS_Type_Cue];
    outcome_data = [INS_PE_Outcome, INS_Color_Outcome, INS_Type_Outcome, INS_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% Insula Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% Insula Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Insula Encoding of Task Variables');

    % SUCCESS MODELS ENCODING
    cueSuccess_data = [INS_VE_CueSuccess, INS_Color_CueSuccess, INS_Type_CueSuccess];
    outcomeSuccess_data = [INS_PE_OutcomeSuccess, INS_Color_OutcomeSuccess, INS_Type_OutcomeSuccess];

    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% Insula Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% Insula Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Insula Encoding of Task Variables');

    % (6) Medial Frontal Encoding:
    MDF_VE_Cue = sum(sum(MedFrontal_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(MedFrontal_Contacts)) * 100;
    MDF_Color_Cue = sum(sum(MedFrontal_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(MedFrontal_Contacts)) * 100;
    MDF_Type_Cue = sum(sum(MedFrontal_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(MedFrontal_Contacts)) * 100;

    MDF_PE_Outcome = sum(sum(MedFrontal_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(MedFrontal_Contacts)) * 100;
    MDF_Color_Outcome = sum(sum(MedFrontal_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(MedFrontal_Contacts)) * 100;
    MDF_Type_Outcome = sum(sum(MedFrontal_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(MedFrontal_Contacts)) * 100;
    MDF_Outcome_Outcome = sum(sum(MedFrontal_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(MedFrontal_Contacts)) * 100;

    MDF_VE_CueSuccess = sum(sum(MedFrontal_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(MedFrontal_Contacts)) * 100;
    MDF_Color_CueSuccess = sum(sum(MedFrontal_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(MedFrontal_Contacts)) * 100;
    MDF_Type_CueSuccess = sum(sum(MedFrontal_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(MedFrontal_Contacts)) * 100;

    MDF_PE_OutcomeSuccess = sum(sum(MedFrontal_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(MedFrontal_Contacts)) * 100;
    MDF_Color_OutcomeSuccess = sum(sum(MedFrontal_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(MedFrontal_Contacts)) * 100;
    MDF_Type_OutcomeSuccess = sum(sum(MedFrontal_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(MedFrontal_Contacts)) * 100;

    % Plotting Medial Frontal encoding of task variables
    % Labels and data
    cue_data = [MDF_VE_Cue, MDF_Color_Cue, MDF_Type_Cue];
    outcome_data = [MDF_PE_Outcome, MDF_Color_Outcome, MDF_Type_Outcome, MDF_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% Medial Frontal Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% Medial Frontal Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Medial Frontal Encoding of Task Variables');

    % SUCCESS MODELS ENCODING
    cueSuccess_data = [MDF_VE_CueSuccess, MDF_Color_CueSuccess, MDF_Type_CueSuccess];
    outcomeSuccess_data = [MDF_PE_OutcomeSuccess, MDF_Color_OutcomeSuccess, MDF_Type_OutcomeSuccess];

    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% Medial Frontal Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% Medial Frontal Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Medial Frontal Encoding of Task Variables');

    % (7) Nucleus Accumbens Encoding:
    NAC_VE_Cue = sum(sum(NAccumbens_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(NAccumbens_Contacts)) * 100;
    NAC_Color_Cue = sum(sum(NAccumbens_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(NAccumbens_Contacts)) * 100;
    NAC_Type_Cue = sum(sum(NAccumbens_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(NAccumbens_Contacts)) * 100;

    NAC_PE_Outcome = sum(sum(NAccumbens_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(NAccumbens_Contacts)) * 100;
    NAC_Color_Outcome = sum(sum(NAccumbens_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(NAccumbens_Contacts)) * 100;
    NAC_Type_Outcome = sum(sum(NAccumbens_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(NAccumbens_Contacts)) * 100;
    NAC_Outcome_Outcome = sum(sum(NAccumbens_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(NAccumbens_Contacts)) * 100;

    NAC_VE_CueSuccess = sum(sum(NAccumbens_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(NAccumbens_Contacts)) * 100;
    NAC_Color_CueSuccess = sum(sum(NAccumbens_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(NAccumbens_Contacts)) * 100;
    NAC_Type_CueSuccess = sum(sum(NAccumbens_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(NAccumbens_Contacts)) * 100;

    NAC_PE_OutcomeSuccess = sum(sum(NAccumbens_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(NAccumbens_Contacts)) * 100;
    NAC_Color_OutcomeSuccess = sum(sum(NAccumbens_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(NAccumbens_Contacts)) * 100;
    NAC_Type_OutcomeSuccess = sum(sum(NAccumbens_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(NAccumbens_Contacts)) * 100;

    % Plotting Nucleus Accumbens encoding of task variables
    % Labels and data
    cue_data = [NAC_VE_Cue, NAC_Color_Cue, NAC_Type_Cue];
    outcome_data = [NAC_PE_Outcome, NAC_Color_Outcome, NAC_Type_Outcome, NAC_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% NAccumbens Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% NAccumbens Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('NAccumbens Encoding of Task Variables');

    % SUCCESS MODELS ENCODING
    cueSuccess_data = [NAC_VE_CueSuccess, NAC_Color_CueSuccess, NAC_Type_CueSuccess];
    outcomeSuccess_data = [NAC_PE_OutcomeSuccess, NAC_Color_OutcomeSuccess, NAC_Type_OutcomeSuccess];

    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% NAccumbens Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% NAccumbens Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('NAccumbens Encoding of Task Variables');

    % (8) Orbital Frontal Encoding:
    OFC_VE_Cue = sum(sum(OrbFrontal_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(OrbFrontal_Contacts)) * 100;
    OFC_Color_Cue = sum(sum(OrbFrontal_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(OrbFrontal_Contacts)) * 100;
    OFC_Type_Cue = sum(sum(OrbFrontal_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(OrbFrontal_Contacts)) * 100;

    OFC_PE_Outcome = sum(sum(OrbFrontal_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(OrbFrontal_Contacts)) * 100;
    OFC_Color_Outcome = sum(sum(OrbFrontal_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(OrbFrontal_Contacts)) * 100;
    OFC_Type_Outcome = sum(sum(OrbFrontal_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(OrbFrontal_Contacts)) * 100;
    OFC_Outcome_Outcome = sum(sum(OrbFrontal_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(OrbFrontal_Contacts)) * 100;

    OFC_VE_CueSuccess = sum(sum(OrbFrontal_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(OrbFrontal_Contacts)) * 100;
    OFC_Color_CueSuccess = sum(sum(OrbFrontal_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(OrbFrontal_Contacts)) * 100;
    OFC_Type_CueSuccess = sum(sum(OrbFrontal_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(OrbFrontal_Contacts)) * 100;

    OFC_PE_OutcomeSuccess = sum(sum(OrbFrontal_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(OrbFrontal_Contacts)) * 100;
    OFC_Color_OutcomeSuccess = sum(sum(OrbFrontal_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(OrbFrontal_Contacts)) * 100;
    OFC_Type_OutcomeSuccess = sum(sum(OrbFrontal_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(OrbFrontal_Contacts)) * 100;

    % Plotting Orbital Frontal encoding of task variables
    % Labels and data
    cue_data = [OFC_VE_Cue, OFC_Color_Cue, OFC_Type_Cue];
    outcome_data = [OFC_PE_Outcome, OFC_Color_Outcome, OFC_Type_Outcome, OFC_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% OrbFrontal Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% OrbFrontal Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('OrbFrontal Encoding of Task Variables');

    % SUCCESS MODELS ENCODING
    cueSuccess_data = [OFC_VE_CueSuccess, OFC_Color_CueSuccess, OFC_Type_CueSuccess];
    outcomeSuccess_data = [OFC_PE_OutcomeSuccess, OFC_Color_OutcomeSuccess, OFC_Type_OutcomeSuccess];

    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% OrbFrontal Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% OrbFrontal Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('OrbFrontal Encoding of Task Variables');


    % (9) Striatum Encoding:
    STR_VE_Cue = sum(sum(striatum_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(striatum_Contacts)) * 100;
    STR_Color_Cue = sum(sum(striatum_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(striatum_Contacts)) * 100;
    STR_Type_Cue = sum(sum(striatum_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(striatum_Contacts)) * 100;

    STR_PE_Outcome = sum(sum(striatum_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(striatum_Contacts)) * 100;
    STR_Color_Outcome = sum(sum(striatum_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(striatum_Contacts)) * 100;
    STR_Type_Outcome = sum(sum(striatum_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(striatum_Contacts)) * 100;
    STR_Outcome_Outcome = sum(sum(striatum_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(striatum_Contacts)) * 100;

    STR_VE_CueSuccess = sum(sum(striatum_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(striatum_Contacts)) * 100;
    STR_Color_CueSuccess = sum(sum(striatum_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(striatum_Contacts)) * 100;
    STR_Type_CueSuccess = sum(sum(striatum_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(striatum_Contacts)) * 100;

    STR_PE_OutcomeSuccess = sum(sum(striatum_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(striatum_Contacts)) * 100;
    STR_Color_OutcomeSuccess = sum(sum(striatum_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(striatum_Contacts)) * 100;
    STR_Type_OutcomeSuccess = sum(sum(striatum_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(striatum_Contacts)) * 100;

    % Plotting Striatum encoding of task variables
    % Labels and data
    cue_data = [STR_VE_Cue, STR_Color_Cue, STR_Type_Cue];
    outcome_data = [STR_PE_Outcome, STR_Color_Outcome, STR_Type_Outcome, STR_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% Striatum Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% Striatum Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Striatum Encoding of Task Variables');

    % SUCCESS MODELS ENCODING
    cueSuccess_data = [STR_VE_CueSuccess, STR_Color_CueSuccess, STR_Type_CueSuccess];
    outcomeSuccess_data = [STR_PE_OutcomeSuccess, STR_Color_OutcomeSuccess, STR_Type_OutcomeSuccess];

    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% Striatum Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% Striatum Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Striatum Encoding of Task Variables');

    % (1) Superior Frontal Encoding:
    SFR_VE_Cue = sum(sum(SupFrontal_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(SupFrontal_Contacts)) * 100;
    SFR_Color_Cue = sum(sum(SupFrontal_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(SupFrontal_Contacts)) * 100;
    SFR_Type_Cue = sum(sum(SupFrontal_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(SupFrontal_Contacts)) * 100;

    SFR_PE_Outcome = sum(sum(SupFrontal_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(SupFrontal_Contacts)) * 100;
    SFR_Color_Outcome = sum(sum(SupFrontal_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(SupFrontal_Contacts)) * 100;
    SFR_Type_Outcome = sum(sum(SupFrontal_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(SupFrontal_Contacts)) * 100;
    SFR_Outcome_Outcome = sum(sum(SupFrontal_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(SupFrontal_Contacts)) * 100;

    SFR_VE_CueSuccess = sum(sum(SupFrontal_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(SupFrontal_Contacts)) * 100;
    SFR_Color_CueSuccess = sum(sum(SupFrontal_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(SupFrontal_Contacts)) * 100;
    SFR_Type_CueSuccess = sum(sum(SupFrontal_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(SupFrontal_Contacts)) * 100;

    SFR_PE_OutcomeSuccess = sum(sum(SupFrontal_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(SupFrontal_Contacts)) * 100;
    SFR_Color_OutcomeSuccess = sum(sum(SupFrontal_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(SupFrontal_Contacts)) * 100;
    SFR_Type_OutcomeSuccess = sum(sum(SupFrontal_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(SupFrontal_Contacts)) * 100;

    % Plotting Superior Frontal encoding of task variables
    % Labels and data
    cue_data = [SFR_VE_Cue, SFR_Color_Cue, SFR_Type_Cue];
    outcome_data = [SFR_PE_Outcome, SFR_Color_Outcome, SFR_Type_Outcome, SFR_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% SupFrontal Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% SupFrontal Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('SupFrontal Encoding of Task Variables');

    % SUCCESS MODELS ENCODING
    cueSuccess_data = [SFR_VE_CueSuccess, SFR_Color_CueSuccess, SFR_Type_CueSuccess];
    outcomeSuccess_data = [SFR_PE_OutcomeSuccess, SFR_Color_OutcomeSuccess, SFR_Type_OutcomeSuccess];

    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% SupFrontal Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% SupFrontal Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('SupFrontal Encoding of Task Variables');


    % (1) Thalamus Encoding:
    % (1) Thalamus Encoding:
    THAL_VE_Cue = sum(sum(thalamus_Contacts & sig_rstd_VE_Cue)) ./  sum(sum(thalamus_Contacts)) * 100;
    THAL_Color_Cue = sum(sum(thalamus_Contacts & sig_rstd_trialColor_Cue)) ./  sum(sum(thalamus_Contacts)) * 100;
    THAL_Type_Cue = sum(sum(thalamus_Contacts & sig_rstd_trialType_Cue)) ./  sum(sum(thalamus_Contacts)) * 100;

    THAL_PE_Outcome = sum(sum(thalamus_Contacts & sig_rstd_PE_Outcome)) ./  sum(sum(thalamus_Contacts)) * 100;
    THAL_Color_Outcome = sum(sum(thalamus_Contacts & sig_rstd_trialColor_Outcome)) ./  sum(sum(thalamus_Contacts)) * 100;
    THAL_Type_Outcome = sum(sum(thalamus_Contacts & sig_rstd_trialType_Outcome)) ./  sum(sum(thalamus_Contacts)) * 100;
    THAL_Outcome_Outcome = sum(sum(thalamus_Contacts & sig_rstd_outcome_Outcome)) ./  sum(sum(thalamus_Contacts)) * 100;

    THAL_VE_CueSuccess = sum(sum(thalamus_Contacts & sig_rstd_VE_CueSuccess)) ./  sum(sum(thalamus_Contacts)) * 100;
    THAL_Color_CueSuccess = sum(sum(thalamus_Contacts & sig_rstd_trialColor_CueSuccess)) ./  sum(sum(thalamus_Contacts)) * 100;
    THAL_Type_CueSuccess = sum(sum(thalamus_Contacts & sig_rstd_trialType_CueSuccess)) ./  sum(sum(thalamus_Contacts)) * 100;

    THAL_PE_OutcomeSuccess = sum(sum(thalamus_Contacts & sig_rstd_PE_OutcomeSuccess)) ./  sum(sum(thalamus_Contacts)) * 100;
    THAL_Color_OutcomeSuccess = sum(sum(thalamus_Contacts & sig_rstd_trialColor_OutcomeSuccess)) ./  sum(sum(thalamus_Contacts)) * 100;
    THAL_Type_OutcomeSuccess = sum(sum(thalamus_Contacts & sig_rstd_trialType_OutcomeSuccess)) ./  sum(sum(thalamus_Contacts)) * 100;

    % Plotting Thalamus encoding of task variables
    % Labels and data
    cue_data = [THAL_VE_Cue, THAL_Color_Cue, THAL_Type_Cue];
    outcome_data = [THAL_PE_Outcome, THAL_Color_Outcome, THAL_Type_Outcome, THAL_Outcome_Outcome];
    figure;
    % Cue Panel
    subplot(2, 2, 1);
    bar(cue_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cue_labels);
    ylabel('% Thalamus Contacts');
    ylim([0 100]);
    title('Cue-Aligned Encoding');
    for i = 1:length(cue_data)
        text(i, cue_data(i) + 1.5, sprintf('%.1f%%', cue_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Panel
    subplot(2, 2, 2);
    bar(outcome_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcome_labels);
    ylabel('% Thalamus Contacts');
    ylim([0 100]);
    title('Outcome-Aligned Encoding');
    for i = 1:length(outcome_data)
        text(i, outcome_data(i) + 1.5, sprintf('%.1f%%', outcome_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Thalamus Encoding of Task Variables');

    % SUCCESS MODELS ENCODING
    cueSuccess_data = [THAL_VE_CueSuccess, THAL_Color_CueSuccess, THAL_Type_CueSuccess];
    outcomeSuccess_data = [THAL_PE_OutcomeSuccess, THAL_Color_OutcomeSuccess, THAL_Type_OutcomeSuccess];

    % Cue Success Panel
    subplot(2, 2, 3);
    bar(cueSuccess_data, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTickLabel', cueSuccess_labels);
    ylabel('% Thalamus Contacts');
    ylim([0 100]);
    title('Success Cue-Aligned Encoding');
    for i = 1:length(cueSuccess_data)
        text(i, cueSuccess_data(i) + 1.5, sprintf('%.1f%%', cueSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % Outcome Success Panel
    subplot(2, 2, 4);
    bar(outcomeSuccess_data, 'FaceColor', [0.8 0.4 0.4]);
    set(gca, 'XTickLabel', outcomeSuccess_labels);
    ylabel('% Thalamus Contacts');
    ylim([0 100]);
    title('Success Outcome-Aligned Encoding');
    for i = 1:length(outcomeSuccess_data)
        text(i, outcomeSuccess_data(i) + 1.5, sprintf('%.1f%%', outcomeSuccess_data(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    sgtitle('Thalamus Encoding of Task Variables');


    %% ~~~~~~~~~~~~~~~~~~~~~~~ Model Variables by Region ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    regions = {'Amygdala', 'Cingulate', 'Hippocampus', 'InfFrontal', 'Insula', ...
        'MedFrontal', 'NAccumbens', 'OrbFrontal', 'Striatum', 'SupFrontal', 'Thalamus'};

    % CUE MODEL REGIONAL ENCODING
    figure;
    VE_Cue_values = [AMY_VE_Cue, CING_VE_Cue, HIPP_VE_Cue, INF_VE_Cue, INS_VE_Cue, ...
        MDF_VE_Cue, NAC_VE_Cue, OFC_VE_Cue, STR_VE_Cue, SFR_VE_Cue, THAL_VE_Cue];
    subplot(2,3,1)
    bar(VE_Cue_values, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('VE Cue Encoding Across Brain Regions');
    for i = 1:length(VE_Cue_values)
        text(i, VE_Cue_values(i) + 1.5, sprintf('%.1f%%', VE_Cue_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    % Define Color_Cue values for each region
    Color_Cue_values = [AMY_Color_Cue, CING_Color_Cue, HIPP_Color_Cue, INF_Color_Cue, INS_Color_Cue, ...
        MDF_Color_Cue, NAC_Color_Cue, OFC_Color_Cue, STR_Color_Cue, SFR_Color_Cue, THAL_Color_Cue];
    % Plot
    subplot(2,3,2)
    bar(Color_Cue_values, 'FaceColor', [0.6 0.5 0.9]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('Color Cue Encoding Across Brain Regions');
    for i = 1:length(Color_Cue_values)
        text(i, Color_Cue_values(i) + 1.5, sprintf('%.1f%%', Color_Cue_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    % Define Type_Cue values for each region
    Type_Cue_values = [AMY_Type_Cue, CING_Type_Cue, HIPP_Type_Cue, INF_Type_Cue, INS_Type_Cue, ...
        MDF_Type_Cue, NAC_Type_Cue, OFC_Type_Cue, STR_Type_Cue, SFR_Type_Cue, THAL_Type_Cue];
    % Plot
    subplot(2,3,3)
    bar(Type_Cue_values, 'FaceColor', [0.85 0.6 0.2]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('Type Cue Encoding Across Brain Regions');
    for i = 1:length(Type_Cue_values)
        text(i, Type_Cue_values(i) + 1.5, sprintf('%.1f%%', Type_Cue_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    % CUE SUCCESS MODEL REGIONAL ENCODING
    VE_CueSuccess_values = [AMY_VE_CueSuccess, CING_VE_CueSuccess, HIPP_VE_CueSuccess, INF_VE_CueSuccess, INS_VE_CueSuccess, ...
        MDF_VE_CueSuccess, NAC_VE_CueSuccess, OFC_VE_CueSuccess, STR_VE_CueSuccess, SFR_VE_CueSuccess, THAL_VE_CueSuccess];
    subplot(2,3,4)
    bar(VE_CueSuccess_values, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('VE Cue Success Encoding Across Brain Regions');
    for i = 1:length(VE_CueSuccess_values)
        text(i, VE_CueSuccess_values(i) + 1.5, sprintf('%.1f%%', VE_CueSuccess_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    % Define Color_Cue values for each region
    Color_CueSuccess_values = [AMY_Color_CueSuccess, CING_Color_CueSuccess, HIPP_Color_CueSuccess, INF_Color_CueSuccess, INS_Color_CueSuccess, ...
        MDF_Color_CueSuccess, NAC_Color_CueSuccess, OFC_Color_CueSuccess, STR_Color_CueSuccess, SFR_Color_CueSuccess, THAL_Color_CueSuccess];
    % Plot
    subplot(2,3,5)
    bar(Color_CueSuccess_values, 'FaceColor', [0.6 0.5 0.9]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('Color Cue Success Encoding Across Brain Regions');
    for i = 1:length(Color_CueSuccess_values)
        text(i, Color_CueSuccess_values(i) + 1.5, sprintf('%.1f%%', Color_CueSuccess_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    % Define Type_Cue values for each region
    Type_CueSuccess_values = [AMY_Type_CueSuccess, CING_Type_CueSuccess, HIPP_Type_CueSuccess, INF_Type_CueSuccess, INS_Type_CueSuccess, ...
        MDF_Type_CueSuccess, NAC_Type_CueSuccess, OFC_Type_CueSuccess, STR_Type_CueSuccess, SFR_Type_CueSuccess, THAL_Type_CueSuccess];
    % Plot
    subplot(2,3,6)
    bar(Type_CueSuccess_values, 'FaceColor', [0.85 0.6 0.2]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('Type Cue Success Encoding Across Brain Regions');
    for i = 1:length(Type_CueSuccess_values)
        text(i, Type_CueSuccess_values(i) + 1.5, sprintf('%.1f%%', Type_CueSuccess_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    % OUTCOME ALIGNED MODELS by REGION
    PE_Outcome_values = [AMY_PE_Outcome, CING_PE_Outcome, HIPP_PE_Outcome, INF_PE_Outcome, INS_PE_Outcome, ...
        MDF_PE_Outcome, NAC_PE_Outcome, OFC_PE_Outcome, STR_PE_Outcome, SFR_PE_Outcome, THAL_PE_Outcome];
    figure;
    subplot(2,4,1)
    bar(PE_Outcome_values, 'FaceColor', [0.4 0.7 0.4]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('PE Outcome Encoding Across Brain Regions');
    for i = 1:length(PE_Outcome_values)
        text(i, PE_Outcome_values(i) + 1.5, sprintf('%.1f%%', PE_Outcome_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    Color_Outcome_values = [AMY_Color_Outcome, CING_Color_Outcome, HIPP_Color_Outcome, INF_Color_Outcome, INS_Color_Outcome, ...
        MDF_Color_Outcome, NAC_Color_Outcome, OFC_Color_Outcome, STR_Color_Outcome, SFR_Color_Outcome, THAL_Color_Outcome];
    subplot(2,4,2)
    bar(Color_Outcome_values, 'FaceColor', [0.6 0.5 0.9]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('Color Outcome Encoding Across Brain Regions');
    for i = 1:length(Color_Outcome_values)
        text(i, Color_Outcome_values(i) + 1.5, sprintf('%.1f%%', Color_Outcome_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    Type_Outcome_values = [AMY_Type_Outcome, CING_Type_Outcome, HIPP_Type_Outcome, INF_Type_Outcome, INS_Type_Outcome, ...
        MDF_Type_Outcome, NAC_Type_Outcome, OFC_Type_Outcome, STR_Type_Outcome, SFR_Type_Outcome, THAL_Type_Outcome];
    subplot(2,4,3)
    bar(Type_Outcome_values, 'FaceColor', [0.85 0.6 0.2]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('Type Outcome Encoding Across Brain Regions');
    for i = 1:length(Type_Outcome_values)
        text(i, Type_Outcome_values(i) + 1.5, sprintf('%.1f%%', Type_Outcome_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    Outcome_Outcome_values = [AMY_Outcome_Outcome, CING_Outcome_Outcome, HIPP_Outcome_Outcome, INF_Outcome_Outcome, INS_Outcome_Outcome, ...
        MDF_Outcome_Outcome, NAC_Outcome_Outcome, OFC_Outcome_Outcome, STR_Outcome_Outcome, SFR_Outcome_Outcome, THAL_Outcome_Outcome];
    subplot(2,4,4)
    bar(Outcome_Outcome_values, 'FaceColor', [0.1 0.2 0.7]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('Outcome Outcome Encoding Across Brain Regions');
    for i = 1:length(Outcome_Outcome_values)
        text(i, Outcome_Outcome_values(i) + 1.5, sprintf('%.1f%%', Outcome_Outcome_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    % OUTCOME SUCCESS MODELS
    PE_OutcomeSuccess_values = [AMY_PE_OutcomeSuccess, CING_PE_OutcomeSuccess, HIPP_PE_OutcomeSuccess, INF_PE_OutcomeSuccess, INS_PE_OutcomeSuccess, ...
        MDF_PE_OutcomeSuccess, NAC_PE_OutcomeSuccess, OFC_PE_OutcomeSuccess, STR_PE_OutcomeSuccess, SFR_PE_OutcomeSuccess, THAL_PE_OutcomeSuccess];
    subplot(2,4,5)
    bar(PE_OutcomeSuccess_values, 'FaceColor', [0.4 0.7 0.4]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('PE Outcome Success Encoding Across Brain Regions');
    for i = 1:length(PE_OutcomeSuccess_values)
        text(i, PE_OutcomeSuccess_values(i) + 1.5, sprintf('%.1f%%', PE_OutcomeSuccess_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    Color_OutcomeSuccess_values = [AMY_Color_OutcomeSuccess, CING_Color_OutcomeSuccess, HIPP_Color_OutcomeSuccess, INF_Color_OutcomeSuccess, INS_Color_OutcomeSuccess, ...
        MDF_Color_OutcomeSuccess, NAC_Color_OutcomeSuccess, OFC_Color_OutcomeSuccess, STR_Color_OutcomeSuccess, SFR_Color_OutcomeSuccess, THAL_Color_OutcomeSuccess];
    subplot(2,4,6)
    bar(Color_OutcomeSuccess_values, 'FaceColor', [0.6 0.5 0.9]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('Color Outcome Success Encoding Across Brain Regions');
    for i = 1:length(Color_OutcomeSuccess_values)
        text(i, Color_OutcomeSuccess_values(i) + 1.5, sprintf('%.1f%%', Color_OutcomeSuccess_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

    Type_OutcomeSuccess_values = [AMY_Type_OutcomeSuccess, CING_Type_OutcomeSuccess, HIPP_Type_OutcomeSuccess, INF_Type_OutcomeSuccess, INS_Type_OutcomeSuccess, ...
        MDF_Type_OutcomeSuccess, NAC_Type_OutcomeSuccess, OFC_Type_OutcomeSuccess, STR_Type_OutcomeSuccess, SFR_Type_OutcomeSuccess, THAL_Type_OutcomeSuccess];
    subplot(2,4,7)
    bar(Type_OutcomeSuccess_values, 'FaceColor', [0.85 0.6 0.2]);
    set(gca, 'XTick', 1:length(regions), 'XTickLabel', regions, 'XTickLabelRotation', 45);
    ylabel('Percentage of Contacts (%)');
    ylim([0 100]);
    title('Type Outcome Success Encoding Across Brain Regions');
    for i = 1:length(Type_OutcomeSuccess_values)
        text(i, Type_OutcomeSuccess_values(i) + 1.5, sprintf('%.1f%%', Type_OutcomeSuccess_values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    grid on;

else
end % if propBrainEncoding plots

%% ~~~~~~~~~~~~~~~ GLASS BRAIN TASK VARIABLE PROJECTIONS ~~~~~~~~~~~~~~~~~~~~ %%

% Load xzy coordinates to project onto the brain
% (2) Load Electrode Projections
final_ElecXYZProj = load('\\155.100.91.44\d\Data\Rhiannon\BART_RLDM_outputs\BrainRegions\final_ElecXYZProj');
ElecXYZProj_All = final_ElecXYZProj.final_ElecXYZProj;
load(fullfile(mainPath,ptID,'Imaging','Registered','Surfaces.mat')); % loading glass brain

% VE CUE PLOT:
figure(10)
subtitle('RSTD VE Cue')
hold on
VE_Cue_data = ElecXYZProj_All(sig_rstd_VE_Cue_logical,:);
myColormap = [0 0 0;    % black
    0 1 0];   % green
s1 = scatter3(VE_Cue_data(:,1), VE_Cue_data(:,2), VE_Cue_data(:,3), ...
    70, [0 1 0], 'filled', 'MarkerEdgeColor', 'k');
%colorbar;
colormap(myColormap);
view(3);
xlabel('Lateral --- Medial --- Lateral');
ylabel('Posterior --- Anterior');
zlabel('Inferior --- Superior');
% plotting actual glass brain....
patchs = patch('Faces', BrainSurfRaw.faces, 'Vertices', BrainSurfRaw.vertices, 'edgecolor', 'none', 'facecolor', 'flat', 'facealpha', .2);
facecolor = repmat([1 1 1],length(BrainSurfRaw.faces),1);
set(patchs, 'FaceVertexCData', facecolor);
camlight
lighting gouraud


% PE OUTCOME PLOT:
figure(11)
subtitle('RSTD PE Outcome')
hold on
PE_Outcome_data = ElecXYZProj_All(sig_rstd_PE_Outcome_logical,:);
s1 = scatter3(PE_Outcome_data(:,1), PE_Outcome_data(:,2), PE_Outcome_data(:,3), ...
    70, [1 0 0], 'filled', 'MarkerEdgeColor', 'k');
%colorbar;
view(3);
xlabel('Lateral --- Medial --- Lateral');
ylabel('Posterior --- Anterior');
zlabel('Inferior --- Superior');
% plotting actual glass brain....
patchs = patch('Faces', BrainSurfRaw.faces, 'Vertices', BrainSurfRaw.vertices, 'edgecolor', 'none', 'facecolor', 'flat', 'facealpha', .2);
facecolor = repmat([1 1 1],length(BrainSurfRaw.faces),1);
set(patchs, 'FaceVertexCData', facecolor);
camlight
lighting gouraud




%% ROIs glass brain plots

% (1) Flatten trodeLabels_All
trodeLabels_All_Flat = [];
for i = 1:size(trodeLabels_All, 1)
    rowLabels = trodeLabels_All(i, :);
    trodeLabels_All_Flat = [trodeLabels_All_Flat; rowLabels(:)];
end
trodeLabels_All_Flat = trodeLabels_All_Flat(~cellfun('isempty', trodeLabels_All_Flat));

% Getting region specific data
amygdala_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "Amygdala"), :);
cingulate_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "Cingulate"), :);
hippocampus_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "Hippocampus"), :);
infFrontal_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "InfFrontal"), :);
insula_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "Insula"), :);
nAccumbens_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "N.Accumbens"), :);
orbFrontal_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "OrbFrontal"), :);
striatum_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "Striatum"), :);
supFrontal_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "SupFrontal"), :);
thalamus_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "Thalamus"), :);

% get logical vars
isAmygdala = strcmp(trodeLabels_All_Flat, "Amygdala");
isCingulate = strcmp(trodeLabels_All_Flat, "Cingulate");
isHippocampus = strcmp(trodeLabels_All_Flat, "Hippocampus");
isInfFrontal = strcmp(trodeLabels_All_Flat, "InfFrontal");
isInsula = strcmp(trodeLabels_All_Flat, "Insula");
isNAccumbens = strcmp(trodeLabels_All_Flat, "N.Accumbens");
isOrbFrontal = strcmp(trodeLabels_All_Flat, "OrbFrontal");
isStriatum = strcmp(trodeLabels_All_Flat, "Striatum");
isSupFrontal = strcmp(trodeLabels_All_Flat, "SupFrontal");
isThalamus = strcmp(trodeLabels_All_Flat, "Thalamus");


% Creating index to plot single hemispheres
% Step 1: Find indices of vertices in the left hemisphere
left_idx = find(BrainSurfRaw.vertices(:,1) < 0);
right_idx = find(BrainSurfRaw.vertices(:,1) > 0);

% Step 2: Create a mask for faces in both hemis
face_mask_left = all(ismember(BrainSurfRaw.faces, left_idx), 2);
face_mask_right = all(ismember(BrainSurfRaw.faces, right_idx), 2);

% Step 3: Extract faces and remap vertex indices
faces_left = BrainSurfRaw.faces(face_mask_left, :);
faces_right = BrainSurfRaw.faces(face_mask_right, :);

% Step 4: Remap subset of vertices
[uniqueLeft_idx, ~, new_faces_left] = unique(faces_left(:));
[uniqueRight_idx, ~, new_faces_right] = unique(faces_right(:));
vertices_left = BrainSurfRaw.vertices(uniqueLeft_idx, :);
vertices_right = BrainSurfRaw.vertices(uniqueRight_idx, :);
faces_left = reshape(new_faces_left, size(faces_left));
faces_right = reshape(new_faces_right, size(faces_right));

% Get R2 values for each electrode
rstd_Cue_R2adj(sig_rstd_VE_Cue_logical)

% GETTING LOGICAL VARIABLES TO PROJECT DIFFERENT SIZES OF ELECTRODES ONTO BRAIN
% Initialize dotSize variable for Cue
dotSize_Cue = zeros(size(flat_rstd_Cue_R2adj_clean));

% Assign increasing dot sizes based on R² range for Cue
dotSize_Cue(flat_rstd_Cue_R2adj_clean <= 0.01) = 10;
dotSize_Cue(flat_rstd_Cue_R2adj_clean > 0.01 & flat_rstd_Cue_R2adj_clean <= 0.1) = 30;
dotSize_Cue(flat_rstd_Cue_R2adj_clean > 0.1 & flat_rstd_Cue_R2adj_clean <= 0.2) = 40;
dotSize_Cue(flat_rstd_Cue_R2adj_clean > 0.2 & flat_rstd_Cue_R2adj_clean <= 0.3) = 50;
dotSize_Cue(flat_rstd_Cue_R2adj_clean > 0.3 & flat_rstd_Cue_R2adj_clean <= 0.4) = 60;
dotSize_Cue(flat_rstd_Cue_R2adj_clean > 0.4 & flat_rstd_Cue_R2adj_clean <= 0.5) = 70;
dotSize_Cue(flat_rstd_Cue_R2adj_clean > 0.5 & flat_rstd_Cue_R2adj_clean <= 0.6) = 80;
dotSize_Cue(flat_rstd_Cue_R2adj_clean > 0.6) = 90;

% Initialize dotSize variable with zeros (same size as flat_rstd_CueSuccess_R2adj_clean)
dotSize_CueSuccess = zeros(size(flat_rstd_CueSuccess_R2adj_clean));

% Assign increasing dot sizes based on R² range
dotSize_CueSuccess(flat_rstd_CueSuccess_R2adj_clean <= 0.01) = 10;
dotSize_CueSuccess(flat_rstd_CueSuccess_R2adj_clean > 0.01 & flat_rstd_CueSuccess_R2adj_clean <= 0.1) = 30;
dotSize_CueSuccess(flat_rstd_CueSuccess_R2adj_clean > 0.1 & flat_rstd_CueSuccess_R2adj_clean <= 0.2) = 40;
dotSize_CueSuccess(flat_rstd_CueSuccess_R2adj_clean > 0.2 & flat_rstd_CueSuccess_R2adj_clean <= 0.3) = 50;
dotSize_CueSuccess(flat_rstd_CueSuccess_R2adj_clean > 0.3 & flat_rstd_CueSuccess_R2adj_clean <= 0.4) = 60;
dotSize_CueSuccess(flat_rstd_CueSuccess_R2adj_clean > 0.4 & flat_rstd_CueSuccess_R2adj_clean <= 0.5) = 70;
dotSize_CueSuccess(flat_rstd_CueSuccess_R2adj_clean > 0.5 & flat_rstd_CueSuccess_R2adj_clean <= 0.6) = 80;
dotSize_CueSuccess(flat_rstd_CueSuccess_R2adj_clean > 0.6) = 90;

% Initialize dotSize variable for Outcome
dotSize_Outcome = zeros(size(flat_rstd_Outcome_R2adj_clean));

% Assign increasing dot sizes based on R² range for Outcome
dotSize_Outcome(flat_rstd_Outcome_R2adj_clean <= 0.01) = 10;
dotSize_Outcome(flat_rstd_Outcome_R2adj_clean > 0.01 & flat_rstd_Outcome_R2adj_clean <= 0.1) = 30;
dotSize_Outcome(flat_rstd_Outcome_R2adj_clean > 0.1 & flat_rstd_Outcome_R2adj_clean <= 0.2) = 40;
dotSize_Outcome(flat_rstd_Outcome_R2adj_clean > 0.2 & flat_rstd_Outcome_R2adj_clean <= 0.3) = 50;
dotSize_Outcome(flat_rstd_Outcome_R2adj_clean > 0.3 & flat_rstd_Outcome_R2adj_clean <= 0.4) = 60;
dotSize_Outcome(flat_rstd_Outcome_R2adj_clean > 0.4 & flat_rstd_Outcome_R2adj_clean <= 0.5) = 70;
dotSize_Outcome(flat_rstd_Outcome_R2adj_clean > 0.5 & flat_rstd_Outcome_R2adj_clean <= 0.6) = 80;
dotSize_Outcome(flat_rstd_Outcome_R2adj_clean > 0.6) = 90;

% Initialize dotSize variable with zeros (same size as flat_rstd_OutcomeSuccess_R2adj_clean)
dotSize_OutcomeSuccess = zeros(size(flat_rstd_OutcomeSuccess_R2adj_clean));

% Assign increasing dot sizes based on R² range
dotSize_OutcomeSuccess(flat_rstd_OutcomeSuccess_R2adj_clean <= 0.01) = 10;
dotSize_OutcomeSuccess(flat_rstd_OutcomeSuccess_R2adj_clean > 0.01 & flat_rstd_OutcomeSuccess_R2adj_clean <= 0.1) = 30;
dotSize_OutcomeSuccess(flat_rstd_OutcomeSuccess_R2adj_clean > 0.1 & flat_rstd_OutcomeSuccess_R2adj_clean <= 0.2) = 40;
dotSize_OutcomeSuccess(flat_rstd_OutcomeSuccess_R2adj_clean > 0.2 & flat_rstd_OutcomeSuccess_R2adj_clean <= 0.3) = 50;
dotSize_OutcomeSuccess(flat_rstd_OutcomeSuccess_R2adj_clean > 0.3 & flat_rstd_OutcomeSuccess_R2adj_clean <= 0.4) = 60;
dotSize_OutcomeSuccess(flat_rstd_OutcomeSuccess_R2adj_clean > 0.4 & flat_rstd_OutcomeSuccess_R2adj_clean <= 0.5) = 70;
dotSize_OutcomeSuccess(flat_rstd_OutcomeSuccess_R2adj_clean > 0.5 & flat_rstd_OutcomeSuccess_R2adj_clean <= 0.6) = 80;
dotSize_OutcomeSuccess(flat_rstd_OutcomeSuccess_R2adj_clean > 0.6) = 90;

% CUE-ALIGNED ENCODING PROJECTIONS
figure(11)
subplot(1,2,1)
subtitle('Striatum Cue')
hold on
% all contacts in black
striatum_data = ElecXYZProj_All(strcmp(trodeLabels_All_Flat, "Striatum"), :);
s1 = scatter3(striatum_data(:,1), striatum_data(:,2), striatum_data(:,3), ...
    dotSize_Cue(isStriatum), [0 0 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.2); % black
% sig contacts for PE
isSigVE_Striatum = isStriatum & sig_rstd_VE_Cue_logical;
VE_CueSig_data = ElecXYZProj_All(isSigVE_Striatum, :);
s2 = scatter3(VE_CueSig_data(:,1), VE_CueSig_data(:,2), VE_CueSig_data(:,3), ...
    dotSize_Cue(isSigVE_Striatum), [1 0 0], 'filled', 'MarkerEdgeColor', 'k'); % red
% sig contacts for trialColor
isSigColor_Striatum = isStriatum & sig_rstd_trialColor_Cue_logical;
trialColor_CueSig_data = ElecXYZProj_All(isSigColor_Striatum, :);
s3 = scatter3(trialColor_CueSig_data(:,1), trialColor_CueSig_data(:,2), trialColor_CueSig_data(:,3), ...
    dotSize_Cue(isSigColor_Striatum), [1 1 0], 'filled', 'MarkerEdgeColor', 'k'); % yellow
% sig contacts for trialType
isSigType_Striatum = isStriatum & sig_rstd_trialType_Cue_logical;
trialType_CueSig_data = ElecXYZProj_All(isSigType_Striatum, :);
s4 = scatter3(trialType_CueSig_data(:,1), trialType_CueSig_data(:,2), trialType_CueSig_data(:,3), ...
    dotSize_Cue(isSigType_Striatum), [0 1 0], 'filled', 'MarkerEdgeColor', 'k'); % green
%colorbar;
%view(3);
%xlabel('Lateral --- Medial --- Lateral');
xlabel('Superior --- Inferior');
ylabel('Anterior  --- Posterior');
% plotting actual glass brain....
patchs = patch('Faces', BrainSurfRaw.faces, 'Vertices', BrainSurfRaw.vertices, 'edgecolor', 'none', 'facecolor', 'flat', 'facealpha', .2);
facecolor = repmat([1 1 1],length(BrainSurfRaw.faces),1);
set(patchs, 'FaceVertexCData', facecolor);
camlight
lighting gouraud
%view([-1 0 0]); % for medial view of right hemi -1
legend([s1, s2, s3, s4], ...
    {'All Contacts', 'VE', 'TrialColor', 'TrialType'}, ...
    'Location', 'bestoutside');

% OUTCOME-ALIGNED ENCODING PROJECTIONS
subplot(1,2,2)
subtitle('Striatum Outcome')
hold on
% all contacts in black
s1 = scatter3(striatum_data(:,1), striatum_data(:,2), striatum_data(:,3), ...
    70, [0 0 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.2); % black
% sig contacts for PE
isSigPE_Striatum = isStriatum & sig_rstd_PE_Outcome_logical;
PE_OutcomeSig_data = ElecXYZProj_All(isSigPE_Striatum, :);
s2 = scatter3(PE_OutcomeSig_data(:,1), PE_OutcomeSig_data(:,2), PE_OutcomeSig_data(:,3), ...
    70, [1 0 0], 'filled', 'MarkerEdgeColor', 'k'); % red
% sig contacts for trialColor
isSigColor_Striatum = isStriatum & sig_rstd_trialColor_Outcome_logical;
trialColor_OutcomeSig_data = ElecXYZProj_All(isSigColor_Striatum, :);
s3 = scatter3(trialColor_OutcomeSig_data(:,1), trialColor_OutcomeSig_data(:,2), trialColor_OutcomeSig_data(:,3), ...
    70, [1 1 0], 'filled', 'MarkerEdgeColor', 'k'); % yellow
% sig contacts for trialType
isSigType_Striatum = isStriatum & sig_rstd_trialType_Outcome_logical;
trialType_OutcomeSig_data = ElecXYZProj_All(isSigType_Striatum, :);
s4 = scatter3(trialType_OutcomeSig_data(:,1), trialType_OutcomeSig_data(:,2), trialType_OutcomeSig_data(:,3), ...
    70, [0 1 0], 'filled', 'MarkerEdgeColor', 'k'); % green
% sig contacts for outcome
isSigOutcome_Striatum = isStriatum & sig_rstd_Outcome_Outcome_logical;
Outcome_OutcomeSig_data = ElecXYZProj_All(isSigOutcome_Striatum, :);
s5 = scatter3(Outcome_OutcomeSig_data(:,1), Outcome_OutcomeSig_data(:,2), Outcome_OutcomeSig_data(:,3), ...
    70, [0 0 1], 'filled', 'MarkerEdgeColor', 'k'); % blue
%colorbar;
%view(3);
xlabel('Lateral --- Medial --- Lateral');
ylabel('Posterior --- Anterior');
zlabel('Inferior --- Superior');
% plotting actual glass brain....
patchs = patch('Faces', BrainSurfRaw.faces, 'Vertices', BrainSurfRaw.vertices, 'edgecolor', 'none', 'facecolor', 'flat', 'facealpha', .2);
facecolor = repmat([1 1 1],length(BrainSurfRaw.faces),1);
set(patchs, 'FaceVertexCData', facecolor);
camlight
lighting gouraud
legend([s1, s2, s3, s4, s5], ...
    {'All Contacts', 'PE', 'TrialColor', 'TrialType', 'Outcome'}, ...
    'Location', 'bestoutside');

% CREATE CUE-ALIGNED COLOR VARIABLES:

% Create overlapping index for variables that are encoding on the same
% contact:
% Preallocate index with zeros (no encoding)
cue_Encoding_Index = zeros(size(sig_rstd_VE_Cue_logical)); % zero == no encoding

% Convert logicals to integer bitmask components
VE_bit    = 1 * sig_rstd_VE_Cue_logical;
Color_bit = 2 * sig_rstd_trialColor_Cue_logical;
Type_bit  = 4 * sig_rstd_trialType_Cue_logical;

% Sum the bits to get unique combination codes
cue_Encoding_Index = VE_bit + Color_bit + Type_bit;

% Define RGB colors for each encoding index
cue_Encoding_ColorMap = [
    0.0, 0.0, 0.0;      % 0 = No encoding (optional, black)
    0.0, 0.6, 0.0;      % 1 = VE (Green)
    0.9, 0.0, 0.0;      % 2 = Color (Red)
    0.9, 0.6, 0.0;      % 3 = VE + Color (Yellow)
    0.0, 0.0, 1.0;      % 4 = Type (Blue)
    0.0, 0.6, 1.0;      % 5 = VE + Type (Cyan)
    0.9, 0.0, 1.0;      % 6 = Color + Type (Magenta)
    0.5, 0.5, 0.5;      % 7 = All three (Gray)
    ];


% CREATE OUTCOME-ALIGNED COLOR VARIABLES:
% Create overlapping index for variables that are encoding on the same
% contact:
% Preallocate index with zeros (no encoding)
outcome_Encoding_Index = zeros(size(sig_rstd_PE_Outcome_logical)); % zero == no encoding

% Convert logicals to integer bitmask components
PE_outcomebit    = 1 * sig_rstd_PE_Outcome_logical;
Color_outcomebit = 2 * sig_rstd_trialColor_Outcome_logical;
Type_outcomebit  = 4 * sig_rstd_trialType_Outcome_logical;
Outcome_outcomebit  = 8 * sig_rstd_Outcome_Outcome_logical;

% Sum the bits to get unique combination codes
outcome_Encoding_Index = PE_outcomebit + Color_outcomebit + Type_outcomebit + Outcome_outcomebit;

% Define RGB colors for each encoding index
outcome_Encoding_ColorMap = [
    0.0, 0.0, 0.0;      % 0 = No encoding
    0.0, 0.6, 0.0;      % 1 = PE only
    0.9, 0.0, 0.0;      % 2 = Color only
    0.9, 0.6, 0.0;      % 3 = PE + Color
    0.0, 0.0, 1.0;      % 4 = Type only
    0.0, 0.6, 1.0;      % 5 = PE + Type
    0.9, 0.0, 1.0;      % 6 = Color + Type
    0.5, 0.5, 0.5;      % 7 = PE + Color + Type
    1.0, 1.0, 0.0;      % 8 = Outcome only (Yellow)
    0.0, 1.0, 0.0;      % 9 = PE + Outcome
    1.0, 0.0, 0.5;      % 10 = Color + Outcome
    0.0, 1.0, 1.0;      % 11 = PE + Color + Outcome
    0.5, 0.0, 1.0;      % 12 = Type + Outcome
    0.4, 0.8, 0.2;      % 13 = PE + Type + Outcome
    1.0, 0.5, 0.2;      % 14 = Color + Type + Outcome
    1.0, 1.0, 1.0;      % 15 = All four variables
    ];

%% HEMISPHERE FIGURE WITH R2
%% Step 2: Plotting
figure(11); clf;
subtitle('Striatum Cue');

%% --- RIGHT HEMISPHERE ---
subplot(2,4,1)
title('Cue-Aligned Task Variable Encoding (rHemi)')
hold on
isRightStriatum = isStriatum & ElecXYZProj_All(:,1) > 0;
cueIdx_R = cue_Encoding_Index(isRightStriatum);
xyz_R = ElecXYZProj_All(isRightStriatum, :);
dot_R = dotSize_Cue(isRightStriatum);
colors_R = cue_Encoding_ColorMap(cueIdx_R + 1, :);
scatter3(xyz_R(:,1), xyz_R(:,2), xyz_R(:,3), ...
    dot_R, colors_R, 'filled', 'MarkerEdgeColor', 'k');
rightHemi = getHemisphereSurface(BrainSurfRaw, 'right');
patchs = patch('Faces', rightHemi.faces, 'Vertices', rightHemi.vertices, ...
    'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.2);
set(patchs, 'FaceVertexCData', repmat([1 1 1], length(rightHemi.faces), 1));
camlight; lighting gouraud;
view([-1 0 0]);
ylabel('Anterior --- Posterior');
zlabel('Inferior --- Superior');

%% --- LEFT HEMISPHERE ---
subplot(2,4,2)
title('Cue-Aligned Task Variable Encoding (lHemi)')
hold on
isLeftStriatum = isStriatum & ElecXYZProj_All(:,1) < 0;
cueIdx_L = cue_Encoding_Index(isLeftStriatum);
xyz_L = ElecXYZProj_All(isLeftStriatum, :);
dot_L = dotSize_Cue(isLeftStriatum);
colors_L = cue_Encoding_ColorMap(cueIdx_L + 1, :);
scatter3(xyz_L(:,1), xyz_L(:,2), xyz_L(:,3), ...
    dot_L, colors_L, 'filled', 'MarkerEdgeColor', 'k');
leftHemi = getHemisphereSurface(BrainSurfRaw, 'left');
patchs = patch('Faces', leftHemi.faces, 'Vertices', leftHemi.vertices, ...
    'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.2);
set(patchs, 'FaceVertexCData', repmat([1 1 1], length(leftHemi.faces), 1));
camlight; lighting gouraud;
view([1 0 0]);
ylabel('Posterior --- Anterior');
zlabel('Inferior --- Superior');

%% --- TOP-DOWN VIEW ---
subplot(2,4,5)
title('Cue-Aligned Encoding (Top View)')
hold on
xyz = ElecXYZProj_All(isStriatum, :);
dot = dotSize_Cue(isStriatum);
cueIdx = cue_Encoding_Index(isStriatum);
colors = cue_Encoding_ColorMap(cueIdx + 1, :);
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), ...
    dot, colors, 'filled', 'MarkerEdgeColor', 'k');
patchs = patch('Faces', BrainSurfRaw.faces, 'Vertices', BrainSurfRaw.vertices, ...
    'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.2);
set(patchs, 'FaceVertexCData', repmat([1 1 1], length(BrainSurfRaw.faces), 1));
camlight; lighting gouraud;
view([0 0 1]);
xlabel('Left --- Right');
ylabel('Posterior --- Anterior');
zlabel('Depth');

%% --- FRONT VIEW ---
subplot(2,4,6)
title('Cue-Aligned Encoding (Front View)')
hold on
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), ...
    dot, colors, 'filled', 'MarkerEdgeColor', 'k');
patchs = patch('Faces', BrainSurfRaw.faces, 'Vertices', BrainSurfRaw.vertices, ...
    'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.2);
set(patchs, 'FaceVertexCData', repmat([1 1 1], length(BrainSurfRaw.faces), 1));
camlight; lighting gouraud;
view([0 1 0]);
xlabel('Left --- Right');
ylabel('Inferior --- Superior');
zlabel('Depth');
%% R2 DOTSIZE
subplot(2,4,3)
axis off
% Define the mapping between R² ranges and dot sizes
r2_labels = {
    'R^2 ≤ 0.01', ...
    '0.01 < R^2 ≤ 0.1', ...
    '0.1 < R^2 ≤ 0.2', ...
    '0.2 < R^2 ≤ 0.3', ...
    '0.3 < R^2 ≤ 0.4', ...
    '0.4 < R^2 ≤ 0.5', ...
    '0.5 < R^2 ≤ 0.6', ...
    'R^2 > 0.6'
    };

dot_sizes = [10, 30, 40, 50, 60, 70, 80, 90];  % Match your scatter sizes
markerSizes = sqrt(dot_sizes);  % Convert to marker sizes (approximate for plot)

% Create dummy handles
dummyHandles = gobjects(length(dot_sizes),1);
hold on
for i = 1:length(dot_sizes)
    dummyHandles(i) = plot(nan, nan, 'ko', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', markerSizes(i), ...
        'DisplayName', r2_labels{i});
end
% Create the legend
lgd = legend(dummyHandles, r2_labels, 'Location', 'best');
title(lgd, 'Adj. R^2 (Cue)');

%% LEGEND FOR COLORS
% Define labels
labels = {'0: None','1: VE', '2: Color', '3: VE + Color', '4: Type', '5: VE + Type', '6: Color + Type', '7: VE + Color + Type'};
subplot(2,4,7)
hold on;
for i = 1:8
    rectangle('Position', [0, 9-i, 1, 1], ...
        'FaceColor', cue_Encoding_ColorMap(i, :), ...
        'EdgeColor', 'k');
    text(1.2, 9.5-i, labels{i}, 'FontSize', 12, 'VerticalAlignment', 'middle');
end
axis equal off;
title('Cue Encoding Color Map Legend', 'FontSize', 14);

%% STRIATUM ENCODING VARIABLES BAR:
% Data organized by variable (columns) and condition (rows)
% Each column = a variable: VE, TrialColor, TrialType
% Each row = a condition: Cue, CueSuccess
data = [STR_VE_Cue, STR_Color_Cue, STR_Type_Cue; STR_VE_CueSuccess, STR_Color_CueSuccess, STR_Type_CueSuccess];
% Plot settings
variable_labels = {'VE', 'TrialColor', 'TrialType'};
condition_labels = {'Cue', 'CueSuccess'};
bar_colors = [0.0, 0.6, 0.0; 0.9, 0.0, 0.0; 0.0, 0.0, 1.0];  % VE: green, Color: red, Type: blue
subplot(2,4,4)
bh = bar(data', 'grouped');
for i = 1:numel(bh)
    bh(i).FaceColor = 'flat';
    for j = 1:length(bh(i).YData)
        bh(i).CData(j,:) = bar_colors(j,:);
    end
    if i == 2
        bh(i).FaceAlpha = 0.3;
    else
        bh(i).FaceAlpha = 1;
    end
end
% Axes and labels
set(gca, 'XTickLabel', variable_labels);
ylabel('% Striatum Contacts');
ylim([0 100]);
title('Striatum Encoding: Cue vs CueSuccess');
% Add percentage labels on bars
for i = 1:numel(bh)
    x = bh(i).XData + bh(i).XOffset;
    y = bh(i).YData;
    for j = 1:length(y)
        text(x(j), y(j)+1.5, sprintf('%.1f%%', y(j)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
end








%% BASIC VENN DIAGRAM I CAN CHANGE IN AI :
subplot(1,2,2)
axis equal off;
hold on;
% Define circle centers and radius
r = 1.5;
center_VE    = [-1, 0];
center_Color = [1, 0];
center_Type  = [0, 1.7];
% Create circular masks using fill (approximate a Venn layout)
theta = linspace(0, 2*pi, 100);
% VE - Green
x1 = r * cos(theta) + center_VE(1);
y1 = r * sin(theta) + center_VE(2);
fill(x1, y1, [0.0, 0.6, 0.0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
% Color - Red
x2 = r * cos(theta) + center_Color(1);
y2 = r * sin(theta) + center_Color(2);
fill(x2, y2, [0.9, 0.0, 0.0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
% Type - Blue
x3 = r * cos(theta) + center_Type(1);
y3 = r * sin(theta) + center_Type(2);
fill(x3, y3, [0.0, 0.0, 1.0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
% Add text labels
%text(center_VE(1)-0.6, center_VE(2), 'VE', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0 0.6 0]);
%text(center_Color(1)+0.6, center_Color(2), 'Color', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.9 0 0]);
%text(center_Type(1), center_Type(2)+1, 'Type', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0 0 1]);
%title('Cue Encoding Venn Diagram', 'FontSize', 16);
hold off;







%% Step 2: Plotting
figure(12); clf;
subtitle('Striatum Outcome');

%% --- RIGHT HEMISPHERE ---
subplot(2,4,1)
title('Outcome-Aligned Task Variable Encoding (rHemi)')
hold on
isRightStriatum = isStriatum & ElecXYZProj_All(:,1) > 0;
outcomeIdx_R = outcome_Encoding_Index(isRightStriatum);
xyz_R = ElecXYZProj_All(isRightStriatum, :);
dot_R = dotSize_Outcome(isRightStriatum);
colors_R = outcome_Encoding_ColorMap(outcomeIdx_R + 1, :);
scatter3(xyz_R(:,1), xyz_R(:,2), xyz_R(:,3), ...
    dot_R, colors_R, 'filled', 'MarkerEdgeColor', 'k');
rightHemi = getHemisphereSurface(BrainSurfRaw, 'right');
patchs = patch('Faces', rightHemi.faces, 'Vertices', rightHemi.vertices, ...
    'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.2); % add more colors to make the figures pop out a bit...
set(patchs, 'FaceVertexCData', repmat([1 1 1], length(rightHemi.faces), 1));
camlight; lighting gouraud;
view([-1 0 0]);
ylabel('Anterior --- Posterior');
zlabel('Inferior --- Superior');

%% --- LEFT HEMISPHERE ---
subplot(2,4,2)
title('Outcome-Aligned Task Variable Encoding (lHemi)')
hold on
isLeftStriatum = isStriatum & ElecXYZProj_All(:,1) < 0;
outcomeIdx_L = outcome_Encoding_Index(isLeftStriatum);
xyz_L = ElecXYZProj_All(isLeftStriatum, :);
dot_L = dotSize_Outcome(isLeftStriatum);
colors_L = outcome_Encoding_ColorMap(outcomeIdx_L + 1, :);
scatter3(xyz_L(:,1), xyz_L(:,2), xyz_L(:,3), ...
    dot_L, colors_L, 'filled', 'MarkerEdgeColor', 'k');
leftHemi = getHemisphereSurface(BrainSurfRaw, 'left');
patchs = patch('Faces', leftHemi.faces, 'Vertices', leftHemi.vertices, ...
    'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.2);
set(patchs, 'FaceVertexCData', repmat([1 1 1], length(leftHemi.faces), 1));
camlight; lighting gouraud;
view([1 0 0]);
ylabel('Posterior --- Anterior');
zlabel('Inferior --- Superior');

%% --- TOP-DOWN VIEW ---
subplot(2,4,5)
title('Outcome-Aligned Encoding (Top View)')
hold on
xyz = ElecXYZProj_All(isStriatum, :);
dot = dotSize_Outcome(isStriatum);
outcomeIdx = outcome_Encoding_Index(isStriatum);
colors = outcome_Encoding_ColorMap(outcomeIdx + 1, :);
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), ...
    dot, colors, 'filled', 'MarkerEdgeColor', 'k');
patchs = patch('Faces', BrainSurfRaw.faces, 'Vertices', BrainSurfRaw.vertices, ...
    'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.2);
set(patchs, 'FaceVertexCData', repmat([1 1 1], length(BrainSurfRaw.faces), 1));
camlight; lighting gouraud;
view([0 0 1]);
xlabel('Left --- Right');
ylabel('Posterior --- Anterior');
zlabel('Depth');

%% --- FRONT VIEW ---
subplot(2,4,6)
title('Outcome-Aligned Encoding (Front View)')
hold on
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), ...
    dot, colors, 'filled', 'MarkerEdgeColor', 'k');
patchs = patch('Faces', BrainSurfRaw.faces, 'Vertices', BrainSurfRaw.vertices, ...
    'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.2);
set(patchs, 'FaceVertexCData', repmat([1 1 1], length(BrainSurfRaw.faces), 1));
camlight; lighting gouraud;
view([0 1 0]);
xlabel('Left --- Right');
ylabel('Inferior --- Superior');
zlabel('Depth');

%% R2 DOTSIZE
subplot(2,4,3)
axis off
% Define the mapping between R² ranges and dot sizes
r2_labels = {
    'R^2 ≤ 0.01', ...
    '0.01 < R^2 ≤ 0.1', ...
    '0.1 < R^2 ≤ 0.2', ...
    '0.2 < R^2 ≤ 0.3', ...
    '0.3 < R^2 ≤ 0.4', ...
    '0.4 < R^2 ≤ 0.5', ...
    '0.5 < R^2 ≤ 0.6', ...
    'R^2 > 0.6'};

dot_sizes = [10, 30, 40, 50, 60, 70, 80, 90];  % Match your scatter sizes
markerSizes = sqrt(dot_sizes);  % Convert to marker sizes (approximate for plot)

% Create dummy handles
dummyHandles = gobjects(length(dot_sizes),1);
hold on
for i = 1:length(dot_sizes)
    dummyHandles(i) = plot(nan, nan, 'ko', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', markerSizes(i), ...
        'DisplayName', r2_labels{i});
end
% Create the legend
lgd = legend(dummyHandles, r2_labels, 'Location', 'best');
title(lgd, 'Adj. R^2 (Outcome)');

%% LEGEND FOR COLORS
% Define labels for all 16 encoding combinations
labels = {'0: None','1: PE','2: Color','3: PE + Color','4: Type','5: PE + Type','6: Color + Type','7: PE + Color + Type',...
    '8: Outcome','9: PE + Outcome','10: Color + Outcome','11: PE + Color + Outcome','12: Type + Outcome','13: PE + Type + Outcome','14: Color + Type + Outcome','15: PE + Color + Type + Outcome'};
subplot(2,4,7)
hold on;
for i = 1:16
    rectangle('Position', [0, 17-i, 1, 1], ...
        'FaceColor', outcome_Encoding_ColorMap(i, :), ...
        'EdgeColor', 'k');
    text(1.2, 17.5-i, labels{i}, 'FontSize', 10, 'VerticalAlignment', 'middle');
end
axis equal off;
title('Outcome Encoding Color Map Legend', 'FontSize', 14);

%% STRIATUM ENCODING VARIABLES BAR:
% Data organized by variable (columns) and condition (rows)
% Each column = a variable: PE, TrialColor, TrialType
% Each row = a condition: Outcome, OutcomeSuccess
data = [STR_PE_Outcome, STR_Color_Outcome, STR_Type_Outcome, STR_Outcome_Outcome; STR_PE_OutcomeSuccess, STR_Color_OutcomeSuccess, STR_Type_OutcomeSuccess];
% Plot settings
variable_labels = {'PE', 'TrialColor', 'TrialType', 'Outcome'};
condition_labels = {'Outcome', 'OutcomeSuccess'};
bar_colors = [0.0, 0.6, 0.0; 0.9, 0.0, 0.0; 0.0, 0.0, 1.0; 1.0, 1.0, 0.0];  % PE: green, Color: red, Type: blue, Outcome: yellow
subplot(2,4,4)
bh = bar(data', 'grouped');
for i = 1:numel(bh)
    bh(i).FaceColor = 'flat';
    for j = 1:length(bh(i).YData)
        bh(i).CData(j,:) = bar_colors(j,:);
    end
    if i == 2
        bh(i).FaceAlpha = 0.3;
    else
        bh(i).FaceAlpha = 1;
    end
end
% Axes and labels
set(gca, 'XTickLabel', variable_labels);
ylabel('% Striatum Contacts');
ylim([0 100]);
title('Striatum Encoding: Outcome vs OutcomeSuccess');
% Add percentage labels on bars
for i = 1:numel(bh)
    x = bh(i).XData + bh(i).XOffset;
    y = bh(i).YData;
    for j = 1:length(y)
        text(x(j), y(j)+1.5, sprintf('%.1f%%', y(j)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
end


























figure(12)
clf
subtitle('Striatum Outcome')

%% RIGHT HEMISPHERE (TOP)
subplot(2,1,1)
title('Outcome-Aligned Task Variable Encoding (rHemi)')
hold on

% Get right hemisphere striatum data
isRightStriatum = isStriatum & ElecXYZProj_All(:,1) > 0; % X > 0 for right hemi
striatum_data_R = ElecXYZProj_All(isRightStriatum, :);
s1 = scatter3(striatum_data_R(:,1), striatum_data_R(:,2), striatum_data_R(:,3), ...
    dotSize_Outcome(isRightStriatum), [0 0 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.2); % black

% Significant VE
isSigPE_R = isRightStriatum & sig_rstd_PE_Outcome_logical;
PE_data_R = ElecXYZProj_All(isSigPE_R, :);
s2 = scatter3(PE_data_R(:,1), PE_data_R(:,2), PE_data_R(:,3), ...
    dotSize_Outcome(isSigPE_R), [1 0 0], 'filled', 'MarkerEdgeColor', 'k'); % red

% Significant Color
isSigColor_R = isRightStriatum & sig_rstd_trialColor_Outcome_logical;
Color_data_R = ElecXYZProj_All(isSigColor_R, :);
s3 = scatter3(Color_data_R(:,1), Color_data_R(:,2), Color_data_R(:,3), ...
    dotSize_Outcome(isSigColor_R), [1 1 0], 'filled', 'MarkerEdgeColor', 'k'); % yellow

% Significant Type
isSigType_R = isRightStriatum & sig_rstd_trialType_Outcome_logical;
Type_data_R = ElecXYZProj_All(isSigType_R, :);
s4 = scatter3(Type_data_R(:,1), Type_data_R(:,2), Type_data_R(:,3), ...
    dotSize_Outcome(isSigType_R), [0 1 0], 'filled', 'MarkerEdgeColor', 'k'); % green

% Significant Type
isSigOutcome_R = isRightStriatum & sig_rstd_Outcome_Outcome_logical;
Outcome_data_R = ElecXYZProj_All(isSigOutcome_R, :);
s5 = scatter3(Outcome_data_R(:,1), Outcome_data_R(:,2), Outcome_data_R(:,3), ...
    dotSize_Outcome(isSigOutcome_R), [0 0 1], 'filled', 'MarkerEdgeColor', 'k'); % blue

% Brain surface
rightHemi = getHemisphereSurface(BrainSurfRaw, 'right');
patchs = patch('Faces', rightHemi.faces, 'Vertices', rightHemi.vertices, ...
    'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.2);
facecolor = repmat([1 1 1], length(rightHemi.faces), 1);
set(patchs, 'FaceVertexCData', facecolor);
camlight; lighting gouraud;
view([-1 0 0]); % Medial view of right hemi
axis square

ylabel('Anterior --- Posterior');
zlabel('Inferior --- Superior');

legend([s1, s2, s3, s4, s5], {'All Contacts', 'VE', 'TrialColor', 'TrialType', 'Outcome'}, ...
    'Location', 'bestoutside');

%% LEFT HEMISPHERE (BOTTOM)
subplot(2,1,2)
title('Outcome-Aligned Task Variable Encoding (lHemi)')
hold on

% Get left hemisphere striatum data
isLeftStriatum = isStriatum & ElecXYZProj_All(:,1) < 0; % X < 0 for left hemi
striatum_data_L = ElecXYZProj_All(isLeftStriatum, :);
s1 = scatter3(striatum_data_L(:,1), striatum_data_L(:,2), striatum_data_L(:,3), ...
    dotSize_Outcome(isLeftStriatum), [0 0 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.2); % black

% Significant VE
isSigPE_L = isLeftStriatum & sig_rstd_PE_Outcome_logical;
PE_data_L = ElecXYZProj_All(isSigPE_L, :);
s2 = scatter3(PE_data_L(:,1), PE_data_L(:,2), PE_data_L(:,3), ...
    dotSize_Outcome(isSigPE_L), [1 0 0], 'filled', 'MarkerEdgeColor', 'k'); % red

% Significant Color
isSigColor_L = isLeftStriatum & sig_rstd_trialColor_Outcome_logical;
Color_data_L = ElecXYZProj_All(isSigColor_L, :);
s3 = scatter3(Color_data_L(:,1), Color_data_L(:,2), Color_data_L(:,3), ...
    dotSize_Outcome(isSigColor_L), [1 1 0], 'filled', 'MarkerEdgeColor', 'k'); % yellow

% Significant Type
isSigType_L = isLeftStriatum & sig_rstd_trialType_Outcome_logical;
Type_data_L = ElecXYZProj_All(isSigType_L, :);
s4 = scatter3(Type_data_L(:,1), Type_data_L(:,2), Type_data_L(:,3), ...
    dotSize_Outcome(isSigType_L), [0 1 0], 'filled', 'MarkerEdgeColor', 'k'); % green

% Significant Outcome
isSigOutcome_L = isLeftStriatum & sig_rstd_Outcome_Outcome_logical;
Outcome_data_L = ElecXYZProj_All(isSigOutcome_L, :);
s5 = scatter3(Outcome_data_L(:,1), Outcome_data_L(:,2), Outcome_data_L(:,3), ...
    dotSize_Outcome(isSigOutcome_L), [0 0 1], 'filled', 'MarkerEdgeColor', 'k'); % blue

% Brain surface
leftHemi = getHemisphereSurface(BrainSurfRaw, 'left');
patchs = patch('Faces', leftHemi.faces, 'Vertices', leftHemi.vertices, ...
    'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.2);

facecolor = repmat([1 1 1], length(leftHemi.faces), 1);
set(patchs, 'FaceVertexCData', facecolor);
camlight; lighting gouraud;
view([1 0 0]); % Medial view of left hemi (mirror of right)

ylabel('Posterior --- Anterior');
zlabel('Inferior --- Superior');

axis square

% Define the mapping between R² ranges and dot sizes
r2_labels = {
    'R^2 ≤ 0.01', ...
    '0.01 < R^2 ≤ 0.1', ...
    '0.1 < R^2 ≤ 0.2', ...
    '0.2 < R^2 ≤ 0.3', ...
    '0.3 < R^2 ≤ 0.4', ...
    '0.4 < R^2 ≤ 0.5', ...
    '0.5 < R^2 ≤ 0.6', ...
    'R^2 > 0.6'
    };
dot_sizes = [10, 30, 40, 50, 60, 70, 80, 90];

% Convert scatter sizes to marker sizes (approximate)
markerSizes = sqrt(dot_sizes);  % Adjust as needed

% Create dummy 2D points for legend with real sizes
dummyHandles = gobjects(length(dot_sizes),1);
for i = 1:length(dot_sizes)
    dummyHandles(i) = plot(nan, nan, 'ko', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', markerSizes(i), ...
        'DisplayName', r2_labels{i});
end

% Create the legend with correct dot sizes
lgd = legend(dummyHandles, r2_labels, 'Location', 'northeastoutside');
title(lgd, 'Adj. R^2 (Outcome)');



