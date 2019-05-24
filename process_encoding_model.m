%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Ben Engelhard, Princeton University (2019).
%
%    This program is provided free without any warranty; you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation.
%    If this code is used, please cite: B Engelhard et al. Specialized coding of sensory, motor, and cognitive variables in VTA dopamine neurons. Nature, 2019.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% process_encoding_model.m
%%%
%%% Description: Process the encoding model on a matrix of predictors belonging to different behavioral variables and a matrix of neural activity traces correspoding to one or more neurons. return the relative
%%% contribution of each behavioral variable to each neuron and the F-statistic derived from a nested model comparison for each behavioral variable and neuron to determine if the behavioral variable is 
%%% significantly represented in the neuron's activity.
%
% arguments: pred_allmat      - cell array correspoding to a matrix of behavioral predictors, each term is a trial and contains a matrix where rows are timepoints and columns are behavioral predictors.
%            pred_inds_cell   - cell array where each term has a vector of indices of the predictors that belong to a specific behavioral variable
%            neural_act_mat   - cell array correspoding to a matrix of activity traces, each term is a trial and contains a matrix where rows are timepoints and columns are traces correspoding to different
%                               neurons. Timepoints where the neuronal activity is not defined (e.g. a the imaging became unstable) are filled with NaNs (Currently it is assumed that before the first NaN
%                               activity was always defined, and after the first NaN activity is never defined).
%            pred_types_cell  - cell array where each terms indicates the type of behavioral variable ('event', 'whole-trial', or 'continuous').
%            approach         - 'norefit': calculate regression weights with the full model, then zero the weights correspoding to the predictors being dropped. 'refit' : calculate regression weights
%                               without the weights correspoding to the predictors being dropped (partial model).
%
% outputs:   relative_contrib - a matrix where rows are behavioral variables and columns are neurons. each term is the relative contribution of the behavioral variable to the neural activity.
%            Fstat_mat        - a matrix where rows are behavioral variables and columns are neurons. each term is the F-statistic associated with the nested model comparison where the predicted variable is the
%                               activity of the neuron, the full model is the model containing the predictors from al behavioral variables, and the partial model is the model containing the predictors from all
%                               variables except the one being tested. The value of this statistic shoudl be compared to a distirbution of staitsitc values obtained by erforing the same oprateion on shuffled 
%                               data, directly using the p-value assocaited with the staitstic is not valid given the autocorrelations in the data.

function [relative_contrib,Fstat_mat] = process_encoding_model(pred_allmat, pred_inds_cell, neural_act_mat, pred_types_cell,approach)

numcells = size(neural_act_mat{1},2);
numtrials_all = length(pred_allmat);

% find for each neuron trials where activity is defined, and also the length of each trial
defined_mat = zeros(numtrials_all,numcells);
trial_length_vec = zeros(numtrials_all,1);

for trctr=1:numtrials_all
    defined_mat(trctr,:) = ~sum(isnan(neural_act_mat{trctr}));
    trial_length_vec (trctr) = size(neural_act_mat{trctr},1);
end


% use crossvalidation to find the best polynomial degree to apply for the continuous variables
num_cv_folds = 5;
max_poly_deg = 3;

% rng(0,'twister')
full_R2_vec = zeros(numcells,1);
partial_R2_vec = zeros(numcells,length(pred_inds_cell));
relative_contrib = zeros(numcells,length(pred_inds_cell));
for cellctr = 1:numcells
    cur_good_trials = 1:find(defined_mat(:,cellctr),1,'last');
    num_trials_per_fold = ceil(length(cur_good_trials)/num_cv_folds);
    
    temp_neural_act = cell2mat(neural_act_mat(cur_good_trials));
    cur_neural_act_mat = mat2cell(temp_neural_act(:, cellctr),trial_length_vec(cur_good_trials),1);
    
    % zscore the predictors
    cur_pred_allmat_z = mat2cell(zscore(cell2mat(pred_allmat(cur_good_trials))),trial_length_vec(cur_good_trials),size(pred_allmat{1},2));
    rng('default')
    cur_random_vector = randperm(length(cur_good_trials));
    
    % get indices of test and train trials for CV
    for foldctr = 1:num_cv_folds
        kf_inds{foldctr} = num_trials_per_fold*(foldctr-1)+1:min(num_trials_per_fold*foldctr,length(cur_good_trials));
        test_trials_folds{foldctr} = cur_random_vector(kf_inds{foldctr});
        train_trials_folds{foldctr} = setdiff(cur_random_vector,test_trials_folds{foldctr});
    end
    
    [~,F_vec] = get_f_pvals_reg(cell2mat(cur_pred_allmat_z),zscore(cell2mat(cur_neural_act_mat)),pred_inds_cell);
    Fstat_mat(cellctr,:) = F_vec;
    
    % make matrix of all possible combinations of polynomial degrees for all continuous variables
    cont_inds = find_non_empty_cells(strfind(pred_types_cell,'continuous'));
    non_cont_inds = setdiff(1:length(pred_types_cell),cont_inds);
    num_cont_inds = length(cont_inds);
    all_degs_mat = [];
    clear cur_cont_preds_inds cur_base_preds
    temp_predmat = cell2mat(cur_pred_allmat_z);
    for cont_var_ctr = 1:num_cont_inds
        all_degs_mat  = ceil([all_degs_mat mod((1:max_poly_deg^num_cont_inds)/(max_poly_deg^(cont_var_ctr-1)),max_poly_deg+.01)']);
        cur_cont_preds_inds{cont_var_ctr} = pred_inds_cell{cont_inds(cont_var_ctr)};
        cur_base_preds{cont_var_ctr} = temp_predmat(:,cur_cont_preds_inds{cont_var_ctr}); % predictors for the current continuous variable before adding any additional polynomial degrees
    end
    all_cont_pred_inds = cell2mat(cur_cont_preds_inds);
    all_non_cont_pred_inds = setdiff(1:size(cur_pred_allmat_z{1},2),all_cont_pred_inds);
    non_cont_predmat = temp_predmat(:,all_non_cont_pred_inds);
    
    
    deg_R2_vec = zeros(1,size(all_degs_mat,1));
    for degctr = 1:size(all_degs_mat,1)
        full_predmat = non_cont_predmat;
        for cont_var_ctr = 1:num_cont_inds
            cur_cont_pred_add = [];
            for curdegctr = 1:all_degs_mat(degctr,cont_var_ctr)
                cur_cont_pred_add = [cur_cont_pred_add cur_base_preds{cont_var_ctr}.^curdegctr];
            end
            full_predmat = [full_predmat cur_cont_pred_add];
        end
        
        %get CV R2 for this predictor matrix
        full_predmat_cell = mat2cell(full_predmat,trial_length_vec(cur_good_trials),size(full_predmat,2));
        deg_R2_vec(degctr) = get_CV_R2(full_predmat_cell,cur_neural_act_mat,test_trials_folds,train_trials_folds,trial_length_vec(cur_good_trials));
    end
    
    [cur_full_R2,best_deg_ind] = max(deg_R2_vec);
    full_R2_vec(cellctr,1) = cur_full_R2;
    
    % now make the matrix with the optimal poly degree for each continuous variable and update the predictor indices as well
    clear pred_inds_cell_new
    full_predmat = [];
    pred_inds_cell_new = {};
    new_pred_inds_ctr = 1;
    for non_cont_var_ctr = 1:length(non_cont_inds)
        cur_predmatvar = temp_predmat(:,cell2mat(pred_inds_cell(non_cont_inds(non_cont_var_ctr))));
        cur_num_preds = size(full_predmat ,2);
        full_predmat = [full_predmat cur_predmatvar];
        pred_inds_cell_new{non_cont_inds(non_cont_var_ctr)} = (1:size(cur_predmatvar,2))+cur_num_preds ;
    end
    for cont_var_ctr = 1:num_cont_inds
        all_cont_pred_add{cont_var_ctr} = [];
        for curdegctr = 1:all_degs_mat(best_deg_ind,cont_var_ctr)
            all_cont_pred_add{cont_var_ctr} = [all_cont_pred_add{cont_var_ctr} cur_base_preds{cont_var_ctr}.^curdegctr];
        end
        cur_num_preds = size(full_predmat ,2);
        full_predmat = [full_predmat zscore(all_cont_pred_add{cont_var_ctr})];
        pred_inds_cell_new{cont_inds(cont_var_ctr)} = (1:size(all_cont_pred_add{cont_var_ctr},2))+cur_num_preds ;
    end
    
    full_predmat_cell = mat2cell(full_predmat,trial_length_vec(cur_good_trials),size(full_predmat,2));
    
    
    % now calculate the relative contributions.  first calculate the R2 when each variable is omitted.
    
    if ~isempty(full_R2_vec(cellctr,1))
        full_R2_vec(cellctr,1) = get_CV_R2(full_predmat_cell,cur_neural_act_mat,test_trials_folds,train_trials_folds,trial_length_vec(cur_good_trials));
    end
    for varctr = 1:length(pred_inds_cell_new)
        partial_R2_vec(cellctr,varctr) =  get_CV_R2(full_predmat_cell,cur_neural_act_mat,test_trials_folds,train_trials_folds,trial_length_vec(cur_good_trials),cell2mat(pred_inds_cell_new(varctr)),approach);
    end
    
    cur_R2_diff = (full_R2_vec(cellctr,1) - partial_R2_vec(cellctr,:))/full_R2_vec(cellctr,1);
    cur_R2_diff(cur_R2_diff<0)=0;
    cur_R2_diff(cur_R2_diff>1)=1;
    relative_contrib(cellctr,:) = cur_R2_diff/sum(cur_R2_diff);
    
    
end



