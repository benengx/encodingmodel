% drop_type: 'norefit' - calculate regression weights with the full model, then zero the weights correspoding to the predictors being dropped
% drop_type: 'refit'   - calculate regression weights without the weights correspoding to the predictors being dropped (partial model)

function [R2,all_predicted,B_all] = get_CV_R2(full_predmat_cell,cur_neural_act_mat,test_trials_folds,train_trials_folds,trial_length_vec,inds_to_drop,drop_type,trial_types_to_match)

if nargin<6
    inds_to_drop=[];
end
if nargin<7
    drop_type='norefit';
end
if nargin<8
    trial_types_to_match = [];
end
if ~isempty(trial_types_to_match)
    for trctr=1:length(cur_neural_act_mat)
        tr_types_cell{trctr,1} = ones(size(cur_neural_act_mat{trctr},1),1)*trial_types_to_match(trctr);
    end
end

all_predicted = cell(size(cur_neural_act_mat));
all_neural_act = cell(size(cur_neural_act_mat));
for k=1:length(test_trials_folds)
    
    cur_Xtrain = cell2mat(full_predmat_cell(train_trials_folds{k}));
    cur_Ytrain = cell2mat(cur_neural_act_mat(train_trials_folds{k}));
    cur_Xtest  = cell2mat(full_predmat_cell(test_trials_folds{k}));
    cur_Ytest  = cell2mat(cur_neural_act_mat(test_trials_folds{k}));
    
    if ~isempty(inds_to_drop)
        switch drop_type
            case 'norefit'
                if isempty(trial_types_to_match)
                    curB = glmfit(cur_Xtrain,cur_Ytrain,'normal','constant','on');
                else
                    
                    cur_tr_types = cell2mat(tr_types_cell(train_trials_folds{k}));
                    all_tr_types = unique(cur_tr_types);
                    tr_types_num = zeros(1,length(all_tr_types));
                    for l=1:length(all_tr_types)
                        tr_types_num(l) = sum(cur_tr_types==all_tr_types(l));
                    end
                    weights_types = prod(tr_types_num)./tr_types_num/sum(tr_types_num);
                    weights_timepoints = zeros(size(cur_tr_types));
                    for l=1:length(all_tr_types)
                        weights_timepoints(cur_tr_types==all_tr_types(l)) = weights_types(l);
                    end
                    
                    curB = glmfit(cur_Xtrain,cur_Ytrain,'normal','constant','on','weights',weights_timepoints);
                end
                curB(inds_to_drop+1)=0;
                cur_Ypred = [ones(size(cur_Xtest,1),1) cur_Xtest]*curB;
            case 'refit'
                inds_to_use = setdiff(1:size(cur_Xtrain,2),inds_to_drop);
                if isempty(trial_types_to_match)
                    curB = glmfit(cur_Xtrain(:,inds_to_use ),cur_Ytrain,'normal','constant','on');
                else
                    
                    cur_tr_types = cell2mat(tr_types_cell(train_trials_folds{k}));
                    all_tr_types = unique(cur_tr_types);
                    tr_types_num = zeros(1,length(all_tr_types));
                    for l=1:length(all_tr_types)
                        tr_types_num(l) = sum(cur_tr_types==all_tr_types(l));
                    end
                    weights_types = prod(tr_types_num)./tr_types_num/sum(tr_types_num);
                    weights_timepoints = zeros(size(cur_tr_types));
                    for l=1:length(all_tr_types)
                        weights_timepoints(cur_tr_types==all_tr_types(l)) = weights_types(l);
                    end
                    
                    curB = glmfit(cur_Xtrain(:,inds_to_use ),cur_Ytrain,'normal','constant','on','weights',weights_timepoints);
                end
                cur_Ypred = [ones(size(cur_Xtest,1),1) cur_Xtest(:,inds_to_use)]*curB;
            otherwise
                error('unknown drop type')
        end
    else
        if isempty(trial_types_to_match)
            
            curB = glmfit(cur_Xtrain,cur_Ytrain,'normal','constant','on');
        else
            
            cur_tr_types = cell2mat(tr_types_cell(train_trials_folds{k}));
            all_tr_types = unique(cur_tr_types);
            tr_types_num = zeros(1,length(all_tr_types));
            for l=1:length(all_tr_types)
                tr_types_num(l) = sum(cur_tr_types==all_tr_types(l));
            end
            weights_types = prod(tr_types_num)./tr_types_num/sum(tr_types_num);
            weights_timepoints = zeros(size(cur_tr_types));
            for l=1:length(all_tr_types)
                weights_timepoints(cur_tr_types==all_tr_types(l)) = weights_types(l);
            end
            curB = glmfit(cur_Xtrain,cur_Ytrain,'normal','constant','on','weights',weights_timepoints);
            
        end
        cur_Ypred = [ones(size(cur_Xtest,1),1) cur_Xtest]*curB;
    end
    all_predicted(test_trials_folds{k},1) = mat2cell(cur_Ypred ,trial_length_vec(test_trials_folds{k}),1);
    all_neural_act(test_trials_folds{k},1) = mat2cell(cur_Ytest ,trial_length_vec(test_trials_folds{k}),1);
end
R2 = corr(cell2mat(all_neural_act),cell2mat(all_predicted)).^2;


% get weights for regression on all data
if nargout>2
    cur_Xall = cell2mat(full_predmat_cell);
    cur_Yall = cell2mat(cur_neural_act_mat);
    
    if isempty(trial_types_to_match)
        B_all = glmfit(cur_Xall ,cur_Yall ,'normal','constant','on');
    else
        
        cur_tr_types = cell2mat(tr_types_cell);
        all_tr_types = unique(cur_tr_types);
        tr_types_num = zeros(1,length(all_tr_types));
        for l=1:length(all_tr_types)
            tr_types_num(l) = sum(cur_tr_types==all_tr_types(l));
        end
        weights_types = prod(tr_types_num)./tr_types_num/sum(tr_types_num);
        weights_timepoints = zeros(size(cur_tr_types));
        for l=1:length(all_tr_types)
            weights_timepoints(cur_tr_types==all_tr_types(l)) = weights_types(l);
        end
        
        B_all = glmfit(cur_Xall ,cur_Yall,'normal','constant','on','weights',weights_timepoints);
    end
    
end

