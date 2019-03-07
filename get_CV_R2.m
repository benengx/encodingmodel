% drop_type: 'norefit' - calculate regression weights with the full model, then zero the weights correspoding to the predictors being dropped
% drop_type: 'refit'   - calculate regression weights without the weights correspoding to the predictors being dropped (partial model)

function R2 = get_CV_R2(full_predmat_cell,cur_neural_act_mat,test_trials_folds,train_trials_folds,trial_length_vec,inds_to_drop,drop_type)

if nargin<6
    inds_to_drop=[];
end
if nargin<7
    drop_type='norefit';
end

all_predicted = cell(size(cur_neural_act_mat));
for k=1:length(test_trials_folds)
    
    cur_Xtrain = cell2mat(full_predmat_cell(train_trials_folds{k}));
    cur_Ytrain = cell2mat(cur_neural_act_mat(train_trials_folds{k}));
    cur_Xtest  = cell2mat(full_predmat_cell(test_trials_folds{k}));
    
    if ~isempty(inds_to_drop)
        switch drop_type
            case 'norefit'
                curB = glmfit(cur_Xtrain,cur_Ytrain,'normal','constant','on');
                curB(inds_to_drop+1)=0;
                cur_Ypred = [ones(size(cur_Xtest,1),1) cur_Xtest]*curB;
            case 'refit'
                inds_to_use = setdiff(1:size(cur_Xtrain,2),inds_to_drop);
                curB = glmfit(cur_Xtrain(:,inds_to_use ),cur_Ytrain,'normal','constant','on');                
                cur_Ypred = [ones(size(cur_Xtest,1),1) cur_Xtest(:,inds_to_use)]*curB;
            otherwise
                error('unknown drop type')
        end
    else
        curB = glmfit(cur_Xtrain,cur_Ytrain,'normal','constant','on');
        cur_Ypred = [ones(size(cur_Xtest,1),1) cur_Xtest]*curB;
    end
    all_predicted(test_trials_folds{k},1) = mat2cell(cur_Ypred ,trial_length_vec(test_trials_folds{k}),1);
end

R2 = corr(cell2mat(cur_neural_act_mat),cell2mat(all_predicted)).^2;







