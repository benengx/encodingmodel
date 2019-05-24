%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Ben Engelhard, Princeton University (2019).
%
%    This program is provided free without any warranty; you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation.
%    If this code is used, please cite: B Engelhard et al. Specialized coding of sensory, motor, and cognitive variables in VTA dopamine neurons. Nature, 2019.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% make_predictor_matrix_generalcase.m
%%%
%%% Description: Make the predictor matrix to be used with the process_encoding_model function
%
% arguments: base_variables - a cell array where each term is a base variable which can be either a cell array or a vector. in both cases the division is by trials in a session. For a cell array base variable, 
%                             each term is a vector that corresponds to one trial in the session, and include the base variable's values in all timepoints of that trial. A vector base variable is necessarily a
%                             'whole-trial' variable where each term is '0' or '1' and denotes the value of that variable for all timepoints in that trial. At least one base variable needs to be a cell array.          
%            var_types      - a cell array with the same length as the 'base_variables' argument where each term denotes the type of base variable. possible types are: 
%                             'event'      : a binary base variable where values of '1' denote the occurence of the event variable.
%                             'whole-trial': a binary bas e variable that takes the same value in all timepoints of a given trial, corrsponding to wether the variable's occurence was valid in that trial. 
%                             'continuous' : a base variable that takes different values at different timepoints.
%            groupings      - a cell array that defines how the base variables are grouped into the behavioral variables used for subsequent processing in the encoding model. each term is a vector of
%                             indices (corresponding to the base_variables argument). All indices in one term of the cell arary are considered to belong to the same variable when processing the encoding model.
%                             Currently, only base variables of the same type can be grouped.
%
% outputs:   pred_allmat       - cell array corresponding to a matrix of predictors, each term corresponds to one trial where rows are timepoints and columns are the different predictors
%            pred_inds_cell    - cell array where each term has a vector of indices of the predictors that belong to a specific behavioral variable
%            grouped_var_types - cell array where each term has the type ('event','whole-trial', or 'continuous') of the corresponding behavioral variable. if no value is given, then each base variable is 
%                                considered to be a behavioral variable

function [pred_allmat,pred_inds_cell,grouped_var_types] = make_predictor_matrix_generalcase(base_variables,var_types,groupings)

if nargin<3
    groupings = mat2cell(1:length(base_variables),1,ones(1,length(base_variables)));
end

num_base_vars = length(base_variables);
numtrials = length(base_variables{1});
load('spline_basis30_int.mat'); % this results in the variable 'spline_basis' which is used for event variables. It is curently built using 7 splines and 30 timepoints. 
                                % change this if you want a different spline basis set.

% find the first cell array base variable to get from it the number of timepoints in each trial
num_points_per_trial = zeros(numtrials,1);
found_cellarraybasevar = 0;
for basevarctr = 1:num_base_vars
    if iscell(base_variables{basevarctr})
        for k=1:numtrials
            num_points_per_trial(k) = length(base_variables{basevarctr}{k});
        end
        found_cellarraybasevar = 1;
        continue
    end
end
if ~found_cellarraybasevar
    error('At least one base variable needs to be a cell array. If necessary, convert a whole trial-variable to a cell array variable')
end

for k=1:numtrials
    clear pred_curmat
    
    % process event variables
    event_vars_inds = find_non_empty_cells(strfind(var_types,'event'));
    
    for varctr = 1:length(event_vars_inds)
        clear cur_event_var
        for spctr = 1:size(spline_basis,2)
            w = conv(base_variables{event_vars_inds(varctr)}{k} ,spline_basis(:,spctr));
            cur_event_var(:,spctr) = w(1:length(base_variables{event_vars_inds(varctr)}{k}));
        end
        pred_curmat{1,event_vars_inds(varctr)} = cur_event_var;
    end

    % process whole-trial variables. they can be either vectors or cell arrays
    wholetrial_vars_inds = find_non_empty_cells(strfind(var_types,'whole-trial'));
    
    for varctr = 1:length(wholetrial_vars_inds)
        if iscell(base_variables{wholetrial_vars_inds(varctr)})
            cur_wt_var = base_variables{wholetrial_vars_inds(varctr)}{k};
        else
            cur_wt_var  = zeros(num_points_per_trial(k),1) + base_variables{wholetrial_vars_inds(varctr)}(k);
        end
        pred_curmat{1,wholetrial_vars_inds(varctr)} = cur_wt_var;
    end

    % process continuous variables. they can be either vectors or cell arrays
    cont_vars_inds = find_non_empty_cells(strfind(var_types,'continuous'));
    
    for varctr = 1:length(cont_vars_inds)
        pred_curmat{1,cont_vars_inds(varctr)} = base_variables{cont_vars_inds(varctr)}{k};
    end

    pred_allmat{k,1} = cell2mat(pred_curmat);
    
    % get the indices of the predictors coresponding to the different variables
    if k==1
        for varctr = 1:size(pred_curmat,2)
           preds_inds_basevars{varctr} = (1:size(pred_curmat{1,varctr},2)) + size(cell2mat(pred_curmat(1,1:varctr-1)),2);
        end
        for grpctr = 1:length(groupings)
            pred_inds_cell{grpctr} = cell2mat(preds_inds_basevars(groupings{grpctr}));
            cur_types = var_types(groupings{grpctr});
            grouped_var_types{grpctr} = cur_types{1};
            
            if ~isempty(setdiff(cur_types,cur_types{1}))
               error('Grouped variables must have the same type') 
            end
        end
    end
    
end





