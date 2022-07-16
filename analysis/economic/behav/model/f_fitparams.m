function  ll = f_fitparams(params, Generate_params)

% We need to temporarilty change the parameter fields to the current
% parameter settings if are to use generate_a_models_data to get performance
it      = 1;
fields  = fieldnames(Generate_params.model(Generate_params.current_model));

% loop through all free parameter indices (except beta)
for field = Generate_params.model(Generate_params.current_model).this_models_free_parameters'
    
    Generate_params.model(Generate_params.current_model).(fields{field}) = params(it);
    it=it+1;
    
end

% and now assign beta too
b = params(end);

% (generate_a_models_data can do multiple subjects but here we want to fit
% one subject at a time and the number of subjects to be run is set before f_fitparams function call in a
% field of Generate_params
[num_samples ranks choiceStop_all choiceCont_all] = generate_a_models_data(Generate_params);

ll = 0;



return