function [LMS_All, var_age, vAll ] = fnc_genMonteCarloObs_UNCensusAgeDist( n_population, fs)
    % Use US Census Data for list_Age

    [t_num, t_txt, t_raw] = xlsread('UN_World_Population_2019.xlsx');

    age_data = t_num(1:end-1,:)

    length_age = length(age_data);
    list_Age = [];
    for x = 1:length_age
        list_Age = [list_Age; randi([age_data(x, 1) age_data(x, 2)], int32(age_data(x, 3) .* n_population ./ 100), 1)];
    end
    
    while(length(list_Age) < n_population)
        list_Age = [list_Age; randi([20 80], 1, 1)];
    end

    [LMS_All, var_age, vAll ] = fnc_genMonteCarloObs( n_population, list_Age, fs);

end

