function Y=fun_full_D_weighted_function_vector(Coeff,state_t_plus_1,state_t,K)

%     tmp=state_t_plus_1-Coeff'*state_t;
    Y=sum((state_t_plus_1-Coeff'*state_t).^2);
    Coeff'*state_t;
    Y=Y/K;

%     Y=y/(size(state_t_plus_1,2));

end


