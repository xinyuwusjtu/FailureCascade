function Y=fun_full_D_weighted_function_vector(Coeff,state_t_plus_1,state_t,K)

    Y=sum((state_t_plus_1-Coeff'*state_t).^2);
    Y=Y/K;

end


