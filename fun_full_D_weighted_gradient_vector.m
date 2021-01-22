function Grad=fun_full_D_weighted_gradient_vector(Coeff,state_t_plus_1_prod,state_t_prod,K)

    Grad=2*(state_t_prod*Coeff-state_t_plus_1_prod);
    
    Grad=Grad/K;

end



