function Epsilon_Logistic_Learning_Entropy_A

    Train_num=20000;

    load(['IEEE118_Capacity_1_2_Flow_Initial_state_A_',num2str(Train_num),'.mat'],'Initial_state')
    load(['IEEE118_Capacity_1_2_Flow_Final_state_A_',num2str(Train_num),'.mat'],'Final_state')
    load(['IEEE118_Capacity_1_2_Flow_Epsilon_opt_A_',num2str(Train_num),'.mat'],'Epsilon_opt')
 
%     load(['IEEE39_Capacity_0_7_Initial_state_A_',num2str(Train_num),'.mat'],'Initial_state')
%     load(['IEEE39_Capacity_0_7_Final_state_A_',num2str(Train_num),'.mat'],'Final_state')
%     load(['IEEE39_Capacity_0_7_Epsilon_opt_A_',num2str(Train_num),'.mat'],'Epsilon_opt')

%     load(['IEEE118_Capacity_1_2_Flow_Initial_state_A_',num2str(Train_num),'_Random_0_8_1_0_Prob_0_7.mat'],'Initial_state')
%     load(['IEEE118_Capacity_1_2_Flow_Final_state_A_',num2str(Train_num),'_Random_0_8_1_0_Prob_0_7.mat'],'Final_state')
%     load(['IEEE118_Capacity_1_2_Flow_Epsilon_opt_A_',num2str(Train_num),'_Random_0_8_1_0_Prob_0_7.mat'],'Epsilon_opt')
    
    M=size(Initial_state,1);
    K=size(Initial_state,2);
    Beta_record=zeros(M+1,M);
    
    for i=1:M
        i
        tic
        options=optimoptions('fmincon','Algorithm','interior-point','GradObj','on','Hessianfcn',@(Coeff,lambda)Hessian_Beta_Entropy(Coeff,lambda,Initial_state(:,1:K),zeros(M+1)));
        [y,fval,exitflag,output]=fmincon(@(Coeff)fun_beta_entropy(Initial_state(:,1:K),Final_state(i,:),Coeff),rand(M+1,1),[],[],[],[],[],[],[],options); 
        fval/K;
        Beta_record(:,i)=y;
        toc
    end
%     save(['IEEE118_Capacity_1_2_Flow_Beta_record_AC_A_',num2str(Train_num),'.mat'],'Beta_record')
    save(['IEEE39_Capacity_0_7_Initial_state_A_',num2str(Train_num),'.mat'],'Beta_record')
%     save(['IEEE118_Capacity_1_2_Flow_Beta_record_A_',num2str(Train_num),'_Random_0_8_1_0_Prob_0_7.mat'],'Beta_record')
    
end



