%% Construct the threshold pool from the training failure cascade samples

Test_num=5000;

load(['IEEE2383_Capacity_1_Flow_1_InitFail_2_D_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'D')
load(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_1_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'A_1')
load(['IEEE2383_Capacity_1_Flow_1_InitFail_2_A_2_SampleNum_',num2str(Test_num),'_Balanced_New.mat'],'A_2')

for i=1:size(A_2,1)
    A_2(i,i)=0; 
end

P_1=A_1';
P_2=A_2';
M=size(D,1);

%Calculate M_est,b_est
M_est=D.*(P_1-P_2);
sum_tmp_1=zeros(M,1);
M_1=D.*P_1;
for i=1:M
    sum_tmp_1(i,1)=sum(M_1(i,:)); 
end
sum_tmp_2=zeros(M,1);
M_2=D.*P_2;
for i=1:M
    sum_tmp_2(i,1)=sum(M_2(i,:)); 
end
b_est=zeros(M,1);
for i=1:M
    b_est(i,1)=D(i,:)*P_2(i,:)'; 
end

load('IEEE2383_Capacity_1_Flow_1_InitFail_2_Train_Balanced_New.mat','cascade_train')
cascade=cascade_train;

Test_num=20000;
K=Test_num;

size_state_vector=zeros(K,1);
weight_vector=zeros(K,1);
for i=1:K
    max_i=max(cascade(:,i));
    tmp_cascade=cascade(:,i);
    tmp_cascade(find(cascade(:,i)==max_i))=0;
    submax_i=max(tmp_cascade);
    if max_i-submax_i>=2
        size_state=submax_i+1; 
    else
        size_state=submax_i;
    end
    size_state_vector(i,1)=size_state;
    weight_vector(i,1)=1;
end

Epsilon_opt=zeros(M,Test_num);
Initial_state=zeros(M,Test_num);
Final_state=zeros(M,Test_num);
Flag_opt=zeros(M,Test_num);

for k=1:Test_num
    real_vector=cascade(:,k);

    real_final_state=zeros(M,1);
    max_real=max(real_vector(:,1));
    tmp_real=real_vector(:,1);
    tmp_real(find(real_vector(:,1)==max_real))=0;
    submax_real=max(tmp_real);
    if max_real-submax_real>=2
        real_final_state(find(real_vector(:,1)==max_real),1)=1;
    end

    Final_state(:,k)=real_final_state;

    %%Find the corresponding best epsilon
    state=ones(M,1);
    tmp_rnd=find(real_vector==1);
    state(tmp_rnd)=0;
    Initial_state(:,k)=state;
    state_1_cum=[];
    state_1_cum=[state_1_cum,state];

    count=0;
    while count<size_state_vector(k,1)
        count=count+1;
        state_1=(M_est*state+b_est);
        state=state_1;
        state_1_cum=[state_1_cum,state_1];
    end

    %Determine the threshold
    for i=1:M
        if real_vector(i)>1 && real_final_state(i)==0 
            Epsilon_opt(i,k)=0.5*state_1_cum(i,real_vector(i))+0.5*state_1_cum(i,real_vector(i)-1);  %Here?
        end
        if real_vector(i)==1
            Epsilon_opt(i,k)=1; 
        end
        if real_final_state(i)==1
            Epsilon_opt(i,k)=state_1_cum(i,real_vector(i))*0.8;
            Flag_opt(i,k)=1;
        end
    end
end

save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_SampleNum_',num2str(Test_num),'_Initial_state_Balanced.mat'],'Initial_state')
save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_SampleNum_',num2str(Test_num),'_Final_state_Balanced.mat'],'Final_state')
save(['IEEE2383_Capacity_1_Flow_1_InitFail_2_SampleNum_',num2str(Test_num),'_Epsilon_opt_Balanced.mat'],'Epsilon_opt')


