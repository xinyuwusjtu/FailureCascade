%% Output the basic information of the tested power systems and the corresponding failure cascade statistics

mpc=case118; %case118;case2383wp;case1888rte;case1354pegase

%     load('IEEE2383_Capacity_1_InitFail_2_New_Train_Tmp_2.mat','cascade_train')
%     cascade_2=cascade_train;    
%     load('IEEE2383_Capacity_1_InitFail_2_New_Test_Tmp_2.mat','cascade_test')
%     cascade_3=cascade_test;
%     cascade_train=[cascade_2,cascade_3];

load('IEEE118_Capacity_1_2_InitFail_3_Train_Balanced_New_AC.mat','cascade_train')

cascade=cascade_train;
Train_num=50000;
M=size(cascade_train,1);

fail_size_record=zeros(Train_num,1);

for k=1:Train_num
    k

    real_vector=cascade(:,k);
    tmp_rnd=cascade(:,k)==1;

    real_final_state=zeros(M,1);
    max_real=max(real_vector(:,1));
    tmp_real=real_vector(:,1);
    tmp_real(real_vector(:,1)==max_real)=0;
    submax_real=max(tmp_real);
    if max_real-submax_real>=2
        real_final_state(real_vector(:,1)==max_real,1)=1;
    end

    fail_size_record(k,1)=size(find(real_final_state==0),1);       

end

size_state_vector=zeros(Train_num,1);
for i=1:Train_num
    max_i=max(cascade(:,i));
    tmp_cascade=cascade(:,i);
    tmp_cascade(cascade(:,i)==max_i)=0;
    submax_i=max(tmp_cascade);
    if max_i-submax_i>=2
        size_state=submax_i+1; 
    else
        size_state=submax_i;
    end
    size_state_vector(i,1)=size_state;
end

fail_size_record_avg=sum(fail_size_record)/Train_num
fail_size_record_max=max(fail_size_record)
fail_size_record_min=min(fail_size_record)

figure('visible','on')
[counts,centers]=hist(fail_size_record,20);
fig_tmp=bar(centers,counts/sum(counts));

num_generator=size(mpc.gen,1);
num_link=M;

failure_exist_flag=zeros(M,1);
for i=1:Train_num 
    tmp_failure=zeros(M,1);
    tmp_failure(cascade(:,i)<size_state_vector(i,1) & cascade(:,i)>1,1)=1;
    failure_exist_flag=failure_exist_flag+tmp_failure;
end

num_eff_link=size(find(failure_exist_flag>0),1);
eff_rate=num_eff_link/num_link;


num_generator
num_link
num_eff_link
eff_rate
fail_size_record_avg
fail_size_record_max
fail_size_record_min


