%% Critical Component Identification

K=10000;
load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_D_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'D')
load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_A_1_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_1')
load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_A_2_SampleNum_',num2str(K),'_Balanced_New_AC.mat'],'A_2')

M=size(D,1);

for i=1:size(D,1)
    A_2(i,i)=0; 
    D(i,i)=0;
end

P_1=A_1';
P_2=A_2';

M_est=D.*(P_1-P_2);
b_est=zeros(M,1);
for i=1:M
    b_est(i,1)=D(i,:)*P_2(i,:)'; 
end    

Link_index=zeros(M,1);
for i=1:M
    Link_index(i,1)=i;
end

M_tmp=M_est;
b_tmp=repmat(b_est,[1,M]);

M_vec=sum((M_tmp),1);

M_profile=[Link_index,M_vec'];
M_vec_sort=sortrows(M_profile,2);

save('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_M_vec_sort_Balanced_New_AC_2.mat','M_vec_sort')    

