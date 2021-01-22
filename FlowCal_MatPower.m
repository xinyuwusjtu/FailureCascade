
%%%%%%%%%%%%%%%%%   pctRunOnAll warning off

    %load('italyNet.mat' //ItalianGrid
    mpc = case1354pegase; %Switch different size of grids //case1888rte, case1354pegase, case2383wp, case2736sp, case2848rte, case2869pegase, case3012wp
    mpc.branch=mpc.branch(:,1:13);
    mpc.bus=mpc.bus(:,1:13);
    mpc.gen=mpc.gen(:,1:21);
    
    beta=1;
    mpc.bus(:,3)=beta*mpc.bus(:,3);
    mpc.bus(:,4)=beta*mpc.bus(:,4);
    mpc.gen(:,2)=beta*mpc.gen(:,2);
    mpc.gen(:,3)=beta*mpc.gen(:,3);
    
    bus = mpc.bus;

    %IEEE1888
    Node_Index_List=zeros(size(bus,1),2);
    for i=1:size(bus,1)
        Node_Index_List(i,1)=i;
        Node_Index_List(i,2)=bus(i,1);
        bus(i,1)=i;
    end
    bus=sortrows(bus,1);
    mpc.bus = bus;
    %//
    gen = mpc.gen;
    for i=1:size(gen,1)
        gen(i,1)=Node_Index_List(find(gen(i,1)==Node_Index_List(:,2)),1); 
    end
    gen=sortrows(gen,1);
    mpc.gen=gen;
    %//
    branch_tmp=mpc.branch;
%     find(branch_tmp(i,1)==Node_Index_List(:,2))
    for i=1:size(branch_tmp,1)
        branch_tmp(i,1)=Node_Index_List(find(branch_tmp(i,1)==Node_Index_List(:,2)),1);
        branch_tmp(i,2)=Node_Index_List(find(branch_tmp(i,2)==Node_Index_List(:,2)),1);
    end
    mpc.branch=branch_tmp;
    
    %Delete the multi-links (Combining the capacity together)
    for i=1:size(mpc.branch,1)
        for j=i+1:size(mpc.branch,1)
            if (mpc.branch(i,1)==mpc.branch(j,1) && mpc.branch(i,2)==mpc.branch(j,2)) || (mpc.branch(i,1)==mpc.branch(j,2) && mpc.branch(i,2)==mpc.branch(j,1))
                mpc.branch(i,6)=mpc.branch(i,6)+mpc.branch(j,6);
                mpc.branch(j,:)=0; 
            end
        end
    end
    
    mpc.branch(find(mpc.branch(:,1)==0),:)=[];
        
    branch = mpc.branch;
    reactance_vector= branch(:,4);
    
    powerInj_consume=bus(:,3);
    powerInj_supply=zeros(size(powerInj_consume,1),1);

    for i=1:size(gen,1)
        powerInj_supply(find(bus(:,1)==gen(i,1)),1)=powerInj_supply(find(bus(:,1)==gen(i,1)),1)+gen(i,2);
    end
    
%     powerInj_supply(100)=powerInj_supply(100)-135.4;
%     
%     index=find(gen(:,1)==100);
%     mpc.gen(index,2)=gen(index,2)-135.4;

    %IEEE57
%     powerInj_supply(1)=powerInj_supply(1)+80.475;
%     powerInj_supply(3)=powerInj_supply(3)+80.475;
%     powerInj_supply(8)=powerInj_supply(8)+80.475;
%     powerInj_supply(12)=powerInj_supply(12)+80.475;
    
%     index=find(gen(:,1)==1);
%     mpc.gen(index,2)=gen(index,2)+80.475;
%     index=find(gen(:,1)==3);
%     mpc.gen(index,2)=gen(index,2)+80.475;
%     index=find(gen(:,1)==8);
%     mpc.gen(index,2)=gen(index,2)+80.475;
%     index=find(gen(:,1)==12);
%     mpc.gen(index,2)=gen(index,2)+80.475;

    powerInj_vector=powerInj_supply-powerInj_consume;
    
    %%%%%%%%%%%%%%%DC should be balanced
    powerInj_consume=powerInj_consume+sum(powerInj_vector)/size(powerInj_vector,1);
    powerInj_vector=powerInj_supply-powerInj_consume;
    bus(:,3)=powerInj_consume;
    mpc.bus=bus;
    
    sum(powerInj_vector)
    
    size(find(mpc.branch(:,6)==0),1)/size(mpc.branch,1)
    
    Adj=zeros(size(bus,1));
    Reactance=zeros(size(bus,1));
    Capacity=zeros(size(bus,1));
    for i=1:size(branch,1)
        Adj(branch(i,1),branch(i,2))=1;
        Adj(branch(i,2),branch(i,1))=1;
        Reactance(branch(i,1),branch(i,2))=branch(i,4);
        Reactance(branch(i,2),branch(i,1))=branch(i,4);
        if branch(i,6)>0
            Capacity(branch(i,1),branch(i,2))=branch(i,6);
            Capacity(branch(i,2),branch(i,1))=branch(i,6);
        else
            Capacity(branch(i,1),branch(i,2))=9900;
            Capacity(branch(i,2),branch(i,1))=9900;
        end
    end
    Capacity_vector=branch(:,6); %IEEE30/IEEE39
    count=0;
    for i=1:size(Capacity_vector,1)
        if Capacity_vector(i,1)==0
            count=count+1;
            Capacity_vector(i,1)=9900;
        end
    end
    branch(:,6)=Capacity_vector;
    mpc.branch=branch;
        
    powerInjected_ori=powerInj_vector';
    adjacency_ori=sparse(Adj);
    Adjacency=adjacency_ori;
    reactance=sparse(Reactance);
    
    M=size(find(triu(full(adjacency_ori))),1);
    
    link_row=branch(:,1);
    link_column=branch(:,2);
    
    z=[link_row,link_column];
    RowColumnIndex=z;
    
    %save('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Adjacency.mat','Adjacency')
    filename = 'IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Adjacency.xlsx';
    xlswrite(filename,Adj);
    %save('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_RowColumnIndex.mat','RowColumnIndex')
    filename = 'IEEE1354_Capacity_1_Flow_1_5_InitFail_2_RowColumnIndex.xlsx';
    xlswrite(filename,RowColumnIndex);
    load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Train_Balanced_New_AC.mat'],'cascade_link_cum')
    cascade_train=cascade_link_cum(:,1:200);
    filename = 'IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Train_Balanced_New_AC.xlsx';
    xlswrite(filename,cascade_train);    
    index_vector=sort(randperm(M),'ascend')';
    mpc.branch=[mpc.branch,index_vector];
    
    %%%
%     mpc.branch(2084,:);

%     mpc.bus(:,3)=0;
%     mpc.bus(:,4)=0;
%     mpc.gen(:,2)=0;
%     mpc.gen(:,3)=0;
    
%     result=rundcpf_me(mpc);
%     result.branch=[result.branch,mpc.branch(:,14)];
%     tempur_vector=abs(result.branch(:,14)); %%%%DC
    
    result=rundcpf_me(mpc);
    result.branch=[result.branch,mpc.branch(:,14)];
    tempur_vector=1/2*(sqrt(result.branch(:,14).^2+result.branch(:,15).^2)+sqrt(result.branch(:,16).^2+result.branch(:,17).^2));  %%%%AC
    
    Delta=[tempur_vector,Capacity_vector];
    Link_index_tmp=find(tempur_vector-Capacity_vector>=0);
    [Link_index_tmp,Delta(Link_index_tmp,:)];
    
%     P_flow=zeros(size(mpc.bus,1));
%     for i=1:M
%         P_flow(result.branch(i,1),result.branch(i,2))=tempur_vector(i,1);
%         P_flow(result.branch(i,2),result.branch(i,1))=tempur_vector(i,1);
%     end
% 
%     Capacity=1.5*abs(P_flow);
%     %Capacity-abs(P_flow)
%     Capacity_vector=zeros(M,1); %IEEE118
%     for i=1:M
%         Capacity_vector(i,1)=Capacity(link_row(i),link_column(i)); 
%     end

    alpha=1;
    K=5000;
    link_set=[link_row,link_column];
    
    for k=2:2
    cascade_link_cum=zeros(M,K);
    k;
    S_vec=zeros(K,1);
    tic
    for i=1:K
        i
%         tic
        rnd=randperm(M);
        tmp_index=(rnd(1:k))';
%         tmp_index=i;
%         tmp_index=[79;827];
        initFail=link_set(tmp_index,:);
        initFail_Index=tmp_index;
%         tic
        [yield, adjacency, history, powerInjected, S] = fCascadeNew_Matpower_Revised_Slack(mpc, powerInjected_ori, adjacency_ori, reactance, Capacity, Capacity_vector, alpha, initFail,tmp_index);
%         toc
        cascade=full(history);
        cascade_link=zeros(M,1);
        for j=1:M
            cascade_link(j,1)=cascade(link_row(j),link_column(j));
        end
        cascade_link_cum(:,i)=cascade_link;
        
        1;
        
%         size(find(cascade_link==2),1);
        S_vec(i,1)=S;
    end
    toc
    
%       save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_' num2str(k) '_Balanced_New.mat'], 'cascade_link_cum')
%     save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_' num2str(k) '_Balanced_New_Retry.mat'], 'cascade_link_cum')
%     save(['IEEE2383_Capacity_1_5_Flow_1_InitFail_' num2str(k) '_Balanced_New_AC.mat'], 'cascade_link_cum')
%     save(['IEEE2383_Capacity_1_Flow_1_InitFail_' num2str(k) '_Balanced_New.mat'], 'cascade_link_cum')
%     save(['IEEE118_Capacity_1_5_InitFail_' num2str(k) '_Balanced.mat'], 'cascade_link_cum')    
    save(['IEEERTS_Capacity_1_Flow_1_InitFail_' num2str(k) '_Balanced_New.mat'], 'cascade_link_cum')
    max(max(cascade_link_cum));
    
    cascade_train=cascade_link_cum(:,1:ceil(0.9*K));
%     index_train=[];
%     for i=1:size(cascade_train,2)
%         for j=i+1:size(cascade_train,2)
%             if sum(abs(cascade_train(:,i)-cascade_train(:,j)))==0
%                cascade_train(:,j)=0;
%                index_train=[index_train,j];
%             end
%         end
%     end
%     cascade_train(:,index_train)=[];
%     size(cascade_train,2)

%     save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_' num2str(k) '_Train_Balanced_New_Retry.mat'], 'cascade_train')
%    save(['IEEE2383_Capacity_1_Flow_1_InitFail_' num2str(k) '_Train_Balanced_New.mat'], 'cascade_train')
%     save(['IEEE1354_Capacity_1_5_Flow_1_InitFail_' num2str(k) '_Train_Balanced_New_AC.mat'], 'cascade_train')
%     save(['IEEE2383_Capacity_1_5_Flow_1_InitFail_' num2str(k) '_Train_Balanced_New_AC.mat'], 'cascade_train')
%     save(['IEEE118_Capacity_1_5_InitFail_' num2str(k) '_Train_Balanced.mat'], 'cascade_train')
    save(['IEEERTS_Capacity_1_Flow_1_InitFail_' num2str(k) '_Train_Balanced_New.mat'], 'cascade_train')
    cascade_test=cascade_link_cum(:,ceil(0.9*K)+1:K);
    
%     index_test=[];
%     for i=1:size(cascade_test,2)
%         for j=i+1:size(cascade_test,2)
%             if sum(abs(cascade_test(:,i)-cascade_test(:,j)))==0
%                cascade_test(:,j)=0;
%                index_test=[index_test,j];
%             end
%         end
%     end

%     save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_' num2str(k) '_Test_Balanced_New_Retry.mat'], 'cascade_test')
%    save(['IEEE2383_Capacity_1_Flow_1_InitFail_' num2str(k) '_Test_Balanced_New.mat'], 'cascade_test')
%     save(['IEEE1354_Capacity_1_5_Flow_1_InitFail_' num2str(k) '_Test_Balanced_New_AC.mat'], 'cascade_test')
%     save(['IEEE2383_Capacity_1_5_Flow_1_InitFail_' num2str(k) '_Test_Balanced_New_AC.mat'], 'cascade_test')
%     save(['IEEE118_Capacity_1_5_InitFail_' num2str(k) '_Test_Balanced.mat'],'cascade_test')
    save(['IEEERTS_Capacity_1_Flow_1_InitFail_' num2str(k) '_Test_Balanced_New.mat'], 'cascade_test')
    
    S_mean=mean(S_vec)
    S_med=median(S_vec)
    S_min=min(S_vec)
    S_max=max(S_vec)
    
    end
    
    
