    %% Flow Calculation (DC/AC)
    mpc = case1354pegase;
    mpc.branch=mpc.branch(:,1:13);
    mpc.bus=mpc.bus(:,1:13);
    mpc.gen=mpc.gen(:,1:21);
    
    beta=1;
    mpc.bus(:,3)=beta*mpc.bus(:,3);
    mpc.bus(:,4)=beta*mpc.bus(:,4);
    mpc.gen(:,2)=beta*mpc.gen(:,2);
    mpc.gen(:,3)=beta*mpc.gen(:,3);
    
    bus = mpc.bus;

    Node_Index_List=zeros(size(bus,1),2);
    for i=1:size(bus,1)
        Node_Index_List(i,1)=i;
        Node_Index_List(i,2)=bus(i,1);
        bus(i,1)=i;
    end
    bus=sortrows(bus,1);
    mpc.bus = bus;

    gen = mpc.gen;
    for i=1:size(gen,1)
        gen(i,1)=Node_Index_List(find(gen(i,1)==Node_Index_List(:,2)),1); 
    end
    gen=sortrows(gen,1);
    mpc.gen=gen;

    branch_tmp=mpc.branch;
    for i=1:size(branch_tmp,1)
        branch_tmp(i,1)=Node_Index_List(find(branch_tmp(i,1)==Node_Index_List(:,2)),1);
        branch_tmp(i,2)=Node_Index_List(find(branch_tmp(i,2)==Node_Index_List(:,2)),1);
    end
    mpc.branch=branch_tmp;
    
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

    powerInj_vector=powerInj_supply-powerInj_consume;
    
    powerInj_consume=powerInj_consume+sum(powerInj_vector)/size(powerInj_vector,1);
    powerInj_vector=powerInj_supply-powerInj_consume;
    bus(:,3)=powerInj_consume;
    mpc.bus=bus;
    
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
    
%     %save('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Adjacency.mat','Adjacency')
%     filename = 'IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Adjacency.xlsx';
%     xlswrite(filename,Adj);
%     %save('IEEE1354_Capacity_1_Flow_1_5_InitFail_2_RowColumnIndex.mat','RowColumnIndex')
%     filename = 'IEEE1354_Capacity_1_Flow_1_5_InitFail_2_RowColumnIndex.xlsx';
%     xlswrite(filename,RowColumnIndex);
%     load(['IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Train_Balanced_New_AC.mat'],'cascade_link_cum')
%     cascade_train=cascade_link_cum(:,1:200);
%     filename = 'IEEE1354_Capacity_1_Flow_1_5_InitFail_2_Train_Balanced_New_AC.xlsx';
%     xlswrite(filename,cascade_train);    
%     index_vector=sort(randperm(M),'ascend')';
%     mpc.branch=[mpc.branch,index_vector];
    
    result=rundcpf_me(mpc);
    result.branch=[result.branch,mpc.branch(:,14)];
    tempur_vector=1/2*(sqrt(result.branch(:,14).^2+result.branch(:,15).^2)+sqrt(result.branch(:,16).^2+result.branch(:,17).^2));  %%%%AC
    
    Delta=[tempur_vector,Capacity_vector];
    Link_index_tmp=find(tempur_vector-Capacity_vector>=0);

    alpha=1;
    K=5000;
    link_set=[link_row,link_column];
    
    for k=2:2
    cascade_link_cum=zeros(M,K);

    S_vec=zeros(K,1);
    for i=1:K
        rnd=randperm(M);
        tmp_index=(rnd(1:k))';
        initFail=link_set(tmp_index,:);
        initFail_Index=tmp_index;

        [yield, adjacency, history, powerInjected, S] = fCascadeNew_Matpower_Revised_Slack_AC(mpc, powerInjected_ori, adjacency_ori, reactance, Capacity, Capacity_vector, alpha, initFail,tmp_index);

        cascade=full(history);
        cascade_link=zeros(M,1);
        for j=1:M
            cascade_link(j,1)=cascade(link_row(j),link_column(j));
        end
        cascade_link_cum(:,i)=cascade_link;
        
        S_vec(i,1)=S;
    end
    toc
    
     save(['IEEE1354_Capacity_1_Flow_1_5_InitFail_' num2str(k) '_Balanced_New_AC.mat'], 'cascade_link_cum')
    cascade_train=cascade_link_cum(:,1:ceil(0.9*K));
     save(['IEEE1354_Capacity_1_5_Flow_1_InitFail_' num2str(k) '_Train_Balanced_New_AC.mat'], 'cascade_train')
    cascade_test=cascade_link_cum(:,ceil(0.9*K)+1:K);
     save(['IEEE1354_Capacity_1_5_Flow_1_InitFail_' num2str(k) '_Test_Balanced_New_AC.mat'], 'cascade_test')

    
    S_mean=mean(S_vec)
    S_med=median(S_vec)
    S_min=min(S_vec)
    S_max=max(S_vec)
    
    end
    
    
