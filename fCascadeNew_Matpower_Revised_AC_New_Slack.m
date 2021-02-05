function [yield, adjacency, history, powerInj, S] = fCascadeNew_Matpower_Revised_AC_New_Slack(mpc, powerInj, adjacency, reactance, capacity, capacity_vector, alpha, initFail,tmp_index)
    
    history = adjacency; 
    shedding_factor=0.8;
    count=0;
    while 1
        count=count+1;
        result=runacpf_me(mpc);
        result.success;
        if ~result.success
            mpc.bus(:,3)=mpc.bus(:,3)*shedding_factor;
            mpc.bus(:,4)=mpc.bus(:,4)*shedding_factor;
            mpc.gen(:,2)=mpc.gen(:,2)*shedding_factor;
            mpc.gen(:,3)=mpc.gen(:,3)*shedding_factor;
        else
            break;
        end
        if count>10
            mpc.bus(:,3)=0;
            mpc.bus(:,4)=0;
            mpc.gen(:,2)=0;
            mpc.gen(:,3)=0;
           break; 
        end
    end
    result.branch=[result.branch,mpc.branch(:,14)];

    initDem = -sum( min(powerInj, 0) );
        
    tempur_vector=1/2*(sqrt(result.branch(:,14).^2+result.branch(:,15).^2)+sqrt(result.branch(:,16).^2+result.branch(:,17).^2));
    
    tempur=zeros(size(result.bus,1));
    for i=size(tempur_vector,1)
        tempur(result.branch(i,1),result.branch(i,2))=tempur_vector(i,1);
        tempur(result.branch(i,2),result.branch(i,1))=tempur_vector(i,1);        
    end
    
    for i=1:size( initFail, 1)
        adjacency( initFail(i,1), initFail(i,2) ) = 0;
        adjacency( initFail(i,2), initFail(i,1) ) = 0;
    end
        
    link_flag_vector=ones(size(result.branch,1),1);
    
    link_flag_vector(tmp_index,1)=0;
    for t=1:size(tmp_index,1)
        mpc.branch(find(mpc.branch(:,14)==tmp_index(t)),:)=[];
    end
    
    link_row=mpc.branch(:,1);
    link_column=mpc.branch(:,2);
    G=digraph(link_row',link_column');
    [C,weak_size] = conncomp(G,'Type','weak');
    S=size(weak_size,2);

    %Load Shedding
    if S>1
        result_tmp=result;
        result_tmp.bus=[];
        result_tmp.gen=[];
        result_tmp.branch=[];
        for i=1:S
            CC_Index=find(C==i);
            CC_Load=mpc.bus(CC_Index,3);
            CC_Gen=[];
            CC_Gen_Index=[];
            for j=1:size(mpc.gen(:,1),1)
                 if size(find(CC_Index==mpc.gen(j,1)),2)==1
                     CC_Gen=[CC_Gen;mpc.gen(j,2)];
                     CC_Gen_Index=[CC_Gen_Index;j];
                 end
            end
            if sum(CC_Load)>sum(CC_Gen)
                if sum(CC_Gen)>0
                   CC_Load=CC_Load*sum(CC_Gen)/sum(CC_Load);
                   mpc.bus(CC_Index,3)=CC_Load;
                else
                   mpc.bus(CC_Index,3)=0;
                end
            else
                if sum(CC_Load)<sum(CC_Gen)
                    if sum(CC_Load)>0
                       CC_Gen=CC_Gen*sum(CC_Load)/sum(CC_Gen);
                       mpc.gen(CC_Gen_Index,2)=CC_Gen;
                    else
                       mpc.gen(CC_Gen_Index,2)=0;
                    end
                end
            end
            
            mpc_i=mpc;
            mpc_i.bus=mpc.bus(CC_Index,:);
            mpc_i.gen=mpc.gen(CC_Gen_Index,:);
            mpc_i.gencost=mpc.gencost(CC_Gen_Index,:);

            if size(find(mpc_i.bus(:,2)==3),1)==0 && size(mpc_i.bus(:,2),1)>1
                if size(find(mpc_i.bus(:,2)==2),1)>0
                    tmp_index_new=find(mpc_i.bus(:,2)==2);
                    index_slack_set=find(mpc_i.bus(tmp_index_new,3)==min(mpc_i.bus(tmp_index_new,3)));
                    index_slack=index_slack_set(1,1);
                    mpc_i.bus(tmp_index_new(index_slack),2)=3;
                end
            end
            
            CC_Branch_Index=[];
            for t=1:size(mpc.branch,1)
                if size(find(mpc_i.bus(:,1)==mpc.branch(t,1)),1)>0 && size(find(mpc_i.bus(:,1)==mpc.branch(t,2)),1)>0
                    CC_Branch_Index=[CC_Branch_Index;t];
                end
            end
            mpc_i.branch=mpc.branch(CC_Branch_Index,:);

            if size(mpc_i.bus,1)~=0 &&  size(mpc_i.gen,1)~=0 && size(mpc_i.branch,1)~=0
                count=0;
                while 1
                    count=count+1;
                    result_i=runacpf_me(mpc_i);
                    result_i.success;
                    if ~result_i.success
                        mpc_i.bus(:,3)=mpc_i.bus(:,3)*shedding_factor;
                        mpc_i.bus(:,4)=mpc_i.bus(:,4)*shedding_factor;
                        mpc_i.gen(:,2)=mpc_i.gen(:,2)*shedding_factor;
                        mpc_i.gen(:,3)=mpc_i.gen(:,3)*shedding_factor;
                        mpc.bus(CC_Index,3)=mpc.bus(CC_Index,3)*shedding_factor;
                        mpc.bus(CC_Index,4)=mpc.bus(CC_Index,4)*shedding_factor;
                        mpc.gen(CC_Gen_Index,2)=mpc.gen(CC_Gen_Index,2)*shedding_factor;
                        mpc.gen(CC_Gen_Index,3)=mpc.gen(CC_Gen_Index,3)*shedding_factor; 
                    else
                        break;
                    end
                    if count>10
                        mpc_i.bus(:,3)=0;
                        mpc_i.bus(:,4)=0;
                        mpc_i.gen(:,2)=0;
                        mpc_i.gen(:,3)=0;
                        mpc.bus(CC_Index,3)=0;
                        mpc.bus(CC_Index,4)=0;
                        mpc.gen(CC_Gen_Index,2)=0;
                        mpc.gen(CC_Gen_Index,3)=0; 
                        break; 
                    end
                end
                result_i.branch=[result_i.branch,mpc_i.branch(:,14)];
                result_tmp.bus=[result_tmp.bus;result_i.bus];
                result_tmp.gen=[result_tmp.gen;result_i.gen];
                result_tmp.branch=[result_tmp.branch;result_i.branch]; 
            else 
                mpc_i.gen(:,2)=0;
                mpc_i.gen(:,3)=0;
                mpc_i.bus(:,3)=0;
                mpc_i.bus(:,4)=0;
                mpc.bus(CC_Index,3)=0;
                mpc.bus(CC_Index,4)=0;
                mpc.gen(CC_Gen_Index,2)=0;
                mpc.gen(CC_Gen_Index,3)=0;               
                result_tmp.bus=[result_tmp.bus;mpc_i.bus];
                result_tmp.gen=[result_tmp.gen;mpc_i.gen];
                result_tmp.branch=[result_tmp.branch;[mpc_i.branch(:,1:size(mpc_i.branch,2)-1),zeros(size(mpc_i.branch,1),4),mpc_i.branch(:,size(mpc_i.branch,2))]];
            end
        end
        result.bus=sortrows(result_tmp.bus,1);
        result.gen=sortrows(result_tmp.gen,1);
        result.branch=sortrows(result_tmp.branch,18); %%%Note
    else
        count=0;
        while 1
            count=count+1;
            result=runacpf_me(mpc);
            result.success;
            if ~result.success
                mpc.bus(:,3)=mpc.bus(:,3)*shedding_factor;
                mpc.bus(:,4)=mpc.bus(:,4)*shedding_factor;
                mpc.gen(:,2)=mpc.gen(:,2)*shedding_factor;
                mpc.gen(:,3)=mpc.gen(:,3)*shedding_factor;
            else
                break;
            end
            if count>10
                mpc.bus(:,3)=0;
                mpc.bus(:,4)=0;
                mpc.gen(:,2)=0;
                mpc.gen(:,3)=0;
               break; 
            end
        end        
       result.branch=[result.branch,mpc.branch(:,14)];
       result.branch=sortrows(result.branch,18);
    end

    for i=1:size(result.branch,1)
        if isnan(result.branch(i,14))
            result.branch(i,14)=0;
            result.branch(i,16)=0;
        end
        if isnan(result.branch(i,15))
            result.branch(i,15)=0;
            result.branch(i,17)=0;                
        end
    end
     
    tempur = tempur .* adjacency;
    history = history + adjacency;

    Node_index=zeros(size(adjacency,1),1);
    for i=1:size(adjacency,1)
        Node_index(i,1)=i; 
    end

    bOverload = 1;
    while bOverload

        flow_vector=1/2*(sqrt(result.branch(:,14).^2+result.branch(:,15).^2)+sqrt(result.branch(:,16).^2+result.branch(:,17).^2));
        powerFlow=zeros(size(mpc.bus,1));
        
        for i=1:size(flow_vector,1)
            powerFlow(result.branch(i,1),result.branch(i,2))=flow_vector(i,1);
            powerFlow(result.branch(i,2),result.branch(i,1))=flow_vector(i,1);        
        end
        
        tempur = alpha.*abs( powerFlow ) + (1-alpha).*tempur;

        bOverload =  max( max( (abs(powerFlow) > capacity).*adjacency ) );
        adjacency = ( tempur <= capacity ).*adjacency;
        
        fail_index = find(flow_vector>capacity_vector(find(link_flag_vector==1)));
        tmp_index=mpc.branch(fail_index,14);

        link_flag_vector(tmp_index,1)=0;
        for t=1:size(tmp_index,1)
            mpc.branch(find(mpc.branch(:,14)==tmp_index(t)),:)=[];
        end
        
        link_row=mpc.branch(:,1);
        link_column=mpc.branch(:,2);
        link_row=[link_row;Node_index];
        link_column=[link_column;Node_index];
        G=digraph(link_row',link_column');
        [C,weak_size] = conncomp(G,'Type','weak');
        S=size(weak_size,2);
        
        if S>1
            result_tmp=result;
            result_tmp.bus=[];
            result_tmp.gen=[];
            result_tmp.branch=[];
            for i=1:S
                CC_Index=find(C==i);
                CC_Load=mpc.bus(CC_Index,3);
                CC_Gen=[];
                CC_Gen_Index=[];
                for j=1:size(mpc.gen(:,1))
                     if size(find(CC_Index==mpc.gen(j,1)),2)==1
                         CC_Gen=[CC_Gen;mpc.gen(j,2)];
                         CC_Gen_Index=[CC_Gen_Index;j];
                     end
                end

                mpc_i=mpc;
                mpc_i.bus=mpc.bus(CC_Index,:);
                mpc_i.gen=mpc.gen(CC_Gen_Index,:);
                mpc_i.gencost=mpc.gencost(CC_Gen_Index,:);

                if size(find(mpc_i.bus(:,2)==3),1)==0 && size(mpc_i.bus(:,2),1)>1
                    if size(find(mpc_i.bus(:,2)==2),1)>0
                        tmp_index_new=find(mpc_i.bus(:,2)==2);
                        index_slack_set=find(mpc_i.bus(tmp_index_new,3)==min(mpc_i.bus(tmp_index_new,3)));
                        index_slack=index_slack_set(1,1);
                        mpc_i.bus(tmp_index_new(index_slack),2)=3;
                    end
                end

                CC_Branch_Index=[];
                for t=1:size(mpc.branch,1)
                    if size(find(mpc_i.bus(:,1)==mpc.branch(t,1)),1)>0 && size(find(mpc_i.bus(:,1)==mpc.branch(t,2)),1)>0
                        CC_Branch_Index=[CC_Branch_Index;t];
                    end
                end
                mpc_i.branch=mpc.branch(CC_Branch_Index,:);
                
                if size(mpc_i.bus,1)~=0 && size(mpc_i.gen,1)~=0 && size(mpc_i.branch,1)~=0
                    count=0;
                    while 1
                        count=count+1;
                        result_i=runacpf_me(mpc_i);
                        result_i.success;
                        if ~result_i.success
                            mpc_i.bus(:,3)=mpc_i.bus(:,3)*shedding_factor;
                            mpc_i.bus(:,4)=mpc_i.bus(:,4)*shedding_factor;
                            mpc_i.gen(:,2)=mpc_i.gen(:,2)*shedding_factor;
                            mpc_i.gen(:,3)=mpc_i.gen(:,3)*shedding_factor;
                            mpc.bus(CC_Index,3)=mpc.bus(CC_Index,3)*shedding_factor;
                            mpc.bus(CC_Index,4)=mpc.bus(CC_Index,4)*shedding_factor;
                            mpc.gen(CC_Gen_Index,2)=mpc.gen(CC_Gen_Index,2)*shedding_factor;
                            mpc.gen(CC_Gen_Index,3)=mpc.gen(CC_Gen_Index,3)*shedding_factor;                            
                        else
                            break;
                        end
                        if count>10
                            mpc_i.bus(:,3)=0;
                            mpc_i.bus(:,4)=0;
                            mpc_i.gen(:,2)=0;
                            mpc_i.gen(:,3)=0;
                            mpc.bus(CC_Index,3)=0;
                            mpc.bus(CC_Index,4)=0;
                            mpc.gen(CC_Gen_Index,2)=0;
                            mpc.gen(CC_Gen_Index,3)=0; 
                            break; 
                        end
                    end
                    result_i.branch=[result_i.branch,mpc_i.branch(:,14)];
                    result_tmp.bus=[result_tmp.bus;result_i.bus];
                    result_tmp.gen=[result_tmp.gen;result_i.gen];
                    result_tmp.branch=[result_tmp.branch;result_i.branch]; 
                else
                    mpc_i.gen(:,2)=0;
                    mpc_i.gen(:,3)=0;
                    mpc_i.bus(:,3)=0;
                    mpc_i.bus(:,4)=0;
                    mpc.bus(CC_Index,3)=0;
                    mpc.bus(CC_Index,4)=0;
                    mpc.gen(CC_Gen_Index,2)=0;
                    mpc.gen(CC_Gen_Index,3)=0; 
                    result_tmp.bus=[result_tmp.bus;mpc_i.bus];
                    result_tmp.gen=[result_tmp.gen;mpc_i.gen];
                    result_tmp.branch=[result_tmp.branch;[mpc_i.branch(:,1:size(mpc_i.branch,2)-1),zeros(size(mpc_i.branch,1),4),mpc_i.branch(:,size(mpc_i.branch,2))]];
                end
            end
            result.bus=sortrows(result_tmp.bus,1);
            result.gen=sortrows(result_tmp.gen,1);
            result.branch=sortrows(result_tmp.branch,18);
        else
            count=0;
            while 1
                count=count+1;
                result=runacpf_me(mpc);
                result.success;
                if ~result.success
                    mpc.bus(:,3)=mpc.bus(:,3)*shedding_factor;
                    mpc.bus(:,4)=mpc.bus(:,4)*shedding_factor;
                    mpc.gen(:,2)=mpc.gen(:,2)*shedding_factor;
                    mpc.gen(:,3)=mpc.gen(:,3)*shedding_factor;
                else
                    break;
                end
                if count>10
                    mpc.bus(:,3)=0;
                    mpc.bus(:,4)=0;
                    mpc.gen(:,2)=0;
                    mpc.gen(:,3)=0;
                   break; 
                end
            end
            result.branch=[result.branch,mpc.branch(:,14)];
            result.branch=sortrows(result.branch,18);
        end
        
        %Add a prevention for NaN
        for i=1:size(result.branch,1)
            if isnan(result.branch(i,14))
                result.branch(i,14)=0;
                result.branch(i,16)=0;
            end
            if isnan(result.branch(i,15))
                result.branch(i,15)=0;
                result.branch(i,17)=0;                
            end
        end
        
        tempur = tempur .* adjacency;
        history = history + adjacency;
        
        1;       
    end

    finalDem = -sum( min(powerInj, 0) );
    
    yield = finalDem/initDem;
end
