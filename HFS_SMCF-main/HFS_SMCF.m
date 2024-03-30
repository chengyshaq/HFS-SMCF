function [feature_slct,W] = HFS_SMCF(X, Y, tree, lambda, treecor,sibcor,alpha,beta,gama,flag,labelsimi)
q=1;ma=[];mi=[];
internalNodes = tree_InternalNodes(tree);
indexRoot = tree_Root(tree);% The root of the tree
noLeafNode =[internalNodes;indexRoot];
eps = 1e-8; % set your own tolerance
maxIte = 10;
[r,c]=size(treecor);
[r2,c2]=size(sibcor);
v1=ones(r,c);
v2=ones(r2,c2);
treecor=v1-treecor;
sibcor=v2-sibcor;
for i = 1:length(noLeafNode)
    children_set = get_children_set(tree, noLeafNode(i));
    m(noLeafNode(i)) = length(children_set);
end
maxm=max(m);
for i = 1:length(noLeafNode)
    [r,c]=size(Y{noLeafNode(i)});
    Ytemp=zeros(r,maxm-c);
    Y{noLeafNode(i)}=[Y{noLeafNode(i)} Ytemp];
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1;
[~,d] = size(X{indexRoot}); % get the number of features
I = eye(maxm);
%% initialize
for j = 1:length(noLeafNode)
    W{noLeafNode(j)} = rand(d, maxm); % initialize W
    %%
    XX{noLeafNode(j)} = X{noLeafNode(j)}' * X{noLeafNode(j)};
    XY{noLeafNode(j)} = X{noLeafNode(j)}' * Y{noLeafNode(j)};
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:maxIte
    %% Update the root node
    %% initialization
    for j = 1:length(noLeafNode)
        D{noLeafNode(j)} = diag(0.5./max(sqrt(sum(W{noLeafNode(j)}.*W{noLeafNode(j)},2)),eps));
    end
    %% Update the root node
    W_current1 = zeros(d,maxm);
    for k = 1:length(noLeafNode)
        if isempty(W{noLeafNode(k)})
            continue
        end
        %W_current =  W_current + R(noLeafNode(j), compareNode(jcj))*W{compareNode(jcj)};
        W_current1 =  W_current1 + treecor(indexRoot, noLeafNode(k))*W{noLeafNode(k)};
    end
    W{indexRoot} = inv(XX{indexRoot} + lambda * D{indexRoot}) * (XY{indexRoot}-alpha*W_current1);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Update the internal nodes
    
    for j = 1:length(internalNodes)
        cur_sib=tree_Sibling(tree,internalNodes(j));
        W_current1 = zeros(d,maxm);
        W_current2= zeros(d,maxm);
        compareNode1 = setdiff(noLeafNode,internalNodes(j));
        compareNode2=setdiff(cur_sib,internalNodes(j));
        for k = 1:length(compareNode1)
            if isempty(W{compareNode1(k)})
                continue
            end
            W_current1 =  W_current1 + treecor(internalNodes(j), compareNode1(k))*W{compareNode1(k)};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        W_current3= zeros(d,maxm);
        [r3,c3]=size(labelsimi);
        v3=ones(r3,c3);
        labelsimi=v3-labelsimi;
        for k = 1:length(compareNode2)
            if isempty(W{compareNode1(k)})
                continue
            end
            W_current3 =  W_current3 + labelsimi(internalNodes(j), compareNode2(k))*W{compareNode1(k)};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        for k = 1:length(compareNode2)
            if isempty(W{compareNode2(k)})
                continue
            end
            W_current2 =  W_current2 + sibcor(internalNodes(j), compareNode2(k))*W{compareNode2(k)};
        end
%            W{internalNodes(j)} = inv(XX{internalNodes(j)} + lambda * D{internalNodes(j)}) * (XY{internalNodes(j)}-alpha*W_current1-beta*W_current2-0.35*W_current3);
          W{internalNodes(j)} = inv(XX{internalNodes(j)} + lambda * D{internalNodes(j)}) * (XY{internalNodes(j)})-(alpha*W_current1+beta*W_current2+gama*W_current3);
      %  W{internalNodes(j)} = inv(XX{internalNodes(j)} + lambda * D{internalNodes(j)}) * (XY{internalNodes(j)}-alpha*W_current1-beta*W_current2);


        
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calculate magnitude gaps
        a1=[];a2=[];b1=[];b2=[];c1=[];c2=[];
        a=inv(XX{internalNodes(j)} + lambda * D{internalNodes(j)}) * (XY{internalNodes(j)});
        b=alpha*W_current1+beta*W_current2+0.1*W_current3;
        c=abs(a./b);
        c_vec = reshape(c, [], 1); %
        
        sorted_c_max = sort(c_vec, 'descend'); 
        top_10_values = sorted_c_max(1:500);  
        average_value_max = mean(top_10_values(:));
        
        sorted_c_min = sort(c_vec,'ascend');
        top_10_smallest_values = sorted_c_min(1:2000);
        average_value_min = mean(top_10_smallest_values(:));
        
        c1=average_value_max;c2=average_value_min;
        mc(q,:)=[c1,c2,c1/c2];
        q=q+1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
        
        
    end
    %% Print the value of object function if flag is 1.
    if (flag ==1)
        obj(i)=norm(X{indexRoot}*W{indexRoot}-Y{indexRoot},'fro')^2+lambda*L21(W{indexRoot});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        objj(i)=norm(W{indexRoot},'fro')^2;
        W_current1 = 0;
        for j = 1:length(noLeafNode)
            if isempty(W{noLeafNode(j)})
                continue
            end
            
            W_current1 =  W_current1 + treecor(indexRoot, noLeafNode(j))*norm(W{indexRoot}'*W{noLeafNode(j)},'fro');
        end
        objj(i)=objj(i)+alpha*W_current1;
        for j = 1:length(internalNodes)
            cur_sib=tree_Sibling(tree,internalNodes(j));
            W_current1 = 0;
            W_current2= 0;
            W_current3= 0;
            compareNode1 = setdiff(noLeafNode,internalNodes(j));
            compareNode2=setdiff(cur_sib,internalNodes(j));
%             objj(i)=objj(i)-(norm(X{internalNodes(j)}*W{internalNodes(j)}-Y{internalNodes(j)},'fro'))^2+lambda*L21(W{internalNodes(j)});
            if(alpha>eps)
                for k = 1:length(compareNode1)
                    if isempty(W{compareNode1(k)})
                        continue
                    end
                    W_current1 =  W_current1 + treecor(internalNodes(j), compareNode1(k))*norm(W{internalNodes(j)}'*W{compareNode1(k)}-I,'fro');
                end
                objj(i)=objj(i)+alpha*W_current1;
            end
            if(beta>eps)
                for k = 1:length(compareNode2)
                    if isempty(W{compareNode2(k)})
                        continue
                    end
                    W_current2 =  W_current2 + sibcor(internalNodes(j), compareNode2(k))*norm(W{internalNodes(j)}'*W{compareNode2(k)}-I,'fro');
                end
                objj(i)=objj(i)+beta*W_current2;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(beta>eps)
                for k = 1:length(compareNode2)
                    if isempty(W{compareNode2(k)})
                        continue
                    end
                    W_current3 =  W_current3 + sibcor(internalNodes(j), compareNode2(k))*norm(W{internalNodes(j)}'*W{compareNode2(k)}-I,'fro');
                end
                objj(i)=objj(i)+gama*W_current3;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%obj
for i = 1: length(noLeafNode)
    W1=W{noLeafNode(i)};
    W{noLeafNode(i)} = W1(:,1:m(noLeafNode(i)));%去掉无用的列
end

% clear W1;
for j = 1: length(noLeafNode)
    tempVector = sum(W{noLeafNode(j)}.^2, 2);
    [atemp, value] = sort(tempVector, 'descend'); % sort tempVecror (W) in a descend order
    xxx{noLeafNode(j)}=X{noLeafNode(j)}(:,value);
    clear tempVector;
    feature_slct{noLeafNode(j)} = value(1:end);
end
if (flag == 1)
    fontsize = 20;
    figure1 = figure('Color',[1 1 1]);
    axes1 = axes('Parent',figure1,'FontSize',fontsize,'FontName','Times New Roman');
    
    plot(obj,'LineWidth',4,'Color',[0 0 1]);
    xlim(axes1,[0.8 10]);
    %     ylim(axes1,[16000,36000]);%Cifar
    % set(gca,'yscale','log')
    set(gca,'FontName','Times New Roman','FontSize',fontsize);
    xlabel('Iteration number');
%     ylabel('F1');
    ylabel('DD-F1');
%     ylabel('F194-F1');
%     ylabel('CLEF-F1');
%     ylabel('VOC-F1');
end
if (flag == 1)
    fontsize = 20;
    figure1 = figure('Color',[1 1 1]);
    axes1 = axes('Parent',figure1,'FontSize',fontsize,'FontName','Times New Roman');
    
    plot(objj,'LineWidth',4,'Color',[0 0 1]);
    xlim(axes1,[0.8 10]);
    %     ylim(axes1,[16000,36000]);%Cifar
    % set(gca,'yscale','log')
    set(gca,'FontName','Times New Roman','FontSize',fontsize);
    xlabel('Iteration number');
%     ylabel('F2');
    ylabel('DD-F2');
%     ylabel('F194-F2');
%     ylabel('CLEF-F2');
%     ylabel('VOC-F2');
end
end



