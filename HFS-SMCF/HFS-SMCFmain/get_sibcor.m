function [sibcor] = get_sibcor(hiercor,Y,tree)
%根据样本特征间的相似度计算兄弟节点间的相似度,不是兄弟用0表示
internalNodes = tree_InternalNodes(tree);
internalNodes(find(internalNodes==-1))=[];
sibcor=zeros(length(Y),length(Y));
for i=1:length(internalNodes)
    if isempty(Y{internalNodes(i)})
        continue
    end
    cur_sib = tree_Sibling(tree,internalNodes(i));%兄弟节点
    cur_par=tree(internalNodes(i),1);%根节点
    cur_index=find(Y{cur_par}==internalNodes(i));%非叶子节点对应实例序号
    cur_cor=hiercor{cur_par}(cur_index,:);%非叶子节点对应的实例相关性矩阵
    for j=1:length(cur_sib)%遍历兄弟节点
        sib_index=find(Y{cur_par}==cur_sib(j));%兄弟节点的实例序号
        if isempty(sib_index)
            continue
        end
        cur_sib_cor=cur_cor(:,sib_index);%兄弟节点相关性矩阵
        sibcor(internalNodes(i),cur_sib(j))=mean(mean(cur_sib_cor,2));%求均值，归结到32个节点上
    end
end
temp=sibcor;
temp_index=find(temp==0);
temp(temp_index)=[];
if length(temp)>0
sibcor=(sibcor-min(temp))./(max(temp)-min(temp));
end
sibcor(temp_index)=0;
end

