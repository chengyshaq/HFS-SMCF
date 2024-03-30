function [clus] = cluster(X,tree,data)
internalNodes = newtree_InternalNodes(tree);%
internalNodes(find(internalNodes==-1))=[];%
indexRoot = tree_Root(tree);% The root of the tree?
noLeafNode =[internalNodes;indexRoot];%
for i = 1:length(noLeafNode)%
     %非叶子节点聚类
     if isempty(X{noLeafNode(i)});
         X{noLeafNode(i)}=[0];
     end
      [ntl,nt0(noLeafNode(i),:)]=kmeans(X{noLeafNode(i)},1,'EmptyAction','singleton','OnlinePhase','off','Display','off');
end
a=tree_LeafNode( tree );
for i = 1:size(a,2)
    indices=data(:,end);
    ntlID = (indices == i);%?
    ntl_data = data(ntlID,1:end-1);%
    %叶子节点聚类
    if isempty(ntl_data);
         ntl_data=[0];
     end
    [tl,nt0(i,:)]=kmeans(ntl_data,1,'EmptyAction','singleton','OnlinePhase','off','Display','off');
end
clus=nt0;