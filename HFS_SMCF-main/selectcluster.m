function [clus2] = selectcluster(X2,tree,data,feature_slct)
internalNodes = newtree_InternalNodes(tree);%
internalNodes(find(internalNodes==-1))=[];%
indexRoot = tree_Root(tree);% The root of the tree?
noLeafNode =[internalNodes;indexRoot];%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m, numFeature] = size(data);
numFeature = numFeature - 1; 
numSeleted = round(numFeature * 0.1);
nt0=[];
for i = 1:length(noLeafNode)%
     %非叶子节点聚类
     if isempty(X2{noLeafNode(i)});
         X2{noLeafNode(i)}=[0];
     end
     selFeature = feature_slct{noLeafNode(i)}(1:numSeleted,:);
%      b=X2{noLeafNode(i)};
     if X2{noLeafNode(i)}==0;
          else
     a=X2{noLeafNode(i)}(:,selFeature);
      [ntl,nt0(noLeafNode(i),:)]=kmeans(X2{noLeafNode(i)}(:,selFeature),1,'EmptyAction','singleton','OnlinePhase','off','Display','off');
      children_set = get_children_set(tree, noLeafNode(i));%得到当前节点的叶子节点
      %叶子节点聚类
     
      for j = 1:length(children_set)
          indices=X2{noLeafNode(i)}(:,end);
          ntlID =  find(indices == children_set(j));%?
          if isempty(ntlID);
          else
          data=X2{noLeafNode(i)}(:,selFeature);
          ntl_data = data(ntlID,1:end);%
          [ntl,nt0(children_set(j),:)]=kmeans(ntl_data,1,'EmptyAction','singleton','OnlinePhase','off','Display','off');
          end
      end
     end
end
% a=tree_LeafNode( tree );
% for i = 1:size(a,2)
%     indices=data(:,end);
%     ntlID = (indices == i);%?
%     ntl_data = data(ntlID,1:end-1);%
%     %叶子节点聚类
%     if isempty(ntl_data);
%          ntl_data=[0];
%      end
%     [tl,nt0(i,:)]=kmeans(ntl_data,1,'EmptyAction','singleton','OnlinePhase','off','Display','off');
% end
clus2=nt0;