function [Yd] = create_hier_distribution(Y,tree,cor,alpha,beta,ga,X,clus,X2)

internalNodes = newtree_InternalNodes(tree);
internalNodes(find(internalNodes==-1))=[];
indexRoot = tree_Root(tree);% The root of the tree
noLeafNode =[internalNodes;indexRoot];
for i = 1:length(noLeafNode)
    m(noLeafNode(i)) = length(find(tree(:,1)==noLeafNode(i)));
end
maxm=max(m);
for i = 1:length(noLeafNode)
    children_set = get_children_set(tree, noLeafNode(i));
    Y{noLeafNode(i)}=myconversionY01(Y{noLeafNode(i)},length(children_set),children_set);%extend 2 to [1 0]
end


for i=1:length(noLeafNode)
    if isempty(Y{noLeafNode(i)})
        continue
    end
    path=tree(noLeafNode(i),2);
    for j=1:size(Y{noLeafNode(i)},1)
        y1=Y{noLeafNode(i)}(j,:);
        y2=Y{noLeafNode(i)}(j,:)+path;
        y2=y2/sum(y2);
        %         cor_sum=sum(cor{noLeafNode(i)}(j,:))-1;
        y3=Y{noLeafNode(i)};
        y3_temp=[];
        cor_temp=cor{noLeafNode(i)}(j,:);
        y3(j,:)=[];
        cor_temp(j)=[];
        if length(cor_temp)>1
            cor_temp=(cor_temp-min(cor_temp))./(max(cor_temp)-min(cor_temp));
        end
        [~,col]=size(y3);
        y3_temp=y3.* repmat(cor_temp',[1,col]);
        
        y3=sum(y3,1);
        y3_temp=sum(y3_temp,1);
        y3=y3_temp./y3;
        y3(find(isnan(y3)==1))=0;
        y3=y3/sum(y3);
        %%
        %%%%%%%%%%%%%%%%%%5
        temp = X{noLeafNode(i)}(j,:);
        A = ones(size(temp,1),size(temp,2));
        B=ones(size(y1,1),size(y1,2))';
        B2=unique(X2{noLeafNode(i)}(:,end));
        if size(B2,1)==noLeafNode(1)-1
%             for p=noLeafNode(1):noLeafNode(end)
%                 dist(p,:)=A-abs(temp-clus(k,:))
%             end
            y4=zeros(size(y1,2),1);
        else
            dist=[];
            for k=1:size(B,1)
                dist(k,:)=A-abs(temp-clus(k,:));
                dist(k,:)=(dist(k,:)-min(dist(k,:)))./(max(dist(k,:))-min(dist(k,:)));
            end
            y4=sum(dist,2);
            y4=(y4-min(y4))./(max(y4)-min(y4));
        end
        y4=y4';
         
        if isempty(y4)
            y4=zeros(size(y1,1),size(y1,2));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Yd{noLeafNode(i)}(j,:)=(1-alpha-beta-ga)*y1+alpha*y2+beta*y3+ga*y4;
    end
end

