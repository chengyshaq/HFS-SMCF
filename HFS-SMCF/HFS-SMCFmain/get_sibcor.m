function [sibcor] = get_sibcor(hiercor,Y,tree)
%������������������ƶȼ����ֵܽڵ������ƶ�,�����ֵ���0��ʾ
internalNodes = tree_InternalNodes(tree);
internalNodes(find(internalNodes==-1))=[];
sibcor=zeros(length(Y),length(Y));
for i=1:length(internalNodes)
    if isempty(Y{internalNodes(i)})
        continue
    end
    cur_sib = tree_Sibling(tree,internalNodes(i));%�ֵܽڵ�
    cur_par=tree(internalNodes(i),1);%���ڵ�
    cur_index=find(Y{cur_par}==internalNodes(i));%��Ҷ�ӽڵ��Ӧʵ�����
    cur_cor=hiercor{cur_par}(cur_index,:);%��Ҷ�ӽڵ��Ӧ��ʵ������Ծ���
    for j=1:length(cur_sib)%�����ֵܽڵ�
        sib_index=find(Y{cur_par}==cur_sib(j));%�ֵܽڵ��ʵ�����
        if isempty(sib_index)
            continue
        end
        cur_sib_cor=cur_cor(:,sib_index);%�ֵܽڵ�����Ծ���
        sibcor(internalNodes(i),cur_sib(j))=mean(mean(cur_sib_cor,2));%���ֵ����ᵽ32���ڵ���
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

