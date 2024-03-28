%% creatSubTable������ģ��
% Written by Yu Wang
% Modified by Hong Zhao
% Modified by Haoyang Liu
%% Creat subtable
function [DataMod,LabelMod,IndxMod]=create_SubTablezh2(dataset, tree)
Data = dataset(:,1:end);%������������
Label =  dataset(:,end);%������ǩ����
[numTrain,~] = size(dataset);%������
internalNodes = tree_InternalNodes(tree);%�������м�ڵ�
internalNodes(find(internalNodes==-1))=[];%�м�ڵ�Ϊ-1���ÿ�
indexRoot = tree_Root(tree);% The root of the tree���ĸ��ڵ�
noLeafNode =[internalNodes;indexRoot];%��Ҷ�ӽڵ�
for i = 1:length(noLeafNode)%ѭ������ÿһ����Ҷ�ӽڵ�
    cur_descendants = tree_Descendant(tree, noLeafNode(i));%��i����Ҷ�ӽڵ���ӽڵ�
    ind_d = 1;  % index for id subscript increment  id�±�����������
    id = [];        % data whose labels belong to the descendants of the current nodes���ǩ���ڵ�ǰ�ڵ�ĺ��������
    for n = 1:numTrain  %ѭ������ÿһ������
        if (ismember(Label(n), cur_descendants) ~= 0)%�жϸ������Ƿ����ڵ�i���ڵ�
            id(ind_d) =  n;
            ind_d = ind_d +1;
        end
    end
    Label_Uni_Sel = Label(id,:);%ѡ����i���м�ڵ��������Ӧ�ı�ǩ
    DataSel = Data(id,:);     %select relative training data for the current classifierΪ��ǰ������ѡ����Ӧ��ѵ������

    numTrainSel = size(Label_Uni_Sel,1);%��i����Ҷ�ӽڵ��Ӧ��������
    LabelUniSelMod = mylabel_modify_MLNP(Label_Uni_Sel, noLeafNode(i), tree);%% �õ���ǩ��ģ��
    % Get the sub-training set containing only relative nodes�õ�ֻ������Խڵ����ѵ����
    ind_tdm = 1;
    index = [];     % data whose labels belong to the children of the current nodes��ǩ���ڵ�ǰ�ڵ���ӽڵ������
    children_set = get_children_set(tree, noLeafNode(i));%�õ���ǰ�ڵ��Ҷ�ӽڵ�
    for ns = 1:numTrainSel%ѭ������ÿһ��ʵ��
        if (ismember(LabelUniSelMod(ns), children_set) ~= 0)
            index(ind_tdm) =  ns;
            ind_tdm = ind_tdm +1;%ͳ��ʵ������
        end
    end
    DataMod{noLeafNode(i)} = DataSel(index, :);   % Find the sub training set of relative to-be-classified nodes
                                                    %�ҵ���Եķ���ڵ����ѵ����
    LabelMod{noLeafNode(i)} = LabelUniSelMod(index, :);%Ҷ�ӽڵ��ǩģ��
    IndxMod{noLeafNode(i)} = id;%���ڵ�i����Ҷ�ӽڵ����������ģ��
end
end