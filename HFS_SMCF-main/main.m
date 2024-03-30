clc; clear;
str1={'DD'};
%%%%%%%%%% Label distribution parameters
a=0.1;
b=0.2;
c=0.15;
%%%%%%%%%%%Coarse to fine split modeling parameters 
alpha=0.2;
beta=0.5;
gama=0.1;
%%%%%%%%%%%%%%
m = length(str1);
rng('default');
for i =1:m
    indx=[];
    filename = [str1{i} 'Train']    
    %     filename1 = [str1{i} 'cor1']
    load (filename);
    %     load (filename1);
    %Pearson similarity
    cor=corr(data_array(:,1:end-1)',data_array(:,1:end-1)','type','pearson');
    [X,Y,~,cor]=create_SubTable(data_array, tree,cor);
    tic;
    TDTime(i)= toc;
    %%%%%%%%%%%%
    [X2,Y2,~]=create_SubTable2(data_array, tree);
    [clus] = cluster(X,tree,data_array);
    %%%%%%%%%%%%
    [Yd] = create_hier_distribution(Y,tree,cor,a,b,c,X,clus,X2);
    [treecor] = get_treecor(tree,clus);
    [sibcor] = get_sibcor(cor,Y,tree);
    [labelsimi] = labelsim(tree,data_array);
    %Feature selectionFHStdclus
    tic;
    [feature_slct,W] = HFS_SMCF(X, Yd, tree, 10, treecor,sibcor,alpha,beta,gama,0,labelsimi);
    %%%%%%%%%%%%%%
    [clus2] = selectcluster(X2,tree,data_array,feature_slct);
    distances = pdist(clus2);
    % 对距离矩阵进行求和
    total_distance = sum(distances);
    %%%%%%%%%%%%%%%%%%%%%%%%
    FSTime(i) =toc;
    count_num=tabulate(data_array(:,end));
    max_num=max(count_num(:,2));
    mon_num=max_num*0.2;
    mon_index{i}=find(count_num(:,2)<=mon_num);
    %Test feature batch
    testFile = [str1{i}, 'Test.mat']
    load (testFile);
    [accuracyMean(i), accuracyStd(i), F_LCAMean(i), FHMean(i), TIEmean(i), TestTime(i),accuracy_l{i},accuracy_mon(i),FHStd(i), TIEStd(i),accuracy_monStd(i),FH{i},TIE{i}] = HierSVMPredictionBatch(data_array, tree, feature_slct,mon_index{i},str1{i});        %
    [t_r,~]=size(data_array);
    tiemean(i)=TIEmean(i)/t_r;
    tieStd(i)=TIEStd(i)./t_r;
    tie{i}=TIE{i}./t_r;
end
