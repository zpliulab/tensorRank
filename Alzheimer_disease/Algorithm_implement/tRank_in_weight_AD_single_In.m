clear; clc;
close all;
%%  data import
Ast = importdata('select_In.csv');
% Ast = importdata('select_In_no_impute.csv');
% Ast = importdata('PPwith_imputaion_dmi_new_small.csv');
Ast_net = Ast.textdata(2:end,:);
Ast_weight = Ast.data(:,1);
%%  unique all nodes
Node1=[];j =1;
% Ast_net=[];cc=1;
for(i =1:size(Ast_net,1))
      Node1{j,1} =Ast_net{i,1};
                       j = j+1;
      Node1{j,1} = Ast_net{i,2};
                       j = j+1;
end
Node=unique(Node1);
Ad1= importdata('AD_disease.txt');
Ad2 =importdata('23_AD.txt');
Ad1=Ad1(1:end);
Ad2=Ad2(1:end);
loc2=[];j = 1;
for(ii = 1:length(Ad1))
     z2=strcmp(Ad1{ii},Node);
    z3=find(z2==1);
    if(~isempty(z3))
        loc2{j,1} = Ad1{ii};
                j = j+1;
    end
end
Ad3= importdata('of_M_100.txt');
 Ad3=Ad3(1:100);
Ad2=Ad2(1:end);
% loc2=[];
for(ii = 1:length(Ad2))
     z2=strcmp(Ad2{ii},Node);
    z3=find(z2==1);
    if(~isempty(z3))
        loc2{j,1} = Ad2{ii};
                j = j+1;
    end
end
for(ii = 1:length(Ad3))
     z2=strcmp(Ad3{ii},Node);
    z3=find(z2==1);
    if(~isempty(z3))
        loc2{j,1} = Ad3{ii};
                j = j+1;
    end
end
loc=unique(loc2);
loc(41)=[];loc(39)=[];loc(25)=[];
loc(12)=[];loc(8)=[];
%% create adjacent matrix
% if it does not exist, we fill 0
restart = 0.85; layer=1;%0.85
num =length(Node);
%% Ast matrix
A1 = zeros(length(Node),length(Node));
A1w = zeros(length(Node),length(Node));
for(i =1:size(Ast_net,1))
    k1 = strcmp(Ast_net{i,1},Node);k11 = find(k1==1);
    k2 = strcmp(Ast_net{i,2},Node);k22 = find(k2==1);
    A1(k2,k1)=1;
    A1w(k2,k1)=Ast_weight(i);
end
%Normalize adjacent matrix
% Adjacent matrix add weight
Acol=sum(A1,1);
Awcol=sum(A1w,1);
% save('degree_Ast.mat','Acol','-v6')
At=A1;
for(i = 1:length(Acol))
    if(Acol(i) >0)
        for(j=1:length(Acol))
        At(j,i) =A1w(j,i)/Awcol(i);
        end
    else
        At(:,i) =1/length(Node);
    end
end

%% 
slice=1;
TT = zeros(length(Node),length(Node),layer);
% Ast Oli In
TT(1:num,1:num,1) =At;
%% 现在趋势是对的
delta_thresh =1e-8;
iter_max =100;iter = 0;
del = 200;
PP1 =1/(layer*num)*ones(layer*num,1,slice);
I_now1 =PP1;del_stro=[];count=1;
while(del > delta_thresh)
     tic
     pre_PR1=I_now1;
     I_now1= restart*TT*I_now1+(1-restart)*PP1;
     del=norm(I_now1-pre_PR1);
    del_stro(count,1)=del;
    count=count+1;
    iter=iter+1;  
    toc
end
%% （1）layer 1 from 3 different single information
Ast_PR=I_now1;
length(unique(Ast_PR))
k1 = sort(unique(Ast_PR),'descend');
Rank =[];
cc=1;
for(i=1:length(Ast_PR))
    zz1=find(Ast_PR(i)==k1);
    if(~isempty(zz1))
        Rank(cc,1)=zz1;
        cc=cc+1;
    end
end

% save('PP_weight.mat','Node','Ast_PR','Rank','-v6')
% save('Beta_weight.mat','Node','Ast_PR','Rank','-v6')
% save('Alpa_weight.mat','Node','Ast_PR','Rank','-v6')
%%
% loc=loc2;

%% disease PR and Rank
D_Rank=[];D_PR=[];cc=1;
for(ii=1:length(loc))
    kk1=strcmp(loc{ii},Node);
    pp1 = find(kk1==1);
    if(~isempty(pp1))
        D_Rank(cc,1)=Rank(pp1);
        D_PR(cc,1)=Ast_PR(pp1);
        cc=cc+1;
    end
end
boxplot(D_Rank)
Nor=setdiff(Node,loc);
P_Auc=[];
R_Rank=[];R_PR=[];
for(count=1:1000)
rd = randperm(length(Nor),length(loc));
rd_sy=[];
for(ii =1:length(rd))
    rd_sy{ii,1}=Nor{rd(ii)};
end
%% rando 的rank 和PR
cc=1;
for(ii =1:length(rd_sy))
    kk1=strcmp(rd_sy{ii},Node);
    zzz1=find(kk1==1);
    if(~isempty(zzz1))
      R_Rank(cc,count)=Rank(zzz1);
      R_PR(cc,count)=Ast_PR(zzz1);
      cc=cc+1;
    end
end

boxplot([D_Rank,R_Rank(:,count)])
[h,p,ci,stats] = ttest2(D_Rank,R_Rank(:,count));
P_Auc(count,1)=p;
a_pr =[D_PR',R_PR(:,count)']';
resp = (1:(2*length(D_Rank)))'<length(D_Rank);
[X,Y,T,AUC] = perfcurve(resp,a_pr,'true');
% AUC
% plot(X,Y,'b')
% xlabel('False positive rate'); ylabel('True positive rate')
P_Auc(count,2)=AUC;
end
boxplot(P_Auc(:,2))
mean(P_Auc(:,2))
std(P_Auc(:,2))
length(find(P_Auc(:,1)<0.05))
save In_weighted.mat  Node Rank D_PR D_Rank R_PR R_Rank -v6
% save Beta_weigthed.mat  D_PR D_Rank R_PR R_Rank -v6
% save 
% save 'Ex_net_big.mat'
% save 'Ex_net_B_imputation_weighted_DMI.mat'
%% 排名靠前的node
% loo =find(Rank<=30);
% Former_node=[];
% for(j =1:length(loo))
%     Former_node{j,1}=Node{loo(j)};
% end
% fid=fopen('AD_In_top_30_weight.txt','a');
% for jj=1:length(Former_node)
%     fprintf(fid,'%s\r\n',Former_node{jj});   %
% %     fprintf(fid,'%.4\t',A(jj)); 
% end
% fclose(fid);
% % save('PR_single_Ast.mat','Ast_PR','-v6')
% Ot_PR=zzza((num+1):(2*num));
% length(unique(Ot_PR))
% It_PR=zzza((2*num+1):(3*num));
% length(unique(It_PR))
% save('Node_name_Ast.mat','Node','-v6')

% save('Ot_PR_single_new4.mat','Ot_PR','-v6')
% save('It_PR_single_new4.mat','It_PR','-v6')
%% (2) layer 2 from 2 different information At and Ot
% Ast_2_PR=zzzb(1:num);
% length(unique(Ast_2_PR))
% Ot_2_PR=zzzb((num+1):(2*num));
% length(unique(Ot_2_PR))
% % 利用几何均值来进行融合
% A_O_PR=[];
% for(i = 1:length(Ast_2_PR))
%     A_O_PR(i,1)=sqrt(Ast_2_PR(i)*Ot_2_PR(i));
% end
% length(unique(A_O_PR))
% % save('A_O_PR_two_new4.mat','A_O_PR','-v6')
% %% (3) layer 3 from 2 different information At and It 
% Ast_3_PR=zzzc(1:num);
% length(unique(Ast_3_PR))
% It_3_PR=zzzc((2*num+1):(3*num));
% length(unique(It_3_PR))
% %利用几何均值来进行融合
% A_I_PR=[];
% for(i = 1:length(Ast_3_PR))
%     A_I_PR(i,1)=sqrt(Ast_3_PR(i)*It_3_PR(i));
% end
% length(unique(A_I_PR))
% % save('A_I_PR_two_new4.mat','A_I_PR','-v6')
% % save('Node_name_new4.mat','Node','-v6')
%% rank index and rank 
% 保存为.mat 文件
% 用R 进行之后的处理

