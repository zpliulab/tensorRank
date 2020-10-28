clear; clc;
close all;
%%  data import
% Ast

Ast = importdata('Selected_Alpa.csv');
% Ast = importdata('Hepat_without_imputaion_dmi_1.csv');
% Ast = importdata('Hepawith_imputaion_dmi.csv');
% Ast = importdata('Hepawith_imputaion_dmi1.csv');
% Ast = importdata('Hepawith_imputaion_cor.csv');
% Ast = importdata('Hepawith_imputaion_cor1.csv');
% Ast = importdata('Plaswith_imputaion_cor.csv');
% Ast = importdata('Hepat_without_imputaion_cor1.csv');
Ast_net=Ast.textdata(2:end,1:2);
% Ast_net = Ast.textdata;
% Ast_weight = Ast.data;
% % In
% In = importdata('used_select_In.csv');
% In_net = In.textdata;
% In_weight = In.data; 
% Oli
% Oli = importdata('used_select_Oli.csv');
% Oli_net = Oli.textdata;
% Oli_weight = Oli.data; 
%%  unique all nodes
Node1=[];j =1;

for(i =1:size(Ast_net,1))
      Node1{j,1} =Ast_net{i,1};
                       j = j+1;
      Node1{j,1} = Ast_net{i,2};
                       j = j+1;
end

Node=unique(Node1);

% Ad2= importdata('unique_liver_disease_genes.txt');
Ad2 =importdata('liver_disease_new12.txt');
% Ad2= importdata('malacards_2.txt');
Ad2=Ad2(1:end);
loc2=[];j = 1;
for(ii = 1:length(Ad2))
     z2=strcmp(Ad2{ii},Node);
    z3=find(z2==1);
    if(~isempty(z3))
        loc2{j,1} = Ad2{ii};
                j = j+1;
    end
end
loc=unique(loc2);
% loc=setdiff(loc2,{'SLC19A2','CDKN2B'});
% loc2(19)=[];
% loc2(13)=[];
% loc=loc2;
% loc(23)=[];
% loc=setdiff(loc2,{'EIF3E','SERPINC1','CYCS','IFNAR1',...
%                           'MAP2K1','YWHAB'});%, 'YWHAE',...
               %       'YWHAH', 'YWHAG', 'YWHAQ', 'YWHAZ'});%setdiff(loc2,{'TF','EIF3E','MALAT1','SERPINC1'});%unique([loc1' loc2']');%
% loc=loc;%(1:27);%(1:41);
% fid=fopen('disease_in_dcor1.txt','a');
% for jj=1:length(loc)
%     fprintf(fid,'%s\r\n',loc{jj});   %
% %     fprintf(fid,'%.4\t',A(jj)); 
% end
% fclose(fid);
%% create adjacent matrix

restart = 0.85; layer=1;
num =length(Node);
%% Ast matrix
A1 = zeros(length(Node),length(Node));

for(i =1:size(Ast_net,1))
    k1 = strcmp(Ast_net{i,1},Node);k11 = find(k1==1);
    k2 = strcmp(Ast_net{i,2},Node);k22 = find(k2==1);
    A1(k2,k1)=1;
end
%Normalize adjacent matrix
% Adjacent matrix add weight
Acol=sum(A1,1);
% save('degree_Ast.mat','Acol','-v6')
At=A1;
for(i = 1:length(Acol))
    if(Acol(i) >0)
        for(j=1:length(Acol))
        At(j,i) =A1(j,i)/Acol(i);
        end
    else
        At(:,i) =1/length(Node);
    end
end
%% 
slice=1;
TT = zeros(layer*length(Node),layer*length(Node),slice);
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
% end

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
% save('PP_origin.mat','Node','Ast_PR','Rank','-v6')
% save('Beta_origin.mat','Node','Ast_PR','Rank','-v6')
% save('Alpa_origin.mat','Node','Ast_PR','Rank','-v6')
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
Nor=setdiff(Node,loc);
P_Auc=[];R_Rank=[];R_PR=[];
for(count=1:100)
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
resp = (1:(2*length(D_PR)))'<(length(D_PR));
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
save alp_origin_new.mat D_PR D_Rank loc Node Rank R_PR R_Rank '-v6'
% save Beta_origin.mat  D_PR D_Rank R_PR R_Rank -v6
%% 排名靠前的node
% load('Ex_net_small.mat')
% loo =find(Rank<=100);
% Former_node=[];
% for(j =1:length(loo))
%     Former_node{j,1}=Node{loo(j)};
% end
% fid=fopen('small_former_node100.txt','a');
% for jj=1:length(Former_node)
%     fprintf(fid,'%s\r\n',Former_node{jj});   %
% %     fprintf(fid,'%.4\t',A(jj)); 
% end
% fclose(fid);