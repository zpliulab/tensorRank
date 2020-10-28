clear; clc;
close all;
%%  data import
% Ast
Alp = importdata('Selected_Alpa.csv');
Alp_net = Alp.textdata(2:end,:);
Alp_weight = Alp.data(:,1);

PP = importdata('Selected_PP.csv');
PP_net = PP.textdata(2:end,:);
PP_weight = PP.data(:,1);

%%  unique all nodes
%% alpha cell
Node_alp=[];j =1;
for(i =1:size(Alp_net,1))
      Node_alp{j,1} =Alp_net{i,1};
                       j = j+1;
      Node_alp{j,1} = Alp_net{i,2};
                       j = j+1;
end
Node_alp=unique(Node_alp);
%% PP cell
Node_PP=[];j =1;
for(i =1:size(PP_net,1))
      Node_PP{j,1} =PP_net{i,1};
                       j = j+1;
      Node_PP{j,1} = PP_net{i,2};
                       j = j+1;
end
Node_PP=unique(Node_PP);


Node=unique([Node_alp',Node_PP']);



% Ad2= importdata('unique_liver_disease_genes.txt');
% Ad2 =importdata('marker_genes.txt');
Ad2= importdata('liver_disease_new12.txt');
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
% fid=fopen('liver_disease.txt','a');
% for jj=1:length(loc)
%     fprintf(fid,'%s\r\n',loc{jj});   %
% %     fprintf(fid,'%.4\t',A(jj)); 
% end
% fclose(fid);

%% create adjacent matrix
%% Alp matrix
Alpm = zeros(length(Node),length(Node));
Alpmw = zeros(length(Node),length(Node));
for(i =1:size(Alp_net,1))
    k1 = strcmp(Alp_net{i,1},Node);k11 = find(k1==1);
    k2 = strcmp(Alp_net{i,2},Node);k22 = find(k2==1);
    Alpm(k2,k1)=1;
    Alpmw(k2,k1)=Alp_weight(i);
end
%Normalize adjacent matrix
% Adjacent matrix add weight
Alpmcol=sum(Alpm,1);
Alpmwcol=sum(Alpmw,1);
% save('degree_Ast.mat','Acol','-v6')
Alpt=Alpm;
for(i = 1:length(Alpmcol))
    if(Alpmcol(i) >0)
        for(j=1:length(Alpmcol))
        Alpt(j,i) =Alpmw(j,i)/Alpmwcol(i);
        %*A1(j,i)/Acol(i);%
        end
    else
        Alpt(:,i) =1/length(Node);
    end
end

%% PP matrix
PPm = zeros(length(Node),length(Node));
PPmw = zeros(length(Node),length(Node));
for(i =1:size(PP_net,1))
    k1 = strcmp(PP_net{i,1},Node);k11 = find(k1==1);
    k2 = strcmp(PP_net{i,2},Node);k22 = find(k2==1);
    PPm(k2,k1)=1;
    PPmw(k2,k1)=PP_weight(i);
end
%Normalize adjacent matrix
% Adjacent matrix add weight
PPmcol=sum(PPm,1);
PPmwcol=sum(PPmw,1);

PPt=PPm;
for(i = 1:length(PPmcol))
    if(PPmcol(i) >0)
        for(j=1:length(PPmcol))
        PPt(j,i) =PPmw(j,i)/PPmwcol(i);
        %*A1(j,i)/Acol(i);%
        end
    else
        PPt(:,i) =1/length(Node);
    end
end
%% 
restart = 0.85; layer=2;%0.85
num =length(Node);
TT = zeros(length(Node),length(Node),layer);
TT(1:num,1:num,1) =Alpt;
TT(1:num,1:num,2)=PPt;
TT=tensor(TT);
%% 现在趋势是对的
delta_thresh =1e-8;
iter_max =100;iter = 0;
del = 200;
PP1 =1/(num)*ones(num,1);
size(PP1)
dd = 1/layer *ones(layer,1);
I_now1 =PP1;
del_stro=[];count=1;
while(del > delta_thresh)
     tic
     pre_PR1=I_now1;
     KK = ttv(ttv(TT,pre_PR1,2),dd,2);
     zz = double(tenmat(KK,1));
     I_now1 = (1-restart)*PP1+restart*zz;
     del=norm(I_now1-pre_PR1) ;
     length(unique(I_now1))
    del_stro(count,1)=del;
    count=count+1;
    iter=iter+1;  
    toc
end
% X = tensor(rand(5,3,4,2));
% A = rand(5,1); B = rand(3,1); C = rand(4,1); D = rand(2,1);
% Y = ttv(X, A, 1) %<-- X times A in mode 1
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
%%
% loc=loc2;
% loc=setdiff(loc2,{'SLC19A2','CDKN2B'});
% loc2(19)=[];
% loc2(13)=[];
loc=unique(loc2);

% loc=setdiff(loc2,{'EIF3E','SERPINC1','CYCS','IFNAR1',...
%                           'MAP2K1','YWHAB'});
                       %, 'YWHAE',...
                    %  'YWHAH', 'YWHAG', 'YWHAQ', 'YWHAZ'});%{'TF','MALAT1',});%(1:33);%unique([loc1' loc2']');
% loc=loc;%(1:27);%(1:41);
% fid=fopen('liver_disease_new12.txt','a');
% for jj=1:length(loc)
%     fprintf(fid,'%s\r\n',loc{jj});   %
% %     fprintf(fid,'%.4\t',A(jj)); 
% end
% fclose(fid);
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
save alp_PP_weighted_new1.mat D_PR D_Rank loc Node Rank R_PR R_Rank  '-v6'
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
% fid=fopen('B_big_former_node30.txt','a');
% for jj=1:length(Former_node)
%     fprintf(fid,'%s\r\n',Former_node{jj});   %
% %     fprintf(fid,'%.4\t',A(jj)); 
% end
% fclose(fid);
% save('PR_single_Ast.mat','Ast_PR','-v6')
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

