% Function: PLS and do bootstrape.
% Updata date: 2022.09.18
% Updata date: 2022.11.5
% Updata date: order the gene symbol!!!!2023.04.12
% Updata date: order the gene symbol!!!!2023.05.19
% Email:luolongcao@163.com
function test2 = Method_PLSbootstrap(X,Y,dim,max_dim,boot,geneSymbol,out_path,fig_name)
% X: the predictors in matrix X
% Y: the responses in matrix Y
% dim:dimension number of PLS.
% boot: bootstrap number.
% geneSymbol: the column name of matrix X.

[~,~,XS,~,~,PCTVAR,~,stats]=plsregress(X,Y,dim);

% temp_dim = find(PCTVAR(2,:) ==max(PCTVAR(2,:)));
temp_dim = max_dim;
%store regions' IDs and weights in descending order of weight for both components:
[R1,~]=corr([XS(:,1),XS(:,2),XS(:,temp_dim)],Y);%XS: predicted X score

%align PLS components with desired direction for interpretability 
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
%     XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
%     XS(:,2)=-1*XS(:,2);
end
if R1(3,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,temp_dim)=-1*stats.W(:,temp_dim);
%     XS(:,temp_dim)=-1*XS(:,temp_dim);
end
%real weight and sorted by weight
% [PLS1w,x1] = sort(stats.W(:,1),'descend');%W: A p-by-ncomp matrix of PLS weights so that XS = X0*W.
% PLS1ids=gene_name(x1);
% [PLS2w,x2] = sort(stats.W(:,2),'descend');
% PLS2ids=gene_name(x2);
[PLS3w,x3] = sort(stats.W(:,temp_dim),'descend');
Symbol = geneSymbol(x3); % order the gene symbol!!!!2023.04.12
% PLS3ids=gene_name(x3);

%define variables for storing the (ordered) weights from all bootstrap runs
% PLS1weights=[];
% PLS2weights=[];
PLS3weights=[];

%start bootstrap
parfor i=1:boot
    myresample = randsample(size(X,1),size(X,1),1);%replacement，有放回抽样
%     res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [~,~,~,~,~,~,~,stats_b]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
%     temp1=stats_b.W(:,1);%extract PLS1 weights
%     newW1=temp1(x1); %order the newly obtained weights the same way as initial PLS 
%     if corr(PLS1w,newW1)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
%         newW1=-1*newW1;
%     end
%     PLS1weights=[PLS1weights,newW1];%store (ordered) weights from this bootstrap run
%     
%     temp2=stats_b.W(:,2);%extract PLS2 weights
%     newW2=temp2(x2); %order the newly obtained weights the same way as initial PLS 
%     if corr(PLS2w,newW2)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
%         newW2=-1*newW2;
%     end
%     PLS2weights=[PLS2weights,newW2];%store (ordered) weights from this bootstrap run
    
    temp3=stats_b.W(:,temp_dim);%extract PLS2 weights
    newW3=temp3(x3); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS3w,newW3)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW3=-1*newW3;
    end
    PLS3weights=[PLS3weights,newW3]; %store (ordered) weights from this bootstrap run    

end

%get standard deviation of weights from bootstrap runs
% PLS1sw=std(PLS1weights');
% PLS2sw=std(PLS2weights');
PLS3sw=std(PLS3weights');

%get bootstrap weights (Z)
% PLS1Z=PLS1w./PLS1sw';
% PLS2Z=PLS2w./PLS2sw';
PLS3Z=(PLS3w)./PLS3sw';
% PLS3Z=(PLS3w-mean(PLS3weights,2))./PLS3sw';
% out_genes_PLS1 = table(geneSymbol, PLS1Z);
% out_genes_PLS2 = table(geneSymbol, PLS2Z);
[~,idx] = sort(x3,'ascend');
test2.gene = Symbol(idx);
test2.Z = PLS3Z(idx);

% test2.P = [];
% for temp = 1:size(idx,1)
%     test2.P = [test2.P,length(find(PLS3weights(temp,:)>abs(PLS3w(temp))))/boot];
% %     test2.P = [test2.P,(1-length(intersect(find(PLS3weights(temp,:)>-abs(PLS3w(temp))),find(PLS3weights(temp,:)<abs(PLS3w(temp)))))/boot)/2];
% end
% test2.P = test2.P(idx);

% out_genes_PLS3 = table(Symbol(idx), PLS3Z(idx));temp_name = [out_path,'PLS_boot_',fig_name,'.csv'];writetable(out_genes_PLS3,temp_name);
% disp([num2str(boot),'bootstrap complated!']);
end