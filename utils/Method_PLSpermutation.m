% Function: PLS and do permutation.
%perform full PLS and plot variance in Y explained by top 15 components
%typically top 2 or 3 components will explain a large part of the variance
%(hopefully!)
% Updata date: 2022.09.18
% Email:luolongcao@163.com
function PLSperm_result = Method_PLSpermutation(X,Y,dim,rep,out_path,fig_name)
% X: the predictors in matrix X
% Y: the responses in matrix Y
% dim: dimension number of PLS.
% rep: permutation number.
% out_path:
% fig_name: figure name for storage.

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
Rsquared_real_y = cumsum(100*PCTVAR(2,1:dim));
% permutation testing to assess significance of PLS result as a function of

Rsq_y = [];
Pvalue_perm_y = ones(1,dim);
for i = 1:dim
    parfor j = 1:rep % automatic parallel calculate
    order=randperm(size(Y,1));
    Xp=X(order,:);
    [~,~,~,~,~,PCTVAR_p,~,~]=plsregress(Xp,Y,dim);

    temp_y=cumsum(100*PCTVAR_p(2,1:dim));
    Rsq_y(i,j) = temp_y(i);
    end
    Pvalue_perm_y(i) = length(find(Rsq_y(i,:)>=Rsquared_real_y(i)))/rep;
end

% %plot 1
% figure;
% subplot(2,1,1)
% plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
% set(gca,'Fontsize',14)
% str_xlab=strcat('PLS components');
% xlabel(str_xlab,'FontSize',14);
% ylabel(strcat('Variance'),'FontSize',14);
% text(1:dim,cumsum(100*PCTVAR(2,1:dim))-3,num2str(cumsum(100*PCTVAR(2,1:dim)).','%.1f'))
% grid on
% %plot 2
% subplot(2,1,2);
% plot(1:dim, Pvalue_perm_y,'ok','MarkerSize',8,'MarkerFaceColor','r');
% xlabel(str_xlab,'FontSize',14);
% ylabel('p-value','FontSize',14);
% text(1:dim,Pvalue_perm_y,num2str(Pvalue_perm_y.','%.3f'))
% grid on
% disp([num2str(rep),'permutation complated!']);
% temp_name = [out_path,datestr(datetime('now'),'yyyy-mm-dd-HH-MM'),'_PLS_permut_',fig_name,'.fig'];
% % savefig(temp_name);
% close all;%×¢ÊÍºó¿É±£Áôfigure
PLSperm_result.PCVAR = PCTVAR(2,:)';
PLSperm_result.Pvalue_perm = Pvalue_perm_y';
PLSperm_result.X=X;
PLSperm_result.Y=Y;
PLSperm_result.dim=dim;
PLSperm_result.rep=rep;
% temp_name = [out_path,datestr(datetime('now'),'yyyy_mm_dd_HH_MM'),'_permut_',fig_name,'_ROIs.mat'];save(temp_name,'PLSperm_result');

end