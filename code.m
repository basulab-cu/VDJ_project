clear;

%% correlate ncRNAs and Ig V genes
load('Ig_expr_proB.mat');
load('Ig_ncRNA_expr_proB.mat');

tag = 'Ighv';
idx = grep(igGenes.gene_name,['^',tag]);
genes = [igGenes(idx,:),fpkm(idx,:)];

idx = myStrfind(genes.gene_id,ncRNA.Geneid);
ncRNA_s = ncRNA(idx,:);

mu1 = [mean(genes{:,[8,10,12]},2),mean(genes{:,[9,11,13]},2)];
mu2 = [mean(ncRNA_s{:,[7,9,11]},2),mean(ncRNA_s{:,[8,10,12]},2)];

log2fc = [log2(mu1(:,2)./mu1(:,1)),log2(mu2(:,2)./mu2(:,1))];
myfigure2;hold on;
plot(log2fc(:,1),log2fc(:,2),'.','MarkerSize',10);
lim = [-6,8];
set(gca,'xlim',lim,'ylim',lim);
set(gca,'FontSize',15);
plot(lim,[0,0],'k--','LineWidth',0.5)
plot([0,0],lim,'k--','LineWidth',0.5)

% one sample proportion test
log2fc_s = log2fc(all(~isnan(log2fc)&~isinf(log2fc),2),:);
n = size(log2fc_s,1);
p = sum(log2fc_s(:,1)<0&log2fc_s(:,2)>0)/n;
p0 = 1/4;
z = (p-p0)/sqrt(p0*(1-p0)/n);
[~,pval] = ztest(z,0,1);
