clearvars;clc

load('brian2.mat')

it = 1;

si(st<it) = [];
st(st<it) = [];

st(si > Ne) = [];
st = st-it;

si(si > Ne) = [];

id_vec = 1:Ne;
id_edges = [(id_vec-0.5)  id_vec(end)+0.5];

mfr_brian = histcounts(si,id_edges) / (T-it);

si_brian = si;
st_brian = st;


% t = 0 : dt : (T-0.2-dt);
% t_edges = [(t - 0.5*dt) t(end)+0.5*dt];   


[yvals, edges] = histcounts(mfr_brian,'Normalization','pdf');
xvals = edges(1 : (end-1) ) + diff(edges)/2;
xvals(yvals==0) = [];
yvals(yvals==0) = [];
mfr_hist_brian = [xvals; yvals];

%%
% load('NETSIM-test1/data/data.mat','st','si')
filepath = './data/00000002spk.bin';
fid = fopen( filepath, 'rb' );
si = fread( fid, inf, 'int32', 8 ); frewind( fid ); fseek(fid,4,'bof');
si = si+1;
st = fread( fid, inf, 'double', 4 );
fclose( fid );

si(st<it) = [];
st(st<it) = [];

st(si > Ne) = [];
si(si > Ne) = [];

st = st - it;

mfr_netsim = histcounts(si,id_edges)/(T-it);

st_netsim = st;
si_netsim = si;

%%
% figure
% histogram(mfr_brian);hold on;histogram(mfr_netsim);legend('brian','netsim');xlabel('firing rate (sp/s)');ylabel('count')
%set(gca,'XScale','log')

% disp('average mean firing rate +/- SD:')
% fprintf('brian: %d +/- %d sp/s \n', mean(mfr_brian), std(mfr_brian))
% fprintf('netsim: %d +/- %d sp/s \n', mean(mfr_netsim), std(mfr_netsim))

%%
[yvals, edges] = histcounts(mfr_netsim,'Normalization','pdf');
xvals = edges(1 : (end-1) ) + diff(edges)/2;
xvals(yvals==0) = [];
yvals(yvals==0) = [];
mfr_hist_netsim = [xvals; yvals];

% figure
% plot(st,si,'.k')
% title('netsim')

% figure
% plot(mfr_hist_brian(1,:),mfr_hist_brian(2,:),mfr_hist_netsim(1,:),mfr_hist_netsim(2,:));xlabel('mean spike rate (sp/s)');ylabel('count');legend('brian','netsim')

%%

mean_brian = mean(mfr_hist_brian(1,:));
sem_brian = std(mfr_hist_brian(1,:)) / sqrt(size(mfr_hist_brian,2));

mean_netsim = mean(mfr_hist_netsim(1,:));
sem_netsim= std(mfr_hist_netsim(1,:)) / sqrt(size(mfr_hist_netsim,2)) ;