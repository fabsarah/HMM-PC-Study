%% Let's correlate some TProb matrices!
%In this script, I extract transitional probability (TP) matrices with music data over time, and use these correlation matrices in PLS

% Step one: get the matrices
data = Rest_hmmPaper;%where I store the variables
vpath = data.vpath;%getting the vpath...see Step_One for how this is extracted
T = 792*20;
Xi = data.Xi;%This can be extracted when you run hmm the first time. I forgot, so am running it now
Gamma = data.Gamma;%The Gamma and hmm arguments are output from the initial hmm run
hmm = data.hmm;

testMask = 1:length(vpath);%Get the mask indices
Mask = reshape(testMask,[],T)';%reshape them task-wise
outMask = num2cell(Mask,2); %convert the Mask indices to a cell
[P,Pi] = getMaskedTransProbMats(vpath,T*(5120/20),hmm,outMask,Gamma,Xi);%get the matrices

%%
time_P = reshape(P,20,792)';%have checked this against the original!
time_Pi = reshape(Pi,20,792)';
%%
temp_P = cell(1,792);
%temp_Pi = cell(1,792);
for i = 1:792
    temptime = nan(144,20);
    for time = 1:20
        tempdata = time_P{i,time};
        tempvec = reshape(tempdata,[],1);
        temptime(:,time) = tempvec;
    end
    temp_P{i} = temptime;
end
%% Getting the matrices into participant*task order
index = Rest_hmmPaper.N_nonempty;
TPT_Ptime = cell(17,53);
%TPT_Pi = cell(17,53);
start = 1;
for part = 1:17
    dex = length(index{part});
    stop = (start-1)+dex;%how many cells to grab
    TPT_Ptime(part,index{part}) = temp_P(start:stop);%put those cells into task order (index) in the new matrix
    %TPT_Pi(part,index{part}) = Pi(start:stop);%same for Pi
    start = stop+1;
end
clear index part i dex time stop start temptime tempdata tempvec
Rest_hmmPaper.TPT_time.TPT_Ptime = TPT_Ptime;
%% It's correlation time!
%first, we need the MIR data
data = Rest_hmmPaper;
TP_Corrs = cell(17,40);
partdata = data.TPT_time.TPT_Ptime;
for piece = 1:40
    temp_MIR = nan(5,20);
    for i = 1:5
        temp_MIR(i,:) = Mir_corrs{i}(piece,:);
    end
    for part = 1:17
        tempdata = partdata{part,piece+6};
        if isempty(tempdata)
            continue
        else
        tempcorr = corr(tempdata',temp_MIR','Rows','all');    
        end
        TP_Corrs{part,piece} = tempcorr;
    end
end
clear data partdata piece temp_MIR i part tempdata tempcorr
Rest_hmmPaper.TPT_time.TP_Corrs = TP_Corrs;
%% Now, let's average these!
data = Rest_hmmPaper.TPT_time.TP_Corrs;
TP_task = cell(17,2);
for part = 1:17
    tempdata = data(part,:);
    ctrl = tempdata([1:5 36:40]);
    exp = tempdata(6:35);
    ctrl = cat(3,ctrl{:});
    exp = cat(3,exp{:});
    ctrl = nanmean(ctrl,3);
    exp = nanmean(exp,3);
    TP_task{part,1} = ctrl;
    TP_task{part,2} = exp;
end
Rest_hmmPaper.TPT_time.TP_task = TP_task;
clear data part tempdata ctrl exp
%% Test plot!
part = 1;
data = Rest_hmmPaper.TPT_time.TP_task;

figure
subplot(1,2,1)
shower_tile_plot(data{part,1})
ylabel('Vectorized TP Matrix')
xlabel('MIR State FO')
colorbar
title(sprintf('Part %d, Control',part),'FontSize',16)

subplot(1,2,2)
shower_tile_plot(data{part,2})
ylabel('Vectorized TP Matrix')
xlabel('MIR State FO')
colorbar
title(sprintf('Part %d, Experiment',part),'FontSize',16)
%%

figure
for i = 1:5
    subplot(2,5,i)
    plotdata = (data{part,1}(:,i));
    plotdata = reshape(plotdata,12,12);
    shower_tile_plot(plotdata);
    pbaspect([1,1,1])
    colorbar
    title(sprintf('Part %d, Control',part),'FontSize',14)
end

for i = 6:10
    subplot(2,5,i)
    plotdata = (data{part,2}(:,i-5));
    plotdata = reshape(plotdata,12,12);
    shower_tile_plot(plotdata);
    pbaspect([1,1,1])
    colorbar
    title(sprintf('Part %d, Experiment',part),'FontSize',14)
end
%% Now, to PLS! Get the data in PLS shape!
addpath(genpath('pls'));
data = Rest_hmmPaper.TPT_time.TP_task;
temp_indata = cell(17,2);
PLS_indata = cell(2,1);
for c = 1:2
    for i = 1:17
    tempdata = data{i,c};
    tempdata = reshape(tempdata,[],1);
    tempdata = tempdata';
    temp_indata{i,c} = tempdata;
    end
    PLS_indata{c} = cell2mat(temp_indata(:,c));
end

Rest_hmmPaper.TPT_time.PLS_indata = PLS_indata;
clear data c i tempdata
%% Do the PLS
clear option
option.method = 1;% mean-centred PLS
option.num_perm = 500;
option.num_boot = 100;
indata = cell2mat(PLS_indata);

tempres = pls_analysis({indata}, 17,2,[option]);%running the PLS
%%
Rest_hmmPaper.TPT_time.TP_res = tempres;
%% Adding the mouse data
Side_2 = Rest_hmmPaper.TPT.Input.Side_2;%just mouse
clear option
option.method = 3;% Behavioural PLS
option.num_perm = 500;
option.num_boot = 100;
option.stacked_behavdata = Side_2{1};% the behaviour matrix
indata = PLS_indata{1};

mousectrlres = pls_analysis({indata}, 17,1,[option]);%running the PLS
%%
Rest_hmmPaper.TPT_time.TP_mouseallres = mouseallres;
Rest_hmmPaper.TPT_time.TP_mousectrlres = mousectrlres;
Rest_hmmPaper.TPT_time.TP_mouseexpres = mouseexpres;
%% Plot it all!
res = Rest_hmmPaper.TPT_time.TP_mouseallres;%Rest_hmmPaper.TPT section
LV = 1;
p = res.perm_result.sprob(LV);
figure

subplot(4,3,[1,4])
z = res.boot_result.orig_corr;
bar(z(1:10,1))
hold on
yneg = res.boot_result.llcorr(1:10,LV);
ypos = res.boot_result.ulcorr(1:10,LV);
errorbar(1:length(z(1:10,LV)),z(1:10,LV),yneg-z(1:10,LV),ypos-z(1:10,LV),'.')
colorbar off
grid on
xlim([0 11])
xticks(1:10)
xticklabels(Rest_hmmPaper.BehavLabs(1:10))
pbaspect([2 3 3])
xtickangle(45)
title(sprintf('Full Model, Orig Corr, p = %d',p),'FontSize',16)
ylabel('Control','FontSize',16)

subplot(4,3,[7,10])
z = res.boot_result.orig_corr;
bar(z(11:20,1))
hold on
yneg = res.boot_result.llcorr(11:20,LV);
ypos = res.boot_result.ulcorr(11:20,LV);
errorbar(1:length(z(11:20,LV)),z(11:20,LV),yneg-z(11:20,LV),ypos-z(11:20,LV),'.')
colorbar off
grid on
xlim([0 11])
xticks(1:10)
xticklabels(Rest_hmmPaper.BehavLabs(11:20))
pbaspect([2 3 3])
xtickangle(45)
%title(sprintf('Full Model Exp, Orig Corr, p = %d',p),'FontSize',16)
ylabel('Experiment','FontSize',16)

x = (res.boot_result.compare_u(:,LV));
x = reshape(x,5,144)';
x(abs(x)<2.3) = nan;
for i = 1:5
    plotdex = [2,3,5,6,8,9];
    subplot(3,3,plotdex(i))
    plotdata = x(:,i);
    plotdata = reshape(plotdata,12,12);
    shower_tile_plot(plotdata);
    caxis([-4 4])
    colorbar
    pbaspect([1 1 1])
    yticks(1:12)
    yticklabels(12:-1:1)
    title(sprintf('TP * MIR FO State %d',i),'FontSize',16)
end

clear s K x y z res
%% Reliability:
%the experimental findings are in at -0.35; the control at
%-0.87, and the combined at -0.62!
x = cell2mat(PLS_indata);
y = cell2mat(Side_2);
[pls_repro,splitflag]=PLSCCA_TestTrain(x,y,500,0.6,0);
%% Plotting the input
plotc = mean(PLS_indata{1});
plote = mean(PLS_indata{2});
plotc = reshape(plotc,144,5);
plote = reshape(plote,144,5);
MIR_labs = {'Texture';'Mid-High Hz';'Low Hz';'Suddenness';'Predictability'};

figure
for i = 1:5
    subplot(2,5,i)
    plotdata = plotc(:,i);
    plotdata(plotdata<0) = -1;
    plotdata(plotdata>0) = 1;
    %plotdata(abs(plotdata)<0.012) = nan;
    shower_tile_plot(reshape(plotdata,12,12));
    yticklabels(12:-2:0)
    pbaspect([1 1 1])
    caxis([-0.03 0.03])
    colorbar
    ylabel('Brain State')
    title(sprintf('Control * %s',MIR_labs{i}),'FontSize',14)

    subplot(2,5,i+5)
    plotdata = plote(:,i);
    plotdata(plotdata<0) = -1;
    plotdata(plotdata>0) = 1;
    %plotdata(abs(plotdata)<0.012) = nan;
    shower_tile_plot(reshape(plotdata,12,12));
    pbaspect([1 1 1])
    yticklabels(12:-2:0)
    caxis([-0.03 0.03])
    colorbar
    ylabel('Brain State')
    title(sprintf('Experiment * %s',MIR_labs{i}),'FontSize',14)

end
%% Check variability (we don't want clusters)
LV = 1;

figure
sc_x = mouseallres.usc(1:17,LV);
sc_y = mouseallres.vsc(1:17,LV);
subplot(1,2,1)
scatter(sc_x, sc_y)
grid on
pbaspect([1 1 1])
title('Control','FontSize',18)

sc_x = mouseallres.usc(18:end,LV);
sc_y = mouseallres.vsc(18:end,LV);
subplot(1,2,2)
scatter(sc_x, sc_y)
grid on
pbaspect([1 1 1])
title('Experiment','FontSize',18)
clear sc_x sc_y LV
%% Plotting Syntax
K = 144;
LV = 1;
%this is plotting the bootstrap coefficient. It shows the weight of channel
%by frequency.

%data = tempres;
data = Rest_hmmPaper.TPT_time.TP_res;
x = (data.boot_result.compare_u(:,LV));
x = reshape(x,144,[])';
%x(abs(x)<4) = 0;
%z = data.boot_result.orig_corr(:,LV);
%yneg = data.boot_result.llcorr(:,LV);
%ypos = data.boot_result.ulcorr(:,LV);


z = data.boot_result.orig_usc(:,LV);
yneg = data.boot_result.llusc(:,LV);
ypos = data.boot_result.ulusc(:,LV);
%%
figure
subplot(1,2,1)
imagesc(x')
%caxis([-15 15])
colorbar
pbaspect([1 1 1])
ylabel('TP Correlations')
xticks(1:5)
xlabel('MIR FO')
title(sprintf('Brain States, LV %d, p = %d',LV,data.perm_result.sprob(LV)),...
    'FontSize',16)

subplot(1,2,2)
bar(z)
grid on
hold on
errorbar(1:length(z),z,yneg-z,ypos-z,'.')
ylabel('Design Scores')
xticklabels({'Rest','Control','Experiment'})
xlim([0,3])
%xticks(1:20)
%xticklabels(BehavLabs)
xtickangle(45)
xlabel('Condition')
title(sprintf('Latent Variable %d, Transitional Probs by Task',LV),'FontSize',16)

hold off
clear x y z K LV yneg ypos
