     close all;
clear;
%% Options
% Labels & Column name
SetCol      =  'Dataset';
MarkCol     =  'Time';
LabCol      =  'Label';
fileCol     =  'File Name';
RepCol      =  'Replicate Label';
ExpCol      =  'Experiment Name';
SamCol      =  'File Name';
% Analytes
DoRatio=false;
DRSel = strtrim({'nC17/Pr',	'nC18/Ph',	'	2MN/1MN',	'	1MP/9MP',	'	2+3MP/9MP',	'	2+3MF/4-MF',	'	1MF/4-MF',	'	2+3MDBT/4MDBT',	'	1MDBT/4MDBT','1,6+1,3-DMN/C2N'});



%  CompSel = strtrim({'Naphthalene','	Acenaphthylene','	Acenaphthene','	Fluorene','	Dibenzothiophene','	Phenanthrene','	Anthracene','	Fluoranthene','	Pyrene','	Benz(a)anthracene','	Chrysene','	Benzo[b]fluoranthene','	Benzo[k]fluoranthene','	Benzo[e]pyrene','	Benzo[a]pyrene','	Perylene','	Indeno(1,2,3-c,d)pyrene','	Dibenz[a,h]anthracene','	Benzo[g,h,i]perylene'});

%  CompSel = strtrim({'	n-C17','	pristane','	n-C18','	Phytane	'});

%  CompSel = strtrim({' C1-N','	C2-N','	C3-N','	C4-N','	C1-BT','	C2-BT','	C3-BT','	C4-BT','	C1-DBT','	C2-DBT','	C3-DBT','	C1-F','	C2-F','	C3-F','	C1-P','	C2-P','	C3-P','	C4-P','	C1-Py','	C2-Py','	C1-Ch','	C2-Ch','	C3-Ch'});

  CompSel = strtrim({'Naphthalene','	Acenaphthylene','	Acenaphthene','	Fluorene','	Dibenzothiophene','	Phenanthrene','	Anthracene','	Fluoranthene','	Pyrene','	Benz(a)anthracene','	Chrysene','	Benzo[b]fluoranthene','	Benzo[k]fluoranthene','	Benzo[e]pyrene','	Benzo[a]pyrene','	Perylene','	Indeno(1,2,3-c,d)pyrene','	Dibenz[a,h]anthracene','	Benzo[g,h,i]perylene','	n-C17','	pristane','	n-C18','	Phytane	',' C1-N','	C2-N','	C3-N','	C4-N','	C1-BT','	C2-BT','	C3-BT','	C4-BT','	C1-DBT','	C2-DBT','	C3-DBT','	C1-F','	C2-F','	C3-F','	C1-P','	C2-P','	C3-P','	C4-P','	C1-Py','	C2-Py','	C1-Ch','	C2-Ch','	C3-Ch','17alfa, 21beta-hopane'});

  % DRE HAS BEEN HERE AND IS STILL
%      CompSel = strtrim({'Naphthalene','d8-naphthalene','	d10-Acenaphthene','	Fluorene','	d10-Fluorene','	d8-Dibenzothiophene','	d10-Phenanthrene','	Pyrene','	d10-Pyrene','	d12-Chrysene','	d12-Benzo(k)fluoranthene','	d12-Benzo[g,h,i]perylene','	d8-Acenaphthylene','	d10-Anthracene','	d10-fluoranthene','	d12-Benz(a)anthracene','	d12-benzo(a)pyrene','	d12-indeno(1,2,3-c,d)pyrene'});
%% Sediment Sample Compounds 
% CompSel = strtrim({'Naphthalene','	Acenaphthylene','	Acenaphthene','	Fluorene','	Dibenzothiophene','	Phenanthrene','	Anthracene','	Fluoranthene','	Pyrene','	Benz(a)anthracene','	Chrysene','	Benzo[b]fluoranthene','	Benzo[k]fluoranthene','	Benzo[e]pyrene','	Benzo[a]pyrene','	Perylene','	Indeno(1,2,3-c,d)pyrene','	Dibenz[a,h]anthracene','	Benzo[g,h,i]perylene','	n-C17','	pristane','	n-C18','	Phytane	',' C1-N','	C2-N','	C3-N','	C4-N','	C1-BT','	C2-BT','	C3-BT','	C4-BT','	C1-DBT','	C2-DBT','	C3-DBT','	C1-F','	C2-F','	C3-F','	C1-P','	C2-P','	C3-P','	C4-P','	C1-Py','	C2-Py','	C1-Ch','	C2-Ch','	C3-Ch','17alfa, 21beta-hopane'});

%% Water Sample Compounds 
% CompSel = strtrim({'Naphthalene','	Acenaphthylene','	Acenaphthene','	Fluorene','	Dibenzothiophene','	Phenanthrene','	Anthracene','	Fluoranthene','	Pyrene','	Benz(a)anthracene','	Chrysene','	Benzo[b]fluoranthene','	Benzo[k]fluoranthene','	Benzo[e]pyrene','	Benzo[a]pyrene','	Perylene','	Indeno(1,2,3-c,d)pyrene','	Dibenz[a,h]anthracene','	Benzo[g,h,i]perylene','	n-C17','	pristane','	n-C18','	Phytane	',' C1-N','	C2-N','	C3-N','	C4-N','	C1-BT','	C2-BT','	C3-BT','	C4-BT','	C1-DBT','	C2-DBT','	C3-DBT','	C1-F','	C2-F','	C3-F','	C1-P','	C2-P','	C3-P','	C4-P','	C1-Py','	C2-Py','	C1-Ch','	C2-Ch','	C3-Ch','17alfa, 21beta-hopane'});


%% Filter Sample Compounds 

% All Compounds
%     CompSel = strtrim({'	Acenaphthylene','	Fluorene','	Dibenzothiophene','	Phenanthrene','	n-C17','	pristane','	n-C18','	Phytane	','Naphthalene',' C1-N','	C2-N','	C3-N','	C4-N','	C1-BT','	C2-BT','	C3-BT','	C1-DBT','	C2-DBT','	C3-DBT','	C1-F','	C2-F','	C3-F','	C1-P','	C2-P','	C3-P','	C4-P','	C1-Py','	C2-Py','	C1-Ch','	C2-Ch','2-MN','	1-MN','	3-MF','	2-MF','	1-MF','	4-MF','	4-MDBT','	2+3-MDBT','	1-MDBT','	3-MP','	2-MP','	9+4-MP','	1-MP','Pyrene'});
  
% PAHs
%        CompSel = strtrim({'Naphthalene',' C1-N','	C2-N','	C3-N','	C4-N',	'C1-BT','	C2-BT','	C3-BT', 'Fluorene','	C1-F','	C2-F','	C3-F','	Phenanthrene','	C1-P','	C2-P','	C3-P','	Dibenzothiophene','	C1-DBT','	C2-DBT','	C3-DBT'});
    % Removed PAHs: 'Pyrene','C1-Py','	C2-Py',' '	C4-P','	Chrysene','	C1-Ch','	C2-Ch'
 

 
% Ratios
%   CompSel = strtrim({'2-MN','	1-MN','	3-MF','	2-MF','	1-MF','	4-MF','	4-MDBT','	2+3-MDBT','	1-MDBT','	3-MP','	2-MP','	9+4-MP','	1-MP'});

% 1st Removed : ,'	Anthracene', ,'	Chrysene' ,'	Benzo[k]fluoranthene' ,' Benzo[a]pyrene', '	Perylene','	Indeno(1,2,3-c,d)pyrene', ,'	C4-BT','	C3-Ch','17alfa, 21beta-hopane'
% 2nd Removed : '	Acenaphthene', '	Dibenz[a,h]anthracene','	Benzo[g,h,i]perylene','	Fluoranthene','	Benz(a)anthracene','	Benzo[b]fluoranthene','	Benzo[e]pyrene',


% Changing of Compound Labels 
CompChange{1} =  strtrim({'Fluorene','	Dibenzothiophene','	Phenanthrene','Naphthalene','Pyrene','	Chrysene'}); % Labels Before
CompChange{2}  = strtrim({'C0-F','C0-DBT','C0-P','C0-N','C0-Py','C0-Ch'});                                       % Labels After

%%
NormComp={'17alfa, 21beta-hopane'};
NormSample={'Normalisation Sample','Experiment Name','File Name'};

DoQCNorm   = false;
DoSampNorm = true;

% Sample Set

leaveout = 'Removed';

% Plotting
AvgPlot =true;

XYplot=false;
Xax = {'BB/AA','29ba/30ab'};
Yax = {'DAS/RS','28ab/30ab'};
ShowConvexHulls=[0 1];
ShowSampleNames=false;

% Reporting
DoReporting=false;
ModelTitle  = {'Filter - PAHs'};
ModelInfo   = {'IS + T0 Normalisation + Averaging'};

Tables      = {};
Titles      = {};

%% Reading Folder
designFile = which('Design_NEG_Filter.xlsx');
destFolder = cd; % Folder where the data will be saved
folder   = fileparts(designFile);
fileList = dir(fullfile(lower(folder),'*.csv'));
fileList = fullfile(folder,{fileList.name})';
%% Reading files
for i=1:numel(fileList)
    [tt{i},nn{i}] = MassHunterImport(fileList{i},200);
    % Data addition
    if(i==1)
        TT=tt{i};
        NN=nn{i};
    else
        if(isequal(TT(1:2,:),tt{i}(1:2,:)))
            TT=[TT ;tt{i}(3:end,:)];
            NN=[NN ;nn{i}(3:end,:)];
        else
            error('Data inconsistent')
        end
    end
end

% Samples
HeadCol     = TT(2,:);
SampleName  = TT(3:end,ismember(HeadCol,'Name'));
FileName_batch    = TT(3:end,ismember(HeadCol,'Data File'));
FileName_batch    = strtok(FileName_batch,'.');
FileName_batch    = cellstr(FileName_batch);

% Compounds
Compounds   = TT(1,2:end);
Compounds   = Compounds(~ismember(Compounds,''));
%  Compounds   = squeeze(split(Compounds,'Results'));
Compounds=regexp(Compounds,'Results','split');
for i=1:length(Compounds), CLab(i,1)=Compounds{i}(1,1);end
Compounds   = CLab;
Compounds   = cellstr(Compounds (:,1));
Compounds   = strtrim(Compounds);
RtCol       = ismember(TT(2,:),'RT');
RespCol     = ismember(TT(2,:),'Resp.');

% Design
[~,~,Design]=xlsread(designFile,'Design');
Colhead            = Design(1,:);
FileName_design    = Design(2:end,ismember(Design(1,:),fileCol));
FileName_design    = strtok(FileName_design,'.');
LOCA               = ismember(FileName_design,FileName_batch);
LOCB               = ismember(FileName_batch,FileName_design);

Design             = Design(2:end,:);
TT                 = TT(3:end,:);
NN                 = NN(3:end,:);

Design             = Design(LOCA,:);
TT                 = TT(LOCB,:);
NN                 = NN(LOCB,:);

FileName_design = FileName_design(LOCA);
FileName_batch = FileName_batch(LOCB);

% Sorting Data according to descending
[~,orddesign] = sort(FileName_design);
[~,ordbatch]  = sort(FileName_batch);

FileName_design= FileName_design(orddesign);
FileName_batch= FileName_batch(ordbatch);
Design        = Design(orddesign,:);
TT            = TT(ordbatch,:);
NN            = NN(ordbatch,:);

if(~isequal(FileName_design,FileName_batch)), error('Design and Batches are not the same order'); end

% % Remove Samples
SamSel= ~ismember(Design(:,strcmp(Colhead,SetCol)),leaveout);

% Data Collection
DesignSel = Design(SamSel,:);
Rt          = NN(SamSel,RtCol);
RespRaw     = NN(SamSel,RespCol);

%% Read Analytes
[~,~,Analytes]=xlsread(designFile,'Analytes');
AnalyteCol = Analytes(2,:);
if(isempty(DRSel)), DRSel=AnalyteCol(1,strcmp(Analytes(1,:),'Diagnostic Ratio')); end
Analytes   = Analytes(3:end,:);
RoILabs = strtrim(Analytes(:,strcmp(AnalyteCol,'label')));
RoINorm = strtrim(Analytes(:,strcmp(AnalyteCol,'Normalisation Compound')));

if(~isequal(Compounds, RoILabs)), error('Compounds do not match Analytes'); end

RoINorm = strtrim(Analytes(:,strcmp(AnalyteCol,'Normalisation Compound')));

if(isempty(RoINorm)), error('Problem with Internal standard definition'); end 
    
DRCol = ismember(AnalyteCol,DRSel);
if(~isequal(numel(DRSel), sum(DRCol))), error('Something wrong with Ratio Def.'); end
DRsplit=regexp(DRSel,'/','split');

for i=1:sum(DRCol)
    RoIs{i,1}     = strcmp(Analytes(:,strcmp(AnalyteCol,DRSel{i})),DRsplit{i}(1,1));
    RoIs{i,2}     = strcmp(Analytes(:,strcmp(AnalyteCol,DRSel{i})),DRsplit{i}(1,2));
end
AnalyteCol = AnalyteCol(DRCol);

%% Data Normalization
if(~isempty(NormComp))
    NormComp=strtrim(NormComp);
    switch NormComp{1}
        case 'Sum'
            RespNorm = bsxfun(@rdivide,RespRaw,nansum(RespRaw,2));
        case 'Euclidian'
            nF = sqrt(sum(RespRaw.^2,2));
            nF = mean(nF)./nF;
            RespNorm = bsxfun(@rdivide,RespRaw,nF);
        case 'Max'
            RespNorm = bsxfun(@rdivide,RespRaw,nanmax(RespRaw'));
            
        case 'Internal Standard'
            [~,index] = ismember(RoINorm,RoILabs);
           RespNorm = nan(size(RespRaw,1),size(RespRaw,2));
          for i=1:numel(index), RespNorm(:,i) = bsxfun(@rdivide ,RespRaw(:,i),RespRaw(:,index(i))); end 
       
        otherwise
            if(numel(NormComp)>1 && isequal(sum(ismember(RoILabs,NormComp)),numel(NormComp)))
                RespNorm = bsxfun(@rdivide,RespRaw,nansum(RespRaw(:,ismember(RoILabs,NormComp)),2));
            elseif(numel(NormComp)==1 && isequal(sum(ismember(RoILabs,NormComp)),numel(NormComp)))
                RespNorm = bsxfun(@rdivide,RespRaw,RespRaw(:,ismember(RoILabs,NormComp)));
            else
                error('Normalisation Compounds not found')
            end
    end
    X=RespNorm;
else
    X=RespRaw;
end

%% QC Normalisation 
if(DoQCNorm)
    [~,index] = ismember(DesignSel(:,strcmpi(Colhead,'QC Sample')),DesignSel(:,strcmpi(Colhead,'File Name')));    
    RespQC = nan(size(X,1),size(X,2));   
    for i=1:numel(index)
        for j=1:size(X,2)
            RespQC(i,j) = bsxfun(@rdivide ,X(i,j),X(index(i),j));
        end
    end    
    X = RespQC;
end
%% Calculation of Diagnostic Ratios
if(DoRatio)
    for i=1:size(X,1)
        for j=1:size(RoIs,1)
            X_Ratio(i,j)= nansum(X(i,RoIs{j,1}))./nansum(X(i,RoIs{j,2}));
        end
    end
end

%%  Averages
if(AvgPlot)
    [a,b,c]     =   unique(DesignSel(:,ismember(Colhead,RepCol)));
    DesignSel   =   DesignSel(b,:);
    nSam=size(a,1);
    if(DoRatio),X=X_Ratio; end
    X(X == 0) = NaN;
    
    for i=1:nSam
        for j=1:size(X,2)
            n_Rep(i,j)=sum(~isnan(X(c==i,j)));
            n_Rep_Label{i,j}=strcat('             n =  ', num2str(sum(~isnan(X(c==i,j)))));
            X_avg(i,j)= nanmean(X(c==i,j),1);
            X_std(i,j)= nanstd(X(c==i,j),1);
            X_Rstd(i,j)= nanstd(X(c==i,j),1)./nanmean(X(c==i,j),1);
            %             X_Rstd_Weighted(i,j)= (nanstd(X(c==i,j),1)./nanmean(X(c==i,j),1))*(n_Rep(i)-1);
            %             X_var(i,j)= nanvar(X(c==i,j),1);
            %             X_var_weighted_Rep(i,j)=(n_Rep(i)-1)*X_var(i,j);
        end
    end
    X=X_avg;
    X_Error=X_std;
    
    n_Rep_Label(n_Rep==3)={''};
    n_Rep_Label=reshape(n_Rep_Label,[size(n_Rep,1),size(n_Rep,2)]);
    
    n_Rep_Label(n_Rep==0)={'             <LOD'};
    n_Rep_Label=reshape(n_Rep_Label,[size(n_Rep,1),size(n_Rep,2)]);
    
    
end% Calculating Mean and Std of all Samples


%% Sample Normalisation
if (DoSampNorm && any(ismember(Colhead,NormSample{1})))
    if(any(cellfun(@isnumeric, DesignSel(:,strcmp(Colhead,NormSample{1}))))),  error('Missing Normalisation Sample from Designfile'); end
for i=1:size(X,2)
        for j=1:size(X,1)
            Xcorr(j,i)= X(j,i)./X(ismember(DesignSel(:,strcmp(Colhead,NormSample{3})),DesignSel(j,strcmp(Colhead,NormSample{1}))),i)*100;
            if(AvgPlot), Xcorr_Error(j,i)= X_Error(j,i)./X(ismember(DesignSel(:,strcmp(Colhead,NormSample{3})),DesignSel(j,strcmp(Colhead,NormSample{1}))),i)*100; end 
            NormFactor(j,i)=X(ismember(DesignSel(:,strcmp(Colhead,NormSample{3})),DesignSel(j,strcmp(Colhead,NormSample{1}))),i);
            NormSampleCheck{j}=DesignSel(j,strcmp(Colhead,NormSample{1}));
        end
    end
    X=Xcorr;
    if(AvgPlot), X_Error=Xcorr_Error; end 
    Ylab='Procent(%) i forhold til uforvitret olie';
else if(DoRatio && ~DoSampNorm),Ylab='Diagnostic Ratio';  else  Ylab='Relative Intensity';
    end  
end
%% Barplot

% [a,b,c]=unique(DesignSel(:,strcmp(Colhead,NormSample{2})));
MarkLabs=DesignSel(:,strcmp(Colhead,MarkCol));
ExpLabs=DesignSel(:,strcmp(Colhead,ExpCol));
SamLab=DesignSel(:,strcmp(Colhead,SamCol));

% Sorting according to T0-T2
[~,ordplot] = sort(MarkLabs);
MarkLabs=MarkLabs(ordplot);
ExpLabs=ExpLabs(ordplot);
SamLab=SamLab(ordplot);
DesignSel=DesignSel(ordplot,:);
X=X(ordplot,:);

if(AvgPlot)
X_Error=X_Error(ordplot,:); 
n_Rep_Label = n_Rep_Label(ordplot,:);
end

[a,b,c]=unique(DesignSel(:,strcmp(Colhead,NormSample{2})));
% end

for i=1:numel(a)
    figH{i}=figure;
    if(DoRatio)
        [XSel,Ind] = ismember(AnalyteCol,DRSel);
        [~,ordcomp] = sort(Ind(XSel));
        CompLab = AnalyteCol(XSel);
        CompLab=CompLab(ordcomp);
        Xlab = 'Diagnostic Ratios';
%          if (~isequal(numel(Compounds(XSel))),numel(CompLab)), error('Problem with  CompSel definition '); end
    else
       [XSel,Ind] =  ismember(Compounds,CompSel);
       [~,ordcomp] = sort(Ind(XSel));
       CompLab = Compounds(XSel);
       CompLab=CompLab(ordcomp);
       Xlab = 'Compounds Groups';   
         if (~isequal(numel(Compounds(XSel)),numel(CompLab))), error('Problem with  CompSel definition '); end
    end
    
    
 % Changing Specific Compound Label
 if (~isempty(CompChange{1}) || ~isempty(CompChange{2}))
    
    [~,Ind] = ismember(CompLab,CompChange{1});
    CompChange{3} = find(Ind>0); 
    CompLab(CompChange{3}) = CompChange{2} (Ind(CompChange{3}));   
 end 
 
    if(AvgPlot)
        
        SamNames{i}=MarkLabs(ismember(c,i));
        
        Ydata{i}  =  X(ismember(c,i),XSel);  
        Ydata{i}  =  Ydata{i} (:,ordcomp);
        Ydata{i}(isnan(Ydata{i}))=0;
        
        Ystd{i}    =  X_Error(ismember(c,i),XSel);
        Ystd{i}   =  Ystd{i} (:,ordcomp);
        Ystd{i}(isnan(Ystd{i}))=0;
       
        Rep_Label =  n_Rep_Label(ismember(c,i),XSel);
        Rep_Label =  Rep_Label(:,ordcomp);
      
        h{i} = bar(Ydata{i}'); hold on;
%         ctr{i}   = bsxfun(@plus, 1:numel(CompLab), [h{i}.XOffset]');
        ctr{i}   = bsxfun(@plus, 1:size(Ydata{i},2), [h{i}.XOffset]');
        
        set(gca,'Xtick',[h{i}(1,1).XData]);
        xticklabels(CompLab);
        set(gca,'XTickLabelRotation',45);
        set(gca,'YTickLabelRotation',45);
        set(gca,'fontsize', 11);
        legend(h{i},MarkLabs(ismember(c,i)),'autoupdate','off','Location','SouthEastOutside');
         errorbar(ctr{i}, Ydata{i},Ystd{i},'.k','linewidth',1);
%         title(ExpLabs(b(i)));
        set(figH{i},'name',ExpLabs{b(i)});
        txHan = text(reshape(ctr{i},1,[]),reshape((Ydata{i}+Ystd{i}),1,[])*1.025, reshape(Rep_Label,1,[]));
        set(txHan,'horizontalalignment','center','verticalalignment','middle','clipping','on');
        set(txHan,'fontsize',9,'Rotation',90);
%         ylim([0,round(max(max(Ydata{i}+Ystd{i}))*1.05)])
        
       
        Rstd_Avg{i} = nanmean(X_Rstd(ismember(c,i),XSel));
        Rstd_Std{i} = nanstd(X_Rstd(ismember(c,i),XSel));
        Rstd_Skw{i} = skewness(X_Rstd(ismember(c,i),XSel),0,1);
        Rstd_Kur{i} = kurtosis(X_Rstd(ismember(c,i),XSel),0,1);
        
        
%         hb{i} = bar(Rstd_Avg{i}'); hold on;
%         ctrb{i}   = bsxfun(@plus, 1:size(Rstd_Avg{i},2), [hb{i}.XOffset]');
%         errorbar(ctrb{i}, Rstd_Avg{i},Rstd_Std{i},'.k','linewidth',1);
%         set(gca,'Xtick',[hb{i}(1,1).XData]);
%         xticklabels(CompLab);
        
    else
        Ydata{i}    = X(ismember(c,i),XSel);
        Ydata{i}=Ydata{i} (:,ordcomp);
        SamNames{i}=strcat(MarkLabs(ismember(c,i)),SamLab(ismember(c,i)));
        h{i} = bar(Ydata{i}'); hold on;
        set(gca,'Xtick',[h{i}(1,1).XData]);
        xticklabels(CompLab);
        set(gca,'XTickLabelRotation',45);
        set(gca,'YTickLabelRotation',45);
        set(gca,'fontsize', 11);
%        legend(h{i},MarkLabs(ismember(c,i)),'autoupdate','off','Location','SouthEastOutside');
        legend(h{i},strcat(MarkLabs(ismember(c,i)),SamLab(ismember(c,i))),'autoupdate','off','Location','SouthEastOutside');
        
%         title(ExpLabs(b(i)));
        set(figH{i},'name',ExpLabs{b(i)});
    end
%     xlabel(Xlab); 
    ylabel(Ylab);
end

% 
% for i=1:numel(a)
%     figHB{i}=figure;
%      
%         hb{i} = bar(Rstd_Avg{i}'); hold on;
%         ctrb{i}   = bsxfun(@plus, 1:size(Rstd_Avg{i},2), [hb{i}.XOffset]');
%         errorbar(ctrb{i}, Rstd_Avg{i},Rstd_Std{i},'.k','linewidth',1);
%         set(gca,'Xtick',[hb{i}(1,1).XData]);
%         xticklabels(CompLab);
% end

%% PW Plot
% [a,b,c]=unique(DesignSel(:,strcmp(Colhead,NormSample{1})));
% MarkLabs=DesignSel(:,strcmp(Colhead,MarkCol));
% for i=1:numel(a)
%     figure; hold on;
%     for j=1:size(Ydata{i},1), plot(categorical(CompLab), Ydata{i}(j,:)','-o'); end
% end
% X-Y Plot
% if(XYplot && ~isempty(Xax) && isequal(numel(Yax),numel(Xax)))
%     MarkLabs=DesignSel(:,strcmp(Colhead,MarkCol));
%     [groupLab,indLeg,groupInd] = unique(MarkLabs);
%     MarkColor =num2cell(MyCM(length(indLeg)),2);
%     %     MarkColor =num2cell(jet(length(indLeg)),2);
%     MarkColor = cat(1,MarkColor{groupInd});
%
%     for i=1:numel(Xax)
%         figure,title([Xax{i}  ' vs ' Yax{i}]);
%         ph{1}=CHScatter(X(:,strcmp(AnalyteCol,Xax{i})),X(:,strcmp(AnalyteCol,Yax{i})),MarkLabs,[],MarkColor,ShowConvexHulls); hold on;
%         xlabel(Xax{i});ylabel(Yax{i});
%         set(ph{1},'marker','o','linestyle','none','MarkerSize',5);
%         legend(cat(1,ph{1}),groupLab{:})
%
%         if (any(ShowConvexHulls))
%             groupLabX = accumarray(groupInd,X(:,strcmp(AnalyteCol,Xax{i})),[length(groupLab) 1],@mean);
%             groupLabY = accumarray(groupInd,X(:,strcmp(AnalyteCol,Yax{i})),[length(groupLab) 1],@mean);
%             txHan = text(groupLabX,groupLabY,char(regexprep(groupLab,'^(.)','  $1')));
%             set(txHan,'horizontalalignment','center','clipping','on');
%             set(txHan,'fontsize',8);
%         end
%
%         if(ShowSampleNames)
%             [a,b,c] = unique(DesignSel(:,strcmp(Colhead,LabCol)));
%             for j=1:length(c), text(X(j,strcmp(AnalyteCol,Xax{i})),X(j,strcmp(AnalyteCol,Yax{i})),a{c(j)},'color','k','Fontsize',8,'FontWeight','bold','VerticalAlignment','top','HorizontalAlignment','center','clipping','on'); end
%         end
%     end
% end
%% figure Size Correction
 figHandles = get(groot, 'Children');
 axHandles = findall(figHandles,'type','axes','tag','');
%  set(figHandles, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 0.5]);
%  set(figHandles, 'Units', 'Normalized', 'OuterPosition', [0.5 0.5 1 1]);
 set(axHandles, 'Units', 'Normalized', 'OuterPosition', [-0.04 0 1.1 1]);
%% Reporting
 %       %Make a long list of some of the methods available in MS-Word
    %      Category='Selection'; % Category='ActiveDocument';
    %      PrintMethods(ActXWord,Category)
    
if(DoReporting)
    
    WordFileName = [datestr(now,1) '_Report.doc'];
    CurDir   = pwd;
    FileSpec = fullfile(CurDir,WordFileName);
    [ActXWord,WordHandle]=StartWord(FileSpec);
    fprintf('Document will be saved in %s\n',FileSpec);
    
    % Section Title
    
    if(isempty(ModelTitle)), ModelTitle  = {''}; end
    ActXWord.Selection.InsertBreak; %pagebreak
    TextString= ['Modeling Report - ' ModelTitle{1} ' - ' datestr(now)];
    WordText(ActXWord,TextString,'Heading 1',[0,1]);
    
    if(~isempty(ModelInfo))
        % Model Information
        ActXWord.Selection.InsertParagraph;
        TextString=['Model Information - ' ModelInfo{1}];
        WordText(ActXWord,TextString,'Heading 2',[0,1]);
    end
    if(~isempty(Tables))
        % Tables
        ActXWord.Selection.InsertParagraph;
        WordText(ActXWord,'Tables','Heading 2',[0,1]);
        for k=1:numel(Tables)
            [NoRows,NoCols]=size(Tables{k});
            WordCreateTable(ActXWord,NoRows,NoCols,Tables{k},1);
            TextString=['Table ' num2str(k) ': ' Titles{k}];
            WordText(ActXWord,TextString,'Normal',[0,1]);
        end
    end   
    % Figures
    ActXWord.Selection.InsertParagraph;
    WordText(ActXWord,'Figures','Heading 2',[0,1]);%enter after text
     figHandles = findobj('Type', 'figure');  % figHandles = get(groot, 'Children');
    for i=1:numel(figHandles)
        FigureIntoWord(ActXWord,figHandles(i));
        TextString= ['Figure ' num2str(i) ': Experiment ' figHandles(i).Name];
        WordText(ActXWord,TextString,'Normal',[0,1]);%enter after text
    end
    
    % Concluding Remarks
        ActXWord.Selection.InsertParagraph;
        TextString=['Concluding Remarks'];
        WordText(ActXWord,TextString,'Heading 2',[0,1]);
        
    % Add pagenumbers (0=not on first page)
    WordPageNumbers(ActXWord,'wdAlignPageNumberRight');

    % Updating Table of Content
    % Selection.GoTo What:=wdGoToField, Which:=wdGoToPrevious, Count:=1, Name:= "TOC"
    WordGoTo(ActXWord,7,3,1,'TOC',1);%%last 1 to delete the object
    WordCreateTOC(ActXWord,1,3);
    CloseWord(ActXWord,WordHandle,FileSpec);
end

function [actx_word,word_handle]=StartWord(word_file_p)
% Start an ActiveX session with Word:
actx_word = actxserver('Word.Application');
actx_word.Visible = true;
trace(actx_word.Visible);
if ~exist(word_file_p,'file')
    % Create new document:
    word_handle = invoke(actx_word.Documents,'Add');
    
    % Table of Content
    style='Heading 1';
    text='Table of Contents';
    WordText(actx_word,text,style,[1,1]);%enter before and after text
    WordCreateTOC(actx_word,1,3);
    %     actx_word.Selection.InsertBreak; %pagebreak
else
    % Open existing document:
    word_handle = invoke(actx_word.Documents,'Open',word_file_p);
    % Find end of document and make it the insertion point:
    end_of_doc = get(actx_word.ActiveDocument.Content,'end');
    set(actx_word.Selection,'Start',end_of_doc);
    set(actx_word.Selection,'End',end_of_doc);
end
end

function WordGoTo(actx_word_p,what_p,which_p,count_p,name_p,delete_p)
%Selection.GoTo(What,Which,Count,Name)
actx_word_p.Selection.GoTo(what_p,which_p,count_p,name_p);
if(delete_p)
    actx_word_p.Selection.Delete;
end

end

function WordCreateTOC(actx_word_p,upper_heading_p,lower_heading_p)
%      With ActiveDocument
%         .TablesOfContents.Add Range:=Selection.Range, RightAlignPageNumbers:= _
%             True, UseHeadingStyles:=True, UpperHeadingLevel:=1, _
%             LowerHeadingLevel:=3, IncludePageNumbers:=True, AddedStyles:="", _
%             UseHyperlinks:=True, HidePageNumbersInWeb:=True, UseOutlineLevels:= _
%             True
%         .TablesOfContents(1).TabLeader = wdTabLeaderDots
%         .TablesOfContents.Format = wdIndexIndent
%     End With
actx_word_p.ActiveDocument.TablesOfContents.Add(actx_word_p.Selection.Range,1,...
    upper_heading_p,lower_heading_p);

actx_word_p.Selection.TypeParagraph; %enter
end

function WordText(actx_word_p,text_p,style_p,enters_p,color_p)
%VB Macro
%Selection.TypeText Text:="Test!"
%in Matlab
%set(word.Selection,'Text','test');
%this also works
%word.Selection.TypeText('This is a test');
if(enters_p(1))
    actx_word_p.Selection.TypeParagraph; %enter
end
actx_word_p.Selection.Style = style_p;
if(nargin == 5)%check to see if color_p is defined
    actx_word_p.Selection.Font.ColorIndex=color_p;
end

if(isnumeric(text_p)), text_p=num2str(text_p); end % By KGP

actx_word_p.Selection.TypeText(text_p);
actx_word_p.Selection.Font.ColorIndex='wdAuto';%set back to default color
for k=1:enters_p(2)
    actx_word_p.Selection.TypeParagraph; %enter
end
end

function WordSymbol(actx_word_p,symbol_int_p)
% symbol_int_p holds an integer representing a symbol,
% the integer can be found in MSWord's insert->symbol window
% 176 = degree symbol
actx_word_p.Selection.InsertSymbol(symbol_int_p);
end

function WordCreateTable(actx_word_p,nr_rows_p,nr_cols_p,data_cell_p,enter_p)
%Add a table which auto fits cell's size to contents
if(enter_p(1))
    actx_word_p.Selection.TypeParagraph; %enter
end
%create the table
%Add = handle Add(handle, handle, int32, int32, Variant(Optional))
actx_word_p.ActiveDocument.Tables.Add(actx_word_p.Selection.Range,nr_rows_p,nr_cols_p,1,1);
%Hard-coded optionals
%first 1 same as DefaultTableBehavior:=wdWord9TableBehavior
%last  1 same as AutoFitBehavior:= wdAutoFitContent

%write the data into the table
for r=1:nr_rows_p
    for c=1:nr_cols_p
        %write data into current cell
        WordText(actx_word_p,data_cell_p{r,c},'Normal',[0,0]);
        
        if(r*c==nr_rows_p*nr_cols_p)
            %we are done, leave the table
            actx_word_p.Selection.MoveDown;
        else%move on to next cell
            actx_word_p.Selection.MoveRight;
        end
    end
end
end

function WordPageNumbers(actx_word_p,align_p)
%make sure the window isn't split
if (~strcmp(actx_word_p.ActiveWindow.View.SplitSpecial,'wdPaneNone'))
    actx_word_p.Panes(2).Close;
end
%make sure we are in printview
if (strcmp(actx_word_p.ActiveWindow.ActivePane.View.Type,'wdNormalView') | ...
        strcmp(actx_word_p.ActiveWindow.ActivePane.View.Type,'wdOutlineView'))
    actx_word_p.ActiveWindow.ActivePane.View.Type ='wdPrintView';
end
%view the headers-footers
actx_word_p.ActiveWindow.ActivePane.View.SeekView='wdSeekCurrentPageHeader';
if actx_word_p.Selection.HeaderFooter.IsHeader
    actx_word_p.ActiveWindow.ActivePane.View.SeekView='wdSeekCurrentPageFooter';
else
    actx_word_p.ActiveWindow.ActivePane.View.SeekView='wdSeekCurrentPageHeader';
end
%now add the pagenumbers 0->don't display any pagenumber on first page
actx_word_p.Selection.HeaderFooter.PageNumbers.Add(align_p,0);
actx_word_p.ActiveWindow.ActivePane.View.SeekView='wdSeekMainDocument';
end

function PrintMethods(actx_word_p,category_p)
style='Heading 3';
text=strcat(category_p,'-methods');
WordText(actx_word_p,text,style,[1,1]);

style='Normal';
text=strcat('Methods called from Matlab as: ActXWord.',category_p,'.MethodName(xxx)');
WordText(actx_word_p,text,style,[0,0]);
text='Ignore the first parameter "handle". ';
WordText(actx_word_p,text,style,[1,3]);

MethodsStruct=eval(['invoke(actx_word_p.' category_p ')']);
MethodsCell=struct2cell(MethodsStruct);
NrOfFcns=length(MethodsCell);
for i=1:NrOfFcns
    MethodString=MethodsCell{i};
    WordText(actx_word_p,MethodString,style,[0,1]);
end
end

function FigureIntoWord(actx_word_p,figHandles)
for i=1:numel(figHandles)
    % Capture current figure/model into clipboard:
    print('-dmeta',figHandles(i))
    fprintf('Figure i: %s\n',figHandles(i).Name);
    % Find end of document and make it the insertion point:
    end_of_doc = get(actx_word_p.activedocument.content,'end');
    set(actx_word_p.application.selection,'Start',end_of_doc);
    set(actx_word_p.application.selection,'End',end_of_doc);
    % Paste the contents of the Clipboard:
    %also works Paste(ActXWord.Selection)
    invoke(actx_word_p.Selection,'Paste');
    actx_word_p.Selection.TypeParagraph; %enter
end
end

function CloseWord(actx_word_p,word_handle_p,word_file_p)
if ~exist(word_file_p,'file')
    % Save file as new:
    invoke(word_handle_p,'SaveAs',word_file_p,1);
else
    % Save existing file:
    invoke(word_handle_p,'Save');
end
% Close the word window:
invoke(word_handle_p,'Close');
% Quit MS Word
invoke(actx_word_p,'Quit');
% Close Word and terminate ActiveX:
delete(actx_word_p);
end
