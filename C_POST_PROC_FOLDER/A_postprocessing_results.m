close all; clear, clc;
tic;
curr_fold = cd;
addpath([curr_fold(1:end-19) '\B_PROC_IDEN_FOLDER']);

single_postprocc = 0; % 1 plot diagramm from single id
n_id = 7; % number of identifications performed

%% SINGLE POSTPROCESSING
id_sing = 1; % single identification indicator

if single_postprocc == 1
    load(sprintf('results_identification_filt_%d.mat',id_sing));
    
    % STABILIZATION AND CLUSTERING GRAPH
    % stabilization
    figure
    plot(ID(ind1,1),ID(ind1,end-4),'x k');
    hold on
    plot(IDRED(:,1),IDRED(:,end-4),'r o');
    xlim([clusterA.f_min clusterA.f_up]);
    ylim([ssi.n_min ssi.n_max]);
    xlabel('Frequency [Hz]');
    ylabel('System order');
    legend ('Spurius Modes','Stable Modes');
    
    %clustering
    figure
    plot(dataF(:,1),dataF(:,2),'ro');
    hold on
    plot (CLUSTER_TOT(:,3),CLUSTER_TOT(:,4),'+b','markersize',13,'linewidth',2);
    hold on;
    plot (CLUSTER_TOT(:,3),CLUSTER_TOT(:,4),'oc','markersize',5,'linewidth',1.5,'MarkerFaceColor','c','MarkerEdgeColor','b');
    xlim([clusterA.f_min clusterA.f_up]);
    ylim([stab.d_min stab.d_max]);
    xlabel('Frequency [Hz]');
    ylabel('Damping ratio [-]');
    legend ('Candidate Modes','Identified Modes');
    
else
    
    %% AGGREGATE POSTPROCESSING
    for id = 1:n_id
        load(sprintf('results_identification_filt_%d.mat',id));
        n_channels = size(data,2)-1;
        FRE_ID{id} = ID_RES(:,3)';
        ZIT_ID{id} = ID_RES(:,4)';
        SHA_ID{id} = ID_RES(:,5:5+n_channels-1)';
    end
    FRE_IDmat = cell2mat(FRE_ID);
    ZIT_IDmat = cell2mat(ZIT_ID);
    SHA_IDmat = cell2mat(SHA_ID);
    
    [FRE_IDmatSort,insort] = sort(FRE_IDmat);
    ZIT_IDmatSort = ZIT_IDmat(insort);
    SHA_IDmatSort = SHA_IDmat(:,insort);
    SHA_IDmatSort = SHA_IDmatSort./max(abs(SHA_IDmatSort));
    for ii = 1:size(SHA_IDmatSort,2)
        if min(SHA_IDmatSort(:,ii)) == -1
            SHA_IDmatSort(:,ii) = -SHA_IDmatSort(:,ii);
        end
    end
    MAC_MATRIX = compute_mac(SHA_IDmatSort,SHA_IDmatSort);
    
    %% MANUAL CLASSIFICATION OF MODES BASED ON FREQUENCY AND MAC
    % This cell "class_mode" must be filled manually since the number of class
    % are unknown at priori. Class with just 1 mode are neglected for safety
    mm = 1;
    class_mode{mm} = [1:5];
    %
    mm = mm+1;
    class_mode{mm} = [6:11];
    % OLD     SHA_IDmatSort(:,9) = -SHA_IDmatSort(:,9);
    % OLD     SHA_IDmatSort(:,11) = -SHA_IDmatSort(:,11);
    %
    mm = mm+1;
    class_mode{mm} = [12:18];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [19];
    %
    mm = mm+1;
    class_mode{mm} = [20:25];
    %
    mm = mm+1;
    class_mode{mm} = [26:36];
    %
    mm = mm+1;
    class_mode{mm} = [37:40];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [41];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [42];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [43];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [44];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [45];
    %
    mm = mm+1;
    class_mode{mm} = [46:51];
    %
    mm = mm+1;
    class_mode{mm} = [52:59];
    %
    mm = mm+1;
    class_mode{mm} = [60:69];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [70];
    %
    mm = mm+1;
    class_mode{mm} = [71:85];
    %
    mm = mm+1;
    class_mode{mm} = [86:87];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [88];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [89];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [90];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [91];
    %
    mm = mm+1;
    class_mode{mm} = [92:93];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [94];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [95];
    %
    mm = mm+1;
    class_mode{mm} = [96:98];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [99];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [100];
    %
    mm = mm+1;
    class_mode{mm} = [101:103];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [104];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [105];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [106];
    %
    mm = mm+1;
    class_mode{mm} = [107:109];
    %
    %     mm = mm+1;
    %     class_mode{mm} = [110];
    
    %% AUTOMATIC AVERAGING
    fnID = [];
    zitID = [];
    shaID = [];
    fnIDstd = [];
    zitIDstd = [];
    shaIDstd = [];
    for ii = 1:size(class_mode,2)
        fnID = [ fnID , nanmean(FRE_IDmatSort(class_mode{ii})) ];
        zitID = [ zitID , nanmean(ZIT_IDmatSort(class_mode{ii})) ];
        %
        fnIDstd = [ fnIDstd , nanstd(FRE_IDmatSort(class_mode{ii})) ];
        zitIDstd = [ zitIDstd , nanstd(ZIT_IDmatSort(class_mode{ii})) ];
        
        if size(class_mode{ii},2) == 1
            shaADJ = SHA_IDmatSort(:,class_mode{ii});
            shaADJstd = 0*SHA_IDmatSort(:,class_mode{ii});
        else
            for jj = 1:size(SHA_IDmatSort,1)
                [fout0,xout]=ksdensity(SHA_IDmatSort(jj,class_mode{ii}));
                [val,ind] = max(fout0);
                xsel = xout(ind);
                xsel = xsel(1);
                shaADJ(jj,1) = xsel;
                shaADJstd(jj,1) = nanstd(SHA_IDmatSort(jj,class_mode{ii}));
            end
        end
        
        shaID = [ shaID , shaADJ ];
        shaIDstd = [ shaIDstd , shaADJstd ];
        shaADJ = [];
        shaADJstd = [];
    end
    
    shaID = shaID./ max(abs(shaID));
    for ii = 1:size(shaID,2)
        if min(shaID(:,ii)) == -1
            shaID(:,ii) = -shaID(:,ii);
        end
    end
    MAC_MATRIXID = compute_mac(shaID,shaID);
    
    ID_RESULTS = [fnID ; zitID ; shaID];
    ID_RESULTS_STD = [fnIDstd ; zitIDstd ; shaIDstd];
    
%                 save([cd '\ID_RESULTS'],'ID_RESULTS','ID_RESULTS_STD');
    
    %% AGGREGATE STABILIZATION AND CLUSTERING GRAPH
    for id = 1:n_id
        load(sprintf('results_identification_filt_%d.mat',id));
        if id ==1
            figure
        end
        
        % stabilization
        subplot(1,2,1)
        hold on
        plot(ID(ind1,1),ID(ind1,end-4),'x k');
        hold on
        plot(IDRED(:,1),IDRED(:,end-4),'r o');
        if id ==n_id
            xlim([clusterA.f_min clusterA.f_up]);
            ylim([ssi.n_min ssi.n_max]);
            xlabel('Frequency [Hz]');
            ylabel('System order');
            legend ('Spurius Modes','Stable Modes');
        end
        
        %clustering
        subplot(1,2,2)
        hold on
        if  id ~= 1
            plot(dataF(:,1),dataF(:,2),'ro','HandleVisibility','off');
        else
            plot(dataF(:,1),dataF(:,2),'ro','HandleVisibility','on');
        end
    end
    subplot(1,2,2)
    hold on
    plot (ID_RESULTS(1,:),ID_RESULTS(2,:),'+b','markersize',13,'linewidth',2,'HandleVisibility','on');
    hold on;
    plot (ID_RESULTS(1,:),ID_RESULTS(2,:),'oc','markersize',5,'linewidth',1.5,'MarkerFaceColor','c','MarkerEdgeColor','b','HandleVisibility','off');
    xlim([clusterA.f_min clusterA.f_up]);
    ylim([stab.d_min stab.d_max]);
    xlabel('Frequency [Hz]');
    ylabel('Damping ratio [-]');
    legend ('Candidate Modes','Identified Modes');
    
end

toc;
