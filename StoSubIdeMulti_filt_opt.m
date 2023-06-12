close all; clear; clc;
curr_fold = cd;
%addpath([curr_fold(1:end-19) '\A_PREP_PROC_FOLDER']);
%load("data_filt_2022_12_20_1.mat");

addpath("C:\Users\rebby\Desktop\TESI\PreprocessedData\Test_2022_12_20");
load("data_filt_2022_12_20_6.mat");

win = 60000; % # of samples of the window: 60000/fs / 60 = 10 min (600 sec) --> da ripetere 7 volte per coprire l'intero segnale
for ide_singola = 2:7
    clearvars -except ide_singola win data curr_fold
    %addpath([curr_fold '\functions']);
    addpath(fullfile(curr_fold, 'functions'));
    tic;
    
    % PARAMETRI DATI
    RR = 1; % desampling factor used in preprocessing
    fs = 100 / RR; % sample frequency [Hz]
    start_idx = ceil(win * (ide_singola - 1) + 1);
    end_idx = min(floor(win * ide_singola), size(data, 1));
    signal = data(start_idx:end_idx, 2:end); % acc data (n.samples x n.channels)
    ch_num = size(signal, 2); % number of channels
    
    % PARAMETRI SSI
    ssi.n_min=20;      % Ordine minimo dell'algoritmo SSI - 
    ssi.n_max=200;      % Ordine massimo dell'algoritmo SSI - 
    ssi.n_step=2;      % Variazione dell'ordine del sistema
    ssi.n_ord=ssi.n_min:ssi.n_step:ssi.n_max; % Numero di ordini analizzati
    
    % PARAMETRI STABILIZATION DIAGRAM
    stab.d_min = 0.5/100;      % Smorzamento minimo del modo del xx%
    stab.d_max = 6/100;    % Smorzamanto massimo dei modo del xx%
    stab.de_f  = 0.025;  % Scarto in frequenza minore di xx Hz
    stab.de_z  = 20/100;    % Scarto in smorzamento minore del xx%
    stab.MAC   = 0.98;   % Valore del MAC
    stab.p_min = 0.25;   % Soglia minima Indice di Occorrenza
    stab.Nmin_stab = ceil(stab.p_min*(ssi.n_max-ssi.n_min)/ssi.n_step); % Numero di modi di soglia indice occorrenza
    
    % PARAMETRI CLUSTER ANALYSIS
    clusterA.wF=1; % Peso da dare alla Frequenza
    clusterA.wMAC=1; % Peso da dare al MAC
    clusterA.dif_minF=0.08; % Differenza massima tra due frequenze [Hz]
    clusterA.dif_minMAC=0.95; % Percentuale minima di MAC tra due modi uguali
    clusterA.linkage='average'; % (average,centroid,single,median,ward,weighted)
    clusterA.f_min = 1; % Hz definisce la zona dove applicare il cluster (comunemente è la min cut-off freq di un filtro)
    clusterA.f_up = 10; % Hz definisce la zona dove applicare il cluster (comunemente è la max cut-off freq di un filtro)
    
    %% IDENTIFICATION OF DIFFERENT ORDERS
    % Preallocate ID to optimize code execution
    ID = cell(length(ssi.n_min:ssi.n_step:ssi.n_max), 1);
    count = 1;
    ssi.n_min
    ssi.n_step
    ssi.n_max
    for i2=ssi.n_min:ssi.n_step:ssi.n_max
        n_block=2*ssi.n_max;  % n_block Hankel matrix
        % {ID} colonne  = ordine  del sistema
        % [ID] righe = modi estratti
        % [ID] colonne  = frequenza, smorzamento, forma modale
        ID{count,1}=SSI_opt(signal,n_block,fs,i2);
        count=count+1;
    end
    
    % Results
    % ID
    
    %% STABILIZATION
    %deleting modes with damping > of d_max
    for ii=1:size(ID,1)
        pos= ID{ii}(:,2)>stab.d_max;
        ID{ii}(pos,:)=[];
    end
    %deleting modes with damping < of d_min
    for ii=1:size(ID,1)
        pos= ID{ii}(:,2)<stab.d_min;
        ID{ii}(pos,:)=[];
    end
    
    % Application of stabilization criteria:
    flag=0;
    cont_empty=0;
    IDENT_AUG=ID;
    for i1=1:size(IDENT_AUG,1)
        if isempty(IDENT_AUG{i1})==0
            if flag==0
                nc=size(IDENT_AUG{i1},2);
                IDENT_AUG{i1}(:,nc+1)=inf;
                IDENT_AUG{i1}(:,nc+2)=inf;
                IDENT_AUG{i1}(:,nc+3)=inf;
                flag=1;
            else
                if cont_empty==0
                    for i2=1:size(IDENT_AUG{i1},1)
                        i3_min=1;
                        de_f_min=inf;
                        f_i2=IDENT_AUG{i1}(i2,1);
                        z_i2=IDENT_AUG{i1}(i2,2);
                        fi_i2=IDENT_AUG{i1}(i2,3:nc)';
                        for i3=1:size(IDENT_AUG{i1-1},1)
                            f_i3=IDENT_AUG{i1-1}(i3,1);
                            if f_i2==f_i3
                                stab.de_f_i23=0;
                            elseif f_i2>f_i3
                                stab.de_f_i23=(f_i2-f_i3);
                            else
                                stab.de_f_i23=(f_i3-f_i2);
                            end
                            if stab.de_f_i23<de_f_min, de_f_min=stab.de_f_i23;
                                i3_min=i3;
                            end
                        end
                        stab.de_f_i23=de_f_min;
                        z_i3=IDENT_AUG{i1-1}(i3_min,2);
                        fi_i3=IDENT_AUG{i1-1}(i3_min,3:nc)';
                        if z_i2==z_i3
                            stab.de_z_i23=0;
                        elseif z_i2>z_i3
                            stab.de_z_i23=(z_i2-z_i3)/z_i3;
                        else
                            stab.de_z_i23=(z_i3-z_i2)/z_i2;
                        end
                        stab.MAC_i23=(abs(fi_i2'*fi_i3))^2/(fi_i2'*fi_i2)/(fi_i3'*fi_i3);
                        
                        IDENT_AUG{i1}(i2,nc+1)=stab.de_f_i23;
                        IDENT_AUG{i1}(i2,nc+2)=stab.de_z_i23;
                        IDENT_AUG{i1}(i2,nc+3)=stab.MAC_i23;
                    end
                else
                    IDENT_AUG{i1}(:,nc+1)=inf;
                    IDENT_AUG{i1}(:,nc+2)=inf;
                    IDENT_AUG{i1}(:,nc+3)=inf;
                end
            end
            cont_empty=0;
        else
            cont_empty=cont_empty+1;
        end
    end
    
    % conditions of stabilization
    ID = zeros(sum(cellfun(@(x) size(x, 1), IDENT_AUG)), nc + 8);
    app = zeros(1, nc + 8);

    idx = 1;

    for i1 = 1:size(IDENT_AUG, 1)
        if ~isempty(IDENT_AUG{i1})
            for i2 = 1:size(IDENT_AUG{i1}, 1)
                app(1:nc) = IDENT_AUG{i1}(i2, :);
                app(nc + 1) = ssi.n_ord(i1);

                % STABILE PER TUTTE E TRE LE CONDIZIONI
                if (app(nc + 1) > stab.de_f) || ...
                   (app(nc + 2) > stab.de_z) || ...
                   (app(nc + 3) < stab.MAC)
                    app(nc + 5) = 0;
                else
                    app(nc + 5) = 1;
                end

                % STABILITA PER FREQUENZA
                if (app(nc + 1) > stab.de_f)
                    app(nc + 6) = 0;
                else
                    app(nc + 6) = 1;
                end

                % STABILITA PER SMORZAMENTO
                if (app(nc + 2) > stab.de_z)
                    app(nc + 7) = 0;
                else
                    app(nc + 7) = 1;
                end

                % STABILITA PER MAC
                if (app(nc + 3) < stab.MAC)
                    app(nc + 8) = 0;
                else
                    app(nc + 8) = 1;
                end

                ID(idx, :) = app;
                idx = idx + 1;
            end
        end
    end

    [row] = find(ID(:,end-3)>0);
    [ins] = find(ID(:,end-3)<1);
    [spurius]=find(ID(:,end)<1 & ID(:,end-1)<1 & ID(:,end-2)<1);
    ind1=[1:size(ID,1)];
    ind1(row)=[];
    IDRED=ID(row,:); % modi stabili
    INRED=ID(ins,:); % modi INstabili
    [staF]=find(INRED(:,end-2)>0);
    [staD]=find(INRED(:,end-1)>0);
    [staMac]=find(INRED(:,end)>0);
    staF=INRED(staF,:);
    staD=INRED(staD,:);
    staMac=INRED(staMac,:);
    spurius=ID(spurius,:);
    
    full_id=IDRED(:,:);
    CLUST =IDRED(:,1:end-8);
    
    % Results
    % ind1
    % ID
    % IDRED
    % full_id
    % CLUST
    
    
    %% CLUSTERING PRETRATTAMENTO
    indices=find(CLUST(:,1)<clusterA.f_up & CLUST(:,1)>clusterA.f_min);   % da impostrare in base a dove voglio applicare la FCM clustering
    mod=size(signal,2)+2;
    dataX=CLUST(indices,1:mod);
    dataF=dataX;  %da cluster rinomino dataF
    IDord = sortrows(dataF,1); %Ordino in modo crescente le frequenze
    nID = ones(size(IDord,1),1);%Aggiungo colonna 1 che servirà forse successivamente per il peso
    IDord=[nID';IDord']'; %|1|frequenze|smorz|forme modali|
    clear nID
    
    IDordMAC = compute_mac(IDord(:, 4:end)', IDord(:, 4:end)');
    IDordMAC_D = ones(size(IDordMAC)) - IDordMAC; % Matrice MAC dissimilità

    % AGGIUNGO RIGA PER AUMENTARE LA DIFFERENZA DI MAC
    posNO = IDordMAC_D > (1 - clusterA.dif_minMAC); % ACCETTO MAC 90%
    IDordMAC_D(posNO) = 1;

    clear posNO;
    
    % Calcolo differenza (valori accettabili tra 0 e 0.05)
    ID_ordFREQ=compute_freq2(IDord(:,2),IDord(:,2));
    posNO=find(ID_ordFREQ>clusterA.dif_minF);
    ID_ordFREQ(posNO)=1;
    
    % ID_ordDAMP=compute_damp(IDord(:,3),IDord(:,3));
    % (Valori accettabili tra 0 e 0,5)
    % CALCOLO DA MATRICE DISSIMILITà TOTALE (alfa*Freq + beta*MAC)
    ID_Dissim=(clusterA.wMAC.*IDordMAC_D)+(clusterA.wF.*ID_ordFREQ); %Matrice distanza (DISSIMILITà)
    dist_min=(clusterA.wMAC*(1-clusterA.dif_minMAC))+(clusterA.wF*clusterA.dif_minF);
    clear IDordMAC MAC_one IDordMAC_D ID_ordFREQ ID_ordDAMP
    
    % Agglomerative hierarchical cluster tree
    VECTOR_DISSIM = squareform (ID_Dissim);
    Z = linkage(VECTOR_DISSIM,clusterA.linkage);
    % linkage(X,method) creates the tree using the specified method, where method describes how to measure the distance between clusters.
    CLUSTZ = cluster(Z,'cutoff',dist_min,'criterion','distance');
    % cluster(Z,'cutoff',c,'criterion',criterion) uses the specified criterion for forming clusters, where criterion is one of the strings 'inconsistent' (default) or 'distance'. The 'distance' criterion uses the distance between the two subnodes merged at a node to measure node height. All leaves at or below a node with height less than c are grouped into a cluster.
    % DISTANCE: average centroid complete median single ward weighted
    % c1 = cluster(Z,'maxclust',n)  %constructs a maximum of n clusters using the 'distance' criterion. cluster finds the smallest height at which a horizontal cut through the tree leaves n or fewer clusters.
    H=dendrogram(Z,0); % H=dendrogram(Z,0,'ColorThreshold','default');
    %RIORDINA IL NUMERO DEL CLUSTER da 1 a n
    clear i2
    contat=1;
    CLUSTZ(1,2)=contat;
    for i2=2:size(CLUSTZ,1)
        pos_col=find(CLUSTZ(i2,1)==CLUSTZ(1:i2-1,1));
        [m,n]=size(pos_col);
        if m>0.9 && n>0.9
            CLUSTZ(i2,2)= CLUSTZ(pos_col(1),2);
        else
            contat=contat+1;
            CLUSTZ(i2,2)=  contat;
        end
    end
    clear contat pos_col
    CLUSTZ(:,1)=[]; %Elimino la numerazione che mi restituiva il linkage
    CLUSTZ=[CLUSTZ';IDord']'; %|n°CLUST|occor 1|freq|damp|forma mod|
    zero=zeros(size(CLUSTZ,1),1);
    CLUSTZ=[CLUSTZ';zero']';
    clear zero
    
    %% TROVO LE FREQUENZE DI RIFERIMENTO CLUST
    clust_i=[]; %CLUST iesimo che si analizza e di cui si fa la media
    clust_i_1=[]; %CLUST iesimo dove son presenti frequenze con le medesime forme modali (1;1;1;1;...;1)
    
    clust_i_sc=[];%CLUST iesimo dove son presenti frequenze scartate (1;1;1;1;...;1)
    num_cluster_fin=1; %numerazione cluster finale
    num_cluster_sc=1;
    
    q_clustz=max(CLUSTZ(:,1)); %numero di cluster trovati con il linkage
    % Determine the maximum number of rows for CLUSTER_TOT
    max_rows = q_clustz;

    % Preallocate CLUSTER_TOT with zeros
    CLUSTER_TOT = zeros(max_rows, ch_num + 8); %Matrice dove saranno contenute tutte le frequenze finali clusterizzate
    CLUSTER_SCART = zeros(max_rows_scart, ch_num + 4);%Matrice dove saranno contenute tutte le frequenze scartate

    for n_cluster = 1:q_clustz
        q = find(CLUSTZ(:,1) == n_cluster); % numero di frequenze presenti nel cluster

        % CASO 1) Se sono presenti più di NminSTAB stabilito all'inzio procedo con la valutazione
        if size(q,1) >= stab.Nmin_stab
            clust_i = CLUSTZ(q,:);

            % Trovo posizione forma modale dalla 5 alla 12 colonna nel clust_i
            pos_fin = 0; % posizione degli 1 e -1
            for ic = 5:size(clust_i,2)-1
                pos_1 = find(clust_i(:,ic) == 1 | clust_i(:,ic) == -1);
                if length(pos_1) >= 1 && length(pos_1) > pos_fin
                    clust_i_1(1:size(pos_1,1),1:size(clust_i,2)) = clust_i(pos_1,:);
                    pos_fin = size(clust_i_1,1);
                end
            end
            clear pos_fin pos_1 q ic

            % Metto 1 nella 13ima colonna per indicare quali frequenze verranno prese per il cluster
            % (Salvo le frequenze scartate)
            for isc = 1:size(clust_i_1,1)
                for iY = 1:size(CLUSTZ,1)
                    if CLUSTZ(iY,3:ch_num+4) == clust_i_1(isc,3:ch_num+4)
                        CLUSTZ(iY,ch_num+5) = 1;
                    end
                end
            end
            clear isc iY

            % TROVO LA FREQUENZA DI RIFERIMENTO COME MEDIA - TIPO A
            OCCUR_i = sum(clust_i_1(:,2));
            OCCUR_iP = (sum(clust_i_1(:,2))) / ((ssi.n_max - ssi.n_min) / ssi.n_step);
            media_fdm = mean(clust_i_1(:,3:ch_num+4), 1);
            ScartF = std(clust_i_1(:,3));
            ScartD = std(clust_i_1(:,4));

            % Calcolo il MAC MEDIO
            MAC_i = compute_mac(clust_i_1(:,5:ch_num+4)', clust_i_1(:,5:ch_num+4)');
            MAC_iDiag = triu(MAC_i, +1);
            k = find(MAC_iDiag == 0); % find: cerca l'indice per cui x è nullo
            MAC_iDiag(k) = NaN;
            mediaMAC = mean(MAC_iDiag, 'omitnan');
            mediaMAC = mean(mediaMAC, 'omitnan');

            % CLUSTER_TOT + STATISTIC
            % |n°|occorrenza|freq|smorz|forma modale|min|max|scart f|scart smo|MACmedio
            CLUSTER_TOT(num_cluster_fin, 1) = num_cluster_fin(1, 1);
            CLUSTER_TOT(num_cluster_fin, 2) = OCCUR_i(1, 1);
            CLUSTER_TOT(num_cluster_fin, 3:ch_num+4) = media_fdm(1, 1:ch_num+2);
            CLUSTER_TOT(num_cluster_fin, ch_num+5) = OCCUR_iP(1, 1);
            CLUSTER_TOT(num_cluster_fin, ch_num+6) = ScartF(1, 1);
            CLUSTER_TOT(num_cluster_fin, ch_num+7) = ScartD(1, 1);
            CLUSTER_TOT(num_cluster_fin, ch_num+8) = mediaMAC(1, 1);

            num_cluster_fin = num_cluster_fin + 1;
            clear OCCUR_i media_fdm ScartF ScartD MAC_i MAC_iDiag k mediaMAC

            % TROVO LA FREQUENZA DI RIFERIMENTO CON LA Kmeans - TIPO B (SCARTATA!)
            % OCCUR_i = sum(clust_i_1(:,2));
            % [idx,C,sumd,D] = kmeans(clust_i_1(:,3:4),1,'Distance','cityblock');
            % maxU_i = min(D); % distanza minore frequenza di riferimento
            % posmaxU_i = find(D == maxU_i); % posizione della frequenza di riferimento
            % CLUSTER_TOT(num_cluster_fin,3:ch_num+4) = clust_i_1(posmaxU_i(1),3:ch_num+4);
            % CLUSTER_TOT(num_cluster_fin,2) = OCCUR_i(1,1);
            % CLUSTER_TOT(num_cluster_fin,1) = num_cluster_fin(1,1);
            %
            % num_cluster_fin = num_cluster_fin + 1;
            % clear OCCUR_i idx C sumd D maxU_i posmaxU_i
        elseif size(q,1) < stab.Nmin_stab && size(q,1) > 1
            clust_i_sc = CLUSTZ(q,:);
            OCCUR_i = sum(clust_i_sc(:,2));
            media_fdm = mean(clust_i_sc(:,3:ch_num+4));
            CLUSTER_SCART(num_cluster_sc,3:ch_num+4) = media_fdm(1,1:ch_num+2);
            CLUSTER_SCART(num_cluster_sc,2) = OCCUR_i(1,1);
            CLUSTER_SCART(num_cluster_sc,1) = num_cluster_sc(1,1);
            num_cluster_sc = num_cluster_sc + 1;
            clear OCCUR_i media_fdm
        elseif size(q,1) == 1
            clust_i_sc = CLUSTZ(q,:);
            OCCUR_i = sum(clust_i_sc(:,2));
            media_fdm = clust_i_sc(:,3:ch_num+4);
            CLUSTER_SCART(num_cluster_sc,3:ch_num+4) = media_fdm(1,1:ch_num+2);
            CLUSTER_SCART(num_cluster_sc,2) = OCCUR_i(1,1);
            CLUSTER_SCART(num_cluster_sc,1) = num_cluster_sc(1,1);
            num_cluster_sc = num_cluster_sc + 1;
            clear OCCUR_i media_fdm
        end
        clear q clust_i clust_i_1 clust_i_sc
    end

    % CLUSTER_TOT(num_cluster_fin:end, :) = []; % Remove excess rows from CLUSTER_TOT
    CLUSTER_SCART(num_cluster_sc:end, :) = []; % Remove excess rows from CLUSTER_SCART

    pos_outliers=find(CLUSTZ(1:size(CLUSTZ,1),ch_num+5)==0);
    outliers(:,1:2)=CLUSTZ(pos_outliers,3:4);
    pos_candidate=find(CLUSTZ(1:size(CLUSTZ,1),ch_num+5)==1);
    candidate(:,1:2)=CLUSTZ(pos_candidate,3:4);
    
    CLUSTZ1=CLUSTZ;
    pos_outliers1=find(CLUSTZ1(1:size(CLUSTZ,1),ch_num+5)==0);
    CLUSTZ1(pos_outliers1,:)=[];
    n_clust_tot=max(CLUSTZ1(:,1));
    col=1;
    for i_nct=1:n_clust_tot
        n_freq_clust=find(CLUSTZ1(:,1)==i_nct);
        if size(n_freq_clust)< stab.Nmin_stab
            CLUSTZ1(n_freq_clust,:)=[];
        else
            CLUSTZ1(n_freq_clust,ch_num+6)=col;
            col=col+1;
        end
    end
    
    nnn=[ssi.n_min;ssi.n_max];
    
    ID_RES=CLUSTER_TOT(:,:);
    
    % Results
    % CLUSTER_TOT
    % ID_RES
    
    %% SAVING DATA
    fprintf('# %d Identificazioni terminate!...Stampa risultati in corso... \n',ide_singola)
    
    save(sprintf('results_identification_filt_%d',ide_singola),'-v7.3');
    
    toc;
    
end

%% STABILIZATION AND CLUSTERING GRAPH
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

% Specifica il percorso completo della cartella di destinazione
folder = "C:\Users\rebby\Desktop\TESI\IdentificatedDataFigures\data_filt_2022_12_20";

% Specifica il nome del file e il formato desiderato (ad esempio 'figura.png')
filename = 'Test_20.12.2023_6.png';

% Crea il percorso completo del file
filepath = fullfile(folder, filename);

% Salva la figura nel percorso specificato
saveas(gcf, filepath);
%%
toc;
