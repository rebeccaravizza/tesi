close all; clear; clc;
tic;
load ID_RESULTS.mat;

savefig = 0;
lambda = 1;
modeSTART = 1;
modeEND = size(ID_RESULTS,2);

% NB: THE COORDINATES REPORTED HEREINAFTER ARE ONLY VALID FOR PLOT
% (GRAPHICAL SCOPE) DO NOT USE THEM FOR MODEL UPDATING !!!!!!!!!!!
% pointscoord = [point #	x	y	z	DoF]
sA = 4;
sa = 1;
quota = 0;
pointscoord = [1	0	sa	quota	2 % Instrumented DoFs -v
    1	0	sa	quota	1
    2	-sA	0	quota	2
    
    3	0	-sa	quota	1 % Not Instrumented DoFs -v
    3	0	-sa	quota	2
    3	0	-sa	quota	3
    4	sA	0	quota	1
    4	sA	0	quota	2
    4	sA	0	quota	3];

%% DEFINING GEOMETRY
% Instrumented DoFs
shaIn = ID_RESULTS(3:end,1);
gdlIn = (1:size(ID_RESULTS,1)-2)';
XIn = pointscoord(1:3,2);
YIn = pointscoord(1:3,3);
ZIn = pointscoord(1:3,4);
diIn = pointscoord(1:3,5);

% Not Instrumented DoFs
shaNI = zeros( length(pointscoord(4:end,2)) ,1);
gdlNI = (gdlIn(end)+1:gdlIn(end)+size(shaNI,1))';
XNI = pointscoord(4:end,2);
YNI = pointscoord(4:end,3);
ZNI = pointscoord(4:end,4);
diNI = pointscoord(4:end,5);

% Total DoFs
gdl = [gdlIn;gdlNI];
X = [XIn;XNI];
Y = [YIn;YNI];
Z = [ZIn;ZNI];
di = [diIn;diNI];

Connecting_lines = [1 2 % point i , point j
    2 3
    3 4
    4 1];

coorduniqueDEFxIn = FindDefShape(gdlIn,XIn,YIn,ZIn,diIn,shaIn,lambda);
NumInstrPoints = length(coorduniqueDEFxIn); clearvars coorduniqueDEFxIn;

%% PLOTTING SHAPES
for mode = modeSTART:modeEND

    % Instrumented DoFs
    shaIn = ID_RESULTS(3:end,mode);

    % Total DoFs
    sha = [shaIn;shaNI];

    [coorduniqueDEFx , coorduniqueDEFy , coorduniqueDEFz , coorduniquex , coorduniquey , coorduniquez , relpointgld] = FindDefShape(gdl,X,Y,Z,di,sha,lambda);

    % PLOT
    fA = figure;
    %

    for ii = 1:size(Connecting_lines,1)
        subplot(2,2,1)
        hold on
        plot3([coorduniquex(Connecting_lines(ii,1));coorduniquex(Connecting_lines(ii,2))] , [coorduniquey(Connecting_lines(ii,1));coorduniquey(Connecting_lines(ii,2))] , [coorduniquez(Connecting_lines(ii,1));coorduniquez(Connecting_lines(ii,2))],'.--','Color',[0.5 0.5 0.5],'LineWidth',2);
        hold on
        plot3([coorduniqueDEFx(Connecting_lines(ii,1));coorduniqueDEFx(Connecting_lines(ii,2))] , [coorduniqueDEFy(Connecting_lines(ii,1));coorduniqueDEFy(Connecting_lines(ii,2))] , [coorduniqueDEFz(Connecting_lines(ii,1));coorduniqueDEFz(Connecting_lines(ii,2))],'.-','Color','k','LineWidth',2);
        %
        subplot(2,2,2)
        hold on
        plot3([coorduniquex(Connecting_lines(ii,1));coorduniquex(Connecting_lines(ii,2))] , [coorduniquey(Connecting_lines(ii,1));coorduniquey(Connecting_lines(ii,2))] , [coorduniquez(Connecting_lines(ii,1));coorduniquez(Connecting_lines(ii,2))],'.--','Color',[0.5 0.5 0.5],'LineWidth',2);
        hold on
        plot3([coorduniqueDEFx(Connecting_lines(ii,1));coorduniqueDEFx(Connecting_lines(ii,2))] , [coorduniqueDEFy(Connecting_lines(ii,1));coorduniqueDEFy(Connecting_lines(ii,2))] , [coorduniqueDEFz(Connecting_lines(ii,1));coorduniqueDEFz(Connecting_lines(ii,2))],'.-','Color','k','LineWidth',2);
        %
        subplot(2,2,3)
        hold on
        plot3([coorduniquex(Connecting_lines(ii,1));coorduniquex(Connecting_lines(ii,2))] , [coorduniquey(Connecting_lines(ii,1));coorduniquey(Connecting_lines(ii,2))] , [coorduniquez(Connecting_lines(ii,1));coorduniquez(Connecting_lines(ii,2))],'.--','Color',[0.5 0.5 0.5],'LineWidth',2);
        hold on
        plot3([coorduniqueDEFx(Connecting_lines(ii,1));coorduniqueDEFx(Connecting_lines(ii,2))] , [coorduniqueDEFy(Connecting_lines(ii,1));coorduniqueDEFy(Connecting_lines(ii,2))] , [coorduniqueDEFz(Connecting_lines(ii,1));coorduniqueDEFz(Connecting_lines(ii,2))],'.-','Color','k','LineWidth',2);
        %
        subplot(2,2,4)
        hold on
        plot3([coorduniquex(Connecting_lines(ii,1));coorduniquex(Connecting_lines(ii,2))] , [coorduniquey(Connecting_lines(ii,1));coorduniquey(Connecting_lines(ii,2))] , [coorduniquez(Connecting_lines(ii,1));coorduniquez(Connecting_lines(ii,2))],'.--','Color',[0.5 0.5 0.5],'LineWidth',2);
        hold on
        plot3([coorduniqueDEFx(Connecting_lines(ii,1));coorduniqueDEFx(Connecting_lines(ii,2))] , [coorduniqueDEFy(Connecting_lines(ii,1));coorduniqueDEFy(Connecting_lines(ii,2))] , [coorduniqueDEFz(Connecting_lines(ii,1));coorduniqueDEFz(Connecting_lines(ii,2))],'.-','Color','k','LineWidth',2);
    end
    subplot(2,2,1)
    hold on
    plot3(coorduniqueDEFx(1:NumInstrPoints),coorduniqueDEFy(1:NumInstrPoints),coorduniqueDEFz(1:NumInstrPoints),'o','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
    hold on
    plot3([coorduniquex,coorduniqueDEFx]', [coorduniquey,coorduniqueDEFy]' , [coorduniquez,coorduniqueDEFz]','-.','Color',[0.75 0.75 0.75],'LineWidth',1);
    title(sprintf('mode: %d - freq.: %.2f +-%.3f [Hz] \ndamp.: %.2f +-%.3f [%%]',mode,ID_RESULTS(1,mode),ID_RESULTS_STD(1,mode),100*ID_RESULTS(2,mode),100*ID_RESULTS_STD(2,mode)));
    grid on
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    pbaspect([max(X)-min(X)+1,max(Y)-min(Y)+1,max(Z)-min(Z)+1]);
    view([-45 45]);
    %
    subplot(2,2,2)
    hold on
    plot3(coorduniqueDEFx(1:NumInstrPoints),coorduniqueDEFy(1:NumInstrPoints),coorduniqueDEFz(1:NumInstrPoints),'o','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
    hold on
    plot3([coorduniquex,coorduniqueDEFx]', [coorduniquey,coorduniqueDEFy]' , [coorduniquez,coorduniqueDEFz]','-.','Color',[0.75 0.75 0.75],'LineWidth',1);
    grid on
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    pbaspect([max(X)-min(X)+1,max(Y)-min(Y)+1,max(Z)-min(Z)+1]);
    view([0 90]);
    %
    subplot(2,2,3)
    hold on
    plot3(coorduniqueDEFx(1:NumInstrPoints),coorduniqueDEFy(1:NumInstrPoints),coorduniqueDEFz(1:NumInstrPoints),'o','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
    hold on
    plot3([coorduniquex,coorduniqueDEFx]', [coorduniquey,coorduniqueDEFy]' , [coorduniquez,coorduniqueDEFz]','-.','Color',[0.75 0.75 0.75],'LineWidth',1);
    grid on
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    pbaspect([max(X)-min(X)+1,max(Y)-min(Y)+1,max(Z)-min(Z)+1]);
    view([90 0]);
    %
    subplot(2,2,4)
    hold on
    plot3(coorduniqueDEFx(1:NumInstrPoints),coorduniqueDEFy(1:NumInstrPoints),coorduniqueDEFz(1:NumInstrPoints),'o','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
    hold on
    plot3([coorduniquex,coorduniqueDEFx]', [coorduniquey,coorduniqueDEFy]' , [coorduniquez,coorduniqueDEFz]','-.','Color',[0.75 0.75 0.75],'LineWidth',1);
    grid on
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    pbaspect([max(X)-min(X)+1,max(Y)-min(Y)+1,max(Z)-min(Z)+1]);
    view([0 0]);

    if savefig == 1
        exportgraphics(fA,[cd,'\IdentModeShapePlot\',sprintf('IdMode#%d-fn%.0f[cHz]',mode,100*ID_RESULTS(1,mode)),'.png'],'Resolution',300);
    end

end

toc;
