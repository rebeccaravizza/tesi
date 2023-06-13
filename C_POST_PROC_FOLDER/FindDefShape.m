function [coorduniqueDEFx , coorduniqueDEFy , coorduniqueDEFz , coorduniquex , coorduniquey , coorduniquez , relpointgld] = FindDefShape(gdl,X,Y,Z,di,sha,lambda)

coord = [X Y Z];
dist = nan(length(X));
ind0 = [];
conta = 0;
point = nan(length(X),1);
for ii = 1:length(X)
    for jj = ii+1:length(X)
        dist(ii,jj) = sqrt( nansum( ( coord(ii,:)-coord(jj,:) ).^2 ) );
    end
    ind0 = [ind0 ; find(dist(ii,:)==0)'];

    if isempty(intersect(ii,ind0))
        conta = conta + 1;
        point([ii;find(dist(ii,:)==0)'],1) = conta;
    end
end
[pointuniq,indunique] = unique(point);
coordunique = coord(indunique,:);
coorduniquex = coordunique(:,1);
coorduniquey = coordunique(:,2);
coorduniquez = coordunique(:,3);
relpointgld = [gdl , point];

% sopra lo script è valido per tutte le forme, da sotto va ripetuto con
% ciclo for ser ogni forma

for ii = 1:length(pointuniq)
    indsame = find(point==pointuniq(ii));
    dipointii = di(indsame);
    indsameX = find(dipointii==1);
    indsameY = find(dipointii==2);
    indsameZ = find(dipointii==3);
    indpointiiX = indsame(indsameX);
    indpointiiY = indsame(indsameY);
    indpointiiZ = indsame(indsameZ);

    if isempty(sha(indpointiiX))
        shaiiX = 0;
    else
        shaiiX = sha(indpointiiX);
    end
    if isempty(sha(indpointiiY))
        shaiiY = 0;
    else
        shaiiY = sha(indpointiiY);
    end
    if isempty(sha(indpointiiZ))
        shaiiZ = 0;
    else
        shaiiZ = sha(indpointiiZ);
    end

    shaii = [shaiiX , shaiiY , shaiiZ];
    coorduniqueDEF(ii,:) = coordunique(ii,:) + lambda*shaii;
end
coorduniqueDEFx = coorduniqueDEF(:,1);
coorduniqueDEFy = coorduniqueDEF(:,2);
coorduniqueDEFz = coorduniqueDEF(:,3);

end