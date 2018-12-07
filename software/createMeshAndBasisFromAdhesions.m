function [myMesh]=createMeshAndBasisFromAdhesions(x_vec,y_vec,paxImage,displField,pixelSize,mask)
% createMeshAndBasisFromAdhesions creates mesh from adhesion information in
% paxImage. x_vec and y_vec serve as mask for mesh nodes.
% Sangyoon Han June 2013
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
%
% This file is part of TFM_Package.
% 
% TFM_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TFM_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TFM_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
disp('Creating force mesh using adhesion information in paxillin image ...')
paxImage = double(paxImage);
% Get Cell mask
pId = paxImage/max(paxImage(:));
% alpha = graythresh(pId);
%estimate the intensity level to use for thresholding the image
pIdCrop = pId(any(mask,2),any(mask,1));
level1 = graythresh(pIdCrop); %Otsu
[~, level2] = cutFirstHistMode(pIdCrop,0); %Rosin
alpha1 = 0.1*level2 + 0.9*level1;
alpha2 = 0.1*level1 + 0.9*level2;
figure,subplot(2,1,1),imshow(pIdCrop,[]),title('original image');
subplot(2,2,3),imshow(im2bw(pIdCrop,alpha1)),title(['Otsu-biased, alpha = ' num2str(alpha1) ]);
subplot(2,2,4),imshow(im2bw(pIdCrop,alpha2)),title(['Rosin-biased, alpha = ' num2str(alpha2)]);
alpha = input('type desired alpha to threshold the image so that it encompass entire cell membrane: ');

bwPI = im2bw(pId,alpha);
bwPI2 = bwmorph(bwPI,'clean');
bwPI3 = bwmorph(bwPI2,'erode',5);
bwPI3 = bwmorph(bwPI3,'dilate',4);
bwPI4 = bwmorph(bwPI3,'close',10);
bwPI4 = bwmorph(bwPI4,'dilate',1);
% bwPI5 = refineEdgeWithSteerableFilterGM(pId,bwPI4);
% In case that there is still islands, pick only the largest chunk
[labelPI,nChunk] = bwlabel(bwPI4);
if nChunk>1
    eachArea = zeros(nChunk,1);
    for k=1:nChunk
        currBWPI = labelPI==k;
        eachArea(k) = sum(currBWPI(:));
    end
    [~,indCellArea] = max(eachArea);
    bwPI4 = labelPI == indCellArea;
end
% Get FA mask
disp('detecting focal adhesions ...')
minSize = round((1000/pixelSize)*(200/pixelSize)); %adhesion limit=1um*.5um
maskFA = blobSegmentThreshold(paxImage,minSize,false,bwPI2);

% For each FA, get points in the area. 
[labelFA,nFA]=bwlabel(maskFA);
FANodeX=[];
FANodeY=[];
[x_mat, y_mat]=meshgrid(1:size(paxImage,2),1:size(paxImage,1));
skelMask = false(size(x_mat));
for k=1:nFA
    eachMaskFA = labelFA==k;
    % try skeleton
    eachMaskFAskel = bwmorph(eachMaskFA,'skel',Inf);
    idxEachFAskel = find(eachMaskFAskel);
    [eachFANodesY,eachFANodesX] = ind2sub(size(paxImage),idxEachFAskel(1:5:end));
    % from this, we can try to find if there is any points apart by 5 pixel
    % in the mask
    for p=1:length(eachFANodesX)
        skelMask = skelMask | sqrt((x_mat-eachFANodesY(p)).^2+(y_mat-eachFANodesX(p)).^2)<5;
    end
    marginMask = eachMaskFA & ~skelMask;
    margNodesY = [];
    margNodesX = [];
    if any(marginMask(:))
        [labelMargin,nMarg] = bwlabel(marginMask);
        for q=1:nMarg
            % get node at the margin skeleton
            eachMaskMarg = labelMargin == q;
            eachMargskel = bwmorph(eachMaskMarg,'skel',Inf);
            idxEachMargskel = find(eachMargskel);
            [eachMargNodesY,eachMargNodesX]=ind2sub(size(paxImage),idxEachMargskel(1:5:end));
            margNodesX = [margNodesX; eachMargNodesX];
            margNodesY = [margNodesY; eachMargNodesY];
        end
    end
    % nodes around FA periphery
    fatFA = bwmorph(eachMaskFA,'dilate',2);
    fatFAB=bwboundaries(fatFA,'noholes');
    fatFA_boundary = fatFAB{1};
    interval = ceil(500/pixelSize);
    FAbdNodesX = fatFA_boundary(1:interval:end,2);
    FAbdNodesY = fatFA_boundary(1:interval:end,1);

    % adding alll nodes regarding this FA
    FANodeX = [FANodeX; eachFANodesX; margNodesX; FAbdNodesX];
    FANodeY = [FANodeY; eachFANodesY; margNodesY; FAbdNodesY];
end
% Get equi-spaced nodes on the cell boundary
B=bwboundaries(bwPI4,'noholes');
boundary = B{1};
% Get rid of points at the image boundary
indAtImBD = boundary(:,1)==1 | boundary(:,1)==size(paxImage,1) ...
            | boundary(:,2)==1 | boundary(:,2)==size(paxImage,2);
boundary(indAtImBD,:)=[];
bdNodesX = boundary(1:9:end,2);
bdNodesY = boundary(1:9:end,1);
% mask for band from edge
iMask = imcomplement(bwPI4);
distFromEdge = bwdist(iMask);
bandwidth = 8; % in um
bandwidth_pix = round(bandwidth*1000/pixelSize);
bandMask = distFromEdge <= bandwidth_pix;
bandMask = bandMask & bwPI4;
band_FA = bandMask & ~maskFA;
% find nascent adhesions from paxImage
disp('detecting nascent adhesions ...')
sigmaPSF_NA = 1.48;
pstruct_NA = pointSourceDetection(paxImage, sigmaPSF_NA,...
    'alpha', 0.05, 'mask', band_FA);
NANodeX = ceil(pstruct_NA.x');
NANodeY = ceil(pstruct_NA.y');
% Add hexagonal grid in the band region
spacing = ceil(2000/pixelSize);
[hexNAX,hexNAY]=createHexGridInMask(round(spacing/4),bandMask);

NAs = [NANodeX NANodeY];
NAs = [NAs; bdNodesX bdNodesY;hexNAX hexNAY];

% nodes around NA
nNA = length(NANodeX);
for k=1:nNA
%     eachMaskNA = false(size(paxImage));
%     eachMaskNA(NANodeY(k),NANodeX(k))=1;
    fatNA = sqrt((x_mat-NANodeX(k)).^2+(y_mat-NANodeY(k)).^2)<6;
    fatNAB=bwboundaries(fatNA,'noholes');
    fatNA_boundary = fatNAB{1};
    interval = ceil(400/pixelSize);
    NAbdNodesX = fatNA_boundary(1:interval:end,2);
    NAbdNodesY = fatNA_boundary(1:interval:end,1);
    NAs = [NAs; NAbdNodesX NAbdNodesY];
end

NANodeX = NAs(:,1);
NANodeY = NAs(:,2);
% Get intNodes in interiorMask
interiorMask = bwPI4 & ~bandMask & ~maskFA;
[pstruct_intNA] = pointSourceDetection(paxImage, sigmaPSF_NA*1.5,...
    'alpha', 0.05, 'mask', interiorMask);
intNAs = [ceil(pstruct_intNA.x') ceil(pstruct_intNA.y')];
% Add hexagonal grid in the interior region
spacing = ceil(2000/pixelSize);
[intHexNAX,intHexNAY]=createHexGridInMask(spacing,interiorMask);
intNAs = [intNAs; intHexNAX,intHexNAY];
% intNAs are separated by 20 pixel 
idx = KDTreeBallQuery(intNAs, intNAs, 20);
valid = true(numel(idx),1);
for i = 1:numel(idx)
    if ~valid(i), continue; end
    neighbors_KD = idx{i}(idx{i}~=i);
    valid(neighbors_KD) = false;
end
intNAs = intNAs(valid, :);
intNANodeX = intNAs(:,1);
intNANodeY = intNAs(:,2);
% Get hexagonal nodes with 2000 nm spacing for bgdMask
bgdMask = ~bwPI4;
[bgdNodeX,bgdNodeY]=createHexGridInMask(spacing,bgdMask);

% For all nodes, get their undeformed postitions using displField
allNodesX=[FANodeX; NANodeX; intNANodeX; bgdNodeX];
allNodesY=[FANodeY; NANodeY; intNANodeY; bgdNodeY];

undNodesX=allNodesX;
undNodesY=allNodesY;

for ii=1:length(allNodesX)
    distToNode = sqrt(sum(([displField.pos(:,1) displField.pos(:,2)]- ...
                            ones(length(displField.pos(:,1)),1)*[allNodesX(ii) allNodesY(ii)]).^2,2));
                        
    [~,closest_ind] = min(distToNode);
    dispVec = [displField.vec(closest_ind,1) displField.vec(closest_ind,2)];
    undNodesX(ii) = allNodesX(ii) - dispVec(1);
    undNodesY(ii) = allNodesY(ii) - dispVec(2);
end

% filter undNodesX and Y with ROI
xmin = min(x_vec);
xmax = max(x_vec);
ymin = min(y_vec);
ymax = max(y_vec);

nPoints = length(allNodesX);
idxROI = false(nPoints,1);
for ii=1:nPoints
    if allNodesX(ii)>=xmin && allNodesX(ii)<=xmax ...
            && allNodesY(ii)>=ymin && allNodesY(ii)<=ymax
        idxROI(ii) = true;
    end
end
% Out-most nodes for background
inUndNodesX = undNodesX(idxROI);
inUndNodesY = undNodesY(idxROI);
% undeformed position of ROI boundary
ROINode = [xmin ymin;xmax ymin;xmax ymax;xmin ymax];
ROINodeX = ROINode(:,1);
ROINodeY = ROINode(:,2);
for ii=1:length(ROINodeX)
    distToNode = sqrt(sum(([displField.pos(:,1) displField.pos(:,2)]- ...
                            ones(length(displField.pos(:,1)),1)*[ROINodeX(ii) ROINodeY(ii)]).^2,2));
                        
    [~,closest_ind] = min(distToNode);
    dispVec = [displField.vec(closest_ind,1) displField.vec(closest_ind,2)];
    undROINodeX(ii) = ROINodeX(ii) - dispVec(1);
    undROINodeY(ii) = ROINodeY(ii) - dispVec(2);
end
undXmin = undROINodeX(1);
undXmax = undROINodeX(2);
undYmin = undROINodeY(1);
undYmax = undROINodeY(4);
% for left edge
bdUndNodesYleft = [undYmin:round(spacing/2):undYmax]';
bdUndNodesXleft = undXmin*ones(size(bdUndNodesYleft));
%for top edge
bdUndNodesXtop = [undXmin:round(spacing/2):undXmax]';
bdUndNodesYtop = undYmax*ones(size(bdUndNodesXtop));
%for right edge
bdUndNodesYright = [undYmax:-round(spacing/2):undYmin]';
bdUndNodesXright = undXmax*ones(size(bdUndNodesYright));
%for bottom edge
bdUndNodesXbottom = [undXmax:-round(spacing/2):undXmin]';
bdUndNodesYbottom = undYmin*ones(size(bdUndNodesXbottom));

bdUndNodesX = [bdUndNodesXleft; bdUndNodesXtop; bdUndNodesXright; bdUndNodesXbottom];
bdUndNodesY = [bdUndNodesYleft; bdUndNodesYtop; bdUndNodesYright; bdUndNodesYbottom];

undNodesX = [inUndNodesX; bdUndNodesX];
undNodesY = [inUndNodesY; bdUndNodesY];

% undNodes are separated at least by 5 pixel 
undNodes = [undNodesX undNodesY];
idx = KDTreeBallQuery(undNodes, undNodes, 5);
valid = true(numel(idx),1);
for i = 1:numel(idx)
    if ~valid(i), continue; end
    neighbors_KD = idx{i}(idx{i}~=i);
    valid(neighbors_KD) = false;
end
undNodes = undNodes(valid, :);
undNodesX = undNodes(:,1);
undNodesY = undNodes(:,2);

xvec=undNodesX;
yvec=undNodesY;
dt=DelaunayTri(xvec,yvec);
delaunay_mesh=dt.Triangulation;

%for each node find all its neighbors:
for n=1:length(xvec)
    candidates=[];
    for k=1:length(delaunay_mesh(:,1))
        if delaunay_mesh(k,1)==n || delaunay_mesh(k,2)==n || delaunay_mesh(k,3)==n
            for m=1:length(delaunay_mesh(1,:))
                if length(candidates)==0
                    if delaunay_mesh(k,m)~=n
                        candidates=horzcat(candidates,delaunay_mesh(k,m));
                    end                    
                elseif delaunay_mesh(k,m)~=n && isfinite(sum(1./(candidates-delaunay_mesh(k,m))))
                    candidates=horzcat(candidates,delaunay_mesh(k,m));
                end                    
            end
        end
    end
    neighbors(n).cand=sort(candidates);
    bounds(n).x=[min(xvec(candidates)) max(xvec(candidates))];
    bounds(n).y=[min(yvec(candidates)) max(yvec(candidates))];   
end

myMesh.p=[xvec,yvec];
myMesh.dt=dt;  % DelaunayTri structure
myMesh.neighbors=neighbors;
myMesh.bounds=bounds;
myMesh.numNodes=length(myMesh.p(:,1));
bdPtsID= convexHull(dt);

% meshBase = struct('f_intp_x',@(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(1).f_disc(:,1),'linear'),...
%                                     'f_intp_y',@(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(myMesh.numNodes+1).f_disc(:,1),'linear'));
% myMesh.base = repmat(meshBase, 2*myMesh.numNodes,1);                        
allPtsID = 1:myMesh.numNodes;
intPtsID = setdiff(allPtsID,bdPtsID);
myMesh.numBasis = length(intPtsID);
myMesh.baseID(myMesh.numNodes) = 0;

base = struct('f_disc',zeros(myMesh.numNodes,2));
base = repmat(base,2*myMesh.numBasis,1);
myMesh.base(2*myMesh.numBasis).f_intp_x = @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(myMesh.numBasis).f_disc(:,1),'linear'); 
myMesh.base(2*myMesh.numBasis).f_intp_y = @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(myMesh.numBasis*2).f_disc(:,1),'linear'); 
%create the basis functions and interpolate them using the Delaunay Triangulation:
p=0;
for j=intPtsID
    p=p+1;
    old_cputime = cputime;
    base(p).f_disc(j,1)=1;
    base(myMesh.numBasis+p).f_disc(j,2)=1;
    
    disp(['Creating ' num2str(p) 'th force base ... (' num2str(cputime-old_cputime) ' sec passed)'])
end
p=0;
for j=intPtsID
    p=p+1;
    old_cputime = cputime;
    myMesh.base(p).f_intp_x= @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(p).f_disc(:,1),'linear');
    myMesh.base(p).f_intp_y= @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(p).f_disc(:,2),'linear'); % only zeros
    myMesh.base(p).testNumber = @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(p).f_disc(:,2),'linear',j); % only zeros
    myMesh.base(myMesh.numBasis+p).f_intp_x= @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(myMesh.numBasis+p).f_disc(:,1),'linear'); % only zeros
    myMesh.base(myMesh.numBasis+p).f_intp_y= @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(myMesh.numBasis+p).f_disc(:,2),'linear'); 
    myMesh.base(p).nodeID = j;
    myMesh.baseID(j) = p; %used when calculating Gram matrix
    myMesh.base(p).node = myMesh.p(j,:);
    disp(['Creating ' num2str(p) 'th basis function ... (' num2str(cputime-old_cputime) ' sec passed)'])
end
    function vOut=nan2zeroTriScatteredInterp(x,y,dtIn,vIn,method,j)
        if nargin<6
            F=TriScatteredInterp(dtIn,vIn,method);
            vOut=F(x,y);
            checkVec=isnan(vOut);
            vOut(checkVec)=0;
        else
            vOut = j;
        end
    end
end
% % plot an example to see if it works correctly
% ind=80;
% if length(xvec)>ind-1
%     xmin=min(xvec);
%     ymin=min(yvec);
%     xmax=max(xvec);
%     ymax=max(yvec);
%     
%     pointsPerEdge=round(sqrt(length(xvec)));
%     [x_fine,y_fine]=meshgrid(linspace(xmin,xmax,10*pointsPerEdge) , linspace(xmin,xmax,10*pointsPerEdge));
% 
%     figure(10)
%     plot(myMesh.p(myMesh.neighbors(ind).cand,1),myMesh.p(myMesh.neighbors(ind).cand,2),'or')
%     hold on
%     plot(myMesh.p(ind,1),myMesh.p(ind,2),'ob')
%     triplot(myMesh.dt);
%     plot([myMesh.bounds(ind).x(1) myMesh.bounds(ind).x(1) myMesh.bounds(ind).x(2) myMesh.bounds(ind).x(2) myMesh.bounds(ind).x(1)],[myMesh.bounds(ind).y(1) myMesh.bounds(ind).y(2) myMesh.bounds(ind).y(2) myMesh.bounds(ind).y(1) myMesh.bounds(ind).y(1)],'k')
%     quiver(x_fine,y_fine,myMesh.base(ind).f_intp_x(x_fine,y_fine),myMesh.base(ind).f_intp_y(x_fine,y_fine),'r')
%     quiver(x_fine,y_fine,myMesh.base(myMesh.numNodes+ind).f_intp_x(x_fine,y_fine),myMesh.base(myMesh.numNodes+ind).f_intp_y(x_fine,y_fine),'g')
%     xlim([xmin xmax])
%     ylim([ymin ymax])
%     hold off
% end