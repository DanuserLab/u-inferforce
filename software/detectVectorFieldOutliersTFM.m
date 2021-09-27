function [outlierIndex, sparseIndex, r,neighborhood_distance] = detectVectorFieldOutliersTFM(data,varargin)
% detectVectorFieldOutliersTFM detect and return the outliers in a vector field
%
% Synopsis:        outlierIdx = detectVectorFieldOutliersTFM(data)
%                  [outlierIdx,r] = detectVectorFieldOutliersTFM(data,2)
%
% This function detects outliers within a vectorial field using an extended
% version of the 'median test' of Westerweel et al. 1994, adapted for PTV.
% After finding neighbors within an average bead distance using KDTree,
% the algorithm calculates the directional fluctuation and norm of vector fluctuation
% with respect to the neighborhood median residual for each vertex. 
% A threshold is then applied to this quantity to extract
% the outliers. Directional fluctuation weights more than vector norm
% fluctuation in usual TFM experiment.
%
% Input:
%      data - a vector field, i.e. a matrix of size nx4 where the two first
%      columns give the positions and the two last the displacement
%
%      threshold (optional) - a threshold for the detection criterion.
%      Usually values are between 2-4 depending on the stringency.
%
%      weighted (optional) - a boolean. If true, neighbors influence is
%      weighted using their relative distance to the central point.
%
% Output
%      outlierIndx - the index of the outlier along the first dimension of
%      data
%
%      r - the values of the normalized fluctuation for each element of the
%      vector field
%
% For more information, see:
% J. Westerweel & F. Scarano, Exp. Fluids, 39 1096-1100, 2005.
% J. Duncan et al., Meas. Sci. Technol., 21 057002, 2010.
%
% Copyright (C) 2021, Danuser Lab - UTSouthwestern 
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

% Sebastien Besson, Aug 2011
% Sangyoon Han, Mar 2015

% Input check
ip=inputParser;
ip.addRequired('data',@(x) size(x,2)==4);
ip.addOptional('threshold',2,@isscalar);
ip.addOptional('weighted' ,1,@isscalar);
ip.addParamValue('epsilon',.1,@isscalar);
ip.parse(data,varargin{:})
threshold=ip.Results.threshold;
weighted=ip.Results.weighted;
epsilon=ip.Results.epsilon;

% Filter out NaN from the initial data (but keep the index for the
% outliers)
ind=find(~isnan(data(:,3)));
dataNoNan=data(ind,:);
    
% Take out duplicate points (Sangyoon)
[dataU,idata,~] = unique(dataNoNan,'rows'); %data2 = data(idata,:),data = data2(iudata,:)
    
% calculate maximum closest distance
distance=zeros(length(dataU),1);
% distanceAll=cell(length(data),1); % the second closest distance
distance2=NaN(length(dataU),1); % the second closest distance
neiBeadsWhole = dataU(:,1:2);
% parfor i=1:length(data) % for each vector, we calculate closest distance and neighboring vectors in 5 times bigger distance vicinity, get the second closest distance
%     neiBeads = neiBeadsWhole;
%     % neiBeads(i,:)=[];
%     [~,distance(i)] = KDTreeClosestPoint(neiBeads,data(i,1:2));
%     [~,curDistanceAll] = KDTreeBallQuery(neiBeads,data(i,1:2),5*distance(i));
%     if length(curDistanceAll{1})>1
%         distance2(i) = curDistanceAll{1}(2);
%     else
%         distance2(i) = NaN;
%     end
% end
for scanDist=1:20 % Get the distance until each point has at least 2 neighbors
    [~,distance] = KDTreeBallQuery(neiBeadsWhole,dataU(:,1:2),scanDist);
    idxBeadsEnoughNeis = cellfun(@length,distance);
    if quantile(idxBeadsEnoughNeis,0.01)>=3
        break
    end
end
% Get the 2nd closest distance
distance2(idxBeadsEnoughNeis>1)=cellfun(@(x) x(2),distance(idxBeadsEnoughNeis>1));

% Discard vectors that are in sparse location
% This takes too much time, and distance2 is not gaussian-mixture
% opts = statset('maxIter', 200);
% objDist = cell(3,1);
% n = 0;
% lastwarn('');
% [~, msgidlast] = lastwarn;
% warning('off','stats:gmdistribution:FailedToConverge')
% warning('off','stats:gmdistribution:MissingData')
% fitError = false;
% model-like. Coverting to use median
% while ~strcmp(msgidlast,'stats:gmdistribution:FailedToConverge') && n<4 && ~fitError
%     n = n+1;
%     try
%         objDist{n} = gmdistribution.fit(distance2, n, 'Options', opts);
%         [~, msgidlast] = lastwarn;
%     catch
%         fitError=true;
%     end
% end
% if n>1
%     objDist = objDist(1:n-1);
%     [~,idx] = min(cellfun(@(i) i.BIC, objDist));
%     objDist = objDist{idx};
%     [mu,idx] = sort(objDist.mu);
%     svec = squeeze(objDist.Sigma(:,:,idx));
%     if length(mu)>1
%         threshDist= max(mu(1)+4*svec(1),mu(2)+4*svec(2));
%     else
%         threshDist= mu+4*svec;
%     end
% elseif n==1 && fitError
%     threshDist = 2*mean(distance2);
% end

% ----------- DeBugged by Waddah Moghram on 4/27/2019
threshDist = nanmean(distance2)+2*nanstd(distance2);
idxCloseVectors = (distance2 < threshDist) & ~isnan(distance2);                      
neighborhood_distance = 5*max(distance2(idxCloseVectors));                           %quantile(distance,0.95);%mean(distance);%size(refFrame,1)*size(refFrame,2)/length(beads); 
if isempty(neighborhood_distance)
   neighborhood_distance = 0;
end
idAwayVectors = (distance2>neighborhood_distance) | isnan(distance2);
idxCloseEnoughVectors = (distance2<neighborhood_distance) & ~isnan(distance2);          % Fixed by Waddah Moghram on 4/26/2019. added parantheses
dataFiltered = dataU(idxCloseEnoughVectors,:);
idCloseEnoughVectors = find(idxCloseEnoughVectors);
%------------------------------------

% Find neighbors and distances
% [idx,neiDist] = KDTreeBallQuery(data(:,1:2), data(:,1:2), neighborhood_distance);
[idx,neiDist] = KDTreeBallQuery(dataFiltered(:,1:2), dataFiltered(:,1:2), neighborhood_distance);

% For those who have too few neighbors, increase distance to include more
% neighbors - SH 170308
numNeis = cellfun(@numel,idx);
try
    warning('off','stats:gmdistribution:FailedToConverge')
    opts = statset('maxIter', 20);
    objDist = gmdistribution.fit(numNeis, 3, 'Options', opts);
    thresNumNeis = min(objDist.mu);
catch
    thresNumNeis = quantile(numNeis,0.1);
end
idxTooFewNeis = numNeis<thresNumNeis;
[idx(idxTooFewNeis),neiDist(idxTooFewNeis)] = KDTreeBallQuery(dataFiltered(:,1:2), dataFiltered(idxTooFewNeis,1:2), 2*neighborhood_distance);

% objDist.Sigma

% Get the triangulation edges and calculate all distances
% edges= tri.edges;
% dp=(data(edges(:,2),1:2)-data(edges(:,1),1:2));
% D= cellfun(@(x) sqrt(sum(x.^2)),neiDist);
% if ~weighted
%     D=ones(size(D));
% end

% Convert the edge list into an adjacency list 
% nodes = unique([edges(:,1)' edges(:,2)']);
% N=cell(numel(nodes),1);
% E=cell(numel(nodes),1);
% for e=1:size(edges,1); 
%     N{edges(e,1)}=[N{edges(e,1)},edges(e,2)]; 
%     E{edges(e,1)}=[E{edges(e,1)},e]; 
%     N{edges(e,2)}=[N{edges(e,2)},edges(e,1)]; 
%     E{edges(e,2)}=[E{edges(e,2)},e];
% end
N = idx;

% Measure weighted local and neighborhood velocities
options = {1:size(dataFiltered,1),'Unif',false};
% d =arrayfun(@(x)D(E{x}),options{:});
d = neiDist;
localVel=arrayfun(@(x)dataFiltered(x,3:4)/(median(d{x})+epsilon),options{:}); % Velocity (or displacement) vector 
% of individual points normalized by median distance to the points.
neighVel=arrayfun(@(x)dataFiltered(N{x},3:4)./repmat(d{x}+epsilon,1,2),options{:}); % Velocities (or displacements)
% of all the neighboring vectors

% Get median weighted neighborhood velocity
medianVel=cellfun(@median,neighVel,'Unif',false);

% Calculate normalized fluctuation using neighborhood residuals
medianRes=arrayfun(@(x) median(abs(neighVel{x}-repmat(medianVel{x},size(neighVel{x},1),1))),options{:});
normFluct = arrayfun(@(x) abs(localVel{x}-medianVel{x})./(medianRes{x}+epsilon),options{:});
r=cellfun(@norm, normFluct);

% Get median weighted magnited and orientation separately -added by SH
% 170307
critAng=pi/8;
rOri = cellfun(@(x,y) acos(x*y'/(norm(x)*norm(y))),medianVel,localVel);

% Filter outliers using threshold
% outlierIndex = ind(idata([idCloseVectors(r.^(1/2).*(rOri/critAng).^1.5>threshold)]));
outlierIndex = ind(idata(idCloseEnoughVectors(r>threshold)));
sparseIndex = ind(idata(idAwayVectors));
