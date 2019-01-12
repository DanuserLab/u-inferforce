function [reg_corner,ireg_corner,kappa,h]=regParamSelecetionLcurve(rho,eta,lambda,init_lambda, varargin)%,dataPath)
% [reg_corner,ireg_corner,kappa]=l_curve_corner(rho,eta,reg_param)
% returns l curve corner estimated using a maximum curvature (kappa) estimation 
% in log-log space
% rho is the misfit and eta is the model norm or seminorm
%
% INPUT
%   rho       - misfit
%   eta       - model norm or seminorm
%   reg_param - the regularization parameter
%   inflection - 0: L-corner, 1: Inflection point before corner (lambda smaller than l-corner), 2:
%   Inflectionpoint after corner  (lambda smaller than l-corner). If
%   the second derivative has a shape of which curvature approaches
%   assymtotically to zero, the algorithm uses strict definition of
%   curvature (instead of second derivative) to find L-corner, and find the
%   optimal lambda (by comparing solution norm difference and noise level
%   in non-cell area)...
%   manualSelection - true if you want to iterate selection process to get
%   better reg parameter. (default: false)
%
% OUTPUT
%   reg_corner  - the value of reg_param with maximum curvature
%   ireg_corner - the index of the value in reg_param with maximum curvature
%   kappa       - the curvature for each reg_param
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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

ip = inputParser;
ip.addRequired('rho',@isnumeric);
ip.addRequired('eta',@isnumeric);
ip.addRequired('lambda',@isnumeric);
ip.addOptional('init_lambda',median(lambda),@isnumeric);
ip.addParamValue('inflection',0,@isnumeric);
ip.addParamValue('manualSelection',false,@islogical);
ip.parse(rho,eta,lambda,init_lambda, varargin{:});
rho=ip.Results.rho;
eta=ip.Results.eta;
lambda=ip.Results.lambda;
init_lambda=ip.Results.init_lambda;
manualSelection=ip.Results.manualSelection;
inflection=ip.Results.inflection;

% transform rho and eta into log-log space
x=log(rho);
y=log(eta);

% calculating slopes by linear regression
x_slope = x(3:end-2);
slope = zeros(size(x_slope));
for k=1:length(x)-4 
    [~,slope(k,1),~] = regression(x(k:k+4)',y(k:k+4)');
end

% calculating kappas by linear regression
x_kappa = x_slope(2:end-1);
kappa = zeros(size(x_kappa));
lambda_cut = lambda(4:end-3);
for k=1:length(x_slope)-2 
    [~,kappa(k),~] = regression(x_slope(k:k+2)',slope(k:k+2)');
end

% kappa2 = diff(diff(y_cut)./diff(x_cut))./diff(x_cut(1:end-1));
% kappadiff = diff(kappa);

% find a local maximum with three sections
nSections = 3;
nPoints = length(kappa);
p=0;
maxKappaCandIdx = [];
for ii=1:nSections
    [~, maxKappaIdx] = max(kappa(floor((ii-1)*nPoints/nSections)+1:floor(ii*nPoints/nSections))); % this is right at the L-corner which is usually over-smoothing
    maxKappaIdx = maxKappaIdx + floor((ii-1)*nPoints/nSections);
    % check if this is truly local maximum
    if maxKappaIdx>1 && maxKappaIdx<nPoints && (kappa(maxKappaIdx)>kappa(maxKappaIdx-1) && kappa(maxKappaIdx)>kappa(maxKappaIdx+1))
        p=p+1;
        maxKappaCandIdx(p) = maxKappaIdx;
    end
end
if length(maxKappaCandIdx)==1
    maxKappaIdx = maxKappaCandIdx(1);
elseif length(maxKappaCandIdx)>1
    % pick the one which is closer to initial lambda
    [~,Idx_close] = min(abs(log(lambda_cut(maxKappaCandIdx))-log(init_lambda)));
    maxKappaIdx = maxKappaCandIdx(Idx_close);
%     [~,tempIndex] = max(kappa(maxKappaCandIdx));% use the first one %max(maxKappaCandIdx);
%     maxKappaIdx = maxKappaCandIdx(tempIndex);
elseif isempty(maxKappaCandIdx)
    disp('There is no local maximum in curvature in the input lambda range.Using global maximum instead ...');
    [~, maxKappaIdx] = max(kappa);
end
if inflection==1 % if inflection point larger than lcorner is to be chosen.
    inflectionIdx = find(kappa<0 & (1:nPoints)'>maxKappaIdx,1,'first');
    ireg_corner= inflectionIdx+3;%round((maxKappaIdx+maxKappaDiffIdx)/2); % thus we choose the mean of those two points.
    reg_corner = lambda_cut(inflectionIdx);
    disp(['L-inflection value (larger than L-corner): ' num2str(reg_corner)])
elseif inflection==2 % if inflection point smaller than lcorner is to be chosen.
    inflectionIdx = find(kappa<0 & (1:nPoints)'<maxKappaIdx,1,'last');
    ireg_corner= inflectionIdx+3;%round((maxKappaIdx+maxKappaDiffIdx)/2); % thus we choose the mean of those two points.
    reg_corner = lambda_cut(inflectionIdx);
    disp(['L-inflection value (smaller than L-corner): ' num2str(reg_corner)])
else
    if maxKappaIdx==1 || sum(kappa>0)/length(kappa)>0.8 % if kappa is assymtotically approaching to zero from large positive...
        [reg_corner,ireg_corner,kappa]=l_curve_corner(rho,eta,lambda);
        kappa = kappa(4:end-3);
    else
        ireg_corner= maxKappaIdx+3;%round((maxKappaIdx+maxKappaDiffIdx)/2); % thus we choose the mean of those two points.
        reg_corner = lambda_cut(maxKappaIdx);
    end
    disp(['L-corner regularization parameter value: ' num2str(reg_corner)])
end
% [~, maxKappaDiffIdx] = max(kappadiff(1:maxKappaIdx)); %  this is steepest point right before L-corner. This is usually too small.
% find an index at kappa = 0 before maxKappaIdx

% ireg_corner = corner(rho,eta,false);

% disp(['Initial L-corner regularization parameter value: ' num2str(reg_corner) '.'])

% % poly-fit version
% firstKIdx = max(1,maxKappaIdx-2);
% lastKIdx = min(length(kappa),maxKappaIdx+2);
% p=polyfit(firstKIdx:lastKIdx,kappa(firstKIdx:lastKIdx)',2);
% if p(1)<0
%     ireg_corner = -p(2)/(2*p(1));
%     q=polyfit(firstKIdx:lastKIdx,lambda_cut(firstKIdx:lastKIdx),1);
%     reg_corner = polyval(q,ireg_corner);
%     disp(['Sub-knot resolution L-corner regularization parameter value: ' num2str(reg_corner) '.'])
% else
%     disp('The corner''''s L-corner does not have positive curvature')
%     reg_corner = lambda_cut(ireg_corner);
% end

if manualSelection
%     numCutPoints = 0;

    % show the l curve and make sure this is fittable with 5th order polynomial
%     x_cut = x((numCutPoints+1:end));
    
    h=figure; set(h,'Position',[1000,100,350,600])
    subplot(3,1,1),plot(x,y,'k')
    xlabel('Residual Norm ||Gm-d||_{2}');
    ylabel('Simi-Norm ||Lm||_{2}');
    subplot(3,1,1), hold on,plot(x(ireg_corner),y(ireg_corner),'ro')
    text(x(ireg_corner),1.01*y(ireg_corner),...
        ['    ',num2str(reg_corner,'%5.3e')]);
    subplot(3,1,2), plot(x_slope,slope), title('slope')
    if ireg_corner-2>0
        subplot(3,1,2), hold on,plot(x_slope(ireg_corner-2),slope(ireg_corner-2),'ro')
    end
    subplot(3,1,3), plot(x_kappa,kappa), title('curvature')
%     subplot(3,1,3), plot(x_cut(3:end-5),diff(kappa)),title('jerk')
    if ireg_corner-3>0
        subplot(3,1,3), hold on,plot(x_kappa(ireg_corner-3),kappa(ireg_corner-3),'ro')
    end
%     subplot(3,1,3), hold on,plot(x_cut(maxKappaIdx),kappadiff(maxKappaIdx),'ro')

%     poly5ivity = input('Is the curve going down with two concaveness (y/n)?','s');
%     while poly5ivity == 'n'
%         numCutPoints = input('how many entry points do you want to eliminate from the beginning?');
%         subplot(3,1,1),plot(x(numCutPoints+1:end),y(numCutPoints+1:end),'k')
%         x_cut = x((numCutPoints+1:end));
%         y_cut = y((numCutPoints+1:end));
% 
%         kappa = diff(diff(y_cut)./diff(x_cut))./diff(x_cut(1:end-1));
%         subplot(3,1,2), plot(x_cut(1:end-2),kappa)
%         kappadiff = diff(kappa);
%         subplot(3,1,3), plot(x_cut(1:end-3),diff(kappa))
% 
%         p=0;
%         maxKappaCandIdx = [];
%         nPoints = length(kappa);
%         for ii=1:nSections
%             [~, maxKappaIdx] = max(kappa(floor((ii-1)*nPoints/nSections)+1:floor(ii*nPoints/nSections))); % this is right at the L-corner which is usually over-smoothing
%             maxKappaIdx = maxKappaIdx + floor((ii-1)*nPoints/nSections);
%             % check if this is truly local maximum
%             if maxKappaIdx>1 && maxKappaIdx<nPoints && (kappa(maxKappaIdx)>kappa(maxKappaIdx-1) && kappa(maxKappaIdx)>kappa(maxKappaIdx+1))
%                 p=p+1;
%                 maxKappaCandIdx(p) = maxKappaIdx;
%             end
%         end
%         if length(maxKappaCandIdx)==1
%             maxKappaIdx = maxKappaCandIdx(1);
%         elseif length(maxKappaCandIdx)>1
%             maxKappaIdx = max(maxKappaCandIdx);
%         elseif isempty(maxKappaCandIdx)
%             error('there is no local maximum in curvature in the input lambda range');
%         end
%     %     [~, maxKappaDiffIdx] = max(kappadiff(1:maxKappaIdx)); %  this is steepest point right before L-corner. This is usually too small.
%         % find an index at kappa = 0 before maxKappaIdx
%         ireg_corner= numCutPoints+maxKappaIdx;%round((maxKappaIdx+maxKappaDiffIdx)/2); % thus we choose the mean of those two points.
%         subplot(3,1,1), hold on,plot(x(ireg_corner),y(ireg_corner),'ro')
%         text(x(ireg_corner),1.1*y(ireg_corner), ['    ',num2str(lambda(ireg_corner),'%5.3e')]);
% 
%         subplot(3,1,2), hold on,plot(x_cut(ireg_corner-numCutPoints),kappa(ireg_corner-numCutPoints),'ro')
%         subplot(3,1,3), hold on,plot(x_cut(ireg_corner-numCutPoints),kappadiff(ireg_corner-numCutPoints),'ro')
% 
%         poly5ivity = input('Is the curve going down with two concaveness (y/n)?','s');
%     end
%     reg_corner = lambda(ireg_corner);
end
% % fit it in 5th order polynomial - fitting with polynomial is dangerous!
% f = fit(x(numCutPoints+1:end), y(numCutPoints+1:end),  'poly5');
% hold on, plot(f)
% [~, fxx] = differentiate(f,x(numCutPoints+1:end));
% 
% % third derivative for finding local peak in curvature
% p = [f.p1 f.p2 f.p3 f.p4 f.p5 f.p6];
% 
% pxxx = polyder(polyder(polyder(p)));
% peaks2 = roots(pxxx);
% peak = peaks2(1);
% [~,peakIdx]=min(abs(x(numCutPoints+1:end)-peak)); % this is right at the L-corner which is over smoothing
% 
% pxxxx = polyder(polyder(polyder(polyder(p))));
% curvCurv = roots(pxxxx); % curvature of curvature
% [~,curvCurvIdx]=min(abs(x(numCutPoints+1:end)-curvCurv)); % this is steepest point right before L-corner. This is usually too small.
% 
% ireg_corner= round((peakIdx+curvCurvIdx)/2); % thus we choose the mean of those two points.
% reg_corner = lambda(numCutPoints+ireg_corner);
% kappa = fxx;
% 
% hold on
% plot(x(ireg_corner+numCutPoints),y(ireg_corner+numCutPoints),'ro')
% % print(h,[dataPath filesep 'Lcurve.eps'],'-depsc')
% close(h)
% find a positive peak in curvature (diff based, discrete)



