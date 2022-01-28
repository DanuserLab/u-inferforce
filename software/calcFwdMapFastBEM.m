function [M]=calcFwdMapFastBEM(x_vec_u, y_vec_u, forceMesh, E,varargin)
% Synoposis:  M=calcFwdMapFastBEM(x_vec_u, y_vec_u, forceMesh, E)
%
% Copyright (C) 2022, Danuser Lab - UTSouthwestern 
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
ip.addRequired('x_vec_u',@isnumeric);
ip.addRequired('y_vec_u',@isnumeric);
ip.addRequired('forceMesh',@isstruct);
ip.addRequired('E',@isscalar);
ip.addOptional('meshPtsFwdSol',[],@isscalar);
ip.addOptional('doPlot',0,@isscalar);
ip.addParameter('basisClassTblPath','basisClassTbl.mat',@ischar);
ip.addParameter('wtBar',-1,@isscalar);
ip.addParameter('imgRows',1024,@isscalar);
ip.addParameter('imgCols',1334,@isscalar);
ip.addParameter('thickness',472,@isscalar); % default assuming 34 um with 72 nm/pix resolution
ip.addParameter('PoissonRatio',0.5,@isscalar); 
ip.parse(x_vec_u, y_vec_u, forceMesh, E,varargin{:})
meshPtsFwdSol=ip.Results.meshPtsFwdSol;
doPlot=ip.Results.doPlot;
basisClassTblPath=ip.Results.basisClassTblPath;
wtBar=ip.Results.wtBar;
imgRows = ip.Results.imgRows;
imgCols = ip.Results.imgCols;
thickness = ip.Results.thickness;    
v = ip.Results.PoissonRatio;

% try to load the lookup table:
try
    basisClassTblData=load(basisClassTblPath);
    basisClassTbl=basisClassTblData.basisClassTbl;
    addAtLeastOneToTbl=0;
catch ME
    basisClassTbl=struct([]) ;
    addAtLeastOneToTbl=1;
end


% for a field of view with 1024x1344 pix and a square grid with grid size
% of 10 pix, the average rel. difference between the solutions with 2^11 and
% 2^12 pts for the fwd solution is less than 0.27% for all positions with
% non-vanishing forces. The maximum error is 4%. These results indicate
% that actually 2^11 points are fine for calculating the basis solutions
% for such a configuration. 2^12 should be fine for sure!
forceSpan=1;
% imgRows=1024;
% imgCols=1344;
    

% transform to column vectors:
x_vec_u=x_vec_u(:);
y_vec_u=y_vec_u(:);


% test if displacement vectors have been measured only at integer
% positions:
diff_x_u = abs(x_vec_u-round(x_vec_u));
diff_y_u = abs(y_vec_u-round(y_vec_u));

% test if the basis function for the force are located only at integer
% positions:
allNodes = vertcat(forceMesh.basis(:).node);
diff_xy_f = abs(allNodes-round(allNodes));

if sum(diff_x_u(:))+sum(diff_y_u(:))<10^(-3) && sum(diff_xy_f(:))<10^(-3)
    method='direct';
else
    method='*cubic';
end

% determine the size of M
numBasis=length(forceMesh.basis);
numPts  =length(x_vec_u);

% Initialize M
M=NaN(2*numPts,2*numBasis);

% here actually each basis function hast to be evaluated twice, which means
% it actually give two values:
ux=zeros(length(x_vec_u),2*numBasis);
uy=zeros(length(y_vec_u),2*numBasis);

% To make sure that the range over which the solution is calculated,
% take the double of the initial x and y ranges:
xmin=min(x_vec_u); xmax=max(x_vec_u);
ymin=min(y_vec_u); ymax=max(y_vec_u);
dx=xmax-xmin;
dy=ymax-ymin;

% The minimum x/y-range over which the basis solution has/had to be
% calculated:
xrangeReq=[-dx dx]';
yrangeReq=[-dy dy]';
    
% calculate the basis solution for all basis classes:
numClass=length(forceMesh.basisClass);
for class=1:numClass
    display(['Work on class: ',num2str(class),' of: ',num2str(numClass),'. Each class [~2x2min]:... ']);
    
    % Integration bounds used in the refine step, in case the fwdSolution
    % has to be calculated:
    xbd_min=min(forceMesh.basisClass(class).neighPos(:,1));
    xbd_max=max(forceMesh.basisClass(class).neighPos(:,1));
    ybd_min=min(forceMesh.basisClass(class).neighPos(:,2));
    ybd_max=max(forceMesh.basisClass(class).neighPos(:,2));

    % try to find the basis class in the table base:
    basisClassIn=forceMesh.basisClass(class);
    [idBestMatch]=findBasisClassInTbl(basisClassTbl,basisClassIn,xrangeReq,yrangeReq,meshPtsFwdSol);
    
    % now run through the x/y-comp. and either pull the fwdSol from the
    % table base or calculate it from scratch. The latter will be saved in
    % the basisClassTbl for later use
    for oneORtwo=1:2        
        if ~isempty(idBestMatch)
            % Then we can use the stored solution.
            % scale the basis solution with the right Youngs modulus. This
            % works for the boussinesq-BCs but might fail for more general BCs:
            scaleE = basisClassTbl(end).uSol.E/E; % for details see Landau Lifschitz p32.
            ux_model_pix = scaleE*double(basisClassTbl(idBestMatch).uSol.comp(oneORtwo).ux);
            uy_model_pix = scaleE*double(basisClassTbl(idBestMatch).uSol.comp(oneORtwo).uy);
            x_model_pix  = double(basisClassTbl(idBestMatch).uSol.x);
            y_model_pix  = double(basisClassTbl(idBestMatch).uSol.y);
        else
            display('Could not find a good match in the tablebase! Have to calculate the solution!')
            % no good match has been found in the table, we have to calculate the
            % solution from scratch. Here we could force xrangeSol>xrangeReq
            % and meshPtsFwdSol>meshPtsFwdSol_min, such that it is more
            % likely that this solution can be used in the future. We will
            % store this solution only for method='direct'!
            if forceSpan==1
                dxSol=max(imgCols,dx);
                dySol=max(imgRows,dy);
                
                xrangeSol=[-dxSol dxSol]';
                yrangeSol=[-dySol dySol]';
            else
                xrangeSol=xrangeReq;
                yrangeSol=yrangeReq;
            end               
            % calculate the solution:
            % Here we use finite thickness for calculating Green's function
            [ux_model, uy_model, x_model, y_model]=fwdSolution(xrangeSol,yrangeSol,E,...
                xbd_min,xbd_max,ybd_min,ybd_max,forceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_x,forceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_y,...
                'fft','noIntp',meshPtsFwdSol,thickness,v); %'fft_finite','noIntp',meshPtsFwdSol,thickness);
            % check if the sampling is fine enough for method 'direct':
            if strcmp(method,'direct') || strcmp(method,'*cubic')
                % This works perfectly for all mesh types as long as the
                % displacment and force nodes are defined at integer positions!
                x_spacing=x_model(2,2)-x_model(1,1);
                y_spacing=y_model(2,2)-y_model(1,1);
                if x_spacing<=1 && y_spacing<=1
                    % Only if the spacing is <1 we have oversampled, interpolate to
                    % integer positions:
                    x_model_pix=x_model(1,1):1:x_model(end,end);
                    y_model_pix=y_model(1,1):1:y_model(end,end);
                    
                    [x_model_pix,y_model_pix]=meshgrid(x_model_pix,y_model_pix);
                    
                    %interpolate the solution to the integer positions:
                    ux_model_pix= interp2(x_model, y_model, ux_model, x_model_pix, y_model_pix); %, 'direct'); There is no such thing as direct, but only 'linear'  %This is ux(:,j)
                    uy_model_pix= interp2(x_model, y_model, uy_model, x_model_pix, y_model_pix); %, 'direct');  %This is uy(:,j)
                else
                    display('Have switched over to *cubic. But is this really necessary? It might well be that even if undersampled the upper search will produce the same result as an interpolation')
                    method='*cubic';
                    pizInterval_x = round(x_spacing);
                    pizInterval_y = round(y_spacing);
                    x_model_pix=x_model(1,1):pizInterval_x:x_model(end,end);
                    y_model_pix=y_model(1,1):pizInterval_y:y_model(end,end);
                    
                    [x_model_pix,y_model_pix]=meshgrid(x_model_pix,y_model_pix);
                    
                    %interpolate the solution to the integer positions:
                    ux_model_pix= interp2(x_model, y_model, ux_model, x_model_pix, y_model_pix, 'cubic');  %This is ux(:,j)
                    uy_model_pix= interp2(x_model, y_model, uy_model, x_model_pix, y_model_pix, 'cubic');  %This is uy(:,j)
                end
            end
        end
        
        toDoBasis=find(vertcat(forceMesh.basis.class)==class)';
        lgthToDoBasis=length(toDoBasis);
        display(['Evaluate ',num2str(lgthToDoBasis),' basis functions'])
        
        
        logMsg = 'Please wait, interpolating basis solutions';
        timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
        tic;
        if ishandle(wtBar)
            wtBar = waitbar(0,wtBar,logMsg);
%         elseif feature('ShowFigureWindows'),
%             wtBar = waitbar(0,logMsg);
        end
        
        for i=1:numel(toDoBasis)
            basisID=toDoBasis(i);
            % lgthToDoBasis=length(toDoBasis);
            % displayText=[num2str(basisID),' of ',num2str(lgthToDoBasis)];
            % progressText(basisID/lgthToDoBasis,displayText);
            % Interpolate the basis-solution:
            xShift = forceMesh.basis(basisID).node(1);
            yShift = forceMesh.basis(basisID).node(2);

            
            if strcmp(method,'direct')
                % Instead of interpolation we can simply search for the
                % correct values in the matrix:
                xmin_pix=x_model_pix(1,1)+xShift;
                ymin_pix=y_model_pix(1,1)+yShift;
                
                % shift the x,y values to indices:
                px=x_vec_u-xmin_pix+1; % +1 because the matrix index starts at 1.
                py=y_vec_u-ymin_pix+1; % +1 because the matrix index starts at 1.
                
                % get the indices:
                [ptInd]=sub2ind(size(x_model_pix),py,px);
                
                if oneORtwo==1
                    M(1:numPts    ,basisID)          = ux_model_pix(ptInd);
                    M(numPts+1:end,basisID)          = uy_model_pix(ptInd);
                elseif oneORtwo==2
                    M(1:numPts    ,basisID+numBasis) = ux_model_pix(ptInd);
                    M(numPts+1:end,basisID+numBasis) = uy_model_pix(ptInd);
                end
                
%                 % This might be a bit more robust but is slower:
%                 forceNearest=0;
%                 if forceNearest==1
%                     if oneORtwo==1
%                         % Then the interpolants of the first function are:
%                         M(1:numPts    ,basisID)          = interp2(x_model_pix+xShift, y_model_pix+yShift, ux_model_pix, x_vec_u, y_vec_u, '*nearest');  %This is ux(:,j)
%                         M(numPts+1:end,basisID)          = interp2(x_model_pix+xShift, y_model_pix+yShift, uy_model_pix, x_vec_u, y_vec_u, '*nearest');  %This is uy(:,j)
%                     elseif oneORtwo==2
%                         % Then the interpolants of the second function are:  (:,j+numBasis)
%                         M(1:numPts    ,basisID+numBasis) = interp2(x_model_pix+xShift, y_model_pix+yShift, ux_model_pix, x_vec_u, y_vec_u, '*nearest');  %This is ux(:,j+numBasis)
%                         M(numPts+1:end,basisID+numBasis) = interp2(x_model_pix+xShift, y_model_pix+yShift, uy_model_pix, x_vec_u, y_vec_u, '*nearest');  %This is uy(:,j+numBasis)
%                     end
%                 end
            elseif strcmp(method,'*cubic') && ~isempty(idBestMatch)
                if oneORtwo==1
                    % Then the interpolants of the first function are:
                    M(1:numPts    ,basisID)          = interp2(x_model_pix+xShift, y_model_pix+yShift, ux_model_pix, x_vec_u, y_vec_u, method);  %This is ux(:,j)
                    M(numPts+1:end,basisID)          = interp2(x_model_pix+xShift, y_model_pix+yShift, uy_model_pix, x_vec_u, y_vec_u, method);  %This is uy(:,j)
                elseif oneORtwo==2
                    % Then the interpolants of the second function are:  (:,j+numBasis)
                    M(1:numPts    ,basisID+numBasis) = interp2(x_model_pix+xShift, y_model_pix+yShift, ux_model_pix, x_vec_u, y_vec_u, method);  %This is ux(:,j+numBasis)
                    M(numPts+1:end,basisID+numBasis) = interp2(x_model_pix+xShift, y_model_pix+yShift, uy_model_pix, x_vec_u, y_vec_u, method);  %This is uy(:,j+numBasis)
                end
            else % we almost do not need this condition! SH
                if oneORtwo==1
                    % Then the interpolants of the first function are:
                    M(1:numPts    ,basisID)          = interp2(x_model+xShift, y_model+yShift, ux_model, x_vec_u, y_vec_u, method);  %This is ux(:,j)
                    M(numPts+1:end,basisID)          = interp2(x_model+xShift, y_model+yShift, uy_model, x_vec_u, y_vec_u, method);  %This is uy(:,j)
                elseif oneORtwo==2
                    % Then the interpolants of the second function are:  (:,j+numBasis)
                    M(1:numPts    ,basisID+numBasis) = interp2(x_model+xShift, y_model+yShift, ux_model, x_vec_u, y_vec_u, method);  %This is ux(:,j+numBasis)
                    M(numPts+1:end,basisID+numBasis) = interp2(x_model+xShift, y_model+yShift, uy_model, x_vec_u, y_vec_u, method);  %This is uy(:,j+numBasis)
                end
            end
            
            % Update the waitbar
            if mod(i,5)==1 && ishandle(wtBar)
                ti=toc;
                waitbar(i/lgthToDoBasis,wtBar,...
                    sprintf([logMsg timeMsg(ti*lgthToDoBasis/i-ti)]));
            end
            
        end
%         % Close waitbar if generated by the function
%         if ishandle(wtBar) && ~ishandle(ip.Results.wtBar), 
%             close(wtBar); 
%         end
%         
        if isempty(idBestMatch) %&& strcmp(method,'direct')
            % Then we have either calculated a previously unknown basis
            % Solution, or we have improved one (by increasing the range or
            % by increasing the meshPtsFwdSol). Enter a new entry only if
            % we are working on the x-comp. The y-comp will be treated in
            % the next loop and will be sorted in into the same
            % basisClassTbl-id:
            if oneORtwo==1
                numClassTbl=length(basisClassTbl);
                currBasisClass=forceMesh.basisClass(class);
                % strip of the basisFunc entry. We don't need this
                % information!
                basisClassTbl(numClassTbl+1).centerPos  = currBasisClass.centerPos;
                basisClassTbl(numClassTbl+1).numNeigh   = currBasisClass.numNeigh;
                basisClassTbl(numClassTbl+1).neighPos   = currBasisClass.neighPos;
                basisClassTbl(numClassTbl+1).dtBaseSup  = currBasisClass.dtBaseSup;
                basisClassTbl(numClassTbl+1).unitVolume = currBasisClass.unitVolume;
            end
            % enter the basis solutions:
            % Scale the basis solution with the right Youngs modulus. This
            % works for the boussinesq-BCs but might fail for more general BCs:
            % To store the forward solution, single precision should be
            % sufficient. x/y positions are integer anyways, store them in
            % int16 format:
            basisClassTbl(end).uSol.comp(oneORtwo).ux = single(ux_model_pix*E); % the factor E is to scale the solution to u(E=1,f)
            basisClassTbl(end).uSol.comp(oneORtwo).uy = single(uy_model_pix*E); % the factor E is to scale the solution to u(E=1,f)
            basisClassTbl(end).uSol.x = int16(x_model_pix);
            basisClassTbl(end).uSol.y = int16(y_model_pix);
            
            % enter parameters:
            basisClassTbl(end).uSol.xrange       = xrangeSol;
            basisClassTbl(end).uSol.yrange       = yrangeSol;
            basisClassTbl(end).uSol.E                = 1; % this could be more general!
            basisClassTbl(end).uSol.method       ='fft';
            basisClassTbl(end).uSol.meshPtsFwdSol= meshPtsFwdSol;
            basisClassTbl(end).uSol.gelHeight   = thickness; % this could be more general!
            addAtLeastOneToTbl=1;
        end
    end    
end
% save the new table base
if addAtLeastOneToTbl
    save(basisClassTblPath, 'basisClassTbl','-v7.3');
end

% plot an example to see if it works correctly
if doPlot==1
    ind=1;
    if forceMesh.numBasis>ind-1
        xmin=min(x_vec_u);
        ymin=min(y_vec_u);
        xmax=max(x_vec_u);
        ymax=max(y_vec_u);
        ux = basisClassTbl(end).uSol.comp(oneORtwo).ux;
        uy = basisClassTbl(end).uSol.comp(oneORtwo).uy;
        figure(11)
        quiver(x_vec_u,y_vec_u,ux(:,ind),uy(:,ind))
        hold on
        quiver(x_vec_u,y_vec_u,ux(:,ind+forceMesh.numBasis),uy(:,ind+forceMesh.numBasis))
        xlim([xmin xmax])
        ylim([ymin ymax])
        hold off
    end
end