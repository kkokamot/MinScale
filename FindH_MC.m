%% Monte Carlo RMS
%Inputs:
    % Z = gridded data to perform calculation oriented so that the
    % slip direction is either perpendicular or parallel to line
    % dx = spacing of gridded data
%Outputs: 
    % H1, H2: RMS Height in the perpendicular and parallel direction
    % L1, L2: Window length for each RMS Height (also in perpendicular and 
    % parallel)
    
function [H1,L1,H2,L2]=FindH_MC(Z,dx)
clearvars -except Z dx
 
sz = size(Z); %size of topographic map
[nrow,ncol] = size(Z);

% preallocate H1, L1, H2, L2
H1=NaN(1,nrow);
L1=H1;
H2=nan(1,ncol);
L2=H2;
for l = 1:2
    [' l:' num2str(l)] %prints 1 at the start of analyzing H1 and prints 2 at the start of H2
    if l ==2 
        Z = rot90(Z); %when running multiple times Z is rotated, switching H1 and H2, possibly switching how the graph is formatted
    end
    flipZ = fliplr(Z); %re-orientates Z

    % This for loop randomly selects sectperrow amount of starting points 
    % from each row to calculate RMS height from
    parfor iwin=5:sz(l) %iwin is window size and ranges from 5 to the size of a row. 
        %Nr is how many rows RMS will be computed for
        %limN is the length of the row analyzed
        sectperrow = 20; %this is how many windows will be analyzed
        S = nan(sectperrow,1); %pre-allocating S, the RMS values for each window analyzed for a specified window length
        if l == 2
            Nr = nrow; 
            S=nan(sectperrow,1);
            limN = ncol;
        else
            Nr = ncol;
            S = nan(sectperrow,1);
            limN = nrow;
        end
        
        %defining GpGinv before the for loop to save computation time while detrending
        G = [[0*(1:iwin)+1]; [1:iwin];]'; 
        Gp=G';
        GpGinv=inv(G'*G);
        
        
        diffIndex = iwin-1; %diffIndex is how many points needed to be added to make a full interval
        
        startrowindex = [1:limN-diffIndex];  %startrowindex chooses what starting points can be used
        
       
        startindexrand = randsample(startrowindex,sectperrow*Nr,true); 
        %startindexrand chooses random starting points (the number depends on sectperrow) for an interval from startrowrowindex
        %(startindexrand will later turn  windowrand of the for loop)
        
        colindex = 1:Nr;
        
        for k = 1:sectperrow
            for i = (k-1)*(Nr)+1:k*(Nr) %goes from 1:Nc k times
                windowrand = [startindexrand(i):startindexrand(i)+diffIndex];
                Zwin = Z(windowrand,colindex(i-((k-1)*Nr))); %change Z to nZ if using flip

                c=GpGinv*(Gp*Zwin);
                fhat = c(1)+c(2)*G(:,2);
                Z3 = Zwin-fhat;

                xhat = sum(Z3)./iwin;
                if l == 1
                   S(i) = sqrt(sum((Z3-xhat).^2)/(length(Z3)-1));
                else 
                   S(i)= sqrt(sum((Z3-xhat).^2)/(length(Z3)-1));
                end
            end
        end
       
        if l == 2
            H2(iwin)=nanmean(S(:));
            L2(iwin)=iwin*dx;
            imagesc(Zwin);
        else
            H1(iwin)=nanmean(S(:));
            L1(iwin)=iwin*dx;
        end
    end
end



