%%Calculate point of parallel-perp. divergence
%data should be cut to only include resolvable data before running function

%H1 must be the perpendicular spectra and H2 should be the parallel spectra
function[div_x, div_y] = FindMinScale(L1, nH1, L2, nH2)

%determine longer dataset
%determine which is perpendicular and which is parallel
if nH1(end) <= nH2(end)
    H1 = nH1;
    H2 = nH2;
else
    H1 = nH2;
    H2 = nH1;
end
    

shortL = min(length(L1),length(L2));
for i = 5:shortL
    j = shortL+1-i;
    %j=i;
    if H2(j) <= 1.1236 *H1(j)
        div_x = L1(j);
        div_y = H1(j);
        return
    else
        if i == shortL
            ['Divergence Wavelength not found']
            div_x = NaN;
            div_y = NaN;
            return
        end
    end
end

% figure(9)
% loglog(L1,H1,L2,H2, '--', div_x, div_y, 'o')