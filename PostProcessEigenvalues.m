%% Define low fidelity sparse grid
S = smolyak_grid_quick_preset(6,3);
Sr = reduce_sparse_grid(S);

% Write to file
%writematrix(Sr.knots, "knotsTest.txt",'')

%% Load LoFi eigenvalue data
params = readmatrix("params.txt");
%bmkeigs = readmatrix("bmkeigs.txt");
scfeigs = readmatrix("scfeigs.txt");

% Define a high fidelity grid to interpolate onto
Shifi = smolyak_grid_quick_preset(6,6);
Shifi_r = reduce_sparse_grid(Shifi);

% Interpolate to hifi grid
hifi_data = interpolate_on_sparse_grid(S,Sr,scfeigs.',Shifi_r.knots);
hifi_params = interpolate_on_sparse_grid(S,Sr,params.',Shifi_r.knots);

eig_diffs = diff(hifi_data);

eig_diff_ratio = eig_diffs(2,:)./eig_diffs(1,:);

costfn = (eig_diff_ratio - 10).^2;

% Find discrete indices that minimise costfn
[minValues,minIndices] = mink(costfn,20);

% Save out parameters
minParams = hifi_params(:,minIndices);
writematrix(minParams,'factor10.txt','Delimiter','space');

%% Plot a realisation of the surrogate

minA = min(params(:,1));
maxA = max(params(:,1));
minR = min(params(:,4));
maxR = max(params(:,4));
minL = min(params(:,5));
maxL = max(params(:,5));

p2k = @(v,minV,maxV) 2*((v-minV)./(maxV-minV) -0.5);

params2knot= @(a,r,ap,rp,l,lp) [p2k(a,minA,maxA),p2k(r,-maxR,-minR),p2k(ap,minA,maxA),p2k(rp,minR,maxR),p2k(l,minL,maxL),p2k(lp,minL,maxL)];

r_of_interest = [-6:0.1:0].';
testPoints = params2knot(5*ones(size(r_of_interest)),...
    r_of_interest,...
    5*ones(size(r_of_interest)),...
    3*ones(size(r_of_interest)),...
    1*ones(size(r_of_interest)),...
    1*ones(size(r_of_interest)));
eigsTest = interpolate_on_sparse_grid(S,Sr,scfeigs.',testPoints.').';
figure
plot(r_of_interest,eigsTest)