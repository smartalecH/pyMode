
if 1
    % FMM
    prkna = [ 10   55.447   0.32065
              15   85.141   0.21469
              20  115.122   0.14044
              25  145.240   0.09312
              30  175.448   0.06134
              35  205.724   0.03998
              40  236.052   0.02578
              45  266.424   0.01645
              50  296.832   0.01040 ];
end

if 0
    % [20]
    prkna = [ 15   85.135   0.21408
              20  115.109   0.14560
              25  145.231   0.10923
              30  175.448   0.08201
              35  205.732   0.05960
              40  236.067   0.04180
              45  266.442   0.02840
              50  296.851   0.01877 ];
end

if 0
    % [14]
    prkna = [ 10   55.437   0.25979
              15   85.131   0.21425
              20  115.114   0.14047
              25  145.235   0.09272
              30  175.444   0.06082
              35  205.721   0.04023
              40  236.051   0.02580
              45  266.424   0.01617
              50  296.832   0.01022 ];
end

%load('bendit03.mat');

for i = 1 : length(Rcurvs)
    for j = 1 : length(Modes{i})
        Neffs(j,i) = Modes{i}{j}.neff;
    end
end

nu = (2*pi*Rcurvs/lambda0) .* Neffs(1,:);
res = [];
for i = length(nu) : -1 : 1
    a = find(abs(prkna(:,1)'-Rcurvs(i)) < 10*eps);
    if ~isempty(a)
        a = a(1);
        relerr_r = abs(real(nu(i)) - prkna(a,2)) / prkna(a,2);
        relerr_i = abs(imag(nu(i)) - prkna(a,3)) / prkna(a,3);
        %        res = vertcat(res, [ PARs(i) real(nu(i)) imag(nu(i)) prkna(a,2) prkna(a,3) log10(relerr_r) log10(relerr_i) ]);
        res = vertcat(res, [ Rcurvs(i) real(nu(i)) imag(nu(i)) ]);
    end
end
format long g
res(:,1)
res(:,2:3)


