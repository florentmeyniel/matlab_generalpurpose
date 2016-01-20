function [p_cor] = MCP_fromPval_fn(p, a, method)
% Compute FWE and FDR p-value at the a level on the provided p values.
% 
%   method       FWE control   adaptative  Assumption
% 'Bonf_pure'      strong         no         none   
% 'Bonf_ind'       strong         no         independant samples
% 'Simes'        weak (FDR)    step-up       PRDS
% 'FDRnoprior'   weak (FDR)    step-up       none
% 'Hochberg'       strong      step-up       PRDS
% 'Sidak'          strong      step-down     Slepian / Dunn-Sidak inequality
% 'Holm'           strong      step-down     none
%
% PRDS: 'Positive Regression Dependency on Subsets'. True with gaussian
% noise + non-null correlation on the null (i.e. true in most imaging data)
%
% Slepian / Dunn-Sidak inequality: (one-side test / two-side test). True for 
% gaussian data with positive correlation. 
%
% Usage: [p_cor] = MCP_fromPval_fn(p, a, method)
%           p: the p-value (any dimension is fine)
%           a: the alpha level for the multiple comparision correction
%      method: the method to use. If not specified, all the corrections
%              are returned in a structure. 
%
% Terminology and formula from:
%   Thomas Nichols & Satoru Hayasaka
%    Controlling the familywise error rate in functional neuroimaing: A
%    comparative review
%    Statistical Methods in Medical Research 2003; 12: 419-446
%
% For more details on FDR:
%   Christopher Genovese, Nicole Lazar, Thomas Nichols
%    Thresholding of statistical maps in functional neuroimaging using the
%    false discovery rate
%    NeuroImage 2003; 15, 870-878
%
% Florent Meyniel 2012-02-09

V = length(p(:));

% DATA INDEPENDANT
% ================
p_cor.bonf    = a / V;
p_cor.bonfind = 1 - (1-a)^(1/V);

% ADAPTATIVE STEP UP
% ==================
p = sort(p(:));
I = (1:V)';
c = sum(1./(1:V));

p_cor.simes      = p(find(p <= a*I/V, 1, 'last'));
p_cor.FDRnoprior = p(find(p <= a*I/(V*c), 1, 'last'));
p_cor.hochberg   = p(find( p<= a./(1 + V - I), 1, 'last'));

% ADAPTATIVE STEP DOWN
% =====================
p = sort(p(:), 1, 'descend');
I = flipud((1:V)');

p_cor.sidak = p(find( p <= (1 - (1-a).^(1./(V-I+1))), 1, 'first'));
p_cor.holm  = p(find( p <= a./(V-I+1), 1, 'first'));
        

if exist('method', 'var') ...
    && any(strcmpi(method, {'Bonf_pure', 'Bonf_ind', 'Simes',...
   'FDRnoprior', 'Hochberg', 'Sidak', 'Holm'}))
    switch lower(method)
        case lower('Bonf_pure')  ; p_cor = p_cor.bonf;     
        case lower('Bonf_ind')   ; p_cor = p_cor.bonfind;
        case lower('Simes')      ; p_cor = p_cor.simes;
        case lower('FDRnoprior') ; p_cor = p_cor.FDRnoprior;
        case lower('Hochberg')   ; p_cor = p_cor.hochberg;
        case lower('Sidak')      ; p_cor = p_cor.sidak;
        case lower('Holm')       ; p_cor = p_cor.holm;
    end
end        
   
end

