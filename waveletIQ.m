function H = waveletIQ(signal, nlevels, wname, nbins, wlen, wstep, entbc)

% Wavelet information quality, per Slobounov, Cao, & Sebastianelli (2009)

% Inputs
%   signal: vector
%   nlevels: number of DWT decomposition levels
%   wname: name of DWT wavelet to use (input to Matlab's wavedec function)
%   nbins: number of bins for entropy calculations
%   slen: sliding window length in units of datapoints
%   step: sliding window step size in units of datapoints
%   entbc: entropy bias correction technique
%       'none': no bias correction (use plug-in estimate)
%       'mm': Miller-Madow
%       'jk': jackknifed
% Output
%   H: wavelet information quality

% Calculate number of windows
nsteps = 1 + ceil((numel(signal) - wlen) / wstep);
% Initialize vector of entropy estimates
H_estimates = NaN(1, nsteps);

for step = 1:nsteps
    % Get current window of data
    wstart = 1 + (step - 1)*wlen;
    wend = min(wstart + wlen - 1, numel(signal));
    subsignal = signal(wstart:wend);

    % Get DWT coefficients
    [coeffs, lens] = wavedec(subsignal, nlevels, wname);
    % Discard approximation coefficients from the final decomposition level
    coeffs = coeffs((lens(1)+1):end);

    % Estimate entropy
    N = numel(coeffs); % Number of datapoints
    switch lower(entbc)
        case 'none' % Plug-in entropy estimate
            curr_H = ent(coeffs, nbins);
        case 'mm' % Miller-Madow entropy estimate
            [curr_H, p] = ent(coeffs, nbins);
            curr_H = curr_H - (nnz(p) - 1)/2/N;
        case 'jk' % Jackknifed entropy estimate
            jk_H = 0;
            for j = 1:numel(coeffs)
                % Remove 1 datapoint
                idx = true(1, numel(coeffs));
                idx(j) = false;
				% Use plug-in estimate 
                jk_H = jk_H + ent(coeffs(idx), nbins);
            end
            curr_H = N*ent(coeffs, nbins) - (N - 1)/N*jk_H;    
    end
	
	% Append current estimate
	H_estimates(step) = curr_H;
end

H = mean(H_estimates);

end

function [H, p] = ent(data, nbins)

% Plug-in entropy estimate

counts = histcounts(data, nbins);
nzidx = counts ~= 0;
p = counts / sum(counts);
H = -sum(p(nzidx).*log(p(nzidx)));

end