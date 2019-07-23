function [Ys] = randn1d(J, Q, fwhm)

% Modified from randomtalk.m by M.Brett and F.Turkheimer (Oct.1999) 
% Downloaded on 2014.06.27 from:
%     http://www.fil.ion.ucl.ac.uk/~wpenny/mbi/index.html


% generate Gaussian kernel:
sd    = fwhm/sqrt(8*log(2));   % sigma for this FWHM
x     = -0.5*(Q-1) : 0.5*(Q-1);
gf    = exp(-(x.*x)/(2*sd*sd));
gf    = gf/sum(sum(gf));

% variance expectation for this kernel
AG    = fft(gf);
Pag   = AG.*conj(AG);     % Power of the noise
COV   = real(ifft2(Pag));
svar  = COV(1,1);
scale = sqrt(1/svar);


% generate random data and convolve with Gaussian kernel:
Y     = randn(J, Q);
Ys  = zeros(J, Q);
for i=1:J
    Ys(i,:) = conv(Y(i,:), gf, 'same');
end
Ys = scale*Ys;
