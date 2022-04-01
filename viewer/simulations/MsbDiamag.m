function varargout = MsbDiamag(Sys, Exp, Opt)

X = linspace(Exp.Range(1), Exp.Range(2), Exp.nPoints).';

if length(Sys.lw)>1

    Y_lrz = Sys.lw(1)/2/pi*( 1./( (X-Sys.IS-Sys.Deq/2).^2+Sys.lw(1)^2/4  )...
                            +1./( (X-Sys.IS+Sys.Deq/2).^2+Sys.lw(1)^2/4  ) );
    if Sys.lw(2)~=0
        glw = Sys.lw(2)/(X(2)-X(1));
        Y_gau = 1/glw/sqrt(2*pi)*( exp(-([1:Exp.nPoints]-Exp.nPoints/2-1).^2/(2*glw^2)));
        
        YG = ifft(fftshift(Y_gau.')); YL = ifft(Y_lrz);
        Y = -abs(fft(YG.*YL));
    else
        Y = Y_lrz;
    end
   
else
    Y = -Sys.lw/2/pi*( 1./( (X-Sys.IS-Sys.Deq/2).^2+Sys.lw^2/4  ) +1./( (X-Sys.IS+Sys.Deq/2).^2+Sys.lw^2/4  ) );
    
end

Y = -Y/sum(Y)/(X(2)-X(1));

if ~nargout
    figure(12); clf; plot(X, Y);
else
    varargout{1} = X;
    varargout{2} = Y;
end
