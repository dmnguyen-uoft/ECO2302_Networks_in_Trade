function out = funscript(x,y,z)

    out1 = x.^2 + .5 .* y.^3 - sqrt(z);
    out2 = x+y;
    
    out.out1 = out1;
    out.out2 = out2;
    
end