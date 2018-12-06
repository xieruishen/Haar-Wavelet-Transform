function signal = InvHaar1D(a,d,level)
    signal = zeros(length(d{1})*2);
    index = length(signal);
    signallength = length(d{1});
    a_next = a;
    d_next = d;
    for i = length(d):(level+1)
        a_next = vertcat(1/sqrt(2)*(a_next - d_next),1/sqrt(2)*(a_next - d_next);
    end
end