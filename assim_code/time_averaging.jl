function avg2(x)
if (typeof(x[1])==Float64) | (typeof(x[1])==Float32)
        return (x[1:2:end] + x[2:2:end])/2;
else
        return x[1:2:end]
end

end


function avg2_2d(x)

return (x[1:2:end,:] + x[2:2:end,:])/2;

end


function avg2_df(x)

ncol = Integer(size(x)[2]);
nrow = Integer(size(x)[1]/2);
ans = deepcopy(x[1:nrow,:]);

for j in 1:ncol
        ans[:,j] = avg2(x[:,j]);
end

return ans;

end

function convert_sim(x)
a = avg2_df(x[1]);
b = avg2_2d(x[2]);
c = avg2_2d(x[3]);
d = avg2_2d(x[4]);

return a,b,c,d
end

