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
d2 = avg2(x[5]);

return a,b,c,d,d2
end


function get_diurnal(x,navg)
	if (typeof(x[1])==Float64) | (typeof(x[1])==Float32)
			return mean(reshape(x,(navg,:)),dims=2)[:,1];
	else
			return x[1:navg]
	end
end


function get_diurnal_summer(x,navg)
	if (typeof(x[1])==Float64) | (typeof(x[1])==Float32)
		day_reshape = reshape(x,(navg,:));
		day_list = collect(1:(size(day_reshape)[2])) .% 365;
		summer_mask = (day_list .>= 150) .* (day_list .< 250);
		summer_only = day_reshape[:,summer_mask];
		return mean(summer_only,dims=2)[:,1];
	else
			return x[1:navg]
	end
end


function get_daily(x,navg)
	if (typeof(x[1])==Float64) | (typeof(x[1])==Float32)
			return mean(reshape(x,(navg,:)),dims=1)[1,:];
	else
			return x[1:navg:end]
	end
end

function get_daily2d(x,navg)
	return mean(reshape(x,(navg,:,size(x)[2])),dims=1)[1,:,:];
end
