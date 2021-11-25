include("./assim_nov24.jl")

#function do_assim(obs_mask, norm_factor, outdir, vodA, vodB, vodC)

#trueTB = get_TB_2(sim_res1[2][:,1], tsoil, tcan, canpot_true, laiM, vodA, vodB, vodC);
#obsH = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
#obsV = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));

datalen = 24*365;
mask1 = zeros(datalen);
mask1[2:(24*3):end] .= 1;
mask1[(2+12):(24*3):end] .= 1;
mask1 = mask1 .== 1;

mask1o = zeros(datalen);
mask1o[2:(24*3):end] .= 1;
mask1o = mask1o .== 1;

mask6 = zeros(datalen);
mask6[7:(24*3):end] .= 1;
mask6[(7+12):(24*3):end] .= 1;
mask6 = mask6 .== 1;

maskAll = ones(datalen) .== 1;

vA, vB, vC = (0.067, 0.82, 0.051);

function run_scenario(scenario, sout)
  #sout = string(scenario,"/")
  if scenario == "o1AMPM"
    do_assim(mask1, 1, sout, vA, vB, vC) 
  

  elseif scenario == "o6AMPM"
    do_assim(mask6, 1, sout, vA, vB, vC)
  

  elseif scenario == "o1AM"
    do_assim(mask1o, 2, sout, vA, vB, vC)
  

  elseif scenario == "oAll"
    do_assim(maskAll, 1.0/36, sout, vA, vB, vC)
  end
end

run_scenario(ARGS[1],ARGS[2])

