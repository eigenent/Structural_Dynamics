function dTdt = dynamics(t, T, u, time_u, A, B)

input_int = interp1(time_u, u, t); % interpolate the data set (time_u, u) at times t

dTdt = A*T + B*input_int; % evaluate ode at times t
end