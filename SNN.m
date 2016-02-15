%Initializes Parameters:
n = 1000; %number of neurons
k_k_inh = rand(n, 1) < 0.2; %the inhibitor nerons
k_exc = ~k_inh; %the excitory neurons
d = 8*k_exc + 2*k_inh; %the d constant for inhibitory and excitory neurons
dt = 0.5; %delta time in ms
a = 0.02; %the a constant for calculating du
b = 0.2; %the b constant for calculating du
c = -65; %resetting the neurons
T = ceil(1000/dt); %The time axis
Tg = 10; %Big Tau - g

n_in = 100;% number of inhibitory fires/0.5*milisecond
frate = 0.002; %firing rate
g_in = zeros(n_in, 1); %g constanat coming in
E_in = zeros(n_in, 1); %E constant coming in
w_in = ones(n, n_in)*0.07; %weightage
w_in(rand(n, n_in)>0.1) = 0; %only 10% of incoming fires are actually affecting the triggered neurons
%the recurrent spiked network
g = zeros(n, 1); 
E = zeros(n, 1);
E(k_inh) = -85;
W = zeros(n, n);
idx = find(rand(n, n) < 0.1);
W(idx) = gamrnd(2,.003,length(idx),1);
W(k_exc, k_inh) = 2*W(k_exc, k_inh);
W = sparse(W);
fn = 0; %neurons fired
%The voltages for each 0.5 ms
v = zeros(n, T);
u = zeros(n, T);
%initial values
v(:, 1) = -70;
u(:, 1) = -14;

%as time passes
for time = 2:T
	if (time*dt >= 200) && (time*dt <= 700) %it starts from 200ms and ends at 700ms
		p = rand(n_in, 1) < frate*dt;
	else
		p = 0;
	end
	%the randomized input form other neurons
	g_in = g_in + p; 
	Iapp = w_in*(g_in.*E_in);
	Iapp = Iapp - (w_in*g_in).*v(:, time-1);
	g_in = (1-dt/Tg)*g_in;
	%The recurrent part of it. Each neuron helps trigger itself
	g = g + fn;
	Isyn = W*(g.*E) - (W*g).*v(:, time-1);
	Iapp = Iapp + Isyn;
	g = (1 - dt/Tg)*g;
	%updating the voltages
	v(:, time) = v(:, time-1) + dt*((0.04*v(:, time-1) + 5).*v(:, time - 1) - u(:, time - 1) + 140 + Iapp);
	u(:, time) = u(:, time-1) + dt*a*(b*v(:, time-1) - u(:, time-1));
	%fired neurons
	fn = v(:, time) >= 35;
	v(fn, time) = c;
	u(fn, time) = u(fn, time - 1) + d(fn);
	v(fn, time-1) = 35;
end
%Plotting it, if it spikes then red, if it's an inhibitory neuron then black
spikes = double(v==35);
spikes(k_inh, :) = 2*spikes(k_inh, :);
clf, hold on;
[X, Y] = meshgrid((0:T-1)*dt, 1:n);
col = 'kr';
for k = 1:2 
	idx = find(spikes == k);
	plot(X(idx), Y(idx), [col(k) '.']);
end
xlim([0, T*dt]); ylim([0, n])
xlabel('Time ms');
ylabel('Unit #');