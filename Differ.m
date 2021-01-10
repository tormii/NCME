function ybar = Differ(y,delta)
% compute the derivative of a discrete time series y
% delta denotes the sampling time interval of y
L = length(y);
ybar = zeros(1,L-2);%
for i = 2 : L-1
  ybar(i-1)=(y(i+1)-y(i-1))/(2*delta);
end
ybar = [(y(2)-y(1))/delta,ybar,(y(end)-y(end-1))/delta];
end