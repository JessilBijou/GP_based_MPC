function indexes = logspacing(n_max,points_in_interval)
%LOGSPACING
%   Prepares a vector that assures that in each zone of the logaritmic plot
%   we have the same number of points:

% Maximum exponent:
n = floor(log10(n_max));

% Initialize:
indexes = [];

for i = 1:n
    if (10^i-10^(i-1))>points_in_interval
        indexes = [indexes linspace(10^(i-1)+1,10^i,points_in_interval)];
    else
        indexes = [indexes 10.^(i-1):10^i];
    end
end

% Turn into integers:
indexes = round(indexes);

end