% Now plot the data like in Fig. 4b of the paper. Complex mode is
% highlighted red.

plotdata1 = nan(5, length(fs));
plotdata2 = nan(5, length(fs));
plotdata3 = nan(5, length(fs));

for i = 1 : size(neffs,1)
    for j = 1 : size(neffs,2)
        n = neffs(i,j);
        if isreal(n)
            % real, propagating mode
            plotdata1(i,j) = n;
        else
            if real(n) == 0
                % imaginary, evanescent mode
                plotdata1(i,j) = -abs(n);
            else
                % complex mode (draw with dashed lines)
                plotdata2(i,j) = +abs(real(n));
                plotdata3(i,j) = -abs(imag(n));
            end
        end
    end
end

clf
plot(fs/1e9, plotdata1, 'k-');
hold on
plot(fs/1e9, plotdata2, 'r--', 'LineWidth', 2);
plot(fs/1e9, plotdata3, 'r--', 'LineWidth', 2);
xlabel 'Frequency / GHz'
ylabel '"Effective indices"'
grid on
