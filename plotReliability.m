function plotReliability(filename)

datab = tdfread(filename);

fcount = 0;
for j=1:1:numel(datab.I1); if datab.I1(j)==datab.I1(1); fcount = fcount + 1; else break; end; end;

color_res = 100;
rb = datab.reliability;
rb = rb .* (rb > 0); % clamp at 0
rb = sqrt(rb);
rb = real(rb/max(rb)) * color_res; % map onto color spectrum

colormap(transpose([1:-1/color_res:0; 1:-1/color_res:0; 1:-1/color_res:0]));
image([datab.freq(1) datab.freq(end)], [datab.I1(1) datab.I1(end)], vec2mat(rb, fcount));
ylabel('I_1');
xlabel('f');
hcb = colorbar;
set(gca, 'YDir', 'normal');
set(hcb,'YTickLabel',0.1:0.1:1);

end