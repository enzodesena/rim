function [s, f_mesh,t_mesh] = plot_spectrogram(rir, spectr_window_length, spectr_overlap, spectr_nfft, FS)

[s,f,t] = spectrogram(rir, spectr_window_length, spectr_overlap, spectr_nfft, FS);
s = s(2:end, :);
s = s./max(max(abs(s)));
[f_mesh, t_mesh] = meshgrid(f(2:end),t);
mesh(f_mesh,t_mesh,20*log10(abs(s')))
view(-90,90)
set(gca,'ydir','reverse')
colormap(gray)
ylabel('Time [s]')
xlabel('Frequency [Hz]')
colorbar
