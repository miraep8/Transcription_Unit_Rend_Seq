fname_3f = 'C:/Users/mirae/Documents/Research/Data/wig Files/Jean_3f.wig';
reads_3f = openWig(fname_3f);
z_score_3f = zScores2(reads_3f);
wig_track_3f = makeVertibiWig(z_score_3f);

num_Peaks = sum(wig_track_3f(2,:) > 5)
new_name = strcat(fname_3f(1:strfind(fname_3f, '.wig')-1), '_peaks.wig');
writeWig(wig_track_3f, new_name);
%writeWig(z_score_3f, 'C:/Users/mirae/Documents/Research/Data/wig Files/Jean_3f_zscore.wig');
fname_3r = 'C:/Users/mirae/Documents/Research/Data/wig Files/Jean_3r.wig';
reads_3r = openWig(fname_3r);
z_score_3r = zScores(reads_3r);
wig_track_3r = makeVertibiWig(z_score_3r);
new_name = strcat(fname_3r(1:strfind(fname_3r, '.wig')-1), '_peaks.wig');
writeWig(wig_track_3r, new_name);

fname_5f = 'C:/Users/mirae/Documents/Research/Data/wig Files/Jean_5f.wig';
reads_5f = openWig(fname_5f);
z_score_5f = zScores(reads_5f);
wig_track_5f = makeVertibiWig(z_score_5f);
new_name = strcat(fname_5f(1:strfind(fname_5f, '.wig')-1), '_peaks.wig');
writeWig(wig_track_5f, new_name);

fname_5r = 'C:/Users/mirae/Documents/Research/Data/wig Files/Jean_5r.wig';
reads_5r = openWig(fname_5r);
z_score_5r = zScores(reads_5r);
wig_track_5r = makeVertibiWig(z_score_5r);
new_name = strcat(fname_5r(1:strfind(fname_5r, '.wig')-1), '_peaks.wig');
writeWig(wig_track_5r, new_name);
