fname_3f = 'C:/Users/mirae/Documents/Research/Data/wig Files/Jean_3f.wig';
reads_3f = openWig(fname_3f);
z_score_3f = zScores(reads_3f);
wig_track_3f = makeVertibiWig(z_score_3f);

num_Peaks = sum(wig_track_3f(2,:) > 5)
new_name = strcat(fname_3f(1:strfind(fname_3f, '.wig')-1), '_peaks.wig');
writeWig(wig_track_3f, fname_3f);
writeWig(z_score_3f, 'C:/Users/mirae/Documents/Research/Data/wig Files/Jean_3f_zscore.wig');
%fname_3r = 'C:/Users/mirae/Documents/Research/Data/wig Files/Jean_3r.wig';
%reads_3r = openWig(fname_3r);
%z_score_3r = zScores(reads_3r);
%wig_track_3r = makeVertibiWig(z_score_3r);
%writeWig(wig_track_3r, fname_3r);
%
%fname_5f = 'C:/Users/mirae/Documents/Research/Data/wig Files/Jean_5f.wig';
%reads_5f = openWig(fname_5f);
%z_score_5f = zScores(reads_5f);
%wig_track_5f = makeVertibiWig(z_score_5f);
%writeWig(wig_track_5f, fname_5f);
%
%fname_5r = 'C:/Users/mirae/Documents/Research/Data/wig Files/Jean_5r.wig';
%reads_5r = openWig(fname_5r);
%z_score_5r = zScores(reads_5r);
%wig_track_5r = makeVertibiWig(z_score_5r);
%writeWig(wig_track_5r, fname_5r);
%










% z_scores = -2:.5:5;
% expected_num = length(z_score_3f(2,:))*(1 - normcdf(z_scores));
% semilogy(z_scores, expected_num)
% hold on
% seen = zeros(1, length(z_scores));
% for i = 1:length(seen)
%     seen(1,i) = sum(z_score_3f(2,:) > z_scores(i));
% end
% plot(z_scores, seen)
% plot(z_scores, length(z_score_3f)*ones(1,length(z_scores)))
% legend('Expected', 'Observed', 'Max Val')
% xlabel('Z score') 
% xticks(min(z_scores):max(z_scores))
% ylabel('Number of Reads with a Higher Z Score') 
% hold off
