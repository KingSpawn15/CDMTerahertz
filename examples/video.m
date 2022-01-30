clearvars;
dir = 'results/combination/exhaustive/';

% video_out = VideoWriter('yourvideo.mp4','MPEG-4'); %create the video object
% video_out.Quality = 10;
% open(video_out);

text_str = cell(3,1);
position = [[0.5,1500];[1000,1500];[0 -20.5]];
box_color = {'red','green','yellow'};


for delay = -55 : 5 : 0
    
    video_out = VideoWriter(['delay=' num2str(delay) '.mp4'],'MPEG-4'); %create the video object
    video_out.Quality = 10;
    open(video_out);
    
    for interaction_gain_factor_rectification = -0.1:0.01:0.15
        for interaction_gain_factor_photodember =-0.1:0.01:0.15
            
            text_str{1} = ['Factor PD: ' num2str(interaction_gain_factor_photodember,'%0.2f')];
            text_str{2} = ['Factor OR: ' num2str(interaction_gain_factor_rectification,'%0.2f')];
            text_str{3} = ['delay: ' num2str(delay)];
            
            str = ['results/combination/exhaustive/','combination_',...
                'pd_gain=',num2str(interaction_gain_factor_photodember),...
                'or_gain=',num2str(interaction_gain_factor_rectification),...
                'delay=',num2str(delay)];
            
            img = imread(strcat(str,'.png')); %read the next image
            
            img = insertText(img,position,text_str,'FontSize',100,'BoxColor',...
                box_color);
            writeVideo(video_out,img);
            
        end
    end
    
    close(video_out);
    
end





