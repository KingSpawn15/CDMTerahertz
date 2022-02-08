clear all;
close all;


parfor pdelay = 1 : 12
    delay = 5 * pdelay;
    videoFReader = VideoReader(['avis/delay=',num2str(delay),'.avi']);
    % Write Video
    videoFWrite = VideoWriter(['avis/delay=',num2str(delay),'.mp4'],'MPEG-4');
    open(videoFWrite);
    for count = 1:abs(videoFReader.Duration*videoFReader.FrameRate)
%         disp(count);
        key_frame = read(videoFReader,count);
        writeVideo(videoFWrite,key_frame);
    end
    %     Release video object
    %     close(videoFReader);
    close(videoFWrite);
    disp('COMPLETED');
    close all;
end