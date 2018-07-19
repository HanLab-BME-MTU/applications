function diaryFlush

    if strcmp( get(0, 'Diary'), 'on' )
        diary off;
        diary on;
    end
    
end