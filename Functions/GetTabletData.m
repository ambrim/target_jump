function tabdata = GetTabletData()
tabdata=[];
 %Fast loop to dump all data from the queue into a variable
 %pkt is 8 element vector. 1=xpos,2=ypos,6=timestamp
    while 1
        pkt = WinTabMex(5);

        if isempty(pkt) %ie the queue is empty
            break;
        else
            tabdata = [tabdata pkt];
        end

        status = uint32(pkt(7));
%         if bitget(status, 1) % Doesnt seem to work
%             disp('Pen is outside the tablet');
%             Screen('DrawText',window_ptr,'Keep Pen on Tablet!',start_x-100,start_y-50,red);
%             Screen('Flip',window_ptr, [], 1);
%         end
        if bitget(status, 2)
            disp('Queue overflow - packets getting lost');
        end
    end
end