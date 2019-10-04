
if exist('Events.nev','file')
    currentLoc = pwd;
    filePath = [currentLoc,'/Events.nev'];
    [eventTimes, eventList] = Nlx2MatEV_v3(filePath,[1 0 0 0 1],0,1,[]);
    save('events.mat')
end