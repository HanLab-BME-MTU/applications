function currentTime = LG_getCurrentTime
%LG_getCurrentTime returns the current time in labelgui2

naviHandles = LG_getNaviHandles;
currentTime = round(get(naviHandles.LG_navi_timepointSlider_sli, 'Value'));