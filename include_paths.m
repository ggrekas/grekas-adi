function chooseMode(enableParaller, enableCompilation)
 
if(enableParaller)
   addpath('mexFiles/paraller/')
else
   addpath('mexFiles/serial/')
end

if(enableCompilation)
   compileMexFiles
end


