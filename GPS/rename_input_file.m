files=dir([run_tag,'_*_kinematic.in']);
% files=dir([run_tag,'_*_',num2str(c2(1)),'_',num2str(c2(2)),'_',num2str(c2(3)),'.log']);
names = {files.name};

for i=1:length(names)
   oldfilename = names{i};
   newfilename = sprintf('%s%s.in',oldfilename(1:end-3),...
       ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_PID',num2str(PID)]);
   copyfile(oldfilename,newfilename);
end