files=dir([run_tag,'_toplog_*_kinematic.tmp']);
% files=dir([run_tag,'_*_',num2str(c2(1)),'_',num2str(c2(2)),'_',num2str(c2(3)),'.log']);
names = {files.name};

for i=1:length(names)
   oldfilename = names{i};
   newfilename = sprintf('%s%s.log',oldfilename(1:end-4),...
       ['_',num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_PID',num2str(PID)]);
   movefile(oldfilename,newfilename);
end