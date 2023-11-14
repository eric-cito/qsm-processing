function [latex_cell] = ScreenShotLatexCell(exam_id, QCfolder)

exam_id = strrep(exam_id,'_','-');

latex_cell = {};

QSMfile = dir(QCfolder);
QSMfile([QSMfile.isdir]) = [];
% QSMfile(cellfun(@isempty,strfind({QSMfile.name},'QSM_'))) = [];
QSMfile(cellfun(@isempty,strfind({QSMfile.name},'.png'))) = [];

latex_cell(end+1) = {['\begin{figure}[!h]'  ]};
for ii = 1:length(QSMfile)
    ImgFile = [QCfolder '/' QSMfile(ii).name];
    QSMname = strrep(QSMfile(ii).name,'.png','');
    QSMname = strrep(QSMname,'QSM_','');
    QSMname = strrep(QSMname,'_','-');
    latex_cell(end+1) = {['\begin{minipage}{.33\textwidth}  ']};
    latex_cell(end+1) = {['\includegraphics[width=0.8\linewidth]{' ImgFile '}  ']};
    latex_cell(end+1) = {['\caption{' exam_id '-' QSMname '}  ']};
    latex_cell(end+1) = {['\end{minipage}  ']};
    
    if mod(ii,3) == 0 && ii ~= length(QSMfile)
        latex_cell(end+1) = {['\end{figure}'  ]};
        latex_cell(end+1) = {['\begin{figure}  ']};
    end
end
latex_cell(end+1) = {['\end{figure}  ']};
latex_cell(end+1) = {['\clearpage  ']};

end