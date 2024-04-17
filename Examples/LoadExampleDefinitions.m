function Examples = LoadExampleDefinitions()

% (c) Sui Tang -2021-05-19

def_files                          = dir('./Examples/*_def.m');   
total_num_defs                     = length(def_files);
Examples                           = cell(1, total_num_defs);
for idx = 1 : total_num_defs
  eval(sprintf('Example = %s();', erase(def_files(idx).name, '.m')));
  Examples{idx}                    = Example;
end

end