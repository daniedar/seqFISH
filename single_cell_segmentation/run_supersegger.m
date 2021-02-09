function run_supersegger(dirname)
    dirname = char(dirname)    
    disp(dirname)
    CONST = loadConstants('60XPa',0)
    CONST.parallel.PARALLEL_FLAG = 0
    CONST.superSeggerOpti.segmenting_fluorescence = 1
    clean_flag = 0
    BatchSuperSeggerOpti(dirname,1,clean_flag,CONST,[2,3])
    close all;

end

