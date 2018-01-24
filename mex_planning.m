function mex_planning()
%MEX_PLANNING Build planning mex libs

cpp_dir = fullfile(fileparts(mfilename('fullpath')), 'private');

for f = {'plan_rays_greedy.cpp'}
    mex('-largeArrayDims', ... %'-v', ...
        'COPTIMFLAGS=-O3 -DNDEBUG -fopenmp', ...
        'CXXOPTIMFLAGS=-O3 -DNDEBUG -fopenmp', ...
        omp_lib(), ... % '-lgomp', '-liomp5'
        '-lboost_system', '-lboost_thread', '-lboost_chrono', ...
        ['-I' cpp_dir], fullfile(cpp_dir, f{1}), '-outdir', cpp_dir);
end

end

function flags = omp_lib()
    [~, cmdout] = system('g++ -liomp5');
    if ~isempty(strfind(cmdout, 'main'))
        flags = '-liomp5';
    else
        flags = '-lgomp';
    end
end
