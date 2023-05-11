cmake_minimum_required(VERSION 3.5)

project(HeuristicFlow_Include LANGUAGES CXX)

unset(heu_include_files)

foreach(md ${Heu_submodules})
    list(APPEND heu_include_files ${md})
    add_subdirectory(src/${md})
endforeach(md ${Heu_submodules})

if(${Heu_do_install})
    install(FILES ${heu_include_files}
        DESTINATION include/HeuristicFlow)
endif()