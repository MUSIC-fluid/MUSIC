add_executable (mpihydro
    main.cpp
    grid.cpp
    eos.cpp
    evolve.cpp
    init.cpp
    reconst.cpp
    minmod.cpp
    advance.cpp
    u_derivative.cpp
    dissipative.cpp
    util.cpp
    grid_info.cpp
    read_in_parameters.cpp
    music.cpp
    cornelius.cpp
    )

target_link_libraries (mpihydro ${GSL_LIBRARIES})

install(TARGETS mpihydro DESTINATION ${CMAKE_HOME_DIRECTORY})
