function(REDSET_ADD_TEST_SERIAL name args outputs)

  # job launcher
  if(NOT DEFINED VELOC_RESOURCE_MANAGER)
    set(test_param mpirun -np 2)
  elseif(${VELOC_RESOURCE_MANAGER} STREQUAL "NONE")
    set(test_param mpirun -np 2)
  elseif(${VELOC_RESOURCE_MANAGER} STREQUAL "LSF")
    set(test_param jsrun -r 1)
  elseif(${VELOC_RESOURCE_MANAGER} STREQUAL "SLURM")
    set(test_param srun -N 1 -n 1)
  endif()

  # Tests
  add_test(NAME ${name} COMMAND ${test_param} ./${name} ${args} )

endfunction(REDSET_ADD_TEST_SERIAL)

