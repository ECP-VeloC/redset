# Header files need not be enumerated
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR})

###########
# LIB REDSET #
###########

# Install header files
LIST(APPEND libredset_install_headers
    redset.h
)
INSTALL(FILES ${libredset_install_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

LIST(APPEND redset_noMPI_srcs
    redset_err_serial.c
    redset_util.c
    redset_io.c
    redset_lofi.c
    redset_partner_serial.c
    redset_xor_serial.c
    redset_reedsolomon_common.c
    redset_reedsolomon_serial.c
)

LIST(APPEND libredset_srcs
    redset_err_mpi.c
    redset_util.c
    redset_util_mpi.c
    redset_io.c
    redset_lofi.c
    redset.c
    redset_single.c
    redset_partner.c
    redset_reedsolomon_common.c
    redset_reedsolomon.c
)

IF(ENABLE_CUDA)
    LIST(APPEND libredset_srcs
        redset_xor.cu
    )
ELSE(ENABLE_CUDA)
    LIST(APPEND libredset_srcs
        redset_xor.c
    )
ENDIF(ENABLE_CUDA)

IF(BUILD_SHARED_LIBS)
   # REDSET Library
   ADD_LIBRARY(redset_o OBJECT ${libredset_srcs})
   TARGET_LINK_LIBRARIES(redset_o PRIVATE ${REDSET_EXTERNAL_LIBS})
   TARGET_COMPILE_DEFINITIONS(redset_o PRIVATE REDSET_ENABLE_MPI)

   ADD_LIBRARY(redset SHARED $<TARGET_OBJECTS:redset_o>)
   ADD_LIBRARY(redset::redset ALIAS redset)
   TARGET_LINK_LIBRARIES(redset PUBLIC ${REDSET_EXTERNAL_LIBS})
   TARGET_INCLUDE_DIRECTORIES(redset PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
   TARGET_COMPILE_DEFINITIONS(redset PUBLIC REDSET_ENABLE_MPI)
   SET_TARGET_PROPERTIES(redset PROPERTIES OUTPUT_NAME redset CLEAN_DIRECT_OUTPUT 1)
   INSTALL(TARGETS redset EXPORT redsetTargets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF()

# REDSET Library
ADD_LIBRARY(redset_static_o OBJECT ${libredset_srcs})
TARGET_LINK_LIBRARIES(redset_static_o PRIVATE ${REDSET_EXTERNAL_STATIC_LIBS})
TARGET_COMPILE_DEFINITIONS(redset_static_o PRIVATE REDSET_ENABLE_MPI)

ADD_LIBRARY(redset-static STATIC $<TARGET_OBJECTS:redset_static_o>)
ADD_LIBRARY(redset::redset-static ALIAS redset-static)
IF(REDSET_LINK_STATIC)
   SET_TARGET_PROPERTIES(redset-static PROPERTIES LINK_SEARCH_START_STATIC 1)
   SET_TARGET_PROPERTIES(redset-static PROPERTIES LINK_SEARCH_END_STATIC 1)
ENDIF(REDSET_LINK_STATIC)
TARGET_LINK_LIBRARIES(redset-static PUBLIC ${REDSET_EXTERNAL_STATIC_LIBS})
TARGET_INCLUDE_DIRECTORIES(redset-static PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
TARGET_COMPILE_DEFINITIONS(redset-static PUBLIC REDSET_ENABLE_MPI)

SET_TARGET_PROPERTIES(redset-static PROPERTIES OUTPUT_NAME redset CLEAN_DIRECT_OUTPUT 1)
INSTALL(TARGETS redset-static EXPORT redsetTargets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Non-MPI library for CLI
if (BUILD_SHARED_LIBS)
   ADD_LIBRARY(redset_noMPI_o OBJECT ${redset_noMPI_srcs})
   TARGET_LINK_LIBRARIES(redset_noMPI_o PRIVATE ${REDSET_SERIAL_LIBS})

   ADD_LIBRARY(redset_base SHARED $<TARGET_OBJECTS:redset_noMPI_o>)
   ADD_LIBRARY(redset::redset_base ALIAS redset_base)
   TARGET_LINK_LIBRARIES(redset_base PUBLIC ${REDSET_SERIAL_LIBS})
   TARGET_INCLUDE_DIRECTORIES(redset_base PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
   INSTALL(TARGETS redset_base EXPORT redsetTargets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

ADD_LIBRARY(redset_noMPI_static_o OBJECT ${redset_noMPI_srcs})
TARGET_LINK_LIBRARIES(redset_noMPI_static_o PRIVATE ${REDSET_SERIAL_STATIC_LIBS})

ADD_LIBRARY(redset_base-static STATIC $<TARGET_OBJECTS:redset_noMPI_static_o>)
ADD_LIBRARY(redset::redset_base-static ALIAS redset_base-static)
TARGET_LINK_LIBRARIES(redset_base-static PUBLIC ${REDSET_SERIAL_STATIC_LIBS})
TARGET_INCLUDE_DIRECTORIES(redset_base-static PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
SET_TARGET_PROPERTIES(redset_base-static PROPERTIES OUTPUT_NAME redset_base CLEAN_DIRECT_OUTPUT 1)
INSTALL(TARGETS redset_base-static EXPORT redsetTargets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

##############
# CLI REDSET #
##############

## Individual binaries generated from a single .c file
#LIST(APPEND redset_c_bins
#	redset_rebuild_xor
#)
#
## Build and install C binaries
#FOREACH(bin IN ITEMS ${redset_c_bins})
#	ADD_EXECUTABLE(${bin} ${bin}.c)
#	TARGET_LINK_LIBRARIES(${bin} redset_base)
#	INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${bin} DESTINATION ${CMAKE_INSTALL_BINDIR})
#ENDFOREACH(bin IN ITEMS ${redset_c_bins})
