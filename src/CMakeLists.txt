add_executable (HBT.e
  main.cpp
  readindata.cpp
  HBT.cpp
  arsenal.cpp
  gauss_quadrature.cpp
  ParameterReader.cpp
  Table.cpp
  )
target_link_libraries (HBT.e ${LIBS})

install(TARGETS HBT.e DESTINATION ${CMAKE_HOME_DIRECTORY})
