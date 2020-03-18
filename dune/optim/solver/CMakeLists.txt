set( HEADERS
  fullydynmatrix.hh
  gauss.hh
  rowrefvector.hh
  gaussjordan.hh
)

install( FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/optim/solver )
