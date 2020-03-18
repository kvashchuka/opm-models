set( HEADERS
  combinedgauss.hh
  gradient.hh
  reducedgauss.hh
)

install( FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/optim/qp )
