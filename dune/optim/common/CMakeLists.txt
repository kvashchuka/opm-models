set( HEADERS
  densesubvector.hh
  fieldmatrix.hh
  integerrange.hh
  matrix.hh
  smallobject.hh
)

install( FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/optim/common )
