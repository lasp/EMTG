# Install script for directory: C:/SwInstalls/gsl-2.4.0

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files/GSL")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/SwInstalls/gsl-2.4.0/build/Debug/gsl.lib")
  elseif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/SwInstalls/gsl-2.4.0/build/Release/gsl.lib")
  elseif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/SwInstalls/gsl-2.4.0/build/MinSizeRel/gsl.lib")
  elseif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/SwInstalls/gsl-2.4.0/build/RelWithDebInfo/gsl.lib")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/SwInstalls/gsl-2.4.0/build/Debug/gslcblas.lib")
  elseif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/SwInstalls/gsl-2.4.0/build/Release/gslcblas.lib")
  elseif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/SwInstalls/gsl-2.4.0/build/MinSizeRel/gslcblas.lib")
  elseif("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/SwInstalls/gsl-2.4.0/build/RelWithDebInfo/gslcblas.lib")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gsl" TYPE FILE FILES
    "C:/SwInstalls/gsl-2.4.0/gsl_inline.h"
    "C:/SwInstalls/gsl-2.4.0/gsl_machine.h"
    "C:/SwInstalls/gsl-2.4.0/gsl_math.h"
    "C:/SwInstalls/gsl-2.4.0/gsl_minmax.h"
    "C:/SwInstalls/gsl-2.4.0/gsl_mode.h"
    "C:/SwInstalls/gsl-2.4.0/gsl_nan.h"
    "C:/SwInstalls/gsl-2.4.0/gsl_pow_int.h"
    "C:/SwInstalls/gsl-2.4.0/gsl_precision.h"
    "C:/SwInstalls/gsl-2.4.0/gsl_types.h"
    "C:/SwInstalls/gsl-2.4.0/gsl_version.h"
    "C:/SwInstalls/gsl-2.4.0/blas/gsl_blas.h"
    "C:/SwInstalls/gsl-2.4.0/blas/gsl_blas_types.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_char.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_complex_double.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_complex_float.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_complex_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_double.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_float.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_int.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_long.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_short.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_uchar.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_uint.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_ulong.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_block_ushort.h"
    "C:/SwInstalls/gsl-2.4.0/block/gsl_check_range.h"
    "C:/SwInstalls/gsl-2.4.0/bspline/gsl_bspline.h"
    "C:/SwInstalls/gsl-2.4.0/cblas/gsl_cblas.h"
    "C:/SwInstalls/gsl-2.4.0/cdf/gsl_cdf.h"
    "C:/SwInstalls/gsl-2.4.0/cheb/gsl_chebyshev.h"
    "C:/SwInstalls/gsl-2.4.0/combination/gsl_combination.h"
    "C:/SwInstalls/gsl-2.4.0/complex/gsl_complex.h"
    "C:/SwInstalls/gsl-2.4.0/complex/gsl_complex_math.h"
    "C:/SwInstalls/gsl-2.4.0/const/gsl_const.h"
    "C:/SwInstalls/gsl-2.4.0/const/gsl_const_cgs.h"
    "C:/SwInstalls/gsl-2.4.0/const/gsl_const_cgsm.h"
    "C:/SwInstalls/gsl-2.4.0/const/gsl_const_mks.h"
    "C:/SwInstalls/gsl-2.4.0/const/gsl_const_mksa.h"
    "C:/SwInstalls/gsl-2.4.0/const/gsl_const_num.h"
    "C:/SwInstalls/gsl-2.4.0/deriv/gsl_deriv.h"
    "C:/SwInstalls/gsl-2.4.0/dht/gsl_dht.h"
    "C:/SwInstalls/gsl-2.4.0/diff/gsl_diff.h"
    "C:/SwInstalls/gsl-2.4.0/eigen/gsl_eigen.h"
    "C:/SwInstalls/gsl-2.4.0/err/gsl_errno.h"
    "C:/SwInstalls/gsl-2.4.0/err/gsl_message.h"
    "C:/SwInstalls/gsl-2.4.0/fft/gsl_dft_complex.h"
    "C:/SwInstalls/gsl-2.4.0/fft/gsl_dft_complex_float.h"
    "C:/SwInstalls/gsl-2.4.0/fft/gsl_fft.h"
    "C:/SwInstalls/gsl-2.4.0/fft/gsl_fft_complex.h"
    "C:/SwInstalls/gsl-2.4.0/fft/gsl_fft_complex_float.h"
    "C:/SwInstalls/gsl-2.4.0/fft/gsl_fft_halfcomplex.h"
    "C:/SwInstalls/gsl-2.4.0/fft/gsl_fft_halfcomplex_float.h"
    "C:/SwInstalls/gsl-2.4.0/fft/gsl_fft_real.h"
    "C:/SwInstalls/gsl-2.4.0/fft/gsl_fft_real_float.h"
    "C:/SwInstalls/gsl-2.4.0/fit/gsl_fit.h"
    "C:/SwInstalls/gsl-2.4.0/histogram/gsl_histogram.h"
    "C:/SwInstalls/gsl-2.4.0/histogram/gsl_histogram2d.h"
    "C:/SwInstalls/gsl-2.4.0/ieee-utils/gsl_ieee_utils.h"
    "C:/SwInstalls/gsl-2.4.0/integration/gsl_integration.h"
    "C:/SwInstalls/gsl-2.4.0/interpolation/gsl_interp.h"
    "C:/SwInstalls/gsl-2.4.0/interpolation/gsl_interp2d.h"
    "C:/SwInstalls/gsl-2.4.0/interpolation/gsl_spline.h"
    "C:/SwInstalls/gsl-2.4.0/interpolation/gsl_spline2d.h"
    "C:/SwInstalls/gsl-2.4.0/linalg/gsl_linalg.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_char.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_complex_double.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_complex_float.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_complex_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_double.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_float.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_int.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_long.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_short.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_uchar.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_uint.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_ulong.h"
    "C:/SwInstalls/gsl-2.4.0/matrix/gsl_matrix_ushort.h"
    "C:/SwInstalls/gsl-2.4.0/min/gsl_min.h"
    "C:/SwInstalls/gsl-2.4.0/monte/gsl_monte.h"
    "C:/SwInstalls/gsl-2.4.0/monte/gsl_monte_miser.h"
    "C:/SwInstalls/gsl-2.4.0/monte/gsl_monte_plain.h"
    "C:/SwInstalls/gsl-2.4.0/monte/gsl_monte_vegas.h"
    "C:/SwInstalls/gsl-2.4.0/multifit/gsl_multifit.h"
    "C:/SwInstalls/gsl-2.4.0/multifit/gsl_multifit_nlin.h"
    "C:/SwInstalls/gsl-2.4.0/multifit_nlinear/gsl_multifit_nlinear.h"
    "C:/SwInstalls/gsl-2.4.0/multilarge/gsl_multilarge.h"
    "C:/SwInstalls/gsl-2.4.0/multilarge_nlinear/gsl_multilarge_nlinear.h"
    "C:/SwInstalls/gsl-2.4.0/multimin/gsl_multimin.h"
    "C:/SwInstalls/gsl-2.4.0/multiroots/gsl_multiroots.h"
    "C:/SwInstalls/gsl-2.4.0/multiset/gsl_multiset.h"
    "C:/SwInstalls/gsl-2.4.0/ntuple/gsl_ntuple.h"
    "C:/SwInstalls/gsl-2.4.0/ode-initval/gsl_odeiv.h"
    "C:/SwInstalls/gsl-2.4.0/ode-initval2/gsl_odeiv2.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permutation.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_char.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_complex_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_complex_float.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_complex_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_float.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_int.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_long.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_char.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_complex_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_complex_float.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_complex_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_float.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_int.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_long.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_short.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_uchar.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_uint.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_ulong.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_matrix_ushort.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_short.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_uchar.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_uint.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_ulong.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_ushort.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_char.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_complex_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_complex_float.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_complex_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_float.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_int.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_long.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_short.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_uchar.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_uint.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_ulong.h"
    "C:/SwInstalls/gsl-2.4.0/permutation/gsl_permute_vector_ushort.h"
    "C:/SwInstalls/gsl-2.4.0/poly/gsl_poly.h"
    "C:/SwInstalls/gsl-2.4.0/qrng/gsl_qrng.h"
    "C:/SwInstalls/gsl-2.4.0/randist/gsl_randist.h"
    "C:/SwInstalls/gsl-2.4.0/rng/gsl_rng.h"
    "C:/SwInstalls/gsl-2.4.0/roots/gsl_roots.h"
    "C:/SwInstalls/gsl-2.4.0/rstat/gsl_rstat.h"
    "C:/SwInstalls/gsl-2.4.0/siman/gsl_siman.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_heapsort.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_char.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_double.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_float.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_int.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_long.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_short.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_uchar.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_uint.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_ulong.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_ushort.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_char.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_double.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_float.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_int.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_long.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_short.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_uchar.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_uint.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_ulong.h"
    "C:/SwInstalls/gsl-2.4.0/sort/gsl_sort_vector_ushort.h"
    "C:/SwInstalls/gsl-2.4.0/spblas/gsl_spblas.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_airy.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_bessel.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_clausen.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_coulomb.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_coupling.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_dawson.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_debye.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_dilog.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_elementary.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_ellint.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_elljac.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_erf.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_exp.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_expint.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_fermi_dirac.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_gamma.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_gegenbauer.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_hermite.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_hyperg.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_laguerre.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_lambert.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_legendre.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_log.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_mathieu.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_pow_int.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_psi.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_result.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_synchrotron.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_transport.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_trig.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_sf_zeta.h"
    "C:/SwInstalls/gsl-2.4.0/specfunc/gsl_specfunc.h"
    "C:/SwInstalls/gsl-2.4.0/splinalg/gsl_splinalg.h"
    "C:/SwInstalls/gsl-2.4.0/spmatrix/gsl_spmatrix.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_char.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_double.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_float.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_int.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_long.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_short.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_uchar.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_uint.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_ulong.h"
    "C:/SwInstalls/gsl-2.4.0/statistics/gsl_statistics_ushort.h"
    "C:/SwInstalls/gsl-2.4.0/sum/gsl_sum.h"
    "C:/SwInstalls/gsl-2.4.0/sys/gsl_sys.h"
    "C:/SwInstalls/gsl-2.4.0/test/gsl_test.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_char.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_complex.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_complex_double.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_complex_float.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_complex_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_double.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_float.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_int.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_long.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_long_double.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_short.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_uchar.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_uint.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_ulong.h"
    "C:/SwInstalls/gsl-2.4.0/vector/gsl_vector_ushort.h"
    "C:/SwInstalls/gsl-2.4.0/wavelet/gsl_wavelet.h"
    "C:/SwInstalls/gsl-2.4.0/wavelet/gsl_wavelet2d.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "C:/SwInstalls/gsl-2.4.0/build/gsl.pc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES "C:/SwInstalls/gsl-2.4.0/build/gsl-config")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "C:/SwInstalls/gsl-2.4.0/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")