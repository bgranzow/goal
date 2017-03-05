#ifndef GOAL_DATA_TYPES_HPP
#define GOAL_DATA_TYPES_HPP

/** \file goal_data_types.hpp */

#include <MatrixMarket_Tpetra.hpp>
#include <Phalanx_KokkosDeviceTypes.hpp>
#include <Sacado_Fad_SLFad.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

/** \brief All goal symbols are contained in this namespace. */
namespace goal {

/** \brief Local ordinal data type. */
typedef int LO;

/** \brief Global ordinal data type. */
typedef long long GO;

/** \brief The scalar data type. */
typedef double ST;

/** \brief Forward automatic differentiation type. */
typedef Sacado::Fad::SLFad<ST, GOAL_FAD_SIZE> FadType;

/** \brief Parallel communication type. */
typedef Teuchos::Comm<int> Comm;

/** \brief Kokkos node device type. */
typedef Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device> KNode;

/** \brief Tpetra map type. */
typedef Tpetra::Map<LO, GO, KNode> Map;

/** \brief Tpetra Crs graph type. */
typedef Tpetra::CrsGraph<LO, GO, KNode> Graph;

/** \brief Tpetra exporter type. */
typedef Tpetra::Export<LO, GO, KNode> Export;

/** \brief Tpetra importer type. */
typedef Tpetra::Import<LO, GO, KNode> Import;

/** \brief Tpetra vector type. */
typedef Tpetra::Vector<ST, LO, GO, KNode> Vector;

/** \brief Tpetra Crs Matrix type. */
typedef Tpetra::CrsMatrix<ST, LO, GO, KNode> Matrix;

/** \brief Tpetra matrix output class type. */
typedef Tpetra::MatrixMarket::Writer<Matrix> MM_Writer;

}  // namespace goal

#endif
