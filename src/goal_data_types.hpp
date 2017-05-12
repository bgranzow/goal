#ifndef goal_data_types_hpp
#define goal_data_types_hpp

/// @file goal_data_types.hpp

#include <MatrixMarket_Tpetra.hpp>
#include <Phalanx_KokkosDeviceTypes.hpp>
#include <Sacado_Fad_SLFad.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

/// @brief All goal symbols are contained in this namespace.
namespace goal {

/// @brief Local ordinal.
using LO = int;

/// @brief Global ordinal.
using GO = long long;

/// @brief Scalar type.
using ST = double;

/// @brief Forward automatic differentiation type.
using FadType = Sacado::Fad::SLFad<ST, Goal_FAD_SIZE>;

/// @brief Teuchos communicator.
using Comm = Teuchos::Comm<int>;

/// @brief Kokkos node.
using KNode =  Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>;

/// @brief Tpetra map.
using Map = Tpetra::Map<LO, GO, KNode>;

/// @brief Tpetra compressed row storage graph.
using Graph = Tpetra::CrsGraph<LO, GO, KNode>;

/// @brief Tpetra importer.
using Import = Tpetra::Import<LO, GO, KNode>;

/// @brief Tpetra exporter.
using Export = Tpetra::Export<LO, GO, KNode>;

/// @brief Tpetra vector.
using Vector = Tpetra::Vector<ST, LO, GO, KNode>;

/// @brief Tpetra multi-vector.
using MultiVector = Tpetra::MultiVector<ST, LO, GO, KNode>;

/// @brief Tpetra compressed row storage graph.
using Matrix = Tpetra::CrsMatrix<ST, LO, GO, KNode>;

/// @brief Tpetra matrix market writer.
using MMWriter = Tpetra::MatrixMarket::Writer<Matrix>;

} // end namespace goal

#endif
