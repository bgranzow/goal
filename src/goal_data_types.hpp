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

/// @brief All goal symbols are contained in this namespace
namespace goal {

using LO = int;
using GO = long long;
using ST = double;
using FadType = Sacado::Fad::SLFad<ST, Goal_FAD_SIZE>;
using Comm = Teuchos::Comm<int>;
using KNode =  Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>;
using Map = Tpetra::Map<LO, GO, KNode>;
using Graph = Tpetra::CrsGraph<LO, GO, KNode>;
using Import = Tpetra::Import<LO, GO, KNode>;
using Export = Tpetra::Export<LO, GO, KNode>;
using Vector = Tpetra::Vector<ST, LO, GO, KNode>;
using MultiVector = Tpetra::MultiVector<ST, LO, GO, KNode>;
using Matrix = Tpetra::CrsMatrix<ST, LO, GO, KNode>;
using MMWriter = Tpetra::MatrixMarket::Writer<Matrix>;

} // end namespace goal

#endif
