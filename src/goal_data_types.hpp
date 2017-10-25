#ifndef goal_data_types_hpp
#define goal_data_types_hpp

#include <MatrixMarket_Tpetra.hpp>
#include <Phalanx_KokkosDeviceTypes.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>

#include "goal_scalar_types.hpp"

namespace goal {

using LO = int;
using GO = long long;
using Comm = Teuchos::Comm<int>;
using KNode =  Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>;
using MapT = Tpetra::Map<LO, GO, KNode>;
using GraphT = Tpetra::CrsGraph<LO, GO, KNode>;
using ImportT = Tpetra::Import<LO, GO, KNode>;
using ExportT = Tpetra::Export<LO, GO, KNode>;
using VectorT = Tpetra::Vector<ST, LO, GO, KNode>;
using MultiVectorT = Tpetra::MultiVector<ST, LO, GO, KNode>;
using MatrixT = Tpetra::CrsMatrix<ST, LO, GO, KNode>;
using MMWriterT = Tpetra::MatrixMarket::Writer<MatrixT>;

}

#endif
