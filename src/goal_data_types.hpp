#ifndef goal_data_types_hpp
#define goal_data_types_hpp

#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Core.hpp>

#include "goal_scalar_types.hpp"

namespace goal {

using LO = int;
using GO = long long;
using Comm = Teuchos::Comm<int>;
using MapT = Tpetra::Map<LO, GO>;
using GraphT = Tpetra::CrsGraph<LO, GO>;
using ImportT = Tpetra::Import<LO, GO>;
using ExportT = Tpetra::Export<LO, GO>;
using VectorT = Tpetra::Vector<ST, LO, GO>;
using MultiVectorT = Tpetra::MultiVector<ST, LO, GO>;
using MatrixT = Tpetra::CrsMatrix<ST, LO, GO>;
using MMWriterT = Tpetra::MatrixMarket::Writer<MatrixT>;

}

#endif
