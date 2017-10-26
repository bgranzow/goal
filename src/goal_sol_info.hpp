#ifndef goal_sol_info_hpp
#define goal_sol_info_hpp

#include "goal_data_types.hpp"

namespace goal {

using Teuchos::RCP;

class Disc;

struct LinearObj {
  RCP<VectorT> R;
  RCP<VectorT> dMdu;
  RCP<MatrixT> dRdu;
};

class SolInfo {
  public:
    SolInfo(Disc* d);
    ~SolInfo();
    Disc* get_disc();
    void gather_R();
    void gather_dMdu();
    void gather_dRdu();
    void gather_all();
    void zero_R();
    void zero_dMdu();
    void zero_dRdu();
    void zero_all();
    void resume_fill();
    void complete_fill();
    LinearObj* owned;
    LinearObj* ghost;
  private:
    Disc* disc;
    RCP<ImportT> importer;
    RCP<ExportT> exporter;
};

SolInfo* create_sol_info(Disc* d);
void destroy_sol_info(SolInfo* s);

}

#endif
