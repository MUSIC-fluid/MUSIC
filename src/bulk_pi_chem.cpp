#include "bulk_pi_chem.h"

namespace {
  BulkPiChem* g_inst = nullptr;

  BulkPiChem* disabled() {
    static ChemBulkConfig cfg;            // chem_bulk_on=0
    static BulkPiChem inst(cfg, nullptr); // safe: PiChem(...) => 0.0
    return &inst;
  }
}

namespace ChemBulk {
  const BulkPiChem& get() { return g_inst ? *g_inst : *disabled(); }
  void set(BulkPiChem* inst) { g_inst = inst; }
}
