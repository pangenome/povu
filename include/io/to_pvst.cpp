#include "./to_pvst.hpp"

namespace povu::io::to_pvst {

void write_pvst(const pvtr::Tree &bt, const std::string &base_name, const core::config &app_config) {
  // TODO: combine and pass as single arg
  std::string bub_file_name = pv_cmp::format("{}/{}.pvst", std::string{app_config.get_output_dir()}, base_name); // file path and name
  std::ofstream bub_file(bub_file_name);

  if (!bub_file.is_open()) {
    std::cerr << "ERROR: could not open file " << bub_file_name << "\n";
    std::exit(1);
  }

  // writer header line
  // ------------------
  bub_file << constants::PVST_HEADER_SYMBOL << pc::COL_SEP << pc::PVST_VERSION
           << pc::COL_SEP << "." << pc::COL_SEP << "." << "\n";

  // write the rest of the PVST
  // --------------------------
  for (std::size_t i {}; i < bt.vtx_count(); ++i) {

    const pvst::VertexBase &v = bt.get_vertex(i);

    // line identifier
    {
    switch (v.get_type()) {
      case pvst::vt_e::slubble:
        bub_file << pc::PVST_CONCEALED_SYMBOL << pc::COL_SEP;
        break;
      case pvst::vt_e::flubble:
        bub_file << pc::PVST_FLUBBLE_SYMBOL << pc::COL_SEP;
        break;
      case pvst::vt_e::dummy:
        bub_file << pc::PVST_DUMMY_SYMBOL << pc::COL_SEP;
        break;
      case pvst::vt_e::tiny:
        bub_file << pc::PVST_TINY_SYMBOL << pc::COL_SEP;
        break;
      case pvst::vt_e::parallel:
        bub_file << pc::PVST_OVERLAP_SYMBOL << pc::COL_SEP;
        break;
      case pvst::vt_e::smothered:
        bub_file << pc::PVST_SMOTHERED_SYMBOL << pc::COL_SEP;
        break;
      case pvst::vt_e::midi:
        bub_file << pc::PVST_MIDI_SYMBOL << pc::COL_SEP;
        break;
      default:
        std::cerr << "ERROR: unknown vertex type in write_bub: " << v.as_str() << "\n";
        std::exit(1);
    }
    }

    // vertex idx
    bub_file << i << pc::COL_SEP;

    // vertex as str
    bub_file << v.as_str() << pc::COL_SEP;

    // children
    if (bt.is_leaf(i)) {
      bub_file << pc::NO_VALUE << "\n";
    }
    else {
      pu::print_with_comma(bub_file, bt.get_children(i), ',');
      bub_file << "\n";
    }
  }

  bub_file.close();
}
} // namespace povu::io::pvst
