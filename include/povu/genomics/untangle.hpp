#ifndef POVU_UNTANGLE_HPP
#define POVU_UNTANGLE_HPP

#include <string_view> // for string_view

#include "povu/genomics/allele.hpp" // for Exp

namespace povu::genomics::untangle
{
inline constexpr std::string_view MODULE = "povu::genomics::untangle";
namespace pga = povu::genomics::allele;

void untangle_ref_walks(pga::Exp &rt);
} // namespace povu::genomics::untangle

#endif // POVU_UNTANGLE_HPP
