#include "./utils.hpp"

namespace povu::utils {
char complement(char nucleotide) {
  switch (nucleotide) {
  case 'A': return 'T';
  case 'T': return 'A';
  case 'C': return 'G';
  case 'G': return 'C';
  // You might want to handle cases like 'N' for unknown nucleotides
  default: return nucleotide;
  }
}

std::string reverse_complement(const std::string& sequence) {
  std::string rc_sequence;

  for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
    rc_sequence += complement(*it);
  }

  return rc_sequence;
}

bool is_numeric_string(const std::string &s) {
#if __cplusplus >= 202002L
  return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
#else
  return !s.empty() && std::all_of(s.begin(), s.end(), [](char c) {
    return std::isdigit(static_cast<unsigned char>(c));
  });
#endif
}

std::string today() {
    // Get the current time
    std::time_t currentTime = std::time(nullptr);

    // Convert the current time to a local time structure
    std::tm* localTime = std::localtime(&currentTime);

    // Extract date and time components
    int year = localTime->tm_year + 1900; // Years since 1900
    int month = localTime->tm_mon + 1;    // Months are 0-based
    int day = localTime->tm_mday;         // Day of the month

    // Create a stringstream to store the formatted date and time
    std::stringstream dateTimeStream;
    dateTimeStream << year
                   << std::setw(2) << std::setfill('0') << month
                   << std::setw(2) << std::setfill('0') << day;

    return dateTimeStream.str();
}

void report_time(std::ostream& os, const std::string &fn_name, const std::string &action, std::chrono::duration<double> period) {
  os << pv_cmp::format("{} INFO Time spent by {}: {:.2f} sec\n", fn_name, action, period.count());
}

std::vector<std::string> immutable_erase(std::vector<std::string> v, std::size_t idx) {
  v.erase(v.begin()+idx, v.begin()+idx+1);
  return v;
}


std::size_t to_bidirected_idx(std::size_t x, bool has_dummy) {
  // we added 1 because we added a dummy start node before bi-edging
  if (has_dummy) { --x; }
  return x % 2 == 0 ? ((x + 2) / 2) - 1 : ((x + 1) / 2) - 1;
}


std::pair<std::size_t, std::size_t> frm_bidirected_idx(std::size_t x, bool has_dummy) {
  auto foo = [&](std::size_t i) { return (2*(i+1)) + (has_dummy ?  1 : 0 ); };
  return { foo(x) - 2, foo(x) - 1 };
}

void split(const std::string &line, char sep, std::vector<std::string> *tokens) {

  // split line at sep and store results in tokens
  std::size_t start = 0;
  std::size_t end = line.find(sep);

  while (end != std::string::npos) {
    tokens->push_back(line.substr(start, end - start));
    start = end + 1;
    end = line.find(sep, start);
  }

  tokens->push_back(line.substr(start, end));
}

} // namespace povu::utils
