#include "config.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

namespace sayram2d::app {
namespace {

std::string trim(const std::string& value) {
  const auto begin = value.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) return {};
  const auto end = value.find_last_not_of(" \t\r\n");
  return value.substr(begin, end - begin + 1);
}

[[noreturn]] void fail(const std::filesystem::path& path, int line,
                       const std::string& message) {
  throw std::runtime_error(path.string() + ":" + std::to_string(line) + ": " +
                           message);
}

double parse_double(const std::string& text, const std::filesystem::path& path,
                    int line, const std::string& key) {
  char* end = nullptr;
  errno = 0;
  const double value = std::strtod(text.c_str(), &end);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE ||
      !std::isfinite(value)) {
    fail(path, line, "invalid number for " + key + ": '" + text + "'");
  }
  return value;
}

int parse_int(const std::string& text, const std::filesystem::path& path,
              int line, const std::string& key) {
  char* end = nullptr;
  errno = 0;
  const long value = std::strtol(text.c_str(), &end, 10);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE || value < 0 ||
      value > std::numeric_limits<int>::max()) {
    fail(path, line, "invalid integer for " + key + ": '" + text + "'");
  }
  return static_cast<int>(value);
}

bool valid_run_id(const std::string& run_id) {
  if (run_id.empty()) return false;
  return std::all_of(run_id.begin(), run_id.end(), [](unsigned char c) {
    return std::isalnum(c) != 0 || c == '_' || c == '-';
  });
}

}  // namespace

Config load_config(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open config file: " + path.string());

  std::map<std::string, std::map<std::string, std::string>> values;
  std::map<std::string, int> lines;
  std::string section;
  std::string line;
  int line_number = 0;

  while (std::getline(input, line)) {
    ++line_number;
    const auto comment = line.find_first_of("#;");
    if (comment != std::string::npos) line.resize(comment);
    line = trim(line);
    if (line.empty()) continue;

    if (line.front() == '[') {
      if (line.back() != ']') fail(path, line_number, "malformed section header");
      section = trim(line.substr(1, line.size() - 2));
      if (section.empty()) fail(path, line_number, "empty section name");
      if (section != "run" && section != "input" && section != "time" &&
          section != "output") {
        fail(path, line_number, "unknown section [" + section + "]");
      }
      continue;
    }

    const auto equals = line.find('=');
    if (equals == std::string::npos) {
      fail(path, line_number, "expected 'key = value'");
    }
    if (section.empty()) fail(path, line_number, "key outside a section");
    const std::string key = trim(line.substr(0, equals));
    const std::string value = trim(line.substr(equals + 1));
    if (key.empty() || value.empty()) {
      fail(path, line_number, "empty key or value");
    }

    const bool known =
        (section == "run" && (key == "run_id" || key == "output_dir")) ||
        (section == "input" && key == "data_file") ||
        (section == "time" && (key == "dt" || key == "nsteps")) ||
        (section == "output" && key == "n_frames");
    if (!known) fail(path, line_number, "unknown key [" + section + "] " + key);

    const std::string qualified = section + "." + key;
    if (values[section].count(key) != 0) {
      fail(path, line_number, "duplicate key " + qualified);
    }
    values[section][key] = value;
    lines[qualified] = line_number;
  }

  auto get = [&](const std::string& sec, const std::string& key,
                 bool required) -> const std::string* {
    const auto section_it = values.find(sec);
    if (section_it != values.end()) {
      const auto key_it = section_it->second.find(key);
      if (key_it != section_it->second.end()) return &key_it->second;
    }
    if (required) {
      throw std::runtime_error(path.string() + ": missing required key [" + sec +
                               "] " + key);
    }
    return nullptr;
  };

  Config config;
  const std::string* run_id = get("run", "run_id", true);
  config.run_id = *run_id;
  if (!valid_run_id(config.run_id)) {
    throw std::runtime_error(
        path.string() +
        ": run_id must contain only letters, digits, '_' or '-'");
  }

  if (const std::string* output_dir = get("run", "output_dir", false)) {
    config.output_dir = *output_dir;
  }
  if (config.output_dir.empty()) {
    throw std::runtime_error(path.string() + ": output_dir must not be empty");
  }

  config.data_file = *get("input", "data_file", true);
  if (config.data_file.empty()) {
    throw std::runtime_error(path.string() + ": data_file must not be empty");
  }

  const std::string* dt_text = get("time", "dt", true);
  config.dt = parse_double(*dt_text, path, lines["time.dt"], "time.dt");
  if (!(config.dt > 0.0)) {
    throw std::runtime_error(path.string() + ": time.dt must be positive");
  }

  const std::string* nsteps_text = get("time", "nsteps", true);
  config.nsteps =
      parse_int(*nsteps_text, path, lines["time.nsteps"], "time.nsteps");
  if (config.nsteps <= 0) {
    throw std::runtime_error(path.string() + ": time.nsteps must be positive");
  }
  if (config.nsteps == std::numeric_limits<int>::max()) {
    throw std::runtime_error(path.string() + ": time.nsteps is too large");
  }

  if (const std::string* n_frames = get("output", "n_frames", false)) {
    config.n_frames =
        parse_int(*n_frames, path, lines["output.n_frames"], "output.n_frames");
  } else {
    config.n_frames = config.nsteps + 1;
  }
  if (config.n_frames < 2) {
    throw std::runtime_error(path.string() + ": output.n_frames must be >= 2");
  }
  if (config.nsteps % (config.n_frames - 1) != 0) {
    throw std::runtime_error(
        path.string() + ": (n_frames - 1) must divide nsteps exactly");
  }

  return config;
}

}  // namespace sayram2d::app
