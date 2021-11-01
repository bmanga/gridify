#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>
#include <string>

namespace {
struct config {
#define X(type, name, default_, ...) type name = {default_};
  ALL_SETTINGS
#undef X
};

config parse_yaml_settings(const std::string &file)
{
  config c;
  YAML::Node y = YAML::LoadFile(file);
#define X(type, name, ...)             \
  if (y[#name]) {                 \
    c.name = y[#name].as<type>(); \
  }

  ALL_SETTINGS

#undef X
  return c;
}

template <class Fn>
config parse_cmd_line_settings(const char *name, const char * desc, int argc, char *argv[],
    const Fn &fn)
{
  auto opts = cxxopts::Options(name, desc);

  opts.add_options()
#define X(type, name, default_, desc) (#name, desc, cxxopts::value<type>()->default_value( \
                                                        std::string(#default_) == "{}" ? "" : #default_                          \
 ))
      ALL_SETTINGS
#undef X
      ("l,load_yaml_defaults", "get option defaults from the specified file", cxxopts::value<std::string>()->default_value(""));

  fn(opts);
  auto parsed_opts = opts.parse(argc, argv);
  config c;

  auto yaml_defaults_file = parsed_opts["load_yaml_defaults"].as<std::string>();
  if (!yaml_defaults_file.empty()) {
    try {
      c = parse_yaml_settings(yaml_defaults_file);
    }
    catch (YAML::BadFile &e) {
      std::cerr << "Error parsing yaml file '" << yaml_defaults_file << "'"
                << std::endl;
      std::cerr << e.what() << std::endl;
      std::exit(-1);
    }
  }

#define X(type, name, ...)                \
  if (parsed_opts.count(#name) != 0) \
    c.name = parsed_opts[#name].as<type>();
  ALL_SETTINGS
#undef X
  return c;
}

void print_configuration(config &config)
{
  fmt::print("Configuration: \n");
#define X(type, name, ...) fmt::print("  -- {}: {}\n", #name, config.name);
  ALL_SETTINGS
#undef X
}

} // namespace