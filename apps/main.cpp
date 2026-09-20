#include <ParseLQCData/Registry.h>
#include <iostream>
#include <vector>
#include <string_view>

int main(int argc, char* argv[]) {
    const auto& tasks = registry::TaskRegistry::instance().all();

    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " [task_name] [args...]\n\n";
        std::cout << "Available tasks:\n";
        for (const auto& [name, entry] : tasks) {
            std::cout << "  - " << name;
            if (!entry.desc.empty()) {
                std::cout << " : " << entry.desc;
            }
            std::cout << "\n";
        }
        return 0;
    }

    std::string_view task_name = argv[1];
    auto it = tasks.find(std::string(task_name));
    if (it == tasks.end()) {
        std::cerr << "Error: Task '" << task_name << "' not found.\n";
        return 1;
    }

    std::vector<std::string_view> args;
    args.reserve(argc > 2 ? argc - 2 : 0);
    for (int i = 2; i < argc; ++i) {
        args.emplace_back(argv[i]);
    }

    return it->second.handler(args);
}
