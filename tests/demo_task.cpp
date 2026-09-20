#include <ParseLQCData/Registry.h>
#include <iostream>

// 基础无参任务
void run_demo() {
    std::cout << "Hello from demo task! The registry works!" << std::endl;
}
REGISTER_TASK("demo", run_demo, "A basic demo task without parameters");

// 标准带参任务
int run_greet(std::span<const std::string_view> args) {
    if (args.empty()) {
        std::cout << "Hello, World! (Pass names as arguments to greet individuals)" << std::endl;
        return 0;
    }
    for (std::string_view name : args) {
        std::cout << "Hello, " << name << "!" << std::endl;
    }
    return 0;
}
REGISTER_TASK("greet", run_greet, "Greet one or more names passed via arguments");
