#pragma once

#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <span>
#include <unordered_map>
#include <functional>
#include <type_traits>

namespace registry {

using TaskHandler = std::function<int(std::span<const std::string_view>)>;

class TaskRegistry {
public:
    struct TaskEntry {
        TaskHandler handler;
        std::string desc;
    };

    static TaskRegistry& instance() {
        static TaskRegistry inst;
        return inst;
    }

    void add(std::string name, TaskHandler handler, std::string desc = "") {
        tasks_[std::move(name)] = {std::move(handler), std::move(desc)};
    }

    int execute(std::string_view name, std::span<const std::string_view> args) const {
        auto it = tasks_.find(std::string(name));
        if (it == tasks_.end()) return -1;
        return it->second.handler(args);
    }

    const auto& all() const { return tasks_; }

private:
    std::unordered_map<std::string, TaskEntry> tasks_;
    TaskRegistry() = default;
};

// 单一泛型包装器抹平返回值与参数差异
template <typename F>
TaskHandler make_handler(F&& f) {
    using DecayedF = std::decay_t<F>;
    return [fn = DecayedF(std::forward<F>(f))](std::span<const std::string_view> args) -> int {
        if constexpr (std::is_invocable_r_v<int, DecayedF, std::span<const std::string_view>>) return fn(args);
        else if constexpr (std::is_invocable_v<DecayedF, std::span<const std::string_view>>) { fn(args); return 0; }
        else if constexpr (std::is_invocable_r_v<int, DecayedF>) return fn();
        else if constexpr (std::is_invocable_v<DecayedF>) { fn(); return 0; }
        else static_assert(!sizeof(DecayedF), "Unsupported function signature");
    };
}

struct AutoRegister {
    template <typename F>
    AutoRegister(const std::string& name, F&& f, const std::string& desc = "") {
        TaskRegistry::instance().add(name, make_handler(std::forward<F>(f)), desc);
    }
};

} // namespace registry

#define REGISTER_TASK(name, func, ...) \
    static const ::registry::AutoRegister reg_##func(name, func, ##__VA_ARGS__)
